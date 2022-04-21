library(tidyverse)
library(janitor)
library(spatialTIME)

load("/Volumes/lab_fridley/Coghill/ProstateStudy/PilotStudy/data/TMA_data.RData")
load("/Volumes/lab_fridley/Coghill/ProstateStudy/PilotStudy/data/patient_sample_data.RData")

example_clinical <- patient_data %>% 
  select(patient_id, age, race, sex, status = hiv_status) %>%
  mutate(race = tolower(race)) %>% 
  mutate(status = stringr::str_remove(status, "HIV")) %>% 
  mutate(status = case_when(
    status == "-" ~ "B",
    TRUE ~ "A"
  )) %>%
  left_join(sample_data_panel2 %>% select(image_tag, patient_id, type), 
            by = "patient_id") %>% 
  filter(!is.na(image_tag)) %>% 
  filter(!is.na(patient_id)) %>% 
  filter(type == "tumor") %>%
  mutate(tissue = gsub("Coghill_P2_", "", image_tag)) %>% 
  mutate(tissue = substring(tissue,1,1)) %>% 
  filter(tissue == "P") %>% 
  select(-tissue)


link_up_file <- example_clinical %>% 
  select(patient_id, image_tag) %>% 
  # mutate(deidentified_id = sample(1:1000, replace = FALSE)) %>% 
  mutate(deidentified_sample = str_replace(image_tag, "Coghill_P2_Anal-Invasive-TMA1", "TMA4")) %>% 
  mutate(deidentified_sample = str_replace(deidentified_sample, "Coghill_P2_Prostate-", ""))

set.seed(8675309)
link_up_file$deidentified_id <- sample(1:1000, nrow(link_up_file), replace = FALSE)

example_clinical <- example_clinical %>% 
  left_join(link_up_file) %>% 
  select(-c(patient_id, image_tag, type))

example_summary <- summary_panel2 %>% 
  left_join(link_up_file %>% select(image_tag, deidentified_sample, deidentified_id), 
            by = c("Image Tag" = "image_tag")) %>% 
  select(-`Image Tag`) %>% 
  select(deidentified_id, deidentified_sample, everything()) %>% 
  filter(!is.na(deidentified_id)) %>% 
  filter(!is.na(deidentified_sample)) %>% 
  filter(deidentified_sample %in% example_clinical$deidentified_sample)

example_spatial <- plyr::llply(tma_panel2, function(x){
  y <- x %>% 
    select(-c(path, Image.Location)) %>%
    left_join(link_up_file %>% select(image_tag, deidentified_sample), 
              by = c("image.tag" = "image_tag")) %>% 
    select(-`image.tag`) %>% 
    select(deidentified_sample, everything())
  
  return(y)
})

names_x <- c()
for(i in 1:length(example_spatial)){
  y <- ifelse(is.na(example_spatial[[i]]$deidentified_sample[[1]]), 
              paste0("missing",i),
              example_spatial[[i]]$deidentified_sample[[1]])
  
  names_x  <- c(names_x, y)
}

names(example_spatial) <- names_x

example_samples <- sample(example_clinical$deidentified_sample, 
                         15, replace = FALSE)


example_spatial_small <- example_spatial[example_samples]
example_clinical_small = example_clinical[example_clinical$deidentified_sample %in% example_samples,]
example_sample_small = example_summary[example_summary$deidentified_sample %in% example_samples,]
colnames(example_spatial_small[[1]]) %>% grep("Positive|CD|Opal", ., value = T) %>% grep("\\Nu|Cyt", ., value = T, invert = T)

example_mif = create_mif(clinical_data = example_clinical_small,
                 sample_data = example_sample_small,
                 spatial_list = example_spatial_small,
                 patient_id = "deidentified_id", 
                 sample_id = "deidentified_sample")
example_mif = compute_metrics(mif = example_mif,
                      mnames = "CD3..Opal.570..Positive",
                      r_range = seq(0, 100, 50),
                      num_permutations = 10,
                      edge_correction = c("translation"),
                      method = c("K"),
                      k_trans = "none",
                      keep_perm_dis = F,
                      workers = 10,
                      overwrite = T,
                      xloc = NULL,
                      yloc = NULL,
                      exhaustive = T)

saveRDS(example_mif,
        file = "/Volumes/lab_fridley/IHC/SMI_2022_Workshop/example_mif_4metrics.rds")

picked_marker_dat = example_mif$derived$univariate_Count %>%
  group_by(deidentified_sample, Marker, r) %>%
  summarise(across(`Theoretical CSR`:`Degree of Clustering Theoretical`, .fns = function(x){mean(x, na.rm=T)})) %>%
  filter(r == 50, Marker %in% c("CD3..Opal.570..Positive"))
new_assignments = picked_marker_dat %>%
  arrange(`Degree of Clustering Theoretical`) %>%
  ungroup() %>%
  mutate(status = c(rep("A", 5), sample(c("A", "B"), 5, replace = T), rep("B", 5))) %>%
  select(deidentified_sample, status)
new_example_clinical_small = example_clinical_small %>%
  select(-status) %>%
  full_join(new_assignments)

write.csv(new_example_clinical_small, 
          file = "/Volumes/Lab_Fridley/IHC/SMI_2022_Workshop/data/deidentified_clinical.csv",
          row.names = FALSE)
write.table(example_sample_small,
            file = "/Volumes/Lab_Fridley/IHC/SMI_2022_Workshop/data/deidentified_summary.csv",
            row.names = FALSE)
save(#example_clinical[example_clinical$deidentified_sample %in% example_samples,],
     # example_summary[example_summary$deidentified_sample %in% example_samples,],
     example_spatial_small, 
     file = "/Volumes/Lab_Fridley/IHC/SMI_2022_Workshop/data/deidentified_example.RData")

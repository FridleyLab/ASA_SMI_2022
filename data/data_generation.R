library(tidyverse)
library(janitor)

load("/Volumes/Lab_Fridley/Coghill/ProstateStudy/PilotStudy/data/TMA_data.RData")
load("/Volumes/Lab_Fridley/Coghill/ProstateStudy/PilotStudy/data/patient_sample_data.RData")

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

write.csv(example_clinical[example_clinical$deidentified_sample %in% example_samples,], 
          file = "/Volumes/Lab_Fridley/IHC/SMI_2022_Workshop/data/deidentified_clinical.csv",
          row.names = FALSE)
write.table(example_summary[example_summary$deidentified_sample %in% example_samples,],
            file = "/Volumes/Lab_Fridley/IHC/SMI_2022_Workshop/data/deidentified_summary.csv",
            row.names = FALSE)
save(#example_clinical[example_clinical$deidentified_sample %in% example_samples,],
     # example_summary[example_summary$deidentified_sample %in% example_samples,],
     example_spatial_small, 
     file = "/Volumes/Lab_Fridley/IHC/SMI_2022_Workshop/data/deidentified_example.RData")

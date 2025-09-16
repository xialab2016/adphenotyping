load("/d_label_new.RData")
load("/dat_label.RData")
load("/ad_incident_registry.RData")

######## Validation of Chart Review and Registry Labels ########

# Function to create 2x2 table and return metrics
create_2x2_table <- function(predicted, actual, rule_name, time_window) {
  # Create confusion matrix
  confusion_matrix <- table(Predicted = predicted, Actual = actual)
  
  # Ensure it's a 2x2 matrix (in case of missing categories)
  if(nrow(confusion_matrix) == 1 | ncol(confusion_matrix) == 1) {
    full_matrix <- matrix(0, nrow = 2, ncol = 2, 
                          dimnames = list(Predicted = c("0", "1"), Actual = c("0", "1")))
    for(i in rownames(confusion_matrix)) {
      for(j in colnames(confusion_matrix)) {
        full_matrix[i, j] <- confusion_matrix[i, j]
      }
    }
    confusion_matrix <- full_matrix
  }
  
  # Extract values
  tn <- confusion_matrix[1,1]  # True Negative
  fp <- confusion_matrix[2,1]  # False Positive  
  fn <- confusion_matrix[1,2]  # False Negative
  tp <- confusion_matrix[2,2]  # True Positive
  
  # Return results as data frame
  return(data.frame(
    Rule = rule_name,
    Time_Window = time_window,
    TP = tp, 
    TN = tn, 
    FP = fp, 
    FN = fn
  ))
}

results_confusion_mat <- rbind(
  create_2x2_table(d.label.new$pred.outcome.1, d.label.new$AD_probable_possible, "UPMC", "Chart Review=AD"),
  create_2x2_table(dat.label$clust.90, dat.label$isAD, "MGB", "Chart Review=AD"),
  create_2x2_table(ad_incident$pred.outcome.1, ad_incident$dx_bin, "UPMC", "Registry=AD"))

results_gold_std_metrics <- results_confusion_mat %>%
  mutate(
    # Calculate performance metrics
    Sensitivity = round(TP / (TP + FN), 3),
    Specificity = round(TN / (TN + FP), 3),
    PPV = round(TP / (TP + FP), 3),
    NPV = round(TN / (TN + FN), 3),
    Accuracy = round((TP + TN) / (TP + TN + FP + FN), 3),
    
    # Handle any potential division by zero
    Sensitivity = ifelse(is.nan(Sensitivity) | is.infinite(Sensitivity), NA, Sensitivity),
    Specificity = ifelse(is.nan(Specificity) | is.infinite(Specificity), NA, Specificity),
    PPV = ifelse(is.nan(PPV) | is.infinite(PPV), NA, PPV),
    NPV = ifelse(is.nan(NPV) | is.infinite(NPV), NA, NPV)
  )

######### Chart Review and ADRC registry Demographics ###############
library(dplyr)
library(data.table)

###### UPMC #########
# Codified data
data.cod = read.csv("/UPMC_AD_2004_to_2023_Codified_rolled_up_data_2023-10-06.csv")

# 1. Demographics
demo_upmc <- read.csv("/UPMC demographics imputed race ethnicity separate.csv")
demo_datefix <- fread("/ad_patient_demographics-datefix.dsv")
demo_upmc$patient_num <- as.character(demo_upmc$PATIENT_STUDY_ID)
## Recoding gender
demo_upmc$GENDER_TITLE <- as.factor(demo_upmc$GENDER_TITLE)
demo_upmc$gender <- case_when(demo_upmc$GENDER_TITLE == "FEMALE" ~ "Female",
                              demo_upmc$GENDER_TITLE == "MALE" ~ "Male",
                              demo_upmc$GENDER_TITLE == "UNKNOWN" ~ NA,
                              demo_upmc$GENDER_TITLE == "" ~ NA,
                              TRUE ~ NA)
## Recoding race and ethnicity
demo_upmc$race_eth_final <- case_when(demo_upmc$race %in% c("White") & demo_upmc$ethnicity == "Non Hispanic" ~ "Non Hispanic White",
                                      demo_upmc$race %in%  c("American Indian or Alaska Native", "Asian", "Black", 
                                                             "Native Hawaiian or Other Pacific Islander", "Other", "White") & 
                                        demo_upmc$ethnicity == "Hispanic" ~ "Hispanic or Latino",
                                      demo_upmc$race %in% c("American Indian or Alaska Native", "Asian", "Black", 
                                                            "Native Hawaiian or Other Pacific Islander", "Other") & 
                                        demo_upmc$ethnicity == "Non Hispanic" ~ demo_upmc$race,
                                      TRUE ~ NA)

# 2. Age at AD diagnosis (first PheCode)
library(eeptools)
## Getting AD diagnosis date
diag.date = data.cod %>% filter(feature_id == "PheCode:290.11" | feature_id == "PheCode:290.1" | feature_id == "PheCode:290") %>% 
  group_by(patient_num) %>% 
  mutate(diag_date = as.Date(start_date, "%Y-%m-%d")) %>%
  arrange(diag_date) %>%
  filter(row_number()==1) %>%
  select(patient_num,diag_date) %>%
  mutate(patient_num = as.character(patient_num))
## Getting birth date
df.bd = demo_datefix %>% mutate(patient_num = as.character(PATIENT_STUDY_ID), 
                                birth_date=base::as.Date(`TO_CHAR(BIRTH_DATE,'MM-DD-YYYY')`, format="%m-%d-%Y")) %>%
  select(patient_num, birth_date)
df.bd <- df.bd %>% left_join(diag.date, by = "patient_num")

# For chart review labels
d.label.new = d.label.new %>% 
  left_join(df.bd, by="patient_num") %>% 
  left_join(demo_upmc[,c("patient_num", "gender", "race_eth_final")], by = "patient_num") %>% 
  filter(!is.na(diag_date)) %>% mutate(age = round(age_calc(birth_date, diag_date)/12))
d.label.new <- d.label.new[,c("patient_num", "AD_probable_possible", "pred.outcome.1", "birth_date", "diag_date", "gender", "race_eth_final", "age")]
save(d.label.new, file = "/d_label_new_for_demographics.RData")
  
# For ADRC registry
ad_incident = ad_incident %>% 
  left_join(df.bd, by="patient_num") %>% 
  left_join(demo_upmc[,c("patient_num", "gender", "race_eth_final")], by = "patient_num") 
ad_incident_linked <- ad_incident %>% filter(!is.na(diag_date)) %>% group_by(patient_num) %>% slice_head(n = 1) %>% ungroup()
ad_incident_linked <- ad_incident_linked %>% filter(!is.na(diag_date)) %>% mutate(age = round(age_calc(birth_date, diag_date)/12))

###### MGB #########

# 1. Demographics
demo_mgb <- read.csv("/MGB demographics imputed race ethnicity separate.csv")[,-1]
demo_mgb$patient_num <- as.character(demo_mgb$PatientNum)
demo_mgb$gender <- case_when(demo_mgb$Gender_Legal_Sex == "Female" ~ "Female",
                             demo_mgb$Gender_Legal_Sex == "Male" ~ "Male",
                             demo_mgb$Gender_Legal_Sex == "Unknown-U" ~ NA,
                             TRUE ~ NA)
demo_mgb$race_eth_final <- case_when(demo_mgb$race %in% c("White") & demo_mgb$ethnicity == "Non Hispanic" ~ "Non Hispanic White",
                                     demo_mgb$race %in%  c("American Indian or Alaska Native", "Asian", "Black", 
                                                           "Native Hawaiian or Other Pacific Islander", "Other", "White") & 
                                       demo_mgb$ethnicity == "Hispanic" ~ "Hispanic or Latino",
                                     demo_mgb$race %in% c("American Indian or Alaska Native", "Asian", "Black", 
                                                          "Native Hawaiian or Other Pacific Islander", "Other") & 
                                       demo_mgb$ethnicity == "Non Hispanic" ~ demo_mgb$race,
                                     TRUE ~ NA)

# 2. Age at AD diagnosis
load("/combined_daily.RData") 
map_file <- read.csv("/codebook_2024-06-12.csv")

## Getting AD diagnosis date
dementia_ad_phecode <- map_file %>% filter(map_file$Group_Code == "PheCode:290.11" | Group_Code == "PheCode:290.1" | Group_Code == "PheCode:290")
diag_df <- combined_data %>%
  semi_join(dementia_ad_phecode, by = "Local_Code")
diag.date = diag_df %>% 
  group_by(PatientNum) %>% 
  mutate(diag_date = as.Date(Start_date, "%Y-%m-%d")) %>%
  arrange(diag_date) %>%
  filter(row_number()==1) %>%
  mutate(patient_num = as.character(PatientNum)) %>%
  ungroup() %>%
  select(patient_num,diag_date) 
## Getting birth date
df.bd = demo_mgb %>% mutate(patient_num = as.character(PatientNum), 
                            birth_date=base::as.Date(Date_of_Birth, format="%Y-%m-%d")) %>%
  select(patient_num, birth_date)
df.bd <- df.bd %>% left_join(diag.date, by = "patient_num")

# For chart review labels
dat.label = dat.label %>% 
  left_join(df.bd, by="patient_num") %>% 
  left_join(demo_mgb[,c("patient_num", "gender", "race_eth_final")], by = "patient_num") %>% 
  filter(!is.na(diag_date)) %>% mutate(age = round(age_calc(birth_date, diag_date)/12))
colnames(dat.label)[2] <- "AD_probable_possible"
colnames(dat.label)[3] <- "pred.outcome.1"
save(dat.label, file = "/dat_label_for_demographics.RData")

###### Combined #########
library(gtsummary)
library(tidyr)

# Combining demographics data
d.label.new$site <- "UPMC"
dat.label$site <- "MGB"
label_com <- rbind(d.label.new, dat.label)

# Overall
table_1 <- select(label_com, c(age, gender, race_eth_final)) %>%
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2,
              missing_text = "Missing",
              label = list(age ~ "Age at AD diagnosis", gender ~ "Gender", race_eth_final ~ "Race and ethnicity")) %>%
  bold_labels() %>%
  modify_header(stat_0 ~ "**Overall**")

# By site
table_2 <- select(label_com, c(age, gender, race_eth_final, site)) %>%
  tbl_summary(
    by = site,
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    missing_text = "Missing",
    label = list(age ~ "Age at AD diagnosis", gender ~ "Gender", race_eth_final ~ "Race and ethnicity")) %>%
  add_p() %>%
  bold_labels() 

# ADRC
table_3 <- select(ad_incident, c(gender, race_eth_final)) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    missing_text = "Missing",
    label = list(gender ~ "Gender", race_eth_final ~ "Race and ethnicity")) %>%
  bold_labels() 

table_4 <- select(ad_incident_linked, c(age)) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 2,
    missing_text = "Missing",
    label = list(age ~ "Age at AD diagnosis")) %>%
  bold_labels() 

# Merge all tables into one
combined_table <- tbl_merge(
  tbls = list(table_1, table_2),
  tab_spanner = c("**Overall**", "**By Site**")
) %>%
  modify_caption("**Table 1. Baseline Characteristics**")

combined_table_csv <- combined_table %>%
  as_tibble()
write.csv(combined_table_csv, "/baseline_characteristics_chartlabels_registry.csv")







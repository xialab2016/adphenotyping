library(lubridate)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(ggplot2)
library(readxl)

##########################################################
# Reading in codified data, NLP data, and mapping file
## combined_daily is the codified data file
## combined_nlp_1 and combined_nlp_2 are the NLP data files
## map_file is used for common ontology codes
##########################################################

load("/combined_daily.RData") 
combined_nlp_1 <- read.csv("/combined_part1.csv")
combined_nlp_2 <- read.csv("/combined_part2.csv")
map_file <- read.csv("/codebook_2024-06-12.csv")

# Output
path.out = "/your_path/"

##########################################################
# Nursing home admission
## nh.cui is the list of nursing home CUIs
##########################################################

## Read in nursing home admission codes
## Revised definition
## Map the nursing home admission codes to common ontology and local codes
nh_def_old_revised <- read_excel("/NH_codes_2023_new_def.xlsx")
nh_def_old_revised$Common_Ontology_Code <- ifelse(grepl("CM", nh_def_old_revised$concept_cd), gsub("CM", "", nh_def_old_revised$concept_cd), nh_def_old_revised$concept_cd)
nh_def_old_revised <- nh_def_old_revised %>% left_join(map_file, by = "Common_Ontology_Code")

## Obtain list of local nursing home admission codes
nh_revised_local <- nh_def_old_revised[,c("Local_Code")] %>% filter(!is.na(Local_Code))

# Get time to nursing home admission
## a. codified
nh_revised_df <- combined_data %>%
  semi_join(nh_revised_local, by = "Local_Code")
nh.revised.date = nh_revised_df %>% mutate(PatientNum = as.character(PatientNum)) %>%
  group_by(PatientNum) %>%
  mutate(nh_date = as.Date(Start_date, "%Y-%m-%d")) %>%
  arrange(nh_date) %>%
  filter(row_number()==1) %>%
  select(PatientNum, nh_date)
save(nh_revised_df, nh.revised.date, file=paste0(path.out, "nursing_home_admission_revised.RData"))

## b. NLP
## Original list of CUIs
nh.cuis = read.table("/nursing_home_cuis_from_dict.dsv", sep="|")$V2
nh_cui_1 <- combined_nlp_1 %>%
  filter(cui %in% nh.cuis) %>%
  distinct(patientnum, date, cui) %>%
  mutate(patientnum = as.character(patientnum)) %>%
  mutate(date = as.Date(date)) 
nh_cui_2 <- combined_nlp_2 %>%
  filter(cui %in% nh.cuis) %>%
  distinct(patientnum, date, cui) %>%
  mutate(patientnum = as.character(patientnum)) %>%
  mutate(date = as.Date(date)) 
nh_cui <- rbind(nh_cui_1, nh_cui_2)
nh.cui.date = nh_cui %>%
  group_by(patientnum) %>%
  mutate(nh_cui_date = as.Date(date, "%Y-%m-%d")) %>%
  arrange(nh_cui_date) %>%
  filter(row_number()==1) %>%
  select(patientnum, nh_cui_date)  %>% 
  mutate(patientnum = as.character(patientnum))
save(nh_cui, nh.cui.date, file=paste0(path.out, "nursing_home_admission_narrative_old_date.RData"))

##########################################################
# MGB imputed demographics
##########################################################

demo <- read.csv("/MGB demographics imputed race ethnicity separate.csv")
## Recoding gender
demo$female <- case_when(demo$Gender_Legal_Sex == "Female" ~ 1,
                         demo$Gender_Legal_Sex == "Male" ~ 0,
                         demo$Gender_Legal_Sex == "Unknown-U" ~ NA,
                         TRUE ~ NA)
d.gender = demo[,c("PatientNum", "female")]
d.gender <- d.gender[complete.cases(d.gender[, c("PatientNum", "female")])]
## Recoding race and ethnicity
demo$race_eth_nhw <- case_when(demo$race %in% c("White") & demo$ethnicity == "Non Hispanic" ~ 1,
                               TRUE ~ NA)
demo$race_eth_other <- case_when(demo$race %in% c("American Indian or Alaska Native", "Asian", "Black", 
                                                  "Native Hawaiian or Other Pacific Islander", "Other") & demo$ethnicity == "Non Hispanic" ~ 1,
                                 demo$race %in% c("American Indian or Alaska Native", "Asian", "Black", 
                                                  "Native Hawaiian or Other Pacific Islander", "Other", "White") & demo$ethnicity == "Hispanic" ~ 2,
                                 TRUE ~ NA)
demo$white <- case_when(demo$race_eth_nhw == 1 ~ 1, demo$race_eth_other %in% c(1, 2) ~ 0, is.na(demo$race_eth_nhw) | is.na(demo$race_eth_other) ~ NA_real_)
d.race = demo[,c("PatientNum", "white")]

## Mortality date
## Converting dates of birth and death into standard format
demo$Date_of_Birth <- as.Date(demo$Date_of_Birth, format = "%Y-%M-%D")
demo$Date_of_Death <- as.Date(demo$Date_of_Death, format = "%m/%d/%Y")
mortality.date = demo %>% mutate(patient_num = as.character(demo$PatientNum), 
                                 death_date = as.Date(Date_of_Death, "%Y-%M-%D")) %>%
  select(patient_num, death_date)
save(mortality.date, file=paste0(path.out, "mortality_date.RData"))

## AD diagnosis date
ad_phecode <- map_file %>% filter(map_file$Group_Code == "PheCode:290.11")
diag_df <- combined_data %>%
  semi_join(ad_phecode, by = "Local_Code")
diag.date = diag_df %>% 
  group_by(PatientNum) %>% 
  mutate(diag_date = as.Date(Start_date, "%Y-%m-%d")) %>%
  arrange(diag_date) %>%
  filter(row_number()==1) %>%
  select(PatientNum,diag_date) %>%
  mutate(PatientNum = as.character(PatientNum))
save(diag.date, file=paste0(path.out, "AD_diagnosis_date.RData"))

# Baseline co-morbidity - get raw icd
## Merge codified data with mapping file
data.icd <- combined_data %>% left_join(map_file, by = c("Local_Code", "Local_Code_Description"))
data.icd <- data.icd[,c("PatientNum", "Common_Ontology_Code", "Start_date")] %>% 
  mutate(PatientNum = as.character(PatientNum))
colnames(data.icd) = c("patient_num", "code", "date")
save(data.icd, file=paste0(path.out, "ICD.RData"))

# AD-related medications
med.list = c("RXNORM:135447", "RXNORM:183379", "RXNORM:4637", "RXNORM:6719", 
             "RXNORM:2557217", "RXNORM:2626143", "RXNORM:1430990")
med.date = data.icd %>% filter(code %in% med.list) %>%
  group_by(patient_num) %>%
  mutate(med_date = as.Date(date, "%Y-%m-%d")) %>%
  arrange(med_date) %>%
  filter(row_number()==1) %>%
  select(patient_num, med_date)
head(med.date)
save(med.date, file=paste0(path.out, "medication_date.RData"))

#################################
## Merge data
#################################

# Get AD diagnosis status as predicted by KOMAP
pred <- read.csv("/prediction_imputed_model1.csv")
pred <- pred[,c("patient_num", "clust.90")]
colnames(d.race) = c("patient_num", "white")
colnames(d.gender) = c("patient_num", "female")
colnames(diag.date) = c("patient_num", "diag_date")
colnames(nh.old.date) = c("patient_num", "nh_date")
colnames(nh.revised.date) = c("patient_num", "nh_date")
colnames(nh.cui.date) = c("patient_num", "nh_cui_date")
colnames(med.date) = c("patient_num", "med_date")
d.race$patient_num = as.character(d.race$patient_num)
d.gender$patient_num = as.character(d.gender$patient_num)
diag.date$patient_num = as.character(diag.date$patient_num)
nh.old.date$patient_num = as.character(nh.old.date$patient_num)
nh.revised.date$patient_num = as.character(nh.revised.date$patient_num)
nh.cui.date$patient_num = as.character(nh.cui.date$patient_num)
med.date$patient_num = as.character(med.date$patient_num)

df_revised = pred %>% mutate(patient_num = as.character(patient_num)) %>%
  filter(clust.90==1) %>% left_join(d.race, by="patient_num")%>% 
  left_join(d.gender, by="patient_num") %>%
  left_join(diag.date, by="patient_num") %>%
  left_join(nh.revised.date, by="patient_num") %>%
  left_join(nh.cui.date, by="patient_num") %>%
  left_join(mortality.date, by="patient_num") %>%
  left_join(med.date, by="patient_num") %>%
  mutate(nh_code_month = mondf(diag_date, nh_date),
         nh_cui_month = mondf(diag_date, nh_cui_date),
         death_month = mondf(diag_date, death_date),
         med_month = mondf(diag_date, med_date))
df_revised <- df_revised %>% filter(is.na(df_revised$female) == FALSE) 
df_revised <- df_revised %>% filter(is.na(df_revised$white) == FALSE) 
save(df_revised, file=paste0(path.out, "survival_analysis_df_revised_1.RData"))

#################################
## Analysis
#################################

monnb = function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01")); lt$year*12 + lt$mon } 
mondf = function(d1, d2) { monnb(d2) - monnb(d1) }

# Time to death
library(survival)
df_revised = df_revised[!is.na(df_revised$diag_date),]
df_revised$death = ifelse(is.na(df_revised$death_month), 0, 1)
df_revised$death_month = ifelse(is.na(df_revised$death_month), mondf(df_revised$diag_date, "2022-12-31"), df_revised$death_month)

# Time to nursing home admission
df_revised$nh_code = ifelse(is.na(df_revised$nh_code_month), 0, 1)
df_revised$nh_cui = ifelse(is.na(df_revised$nh_cui_month), 0, 1)
df_revised$nh = ifelse(df_revised$nh_code+df_revised$nh_cui>=1, 1, 0)

df_revised$nh_code_month = ifelse(is.na(df_revised$nh_code_month), mondf(df_revised$diag_date, "2022-12-31"), df_revised$nh_code_month)
df_revised$nh_cui_month = ifelse(is.na(df_revised$nh_cui_month), mondf(df_revised$diag_date, "2022-12-31"), df_revised$nh_cui_month)
df_revised$nh_month = pmin(df_revised$nh_code_month, df_revised$nh_cui_month)

df.nh = df_revised
df.nh$event_code_month = pmin(df.nh$nh_code_month, df.nh$death_month)
df.nh$event_month = pmin(df.nh$nh_month, df.nh$death_month)

df.nh$event_code = ifelse(df.nh$nh_code + df.nh$death == 0, 0, 
                          ifelse(df.nh$nh_code_month<=df.nh$death_month, 1, 2))
df.nh$event_code = factor(df.nh$event_code, levels = c(0,1,2), labels = c("No","NursingHome","Death"))
df.nh$event = ifelse(df.nh$nh + df.nh$death == 0, 0, 
                     ifelse(df.nh$nh_month<=df.nh$death_month, 1, 2))
df.nh$event = factor(df.nh$event, levels = c(0,1,2), labels = c("No","NursingHome","Death"))
save(df.nh, file=paste0(path.out, file="survival_analysis_df_nh_revised_1.RData"))

# Baseline covariates
# 1. Age
library(eeptools)
df.bd = demo %>% mutate(patient_num = as.character(PatientNum), 
                        birth_date=base::as.Date(Date_of_Birth, format="%Y-%m-%d")) %>%
  select(patient_num, birth_date)
df.nh = df.nh %>% left_join(df.bd, by="patient_num") %>% mutate(age = round(age_calc(birth_date, diag_date)/12))
save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_revised_baseline_1.RData"))

# 2. co-morbidity index
library(comorbidity)
library(stringr)
library(zoo)

df.nh$baseline_start_date_1yr = as.Date(as.yearmon(df.nh$diag_date)-1, frac=1) #one-year
df.nh$baseline_start_date_2yr = as.Date(as.yearmon(df.nh$diag_date)-2, frac=1) #two-year
data.icd.baseline.1yr = data.icd %>% 
  left_join(df.nh%>%select(patient_num, diag_date, baseline_start_date_1yr), by="patient_num") %>%
  filter((date >= baseline_start_date_1yr) & (date <= diag_date))
data.icd.baseline.2yr = data.icd %>% 
  left_join(df.nh%>%select(patient_num, diag_date, baseline_start_date_2yr), by="patient_num") %>%
  filter((date >= baseline_start_date_2yr) & (date <= diag_date))

# One-year
data.icd9.baseline.1yr = data.icd.baseline.1yr %>% filter(grepl("ICD9", code)) %>% mutate(code = gsub("ICD9:|\\.","",code))
comorbidity.icd9.1yr = comorbidity(data.icd9.baseline.1yr, "patient_num", "code", map = "elixhauser_icd9_quan", assign0 = FALSE)
data.icd10.baseline.1yr = data.icd.baseline.1yr %>% filter(grepl("ICD10", code)) %>% mutate(code = str_pad(gsub("ICD10:|\\.","",code), 7, "right"))
comorbidity.icd10.1yr = comorbidity(data.icd10.baseline.1yr, "patient_num", "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
comorbidity.all.1yr = comorbidity.icd9.1yr %>% full_join(comorbidity.icd10.1yr, by="patient_num", suffix = c("_1yr_9", "_1yr_10")) %>% replace(is.na(.), 0)
exli.nm = colnames(comorbidity.icd9.1yr)[-1]
for (nm in exli.nm) {
  comorbidity.all.1yr[,paste0(nm,"_1yr")] = ifelse(comorbidity.all.1yr[,paste0(nm,"_1yr_9")] + comorbidity.all.1yr[,paste0(nm,"_1yr_10")] >=1, 1, 0)
}
# Two-year
data.icd9.baseline.2yr = data.icd.baseline.2yr %>% filter(grepl("ICD9", code)) %>% mutate(code = gsub("ICD9:|\\.","",code))
comorbidity.icd9.2yr = comorbidity(data.icd9.baseline.2yr, "patient_num", "code", map = "elixhauser_icd9_quan", assign0 = FALSE)
data.icd10.baseline.2yr = data.icd.baseline.2yr %>% filter(grepl("ICD10", code)) %>% mutate(code = str_pad(gsub("ICD10:|\\.","",code), 7, "right"))
comorbidity.icd10.2yr = comorbidity(data.icd10.baseline.2yr, "patient_num", "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
comorbidity.all.2yr = comorbidity.icd9.2yr %>% full_join(comorbidity.icd10.2yr, by="patient_num", suffix = c("_2yr_9", "_2yr_10")) %>% replace(is.na(.), 0)
exli.nm = colnames(comorbidity.icd9.2yr)[-1]
for (nm in exli.nm) {
  comorbidity.all.2yr[,paste0(nm,"_2yr")] = ifelse(comorbidity.all.2yr[,paste0(nm,"_2yr_9")] + comorbidity.all.2yr[,paste0(nm,"_2yr_10")] >=1, 1, 0)
}
# Merge with df.nh
df.nh = df.nh %>% left_join(comorbidity.all.1yr %>% select(c(patient_num, paste0(exli.nm,"_1yr"))), by="patient_num") %>%
  left_join(comorbidity.all.2yr %>% select(c(patient_num, paste0(exli.nm,"_2yr"))), by="patient_num")
save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_revised_baseline_1.RData"))

# 3. Healthcare utilization
df.temp = combined_data %>% 
  distinct(PatientNum, Start_date, Local_Code) %>%				## unique features each date per patient
  group_by(PatientNum, Start_date) %>% 					## days of each feature per patient
  summarize(count = n()) %>% ungroup() 
df.temp <- df.temp %>% 
  mutate(patient_num = as.character(PatientNum)) %>%
  mutate(start_date = as.Date(Start_date, format = "%Y-%m-%d")) %>%
  inner_join(diag.date, by="patient_num") %>%
  mutate(num_month = mondf(diag_date, start_date)) %>%
  select(patient_num, start_date, diag_date, count, num_month)

df.codect.1yr = df.temp %>% filter((num_month>=-12) & (num_month <=0)) %>%
  group_by(patient_num) %>%
  summarise(codect_1yr = sum(count))
df.codect.2yr = df.temp %>% filter((num_month>=-24) & (num_month <=0)) %>%
  group_by(patient_num) %>%
  summarise(codect_2yr = sum(count))

df.nh = df.nh %>% 
  left_join(df.codect.1yr, by="patient_num") %>%
  left_join(df.codect.2yr, by="patient_num")
save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_revised_baseline_util.RData"))

# 4. First date of EHR encounter (for excluding those with insufficient pre-index follow-up duration later)
df.first = df.temp %>% group_by(patient_num) %>% summarise(first_month=min(num_month))
df.nh = df.nh %>% left_join(df.first, by="patient_num")
save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_revised_baseline_1.RData"))

#################################
# Creating survival analysis matrices
#################################

library(survival)
library(survminer)
library(cmprsk)

# Names of features
features.nm = c("white", "female", "age", "log(util+1)", "chf", "carit", "valv", "pcd", "pvd", "hypunc", "hypc", "para", "ond",
                "cpd", "diabunc", "diabc", "hypothy", "rf", "ld", "pud", "aids", "lymph",
                "metacanc", "solidtum", "rheumd", "coag", "obes", "wloss", "fed", "blane", 
                "dane", "alcohol", "drug", "psycho", "depre")

# Nursing home codified only
df.code = df.nh[(df.nh$event_code_month>0),]; dim(df.code)
df.code$white <- ifelse(is.na(df.code$white), 0, ifelse(df.code$white == 1, 1, 0))
df.code.1yr = df.code[df.code$first_month<=-12,]; dim(df.code.1yr)
df.code.2yr = df.code[df.code$first_month<=-24,]; dim(df.code.2yr)
df.code.1yr.clean <- df.code.1yr[complete.cases(df.code.1yr[, c("white", "female", "age", "codect_1yr", "chf_1yr", "carit_1yr", 
                                                                "valv_1yr", "pcd_1yr", "pvd_1yr", "hypunc_1yr", "hypc_1yr", 
                                                                "para_1yr", "ond_1yr", "cpd_1yr", "diabunc_1yr", "diabc_1yr", 
                                                                "hypothy_1yr", "rf_1yr", "ld_1yr", "pud_1yr", "aids_1yr",
                                                                "lymph_1yr", "metacanc_1yr", "solidtum_1yr", "rheumd_1yr", 
                                                                "coag_1yr", "obes_1yr", "wloss_1yr", "fed_1yr", "blane_1yr", 
                                                                "dane_1yr", "alcohol_1yr", "drug_1yr", "psycho_1yr", "depre_1yr")]), ]
df.code.1yr.mat <- model.matrix(object = ~ white + female + age + log(codect_1yr + 1) + 
                                  chf_1yr + carit_1yr + valv_1yr + pcd_1yr + pvd_1yr + hypunc_1yr + hypc_1yr + para_1yr + ond_1yr + 
                                  cpd_1yr + diabunc_1yr + diabc_1yr + hypothy_1yr + rf_1yr + ld_1yr + pud_1yr + aids_1yr + lymph_1yr +
                                  metacanc_1yr + solidtum_1yr + rheumd_1yr + coag_1yr + obes_1yr + wloss_1yr + fed_1yr + blane_1yr +
                                  dane_1yr + alcohol_1yr + drug_1yr + psycho_1yr + depre_1yr, data = df.code.1yr.clean)[,-1]
df.code.1yr.mat.s = cbind(df.code.1yr.mat[,c(1,2,3,4)], rowSums(df.code.1yr.mat[,-c(1,2,3,4)])); colnames(df.code.1yr.mat.s)[5] = "score"
df.code.2yr.clean <- df.code.2yr[complete.cases(df.code.2yr[, c("white", "female", "age", "codect_2yr", "chf_2yr", "carit_2yr", 
                                                                "valv_2yr", "pcd_2yr", "pvd_2yr", "hypunc_2yr", "hypc_2yr", 
                                                                "para_2yr", "ond_2yr", "cpd_2yr", "diabunc_2yr", "diabc_2yr", 
                                                                "hypothy_2yr", "rf_2yr", "ld_2yr", "pud_2yr", "aids_2yr",
                                                                "lymph_2yr", "metacanc_2yr", "solidtum_2yr", "rheumd_2yr", 
                                                                "coag_2yr", "obes_2yr", "wloss_2yr", "fed_2yr", "blane_2yr", 
                                                                "dane_2yr", "alcohol_2yr", "drug_2yr", "psycho_2yr", "depre_2yr")]), ]
df.code.2yr.mat = model.matrix(object = ~ white + female + age + log(codect_2yr+1) + 
                                 chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+ 
                                 cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                                 metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                                 dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data = df.code.2yr.clean)[,-1]
df.code.2yr.mat.s = cbind(df.code.2yr.mat[,c(1,2,3,4)], rowSums(df.code.2yr.mat[,-c(1,2,3,4)])); colnames(df.code.2yr.mat.s)[5] = "score"

# Nursing home codified or NLP
df.all = df.nh[(df.nh$event_month>0),]; dim(df.all)
df.all.1yr = df.all[df.all$first_month<=-12,]; dim(df.all.1yr)
df.all.2yr = df.all[df.all$first_month<=-24,]; dim(df.all.2yr)
df.all.1yr.clean <- df.all.1yr[complete.cases(df.all.1yr[, c("white", "female", "age", "codect_1yr", "chf_1yr", "carit_1yr", 
                                                             "valv_1yr", "pcd_1yr", "pvd_1yr", "hypunc_1yr", "hypc_1yr", 
                                                             "para_1yr", "ond_1yr", "cpd_1yr", "diabunc_1yr", "diabc_1yr", 
                                                             "hypothy_1yr", "rf_1yr", "ld_1yr", "pud_1yr", "aids_1yr",  
                                                             "lymph_1yr", "metacanc_1yr", "solidtum_1yr", "rheumd_1yr", 
                                                             "coag_1yr", "obes_1yr", "wloss_1yr", "fed_1yr", "blane_1yr", 
                                                             "dane_1yr", "alcohol_1yr", "drug_1yr", "psycho_1yr", "depre_1yr")]), ]
df.all.1yr.mat = model.matrix(object = ~ white + female + age + log(codect_1yr+1) + 
                                chf_1yr + carit_1yr + valv_1yr + pcd_1yr + pvd_1yr + hypunc_1yr + hypc_1yr + para_1yr + ond_1yr+ 
                                cpd_1yr + diabunc_1yr + diabc_1yr + hypothy_1yr + rf_1yr + ld_1yr + pud_1yr + aids_1yr + lymph_1yr +
                                metacanc_1yr + solidtum_1yr + rheumd_1yr + coag_1yr + obes_1yr + wloss_1yr + fed_1yr + blane_1yr +
                                dane_1yr + alcohol_1yr + drug_1yr + psycho_1yr + depre_1yr, data = df.all.1yr)[,-1]
df.all.1yr.mat.s = cbind(df.all.1yr.mat[,c(1,2,3,4)], rowSums(df.all.1yr.mat[,-c(1,2,3,4)])); colnames(df.all.1yr.mat.s)[5] = "score"
df.all.2yr.clean <- df.all.2yr[complete.cases(df.all.2yr[, c("white", "female", "age", "codect_2yr", "chf_2yr", "carit_2yr", 
                                                             "valv_2yr", "pcd_2yr", "pvd_2yr", "hypunc_2yr", "hypc_2yr", 
                                                             "para_2yr", "ond_2yr", "cpd_2yr", "diabunc_2yr", "diabc_2yr", 
                                                             "hypothy_2yr", "rf_2yr", "ld_2yr", "pud_2yr",  "aids_2yr",
                                                             "lymph_2yr", "metacanc_2yr", "solidtum_2yr", "rheumd_2yr", 
                                                             "coag_2yr", "obes_2yr", "wloss_2yr", "fed_2yr", "blane_2yr", 
                                                             "dane_2yr", "alcohol_2yr", "drug_2yr", "psycho_2yr", "depre_2yr")]), ]
df.all.2yr.mat = model.matrix(object = ~ white + female + age + log(codect_2yr+1) + 
                                chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+ 
                                cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                                metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                                dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data = df.all.2yr)[,-1]
df.all.2yr.mat.s = cbind(df.all.2yr.mat[,c(1,2,3,4)], rowSums(df.all.2yr.mat[,-c(1,2,3,4)])); colnames(df.all.2yr.mat.s)[5] = "score"
save(df.all.2yr.clean, file = paste0(path.out, "df_all_2yr_clean_mgb.RData"))

#################################
# Competing risk Cox PH model
# Uncomment outcome=NursingHome or Death
#################################

#outcome = "NursingHome"
outcome = "Death"
# Nursing home codified only, util codified only, 1yr
CRR.nh = crr(ftime = df.code.1yr.clean$event_code_month, fstatus = df.code.1yr.clean$event_code, cov1 = df.code.1yr.mat, failcode = outcome)
HR.code.code.1yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.code.code.1yr$nh_type = "codified only"
HR.code.code.1yr$base_type = "12 months"
HR.code.code.1yr$features = features.nm

# Nursing home codified only, util codified only, 2yr
CRR.nh = crr(ftime = df.code.2yr.clean$event_code_month, fstatus = df.code.2yr.clean$event_code, cov1 = df.code.2yr.mat, failcode = outcome)
HR.code.code.2yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.code.code.2yr$nh_type = "codified only"
HR.code.code.2yr$base_type = "24 months"
HR.code.code.2yr$features = features.nm

# Nursing home codified+NLP, util codified only, 1yr
CRR.nh = crr(ftime = df.all.1yr.clean$event_month, fstatus = df.all.1yr.clean$event, cov1 = df.all.1yr.mat, failcode = outcome)
HR.both.code.1yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.both.code.1yr$nh_type = "codified & NLP"
HR.both.code.1yr$base_type = "12 months"
HR.both.code.1yr$features = features.nm

# Nursing home codified+NLP, util codified only, 2yr
CRR.nh = crr(ftime = df.all.2yr.clean$event_month, fstatus = df.all.2yr.clean$event, cov1 = df.all.2yr.mat, failcode = outcome)
HR.both.code.2yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.both.code.2yr$nh_type = "codified & NLP"
HR.both.code.2yr$base_type = "24 months"
HR.both.code.2yr$features = features.nm

#################################
# Obtaining adjusted hazard ratios
#################################

pval_sig = Vectorize(function(p) {
  if ((p <= 0.05)& (p >0.01)) {return ("*")}
  else if ((p <= 0.01)& (p >0.001)) {return ("**")}
  else if ((p <= 0.001)& (p >0.0001)) {return("***")}
  else if (p <= 0.0001) {return("****")}
  else {return("")}
})

HR.all = data.frame(rbind(HR.code.code.1yr, HR.code.code.2yr, 
                          HR.both.code.1yr, HR.both.code.2yr))
colnames(HR.all) = c("HR", "HR_lower", "HR_upper", "pval", "nh_type", "base_type", "features")
HR.all = HR.all %>% mutate(pval.sig = pval_sig(pval), features = factor(features, levels=features.nm), 
                           text = paste0(round(HR, 3), " (", round(HR_lower, 3), ",", round(HR_upper, 3), ")"))
HR.all$features = factor(HR.all$features, levels = HR.all$features[1:35], labels = c("race/ethnicity: non-Hispanic white", "gender: women", "age at AD diagnosis", "healthcare utilization: log (util + 1)", "congestive heart failure", "cardiac arrhythmias", "valvular disease", "pulmonary circulation disorders", "peripheral vascular disorder", "hypertension, uncomplicated", "hypertension, complicated", "paralysis", "other neurological disorders", "chronic pulmonary disease", "diabetes, uncontrolled", "diabetes, controlled", "hypothyroidism", "renal failure", "liver disease", "peptic ulcer disease", "AIDS/HIV", "lymphoma", "metastatic cancer", "solid tumor, without metastasis", "rheumatoid arthritis", "coagulopathy", "obesity", "weight loss", "fluid and electrolyte disorders", "blood loss anemia", "deficiency anemia", "alcohol abuse", "drug abuse", "psychoses", "depression"))
HR.all$sig = as.factor(ifelse(HR.all$pval<=0.05, 1, 0))

plot.HR = ggplot(HR.all[(HR.all$nh_type == "codified & NLP") & 
                          (HR.all$base_type == "24 months"),], aes(x=forcats::fct_rev(features), y=HR, color=sig)) +
  scale_color_manual(values=c("black","red")) +
  guides(color="none") +
  geom_point(position=position_dodge(width = 0.5)) +
  geom_errorbar(width=.1, position=position_dodge(width = 0.5), aes(ymin=HR_lower, ymax=HR_upper)) +
  geom_hline(yintercept=1, linetype="dashed") + 
  geom_text(aes(y=3.2, label=pval.sig), color="red", size=6) +   
  geom_text(aes(y=2.7, label=text), nudge_x = 0, size=5) +      
  ylab("Adjusted Hazard Ratio") +
  xlab("Baseline Features") +
  coord_flip() +
  scale_y_continuous(limits=c(0,3.3), expand=c(0,0)) +          
  theme_bw() + 
  ggtitle("Death") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18)
  )
plot.HR

write.csv(HR.all, file=paste0(path.out, "adjusted_HR_death_revised.csv"))
ggsave(plot.HR, file=paste0(path.out, "HR_death_revised.png"), width=7.5 ,height=7.5, dpi=300)

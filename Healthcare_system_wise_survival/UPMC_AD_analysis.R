library(lubridate)
library(dplyr)
library(data.table)
library(tidyverse)
library(tidyr)
library(stringr)
library(ggplot2)
library(readxl)

##########################################################
# Reading in codified data and NLP data
## data.cod is the codified data file
## data is the NLP data file
## freq is the frequency summary file
## features is the list of NLP features
##########################################################

data.cod = read.csv("/UPMC_AD_2004_to_2023_Codified_rolled_up_data_2023-10-06.csv")
data = read.csv("/UPMC_AD_2011_to_2021_Codified+NLP_data_aggregated_2023-10-20_Lin.csv")
freq = read.csv("/UPMC_AD_2011_to_2021_Codified+NLP_frequency_summary_with_desc_2023-03-15.csv")
features = c(freq[freq$freq.code > 5, "feature_id"], "AB0000AD")

# Output
path.out = "/your_path/"

##########################################################
# Nursing home admission
## nh is the list of nursing home diagnosis codes in UPMC data
## nh.cui is the list of nursing home CUIs
##########################################################

## Read in nursing home admission codes
## Revised definition
nh_def_old_revised <- read_excel("/NH_codes_2023_new_def.xlsx")

# Get time to nursing home admission
nh = read.csv("/UPMC_AD_2004_to_2023_Codified_Diagnostic_Codes_Nursing_Home_Admission_2023-10-06.csv")
nh_revised_df <- nh %>%
  filter(concept_cd %in% nh_def_old_revised$concept_cd) 
nh.revised.date = nh_revised_df %>% 
  mutate(patient_num = as.character(patient_num)) %>%
  group_by(patient_num) %>%
  mutate(nh_date = as.Date(start_date, "%Y-%m-%d")) %>%
  arrange(nh_date) %>%
  filter(row_number()==1) %>%
  select(patient_num, nh_date)

## List of CUIs
nh.cuis = read.table("/nursing_home_cuis_from_dict.dsv", sep="|")$V2
nh_cui <- data %>%
  filter(feature_id %in% nh.cuis) %>%
  distinct(patient_num, start_date, feature_id) 
nh.cui.date = data %>% filter(feature_id %in% nh.cuis) %>%
  group_by(patient_num) %>%
  mutate(nh_cui_date = as.Date(start_date, "%Y-%m-%d")) %>%
  arrange(nh_cui_date) %>%
  filter(row_number()==1) %>%
  select(patient_num, nh_cui_date)  %>% 
  mutate(patient_num = as.character(patient_num))

save(nh.revised.date, nh.cui.date, file=paste0(path.out, "nursing_home_admission_date.RData"))

##########################################################
# UPMC imputed demographics
##########################################################

demo <- read.csv("/UPMC demographics imputed race ethnicity separate.csv")
demo_datefix <- fread("/ad_patient_demographics-datefix.dsv")

## Recoding gender
demo$GENDER_TITLE <- as.factor(demo$GENDER_TITLE)
demo$gender <- case_when(demo$GENDER_TITLE == "FEMALE" ~ 1,
                         demo$GENDER_TITLE == "MALE" ~ 0,
                         demo$GENDER_TITLE == "UNKNOWN" ~ NA,
                         demo$GENDER_TITLE == "" ~ NA,
                         TRUE ~ NA)
d.gender = demo[,c("PATIENT_STUDY_ID", "gender")]
names(d.gender) <- c("patient_num", "female")
## Recoding race and ethnicity
demo$race <- as.factor(demo$race)
demo$ethnicity <- as.factor(demo$ethnicity)
demo$white = ifelse((demo$race == "White") & (demo$ethnicity != "Hispanic"), 1, 0)
d.race = demo[,c("PATIENT_STUDY_ID", "white")]
names(d.race) <- c("patient_num", "white")

## Mortality date
demo_datefix$`TO_CHAR(DEATH_DATE,'MM-DD-YYYY')` <- as.Date(demo_datefix$`TO_CHAR(DEATH_DATE,'MM-DD-YYYY')`, "%m-%d-%Y")
mortality.date = demo_datefix %>% mutate(patient_num = as.character(demo_datefix$PATIENT_STUDY_ID), 
                                         death_date = as.Date(demo_datefix$`TO_CHAR(DEATH_DATE,'MM-DD-YYYY')`)) %>%
  select(patient_num, death_date)
save(mortality.date, file=paste0(path.out, "mortality_date.RData"))

## AD diagnosis date
diag.date = data.cod %>% filter(feature_id == "PheCode:290.11") %>% 
  group_by(patient_num) %>% 
  mutate(diag_date = as.Date(start_date, "%Y-%m-%d")) %>%
  arrange(diag_date) %>%
  filter(row_number()==1) %>%
  select(patient_num,diag_date) %>%
  mutate(patient_num = as.character(patient_num))
save(diag.date, file=paste0(path.out, "AD_diagnosis_date.RData"))

# Baseline co-morbidity - get raw icd
path.icd = "/path_to_icd/"
filenames = Sys.glob(paste0(path.icd, "*.csv"))
data.icd = NULL
for (f in filenames) {
  print(f)
  data.icd = rbind(data.icd, read.csv(f, header=F)[,1:5] %>% filter(grepl("ICD", V2)) %>% select(V1, V2, V3))
}
colnames(data.icd) = c("patient_num", "code", "date")
data.icd$patient_num = as.character(data.icd$patient_num)
save(data.icd, file=paste0(path.out, "ICD.RData"))

# AD-related medications
med.list = c("RXNORM:135447", "RXNORM:183379", "RXNORM:4637", "RXNORM:6719", 
             "RXNORM:2557217", "RXNORM:2626143", "RXNORM:1430990")
med.date = data.cod %>% filter(feature_id %in% med.list) %>%
  mutate(patient_num = as.character(patient_num)) %>%
  group_by(patient_num) %>%
  mutate(med_date = as.Date(start_date, "%Y-%m-%d")) %>%
  arrange(med_date) %>%
  filter(row_number()==1) %>%
  select(patient_num, med_date)
save(med.date, file=paste0(path.out, "medication_date.RData"))

#################################
## Merge data
#################################

# Get AD diagnosis status as predicted by KOMAP
pred <- read.csv("/komap_pred_imputed.csv")
pred <- pred[,c("patient_num", "pred.outcome.1")]
d.race$patient_num = as.character(d.race$patient_num)
d.gender$patient_num = as.character(d.gender$patient_num)
diag.date$patient_num = as.character(diag.date$patient_num)
nh.revised.date$patient_num = as.character(nh.revised.date$patient_num)
nh.cui.date$patient_num = as.character(nh.cui.date$patient_num)
med.date$patient_num = as.character(med.date$patient_num)

df = pred %>% mutate(patient_num = as.character(patient_num)) %>%
  filter(pred.outcome.1==1) %>% left_join(d.race, by="patient_num")%>% 
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
save(df, file=paste0(path.out, "survival_analysis_df_revised_1.RData"))

#################################
# Data aggregation for survival analysis
#################################

monnb = function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01")); lt$year*12 + lt$mon } 
mondf = function(d1, d2) { monnb(d2) - monnb(d1) }

# Time to death
library(survival)
df = df[!is.na(df$diag_date),]
df$death = ifelse(is.na(df$death_month), 0, 1)
df$death_month = ifelse(is.na(df$death_month), mondf(df$diag_date, "2022-12-31"), df$death_month)

# Time to nursing home admission
df$nh_code = ifelse(is.na(df$nh_code_month), 0, 1)
df$nh_cui = ifelse(is.na(df$nh_cui_month), 0, 1)
df$nh = ifelse(df$nh_code+df$nh_cui>=1, 1, 0)

df$nh_code_month = ifelse(is.na(df$nh_code_month), mondf(df$diag_date, "2022-12-31"), df$nh_code_month)
df$nh_cui_month = ifelse(is.na(df$nh_cui_month), mondf(df$diag_date, "2022-12-31"), df$nh_cui_month)
df$nh_month = pmin(df$nh_code_month, df$nh_cui_month)

df.nh = df
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
df.bd = demo_datefix %>% mutate(patient_num = as.character(PATIENT_STUDY_ID), 
                        birth_date=base::as.Date(`TO_CHAR(BIRTH_DATE,'MM-DD-YYYY')`, format="%m-%d-%Y")) %>%
  select(patient_num, birth_date)
df.nh = df.nh %>% left_join(df.bd, by="patient_num") %>% mutate(age = round(age_calc(birth_date, diag_date)/12))
save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_baseline_1.RData"))

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
data.icd9.baseline.1yr = data.icd.baseline.1yr %>% filter(grepl("ICD9CM", code)) %>% mutate(code = gsub("ICD9CM:|\\.","",code))
comorbidity.icd9.1yr = comorbidity(data.icd9.baseline.1yr, "patient_num", "code", map = "elixhauser_icd9_quan", assign0 = FALSE)
data.icd10.baseline.1yr = data.icd.baseline.1yr %>% filter(grepl("ICD10CM", code)) %>% mutate(code = str_pad(gsub("ICD10CM:|\\.","",code), 7, "right"))
comorbidity.icd10.1yr = comorbidity(data.icd10.baseline.1yr, "patient_num", "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
comorbidity.all.1yr = comorbidity.icd9.1yr %>% full_join(comorbidity.icd10.1yr, by="patient_num", suffix = c("_1yr_9", "_1yr_10")) %>% replace(is.na(.), 0)
exli.nm = colnames(comorbidity.icd9.1yr)[-1]
for (nm in exli.nm) {
  comorbidity.all.1yr[,paste0(nm,"_1yr")] = ifelse(comorbidity.all.1yr[,paste0(nm,"_1yr_9")] + comorbidity.all.1yr[,paste0(nm,"_1yr_10")] >=1, 1, 0)
}
# Two-year
data.icd9.baseline.2yr = data.icd.baseline.2yr %>% filter(grepl("ICD9CM", code)) %>% mutate(code = gsub("ICD9CM:|\\.","",code))
comorbidity.icd9.2yr = comorbidity(data.icd9.baseline.2yr, "patient_num", "code", map = "elixhauser_icd9_quan", assign0 = FALSE)
data.icd10.baseline.2yr = data.icd.baseline.2yr %>% filter(grepl("ICD10CM", code)) %>% mutate(code = str_pad(gsub("ICD10CM:|\\.","",code), 7, "right"))
comorbidity.icd10.2yr = comorbidity(data.icd10.baseline.2yr, "patient_num", "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
comorbidity.all.2yr = comorbidity.icd9.2yr %>% full_join(comorbidity.icd10.2yr, by="patient_num", suffix = c("_2yr_9", "_2yr_10")) %>% replace(is.na(.), 0)
exli.nm = colnames(comorbidity.icd9.2yr)[-1]
for (nm in exli.nm) {
  comorbidity.all.2yr[,paste0(nm,"_2yr")] = ifelse(comorbidity.all.2yr[,paste0(nm,"_2yr_9")] + comorbidity.all.2yr[,paste0(nm,"_2yr_10")] >=1, 1, 0)
}
# Merge with df.nh
df.nh = df.nh %>% left_join(comorbidity.all.1yr %>% select(c(patient_num, paste0(exli.nm,"_1yr"))), by="patient_num") %>%
  left_join(comorbidity.all.2yr %>% select(c(patient_num, paste0(exli.nm,"_2yr"))), by="patient_num")

save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_baseline_1.RData"))

# 3. Healthcare utilization
df.temp = data.cod %>% filter(feature_id %in% features) %>%
  distinct(patient_num, start_date, feature_id) %>%				## unique features each date per patient
  group_by(patient_num, start_date) %>% 					## days of each feature per patient
  summarize(count = n()) %>% ungroup() %>% 
  mutate(patient_num = as.character(patient_num)) %>%
  inner_join(df.nh[c("patient_num", "diag_date")], by="patient_num") %>%
  mutate(num_month = mondf(diag_date, start_date))

df.codect.1yr = df.temp %>% filter((num_month>=-12) & (num_month <=0)) %>%
  group_by(patient_num) %>%
  summarise(codect_1yr = sum(count))
df.codect.2yr = df.temp %>% filter((num_month>=-24) & (num_month <=0)) %>%
  group_by(patient_num) %>%
  summarise(codect_2yr = sum(count))

df.nh = df.nh %>% 
  left_join(df.codect.1yr, by="patient_num") %>%
  left_join(df.codect.2yr, by="patient_num")

save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_baseline_1.RData"))

# 4. First date of EHR encounter (for excluding those with insufficient pre-index follow-up duration later)
df.first = df.temp %>% group_by(patient_num) %>% summarise(first_month=min(num_month))
df.nh = df.nh %>% left_join(df.first, by="patient_num")
save(df.nh, file=paste0(path.out, "survival_analysis_df_nh_baseline_1.RData"))

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
df.code.1yr = df.code[df.code$first_month<=-12,]; dim(df.code.1yr)
df.code.2yr = df.code[df.code$first_month<=-24,]; dim(df.code.2yr)
df.code.1yr.mat = model.matrix(object = ~ white + female + age + log(codect_1yr+1) + 
                                 chf_1yr + carit_1yr + valv_1yr + pcd_1yr + pvd_1yr + hypunc_1yr + hypc_1yr + para_1yr + ond_1yr+ 
                                 cpd_1yr + diabunc_1yr + diabc_1yr + hypothy_1yr + rf_1yr + ld_1yr + pud_1yr + aids_1yr + lymph_1yr +
                                 metacanc_1yr + solidtum_1yr + rheumd_1yr + coag_1yr + obes_1yr + wloss_1yr + fed_1yr + blane_1yr +
                                 dane_1yr + alcohol_1yr + drug_1yr + psycho_1yr + depre_1yr, data = df.code.1yr)[,-1]
df.code.1yr.mat.s = cbind(df.code.1yr.mat[,c(1,2,3,4)], rowSums(df.code.1yr.mat[,-c(1,2,3,4)])); colnames(df.code.1yr.mat.s)[5] = "score"
df.code.2yr.mat = model.matrix(object = ~ white + female + age + log(codect_2yr+1) + 
                                 chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+ 
                                 cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                                 metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                                 dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data = df.code.2yr)[,-1]
df.code.2yr.mat.s = cbind(df.code.2yr.mat[,c(1,2,3,4)], rowSums(df.code.2yr.mat[,-c(1,2,3,4)])); colnames(df.code.2yr.mat.s)[5] = "score"

# Nursing home codified or NLP
df.all = df.nh[(df.nh$event_month>0),]; dim(df.all)
df.all.1yr = df.all[df.all$first_month<=-12,]; dim(df.all.1yr)
df.all.2yr = df.all[df.all$first_month<=-24,]; dim(df.all.2yr)
df.all.1yr.mat = model.matrix(object = ~ white + female + age + log(codect_1yr+1) + 
                                chf_1yr + carit_1yr + valv_1yr + pcd_1yr + pvd_1yr + hypunc_1yr + hypc_1yr + para_1yr + ond_1yr+ 
                                cpd_1yr + diabunc_1yr + diabc_1yr + hypothy_1yr + rf_1yr + ld_1yr + pud_1yr + aids_1yr + lymph_1yr +
                                metacanc_1yr + solidtum_1yr + rheumd_1yr + coag_1yr + obes_1yr + wloss_1yr + fed_1yr + blane_1yr +
                                dane_1yr + alcohol_1yr + drug_1yr + psycho_1yr + depre_1yr, data = df.all.1yr)[,-1]
df.all.1yr.mat.s = cbind(df.all.1yr.mat[,c(1,2,3,4)], rowSums(df.all.1yr.mat[,-c(1,2,3,4)])); colnames(df.all.1yr.mat.s)[5] = "score"
df.all.2yr.mat = model.matrix(object = ~ white + female + age + log(codect_2yr+1) + 
                                chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+ 
                                cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                                metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                                dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data = df.all.2yr)[,-1]
df.all.2yr.mat.s = cbind(df.all.2yr.mat[,c(1,2,3,4)], rowSums(df.all.2yr.mat[,-c(1,2,3,4)])); colnames(df.all.2yr.mat.s)[5] = "score"
save(df.all.2yr, file = paste0(path.out, "df_all_2yr_clean_upmc.RData"))

#################################
# Competing risk Cox PH model
# Uncomment outcome=NursingHome or Death
#################################

outcome = "Death"
#outcome = "NursingHome"
# Nursing home codified only, util codified only, 1yr
CRR.nh = crr(ftime = df.code.1yr$event_code_month, fstatus = df.code.1yr$event_code, cov1 = df.code.1yr.mat, failcode = outcome)
HR.code.code.1yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.code.code.1yr$nh_type = "codified only"
HR.code.code.1yr$base_type = "12 months"
HR.code.code.1yr$features = features.nm

# Nursing home codified only, util codified only, 2yr
CRR.nh = crr(ftime = df.code.2yr$event_code_month, fstatus = df.code.2yr$event_code, cov1 = df.code.2yr.mat, failcode = outcome)
HR.code.code.2yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.code.code.2yr$nh_type = "codified only"
HR.code.code.2yr$base_type = "24 months"
HR.code.code.2yr$features = features.nm

# Nursing home codified+NLP, util codified only, 1yr
CRR.nh = crr(ftime = df.all.1yr$event_month, fstatus = df.all.1yr$event, cov1 = df.all.1yr.mat, failcode = outcome)
HR.both.code.1yr = data.frame(cbind(summary(CRR.nh)$conf.int[,c("exp(coef)", "2.5%", "97.5%")], summary(CRR.nh)$coef[,"p-value"]))
HR.both.code.1yr$nh_type = "codified & NLP"
HR.both.code.1yr$base_type = "12 months"
HR.both.code.1yr$features = features.nm

# Nursing home codified+NLP, util codified only, 2yr
CRR.nh = crr(ftime = df.all.2yr$event_month, fstatus = df.all.2yr$event, cov1 = df.all.2yr.mat, failcode = outcome)
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

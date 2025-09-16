library(gtsummary)
library(dplyr)
library(data.table)
library(tidyr)
library(broom)
library(survival)

################################################################################
# SURVIVAL ANALYSIS SITEWISE DATA
################################################################################

# Getting UPMC survival analysis data
load("/df_all_2yr_clean_upmc.RData")
upmc_impute_demo <- df.all.2yr
upmc_impute_demo$site <- "UPMC"

# Getting MGB survival analysis data
load("/df_all_2yr_clean_mgb.RData")
mgb_impute_demo <- df.all.2yr.clean
mgb_impute_demo$site <- "MGB"
colnames(mgb_impute_demo)[2] <- "pred.outcome.1"

################################################################################
# IMPUTED 
################################################################################

# UPMC imputed demographics
demo_upmc <- read.csv("/UPMC demographics imputed race ethnicity separate.csv")
demo_upmc$patient_num <- as.character(demo_upmc$PATIENT_STUDY_ID)
demo_upmc_race_eth <- demo_upmc[,c("patient_num", "race", "ethnicity")]
upmc_impute_demo <- upmc_impute_demo %>% left_join(demo_upmc_race_eth, by = "patient_num")

# MGB imputed demographics
demo_mgb <- read.csv("/MGB demographics imputed race ethnicity separate.csv")
demo_mgb$patient_num <- as.character(demo_mgb$PatientNum)
demo_mgb_race_eth <- demo_mgb[,c("patient_num", "race", "ethnicity")]
mgb_impute_demo <- mgb_impute_demo %>% left_join(demo_mgb_race_eth, by = "patient_num")

# Combining demographics data
imputed_com <- rbind(upmc_impute_demo, mgb_impute_demo)
imputed_com$race_eth_final <- case_when(imputed_com$race %in% c("White") & imputed_com$ethnicity == "Non Hispanic" ~ "Non Hispanic White", 
                                        imputed_com$race %in%  c("American Indian or Alaska Native", "Asian", "Black", "Native Hawaiian or Other Pacific Islander", "Other", "White") & imputed_com$ethnicity == "Hispanic" ~ "Hispanic or Latino", 
                                        imputed_com$race %in% c("American Indian or Alaska Native", "Asian", "Black",  "Native Hawaiian or Other Pacific Islander", "Other") & imputed_com$ethnicity == "Non Hispanic" ~ imputed_com$race, 
                                        TRUE ~ NA)

# Comorbidity calculation
imputed_com <- imputed_com %>%
  rowwise() %>%
  mutate(total_score = sum(c_across(chf_1yr:depre_2yr), na.rm = TRUE)) %>%
  ungroup()

weights <- c(chf = 7, carit = 5, valv = -1, pcd = 4, pvd = 2, hypunc = 0, hypc = 0, para = 7, ond = 6, cpd = 3, diabunc = 0, diabc = 0, hypothy = 0,
    rf = 5, ld = 11, pud = 0, aids = 0, lymph = 9, metacanc = 12, solidtum = 4, rheumd = 0, coag = 3, obes = -4, wloss = 6, fed = 5, blane = -2, dane = -2,
    alcohol = 0, drug = -7, psycho = 0, depre = -3)
imputed_comorbidity_df <- imputed_com[,c("patient_num", "chf_2yr", "carit_2yr", "valv_2yr", "pcd_2yr", "pvd_2yr",
                                         "hypunc_2yr", "hypc_2yr", "para_2yr", "ond_2yr", "cpd_2yr",
                                         "diabunc_2yr", "diabc_2yr", "hypothy_2yr", "rf_2yr", "ld_2yr",
                                         "pud_2yr", "aids_2yr", "lymph_2yr", "metacanc_2yr", "solidtum_2yr",
                                         "rheumd_2yr", "coag_2yr", "obes_2yr", "wloss_2yr", "fed_2yr",
                                         "blane_2yr", "dane_2yr", "alcohol_2yr", "drug_2yr", "psycho_2yr",
                                         "depre_2yr")]
names(imputed_comorbidity_df) <- sub("_2yr$", "", names(imputed_comorbidity_df))
imputed_comorbidity_df <- imputed_comorbidity_df %>%
  rowwise() %>%
  mutate(weighted_score = sum(c_across(all_of(names(weights))) * weights[names(weights)]))
imputed_com <- imputed_com %>%
  left_join(imputed_comorbidity_df[,c("patient_num", "weighted_score")], by = "patient_num")
imputed_com$elixhauser <- case_when(imputed_com$weighted_score < 0 ~ "<0",
                                    imputed_com$weighted_score == 0 ~ "0",
                                    imputed_com$weighted_score > 0 & imputed_com$weighted_score <= 5 ~ "1-5",
                                    imputed_com$weighted_score > 5 & imputed_com$weighted_score <= 13 ~ "6-13",
                                    imputed_com$weighted_score > 13 ~ "14+",
                                    TRUE ~ NA)

# Outcome
imputed_com$survival <- ifelse(imputed_com$death == 0, "Alive", "Death")
imputed_com$nursing_home <- ifelse(imputed_com$nh == 1, "Admitted", "Not Admitted")
# AD related medication
imputed_com$med <- ifelse(is.na(imputed_com$med_month), "Not prescribed", "Prescribed")
imputed_com$med_diff <- ifelse(imputed_com$med_month >= 0, "On or after AD diagnosis", "Before AD diagnosis")
save(imputed_com, file = "/imputed_com.RData")

################################################################################
# NOT IMPUTED 
################################################################################

# UPMC not imputed demographics
demo_datefix <- fread("/ad_patient_demographics-datefix.dsv")
## Recoding gender
demo_datefix <- demo_datefix %>% mutate(patient_num = as.character(PATIENT_STUDY_ID)) %>% mutate(patient_num = trimws(tolower(patient_num)))
demo_datefix$GENDER_TITLE <- as.factor(demo_datefix$GENDER_TITLE)
demo_datefix$gender <- case_when(demo_datefix$GENDER_TITLE == "FEMALE" ~ 1,
                                 demo_datefix$GENDER_TITLE == "MALE" ~ 0,
                                 demo_datefix$GENDER_TITLE == "UNKNOWN" ~ NA,
                                 demo_datefix$GENDER_TITLE == "" ~ NA,
                                 TRUE ~ NA)
## Recoding race and ethnicity
demo_datefix$RACE_TITLE <- as.factor(demo_datefix$RACE_TITLE)
demo_datefix$race <- factor(demo_datefix$RACE_TITLE, levels =  c("", "ALASKA NATIVE", "AMERICAN INDIAN/ALASKA NATIVE", "ASIAN INDIAN", "BLACK/AFRICAN AMERICAN", "CHINESE", "FILIPINO", "GUAMANIAN", "JAPANESE", "KOREAN", "NATIVE HAWAIIAN", "NOT SPECIFIED", "OTHER", "OTHER ASIAN", "OTHER PACIFIC ISLANDER", "SAMOAN", "UNREPORTED,CHOSE NOT TO DISCLOSE RACE", "VIETNAMESE", "WHITE"), labels = c("Unknown", "American Indian or Alaska Native", "American Indian or Alaska Native", "Asian", "Black", "Asian", "Asian", "Native Hawaiian or Other Pacific Islander","Asian", "Asian", "Native Hawaiian or Other Pacific Islander", "Unknown", "Other", "Asian", "Native Hawaiian or Other Pacific Islander", "Native Hawaiian or Other Pacific Islander", "Unknown", "Asian", "White"))
demo_datefix$ETHNIC_TITLE <- as.factor(demo_datefix$ETHNIC_TITLE)
demo_datefix$ethnicity <- factor(demo_datefix$ETHNIC_TITLE, levels = c("", "ANOTHER HISPANIC, LATINO/A, OR SPANISH ORIGIN", "MEXICAN, MEXICAN AMERICAN, CHICANO/A", "NON-HISPANIC OR LATINO/A", "NOT SPECIFIED", "PUERTO RICAN", "UNREPORTED/CHOSE NOT TO DISCLOSE"), labels = c("Unknown", "Hispanic", "Hispanic", "Non Hispanic", "Unknown", "Hispanic", "Unknown"))
demo_datefix <- demo_datefix %>% filter(race != "Unknown") %>% filter(ethnicity != "Unknown")
demo_upmc_not_imputed <- demo_datefix[,c("patient_num", "gender", "race", "ethnicity")] 
upmc_not_impute_demo <- merge(upmc_impute_demo, demo_upmc_not_imputed, by = "patient_num")

# MGB not imputed demographics
demo_mgb_noimpute <- read.csv("/MGB demographics deidentified.csv")
demo_mgb_noimpute <- demo_mgb_noimpute %>% mutate(patient_num = as.character(PatientNum)) %>% mutate(patient_num = trimws(tolower(patient_num)))
## Recoding gender
demo_mgb_noimpute$gender <- case_when(demo_mgb_noimpute$Gender_Legal_Sex == "Female" ~ 1,
                                      demo_mgb_noimpute$Gender_Legal_Sex == "Male" ~ 0,
                                      demo_mgb_noimpute$Gender_Legal_Sex == "Unknown-U" ~ NaN,
                                      TRUE ~ NA_real_)
## Recoding race and ethnicity
demo_mgb_noimpute$race <- case_when(demo_mgb_noimpute$Race_Group == "American Indian or Alaska Native" ~ 1,
                                    demo_mgb_noimpute$Race_Group == "Asian" ~ 2,
                                    demo_mgb_noimpute$Race_Group == "Black" ~ 3,
                                    demo_mgb_noimpute$Race_Group == "Native Hawaiian or Other Pacific Islander" ~ 4,
                                    demo_mgb_noimpute$Race_Group == "White" ~ 5,
                                    demo_mgb_noimpute$Race_Group == "Other" ~ 6,
                                    demo_mgb_noimpute$Race_Group == "Two or More" ~ 6,
                                    demo_mgb_noimpute$Race_Group == "Unknown/Missing" ~ 7,
                                    demo_mgb_noimpute$Race_Group == "Declined" ~ 7,
                                    TRUE ~ NA_real_)
demo_mgb_noimpute$race <- factor(demo_mgb_noimpute$race, levels = c(1:7), labels = c("American Indian or Alaska Native", "Asian",
                                                                                     "Black", "Native Hawaiian or Other Pacific Islander",
                                                                                     "White", "Other", "Unknown"))
demo_mgb_noimpute$ethnicity <- case_when(demo_mgb_noimpute$Ethnic_Group == "HISPANIC" ~ 1,
                                         demo_mgb_noimpute$Ethnic_Group == "Non Hispanic" ~ 2,
                                         demo_mgb_noimpute$Ethnic_Group == "Unknown/Missing" ~ 3,
                                         demo_mgb_noimpute$Ethnic_Group == "DECLINED" ~ 3,
                                         TRUE ~ NaN)
demo_mgb_noimpute$ethnicity <- factor(demo_mgb_noimpute$ethnicity, levels = c(1:3), labels = c("Hispanic", "Non Hispanic", "Unknown"))
demo_mgb_noimpute <- demo_mgb_noimpute %>% filter(race != "Unknown") %>% filter(ethnicity != "Unknown")
demo_mgb_not_imputed <- demo_mgb_noimpute[,c("patient_num", "gender", "race", "ethnicity")] 
mgb_not_impute_demo <- merge(mgb_impute_demo, demo_mgb_not_imputed, by = "patient_num")

# Combining demographics data
not_imputed_com <- rbind(upmc_not_impute_demo, mgb_not_impute_demo)
not_imputed_com$race_eth_final <- case_when(not_imputed_com$race %in% c("White") & not_imputed_com$ethnicity == "Non Hispanic" ~ "Non Hispanic White", 
                                        not_imputed_com$race %in%  c("American Indian or Alaska Native", "Asian", "Black", "Native Hawaiian or Other Pacific Islander", "Other", "White") & not_imputed_com$ethnicity == "Hispanic" ~ "Hispanic or Latino", 
                                        not_imputed_com$race %in% c("American Indian or Alaska Native", "Asian", "Black",  "Native Hawaiian or Other Pacific Islander", "Other") & not_imputed_com$ethnicity == "Non Hispanic" ~ not_imputed_com$race, 
                                        TRUE ~ NA)

# Comorbidity calculation
not_imputed_com <- not_imputed_com %>%
  rowwise() %>%
  mutate(total_score = sum(c_across(chf_1yr:depre_2yr), na.rm = TRUE)) %>%
  ungroup()
not_imputed_comorbidity_df <- not_imputed_com[,c("patient_num", "chf_2yr", "carit_2yr", "valv_2yr", "pcd_2yr", "pvd_2yr",
                                         "hypunc_2yr", "hypc_2yr", "para_2yr", "ond_2yr", "cpd_2yr",
                                         "diabunc_2yr", "diabc_2yr", "hypothy_2yr", "rf_2yr", "ld_2yr",
                                         "pud_2yr", "aids_2yr", "lymph_2yr", "metacanc_2yr", "solidtum_2yr",
                                         "rheumd_2yr", "coag_2yr", "obes_2yr", "wloss_2yr", "fed_2yr",
                                         "blane_2yr", "dane_2yr", "alcohol_2yr", "drug_2yr", "psycho_2yr",
                                         "depre_2yr")]
names(not_imputed_comorbidity_df) <- sub("_2yr$", "", names(not_imputed_comorbidity_df))
not_imputed_comorbidity_df <- not_imputed_comorbidity_df %>%
  rowwise() %>%
  mutate(weighted_score = sum(c_across(all_of(names(weights))) * weights[names(weights)]))

not_imputed_com <- not_imputed_com %>%
  left_join(not_imputed_comorbidity_df[,c("patient_num", "weighted_score")], by = "patient_num")
not_imputed_com$elixhauser <- case_when(not_imputed_com$weighted_score < 0 ~ "<0",
                                    not_imputed_com$weighted_score == 0 ~ "0",
                                    not_imputed_com$weighted_score > 0 & not_imputed_com$weighted_score <= 5 ~ "1-5",
                                    not_imputed_com$weighted_score > 5 & not_imputed_com$weighted_score <= 13 ~ "6-13",
                                    not_imputed_com$weighted_score > 13 ~ "14+",
                                    TRUE ~ NA)

# Outcome
not_imputed_com$survival <- ifelse(not_imputed_com$death == 0, "Alive", "Death")
not_imputed_com$nursing_home <- ifelse(not_imputed_com$nh == 1, "Admitted", "Not Admitted")
# AD related medication
not_imputed_com$med <- ifelse(is.na(not_imputed_com$med_month), "Not prescribed", "Prescribed")
not_imputed_com$med_diff <- ifelse(not_imputed_com$med_month >= 0, "On or after AD diagnosis", "Before AD diagnosis")
save(not_imputed_com, file = "/not_imputed_com.RData")

################################################################################
# TABLES
################################################################################

# Run the code below with imputed_com or not_imputed_com data frames

# Overall
table_1 <- select(not_imputed_com, c(age, female, race_eth_final, codect_2yr, total_score, hypunc_2yr, carit_2yr, survival, nursing_home, med, med_diff)) %>%
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
              digits = all_continuous() ~ 2,
              missing_text = "Missing",
              label = list(age ~ "Age at AD diagnosis", female ~ "Gender", race_eth_final ~ "Race and ethnicity", 
                           codect_2yr ~ "Baseline Healthcare Utilization", total_score ~ "Baseline Elixhauser comorbidity index", 
                           hypunc_2yr ~ "Hypertension, uncontrolled", carit_2yr ~ "Cardiac arrhythmia", 
                           survival ~ "Survival", nursing_home ~ "Nursing Home Admission", med ~ "AD-related medication", 
                           med_diff ~ "AD-related medication prescription")) %>%
  bold_labels() %>%
  modify_header(stat_0 ~ "**Overall**")

# Time to nursing home and death
time_nh_death <- not_imputed_com %>%
  group_by(event) %>%
  summarise(
    median_event_month = median(event_month, na.rm = TRUE),
    Q1_event_month = quantile(event_month, 0.25, na.rm = TRUE),
    Q3_event_month = quantile(event_month, 0.75, na.rm = TRUE)
  )
print(time_nh_death)

# By site
table_2 <- select(not_imputed_com, c(age, female, race_eth_final, codect_2yr, total_score, hypunc_2yr, carit_2yr, survival, nursing_home, med, med_diff, site)) %>%
  tbl_summary(
    by = site,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing_text = "Missing",
    label = list(age ~ "Age at AD diagnosis", female ~ "Gender", race_eth_final ~ "Race and ethnicity", 
                 codect_2yr ~ "Baseline Healthcare Utilization", total_score ~ "Baseline Elixhauser comorbidity index", 
                 hypunc_2yr ~ "Hypertension, uncontrolled", carit_2yr ~ "Cardiac arrhythmia", 
                 survival ~ "Survival", nursing_home ~ "Nursing Home Admission", med ~ "AD-related medication", 
                 med_diff ~ "AD-related medication prescription")) %>%
  add_p() %>%
  bold_labels() 

time_nh_death_site <- not_imputed_com %>%
  group_by(event, site) %>%
  summarise(
    median_event_month = median(event_month, na.rm = TRUE),
    Q1_event_month = quantile(event_month, 0.25, na.rm = TRUE),
    Q3_event_month = quantile(event_month, 0.75, na.rm = TRUE)
  )
print(time_nh_death_site)

p_values_site <- not_imputed_com %>%
  filter(event != "No") %>%            # remove "No"
  group_by(event) %>%                  # group by event
  summarise(
    p.value = if(n_distinct(site) == 2) {
      wilcox.test(event_month ~ site)$p.value
    } else {
      NA_real_
    }
  )
print(p_values_site)

# By race / ethnicity
table_3 <- select(not_imputed_com, c(age, female, white, race_eth_final, codect_2yr, total_score, hypunc_2yr, carit_2yr, survival, nursing_home, med, med_diff)) %>%
  tbl_summary(
    by = white,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing_text = "Missing",
    label = list(age ~ "Age at AD diagnosis", female ~ "Gender", race_eth_final ~ "Race and ethnicity", 
                 codect_2yr ~ "Baseline Healthcare Utilization", total_score ~ "Baseline Elixhauser comorbidity index", 
                 hypunc_2yr ~ "Hypertension, uncontrolled", carit_2yr ~ "Cardiac arrhythmia", 
                 survival ~ "Survival", nursing_home ~ "Nursing Home Admission", med ~ "AD-related medication", 
                 med_diff ~ "AD-related medication prescription")) %>%
  add_p() %>%
  bold_labels()

time_nh_death_race_eth <- not_imputed_com %>%
  group_by(event, white) %>%
  summarise(
    median_event_month = median(event_month, na.rm = TRUE),
    Q1_event_month = quantile(event_month, 0.25, na.rm = TRUE),
    Q3_event_month = quantile(event_month, 0.75, na.rm = TRUE)
  )
print(time_nh_death_race_eth)

p_values_race_eth <- not_imputed_com %>%
  filter(event != "No") %>%            # remove "No"
  group_by(event) %>%                  # group by event
  summarise(
    p.value = if(n_distinct(site) == 2) {
      wilcox.test(event_month ~ white)$p.value
    } else {
      NA_real_
    }
  )
print(p_values_race_eth)

# By gender
table_4 <- select(not_imputed_com, c(age, female, race_eth_final, codect_2yr, total_score, hypunc_2yr, carit_2yr, survival, nursing_home, med, med_diff)) %>%
  tbl_summary(
    by = female,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 2,
    missing_text = "Missing",
    label = list(age ~ "Age at AD diagnosis", race_eth_final ~ "Race and ethnicity", 
                 codect_2yr ~ "Baseline Healthcare Utilization", total_score ~ "Baseline Elixhauser comorbidity index", 
                 hypunc_2yr ~ "Hypertension, uncontrolled", carit_2yr ~ "Cardiac arrhythmia", 
                 survival ~ "Survival", nursing_home ~ "Nursing Home Admission", med ~ "AD-related medication", 
                 med_diff ~ "AD-related medication prescription")) %>%
  add_p() %>%
  bold_labels()

time_nh_death_gender <- not_imputed_com %>%
  group_by(event, female) %>%
  summarise(
    median_event_month = median(event_month, na.rm = TRUE),
    Q1_event_month = quantile(event_month, 0.25, na.rm = TRUE),
    Q3_event_month = quantile(event_month, 0.75, na.rm = TRUE)
  )
print(time_nh_death_gender)

p_values_gender <- not_imputed_com %>%
  filter(event != "No") %>%            # remove "No"
  group_by(event) %>%                  # group by event
  summarise(
    p.value = if(n_distinct(site) == 2) {
      wilcox.test(event_month ~ female)$p.value
    } else {
      NA_real_
    }
  )
print(p_values_gender)

# Merge all tables
combined_table <- tbl_merge(
  tbls = list(table_1, table_2, table_3, table_4),
  tab_spanner = c("**Overall**", "**By Site**", "**By Race/Ethnicity**", "**By Gender**")
) %>%
  modify_caption("**Table 1. Baseline Characteristics**")
combined_table_csv <- combined_table %>% as_tibble()
write.csv(combined_table_csv, "/baseline_characteristics_table_not_imputed.csv")

# Follow up time
# To get median time of followup from reverse KM estimator
fit.censor.impute = survfit(Surv(not_imputed_com$event_month, not_imputed_com$event_code=="No") ~ 1)
fup = quantile(fit.censor.impute, probs = c(0.25, 0.5, 0.75))$quantile
paste0("Median time of follow-up (IQR): ", fup[2], " [", fup[1], ",", fup[3], "]")

# by gender
fit.censor.gender = survfit(Surv(not_imputed_com$event_month, not_imputed_com$event_code=="No") ~ not_imputed_com$female)
fup.gender = quantile(fit.censor.gender, probs = c(0.25, 0.5, 0.75))$quantile; fup.gender
surv_diff_gender <- survdiff(Surv(event_month, event_code == "No") ~ female, data = not_imputed_com)

# by race and ethnicity
fit.censor.white = survfit(Surv(not_imputed_com$event_month, not_imputed_com$event_code=="No") ~ not_imputed_com$white)
fup.white = quantile(fit.censor.white, probs = c(0.25, 0.5, 0.75))$quantile; fup.white
surv_diff_white <- survdiff(Surv(event_month, event_code == "No") ~ white, data = not_imputed_com)

# by site
fit.censor.site = survfit(Surv(not_imputed_com$event_month, not_imputed_com$event_code=="No") ~ not_imputed_com$site)
fup.site = quantile(fit.censor.site, probs = c(0.25, 0.5, 0.75))$quantile; fup.site
surv_diff_site <- survdiff(Surv(event_month, event_code == "No") ~ site, data = not_imputed_com)

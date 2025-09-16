library(data.table)
library(dplyr)
library(mice)
library(lubridate)
library(eeptools)

# Reading in UPMC demographic data
demo <- fread("/ad_patient_demographics-datefix.dsv")
demo_old <- fread("/ad_demographics.dsv")
demo_diff <- demo[demo$PATIENT_STUDY_ID %in% demo_old$PATIENT_STUDY_ID]

## Obtaining age
demo_diff$dob <- format(as.Date(demo_diff$`TO_CHAR(BIRTH_DATE,'MM-DD-YYYY')`, format = "%m-%d-%Y"), "%Y-%m-%d")
demo_diff$age <- as.numeric(difftime(as.Date("2022-12-31"), demo_diff$dob, units = "weeks")) %/% 52.25

## Recoding gender
demo_diff$GENDER_TITLE <- as.factor(demo_diff$GENDER_TITLE)
demo_diff$gender <- case_when(demo_diff$GENDER_TITLE == "FEMALE" ~ 1,
                         demo_diff$GENDER_TITLE == "MALE" ~ 0,
                         demo_diff$GENDER_TITLE == "UNKNOWN" ~ NaN,
                         demo_diff$GENDER_TITLE == "" ~ NaN,
                         TRUE ~ NaN)
## Recoding race
demo_diff$RACE_TITLE <- as.factor(demo_diff$RACE_TITLE)
demo_diff$race <- factor(demo_diff$RACE_TITLE, levels =  c("", "ALASKA NATIVE", "AMERICAN INDIAN/ALASKA NATIVE", "ASIAN INDIAN", "BLACK/AFRICAN AMERICAN", "CHINESE", "FILIPINO", "GUAMANIAN", "JAPANESE", "KOREAN", "NATIVE HAWAIIAN", "NOT SPECIFIED", "OTHER", "OTHER ASIAN", "OTHER PACIFIC ISLANDER", "SAMOAN", "UNREPORTED,CHOSE NOT TO DISCLOSE RACE", "VIETNAMESE", "WHITE"),labels = c("Unknown", "American Indian or Alaska Native", "American Indian or Alaska Native", "Asian", "Black", "Asian", "Asian", "Native Hawaiian or Other Pacific Islander", "Asian", "Asian", "Native Hawaiian or Other Pacific Islander", "Unknown", "Other", "Asian", "Native Hawaiian or Other Pacific Islander", "Native Hawaiian or Other Pacific Islander", "Unknown", "Asian", "White"))
## Recoding ethnicity
demo_diff$ETHNIC_TITLE <- as.factor(demo_diff$ETHNIC_TITLE)
demo_diff$ethnicity <- factor(demo_diff$ETHNIC_TITLE, levels = c("", "ANOTHER HISPANIC, LATINO/A, OR SPANISH ORIGIN", "MEXICAN, MEXICAN AMERICAN, CHICANO/A", "NON-HISPANIC OR LATINO/A", "NOT SPECIFIED", "PUERTO RICAN", "UNREPORTED/CHOSE NOT TO DISCLOSE"), labels = c("Unknown", "Hispanic", "Hispanic", "Non Hispanic", "Unknown", "Hispanic", "Unknown"))

# Imputing race as a multinomial variable and ethnicity as a binary variable (Models 1 and 2)

## Including the following demographics for imputation (Model 1): age, gender
demo_impute_race_multinomial <- demo_diff[,c("PATIENT_STUDY_ID", "age", "gender", "race", "ethnicity")]
## Converting race and ethnicity to factors for imputation
demo_impute_race_multinomial$race <- case_when(demo_impute_race_multinomial$race == "Unknown" ~ NA,
                                               demo_impute_race_multinomial$race == "American Indian or Alaska Native" ~ 1, 
                                               demo_impute_race_multinomial$race == "Asian" ~ 2,
                                               demo_impute_race_multinomial$race == "Black" ~ 3,
                                               demo_impute_race_multinomial$race == "Native Hawaiian or Other Pacific Islander" ~ 4,
                                               demo_impute_race_multinomial$race == "White" ~ 5,
                                               demo_impute_race_multinomial$race == "Other" ~ 6,
                                               TRUE ~ NA_real_)
demo_impute_race_multinomial$race <- as.factor(demo_impute_race_multinomial$race)
demo_impute_race_multinomial$ethnicity <- case_when(demo_impute_race_multinomial$ethnicity == "Unknown" ~ NA,
                                                    demo_impute_race_multinomial$ethnicity == "Hispanic" ~ 0,
                                                    demo_impute_race_multinomial$ethnicity == "Non Hispanic" ~ 1,
                                                    TRUE ~ NA_real_)
demo_impute_race_multinomial$ethnicity <- as.factor(demo_impute_race_multinomial$ethnicity)

## Perform imputation
imp_model <- mice(demo_impute_race_multinomial[,-1], m = 5, maxit = 50, method = c("", "", "polyreg", "logreg"))
imputed_data <- complete(imp_model, 5)

## Including the following demographics for imputation (Model 2): age, gender, marital status, zip code
## Uncomment the following lines for alternative imputation model
# demo_impute_race_multinomial <- demo_diff[,c("PATIENT_STUDY_ID", "age", "gender", "race", "ethnicity", "MARITAL_STATUS_TITLE", "ZIP_CODE")]
# demo_impute_race_multinomial$race <- case_when(demo_impute_race_multinomial$race == "Unknown" ~ NA,
#                                                demo_impute_race_multinomial$race == "American Indian or Alaska Native" ~ 1, 
#                                                demo_impute_race_multinomial$race == "Asian" ~ 2,
#                                                demo_impute_race_multinomial$race == "Black" ~ 3,
#                                                demo_impute_race_multinomial$race == "Native Hawaiian or Other Pacific Islander" ~ 4,
#                                                demo_impute_race_multinomial$race == "White" ~ 5,
#                                                demo_impute_race_multinomial$race == "Other" ~ 6,
#                                                TRUE ~ NA_real_)
# demo_impute_race_multinomial$race <- as.factor(demo_impute_race_multinomial$race)
# demo_impute_race_multinomial$ethnicity <- case_when(demo_impute_race_multinomial$ethnicity == "Unknown" ~ NA,
#                                                     demo_impute_race_multinomial$ethnicity == "Hispanic" ~ 0,
#                                                     demo_impute_race_multinomial$ethnicity == "Non Hispanic" ~ 1,
#                                                     TRUE ~ NA_real_)
# demo_impute_race_multinomial$ethnicity <- as.factor(demo_impute_race_multinomial$ethnicity)
# imp_model <- mice(demo_impute_race_multinomial[,-1], m = 5, maxit = 50, method = c("", "", "polyreg", "logreg", "", ""))

## Obtaining output in desired format
imputed_data$`PATIENT_STUDY_ID` <- demo_impute_race_multinomial$PATIENT_STUDY_ID
imputed_data$race <- factor(imputed_data$race, levels = c(1:6), labels = c("American Indian or Alaska Native", "Asian",
                                                                           "Black", "Native Hawaiian or Other Pacific Islander",
                                                                           "White", "Other"))
imputed_data$ethnicity <- factor(imputed_data$ethnicity, levels = c(0:1), labels = c("Hispanic", "Non Hispanic"))
imputed_data <- imputed_data[,c("PATIENT_STUDY_ID", "race", "ethnicity")]

## Export the imputed data
demo_upmc_impute <- demo_diff[,-c("race", "ethnicity")] %>% left_join(imputed_data, by = "PATIENT_STUDY_ID")
write.csv(demo_upmc_impute, file = "/UPMC demographics imputed race ethnicity separate.csv")
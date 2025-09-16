library(data.table)
library(dplyr)
library(mice)

# MGB not imputed demographics
demo_mgb_noimpute <- fread("/MGB demographics deidentified.csv")
## Recoding gender
demo_mgb_noimpute$female <- case_when(demo_mgb_noimpute$Gender_Legal_Sex == "Female" ~ 1,
                         demo_mgb_noimpute$Gender_Legal_Sex == "Male" ~ 0,
                         demo_mgb_noimpute$Gender_Legal_Sex == "Unknown-U" ~ NaN,
                         TRUE ~ NaN)
## Recoding race
demo_mgb_noimpute$race <- case_when(demo_mgb_noimpute$Race_Group == "American Indian or Alaska Native" ~ 1,
                                demo_mgb_noimpute$Race_Group == "Asian" ~ 2,
                                demo_mgb_noimpute$Race_Group == "Black" ~ 3,
                                demo_mgb_noimpute$Race_Group == "Native Hawaiian or Other Pacific Islander" ~ 4,
                                demo_mgb_noimpute$Race_Group == "White" ~ 5,
                                demo_mgb_noimpute$Race_Group == "Other" ~ 6,
                                demo_mgb_noimpute$Race_Group == "Two or More" ~ 6,
                                demo_mgb_noimpute$Race_Group == "Unknown/Missing" ~ 7,
                                demo_mgb_noimpute$Race_Group == "Declined" ~ 7,
                                TRUE ~ NaN)
demo_mgb_noimpute$race <- factor(demo_mgb_noimpute$race, levels = c(1:7), labels = c("American Indian or Alaska Native", "Asian",
                                                           "Black", "Native Hawaiian or Other Pacific Islander",
                                                           "White", "Other", "Unknown"))
## Recoding ethicity
demo_mgb_noimpute$ethnicity <- case_when(demo_mgb_noimpute$Ethnic_Group == "HISPANIC" ~ 1,
                                         demo_mgb_noimpute$Ethnic_Group == "Non Hispanic" ~ 2,
                                         demo_mgb_noimpute$Ethnic_Group == "Unknown/Missing" ~ 3,
                                         demo_mgb_noimpute$Ethnic_Group == "DECLINED" ~ 3,
                                         TRUE ~ NaN)
demo_mgb_noimpute$ethnicity <- factor(demo_mgb_noimpute$ethnicity, levels = c(1:3), labels = c("Hispanic", "Non Hispanic", "Unknown"))

# Create binary race variables, setting Unknown to NA
demo_mgb_noimpute <- demo_mgb_noimpute %>%
  mutate(
    race_AmIndian_AlaskaNative = case_when(
      race == "American Indian or Alaska Native" ~ "Yes",
      race == "Unknown" ~ NA_character_,
      TRUE ~ "No"
    ),
    race_Asian = case_when(
      race == "Asian" ~ "Yes",
      race == "Unknown" ~ NA_character_,
      TRUE ~ "No"
    ),
    race_Black = case_when(
      race == "Black" ~ "Yes",
      race == "Unknown" ~ NA_character_,
      TRUE ~ "No"
    ),
    race_Hawaiian_PacificIslander = case_when(
      race == "Native Hawaiian or Other Pacific Islander" ~ "Yes",
      race == "Unknown" ~ NA_character_,
      TRUE ~ "No"
    ),
    race_White = case_when(
      race == "White" ~ "Yes",
      race == "Unknown" ~ NA_character_,
      TRUE ~ "No"
    ),
    race_Other = case_when(
      race == "Other" ~ "Yes",
      race == "Unknown" ~ NA_character_,
      TRUE ~ "No"
    )
  )

# Create binary ethnicity variables, setting Unknown to NA
demo_mgb_noimpute <- demo_mgb_noimpute %>%
  mutate(
    eth_Hispanic = case_when(ethnicity == "Hispanic" ~ "Yes", 
                             ethnicity == "Unknown" ~ NA_character_,
                             TRUE ~ "No"),
    eth_NonHispanic = case_when(ethnicity == "Non Hispanic" ~ "Yes", 
                                ethnicity == "Unknown" ~ NA_character_,
                                TRUE ~ "No")
    )

# Imputing race as a multinomial variable and ethnicity as a binary variable (Models 1 and 2)

## Including the following demographics for imputation (Model 1): age, gender
demo_impute_race_multinomial <- demo_mgb_noimpute[,c("PatientNum", "Age", "Gender_Legal_Sex", "race", "ethnicity")]
demo_impute_race_multinomial$race <- ifelse(demo_impute_race_multinomial$race == "Unknown", NA, demo_impute_race_multinomial$race)
demo_impute_race_multinomial$ethnicity <- ifelse(demo_impute_race_multinomial$ethnicity == "Unknown", NA, demo_impute_race_multinomial$ethnicity)

## Converting race and ethnicity to factors for imputation
demo_impute_race_multinomial$race <- as.factor(demo_impute_race_multinomial$race)
demo_impute_race_multinomial$ethnicity <- as.factor(demo_impute_race_multinomial$ethnicity)

## Perform imputation
imp_model <- mice(demo_impute_race_multinomial[,-1], m = 5, maxit = 50, method = c("", "", "polyreg", "logreg"))
imputed_data <- complete(imp_model, 5)

## Including the following demographics for imputation (Model 2): age, gender, marital status, religion, veteran, country
## Uncomment the following lines for alternative imputation model
# demo_impute_race_multinomial <- demo_mgb_noimpute[,c("PatientNum", "Age", "Gender_Legal_Sex", "race", "ethnicity", "Marital_status", "Religion", "Is_a_veteran", "Country")]
# demo_impute_race_multinomial$race <- ifelse(demo_impute_race_multinomial$race == "Unknown", NA, demo_impute_race_multinomial$race)
# demo_impute_race_multinomial$ethnicity <- ifelse(demo_impute_race_multinomial$ethnicity == "Unknown", NA, demo_impute_race_multinomial$ethnicity)
# imp_model <- mice(demo_impute_race_multinomial[,-1], m = 5, maxit = 50, method = c("", "", "polyreg", "logreg", "", "", "", ""))

## Obtaining output in desired format
imputed_data$PatientNum <- demo_impute_race_multinomial$PatientNum
imputed_data$race <- factor(imputed_data$race, levels = c(1:6), labels = c("American Indian or Alaska Native", "Asian",
                                                                   "Black", "Native Hawaiian or Other Pacific Islander",
                                                                   "White", "Other"))
imputed_data$ethnicity <- factor(imputed_data$ethnicity, levels = c(1:2), labels = c("Hispanic", "Non Hispanic"))
imputed_data <- imputed_data[,c("PatientNum", "race", "ethnicity")]

## Export the imputed data
demo_mgb_impute <- demo_mgb_noimpute[,-c("race", "ethnicity")] %>% left_join(imputed_data, by = "PatientNum")
write.csv(demo_mgb_impute, file = "/MGB demographics imputed race ethnicity separate.csv")

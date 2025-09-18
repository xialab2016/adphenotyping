library(data.table)
library(dplyr)
library(meta)
library(ggplot2)
library(adjustedCurves)
library(patchwork)
library(survminer)
library(survival)

# Read in imputed data (or single comorbidity score, not imputed data)
upmc_nh <- read.csv("/adjusted_HR_nh_revised.csv")
upmc_death <- read.csv("/adjusted_HR_death_revised.csv")
mgb_nh <- read.csv("/adjusted_HR_nh_revised.csv")
mgb_death <- read.csv("/adjusted_HR_death_revised.csv")

# Keeping only 24 month data from codified or NLP prediction
upmc_nh <- upmc_nh %>% filter(base_type == "24 months", nh_type == "codified & NLP")
upmc_death <- upmc_death %>% filter(base_type == "24 months", nh_type == "codified & NLP")
mgb_nh <- mgb_nh %>% filter(base_type == "24 months", nh_type == "codified & NLP")
mgb_death <- mgb_death %>% filter(base_type == "24 months", nh_type == "codified & NLP")

# Adding in sample sizes
upmc_nh$N <- 13408
upmc_death$N <- 13408
mgb_nh$N <- 15854
mgb_death$N <- 15854

# Adding in site IDs
upmc_nh$study <- "UPMC"
upmc_death$study <- "UPMC"
mgb_nh$study <- "MGB"
mgb_death$study <- "MGB"

# Adding in feature names
upmc_nh$feature_names <- c("race/ethnicity: non-Hispanic white", "gender: women", "age at AD diagnosis", "healthcare utilization: log (util + 1)", 
                           "congestive heart failure", "cardiac arrhythmias", "valvular disease", "pulmonary circulation disorders", 
                           "peripheral vascular disorder", "hypertension, uncomplicated", "hypertension, complicated",
                           "paralysis", "other neurological disorders", "chronic pulmonary disease", "diabetes, uncontrolled", 
                           "diabetes, controlled", "hypothyroidism", "renal failure", "liver disease", 
                           "peptic ulcer disease", "AIDS/HIV", "lymphoma", "metastatic cancer", "solid tumor, without metastasis", 
                           "rheumatoid arthritis", "coagulopathy", "obesity", "weight loss", "fluid and electrolyte disorders", 
                           "blood loss anemia", "deficiency anemia", "alcohol abuse", "drug abuse", "psychoses", "depression")
upmc_death$feature_names <- c("race/ethnicity: non-Hispanic white", "gender: women", "age at AD diagnosis", "healthcare utilization: log (util + 1)", 
                              "congestive heart failure", "cardiac arrhythmias", "valvular disease", "pulmonary circulation disorders", 
                              "peripheral vascular disorder", "hypertension, uncomplicated", "hypertension, complicated",
                              "paralysis", "other neurological disorders", "chronic pulmonary disease", "diabetes, uncontrolled", 
                              "diabetes, controlled", "hypothyroidism", "renal failure", "liver disease",
                              "peptic ulcer disease", "AIDS/HIV", "lymphoma", "metastatic cancer", "solid tumor, without metastasis", 
                              "rheumatoid arthritis", "coagulopathy", "obesity", "weight loss", "fluid and electrolyte disorders", 
                              "blood loss anemia", "deficiency anemia", "alcohol abuse", "drug abuse", "psychoses", "depression")
mgb_nh$feature_names <- c("race/ethnicity: non-Hispanic white", "gender: women", "age at AD diagnosis", "healthcare utilization: log (util + 1)", 
                          "congestive heart failure", "cardiac arrhythmias", "valvular disease", "pulmonary circulation disorders", 
                          "peripheral vascular disorder", "hypertension, uncomplicated", "hypertension, complicated",
                          "paralysis", "other neurological disorders", "chronic pulmonary disease", "diabetes, uncontrolled", 
                          "diabetes, controlled", "hypothyroidism", "renal failure", "liver disease", 
                          "peptic ulcer disease", "AIDS/HIV", "lymphoma", "metastatic cancer", "solid tumor, without metastasis", 
                          "rheumatoid arthritis", "coagulopathy", "obesity", "weight loss", "fluid and electrolyte disorders", 
                          "blood loss anemia", "deficiency anemia", "alcohol abuse", "drug abuse", "psychoses", "depression")
mgb_death$feature_names <- c("race/ethnicity: non-Hispanic white", "gender: women", "age at AD diagnosis", "healthcare utilization: log (util + 1)", 
                             "congestive heart failure", "cardiac arrhythmias", "valvular disease", "pulmonary circulation disorders", 
                             "peripheral vascular disorder", "hypertension, uncomplicated", "hypertension, complicated",
                             "paralysis", "other neurological disorders", "chronic pulmonary disease", "diabetes, uncontrolled", 
                             "diabetes, controlled", "hypothyroidism", "renal failure", "liver disease", 
                             "peptic ulcer disease", "AIDS/HIV", "lymphoma", "metastatic cancer", "solid tumor, without metastasis", 
                             "rheumatoid arthritis", "coagulopathy", "obesity", "weight loss", "fluid and electrolyte disorders", 
                             "blood loss anemia", "deficiency anemia", "alcohol abuse", "drug abuse", "psychoses", "depression")

# Calculating standard errors of HR
upmc_nh$se <- (log(upmc_nh$HR_upper) - log(upmc_nh$HR_lower)) / (2 * 1.96)
upmc_death$se <- (log(upmc_death$HR_upper) - log(upmc_death$HR_lower)) / (2 * 1.96)
mgb_nh$se <- (log(mgb_nh$HR_upper) - log(mgb_nh$HR_lower)) / (2 * 1.96)
mgb_death$se <- (log(mgb_death$HR_upper) - log(mgb_death$HR_lower)) / (2 * 1.96)

# Combining data frames
nh_HR_all <- rbind(upmc_nh, mgb_nh)
death_HR_all <- rbind(upmc_death, mgb_death)

# Fixed effects meta-analysis using metagen

# Group by feature names
## Run this with nh_HR_all or death_HR_all
meta_results <- nh_HR_all %>%
  group_by(feature_names) %>%
  group_modify(~ {
    meta_analysis <- metagen(
      TE = log(.x$HR),                      # Log-transformed hazard ratios
      seTE = .x$se,                         # Standard errors
      studlab = .x$study, # Labels
      com.fixed = T,                       # Fixed-effects model
      random = F,                   
      sm = "HR"                             # Using HR
    )
    tibble(meta_analysis = list(meta_analysis))
  })

# Sort meta_results
## Run this with nh_HR_all or death_HR_all
sorted_meta_results <- nh_HR_all %>%
  dplyr::select(feature_names) %>%
  distinct() %>%
  left_join(meta_results, by = "feature_names")

# Extract the relevant statistics from each meta-analysis and compile into a data frame
meta_summary <- sorted_meta_results %>%
  rowwise() %>%
  mutate(
    Fixed_Effects_HR = exp(meta_analysis$TE.fixed), 
    CI_lower = exp(meta_analysis$lower.fixed),
    CI_upper = exp(meta_analysis$upper.fixed),
    p_value = meta_analysis$pval.fixed,
    Q_statistic = meta_analysis$Q,
    Q_df = meta_analysis$df.Q,
    Q_p_value = meta_analysis$pval.Q,
    tau2 = meta_analysis$tau^2,
    tau = meta_analysis$tau,
    I2 = meta_analysis$I2
  ) %>%
  dplyr::select(feature_names, Fixed_Effects_HR, CI_lower, CI_upper, p_value, Q_statistic, Q_df, Q_p_value, tau2, tau, I2)

pval_sig = Vectorize(function(p) {
  if ((p <= 0.05)& (p >0.01)) {return ("*")}
  else if ((p <= 0.01)& (p >0.001)) {return ("**")}
  else if ((p <= 0.001)& (p >0.0001)) {return("***")}
  else if (p <= 0.0001) {return("****")}
  else {return("")}
})
pval_sig_level = Vectorize(function(p) {
  if ((p <= 0.05)& (p >0.01)) {return ("4")}
  else if ((p <= 0.01)& (p >0.001)) {return ("3")}
  else if ((p <= 0.001)& (p >0.0001)) {return("2")}
  else if (p <= 0.0001) {return("1")}
  else {return("")}
})

meta_summary$sig = as.factor(ifelse(meta_summary$p_value<=0.05, 1, 0))

meta_summary <- meta_summary %>%
  mutate(
    pval.sig = pval_sig(p_value),
    pval.sig.level = pval_sig_level(p_value),
    feature_names = factor(feature_names, levels = feature_names),
    text = paste0(
      formatC(Fixed_Effects_HR, format = "f", digits = 3), " (",
      formatC(CI_lower, format = "f", digits = 3), ", ",
      formatC(CI_upper, format = "f", digits = 3), ")"
    )
  )
meta_summary <- meta_summary %>%
  mutate(pval.sig.level = as.numeric(pval.sig.level)) %>%
  arrange(pval.sig.level, desc(Fixed_Effects_HR)) %>%
  # Update feature_names factor levels to match row order
  mutate(
    feature_names = factor(feature_names, levels = feature_names),
    pval.sig = pval_sig(p_value),  # keep your existing pval.sig function
    text = sprintf("%.3f (%.3f, %.3f)", Fixed_Effects_HR, CI_lower, CI_upper)
  )

plot.HR = ggplot(meta_summary, aes(x=forcats::fct_rev(feature_names), y=Fixed_Effects_HR, color=sig)) +
  scale_color_manual(values = c("0" = "black", "1" = "red")) +
  guides(color = "none") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    width = .1,
    position = position_dodge(width = 0.5),
    aes(ymin=CI_lower, ymax=CI_upper)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  geom_text(aes(y = 3.3, label = pval.sig), color = "red", size = 6) +   
  geom_text(aes(y = 2.8, label = text), nudge_x = 0, size = 5) +      
  ylab("Adjusted Hazard Ratio") +
  xlab("Baseline Features") +
  coord_flip() +
  scale_y_continuous(limits = c(0,3.4), expand = c(0,0)) +          
  theme_bw() + 
  ggtitle("Nursing Home Admission") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18)
  )
plot.HR

save(sorted_meta_results, file = "/sorted_meta_results_nh.RData")
write.csv(meta_summary, file="/fixed_eff_HR_nh.csv")
ggsave(plot.HR, file="/fixed_eff_HR_nh.png", width=11, height=9, dpi=300)

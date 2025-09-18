library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# Custom functions for p-value significance
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


# Load combined patient-level survival analysis data
## Load either imputed or not imputed data
load("/imputed_com.RData") # From Combined_demographics.R

# Comorbidity burden
features.nm = c("white", "female", "age", "log(util+1)", "chf", "carit", "valv", "pcd", "pvd", "hypunc", "hypc", "para", "ond",
                "cpd", "diabunc", "diabc", "hypothy", "rf", "ld", "pud", "aids", "lymph",
                "metacanc", "solidtum", "rheumd", "coag", "obes", "wloss", "fed", "blane", 
                "dane", "alcohol", "drug", "psycho", "depre")
features.df = data.frame(ac = features.nm, nm = c("race/ethnicity: non-Hispanic white", "gender: women", "age at AD diagnosis", "healthcare utilization: log (util + 1)", 
                                                  "congestive heart failure", "cardiac arrhythmias", "valvular disease", "pulmonary circulation disorders", 
                                                  "peripheral vascular disorder", "hypertension, uncomplicated", "hypertension, complicated",
                                                  "paralysis", "other neurological disorders", "chronic pulmonary disease", "diabetes, uncontrolled", 
                                                  "diabetes, controlled", "hypothyroidism", "renal failure", "liver disease", 
                                                  "peptic ulcer disease", "AIDS/HIV", "lymphoma", "metastatic cancer", "solid tumor, without metastasis", 
                                                  "rheumatoid arthritis", "coagulopathy", "obesity", "weight loss", "fluid and electrolyte disorders", 
                                                  "blood loss anemia", "deficiency anemia", "alcohol abuse", "drug abuse", "psychoses", "depression"))
# Pre-existing commorbidities and risk of AD diagnosis (calculate RR for all, by race, by sex)
race.prop = imputed_com %>% group_by(white) %>% summarise(n = n(), across(chf_2yr:depre_2yr, ~sum(.))) %>% gather(features, count, chf_2yr:depre_2yr)
race.prop = race.prop %>% mutate(prop = count / n, 
                                 lower = prop-1.96*sqrt(prop*(1-prop)/n), 
                                 upper = prop+1.96*sqrt(prop*(1-prop)/n))

race.pval = race.prop %>% group_by(features) %>% summarise(pval = prop.test(x=count, n=n)$p.value) %>% mutate(pval.sig = pval_sig(pval))
race.prop = race.prop %>% left_join(race.pval, by="features") %>% 
  mutate(race = ifelse(white==1, "Non-Hispanic White", "Others"),
         features = gsub("_2yr", "", features)) %>% 
  left_join(features.df, by=c("features"="ac"))
race.prop
base.race = ggplot(race.prop, aes(x=reorder(nm,prop), y=prop)) +
  geom_point(aes(color=as.factor(race), group=as.factor(race)), position=position_dodge(width = 0.5)) +
  geom_errorbar(width=.1, position=position_dodge(width = 0.5), aes(color=as.factor(race), group=as.factor(race), ymin=lower, ymax=upper)) +
  geom_text(aes(x=nm, y=0.95, label=pval.sig),color="red") +
  ylim(0, 1) +
  xlab("Comorbidity") + ylab("% with comorbidity") +
  labs(color = "Race") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Baseline Comorbidities by Race/Ethnicity") +
  scale_color_manual(name = "Race/Ethnicity", values=c("#2E9FDF","#E7B800"))

base.race.age = imputed_com %>% group_by(white) %>% summarise(n = n(), mean_age = mean(age), sd_age = sd(age), race = ifelse(white==1, "Non-Hisp White", "Others")) %>% 
  ggplot() + geom_point(aes(x=as.factor(race), y=mean_age, color=as.factor(race), group=as.factor(race)), position=position_dodge(width = 0.5)) +
  geom_errorbar(width=.1, position=position_dodge(width = 0.5), aes(x=as.factor(race), color=as.factor(race), group=as.factor(race), ymin=mean_age-1.96*sd_age/sqrt(n), 
                                                                    ymax=mean_age+1.96*sd_age/sqrt(n))) +
  xlab("age") + ylab("years") +
  coord_flip() +
  theme_bw() +
  theme(legend.position="none") +
  scale_color_manual(values=c("#2E9FDF","#E7B800"))

gender.prop = imputed_com %>% group_by(female) %>% summarise(n = n(), across(chf_2yr:depre_2yr, ~sum(.))) %>% gather(features, count, chf_2yr:depre_2yr)
gender.prop = gender.prop %>% mutate(prop = count / n, 
                                     lower = prop-1.96*sqrt(prop*(1-prop)/n), 
                                     upper = prop+1.96*sqrt(prop*(1-prop)/n))

gender.pval = gender.prop %>% group_by(features) %>% summarise(pval = prop.test(x=count, n=n)$p.value) %>% mutate(pval.sig = pval_sig(pval))
gender.prop = gender.prop %>% left_join(gender.pval, by="features")%>% 
  mutate(gender = factor(ifelse(female==1, "Women", "Men"), levels=c("Women", "Men")),
         features = gsub("_2yr", "", features)) %>% 
  left_join(features.df, by=c("features"="ac"))
base.gender = ggplot(gender.prop, aes(x=reorder(nm,prop), y=prop)) +
  geom_point(aes(color=as.factor(gender), group=as.factor(gender)), position=position_dodge(width = 0.5)) +
  geom_errorbar(width=.1, position=position_dodge(width = 0.5), aes(color=as.factor(gender), group=as.factor(gender), ymin=lower, ymax=upper)) +
  geom_text(aes(x=nm, y=0.85, label=pval.sig),color="red") +
  ylim(0, 0.9) +
  xlab("Comorbidity") + ylab("% with comorbidity") +
  labs(color = "Gender") +
  coord_flip() +
  theme_bw() + 
  ggtitle("Baseline Comorbidities by Gender") 

base.plot = ggarrange(base.race, base.gender, widths = c(0.54, 0.46), nrow=1)
base.plot
ggsave(base.plot, width=12, height=6, dpi=300, file="/baseline_comorbidities_imputed.png")

# Adjusted survival curve (adjusting for both institutions)

## Load either imputed or not imputed data
survival_df <- imputed_com[,c("patient_num", "nh_month", "nh", "death_month", "death",
                              "white", "female", "age", "codect_2yr",
                              "chf_2yr", "carit_2yr", "valv_2yr", "pcd_2yr", "pvd_2yr",
                              "hypunc_2yr", "hypc_2yr", "para_2yr", "ond_2yr", "cpd_2yr",
                              "diabunc_2yr", "diabc_2yr", "hypothy_2yr", "rf_2yr", "ld_2yr",
                              "pud_2yr", "aids_2yr", "lymph_2yr", "metacanc_2yr", "solidtum_2yr",
                              "rheumd_2yr", "coag_2yr", "obes_2yr", "wloss_2yr", "fed_2yr",
                              "blane_2yr", "dane_2yr", "alcohol_2yr", "drug_2yr", "psycho_2yr",
                              "depre_2yr")]

# Mean of all covariates from both institutions
df.mean.race = survival_df %>% group_by(white) %>% summarise(across(female:depre_2yr, ~ mean(.x, na.rm = TRUE)))
df.mean.race

df.mean.gender = survival_df %>% group_by(female) %>% summarise(across(white:depre_2yr, ~ mean(.x, na.rm = TRUE)))
df.mean.gender

# Sample size for imputed
upmc.n = 13408
mgb.n = 15854

# Nursing home admission

load("/sorted_meta_results_nh.RData") # from Fixed_effect_meta-analysis_adjusted_HR.R

# weights for race
survival_df[1:upmc.n, "weight_race"] = weights(sorted_meta_results$meta_analysis[[1]])$p.common[1] / 100
survival_df[(upmc.n+1):(upmc.n+mgb.n), "weight_race"] = weights(sorted_meta_results$meta_analysis[[1]])$p.common[2] / 100

# weights for gender
survival_df[1:upmc.n, "weight_gender"] = weights(sorted_meta_results$meta_analysis[[2]])$p.common[1] / 100
survival_df[(upmc.n+1):(upmc.n+mgb.n), "weight_gender"] = weights(sorted_meta_results$meta_analysis[[2]])$p.common[2] / 100

# a. by race
CS.nh.race = coxph(Surv(nh_month, nh) ~ strata(white) + female + age + log(codect_2yr+1) +
                     chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+
                     cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                     metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                     dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data=survival_df, x=TRUE, weights = survival_df$weight_race)
cox.nh.race = survfit(CS.nh.race, df.mean.race)

g1p = ggsurvplot(cox.nh.race, data = df.mean.race, size = 1, palette = c("#E7B800", "#2E9FDF"),
                 conf.int = TRUE, pval = F, xlim=c(-10, 160),
                 legend.labs = c("Others", "Non-Hispanic White"),
                 legend.title = "",  # Remove title
                 risk.table = F,
                 ggtheme = theme_bw(), xlab="Months after AD diagnosis", break.time.by = 50, ylab="Adjusted\nSurvival\nProbability",
                 title="Nursing Home Admission (By Race/Ethnicity)")$plot + 
  theme(axis.title.y = element_text(angle = 0, vjust = 1.05, hjust = 1.1, size = 14), # Y-axis title size
        axis.text.x = element_text(size = 14),  # X-axis tick labels
        axis.text.y = element_text(size = 14),  # Y-axis tick labels
        axis.title.x = element_text(size = 14), # X-axis title
        plot.title = element_text(size = 18), # Plot title size
        legend.position = "inside",
        legend.position.inside = c(0.72, 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.margin = margin(5, 5, 5, 5),
        legend.box.spacing = unit(0, "pt"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))
fit.nh.race = survfit(Surv(time = nh_month, nh) ~ white, data=survival_df)
g1t = ggsurvtable(fit.nh.race, risk.table.type = "nrisk_cumevents", data=survival_df,color = "strata", palette = c("#E7B800", "#2E9FDF"), y.text=F,
                  legend="none", 
                  ggtheme = theme_bw() + 
                    theme(axis.text.x = element_text(size = 14),  # X-axis tick labels
                          axis.text.y = element_text(size = 14),  # Y-axis tick labels
                          plot.title = element_text(size = 16)), # Title
                  break.time.by =50, xlim=c(-10, 160),
                  fontsize = 4.5,
                  title = "Number of events")$risk.table +
  xlab("") + ylab("") + 
  theme(legend.position = "none")
g1.all = g1p + g1t + plot_layout(ncol = 1, heights=c(4,1)); g1.all

# b. by gender
CS.nh.gender = coxph(Surv(nh_month,nh) ~ white + strata(female) + age + log(codect_2yr+1) +
                       chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+
                       cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                       metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                       dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data=survival_df, x=TRUE, weights = survival_df$weight_gender)
cox.nh.gender = survfit(CS.nh.gender, df.mean.gender)

g2p = ggsurvplot(cox.nh.gender, data = df.mean.gender, size = 1, palette = c("#00BFC4", "#F8766D"),
                 conf.int = TRUE, pval = F, xlim=c(-10, 160),
                 legend.labs = c("Men", "Women"),
                 legend.title = "",  # Remove title
                 risk.table = F,
                 ggtheme = theme_bw(), xlab="Months after AD diagnosis", break.time.by = 50, ylab="Adjusted\nSurvival\nProbability",
                 title="Nursing Home Admission (By Gender)")$plot + 
  theme(axis.title.y = element_text(angle = 0, vjust = 1.05, hjust = 1.1, size = 14), # Y-axis title size
        axis.text.x = element_text(size = 14),  # X-axis tick labels
        axis.text.y = element_text(size = 14),  # Y-axis tick labels
        axis.title.x = element_text(size = 14), # X-axis title
        plot.title = element_text(size = 18), # Plot title size
        legend.position = "inside",
        legend.position.inside = c(0.72, 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.margin = margin(5, 5, 5, 5),
        legend.box.spacing = unit(0, "pt"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))
fit.nh.gender = survfit(Surv(time = nh_month, nh) ~ female, data=survival_df)
g2t = ggsurvtable(fit.nh.gender, risk.table.type = "nrisk_cumevents", data=survival_df,color = "strata", palette = c("#00BFC4", "#F8766D"), y.text=F,
                  legend="none", 
                  ggtheme = theme_bw() + 
                    theme(axis.text.x = element_text(size = 14),  # X-axis tick labels
                          axis.text.y = element_text(size = 14),  # Y-axis tick labels
                          plot.title = element_text(size = 16)), # Title
                  break.time.by =50, xlim=c(-10, 160),
                  fontsize = 4.5,
                  title = "Number of events")$risk.table +
  xlab("") + ylab("") + 
  theme(legend.position = "none")
g2.all = g2p + g2t + plot_layout(ncol = 1, heights=c(4,1)); g2.all


# Death

load("/sorted_meta_results_death.RData") # from Fixed_effect_meta-analysis_adjusted_HR.R

# weights for race
survival_df[1:upmc.n, "weight_race"] = weights(sorted_meta_results$meta_analysis[[1]])$p.common[1] / 100
survival_df[(upmc.n+1):(upmc.n+mgb.n), "weight_race"] = weights(sorted_meta_results$meta_analysis[[1]])$p.common[2] / 100

# weights for gender
survival_df[1:upmc.n, "weight_gender"] = weights(sorted_meta_results$meta_analysis[[2]])$p.common[1] / 100
survival_df[(upmc.n+1):(upmc.n+mgb.n), "weight_gender"] = weights(sorted_meta_results$meta_analysis[[2]])$p.common[2] / 100

# c. by race
CS.death.race = coxph(Surv(death_month,death) ~ strata(white) + female + age + log(codect_2yr+1) +
                        chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+
                        cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                        metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                        dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data=survival_df, x=TRUE, weights = survival_df$weight_race)
cox.death.race = survfit(CS.death.race, df.mean.race)

g3p = ggsurvplot(cox.death.race, data = df.mean.race, size = 1, palette = c("#E7B800", "#2E9FDF"),
                 conf.int = TRUE, pval = F, xlim=c(-10, 160),
                 legend.labs = c("Others", "Non-Hispanic White"),
                 legend.title = "",  # Remove title
                 risk.table = F,
                 ggtheme = theme_bw(), xlab="Months after AD diagnosis", break.time.by = 50, ylab="Adjusted\nSurvival\nProbability",
                 title="Death (By Race/Ethnicity)")$plot + 
  theme(axis.title.y = element_text(angle = 0, vjust = 1.05, hjust = 1.1, size = 14), # Y-axis title size
        axis.text.x = element_text(size = 14),  # X-axis tick labels
        axis.text.y = element_text(size = 14),  # Y-axis tick labels
        axis.title.x = element_text(size = 14), # X-axis title
        plot.title = element_text(size = 18), # Plot title size
        legend.position = "inside",
        legend.position.inside = c(0.72, 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.margin = margin(5, 5, 5, 5),
        legend.box.spacing = unit(0, "pt"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))
fit.death.race = survfit(Surv(time = death_month, death) ~ white, data=survival_df)
g3t = ggsurvtable(fit.death.race, risk.table.type = "nrisk_cumevents", data=survival_df,color = "strata", palette = c("#E7B800", "#2E9FDF"), y.text=F,
                  legend="none", 
                  ggtheme = theme_bw() + 
                    theme(axis.text.x = element_text(size = 14),  # X-axis tick labels
                          axis.text.y = element_text(size = 14),  # Y-axis tick labels
                          plot.title = element_text(size = 16)), # Title
                  break.time.by =50, xlim=c(-10, 160),
                  fontsize = 4.5,
                  title = "Number of events")$risk.table +
  xlab("") + ylab("") + 
  theme(legend.position = "none")
g3.all = g3p + g3t + plot_layout(ncol = 1, heights=c(4,1)); g3.all


# d. by gender
CS.death.gender = coxph(Surv(death_month,death) ~ white + strata(female) + age + log(codect_2yr+1) +
                          chf_2yr + carit_2yr + valv_2yr + pcd_2yr + pvd_2yr + hypunc_2yr + hypc_2yr + para_2yr + ond_2yr+
                          cpd_2yr + diabunc_2yr + diabc_2yr + hypothy_2yr + rf_2yr + ld_2yr + pud_2yr + aids_2yr + lymph_2yr +
                          metacanc_2yr + solidtum_2yr + rheumd_2yr + coag_2yr + obes_2yr + wloss_2yr + fed_2yr + blane_2yr +
                          dane_2yr + alcohol_2yr + drug_2yr + psycho_2yr + depre_2yr, data=survival_df, x=TRUE, weights = survival_df$weight_gender)
cox.death.gender = survfit(CS.death.gender, df.mean.gender)

g4p = ggsurvplot(cox.death.gender, data = df.mean.gender, size = 1, palette = c("#00BFC4", "#F8766D"),
                 conf.int = TRUE, pval = F, xlim=c(-10, 160),
                 legend.labs = c("Men", "Women"),
                 legend.title = "",  # Remove title
                 risk.table = F,
                 ggtheme = theme_bw(), xlab="Months after AD diagnosis", break.time.by = 50, ylab="Adjusted\nSurvival\nProbability",
                 title="Death (By Gender)")$plot + 
  theme(axis.title.y = element_text(angle = 0, vjust = 1.05, hjust = 1.1, size = 14), # Y-axis title size
        axis.text.x = element_text(size = 14),  # X-axis tick labels
        axis.text.y = element_text(size = 14),  # Y-axis tick labels
        axis.title.x = element_text(size = 14), # X-axis title
        plot.title = element_text(size = 18), # Plot title size
        legend.position = "inside",
        legend.position.inside = c(0.72, 0.8),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.key.size = unit(1.2, "cm"),
        legend.spacing.x = unit(0.3, "cm"),
        legend.margin = margin(5, 5, 5, 5),
        legend.box.spacing = unit(0, "pt"),
        legend.background = element_rect(fill = "white", color = "black", size = 0.5))

fit.death.gender = survfit(Surv(time = death_month, death) ~ female, data=survival_df)
g4t = ggsurvtable(fit.death.gender, risk.table.type = "nrisk_cumevents", data=survival_df,color = "strata", palette = c("#00BFC4", "#F8766D"), y.text=F,
                  legend="none", 
                  ggtheme = theme_bw() + 
                    theme(axis.text.x = element_text(size = 14),  # X-axis tick labels
                          axis.text.y = element_text(size = 14),  # Y-axis tick labels
                          plot.title = element_text(size = 16)), # Title
                  break.time.by =50, xlim=c(-10, 160),
                  fontsize = 4.5,
                  title = "Number of events")$risk.table +
  xlab("") + ylab("") + 
  theme(legend.position = "none")
g4.all = g4p + g4t + plot_layout(ncol = 1, heights=c(4,1)); g4.all

survival.plot = ggarrange(g1.all, g2.all, g3.all, g4.all, widths = c(1.5, 1.5), heights = c(1.7, 1.7), ncol=2, nrow=2)
survival.plot

ggsave("/adjKM_nh_death.png", survival.plot, width = 15, height = 17, dpi = 300)
ggsave(g1.all, file=paste0("/adjKM_nh_race.png"), width=7, height=8, dpi=300, device="png")
ggsave(g2.all, file=paste0("/adjKM_nh_gender.png"), width=7, height=8, dpi=300, device="png")
ggsave(g3.all, file=paste0("/adjKM_death_race.png"), width=7, height=8, dpi=300, device="png")
ggsave(g4.all, file=paste0("/adjKM_death_gender.png"), width=7, height=8, dpi=300, device="png")





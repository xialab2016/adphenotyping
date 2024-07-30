library(openxlsx)
library(glmnet)
library(PRROC)
library(ranger)
library(MAP)
library(data.table)
library(stringr)
library(KOMAP)
library(PheNorm)
library(dplyr)

setwd("your_directory") 

cod_ONCE = read.csv("\ONCE\\ONCE_AD_cod_upmc.csv", sep="|")
cod_ONCE = cod_ONCE[cod_ONCE$high_confidence_level==1,]$Variable
codify_select = gsub(cod_ONCE, pattern = "\\:", replacement = ".")

NLP_ONCE = read.csv("\ONCE\\ONCE_AD_nlp_upmc.csv", sep="|")
cui_select_weight = NLP_ONCE[NLP_ONCE$high_confidence_level==1,]$cui

all_dat = read.csv("\data quality check\\AD_newdat_0708.csv")[,-1]



dat_log_wide = all_dat
dat_log_wide[,-1] <- log(all_dat[, -1] + 1)

Pheno_code_lst <- c("PheCode.290.11")
main_CUI_lst <- c('C0002395')
code.label <- c('AD')
full_nm_lst <- c('alzheimer\'s disease')


auc_diff_table = auc_boot_mean_table = c()
p_auc_diff_table = c()
auc_diff_lst = rep(list(c()), 1)
auc_all_lst = rep(list(c()), 1)
auc_table <- c()
F_table <- c()
PRAUC_table <- c()
prev_vec <- c()
label_num_vec <- c()
coef.lst <- c()

lambda <- 0.15
rep <- 50
sparse = TRUE
nm.utl = 'utils'
args = commandArgs(trailingOnly=TRUE)
k_vec = as.numeric(args[1])

auc_table_w = auc_table_b = c()

main_code <- Pheno_code_lst[1] 
main_CUI <- main_CUI_lst[1]
phename <- code.label[1]
disease = full_nm_lst[1]

dat_corrupt = dat_log_wide

dat_corrupt$corrupt_mainICD = dat_corrupt[, main_code]
dat_corrupt$corrupt_mainNLP = dat_corrupt[, main_CUI]
dat_corrupt$corrupt_mainICD[sample(1:nrow(dat_corrupt), round(nrow(dat_corrupt) * 0.2), replace = FALSE)] = mean(as.matrix(dat_corrupt[, main_code]))

id.train = sample(1:nrow(dat_corrupt), round(nrow(dat_corrupt) / 3 * 2))
id.valid = setdiff(1:nrow(dat_corrupt), id.train)
dat.cov.train = cov(dat_corrupt[id.train, ])
dat.cov.valid = cov(dat_corrupt[id.valid, ])
dat.cov = cov(dat_corrupt[,colnames(dat_corrupt) %in% c(c(main_code, main_CUI, "corrupt_mainICD", "corrupt_mainNLP", codify_select, cui_select_weight), "utils")])


##-------------------------------------------------------------------------------------------------------------
# gold labels

dat.label <- read.csv("AD\\AD_phenotyping_yunqing\\label_bin.csv")[,-1]

input.cov.train <- dat.cov.train
input.cov.valid <- dat.cov.valid
input.cov <- dat.cov
nm.corrupt.code = 'corrupt_mainICD'
nm.corrupt.cui = 'corrupt_mainNLP'


dat.label$pi = 1
colnames(dat.label) <- c('PatientNum','Y', 'pi')
dat.label$PatientNum = as.character(dat.label$PatientNum)


dat.merge = merge(dat.label, dat_log_wide, by="PatientNum")
colnames(dat.merge)[which(names(dat.merge)=='PatientNum')] = "patient_num"


#prediction on labeled data
out.corrupt = KOMAP_corrupt(input.cov.train, input.cov.valid, is.wide = TRUE,
                            target.code = main_code, target.cui = main_CUI, nm.disease = phename,
                            nm.utl, nm.corrupt.code, nm.corrupt.cui, nm.multi = NULL, 
                            codify.feature = codify_select, nlp.feature = cui_select_weight,
                            pred = T, eval.real = FALSE, eval.sim = FALSE,
                            # mu0 = NULL, mu1 = NULL, var0 = NULL, var1 = NULL, prev_Y = NULL, B = 10000,
                            dat.part = dat.merge, nm.id = 'patient_num', nm.pi = 'pi', nm.y = 'Y'
)


#prediction on full data

out.corrupt.pred = KOMAP_corrupt(input.cov.train, input.cov.valid, is.wide = TRUE,
                            target.code = main_code, target.cui = main_CUI, nm.disease = phename,
                            nm.utl, nm.corrupt.code, nm.corrupt.cui, nm.multi = NULL, 
                            codify.feature = codify_select, nlp.feature = cui_select_weight,
                            pred = T, eval.real = F, eval.sim = FALSE,
                            dat.part = dat_log_wide, nm.id = 'PatientNum'#, nm.pi = 'pi', nm.y = 'Y'
)

#patient count with 90% and 95% specificity

scores = out.corrupt.pred$pred_prob$pred.score
scores$clust.90 = 0
scores$clust.90[scores$`mainICDNLP + codify & NLP`>=0.1941] = 1 

scores$clust.95 = 0
scores$clust.95[scores$`mainICDNLP + codify & NLP`>=0.299] = 1

sum(scores$clust.90)
sum(scores$clust.95)
write.csv(scores, "Y:\\AD\\AD_phenotyping_yunqing\\pred_on_full_data0710.csv")

##--------------------------------------------------------------------------------------------------------------
#evaluation matrices and validation by group

#re-scale the KOMAP score
score_MAIN = out.corrupt$pred_prob$pred.score
score = transform(
  score_MAIN,
  scale = `mainICDNLP + codify & NLP` / abs(ifelse(`mainICDNLP + codify & NLP` > 0, max(`mainICDNLP + codify & NLP`), min(`mainICDNLP + codify & NLP`)))
  
)

#re-scale log count of main CUI and main code
feature_score = transform(
  dat.merge,
  cui_scale= C0002395 / abs(ifelse(C0002395 > 0, max(C0002395), min(C0002395))),
  phe1 = PheCode.290.11 / abs(ifelse(PheCode.290.11 > 0, max(PheCode.290.11), min(PheCode.290.11))),
  phe2 = PheCode.290.1 / abs(ifelse(PheCode.290.1 > 0, max(PheCode.290.1), min(PheCode.290.1))),
  phe3 = PheCode.290 / abs(ifelse(PheCode.290 > 0, max(PheCode.290), min(PheCode.290)))
  )

feature_score[is.na(feature_score)] = 0


## run the KOMAP_ROC.R before the folling code

get_eval = function(y, s){
  eval = ROC(y, s)
  auc = COMP_AUC(eval$FPR, eval$TPR)
  prauc = PR_AUC(eval$PPV, eval$TPR)
  fpr = eval$FPR[which(eval$FPR<=0.1)[1]]
  tpr = eval$TPR[which(eval$FPR<=0.1)[1]] 
  ppv = eval$PPV[which(eval$FPR<=0.1)[1]] 
  npv = (1-fpr)*(1-eval$prev)/((1-fpr)*(1-eval$prev)+(1-tpr)*eval$prev)
  print(c("auc","prauc","tpr","ppv","npv"))
  print(c(auc,prauc,tpr,ppv,npv))
  #print(eval$PPV)
}


# stratified data: f=female, m=male, w=non Hispanic white, o=other
demo = read.csv("AD\\AD_phenotyping_yunqing\\demo.csv", sep="|")[,c(1,6,13,16)]
colnames(demo) = c("PatientNum", "Gender_Legal_Sex", "Race", "Ethnic_Group")
all_dat_strat = merge(dat_log_wide, demo, by="PatientNum")

dat_log_wide_f = all_dat_strat[all_dat_strat$Gender_Legal_Sex=="Female",]
dat_log_wide_m = all_dat_strat[all_dat_strat$Gender_Legal_Sex=="Male",]

dat_log_wide_w = all_dat_strat[all_dat_strat$Ethnic_Group=="Non Hispanic" & all_dat_strat$Race=="White",]
dat_log_wide_o = all_dat_strat[(all_dat_strat$Race %in% c("Black","Other","Asian","Hawaiian","American Indian",
                                                          "Native Hawaiian or Other Pacific Islander",
                                                          "American Indian or Alaska Native")
                                & all_dat_strat$Ethnic_Group=="HISPANIC")
                               | (all_dat_strat$Race %in% c("Black","Other","Asian","Hawaiian","American Indian",
                                                            "Native Hawaiian or Other Pacific Islander",
                                                            "American Indian or Alaska Native")
                                  & all_dat_strat$Ethnic_Group=="Non Hispanic"),]

unique(all_dat_strat$Race)

#label
yf = dat.merge[dat.merge$patient_num%in%dat_log_wide_f$PatientNum,]$Y
ym = dat.merge[dat.merge$patient_num%in%dat_log_wide_m$PatientNum,]$Y
yw = dat.merge[dat.merge$patient_num%in%dat_log_wide_w$PatientNum,]$Y
yo = dat.merge[dat.merge$patient_num%in%dat_log_wide_o$PatientNum,]$Y
#KOMAP scaled score
sf = score[score$patient_num%in%dat_log_wide_f$PatientNum,]$scale
sm = score[score$patient_num%in%dat_log_wide_m$PatientNum,]$scale
sw = score[score$patient_num%in%dat_log_wide_w$PatientNum,]$scale
so = score[score$patient_num%in%dat_log_wide_o$PatientNum,]$scale
#main CUI/code score
#use main CUI as an example (replace the variable name for main code)
feature_sf = feature_score[feature_score$patient_num%in%dat_log_wide_f$PatientNum,]$phe3#[c("cui_scale", "phe1", "phe2", "phe3")]
feature_sm = feature_score[feature_score$patient_num%in%dat_log_wide_m$PatientNum,]$phe3#[c("cui_scale", "phe1", "phe2", "phe3")]
feature_sw = feature_score[feature_score$patient_num%in%dat_log_wide_w$PatientNum,]$phe3#[c("cui_scale", "phe1", "phe2", "phe3")]
feature_so = feature_score[feature_score$patient_num%in%dat_log_wide_o$PatientNum,]$phe3

cui = feature_score$cui_scale
phe290.11 = feature_score$phe1
phe290.1 = feature_score$phe2
phe290 = feature_score$phe3

get_eval(dat.merge$Y, score$scale)
get_eval(dat.merge$Y, phe290.11)
get_eval(dat.merge$Y, phe290.1)
get_eval(dat.merge$Y, phe290)
get_eval(dat.merge$Y, cui)


##--------------------------------------------------------------------------------------------------------
#trained on full dataset (KOMAP seperate)
get_eval(yw, sw)
get_eval(yo, so)
get_eval(ym, sm)
get_eval(yf, sf)


## without komap
get_eval(yw, feature_sw)
get_eval(yo, feature_so)
get_eval(ym, feature_sm)
get_eval(yf, feature_sf)







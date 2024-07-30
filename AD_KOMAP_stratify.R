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


dat_f = dat_log_wide_f[,-c(105,106,107)]
dat_m = dat_log_wide_m[,-c(105,106,107)]
dat_w = dat_log_wide_w[,-c(105,106,107)]
dat_o = dat_log_wide_o[,-c(105,106,107)]


KOMAP_run = function(dat_log_wide){
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
  
  ##--------------------------------------------------------------
  
  # gold labels
  
  dat.label <- read.csv("AD\\AD_phenotyping_yunqing\\label_bin.csv")[,-1]
  
  ## KOMAP
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
  
  out.corrupt = KOMAP_corrupt(input.cov.train, input.cov.valid, is.wide = TRUE,
                              target.code = main_code, target.cui = main_CUI, nm.disease = phename,
                              nm.utl, nm.corrupt.code, nm.corrupt.cui, nm.multi = NULL, 
                              codify.feature = codify_select, nlp.feature = cui_select_weight,
                              pred = FALSE, eval.real = TRUE, eval.sim = FALSE,
                              # mu0 = NULL, mu1 = NULL, var0 = NULL, var1 = NULL, prev_Y = NULL, B = 10000,
                              dat.part = dat.merge, nm.id = 'patient_num', nm.pi = 'pi', nm.y = 'Y'
  )

  
  score = out.corrupt$pred_prob$pred.score
  score_scale = transform(
    score,
    score_scale = `mainICDNLP + codify & NLP` / abs(ifelse(`mainICDNLP + codify & NLP` > 0, max(`mainICDNLP + codify & NLP`), min(`mainICDNLP + codify & NLP`)))
  )
  
  
  eval = ROC(dat.merge$Y, score_scale$score_scale)
  print(COMP_AUC(eval$FPR, eval$TPR))
  print(PR_AUC(eval$PPV, eval$TPR))
  
  fpr = eval$FPR[which(eval$FPR<=0.1)[1]]
  tpr = eval$TPR[which(eval$FPR<=0.1)[1]] 
  ppv = eval$PPV[which(eval$FPR<=0.1)[1]] 
  npv = (1-fpr)*(1-eval$prev)/((1-fpr)*(1-eval$prev)+(1-tpr)*eval$prev)
  cut = eval$cut[which(eval$FPR<=0.1)[1]]
  print(c(tpr, ppv, npv, cut))
  
  
  
  
}

KOMAP_run(dat_f)
KOMAP_run(dat_m)
KOMAP_run(dat_w)
KOMAP_run(dat_o)

KOMAP_run(dat_log_wide)



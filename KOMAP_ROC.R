ROC <- function(true, pred, cutoffs=seq(-1,1,0.0001)) {
  TPR = FPR = PPV = rep(NA, length(cutoffs))
  prev = mean(true)
  for (i in 1:length(cutoffs)){
    c = cutoffs[i]
    TPR[i] = sum((pred >= c) * true) / sum(true)
    FPR[i] = sum((pred >= c) * (1-true)) / sum(1-true)
    PPV[i] = ifelse( FPR[i]+TPR[i] == 0, 1, TPR[i] * prev / (FPR[i] * (1-prev) + TPR[i] * prev))
  
  }

  return (list(TPR=TPR, FPR=FPR, PPV=PPV, prev=prev, cut=cutoffs))
}

COMP_AUC <- function(FPR, TPR) {
  if (length(FPR) != length(TPR)) {print("Lengths don't match")}
  sums = 0
  for (i in 2:length(FPR)) {
    sums = sums + (TPR[i-1]+TPR[i]) * abs(FPR[i-1] - FPR[i])/2
  }
  return(sums)
}

PR_AUC <- function(PPV, TPR) {
  if (length(PPV) != length(TPR)) {print("Lengths don't match")}
  sums = 0
  for (i in 2:length(TPR)) {
    sums = sums + abs(TPR[i]-TPR[i-1]) * (PPV[i-1] + PPV[i])/2
  }
  return(sums)
}

##-----------------------------------------------------------------
#rm(out.corrupt)
#load('Y:\\AD\\AD_phenotyping_yunqing\\AD_out_corrupt_datnew_female.RData')

# outcome = out.corrupt.pred$pred_prob$pred.score
# outcome$label1 = 0
# outcome$label1[outcome$`mainICD + codify`=="disease"] = 1
# outcome$label2 = 0
# outcome$label2[outcome$`mainICDNLP + codify & NLP`=="disease"] = 1
# View(outcome)
# 
# score = out.corrupt$pred_prob$pred.score
# score1 = transform(
#   score,
#   scale1 = `mainICD + codify` / abs(ifelse(`mainICD + codify` > 0, max(`mainICD + codify`), min(`mainICD + codify`)))
# )
# 
# score2 = transform(
#   score,
#   scale2 = `mainICDNLP + codify & NLP` / abs(ifelse(`mainICDNLP + codify & NLP` > 0, max(`mainICDNLP + codify & NLP`), min(`mainICDNLP + codify & NLP`)))
# )
# 
# 
# 
# dat_female = read.csv("Y:\\AD\\AD_phenotyping_yunqing\\KOMAP_stratified\\dat_female.csv")
# View(dat_female)
# 
# 
# 
# cui_score = transform(
#   dat_female,
#   test.scale= C0002395 / abs(ifelse(C0002395 > 0, max(C0002395), min(C0002395))))
# cui_score[is.na(cui_score$test.scale)] = 0
# View(cui_score)
# cui_eval = ROC(dat_female$is_AD, cui_score$test.scale)
# ppv = cui_eval$PPV[which(cui_eval$FPR<=0.1)[1]]
# 
# cod_score = transform(
#   dat_female,
#   test.scale= PheCode.290 / abs(ifelse(PheCode.290 > 0, max(PheCode.290), min(PheCode.290))))
# cod_score$test.scale[is.na(cod_score$test.scale)] = 0
# View(cod_score)
# cod_eval = ROC(dat_female$is_AD, cod_score$test.scale)
# ppv = cod_eval$PPV[which(cod_eval$FPR<=0.1)[1]]
# 
# 
# 
# 
# eval1 = ROC(dat.merge$Y, score1$scale1)
# COMP_AUC(eval1$FPR, eval1$TPR)
# PR_AUC(eval1$PPV, eval1$TPR)
# #plot(eval1$FPR, eval1$TPR)
# 
# eval2 = ROC(dat.merge$Y, score2$scale2)
# COMP_AUC(eval2$FPR, eval2$TPR)
# PR_AUC(eval2$PPV, eval2$TPR)
# plot(eval2$FPR, eval2$TPR)
# 
# PRAUC(score2$scale2, dat.merge$Y)
# 
# 
# fpr = eval1$FPR[which(eval1$FPR<=0.1)[1]]
# tpr = eval1$TPR[which(eval1$FPR<=0.1)[1]] #first tpr when fpr<=0.05
# ppv = eval1$PPV[which(eval1$FPR<=0.1)[1]] 
# npv = (1-fpr)*(1-eval2$prev)/((1-fpr)*(1-eval2$prev)+(1-tpr)*eval2$prev)
# npv
# tpr
# ppv
# 
# 
# fpr = eval2$FPR[which(eval2$FPR<=0.1)[1]]
# tpr = eval2$TPR[which(eval2$FPR<=0.1)[1]] #first tpr when fpr<=0.05
# ppv = eval2$PPV[which(eval2$FPR<=0.1)[1]] 
# npv = (1-fpr)*(1-eval2$prev)/((1-fpr)*(1-eval2$prev)+(1-tpr)*eval2$prev)
# cut1 = eval2$cut[which(eval2$FPR<=0.05)[1]]
# cut2 = eval2$cut[which(eval2$FPR<=0.1)[1]]
# cut
# npv
# tpr
# ppv
# 
# 
# ##-------------------------------------------------------
# 
# roc.test = ROC(dat.merge$Y, test2$test.scale2)
# comp = COMP_AUC(roc.test$FPR,roc.test$TPR)
# comp
# 
# prauc = COMP_AUC(roc.test$TPR,roc.test$PPV)
# prauc
# #roc(dat.merge$Y, test$test.scale)
# 
# fpr = roc.test$FPR[which(roc.test$FPR<=0.1)[1]]
# tpr = roc.test$TPR[which(roc.test$FPR<=0.1)[1]] #first tpr when fpr<=0.05
# ppv = roc.test$PPV[which(roc.test$FPR<=0.1)[1]] 
# npv = (1-fpr)*(1-roc.test$prev)/((1-fpr)*(1-roc.test$prev)+(1-tpr)*roc.test$prev)
# npv
# tpr
# ppv
# 
# fpr = roc.test$FPR[which(roc.test$FPR<=0.05)[1]]
# tpr = roc.test$TPR[which(roc.test$FPR<=0.05)[1]] #first tpr when fpr<=0.05
# ppv = roc.test$PPV[which(roc.test$FPR<=0.05)[1]] 
# npv = (1-fpr)*(1-roc.test$prev)/((1-fpr)*(1-roc.test$prev)+(1-tpr)*roc.test$prev)
# npv
# tpr
# ppv
# ##control for .05, .1
# 
# 
# 
# 

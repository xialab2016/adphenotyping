library(tidyverse)
library(glmnet)
source("utility_20221127.R")

#############################################################################################
## 1. 	Read in codifed data, 
##      Create patient-level codified/utilization and nlp data 
##      Take log(x+1) of counts and save to "AD_pt_cod_data.Rdata"
#############################################################################################
data.dir = "/your_working_directory/"
data = read.csv(paste0(data.dir, "/Codified+NLP_data_aggregated.csv"))

unq.features = unique(data$feature_id)

freq = read.csv(paste0(data.dir, "/Codified+NLP_frequency_summary.csv"))
features = c(freq[freq$freq.code > 5, "feature_id"], "AB0000AD")

f.counts.patient = function(dat, v.date, v.nm){
  dat = dat %>% rename("v.D" = matches(v.date), "v.Y" = matches(v.nm))
  dat %>% filter(v.Y %in% features) %>%
    distinct(patient_num, v.D, v.Y) %>%				## unique features each date per patient
    group_by(patient_num, v.Y) %>% 					## days of each feature per patient
    summarize(count = n()) %>% ungroup() %>%
    spread(v.Y, count) %>% replace(is.na(.), 0)		## transpose the data from long to wide
}
d.pt.codnlp = f.counts.patient(data, v.date="start_date", v.nm="feature_id")
utils = rowSums(d.pt.codnlp[,-1])
d.pt.codnlp$utils = utils
d.pt.codnlp[,-1] = log(d.pt.codnlp[,-1]+1)
# Replace special characters in column names by '.'
colnames(d.pt.codnlp) = gsub(x=colnames(d.pt.codnlp), pattern="\\:", replacement=".")
colnames(d.pt.codnlp) = gsub(x=colnames(d.pt.codnlp), pattern="\\-", replacement=".")
colnames(d.pt.codnlp) = gsub(x=colnames(d.pt.codnlp), pattern=" ", replacement="")
save(d.pt.codnlp, file="/n/data1/hsph/biostat/celehs/lab/SHARE/UPMC/AD/data_table/data_processed/Lin/AD_pt_codnlp_data_10252023_uniq.Rdata")

#############################################################################################
## 2. 	Read in patient-level data from 1, chart review labels, and demographics, 
##      Run KOMAP and get validation metrics
#############################################################################################
#------------------- Load data -------------------#
load("AD_pt_codnlp_data_10252023_uniq.Rdata")  #d.pt.codnlp
load("AD_label_data_chart_all.rda")  # d.label.new

# Read in demographic data and create variable "white" and "female"
demo = read.table("ad_demographics.dsv", 
                  header=T,
                  sep="\t", quote="")
d = demo[,c("PATIENT_STUDY_ID", "RACE_TITLE", "ETHNIC_TITLE", "GENDER_CODE")]
d$white = ifelse((d$RACE_TITLE=="WHITE") & (d$ETHNIC_TITLE != "HISPANIC OR LATINO"), 1, 0)
d$female = ifelse(d$GENDER_CODE==1, 1, 0)

#---------------- Path out and filtered date ---------------#  *CHANGE
path.out = "results_final/"
date = "2011-01-01"

#--------------- Basic Info of datasets -----------#
dim(d.pt.codnlp)  # 78286 x 23215
dim(d.label.new)  # 204 x 2
d.label = d.label.new
d.pt = d.pt.codnlp

#----------------- code names -----------------#
# Specify main code and outcome variables
main.code = "PheCode.290.11"
parent.code = c("PheCode.290.1", "PheCode.290")
outcome = "AD_probable_possible"
main.cui = "C0002395"
ac.cui = "AB0000AD"

#----------------- merge data -----------------#
d.pt.y = merge(d.pt, 
               d.label[,c("patient_num", outcome)], 
               by="patient_num", 
               all.x=T)
# d.pt.y[,"PheCode.290.1x"] = log(exp(d.pt.y[,"PheCode.290.1"]) - 
#                                   exp(d.pt.y[,"PheCode.290.11"]) + 1)

# Combine race and gender info with others
d.pt.y.group = merge(d.pt.y, 
                     d[,c("PATIENT_STUDY_ID", "white", "female")], 
                     by.x="patient_num", 
                     by.y="PATIENT_STUDY_ID",
                     all.x=T)

# Labeled data
d.labeled = d.pt.y.group[!is.na(d.pt.y.group[outcome]),]
prev = mean(d.labeled[,outcome])

#-------------- Performance of individual features ---------------#
FPR = 0.1
features = c("C0002395", "PheCode.290", "PheCode.290.1", "PheCode.290.11")
group.nm = "white" # group.nm = "female"
  
label.1 = d.labeled[d.labeled[,group.nm]==1,]; prev.1 = mean(label.1[,outcome])
label.0 = d.labeled[d.labeled[,group.nm]==0,]; prev.0 = mean(label.0[,outcome])

results.mat = results.1.mat = results.0.mat = matrix(NA, ncol=5, nrow = length(features))
colnames(results.mat) = colnames(results.1.mat) =colnames(results.0.mat) =c("AUC", "PRAUC", "TPR", "PPV", "NPV")
rownames(results.mat) = rownames(results.1.mat) = rownames(results.0.mat) = features

for (f in features){
  print(paste0("======",f,"======"))
  roc = ROC(d.labeled[,outcome],rank(d.labeled[,f])/nrow(d.labeled))
  print("------All------")
  results.mat[f,"AUC"] = COMP_AUC(roc$FPR, roc$TPR)
  results.mat[f,"PRAUC"] = COMP_PRAUC(roc$TPR, roc$PPV)
  results.mat[f,"TPR"] = TPR = roc$TPR[which.min(abs(roc$FPR-FPR))]
  results.mat[f,"PPV"] = TPR * prev / (FPR * (1-prev) + TPR * prev)
  results.mat[f,"NPV"] = (1-FPR) * (1-prev) / ((1-FPR) * (1-prev) + (1-TPR) * prev)
  
  roc.1 = ROC(label.1[,outcome],rank(label.1[,f])/nrow(label.1))
  roc.0 = ROC(label.0[,outcome],rank(label.0[,f])/nrow(label.0))
  print("------Non-hisp white------")
  results.1.mat[f,"AUC"] = COMP_AUC(roc.1$FPR, roc.1$TPR)
  results.1.mat[f,"PRAUC"] = COMP_PRAUC(roc.1$TPR, roc.1$PPV)
  results.1.mat[f,"TPR"] = TPR.1 = roc.1$TPR[which.min(abs(roc.1$FPR-FPR))]
  results.1.mat[f,"PPV"] = TPR.1 * prev.1 / (FPR * (1-prev.1) + TPR.1 * prev.1)
  results.1.mat[f,"NPV"] = (1-FPR) * (1-prev.1) / ((1-FPR) * (1-prev.1) + (1-TPR.1) * prev.1)
  
  print("------others------")
  results.0.mat[f,"AUC"] = COMP_AUC(roc.0$FPR, roc.0$TPR)
  results.0.mat[f,"PRAUC"] = COMP_PRAUC(roc.0$TPR, roc.0$PPV)
  results.0.mat[f,"TPR"] = TPR.0 = roc.0$TPR[which.min(abs(roc.0$FPR-FPR))]
  results.0.mat[f,"PPV"] = TPR.0 * prev.0 / (FPR * (1-prev.0) + TPR.0 * prev.0)
  results.0.mat[f,"NPV"] = (1-FPR) * (1-prev.0) / ((1-FPR) * (1-prev.0) + (1-TPR.0) * prev.0)
}

results.mat
results.1.mat
results.0.mat

#-------------- Performance of KOMAP ---------------#
once.cod = c(read.csv(file="ONCE_AD_cod.csv")  %>% select(Variable))$Variable
once.nlp = c(read.csv(file="ONCE_AD_nlp.csv") %>% select(cui))$cui
features.once = gsub(c(once.cod, once.nlp), pattern = "\\:", replacement = ".")  

library(KOMAP)
d.pt.komap = d.pt[,colnames(d.pt) %in% c("patient_num", features.once, "utils")]
d.pt.cov = cov(d.pt[,colnames(d.pt) %in% c(features.once, "utils")])

# Overall KOMAP
komap = KOMAP(d.pt.cov , is.wide = TRUE, main.code, main.cui, "utils",
              pred = FALSE, eval.real = FALSE, eval.sim = FALSE)

# Seperate KOMAP
d.pt.1 = d.pt.y.group[d.pt.y.group[,group.nm]==1,]
d.pt.0 = d.pt.y.group[d.pt.y.group[,group.nm]==0,]
group1.cov = cov(na.omit(d.pt.1[,colnames(d.pt.1) %in% c(features.once, "utils")]))
group0.cov = cov(na.omit(d.pt.0[,colnames(d.pt.0) %in% c(features.once, "utils")]))
komap.1 = KOMAP(group1.cov , is.wide = TRUE, main.code, main.cui, "utils",
                pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
komap.0 = KOMAP(group0.cov , is.wide = TRUE, main.code, main.cui, "utils",
                pred = FALSE, eval.real = FALSE, eval.sim = FALSE)

# KOMAP performance
pred.icdnlp = as.matrix(d.labeled[,komap$est$lst$`mainICDNLP + allfeature`$beta$feat]) %*% c(komap$est$lst$`mainICDNLP + allfeature`$beta$theta)
pred.icdnlp = rank(pred.icdnlp) / nrow(d.labeled)

## 1. Overall KOMAP
### 1a. Overall
roc = ROC(d.labeled[,outcome],pred.icdnlp)
AUC = COMP_AUC(roc$FPR, roc$TPR)
PRAUC = COMP_PRAUC(roc$TPR, roc$PPV)
TPR = roc$TPR[which.min(abs(roc$FPR-FPR))]
PPV = TPR * prev / (FPR * (1-prev) + TPR * prev)
NPV = (1-FPR) * (1-prev) / ((1-FPR) * (1-prev) + (1-TPR) * prev)

### 1b. In group 1
roc.1 = ROC(label.1[,outcome],pred.icdnlp[d.labeled[,group.nm]==1])
AUC.1 = COMP_AUC(roc.1$FPR, roc.1$TPR)
PRAUC.1 = COMP_PRAUC(roc.1$TPR, roc.1$PPV)
TPR.1 = roc.1$TPR[which.min(abs(roc.1$FPR-FPR))]
PPV.1 = TPR.1 * prev.1 / (FPR * (1-prev.1) + TPR.1 * prev.1)
NPV.1 = (1-FPR) * (1-prev.1) / ((1-FPR) * (1-prev.1) + (1-TPR.1) * prev.1)

### 1c. In group 0
roc.0 = ROC(label.0[,outcome],pred.icdnlp[d.labeled[,group.nm]==0])
AUC.0 = COMP_AUC(roc.0$FPR, roc.0$TPR)
PRAUC.0 = COMP_PRAUC(roc.0$TPR, roc.0$PPV)
TPR.0 = roc.0$TPR[which.min(abs(roc.0$FPR-FPR))]
PPV.0 = TPR.0 * prev.0 / (FPR * (1-prev.0) + TPR.0 * prev.0)
NPV.0 = (1-FPR) * (1-prev.0) / ((1-FPR) * (1-prev.0) + (1-TPR.0) * prev.0)

### Save results
komap.all.perf = cbind(c(AUC, PRAUC, TPR, PPV, NPV),
                       c(AUC.1, PRAUC.1, TPR.1, PPV.1, NPV.1),
                       c(AUC.0, PRAUC.0, TPR.0, PPV.0, NPV.0))
rownames(komap.all.perf) = c("AUC", "PRAUC", "TPR", "PPV", "NPV")
colnames(komap.all.perf) = c("Overall", "Group1", "Group0")
komap.all.perf
#write.csv(komap.all.perf, file="output/komap_all_perf.csv")

## 2.Separately trained KOMAP
### 2a. Overall
features.ind.1 = which(komap.1$est$lst$`mainICD + allfeature`$beta$feat %in% colnames(d.pt.1))
pred.1.icdnlp = as.matrix(label.1[,komap.1$est$lst$`mainICDNLP + allfeature`$beta$feat[features.ind.1]]) %*% 
  c(komap.1$est$lst$`mainICDNLP + allfeature`$beta$theta[features.ind.1])

features.ind.0 = which(komap.0$est$lst$`mainICD + allfeature`$beta$feat %in% colnames(d.pt.0))
pred.0.icdnlp = as.matrix(label.0[,komap.0$est$lst$`mainICDNLP + allfeature`$beta$feat[features.ind.0]]) %*% 
  c(komap.0$est$lst$`mainICDNLP + allfeature`$beta$theta[features.ind.0])
pred.icdnlp = rank(c(pred.1.icdnlp, pred.0.icdnlp)) / nrow(d.labeled)

roc = ROC(c(d.labeled[d.labeled[,group.nm]==1,outcome], d.labeled[d.labeled[,group.nm]==0,outcome]), pred.icdnlp)
AUC = COMP_AUC(roc$FPR, roc$TPR)
PRAUC = COMP_PRAUC(roc$TPR, roc$PPV)
TPR = roc$TPR[which.min(abs(roc$FPR-FPR))]
PPV = TPR * prev / (FPR * (1-prev) + TPR * prev)
NPV = (1-FPR) * (1-prev) / ((1-FPR) * (1-prev) + (1-TPR) * prev)

### 2b. In group 1
pred.1.icdnlp = rank(pred.1.icdnlp) / nrow(label.1)
roc.1 = ROC(label.1[,outcome],pred.1.icdnlp)
AUC.1 = COMP_AUC(roc.1$FPR, roc.1$TPR)
PRAUC.1 = COMP_PRAUC(roc.1$TPR, roc.1$PPV)
TPR.1 = roc.1$TPR[which.min(abs(roc.1$FPR-FPR))]
PPV.1 = TPR.1 * prev.1 / (FPR * (1-prev.1) + TPR.1 * prev.1)
NPV.1 = (1-FPR) * (1-prev.1) / ((1-FPR) * (1-prev.1) + (1-TPR.1) * prev.1)

### 2c. In group 0
pred.0.icdnlp = rank(pred.0.icdnlp) / nrow(label.0)
roc.0 = ROC(label.0[,outcome],pred.0.icdnlp)
AUC.0 = COMP_AUC(roc.0$FPR, roc.0$TPR)
PRAUC.0 = COMP_PRAUC(roc.0$TPR, roc.0$PPV)
TPR.0 = roc.0$TPR[which.min(abs(roc.0$FPR-FPR))]
PPV.0 = TPR.0 * prev.0 / (FPR * (1-prev.0) + TPR.0 * prev.0)
NPV.0 = (1-FPR) * (1-prev.0) / ((1-FPR) * (1-prev.0) + (1-TPR.0) * prev.0)

### Save results
komap.sep.perf = cbind(c(AUC, PRAUC, TPR, PPV, NPV),
                       c(AUC.1, PRAUC.1, TPR.1, PPV.1, NPV.1),
                       c(AUC.0, PRAUC.0, TPR.0, PPV.0, NPV.0))
rownames(komap.sep.perf) = c("AUC", "PRAUC", "TPR", "PPV", "NPV")
colnames(komap.sep.perf) = c("Overall", "Group1", "Group0")
komap.sep.perf
#write.csv(komap.sep.perf, file="output/komap_sep_perf.csv")

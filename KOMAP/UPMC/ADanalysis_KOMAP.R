library(openxlsx)
library(glmnet)
library(PRROC)
library(ranger)
library(MAP)
library(data.table)
library(stringr)
library(KOMAP)
library(caret)

source('/n/data1/hsph/biostat/celehs/lab/xix636/KOMAP/new_disease_whole_process/gen_paper_plot/sparse_model_0906/functions.R')
path.function = ""
source(paste0(path.function,"prior_function/PheNorm/FUN_PheNorm_Publish_ZH_Yuri.R"))
source(paste0(path.function,"prior_function/PheNorm/library_v2.r"))
source(paste0(path.function,"prior_function/PheNorm/library_v3.r")) 
source(paste0(path.function,'prior_function/library_v2_auc.r'))
source(paste0(path.function,'prior_function/Library_auc.R'))
source(paste0(path.function,'prior_function/function.R'))
source(paste0(path.function,'prior_function/myfunction_new.R'))
source(paste0(path.function,'anchor.R'))

source('corrupt_function/helper_corrupt.R')
#source('corrupt_function/helper.R')
source('corrupt_function/helper_check.R')
source('corrupt_function/other_functions.R')

#############################################################################################
## 1. 	Read in data, 
##      Create patient-level codified/utilization and nlp data 
##      Take log(x+1) of counts and save to "AD_pt_cod_data.Rdata"
#############################################################################################
data.dir = "/n/data1/hsph/biostat/celehs/lab/SHARE/UPMC/AD/data_table/"
data = read.csv(paste0(data.dir, "data_processed/UPMC_AD_2011_to_2021_Codified+NLP_data_aggregated_2023-03-13.csv"))

unq.features = unique(data$feature_id)
"AB0000AD" %in% unq.features

freq = read.csv(paste0(data.dir, "data_summary/UPMC_AD_2011_to_2021_Codified+NLP_frequency_summary_with_desc_2023-03-15.csv"))
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
save(d.pt.codnlp, file="AD_pt_codnlp_data_04272023_uniq.Rdata")

#############################################################################################
## New labels
#############################################################################################
#library(readxl)
lab.new = read.csv("AD_Annotation_Tracker_05162023.csv")
colnames(lab.new) = c("patient_num", "Date_diagnosis", "AD", "Date_death", "After_2011")
lab.new = lab.new[!is.na(lab.new$AD) & lab.new$AD!="",]; nrow(lab.new)
lab.new$label = ifelse(lab.new$AD %in% c("Probable", "Possible", "probable", "possible"), 1, 0)
d.label.new = lab.new[,c("patient_num", "label")]
colnames(d.label.new) = c("patient_num", "AD_probable_possible")
save(d.label.new, file="AD_label_data_chart_all.rda")


#############################################################################################
## 2.    KOMAP performance
#############################################################################################
# Load the Codfied + NLP data
load("AD_pt_codnlp_data_04272023_uniq.Rdata")
dat_lst <- vector('list', 1)
#dat_lst[[1]] <- d.pt.codnlp

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
#dict = readRDS('shiny_0711/dict.RDS')
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

# # Read in demographic data and create variable "white"
# demo = read.table("/n/data1/hsph/biostat/celehs/lab/SHARE/UPMC/AD/data_raw/ad_demographics.dsv", 
#                   header=T,
#                   sep="\t", quote="")
# d.race = demo[,c("PATIENT_STUDY_ID", "RACE_TITLE")]
# d.race$white = ifelse((d.race$RACE_TITLE=="WHITE") & (demo$ETHNIC_TITLE != "HISPANIC OR LATINO"), 1, 0)

# Combine race info with others
# dat_lst[[1]] = merge(d.pt.codnlp, 
#                     d.race[,c("PATIENT_STUDY_ID", "white")], 
#                     by.x="patient_num", 
#                     by.y="PATIENT_STUDY_ID",
#                     all.x=T)
dat_lst[[1]] = d.pt.codnlp

# gold labels
load("AD_label_data_chart_all.rda")
dat.label <- d.label.new

#### codified features selected by KESER 
once.cod = read.csv(file="ONCE_AD_cod.csv")
#codify_select <- gsub(once.cod[once.cod$high_confidence_level==1,"Variable"], pattern = "\\:", replacement = ".")
codify_select <- gsub(once.cod[,"Variable"], pattern = "\\:", replacement = ".")

#### NLP features selected
once.nlp = read.csv(file="ONCE_AD_nlp.csv")
#cui_select_weight <- gsub(once.nlp[once.nlp$high_confidence_level==1,"cui"], pattern = "\\:", replacement = ".")
cui_select_weight <- gsub(once.nlp[,"cui"], pattern = "\\:", replacement = ".")
#cui_select_major <- gsub(once.nlp[once.nlp$high_confidence_level==1,"cui"], pattern = "\\:", replacement = ".")
cui_select_major <- gsub(once.nlp[,"cui"], pattern = "\\:", replacement = ".")


# codified + nlp data
dat.all <- dat_lst[[1]]
#dat.all[is.na(dat.all)] <- 0
dat.filter_pos <- dat.all
# dat.filter_pos <- dat.filter_pos[dat.filter_pos$white == 0, ]
# dat.filter_pos <- na.omit(dat.filter_pos)
# dat.filter_pos <- dat.all[which(dat.all[[main_code]] > 0),]

# zero.rate <- apply(dat.filter_pos, 2, function(x) sum(x==0)/length(x))
# feature.keep.unsup <- colnames(dat.filter_pos)[zero.rate < 0.95]
# dat.filter_pos <- dat.filter_pos[,which(colnames(dat.filter_pos) %in% feature.keep.unsup)]
dat.filter_prev = dat.filter_pos
rownames(dat.filter_prev) = dat.filter_prev$patient_num
dat.filter_prev$patient_num = NULL
utl.index = which(colnames(dat.filter_prev) == "utils")
dat.filter_prev = dat.filter_prev[, c("utils", colnames(dat.filter_prev)[-utl.index])]
# dat.filter_prev = log(dat.filter_prev + 1)

### 06/24/2023 new: add corrupted main surrogates
dat.filter_prev_corrupt = dat.filter_prev %>% as.data.frame()
dat.filter_prev_corrupt$corrupt_mainICD = dat.filter_prev_corrupt[, main_code]
dat.filter_prev_corrupt$corrupt_mainICD[sample(1:nrow(dat.filter_prev_corrupt), round(nrow(dat.filter_prev_corrupt) * 0.2), replace = FALSE)] =
  mean(dat.filter_prev_corrupt[, main_code])
dat.filter_prev_corrupt$corrupt_mainNLP = dat.filter_prev_corrupt[, main_CUI]
dat.filter_prev_corrupt$corrupt_mainNLP[sample(1:nrow(dat.filter_prev_corrupt), round(nrow(dat.filter_prev_corrupt) * 0.2), replace = FALSE)] =
  mean(dat.filter_prev_corrupt[, main_CUI])
id.train = sample(1:nrow(dat.filter_prev_corrupt), round(nrow(dat.filter_prev_corrupt) / 3 * 2))
id.valid = setdiff(1:nrow(dat.filter_prev_corrupt), id.train)
#dat.cov.train = cov(dat.filter_prev_corrupt[id.train, ])
dat.cov.train = cov(dat.filter_prev_corrupt[id.train, colnames(dat.filter_prev_corrupt) %in% c(c(main_code, main_CUI, "corrupt_mainICD", "corrupt_mainNLP", codify_select, cui_select_weight), "utils")])
dat.cov.valid = cov(dat.filter_prev_corrupt[id.valid, colnames(dat.filter_prev_corrupt) %in% c(c(main_code, main_CUI, "corrupt_mainICD", "corrupt_mainNLP", codify_select, cui_select_weight), "utils")])
dat.cov = cov(dat.filter_prev_corrupt[,colnames(dat.filter_prev_corrupt) %in% c(c(main_code, main_CUI, "corrupt_mainICD", "corrupt_mainNLP", codify_select, cui_select_weight), "utils")])

save(dat.cov.train, dat.cov.valid, dat.cov, file="KOMAP_cov_20240515.RData")

### Run MAP
ICD <- exp(dat.filter_pos[[main_code]])-1
NLP <- exp(dat.filter_pos[[main_CUI]])-1
note <- exp(dat.filter_pos$utils)-1
mat = Matrix(data=cbind(ICD,NLP),sparse = TRUE)
note = Matrix(note,ncol=1,sparse = TRUE)
res = MAP(mat = mat, note=note)
dat.filter_pos$map <- as.vector(res$scores)

nm.X <- setdiff(colnames(dat.filter_pos), c(main_code, main_CUI, 'patient_num', 'utils'))
nm.utl <- "utils"

#### 05092024 newly added: Anchor methods
### Anchor
S_ICD <- dat.filter_pos[[main_code]]
threshold <- quantile(S_ICD, 0.8)
S <- ifelse(S_ICD >= threshold, 1, 0)

# Logistic lasso & RF:
# keser 
X <- dat.filter_pos[,union(intersect(nm.X, codify_select), c(main_CUI, nm.utl))]
dat.filter_pos$anchor.codi.lasso <- anchor(S, X, model = 'lasso')
dat.filter_pos$anchor.codi.rf <- anchor(S, X, model = 'RF')

# NLP final weight
X <- dat.filter_pos[,union(intersect(nm.X, cui_select_weight), c(main_CUI, nm.utl))]
dat.filter_pos$anchor.nlp.lasso <- anchor(S, X, model = 'lasso')
dat.filter_pos$anchor.nlp.rf <- anchor(S, X, model = 'RF')

# keser + NLP (weight)
X <- dat.filter_pos[,union(intersect(nm.X, c(cui_select_weight, codify_select)), c(main_CUI, nm.utl))]
dat.filter_pos$anchor.both.lasso <- anchor(S, X, model = 'lasso')
dat.filter_pos$anchor.both.rf <- anchor(S, X, model = 'RF')

save(dat.filter_pos, file="temp.RData")

### Run PheNorm
coef.all <- as.data.frame(matrix(0, nrow = length(unique(c(main_code, main_CUI, codify_select,
                                                           cui_select_weight))) + 1,
                                 ncol = 0))
rownames(coef.all) = unique(c(main_code, main_CUI, codify_select,
                              cui_select_weight, nm.utl))
coef.all = coef.all[rownames(coef.all) %in% colnames(dat.filter_pos),]

# main ICD code + keser (only codified features)
fit.phenorm <- PheNorm.Prob(main_code, "utils",
                            dat.filter_pos[,c(nm.X, main_code, 'utils')], 
                            nm.X = intersect(nm.X, codify_select), corrupt.rate = 0.3, train.size = 30000)
dat.filter_pos$phenorm.keser <- as.vector(fit.phenorm$probs)
coef.all$phenorm_icd = NA
var.nm = str_remove(names(fit.phenorm$betas), 'SX\\.norm\\.corrupt')
coef.all$phenorm_icd[match(var.nm, rownames(coef.all))] = fit.phenorm$betas

# main ICD + main NLP + keser + NLP (weight)
fit.phenorm <- PheNorm.Prob(c(main_code, main_CUI), "utils",
                            dat.filter_pos[,c(nm.X, main_code, main_CUI, 'utils')], 
                            nm.X = intersect(nm.X, c(codify_select, cui_select_weight)), 
                            corrupt.rate = 0.2, train.size = 20000)
dat.filter_pos$phenorm.both.weight <- as.vector(fit.phenorm$probs)
coef.all$phenorm_ICDNLP_icd = coef.all$phenorm_ICDNLP_nlp = NA
var.nm = str_remove(names(fit.phenorm$betas), 'SX\\.norm\\.corrupt')
coef.all$phenorm_ICDNLP_icd[match(var.nm, rownames(coef.all))] = fit.phenorm$betas[,1]
coef.all$phenorm_ICDNLP_nlp[match(var.nm, rownames(coef.all))] = fit.phenorm$betas[,2] 

# Komap
# library(KOMAP)
# d.pt.cov = cov(dat.filter_pos[,colnames(dat.filter_pos) %in% c(c(main_code, main_CUI, codify_select, cui_select_weight), "utils")])
# out = KOMAP(d.pt.cov, is.wide = TRUE,
#             main_code, main_CUI, nm.utl, nm.multi = NULL,
#             codify.feature=codify_select, nlp.feature=cui_select_weight,
#             pred = FALSE, eval.real = FALSE, eval.sim = FALSE)
#                 
# coef.all$komap_icd = NA
# coef.all$komap_icd[match(out$est$lst$`mainICD + codify`$beta$feat, rownames(coef.all))] = 
#   as.vector(out$est$lst$`mainICD + codify`$beta$theta)
# 
# coef.all$komap_icdnlp = NA
# coef.all$komap_icdnlp[match(out$est$lst$`mainICDNLP + codify & NLP`$beta$feat, rownames(coef.all))] = 
#   as.vector(out$est$lst$`mainICDNLP + codify & NLP`$beta$theta)
input.cov.train <- dat.cov.train
input.cov.valid <- dat.cov.valid
input.cov <- dat.cov
nm.corrupt.code = 'corrupt_mainICD'
nm.corrupt.cui = 'corrupt_mainNLP'

dat.label$Y <- dat.label$AD_probable_possible
dat.label$pi = dat.label$w = 1
dat.label <- dat.label[,c('Y', 'pi', 'patient_num')]
dat.label$patient_num = as.character(dat.label$patient_num)
dat.merge = merge(dat.filter_pos, dat.label, by = 'patient_num')
out.corrupt = KOMAP_corrupt(input.cov.train, input.cov.valid, is.wide = TRUE,
                            target.code = main_code, target.cui = main_CUI, nm.disease = phename,
                            nm.utl, nm.corrupt.code, nm.corrupt.cui, nm.multi = NULL, 
                            codify.feature = codify_select, nlp.feature = cui_select_weight,
                            pred = FALSE, eval.real = TRUE, eval.sim = FALSE,
                            # mu0 = NULL, mu1 = NULL, var0 = NULL, var1 = NULL, prev_Y = NULL, B = 10000,
                            dat.part = dat.merge, nm.id = 'patient_num', nm.pi = 'pi', nm.y = 'Y'
)
out = out.corrupt
save(out.corrupt, file="AD_out_corrupt.RData")
coef.all$komap_icd = NA
coef.all$komap_icd[match(out$est$lst$`mainICD + codify`$beta$feat, rownames(coef.all))] =
  as.vector(out$est$lst$`mainICD + codify`$beta$theta)

coef.all$komap_icdnlp = NA
coef.all$komap_icdnlp[match(out$est$lst$`mainICDNLP + codify & NLP`$beta$feat, rownames(coef.all))] =
  as.vector(out$est$lst$`mainICDNLP + codify & NLP`$beta$theta)


### Supervised method & eval at the same time
#dat.label$Y = dat.label$AD_probable_possible
#dat.label$pi = 1
dat.merge <- merge(dat.filter_pos, dat.label, by = 'patient_num')
dat.merge$komap.icd = as.matrix(dat.merge[, rownames(coef.all)[!is.na(coef.all$komap_icd)]]) %*%
  as.matrix(coef.all$komap_icd[!is.na(coef.all$komap_icd)])
dat.merge$komap.both = as.matrix(dat.merge[, rownames(coef.all)[!is.na(coef.all$komap_icdnlp)]]) %*%
  as.matrix(coef.all$komap_icdnlp[!is.na(coef.all$komap_icdnlp)])

dat.sup <- dat.merge[complete.cases(dat.merge),
                     c('Y', 'pi', intersect(nm.X, c(codify_select, cui_select_weight)), 
                       main_code, main_CUI, 'utils')]
dat.merge <- dat.merge[complete.cases(dat.merge), c('Y', 'pi', 'komap.both','komap.icd',
                                                    'phenorm.both.weight', 'phenorm.keser',
                                                    main_code, main_CUI, 'map',
                                                    'anchor.codi.lasso', 'anchor.codi.rf',
                                                    'anchor.nlp.lasso', 'anchor.nlp.rf',
                                                    'anchor.both.lasso', 'anchor.both.rf')]
# method_nm = c('komap.both','komap.icd',
#               'phenorm.both.weight', 'phenorm.keser',
#               main_code, main_CUI, 'map', 'super.both', 'super.icd')
method_nm = c('komap.both','komap.icd',
              'phenorm.both.weight', 'phenorm.keser',
              main_code, main_CUI, 'map', 
              'anchor.codi.lasso', 'anchor.codi.rf',
              'anchor.nlp.lasso', 'anchor.nlp.rf',
              'anchor.both.lasso', 'anchor.both.rf',
              'super.both', 'super.icd')

# main ICD + KESER
set.seed(1)
Y <- dat.sup$Y
X_codi <- dat.sup[,c(intersect(nm.X, codify_select), main_code, 'utils')]
X_codi <- as.matrix(X_codi)
X_both <- dat.sup[,c(intersect(nm.X, c(codify_select, cui_select_weight)), 
                     main_code, main_CUI, 'utils')]
X_both <- as.matrix(X_both)
w <- 1 / dat.sup$pi

# dat.w = dat.sup[dat.sup$white==1, ]
# Y_w <- dat.w$Y
# X_codi_w <- dat.w[,c(intersect(nm.X, codify_select), main_code, 'utils')]
# X_codi_w <- as.matrix(X_codi_w)
# X_both_w <- dat.w[,c(intersect(nm.X, c(codify_select, cui_select_weight)), 
#                      main_code, main_CUI, 'utils')]
# X_both_w <- as.matrix(X_both_w)
# w.w <- 1 / dat.w$pi
# 
# dat.b = dat.sup[dat.sup$white==0, ]
# Y_b <- dat.b$Y
# X_codi_b <- dat.b[,c(intersect(nm.X, codify_select), main_code, 'utils')]
# X_codi_b <- as.matrix(X_codi_b)
# X_both_b <- dat.b[,c(intersect(nm.X, c(codify_select, cui_select_weight)), 
#                      main_code, main_CUI, 'utils')]
# X_both_b <- as.matrix(X_both_b)
# w.b <- 1 / dat.b$pi

prev_vec <- c(prev_vec, mean(Y))
label_num_vec <- c(label_num_vec, length(Y))

num_method <- ncol(dat.merge)
p_auc_diff = auc_diff = auc_vec = auc_vec_w = auc_vec_b <- rep(0, num_method)
F_vec = F_vec_b = F_vec_w <- rep(0, num_method)
PRAUC_vec = PRAUC_vec_b = PRAUC_vec_w <- rep(0, num_method)
auc_vec_sup = F_vec_sup = PRAUC_vec_sup = auc_vec_boot_mean = c()
roc_lst = roc_lst_w = roc_lst_b <- vector('list', num_method)
pr_lst = c()
roc_pac_lst = c()
ttt = proc.time()[1]
k = 1

COMP_PRAUC <- function(TPR, PPV) {
  if (length(TPR) != length(PPV)) {print("Lengths don't match")}
  sums = 0
  for (i in 2:length(TPR)) {
    sums = sums + (PPV[i]+PPV[i-1]) * abs(TPR[i] - TPR[i-1])/2
  }
  return(sums)
}

set.seed(0)
for (t in 1:num_method){
  if(t <= num_method-2){
    print(paste0(method_nm[t],'...'))
    dat.auc = data.frame(`Y` = Y, `score` = as.vector(dat.merge[,t + 2]),
                         `w` = 1 / dat.merge$pi)
    auc_vec[t] <- AUC(Y, as.vector(dat.merge[,t + 2]), wgt=w)
    print(paste0('(full) AUC: ', auc_vec[t]))
    # ROC table
    roc_lst[[t]] <- ROC(Y, as.vector(dat.merge[,t + 2]), wgti=w, seq = seq(.01,.99,by=1e-5))[-1,]
    # Maximum F-score
    F_score_all <- 2 / (1 / roc_lst[[t]][,4] + 1 / roc_lst[[t]][,5])
    F_vec[t] <- max(F_score_all)
    
    # pr = PRROC::pr.curve(scores.class0 = dat.merge[which(Y == 1), t + 2], scores.class1 = dat.merge[which(Y == 0), t + 2],
    #                      weights.class0 = w[which(Y == 1)], weights.class1 = w[which(Y==0)], curve = TRUE,
    #                      max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
    # pr_lst = c(pr_lst, pr$auc.integral)
    # print(paste0('(package) AUPRC: ', roc$auc))
    roc = PRROC::roc.curve(scores.class0 = dat.merge[which(Y == 1), t + 2], scores.class1 = dat.merge[which(Y == 0), t + 2],
                           weights.class0 = w[which(Y == 1)], weights.class1 = w[which(Y==0)], curve = TRUE,
                           max.compute = TRUE, min.compute = TRUE, rand.compute = TRUE)
    roc_pac_lst = c(roc_pac_lst, roc$auc)
    print(paste0('(package) AUROC: ', roc$auc))
  }
  
  if(t == 1){
    # AUC
    ## baseline comparison
    auc_vec_boot = c()
    for(ii in 1:1000){
      set.seed(ii)
      iid = sample(1:nrow(dat.auc), size = nrow(dat.auc), replace = TRUE)
      dat.auc.boot = dat.auc[iid, ]
      while(sum(dat.auc.boot$Y == 1) == 0){
        iid = sample(1:nrow(dat.auc), size = nrow(dat.auc), replace = TRUE)
        dat.auc.boot = dat.auc[iid, ]
      }
      auc_vec_boot = c(auc_vec_boot, AUC(dat.auc.boot$Y,
                                         dat.auc.boot$score,
                                         wgt=dat.auc.boot$w))
    }
    auc_vec_boot_arch = auc_vec_boot
    print(paste0('(boot) Mean AUC: ', mean(auc_vec_boot)))
    auc_vec_boot_mean = c(auc_vec_boot_mean, mean(auc_vec_boot))
    auc_vec_boot_diff = rep(0, 1000)
  }else{
    if(t <= num_method-2){
      ### non-supervised
      auc_vec_boot_diff = auc_vec_boot = c()
      for(ii in 1:1000){
        set.seed(ii)
        iid = sample(1:nrow(dat.auc), size = nrow(dat.auc), replace = TRUE)
        dat.auc.boot = dat.auc[iid, ]
        while(sum(dat.auc.boot$Y == 1) == 0){
          iid = sample(1:nrow(dat.auc), size = nrow(dat.auc), replace = TRUE)
          dat.auc.boot = dat.auc[iid, ]
        }
        auc_boot = AUC(dat.auc.boot$Y,
                       dat.auc.boot$score,
                       wgt=dat.auc.boot$w)
        auc_vec_boot = c(auc_vec_boot, auc_boot)
      }
      print(paste0('(boot) Mean AUC: ', mean(auc_vec_boot)))
      auc_vec_boot_diff = auc_vec_boot_arch - auc_vec_boot
      auc_vec_boot_mean = c(auc_vec_boot_mean, mean(auc_vec_boot))
    }else{
      if(t == num_method-1){
        ### supervised
        print('Supervised, ICD+NLP...')
        ### both icd + nlp
        eval_results <- ROC.FUN.ALASSO.CV(cbind(Y, X_both), wgti=w,
                                          yes.CV=T, yes.seed=F, rep=rep, regularize=T,
                                          train_size = 50)
        model.coef = eval_results$beta
        model.coef = model.coef[str_detect(names(model.coef), '^b\\.') & !str_detect(names(model.coef), '^b\\.all')]
        names(model.coef) = str_remove(names(model.coef), '^b\\.')
        coef.all$super_icdnlp_50 = NA
        model.coef = model.coef[!names(model.coef) %in% c('utl','Intercept')]
        coef.all$super_icdnlp_50[match(names(model.coef), rownames(coef.all))] = as.vector(model.coef)
        auc_vec[t] <- eval_results$roc.cv[1]
        print(paste0('(alg) AUC: ', eval_results$roc.cv[1]))
        F_score <- max(2 / (1 / eval_results$roc.cv[2:101] + 1 / eval_results$roc.cv[102:201]))
        F_vec_sup <- c(F_vec_sup, F_score)
        F_vec[t] <- F_score
        # AUC
        ### supervised
        auc_vec_boot_diff = auc_vec_boot = c()
        for(ii in 1:200){
          if(ii %% 100 == 0) print(paste0(ii, ' th round..'))
          set.seed(ii)
          iid = sample(1:nrow(dat.auc), size = nrow(dat.auc), replace = TRUE)
          eval_results_ii <- ROC.FUN.ALASSO.CV(cbind(Y[iid], X_both[iid, ]), wgti=w[iid],
                                               yes.CV=T, yes.seed=F, rep=rep, regularize=T,
                                               train_size = 50)
          print(paste0('AUC: ', eval_results_ii$roc.cv[1]))
          auc.boot = eval_results_ii$roc.cv[1]
          auc_vec_boot = c(auc_vec_boot, auc.boot)
          auc_vec_boot_diff = c(auc_vec_boot_diff, auc_vec_boot_arch[ii] - auc_boot)
        }
        auc_vec_boot_mean = c(auc_vec_boot_mean, mean(auc_vec_boot))
        
      }else{
        ### supervised
        print('Supervised, ICD only...')
        ### ICD ONLY
        eval_results <- ROC.FUN.ALASSO.CV(cbind(Y, X_codi), wgti=w,
                                          yes.CV=T, yes.seed=F, rep=rep, regularize=T,
                                          train_size = 50)
        model.coef = eval_results$beta
        model.coef = model.coef[str_detect(names(model.coef), '^b\\.') & !str_detect(names(model.coef), '^b\\.all')]
        names(model.coef) = str_remove(names(model.coef), '^b\\.')
        coef.all$super_icd_50 = NA
        model.coef = model.coef[!names(model.coef) %in% c('utl','Intercept')]
        coef.all$super_icd_50[match(names(model.coef), rownames(coef.all))] = as.vector(model.coef)
        auc_vec[t] <- eval_results$roc.cv[1]
        print(paste0('(alg) AUC: ', eval_results$roc.cv[1]))
        F_score <- max(2 / (1 / eval_results$roc.cv[2:101] + 1 / eval_results$roc.cv[102:201]))
        F_vec_sup <- c(F_vec_sup, F_score)
        F_vec[t] <- F_score
        # AUC
        ### supervised
        auc_vec_boot_diff = auc_vec_boot = c()
        for(ii in 1:200){
          if(ii %% 100 == 0) print(paste0(ii, ' th round..'))
          set.seed(ii)
          iid = sample(1:nrow(dat.auc), size = nrow(dat.auc), replace = TRUE)
          eval_results_ii <- ROC.FUN.ALASSO.CV(cbind(Y[iid], X_codi[iid, ]), wgti=w[iid],
                                               yes.CV=T, yes.seed=F, rep=rep, regularize=T,
                                               train_size = 70)
          print(paste0('AUC: ', eval_results_ii$roc.cv[1]))
          auc.boot = eval_results_ii$roc.cv[1]
          auc_vec_boot = c(auc_vec_boot, auc.boot)
          auc_vec_boot_diff = c(auc_vec_boot_diff, auc_vec_boot_arch[ii] - auc_boot)
        }
        auc_vec_boot_mean = c(auc_vec_boot_mean, mean(auc_vec_boot))
      }
    }
  }
  auc_diff_lst[[k]] = rbind(auc_diff_lst[[k]], auc_vec_boot_diff)
  auc_all_lst[[k]] = rbind(auc_all_lst[[k]], auc_vec_boot)
  auc_diff[t] = mean(auc_vec_boot_diff); sd_auc_diff = sd(auc_vec_boot_diff)
  print(paste0('AUC diff: ', auc_diff[t]))
  p_auc_diff[t] = pnorm(mean(auc_vec_boot_diff) / sd(auc_vec_boot_diff),
                        lower.tail = FALSE)
  print(paste0('AUC diff pval: ', p_auc_diff[t]))
}

ttt1 = proc.time()[1]
print(round((ttt1 - ttt)/60, 3))

auc_table <- rbind(auc_table, auc_vec)
p_auc_diff_table = rbind(p_auc_diff_table, p_auc_diff)
auc_diff_table = rbind(auc_diff_table, auc_diff)
auc_boot_mean_table = rbind(auc_boot_mean_table, auc_vec_boot_mean)
F_table <- rbind(F_table, F_vec)

print(k)

#coef.all[is.na(coef.all)] = 0
coef.lst = c(coef.lst, list(coef.all))
colnames(auc_table) = colnames(p_auc_diff_table) =
  colnames(auc_diff_table) = colnames(auc_boot_mean_table) =  method_nm

#save(coef.lst, file = paste0('KOMAP_coef_UPMC_AD_0716.RData'))
save(auc_table, p_auc_diff_table, auc_diff_table, auc_boot_mean_table, F_table,
     coef.lst, auc_diff_lst, auc_all_lst,
     file = paste0('KOMAP_UPMC_AD_20240517.RData'))


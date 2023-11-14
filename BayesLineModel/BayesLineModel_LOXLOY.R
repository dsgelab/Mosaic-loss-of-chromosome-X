
line_dir <- "/Users/aoxliu/Documents/Project3_Finngen_mCA/Analysis_GWAS_LOX/Main_analysis/CompareLOXLOY_ShiftGWAS/CompareLOXLOY/Bayes_Line_Models/"
source(paste0(line_dir, "codes_by_Matti/line_models_functions.R") )


######################################################################################################
#                                     Read in data                                                   #
######################################################################################################

dat <- read.table(paste0(line_dir, "data/20231024_sig_LOXmeta_LOYukbJohn.BayesLineModelInput.LDQCed.tsv"), header=T, sep="\t")



######################################################################################################
#                                Run line model for 3 models                                         #
######################################################################################################

# Use 3 models
## NULL (tau = 0, slope = 1, rho = 0) NOTE: slope and rho can be set arbitrarily when tau = 0
# EFF1 only (tau = 0.15, slope = 0, rho = 0.995)
# EFF2 only (tau = 0.15, slope = Inf, rho = 0.995)
# SAME (tau = 0.15, slope = 1, rho = 0.995)

# tau = 0.15 means that we expect that 95% effects are below 2*0.15 = 0.3
# This is reasonable prior in typical GWAS but if you have larger effects, 
# adjust tau accordingly.


## the slope when effects are shared?
lm_fit <- lm(eff_lox_ivw ~ eff_loy_bolt, data=dat %>% filter(trait=="mLOX and mLOY (same variant)"))
summary(lm_fit)
# old
# (Intercept)  0.007261   0.001297   5.596  0.00139 ** 
# eff_loy_bolt 0.254346   0.004304  59.096 1.58e-09 ***
# new 
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0094795  0.0007116   13.32 1.11e-05 ***
# eff_loy_bolt 0.2573508  0.0023604  109.03 4.01e-11 ***



lm_fit <- lm(eff_lox_ivw ~ eff_loy_bolt, data=dat %>% filter(trait %in% c("mLOX and mLOY (same variant)", "mLOX and mLOY (r2>0.6)") ))
summary(lm_fit)
# old 
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.010364   0.002464   4.206  0.00181 ** 
# eff_loy_bolt 0.249937   0.009865  25.335  2.1e-10 ***
# new 
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.011833   0.002061   5.742 0.000187 ***
# eff_loy_bolt 0.253801   0.008250  30.765 3.09e-11 ***


dat %>% filter(trait=="mLOX and mLOY (same variant)") %>% 
      select(eff_loy_bolt, eff_lox_ivw) %>% 
      mutate(sl=eff_lox_ivw/eff_loy_bolt) %>% summary()  # old: Median :0.2758, Mean   :0.2941; new: Median :0.3056, Mean   :0.3110 


dat %>% filter(trait %in% c("mLOX and mLOY (same variant)", "mLOX and mLOY (r2>0.6)") ) %>% 
      select(eff_loy_bolt, eff_lox_ivw) %>% 
      mutate(sl=eff_lox_ivw/eff_loy_bolt) %>% summary()  # old, Median :0.2842, Mean   :0.3753; new: Median :0.3385, Mean   :0.3880  
freq_loy <- 41791/(205011)   # 0.2038476
freq_lox <- 86093/(818431)   # 0.1051927


taus=c(0.15, 0.15, 0.15)    # prior SD for the larger effect 
# slopes=c(0, Inf, 0.35, -0.2)    # slopes for models
# slopes=c(0, Inf, 0.35)    # slopes for models
# slopes=c(0, Inf, 0.2941)    # slopes for models
slopes=c(0, Inf, 0.30)    # slopes for models
# slopes=c(0, Inf, 1)    # slopes for models
rhos=c(0.995, 0.995, 0.995) # rhos for models
model.names=c("mLOY", "mLOX", "Same")
K=length(taus)   #number of models
model.cols = c("limegreen", "red", "orange")

#Visualize the models:
visualize.line.models(scales=taus, slopes=slopes, cors=rhos, 
                      model.names=model.names, 
                      model.cols=model.cols,
                      legend.position = "bottomright")
#NOTE: set legend.position = NULL to remove legend.


dat_line <- dat %>% filter(!is.na(eff_loy_bolt) & !is.na(eff_lox_ivw)) 
dim(dat_line)   # 199  22


#Assume that these two effects are from independent data sets:
r.lkhood = 0 # 0 means no correlation between estimators of EFF1 and EFF2
#If data come from overlapping samples, then this must be set accordingly.

#Estimate the model probabilities for each data point separately using
#uniform prior on models (here prior is 25% on each of the four models).
separate.prob = line.models(X = dat_line[,c("eff_loy_bolt","eff_lox_ivw")], SE = dat_line[,c("se_loy_bolt","se_lox_ivw")], 
                            scales=taus, slopes=slopes, cors=rhos, 
                            model.priors = rep(1/K,K), 
                            model.names = model.names, r.lkhood = r.lkhood)

res_separate <- cbind(dat_line, separate.prob) %>% 
      rename(mLOY.separate="mLOY", mLOX.separate="mLOX", Same.separate="Same") %>% 
      mutate(grp_95.separate=ifelse(mLOY.separate>=0.95, "mLOY specific", NA),
             grp_95.separate=ifelse(mLOX.separate>=0.95, "mLOX specific", grp_95.separate),
             grp_95.separate=ifelse(Same.separate>=0.95, "shared by mLOX and mLOY", grp_95.separate) ) 
dim(res_separate)   # 199 26 

# Estimate the model probabilities for all data point jointly 
# together with the proportion of data points coming from each model.
dat_line_uniq <- dat_line %>% distinct(eff_loy_bolt, eff_lox_ivw, se_loy_bolt, se_lox_ivw, trait, SNP)
dim(dat_line_uniq)   # 192  6 
joint.prob = line.models.with.proportions(X = dat_line_uniq[,c("eff_loy_bolt","eff_lox_ivw")], SE = dat_line_uniq[,c("se_loy_bolt","se_lox_ivw")],
                                          scales = taus, slopes = slopes, cors = rhos,
                                          model.names = model.names, r.lkhood = r.lkhood,
                                          n.iter = 10000, n.burnin = 1000)
head(joint.prob$groups)   #probabilities of each data.point (increas n.iter to get more accuracy)
joint.prob$params   # posterior of proportions in each group


res_join <- cbind(dat_line_uniq, joint.prob$groups) %>% 
      rename(mLOY.joint="mLOY", mLOX.joint="mLOX", Same.joint="Same") %>% 
      mutate(grp_95.joint=ifelse(mLOY.joint>=0.95, "mLOY specific", NA),
             grp_95.joint=ifelse(mLOX.joint>=0.95, "mLOX specific", grp_95.joint),
             grp_95.joint=ifelse(Same.joint>=0.95, "shared by mLOX and mLOY", grp_95.joint) ) %>% 
      select(SNP, trait, mLOY.joint, mLOX.joint, Same.joint, grp_95.joint)


## format the outputs ---------------
dat_res <- inner_join(res_separate, res_join, by=c("SNP","trait"))
dat_res <- dat_res %>% select(Gene, SNP, CHR, BP, REF, ALT, trait, 
                   eff_lox_ivw, se_lox_ivw, pval_lox_ivw, pval_lox_zscore, 
                   eff_loy_bolt, se_loy_bolt, pval_loy_bolt, 
                   mLOY.separate, mLOX.separate, Same.separate, grp_95.separate, 
                   mLOY.joint, mLOX.joint, Same.joint, grp_95.joint)

dat_res <- dat_res %>% mutate(mLOY.separate=ifelse(round(mLOY.separate,3)<0.001, formatC(mLOY.separate, format="e", digits=2), round(mLOY.separate,3)), 
                              mLOX.separate=ifelse(round(mLOX.separate,3)<0.001, formatC(mLOX.separate, format="e", digits=2), round(mLOX.separate,3)),
                              Same.separate=ifelse(round(Same.separate,3)<0.001, formatC(Same.separate, format="e", digits=2), round(Same.separate,3)),
                              
                              mLOY.joint=ifelse(round(mLOY.joint,3)<0.001, formatC(mLOY.joint, format="e", digits=2), round(mLOY.joint,3)), 
                              mLOX.joint=ifelse(round(mLOX.joint,3)<0.001, formatC(mLOX.joint, format="e", digits=2), round(mLOX.joint,3)),
                              Same.joint=ifelse(round(Same.joint,3)<0.001, formatC(Same.joint, format="e", digits=2), round(Same.joint,3)))

dat_miss <- dat %>% filter(is.na(eff_loy_bolt) | is.na(eff_lox_ivw)) %>% 
      mutate(mLOY.separate=NA, mLOX.separate=NA, Same.separate=NA, grp_95.separate=NA, 
             mLOY.joint=NA, mLOX.joint=NA, Same.joint=NA, grp_95.joint=NA) %>% 
      select(Gene, SNP, CHR, BP, REF, ALT, trait, 
                   eff_lox_ivw, se_lox_ivw, pval_lox_ivw, pval_lox_zscore, 
                   eff_loy_bolt, se_loy_bolt, pval_loy_bolt, 
                   mLOY.separate, mLOX.separate, Same.separate, grp_95.separate, 
                   mLOY.joint, mLOX.joint, Same.joint, grp_95.joint)
dim(dat_miss)   # 13 22

dat_output <- rbind(dat_res, dat_miss) %>% arrange(trait, CHR, BP)


## Add information about extended MHC ---------------
dat_output <- dat_output %>% mutate(MHC_extended=ifelse(CHR==6 & BP>=25726063 & BP<=33400644, "Yes", "No"))  


## Add P for effect size difference from a two-sided t-test  ---------------
dat_output <- dat_output %>% 
      mutate(beta_diff=eff_lox_ivw-eff_loy_bolt, 
             var_diff=se_lox_ivw^2+se_lox_ivw^2,
             p_diff_2sidettest=2*pnorm(abs(beta_diff/sqrt(var_diff)),lower.tail=F), 
             p_diff_2sidettest=ifelse(round(p_diff_2sidettest,2)<0.01, formatC(p_diff_2sidettest, format="e", digits=2), round(p_diff_2sidettest,2))) %>% 
      select(-beta_diff, -var_diff)
# write.table(dat_output, "FormattedTables/TableS12_mLOX_mLOY_AssignGroups.BayesLineModel.LDQCed.slope030.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)
write.table(dat_output, "FormattedTables/TableS12_mLOX_mLOY_AssignGroups.BayesLineModel.LDQCed.slope1.tsv", append=F, quote=F, sep="\t", row.names=F, col.names=T)



dat_res %>% group_by(grp_95.separate) %>% count()
#   grp_95.separate             n
# 1 mLOX specific              36 -> 34
# 2 mLOY specific              24 -> 25
# 3 shared by mLOX and mLOY    23 -> 21
# 4 NA                        113 -> 119


dat_res %>% group_by(grp_95.joint) %>% count()
#   grp_95.joint                n
# 1 mLOX specific              36 -> 31
# 2 mLOY specific              28 -> 35
# 3 shared by mLOX and mLOY    19 -> 19
# 4 NA                        113 -> 114



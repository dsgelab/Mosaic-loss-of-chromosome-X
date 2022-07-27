### For 49 independent loci identified by mLOX GWAS meta-analysis, 
### we investigated whether the effects were the same for cell fraction below 5% and above 5%.


##################################################################################################
#        Run statified GWAS in FinnGen for cell fraction below 5% and above 5%                   #
##################################################################################################

# the same control group was used for two statified analyses.



##################################################################################################
#                 Information of 49 independent loci from GWAS meta-analysis                     #
##################################################################################################

setwd("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_GWAS_LOX/Manuscript/Plots/figures")
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
'%!in%' <- function(x,y)!('%in%'(x,y))


## Read in the liftover list of 49 independent loci -----------
liftover_lst <- read.table("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_GWAS_LOX/Main_analysis/CompareLOXLOY_ShiftGWAS/CompareLOXLOY/Bayes_Line_Models/data/G2G_LOX_MA_GeneScores_Signals_31032022_LOX_Signals.tsv", sep="\t", header=T)
dim(liftover_lst)   # 49 7

lift38 <- read.table("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_GWAS_LOX/Main_analysis/CompareLOXLOY_ShiftGWAS/CompareLOXLOY/Bayes_Line_Models/data/hglft_genome_279aa_df5d00.bed", sep="\t", header=F)
dim(lift38)   # 49 1

sig_LOX <- cbind(liftover_lst, lift38) %>% 
                separate(col=V1, into=c("left", "postion_hg38"), sep="-") %>% 
                mutate(Variant=paste0(CHR,":",postion_hg38,":",ALLELE1,":",ALLELE0), trait="LOX") %>% 
                rename(rsid="SNP") %>% 
                select(Variant, rsid, Top_gene, trait) 

ivw_signals <- read.table("/Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_GWAS_LOX/Main_analysis/CompareLOXLOY_ShiftGWAS/CompareLOXLOY/Bayes_Line_Models/data/G2G_LOX_MA_GeneScores_Signals_31032022_LOX_Zweighted_Signals.tsv", header=T)
ivw_signals <- ivw_signals %>% rename(Top_gene_forcheck="Top_gene")
dim(ivw_signals)   # 49 11

sig_LOX_lst <- cbind(sig_LOX, ivw_signals) %>% filter(Top_gene_forcheck==Top_gene & SNP==rsid) %>% 
      select(CHR, BP, ALLELE1, ALLELE0, Variant, rsid, Top_gene, A1FREQ, BETA, SE, P, Top_gene_summary) %>% 
      arrange(desc(CHR, BP)) 

write.table(sig_LOX_lst %>% select(-CHR, -BP, -ALLELE1, -ALLELE0), "../input/G2G_LOX_GeneScores_49Signals_Signals.hg38.tsv", append=F, quote=F, sep="\t", row.names=F)
# system("gsutil cp ../input/G2G_LOX_GeneScores_49Signals_Signals.hg38.tsv  gs://dsge-aoxing/mocha/LOX_GWAS/Meta_mLOX/EachBiobank/")
# system("gsutil cp  gs://dsge-aoxing/mocha/LOX_GWAS/Meta_mLOX/EachBiobank/mLOX_FinnGenr9_AgeStrata.49loci.20220527.tsv  /Users/aoxliu/Documents/Project2_Finngen_mCA/Analysis_GWAS_LOX/Manuscript/Plots/input")





##############################################################################################
#         Bar plot with P-value from pair-wise comparisons of age-stratified analysis        #
##############################################################################################

## Read in and format age-stratified results -----------
stats_age <- read.table("../input/mLOX_FinnGenr9_AgeStrata.49loci.20220527.tsv", header=T)
dim(stats_age)    # 215  18

stats_age %>% group_by(chrom, pos, ref, alt) %>% count() %>% data.frame()  # 43
stats_age %>% inner_join(liftover_lst, by=c("chrom"="CHR", "pos"="BP")) %>% dim()
stats_age %>% mutate(Variant=paste0(chrom,":",pos,":",ref,":",alt)) %>% filter(Variant %in% c(sig_LOX$Variant,"19:54609641:C:T")) %>% dim()  # 200+5
stats_age %>% mutate(Variant=paste0(chrom,":",pos,":",ref,":",alt)) %>% filter(Variant %!in% sig_LOX$Variant) %>% dim()  # 15
stats_age <- stats_age %>% mutate(Variant=paste0(chrom,":",pos,":",ref,":",alt))
dim(stats_age)

n_variant <- nrow(stats_age)/5
stats_compare <- matrix(NA, nrow=n_variant*10, ncol=11)
colnames(stats_compare) <- c("Variant", "age_1", "age_2", "beta_age_1", "beta_age_2", "se_age_1", "se_age_2",  "p_age_1", "p_age_2", "p_diff", "sig_diff")

stats_compare[,"Variant"] <- rep(unique(stats_age$Variant), each=10)
age_grp1 <- c("00_40", "00_40", "00_40", "00_40", 
              "40_50", "40_50", "40_50", 
              "50_60", "50_60", 
              "60_70")

age_grp2 <- c("40_50", "50_60", "60_70", "70_105",
              "50_60", "60_70", "70_105",
              "60_70", "70_105",
              "70_105")

stats_compare[,"Variant"] <- rep(unique(stats_age$Variant), each=10)
stats_compare[,"age_1"] <- age_grp1
stats_compare[,"age_2"] <- age_grp2
stats_compare <- data.frame(stats_compare)

for (i in 1:nrow(stats_compare)) {
	print(i)
	variant <- stats_compare[i,"Variant"]
	age1 <- stats_compare[i,"age_1"]
	age2 <- stats_compare[i,"age_2"]
	
	d_age1 <- stats_age %>% filter(Variant==variant & age==age1)
	d_age2 <- stats_age %>% filter(Variant==variant & age==age2)
	
	stats_compare[i, "beta_age_1"] <- d_age1[1,"beta"]
	stats_compare[i, "beta_age_2"] <- d_age2[1,"beta"]
	
	stats_compare[i, "se_age_1"] <- d_age1[1,"sebeta"]
	stats_compare[i, "se_age_2"] <- d_age2[1,"sebeta"]
	
	stats_compare[i, "p_age_1"] <- d_age1[1,"pval"]
	stats_compare[i, "p_age_2"] <- d_age2[1,"pval"]
	
	as.numeric(as.character(stats_compare[i, "beta_age_1"]))
	beta_diff <- as.numeric(as.character(stats_compare[i, "beta_age_1"])) -  as.numeric(as.character(stats_compare[i, "beta_age_2"]))
	var_diff <- as.numeric(as.character(stats_compare[i, "se_age_1"]))^2 + as.numeric(as.character(stats_compare[i, "se_age_2"]))^2
	
	zp <- beta_diff/sqrt(var_diff)
	P_val <- 2*pnorm(abs(zp),lower.tail=F)
	sig <- ifelse(P_val<0.05,"T","F")
	stats_compare[i, "p_diff"] <- P_val
	stats_compare[i, "sig_diff"] <- sig
}

stats_compare <- stats_compare %>% inner_join(sig_LOX_lst[,c("Variant", "rsid", "Top_gene", "CHR", "BP")], by="Variant") %>%
      mutate(beta_age_1=round(as.numeric(as.character(beta_age_1)), 3), 
             beta_age_2=round(as.numeric(as.character(beta_age_2)), 3), 
             se_age_1=round(as.numeric(as.character(se_age_1)), 3),
             se_age_2=round(as.numeric(as.character(se_age_2)), 3),
             
             age_1=sub("_","-",age_1), 
             age_1=sub("00-40","<40",age_1), 
             age_1=sub("70-105",">70",age_1), 
             age_2=sub("_","-",age_2), 
             age_2=sub("00-40","<40",age_2), 
             age_2=sub("70-105",">70",age_2), 
            
             p_age_1=as.numeric(as.character(p_age_1)),
             p_age_2=as.numeric(as.character(p_age_2)),
             p_diff=as.numeric(as.character(p_diff)),
             p_age_1=ifelse(round(p_age_1,2)<0.01, formatC(p_age_1, format="e", digits=2), round(p_age_1,2)),
             p_age_2=ifelse(round(p_age_2,2)<0.01, formatC(p_age_2, format="e", digits=2), round(p_age_2,2)),
             p_diff=ifelse(round(p_diff,2)<0.01, formatC(p_diff, format="e", digits=2), round(p_diff,2))) %>% 
      arrange(desc(-CHR), desc(-BP)) %>% 
      select(Variant, Top_gene, 
             age_1, beta_age_1, se_age_1, p_age_1, 
             age_2, beta_age_2, se_age_2, p_age_2, p_diff)
write.table(stats_compare, "../input/G2G_LOX_GeneScores_49Signals_Signals.hg38.tsv", append=F, quote=F, sep="\t", row.names=F)


### function to plot pair-wise comparions -------------------
plot_onset <- function(SNP){
	GENE <- unlist(sig_LOX_lst %>% filter(Variant==SNP) %>% select(Top_gene))
	
	## check the beta from meta-analysis, if is negative, change it to the effect of risk allele --------------
	# risk_dir <- sign(unlist(sig_LOX_lst %>% filter(Variant==SNP) %>% select(BETA)))
	stats_age_risk <- stats_age %>% filter(Variant==SNP) %>% 
	            mutate(sig=ifelse(pval<0.05, "T", "F")) %>% 
	            mutate(age=sub("_","-",age), 
	                   age=sub("00-40","<40",age), 
	                   age=sub("70-105",">70",age),
	                   age=factor(age, levels=c("<40", "40-50", "50-60", "60-70", ">70")) ) %>% 
	            arrange(-desc(age))
	
	BETA <- unlist(sig_LOX_lst %>% filter(Variant==SNP) %>% select(BETA))
	if (BETA>0) {
		stats_age_risk <- stats_age_risk %>% mutate(beta_risk=-beta)
	} else {
		stats_age_risk <- stats_age_risk %>% mutate(beta_risk=beta)
	}

	pair_test <- stats_compare %>% filter(Variant==SNP) %>% 
	        mutate(group1=sub("_","-",age_1), group1=sub("00-40","<40",group1), group1=sub("70-105",">70",group1), 
	               group2=sub("_","-",age_2), group2=sub("00-40","<40",group2), group2=sub("70-105",">70",group2) ) %>% 
	        mutate(p_diff=as.numeric(as.character(p_diff)), 
	               p_diff_label=ifelse(p_diff<0.05,"*", ""), 
	               p_diff_label=ifelse(p_diff<0.05/10,"**", p_diff_label), 
	               p_diff_label=ifelse(p_diff<0.05/100,"***", p_diff_label)) %>% 
	        filter(p_diff_label!="")
	
	beta_max <- max(stats_age_risk$beta_risk + stats_age_risk$sebeta)
	if (nrow(pair_test)==1) {pair_test[,"y.position"] <- beta_max+0.05 }
	if (nrow(pair_test)==2) {pair_test[,"y.position"] <- c(beta_max+0.05, beta_max+0.1) }
	if (nrow(pair_test)==3) {pair_test[,"y.position"] <- c(beta_max+0.05, beta_max+0.1, beta_max+0.15) }
	if (nrow(pair_test)==4) {pair_test[,"y.position"] <- c(beta_max+0.05, beta_max+0.1, beta_max+0.15, beta_max+0.2) }
	
	tit <- paste0(GENE, " (", SNP, ")")
	
	p_onset_effect_test <- ggplot(stats_age_risk) + 
		geom_bar(aes(x=age, y=beta_risk, fill=sig), size=0.5, width=0.5, stat="identity", alpha=0.6) + 
		geom_errorbar(aes(x=age, ymin=beta_risk-sebeta, ymax=beta_risk+sebeta, color=sig), width=0.1, size=0.5, linetype=1) + 
		geom_hline(yintercept=0, col="grey", size=0.5) + 
		labs(y="log(OR) of LOX", x="Age", title=tit) +
		scale_color_manual(breaks=c("T", "F"), values=c("#E41A1C", "black")) +
		scale_fill_manual(breaks=c("T", "F"), values=c("#E41A1C", "black")) +
		scale_size_manual(breaks=c("T", "F"), values=c(2,1)) +
		theme_classic() +
		theme(axis.text=element_text(size=9, face="bold"), axis.title=element_text(size=10, face="bold"), plot.title=element_text(size=12, face="bold")) + 
		theme(strip.background = element_blank(), strip.text=element_text(size=9, face="bold")) + 
		theme(legend.position="none")
	
	if (nrow(pair_test)>=1) {
		p_onset_effect_test <- p_onset_effect_test + stat_pvalue_manual(
			pair_test,
			y.position="y.position",
			xmin="group1",
			xmax="group2",
			label="p_diff_label",
			size=8,
			bracket.size=0.4
		)
	}
	return(p_onset_effect_test)
	# ggsave(paste0("figures/SuppFigure2_",SNP,".",GENE,".tiff"), p_onset_effect_test, width=5, height=5)
}


pdf("figures/SuppFigure2_AgeOnsetStratifiedEffect.pdf")
for (variant in unique(stats_compare$Variant)){
	print(variant)
	print(plot_onset(variant))
}
dev.off()



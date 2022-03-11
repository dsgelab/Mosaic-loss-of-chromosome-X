
# Meta-analysis for GWAS 

1. GWAS, https://github.com/freeseek/mocha/blob/master/wdl/assoc.wdl (applying regenie to MoChA output) is suggested to be used. 

2. QC for imputation info (>0.6) and MAF (>0.001) for each cohort.

3. liftover to 38 build (https://github.com/FINNGEN/META_ANALYSIS/blob/master/wdl/lift.wdl)

4. munge (harmonize with GnomAD).
   https://github.com/FINNGEN/META_ANALYSIS/blob/master/wdl/munge_wo_lift.wdl

5. Meta-analysis 
   Regarding the GWAS outcome, all cohorts defined mLOX as binary outcome with MoChA, except for UKB. For UKB, two measures are used, one as binary (with MoChA), the other as continuous by combing 3 ways of mLOX calling.
   IVW for those used logsitic regression and weighted z-score for aggregating summary stats 
   Regarding the ancestry, Biobank Japan (BBJ) is the only cohort which is not from European ancestry.
   Therefore, we performed 4 meta-analyses considering both outcome scale and ancestry: 
   (1) IVW for biobanks with European ancestry (LOX in UKB as binary outcome). 
   (2) IVW for all biobanks including BBJ (LOX in UKB as binary outcome).
   (3) Weighted z-score for biobanks with European ancestry (LOX in UKB as continuous outcome). 
   (4) Weighted z-score for biobanks all biobanks including BBJ(LOX in UKB as continuous outcome). 
   For (1) and (2) which applied IVW method, meta-analysis was done using https://github.com/FINNGEN/META_ANALYSIS/blob/master/wdl/meta.wdl.
   For (3) and (4) which applied weighted z-score, we used METAL.
   
6. Add N of cases, N of controls, total effective sample size, and effect allele frequency.

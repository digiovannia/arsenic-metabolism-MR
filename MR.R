# Executes IVW and Maximum-likelihood Mendelian Randomization analysis for
# each CpG, storing the following results:
# (1) sig: flag of '***' for estimates of P < 0.05, else '...'
# (2) [predictor]: indicates the predictor name
# (3) cpg: identifier of the CpG site in the outcome
# (4) effect: effect estimate for the primary predictor
# (5) se: standard error
# (6) CI_lower: lower bound for 95% CI
# (7) CI_upper: upper bound for 95% CI
# (8) pval: P-value for estimate for the primary predictor

library(ggplot2)
library(MendelianRandomization)

beta_DMA=c(3.59526,2.21788,5.08819)
se_DMA=c(0.47,0.34,0.51)
snp=c('rs9527','rs11191527','rs61735836')
coh <- c('heals_epic',"heals_meta",'best_meta','meta')

for(j in coh){
  medfile <- read.table(sprintf("RegFiles/snp_%s_bn_lowef.txt", j),
                        header=TRUE,as.is=TRUE)
  write('sig cpg effect se CI_lower CI_upper pval',
        file=sprintf("RegFiles/MR_%s_bn_lowef.txt", j),append=TRUE)
  write('sig cpg effect se CI_lower CI_upper pval',
        file=sprintf("RegFiles/MR_maxlik_%s_bn_lowef.txt", j),append=TRUE)
  for(c in unique(medfile$cpg)){
    subfile <- medfile[medfile$cpg == c,]
    beta_CpG <- subfile$effect
    se_CpG <- subfile$se
    mri_CpG <- mr_ivw(mr_input(bx=beta_DMA, bxse=se_DMA, by=beta_CpG,
                               byse=se_CpG, snp=snp))
    mrml_CpG <- mr_maxlik(mr_input(bx=beta_DMA, bxse=se_DMA, by=beta_CpG,
                                   byse=se_CpG, snp=snp))
    sigi <- if(mri_CpG$Pvalue <= 0.05) {"***"} else {"..."} 
    write(paste(sigi,c,mri_CpG$Estimate,mri_CpG$StdError,mri_CpG$CILower,
                mri_CpG$CIUpper,mri_CpG$Pvalue),
          file=sprintf("RegFiles/MR_%s_bn_lowef.txt", j),append=TRUE)
    sigi <- if(mrml_CpG$Pvalue <= 0.05) {"***"} else {"..."} 
    write(paste(sigi,c,mrml_CpG$Estimate,mrml_CpG$StdError,mrml_CpG$CILower,
                mrml_CpG$CIUpper,mrml_CpG$Pvalue),
          file=sprintf("RegFiles/MR_maxlik_%s_bn_lowef.txt", j),append=TRUE)
  }
}
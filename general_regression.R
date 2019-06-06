# Executes linear regressions for HEALS and BEST cohort phenotypes and SNPs,
# including CpG methylation, and stores results in a text file that includes:
# (1) sig: flag of '***' for estimates of P < 0.05, else '...'
# (2) [predictor]: indicates the predictor name
# (3) cpg: identifier of the CpG site in the outcome
# (4) effect: effect estimate for the primary predictor
# (5) se: standard error
# (6) CI_lower: lower bound for 95% CI
# (7) CI_upper: upper bound for 95% CI
# (8) pval: P-value for estimate for the primary predictor
# The combined SNP analysis file is included for convenience in MR analyses.

library(stringi)

snplist = c("X10.104623578","X10.104795134","exm1580829")
default_covars = c("Age","Sex","cig_curt","cig_former","bmi","Sample_Plate")
predictors = c(c("DMA_pct", "lnDMA", "snp", "recoded_bin_score", "wt_score"),
               snplist)
outcome = "cpg"

gen_lm = function(data,x,y,j,covar=default_covars){
  # Generates the linear model object corresponding to the desired analysis.
  if(j == "meta"){                        
    #combined meta-analysis cohorts must include adjustment for cohort
    covar = c(covar, "cohort")
  }
  if(x == "X10.104623578" | x == "bin_score" | x == "recoded_bin_score"){
    #accounts for LD
    covar = c(covar, "X10.104795134")
  }
  if(x == "X10.104795134"){                                                   
    #accounts for LD
    covar = c(covar, "X10.104623578")
  }
  if(x == "DMA_pct" | x == "lnDMA"){            
    #BMI and education covariate recommended for DMA analyses
    covar = c(covar, "bmi", "Education_year", "lnWArsenic")
  }
  if(x == 'lnDMA'){                           
    #adjust for arsenic exposure confounding in raw DMA measure
    covar = c(covar, 'lnWArsenic', 'lnUrineAs', 'lnUrineCreat')
  }
  covs = paste(covar,collapse=" + ")
  formstring = sprintf("%s ~ %s + %s",y,x,covs)
  return(lm(formula(formstring),data=data))
}

generate_regression = function(data_file,m,y,j,outcome,predictor){
  model = gen_lm(data_file, m, y, j)
  pval = summary(model)$coef[m,4]
  effect = summary(model)$coef[m,1]
  se = summary(model)$coef[m,2]
  CI_lower = confint(model, m, level = 0.95)[1]
  CI_upper = confint(model, m, level = 0.95)[2]
  sig = if(pval <= 0.05) {"***"} else {"..."} 
  write(paste(sig,m,y,effect,se,CI_lower,CI_upper,pval),
        file=sprintf("RegFiles/%s_%s_bn_lowef.txt",
                     predictor,j),append=TRUE)
}

for(predictor in predictors){
  if(predictor == 'DMA_pct' | predictor == 'lnDMA'){
    coh = c("heals_epic","heals_meta")
  } else {
    coh = c("heals_epic","heals_meta","best_meta","meta")
  }
  for(j in coh) { #Separate regressions for the four cohort-array combinations
    data_file = read.table(sprintf("lowef_med_%s_bn.txt",
                                   j),header=TRUE,as.is=TRUE)
    cpglist = names(data_file)[substr(names(data_file),1,2) == "cg"]
    write(paste("sig",predictor,outcome,"effect","se","CI_lower","CI_upper",
                "pval"),
          file=sprintf("RegFiles/%s_%s_bn_lowef.txt",
                       predictor,j),append=TRUE)
    if(predictor == 'snp'){
      for(m in snplist){
        for(y in cpglist){
          generate_regression(data_file,m,y,j,outcome,predictor)
        }
      }
    } else {
      for(y in cpglist){
        generate_regression(data_file,predictor,y,j,outcome,predictor)
      }
    }
  }
}
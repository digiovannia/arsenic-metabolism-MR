# Computes summary statistics for the HEALS and BEST study populations

library(scales)

heals = read.table(sprintf("lowef_med_%s_bn.txt",
                            "heals_meta"),header=TRUE,as.is=TRUE)
healssnp = read.table(sprintf("med_%s_bn.txt",
                               "heals_meta"),header=TRUE,as.is=TRUE)
best = read.table(sprintf("lowef_med_%s_bn.txt",
                           "best_meta"),header=TRUE,as.is=TRUE)
bestsnp = read.table(sprintf("med_%s_bn.txt",
                              "best_meta"),header=TRUE,as.is=TRUE)

for(i in 2:4){
  names(healssnp)[i] = paste0('snp', i-1)
  names(bestsnp)[i] = paste0('snp', i-1)
}

healssnp = healssnp[,1:4]
bestsnp = bestsnp[,1:4]

features = c("id","X10.104623578","X10.104795134","exm1580829","wt_score",
             "recoded_bin_score","Sex","Age","bmi","cig_curt","cig_ever",
             "cig_former","Education_year","lnDMA","DMA_pct","lnUrineCreat",
             "lnUrineAs","WArsenic","UrAsgmCr")
best_features = c("id","X10.104623578","X10.104795134","exm1580829",
                  "wt_score","recoded_bin_score","Sex","Age","bmi",
                  "cig_curt","cig_ever","cig_former","Education_year",
                  "lnUrineCreat","lnUrineAs",'UrAsgmCr')

heals = heals[,features]
best = best[,best_features]

drops = c("exm1580829","wt_score",'recoded_bin_score')

heals_excl_snp3 = heals[,!(names(heals) %in% drops)]
best_excl_snp3 = best[,!(names(best) %in% drops)]

heals = merge(heals, healssnp, by="id")
best = merge(best, bestsnp, by="id")

heals = heals[complete.cases(heals_excl_snp3),]
best = best[complete.cases(best_excl_snp3),]

Male_heals = as.numeric(heals$Sex == 1)
Female_heals = as.numeric(heals$Sex == 2)
Male_best = as.numeric(best$Sex == 1)
Female_best = as.numeric(best$Sex == 2)
Cignever_heals = as.numeric(heals$cig_ever == 0)
Cignever_best = as.numeric(best$cig_ever == 0)
BMIless18pt5_heals = as.numeric(heals$bmi < 18.5)
BMIless18pt5_best = as.numeric(best$bmi < 18.5)
BMI18pt5to25_heals = as.numeric(heals$bmi >= 18.5 & heals$bmi < 25)
BMI18pt5to25_best = as.numeric(best$bmi >= 18.5 & best$bmi < 25)
BMIover25_heals = as.numeric(heals$bmi >= 25)
BMIover25_best = as.numeric(best$bmi >= 25)
snp1_0_heals = as.numeric(heals$snp1 == 0)
snp1_1_heals = as.numeric(heals$snp1 == 1)
snp1_2_heals = as.numeric(heals$snp1 == 2)
snp2_0_heals = as.numeric(heals$snp2 == 0)
snp2_1_heals = as.numeric(heals$snp2 == 1)
snp2_2_heals = as.numeric(heals$snp2 == 2)
snp3_0_heals = as.numeric(heals$snp3 == 0)
snp3_1_heals = as.numeric(heals$snp3 == 1)
snp3_2_heals = as.numeric(heals$snp3 == 2)
snp1_0_best = as.numeric(best$snp1 == 0)
snp1_1_best = as.numeric(best$snp1 == 1)
snp1_2_best = as.numeric(best$snp1 == 2)
snp2_0_best = as.numeric(best$snp2 == 0)
snp2_1_best = as.numeric(best$snp2 == 1)
snp2_2_best = as.numeric(best$snp2 == 2)
snp3_0_best = as.numeric(best$snp3 == 0)
snp3_1_best = as.numeric(best$snp3 == 1)
snp3_2_best = as.numeric(best$snp3 == 2)
Educ0_heals = as.numeric(heals$Education_year == 0)
Educ0_best = as.numeric(best$Education_year == 0)
range1to5_heals = heals$Education_year >= 1 & heals$Education_year <= 5
range1to5_best = best$Education_year >= 1 & best$Education_year <= 5
Educ1to5_heals = as.numeric(range1to5_heals)
Educ1to5_best = as.numeric(range1to5_best)
Educ6up_heals = as.numeric(heals$Education_year >= 6)
Educ6up_best = as.numeric(best$Education_year >= 6)

heals_ct = cbind(heals, Male_heals, Female_heals, Cignever_heals,
                 heals$cig_former, heals$cig_curt, BMIless18pt5_heals,
                 BMI18pt5to25_heals, BMIover25_heals, snp1_0_heals,
                 snp1_1_heals, snp1_2_heals, snp2_0_heals, snp2_1_heals,
                 snp2_2_heals, snp3_0_heals,snp3_1_heals,snp3_2_heals,
                 Educ0_heals, Educ1to5_heals, Educ6up_heals)
best_ct = cbind(best, Male_best, Female_best, Cignever_best,
                best$cig_former, best$cig_curt, BMIless18pt5_best,
                BMI18pt5to25_best, BMIover25_best, snp1_0_best, snp1_1_best,
                snp1_2_best, snp2_0_best, snp2_1_best, snp2_2_best,
                snp3_0_best,snp3_1_best,snp3_2_best, Educ0_best,
                Educ1to5_best, Educ6up_best)
heals_ct = heals_ct[,22:dim(heals_ct)[2]]
heals_ct = colSums(heals_ct, na.rm = TRUE)
heals_perc = percent(heals_ct/dim(heals)[1])
best_ct = best_ct[,19:dim(best_ct)[2]]
best_ct = colSums(best_ct, na.rm = TRUE)
best_perc = percent(best_ct/dim(best)[1])

Age_heals = summary(heals$Age)
Age_best = summary(best$Age)
lnWAs_heals = summary(heals$WArsenic)
lnUrAsgmCr_heals = summary(heals$UrAsgmCr)
lnUrAsgmCr_best = summary(best$UrAsgmCr)
wt_score_heals = summary(heals$wt_score)
wt_score_best = summary(best$wt_score)
DMApct_heals = summary(heals$DMA_pct)
lnDMA_heals = summary(heals$lnDMA)

heals_ct = c(heals_ct, Age_heals[[3]], lnWAs_heals[[3]],
             lnUrAsgmCr_heals[[3]], DMApct_heals[[3]],
             lnDMA_heals[[3]], wt_score_heals[[3]])
# Note that slashes are used in IQRs instead of commas, due to CSV format.
heals_perc = c(heals_perc, sprintf("%s/%s",Age_heals[[2]], Age_heals[[5]]),
               sprintf("%s/%s",lnWAs_heals[[2]], lnWAs_heals[[5]]),
               sprintf("%s/%s",lnUrAsgmCr_heals[[2]], 
                       lnUrAsgmCr_heals[[5]]),
               sprintf("%s/%s", DMApct_heals[[2]], DMApct_heals[[5]]),
               sprintf("%s/%s", lnDMA_heals[[2]], lnDMA_heals[[5]]),
               sprintf("%s/%s", wt_score_heals[[2]], wt_score_heals[[5]]))
best_ct = c(best_ct, Age_best[[3]], NA, lnUrAsgmCr_best[[3]], NA, NA, 
            wt_score_best[[3]])
best_perc = c(best_perc, sprintf("%s/%s", Age_best[[2]], Age_best[[5]]), NA,
              sprintf("%s/%s", lnUrAsgmCr_best[[2]], lnUrAsgmCr_best[[5]]),
              NA, NA, sprintf("%s/%s", wt_score_best[[2]], 
                              wt_score_best[[5]]))
test = rbind(heals_ct,heals_perc,best_ct,best_perc)
heals_final = paste(test[1,]," (",test[2,],")",sep="")
best_final = paste(test[3,]," (",test[4,],")",sep="")
data = rbind(test,heals_final, best_final)[5:6,]
data = t(data)
data = data[2:dim(data)[1],]
write.csv(data,file="describe.csv", quote=FALSE)

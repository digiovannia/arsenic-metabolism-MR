# Given regression and MR results for each analysis, this script computes the
# number and percent of CpGs with association/effect estimates that are 
# directionally consistent with EWAS arsenic associations, as well as binomial
# P-values. Pie charts are generated for each analysis as well, to be
# manually nested in forest plots.

library(ggplot2)

predictors = c('DMA_pct','lnDMA','recoded_bin_score','wt_score','MR',
              'MR_maxlik','X10.104623578','X10.104795134','exm1580829')

for(bonf in c(TRUE, FALSE)){
  for(k in predictors){
    if(k %in% c('DMA_pct', 'lnDMA')){
      coh = c("heals_epic","heals_meta")
    } else {
      coh = c("heals_epic","heals_meta",'best_meta','meta')
    }
    for(j in coh){
      file = read.table(sprintf("RegFiles/%s_%s_bn_lowef.txt",k,j),
                        header=TRUE,as.is=TRUE, stringsAsFactors = FALSE)
      if(j != 'meta'){
        ref = read.csv(sprintf("%s_fdr_results_betascale_2-19-2019.csv",
                               j),stringsAsFactors = FALSE)
      } else {
        ref = read.csv("heals_meta_fdr_results_betascale_2-19-2019.csv",
                       stringsAsFactors = FALSE)
      }
      if(bonf & j != 'heals_epic'){
        #Separate analysis for the Bonferroni threshold CpGs
        ref = ref[!is.na(ref$BF),]
      }
      newref = ref[,c("Name","logFC","P.Value")]
      names(newref)[1] = "cpg"
      if(bonf & j != 'heals_epic'){
        test = merge(newref,file,by='cpg')
      } else {
        test = merge(file,newref,by='cpg',all.y=TRUE)
      }
      # Opposite signs are considered consistent with causal interpretation
      Match = sign(test$effect) != sign(test$logFC)
      final = cbind(test,Match)
      if(bonf){
        write.table(final,
                    file = sprintf("consist/BONF_con_%s_%s_bn_lowef.txt",k,
                                   j), quote=F, row.names = F)
      } else {
        write.table(final, file = sprintf("consist/con_%s_%s_bn_lowef.txt",
                                          k,j), quote=F, row.names = F)
      }
      yes = sum(final$Match == TRUE)
      no = sum(final$Match == FALSE)
      df = data.frame(match = c("Consistent", "Inconsistent"),
                      count = c(yes, no))
      fulldf = df
      fulldf[3,1] = 'Consistent'
      fulldf[3,2] = df[1,2]/(df[1,2] + df[2,2])
      fulldf[4,1] = 'Consistent'
      fulldf[5,1] = 'Consistent'
      onesided_test = binom.test(df[1,2], df[1,2] + df[2,2], 0.5, 'greater')
      twosided_test = binom.test(df[1,2], df[1,2] + df[2,2], 0.5)
      fulldf[4,2] = onesided_test$p.value
      fulldf[5,2] = twosided_test$p.value
      if(bonf & j != 'heals_epic'){
        write.csv(fulldf, file=sprintf('consist/BONF_%s_%s_matches.csv',k,j))
      } else {
        write.csv(fulldf, file=sprintf('consist/%s_%s_matches.csv',k,j))
      }
      if(!bonf){
        # Generate pie charts only for the full analyses
        bp = ggplot(df, aes(x="",y=count, fill=match)) + theme_bw() +
          geom_bar(width=1,stat="identity",color='black',
                   show.legend = FALSE) +
          theme(axis.text.x=element_blank(), panel.grid = element_blank(),
                axis.ticks = element_blank()) +
          scale_fill_manual(values = c("#56B4E9","#D55E00")) +
          theme(axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.border = element_blank()) +
          labs(fill = "") + theme(legend.title=element_text(size=16)) +
          theme(legend.text=element_text(size=14)) +
          theme(plot.title = element_text(face = "bold", size = 20,
                                          hjust = 0.5))
        pie = bp + coord_polar("y", start=0)
        pdf(file = sprintf("plots/%s_cpg_by_%s_pie.pdf",j,k),
            width=6,height=6)
        print(pie)
        dev.off()
      }
    }
  }
}
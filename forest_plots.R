# Based on consistency results for each analysis, this script generates
# forest plots displaying the effect estimates and consistency (by color) of
# each CpG with the causal interpretation.

library(hash)
library(ggplot2)

predictors = c("DMA_pct","lnDMA","recoded_bin_score","wt_score","MR",
              "X10.104623578","X10.104795134","exm1580829")

titles = hash(predictors, c("DMA%", "ln(DMA)", "Binary GP-DMA%",
                             "Weighted GP-DMA%", "Mendelian Randomization",
                             "rs9527", "rs11191527", "rs61735836"))

for(k in predictors){
  if(k %in% c('DMA_pct', 'lnDMA')){
    coh = c("heals_epic","heals_meta")
  } else {
    coh = c("heals_epic","heals_meta",'best_meta','meta')
  }
  for(j in coh){
    file = read.table(sprintf("consist/con_%s_%s_bn_lowef.txt",k,j),
                      header=TRUE,as.is=TRUE)
    file = file[order(file$effect),] #rank CpGs by association estimates
    file$Match = ifelse(file$Match == TRUE, 'Consistent', 'Inconsistent')
    x = factor(file$cpg, levels=file$cpg[order(file$effect)])
    y = file$effect
    ylo = file$CI_lower
    yhi = file$CI_upper
    mat = file$Match
    colcod = ifelse(file$Match == 'Consistent', "#56B4E9", "#D55E00")
    coded_df = data.frame(x,y,ylo,yhi,mat,colcod)
    ylabel = ifelse(k == 'MR', 'MR Causal Association Estimate',
                    sprintf('%s β', titles[[k]]))
    pdf(file = sprintf("plots/%s_forest_%s.pdf",k,j), width=8.5,
        height=3.4)
    if(j == 'heals_epic'){
      coordx = 0.96
      coordy = 0.03
      legsize = 12
    } else {
      coordx =  0.95
      coordy = 0.05
      legsize = 12
    }
    ysize = ifelse(k == 'MR', 12, 14)
    forest = ggplot(coded_df, aes(x=x, y=y, ymin=ylo, ymax=yhi, colour=mat)) +
      scale_colour_manual(breaks = mat, 
                          values = c("#56B4E9", "#D55E00")) +
      theme_bw() +
      geom_pointrange() + geom_hline(yintercept = 0, linetype = 2) +
      xlab('CpG Site Ordered by Association Estimate') + ylab(ylabel) +
      labs(colour = "") +
      theme(plot.title = element_text(face = "bold", size = 18,
                                      hjust = 0.5)) +
      theme(axis.title.y=element_text(size=ysize)) +
      theme(axis.title.x=element_text(size=14)) +
      theme(legend.title=element_blank()) +
      theme(legend.text=element_text(size=legsize)) +
      theme(panel.grid.major = element_blank()) +
      theme(panel.grid.minor = element_blank()) +
      theme(legend.justification=c(1,0), legend.position=c(coordx,coordy)) +
      theme(axis.text.x = element_blank()) +
      theme(axis.ticks = element_blank())
    print(forest)
    dev.off()
  }
}  
  
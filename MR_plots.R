# Creates Mendelian randomization plots and forest plots for the three CpGs
# with the strongest causal effect estimates.

library(ggplot2)
library(hash)
library(MendelianRandomization)
options(scipen=999)

beta_DMA=c(3.59526,2.21788,5.08819)
se_DMA=c(0.47,0.34,0.51)
snp=c('rs9527','rs11191527','rs61735836')

# The following four hash maps define the positions of labels for the plots
xoff = hash(c('cg18050715', 'cg04826368', 'cg14583825'), c(0.3, 0.3, 0.3))
yoff = hash(c('cg18050715', 'cg04826368', 'cg14583825'), 
            c(0.0008, 0.0016, 0.0004))
labxoff = hash(c('cg18050715', 'cg04826368', 'cg14583825'), c(1.2, 1.2, 1.2))
labyoff = hash(c('cg18050715', 'cg04826368', 'cg14583825'),
               c(0.011, 0.023, 0.006))

j = 'meta'

medfile = read.table(sprintf("RegFiles/snp_%s_bn_lowef.txt", j),
                     header=TRUE,as.is=TRUE)
results = read.table('consist/con_MR_meta_bn_lowef.txt',
                     header=TRUE, as.is=TRUE)
results = results[results$sig == '***',]
for(c in c('cg18050715', 'cg04826368', 'cg14583825')){
  subfile = medfile[medfile$cpg == c,]
  beta_CpG = subfile$effect
  se_CpG = subfile$se
  mri_CpG = mr_ivw(mr_input(bx=beta_DMA, bxse=se_DMA, by=beta_CpG,
                            byse=se_CpG, snp=snp))
  ht = mri_CpG
  pdf(file = sprintf("plots/%s_MR_forest.pdf",c), width=5,height=4)
  forestplot_cpg = ggplot(data=NULL, aes(x=factor(snp), y=beta_CpG,
                                         ymin=(beta_CpG - 1.96*se_CpG),
                                         ymax=(beta_CpG + 1.96*se_CpG))) +
    geom_pointrange(lwd=.5) + coord_flip() + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + geom_hline(yintercept=0,
                                                           lty=1, color="red",
                                                           lwd=.4) +
    xlab('SNP') + ylab(sprintf('Per-allele Association with %s Methylation',
                               c)) +
    theme(axis.text.x=element_text(size=12, color="black")) +
    theme(axis.text.y=element_text(size=12, color="black"))
  print(forestplot_cpg)
  dev.off()
  estlabel = ifelse(c == 'cg14583825', 0.0006, signif(ht$Estimate, 2))
  pdf(file = sprintf("plots/%s_MR_plot.pdf",c), width=5,height=4)
  ht_plot = mr_plot(mr_input(bx=beta_DMA, bxse=se_DMA, by=beta_CpG,
                             byse=se_CpG, snp=snp),
                    line='ivw', interactive=F, labels=F)
  fullplot = ht_plot + theme_bw() + geom_abline(mapping=aes(slope=ht$CIUpper,
                                                            intercept=0), 
                                                linetype=2, colour="blue") +
    geom_abline(mapping=aes(slope=ht$CILower,intercept=0),
                linetype=2, colour="blue") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_text(aes(x=beta_DMA-xoff[[c]], y=beta_CpG-yoff[[c]], label=snp)) +
    geom_text(aes(x=labxoff[[c]], y=labyoff[[c]]),
              label = sprintf('Estimate = %s\nP = %s', estlabel,
                              signif(ht$Pvalue, 2)), size=4) +
    xlab("Effect Allele Association with DMA% (β)") +
    ylab(sprintf("Effect Allele Association with %s (β)",c)) +
    theme(axis.title.y=element_text(size=12)) +
    theme(axis.text.y=element_text(size=12, color="black")) +
    theme(axis.title.x=element_text(size=12)) +
    theme(axis.text.x=element_text(size=12, color="black"))
  print(fullplot)
  dev.off()
}
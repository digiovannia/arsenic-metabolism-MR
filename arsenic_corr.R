# Creates scatterplots of DMA% and log-transformed DMA against log-transformed
# creatinine-adjusted urinary arsenic.

library(ggplot2)
library(scales)
library(hash)

phenos = c('DMA_pct','lnDMA')
file = read.table("lowef_med_heals_epic_bn.txt",header=TRUE,
                  as.is=TRUE)

alignmap = hash(phenos, list(c(98, 96, 94), c(2^(9.643856), 2^(8.643856+0.67),
                                              2^(7.643856+1.34))))
phen_val_map = hash(phenos, c('DMA_pct', 'DMA'))
pheno_map = hash(phenos, c('DMA%', 'DMA'))

for(phen in phenos){
  R = round(cor(file[['lnUrAsgmCr']], file[[phen]], use='complete.obs'), 2)
  N = 383
  options(scipen = 1)
  options(digits = 2)
  P = summary(lm(file[[phen]] ~ file[['lnUrAsgmCr']]))$coef[2,4]
  
  pdf(file = sprintf("plots/arsenic_%s_corr.pdf", phen), width=5,height=6)
  corr = ggplot(data = file, aes_string(x='UrAsgmCr',y=phen_val_map[[phen]]))+
         theme_bw() + geom_point() + annotate('text', x=40,
                                              y=alignmap[[phen]][1],
                                              label=sprintf('R = %s', R),
                                              size=5) +
         annotate('text', x=40, y=alignmap[[phen]][2],
                  label=sprintf('P = %s', scientific(P, digits=2)), size=5) +
         annotate('text', x=40, y=alignmap[[phen]][3],
                  label=sprintf('N = 383'), size=5) +
         xlab("Urinary Arsenic (µg/g creatinine)") + ylab(pheno_map[[phen]]) +
         scale_x_continuous(trans = 'log2', limits=c(20,1800),
                            breaks=c(25, 50, 100, 200, 400, 800, 1600),
                            labels=c('25', '50', '100', '200', '400', '800',
                                     '1600')) +
         theme(axis.text.x=element_text(size=12, color="black")) +
         theme(axis.text.y=element_text(size=12, color="black")) +
         theme(axis.title.y=element_text(size=14)) +
         theme(axis.title.x=element_text(size=14)) +
         theme(panel.grid.major = element_blank()) +
         theme(panel.grid.minor = element_blank())
  if(phen == 'lnDMA'){
    ggplot_add(scale_y_continuous(trans = 'log2', limits=c(2,1000),
                                  breaks=c(2,4,8,16,32,64,128,256,512,1024),
                                  labels=c('2','4','8','16','32','64',
                                           '128','256','512','1024')), corr)
  }
  print(corr)
  dev.off()
}
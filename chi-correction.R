#!/usr/bin/env Rscript
library(snpStats)

##########Read Arguments###################################
args = commandArgs(trailingOnly=TRUE)
HKAtest<- read.delim(args[1])
file<- args[2] #name of output
min_sites<- args[3]

##########Prepare Data###################################
HKAtest<-HKAtest[HKAtest$poly_gene+HKAtest$fix_gene>min_sites,]

png(paste0("qqplot_chiHKA",file,".png"), width=600, height=600)
qq.chisq(HKAtest$chi2,df=1)
dev.off()
qqresult<-qq.chisq(HKAtest$chi2,df=1)

HKAtest$corrected<-HKAtest$chi2/qqresult["lambda"]

png(paste0("qqplot_chiHKA_corrected",file,".png"), width=600, height=600)
qq.chisq(HKAtest$corrected,df=1)
dev.off()

HKAtest$pval_corr<-pchisq(HKAtest$corrected, df=1, lower.tail=FALSE)
HKAtest$pval2<-pchisq(HKAtest$chi2, df=1, lower.tail=FALSE)

write.table(HKAtest,paste0(args[1],"_corrected.tab"),
            quote = FALSE, row.names = FALSE, col.names = TRUE,
            sep = "\t")

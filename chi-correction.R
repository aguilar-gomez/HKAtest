#!/usr/bin/env Rscript
library(snpStats)

##########Read Arguments###################################
args = commandArgs(trailingOnly=TRUE)
HKAtest_raw<- read.delim(args[1])
file<- args[2] #name of output
min_sites<- as.numeric(args[3])

##########Prepare Data###################################
HKAtest<-HKAtest_raw[HKAtest_raw$poly_gene+HKAtest_raw$fix_gene>min_sites,]
HKAtest$sites<-HKAtest$poly_gene+HKAtest$fix_gene

png(paste0("qqplot_chiHKA",file,".png"), width=600, height=600)
qq.chisq(HKAtest$chi2,df=1)
dev.off()
qqresult<-qq.chisq(HKAtest$chi2,df=1)

HKAtest$corrected<-HKAtest$chi2/qqresult["lambda"]

png(paste0("qqplot_chiHKA_corrected",file,".png"), width=600, height=600)
qq.chisq(HKAtest$corrected,df=1)
dev.off()

HKAtest$pval_corr<-pchisq(HKAtest$corrected, df=1, lower.tail=FALSE)
#HKAtest$pval2<-pchisq(HKAtest$chi2, df=1, lower.tail=FALSE) Sanity check: uncomment if you want to see the original p-value calculated again with this method (R vs python)

write.table(HKAtest,paste0(args[1],"_corrected.tab"),
            quote = FALSE, row.names = FALSE, col.names = TRUE,
            sep = "\t")

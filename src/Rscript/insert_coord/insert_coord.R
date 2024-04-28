#!/usr/lib/R/bin/Rscript

source('/home/kgagalova/sageflow/rlibaries/insert_plots/R/insert_plot.R')

args<-commandArgs(TRUE)

gff <- args[1]
out <- args[2]
out_summary <- args[3]

summary_all <- assignGroups(gff)

#summarize counts
summary <- data.frame(table(summary_all$Name))
colnames(summary) <- c("Name","counts")

write.table(summary_all,out,sep='\t',row.names=F,quote=F)
write.table(summary,out_summary,sep='\t',row.names=F,quote=F)

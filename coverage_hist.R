setwd("/Users/fco/Desktop/BREMEN_OP/chapter1_2/")



coverage <- read.table("coverage_table",sep = "")

barplot(coverage$V2, names.arg = coverage$V1, las=2,cex.names=0.3,ylab="mean coverage/base",main="Average coverages, 117 hamlet samples")

mean(coverage$V2)

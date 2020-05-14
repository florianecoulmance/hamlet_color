rm(list = ls())

library(ggplot2)
library(gridExtra)


setwd("/Users/fco/Desktop/BREMEN_OP/chapter1_2/")

metrics <- read.table("raw_var_sites.table.txt", sep='', header = TRUE)


mq <- ggplot(metrics, aes(x=MQ)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(50, 70) + geom_vline(xintercept=c(57.5,62.2), linetype="dotted")
#mq

#m <- metrics[which.max(metrics$MQ),]
#print(m)

qd <- ggplot(metrics, aes(x=QD)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(0, 50) + geom_vline(xintercept=c(4), linetype="dotted")

#qd

logFS <- log(metrics$FS)
fs <- ggplot(metrics, aes(x=logFS)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(-5, 15) + geom_vline(xintercept=c(1.77), linetype="dotted")
#fs

mqrs <- ggplot(metrics, aes(x=MQRankSum)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(-1, +1) + geom_vline(xintercept=c(-0.2,0.2), linetype="dotted")
#mqrs

rprs <- ggplot(metrics, aes(x=ReadPosRankSum)) + 
  geom_density(color="darkblue", fill="lightblue") + xlim(-5, +5) + geom_vline(xintercept=c(-2,2), linetype="dotted")
#rprs

figure <- grid.arrange(mq, qd, fs, mqrs, rprs, ncol = 1, nrow = 5)

ggsave(file="filtering_thresholds.pdf", figure)

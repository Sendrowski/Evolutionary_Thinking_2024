setwd("/Users/au732936/PycharmProjects/Evolutionary_Thinking_2024/week39/Friday")

source("simulateWF.R")

WF_twoalleles(5,15)

bluecounts <- WF_twoalleles(5,15)
bluefreq <- bluecounts / (2 * 5)
plot(bluefreq,ylim=c(0,1),type="b",col="blue",pch=19,xlab="generations",ylab="Blue frequency")

WF_manyalleles(5,15)

library(tidyverse)
source("geneticdrift_moi.R")
gf <- genetic_drift(N = 2*500, f_0 = 1/2, G = 1000, R = 50)
plot_genetic_drift(gf)

gf_stats(gf)

track_lineages(N.vec=rep(10,20), n.iter=1, num.tracked=3)

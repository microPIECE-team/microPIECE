#! /usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(reshape2)
input <- args[1]
output <- args[2]
svg(output)

data <- read.csv(input,header=TRUE,sep=";")


ggplot(data,aes(x=length, y=rpm, fill=type))+geom_bar(stat = "identity")+ expand_limits(y=c(0,100000))+theme_bw()

garbage <- dev.off()

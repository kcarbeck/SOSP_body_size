# katherine carbeck
# alaska song sparrow - summary stats
# july 13 2022

library(scales)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(qqman)

options(scipen = 999)

setwd("~/Desktop/PhD/Courses/Bioinformatics/Cornell - WGS Club")

#import all files

fst <- read.csv("AK/satsuma/satsuma_mean_fst_manhattan_plot.csv") # from chromosemble_manhattan_plot.R

taj_sanaka<-read.delim("AK/tajima_pi/sanaka_TD_25kb.Tajima.D",header=T)
taj_maxima<-read.delim("AK/tajima_pi/maxima_TD_25kb.Tajima.D",header=T)
taj_merrilli<-read.delim("AK/tajima_pi/merrilli_TD_25kb.Tajima.D",header=T)
taj_rufina<-read.delim("AK/tajima_pi/rufina_TD_25kb.Tajima.D",header=T)

pi_sanaka<-read.delim("AK/tajima_pi/sanaka_nucleotide_diversity_25kb.windowed.pi",header=T)
pi_maxima<-read.delim("AK/tajima_pi/maxima_nucleotide_diversity_25kb.windowed.pi",header=T)
pi_merrilli<-read.delim("AK/tajima_pi/merrilli_nucleotide_diversity_25kb.windowed.pi",header=T)
pi_rufina<-read.delim("AK/tajima_pi/rufina_nucleotide_diversity_25kb.windowed.pi",header=T)


#remove windows were variant<10
taj_sanaka<-taj_sanaka %>% filter(N_SNPS>9)
taj_maxima<-taj_maxima %>% filter(N_SNPS>9)
taj_merrilli<-taj_merrilli %>% filter(N_SNPS>9)
taj_rufina<-taj_rufina %>% filter(N_SNPS>9)

pi_sanaka<-pi_sanaka %>% filter(N_VARIANTS>9)
pi_maxima<-pi_maxima %>% filter(N_VARIANTS>9)
pi_merrilli<-pi_merrilli %>% filter(N_VARIANTS>9)
pi_rufina<-pi_rufina %>% filter(N_VARIANTS>9)


#for taj, add +1 to Binstart so they match
taj_sanaka$BIN_START <-taj_sanaka$BIN_START+1
taj_maxima$BIN_START <-taj_maxima$BIN_START+1
taj_merrilli$BIN_START <-taj_merrilli$BIN_START+1
taj_rufina$BIN_START <-taj_rufina$BIN_START+1

#colnames
colnames(fst_san_ruf)[colnames(fst_san_ruf)=="MEAN_FST"]<-"fst_san_ruf"
colnames(fst_san_mer)[colnames(fst_san_mer)=="MEAN_FST"]<-"fst_san_mer"
colnames(fst_max_ruf)[colnames(fst_max_ruf)=="MEAN_FST"]<-"fst_max_ruf"
colnames(fst_max_mer)[colnames(fst_max_mer)=="MEAN_FST"]<-"fst_max_mer"
colnames(fst_max_san)[colnames(fst_max_san)=="MEAN_FST"]<-"fst_max_san"
colnames(fst_mer_ruf)[colnames(fst_mer_ruf)=="MEAN_FST"]<-"fst_mer_ruf"

colnames(taj_sanaka)[colnames(taj_sanaka)=="TajimaD"]<-"taj_sanaka"
colnames(taj_maxima)[colnames(taj_maxima)=="TajimaD"]<-"taj_maxima"
colnames(taj_merrilli)[colnames(taj_merrilli)=="TajimaD"]<-"taj_merrilli"
colnames(taj_rufina)[colnames(taj_rufina)=="TajimaD"]<-"taj_rufina"

colnames(pi_sanaka)[colnames(pi_sanaka)=="PI"] <- "pi_sanaka"
colnames(pi_maxima)[colnames(pi_maxima)=="PI"] <- "pi_maxima"
colnames(pi_merrilli)[colnames(pi_merrilli)=="PI"] <- "pi_merrilli"
colnames(pi_rufina)[colnames(pi_rufina)=="PI"] <- "pi_rufina"


#keep only chrom, bin start, and stat value
fst_san_ruf_sub<-fst_san_ruf[c(1,2,6)]
fst_san_mer_sub<-fst_san_mer[c(1,2,6)]
fst_max_ruf_sub<-fst_max_ruf[c(1,2,6)]
fst_max_mer_sub<-fst_max_mer[c(1,2,6)]
fst_max_san_sub<-fst_max_san[c(1,2,6)]
fst_mer_ruf_sub<-fst_mer_ruf[c(1,2,6)]


taj_sanaka_sub<-taj_sanaka[c(1,2,4)]
taj_maxima_sub<-taj_maxima[c(1,2,4)]
taj_merrilli_sub<-taj_merrilli[c(1,2,4)]
taj_rufina_sub<-taj_rufina[c(1,2,4)]
  
  
pi_sanaka_sub<-pi_sanaka[c(1,2,5)]
pi_maxima_sub<-pi_maxima[c(1,2,5)]
pi_merrilli_sub<-pi_merrilli[c(1,2,5)]
pi_rufina_sub<-pi_rufina[c(1,2,5)]


#merge 
fst<-merge(fst_san_ruf_sub,fst_san_mer_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
fst<-merge(fst,fst_max_ruf_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
fst<-merge(fst,fst_max_mer_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
fst<-merge(fst,fst_max_san_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
fst<-merge(fst,fst_mer_ruf_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
quantile(fst$fst_mer_ruf,.99, na.rm=T) #0.1251543   


taj<-merge(taj_sanaka_sub,taj_maxima_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
taj<-merge(taj,taj_merrilli_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
taj<-merge(taj,taj_rufina_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)


pi<-merge(pi_sanaka_sub,pi_maxima_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
pi<-merge(pi,pi_merrilli_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)
pi<-merge(pi,pi_rufina_sub,by=c("CHROM","BIN_START"),all.x=T, all.y=T)


stats<-merge(fst,taj,by=c("CHROM","BIN_START"), all.x=T, all.y=T)  
stats<-merge(stats,pi,by=c("CHROM","BIN_START"), all.x=T, all.y=T)   
#24225 windows


#need this to make sure the BIN_STARTs order correctly
stats <- stats[order(stats$BIN_START),]

write.csv(stats, "AK/FST/stats_for_plot.csv")

#remove scaffolds that have less than 3 windows
x.tab <- table(stats$CHROM)
stats_clean<-stats[stats$CHROM %in% names(x.tab[x.tab>=4]),]
#154014 windows

write.csv(stats_clean, "AK/satsuma/stats_clean_for_plot.csv")

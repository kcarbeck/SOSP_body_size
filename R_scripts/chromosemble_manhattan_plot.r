# katherine carbeck
# jul 8 2022


library(ggplot2)
library(qqman)
library(dplyr)
library(data.table)
library(tidyverse)


#load chain file from satsuma chromosemble output
sosp.zefi.chain <- fread("satsuma_summary.chained2.out") #use the xcorr_first output 
names(sosp.zefi.chain) <- c("sosp_scaff", "sosp_start", "sosp_end", "zefi_chr", "zefi_start", "zefi_end", "V7", "dir")
sosp.zefi.chain2 <- sosp.zefi.chain %>% 
  filter(zefi_chr!="") %>% #remove empty alignments 
  mutate(chr1 = gsub("_.*", "", zefi_chr)) %>%
  mutate(chr = gsub("RRCB.*", "UN", chr1))

#see whats going on
head(sosp.zefi.chain2)
t<-sosp.zefi.chain2 %>% 
  count(chr)
# sosp_scaff sosp_start sosp_end                                                                                          zefi_chr zefi_start zefi_end V7 dir
# 1: Contig0_pilon       1504     1529 RRCB01000109.1_dna:primary_assembly_primary_assembly:bTaeGut1_v1.p:RRCB01000109.1:1:1584006:1_REF     595312   595337  1   -
# 2: Contig0_pilon      63317    63346 RRCB01000109.1_dna:primary_assembly_primary_assembly:bTaeGut1_v1.p:RRCB01000109.1:1:1584006:1_REF     968388   968417  1   -
# 3: Contig0_pilon      66490    66526 RRCB01000109.1_dna:primary_assembly_primary_assembly:bTaeGut1_v1.p:RRCB01000109.1:1:1584006:1_REF     968388   968424  1   -
# 4: Contig0_pilon      66626    66651 RRCB01000109.1_dna:primary_assembly_primary_assembly:bTaeGut1_v1.p:RRCB01000109.1:1:1584006:1_REF     968263   968288  1   -
# 5: Contig1_pilon          8       35                         25_dna:primary_assembly_primary_assembly:bTaeGut1_v1.p:25:1:3022303:1_REF      54008    54035  1   +
# 6: Contig1_pilon       1933     1957                         25_dna:primary_assembly_primary_assembly:bTaeGut1_v1.p:25:1:3022303:1_REF      63849    63873  1   +
#   chr1 chr
# 1: RRCB01000109.1  UN
# 2: RRCB01000109.1  UN
# 3: RRCB01000109.1  UN
# 4: RRCB01000109.1  UN
# 5:             25  25
# 6:             25  25

fwrite(sosp.zefi.chain2, "FORMATTED_satsuma_summary.chained.out")

sosp.zefi.chain2 <- sosp.zefi.chain2 %>% 
mutate(sosp_length = sosp_end - sosp_start)

#maxmatch summary
sosp.chain.df <- ddply(sosp.zefi.chain2,
                .(sosp_scaff,chr),
                summarize,
                MaxMatch=sum(sosp_length),
                zefi_start=min(zefi_start))

sosp.chain <- arrange(sosp.chain.df,chr,zefi_start,MaxMatch) %>% 
  mutate(chr = ifelse(MaxMatch < 5000, "NA", chr))

#select only the place where the whole contig matches best
sosp.chain <- ddply(sosp.chain,.(sosp_scaff),function(e){                
  a <- subset(e,chr==e$chr[e$MaxMatch==max(MaxMatch)])
  if(nrow(a)==1){
    a
  }
})

count(unique(sosp.zefi.chain2$chr))


gr1 <- GRanges(sosp.zefi.chain2$sosp_scaff, IRanges(sosp.zefi.chain2$sosp_start, sosp.zefi.chain2$sosp_end), zefi = sosp.zefi.chain2$chr, zefi_start = sosp.zefi.chain2$zefi_start, zefi_end = sosp.zefi.chain2$zefi_end)
head(gr1)

#
chr_order <- c("1", "1A", "2", "3", "4", "4A", "5", "6","7",  "8", 
               "9", "10", "11",  "12", "13", "14", "15", "16", "17", "18", 
               "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", 
               "29", "30","Z", "UN")

# -------------------------------  fst  -----------------------------------
max_mer <-read.delim("/Users/katherine/Desktop/PhD/Cornell - WGS Club/AK/FST/fst/50kb/FST_maxima_merrilli_50kb_mpileup.windowed.weir.fst",header=T)
max_ruf <-read.delim("/Users/katherine/Desktop/PhD/Cornell - WGS Club/AK/FST/fst/50kb/FST_maxima_rufina_50kb_mpileup.windowed.weir.fst",header=T)
san_mer <-read.delim("/Users/katherine/Desktop/PhD/Cornell - WGS Club/AK/FST/fst/50kb/FST_sanaka_merrilli_50kb_mpileup.windowed.weir.fst",header=T)
san_ruf <-read.delim("/Users/katherine/Desktop/PhD/Cornell - WGS Club/AK/FST/fst/50kb/FST_sanaka_rufina_50kb_mpileup.windowed.weir.fst",header=T)

#remove windows were variant<10
fst_max_mer <- max_mer %>% filter(N_VARIANTS>9)
fst_max_ruf <- max_ruf %>% filter(N_VARIANTS>9)
fst_san_mer <- san_mer %>% filter(N_VARIANTS>9)
fst_san_ruf <- san_ruf %>% filter(N_VARIANTS>9)

#colnames
colnames(fst_san_ruf)[colnames(fst_san_ruf)=="MEAN_FST"]<-"fst_san_ruf"
colnames(fst_san_mer)[colnames(fst_san_mer)=="MEAN_FST"]<-"fst_san_mer"
colnames(fst_max_ruf)[colnames(fst_max_ruf)=="MEAN_FST"]<-"fst_max_ruf"
colnames(fst_max_mer)[colnames(fst_max_mer)=="MEAN_FST"]<-"fst_max_mer"

#keep only chrom, bin start, bin end, and fst
fst_san_ruf_sub<-fst_san_ruf[c(1,2,3,6)]
fst_san_mer_sub<-fst_san_mer[c(1,2,3,6)]
fst_max_ruf_sub<-fst_max_ruf[c(1,2,3,6)]
fst_max_mer_sub<-fst_max_mer[c(1,2,3,6)]

#merge 
fst<-merge(fst_san_ruf_sub,fst_san_mer_sub,by=c("CHROM","BIN_START","BIN_END"),all.x=T, all.y=T)
fst<-merge(fst,fst_max_ruf_sub,by=c("CHROM","BIN_START","BIN_END"),all.x=T, all.y=T)
fst<-merge(fst,fst_max_mer_sub,by=c("CHROM","BIN_START","BIN_END"),all.x=T, all.y=T)
quantile(fst$fst_max_ruf,.99, na.rm=T) #0.3374189 


#####need this to make sure the BIN_STARTs order correctly
stats <- fst[order(fst$BIN_START),]

#remove scaffolds that have less than 3 windows
x.tab <- table(stats$CHROM)
stats_clean<-stats[stats$CHROM %in% names(x.tab[x.tab>=4]),]
#20380 windows

pops_df<-stats_clean[complete.cases(stats_clean),] #remove NAs
# join fst and chain dfs
pops.chr <- suppressWarnings(left_join(pops_df, sosp.chain, by = c("CHROM" = "sosp_scaff")))
pops.chr2 <- pops.chr %>% mutate(chr = ifelse(!chr %in% chr_order, "UN", as.character(chr))) 
unique(pops.chr2$chr)

pops.chr2 <- pops.chr2 %>% 
    mutate(chr_ordered = factor(chr, levels = chr_order)) %>% 
    arrange(chr, zefi_start) 
pops_comb<-pops.chr2
head(pops_comb)



gr1 <- GRanges(sosp.zefi.chain2$sosp_scaff, IRanges(sosp.zefi.chain2$sosp_start, sosp.zefi.chain2$sosp_end), zefi = sosp.zefi.chain2$chr, zefi_start = sosp.zefi.chain2$zefi_start, zefi_end = sosp.zefi.chain2$zefi_end)


##### find overlap with zefi coords
gr2 <- GRanges(pops_comb$CHROM, ranges = IRanges(pops_comb$BIN_START, pops_comb$BIN_START))
fo <- suppressWarnings(findOverlaps(gr1, gr2, maxgap = 650))
overlap <- cbind(sosp.zefi.chain2[queryHits(fo), ], pops_comb[subjectHits(fo), ])
overlap_srt <- overlap %>% arrange(chr, zefi_start)

overlap_srt2 <- overlap_srt %>% 
  mutate(SNP = paste(CHROM, sosp_start, sep = "_")) %>% 
  distinct(SNP, .keep_all = TRUE)
head(overlap_srt2)

chr_pops_comb <- overlap_srt2 



#get df all situated and ready to plot
fst <- chr_pops_comb %>% arrange(chr_ordered, zefi_start)
fst$row<-1:nrow(fst)
fst$rollmean <- zoo::rollmean(fst$fst_san_ruf,50,fill=NA)

chr_breaks <- fst %>% group_by(chr_ordered) %>% 
  dplyr::summarise(chr_breaks = mean(row)) %>% 
  mutate(chr_ordered = as.factor(chr_ordered),
         chr_labels = ifelse(chr_ordered %in% chr_order, as.character(chr_ordered), ".")) %>% 
  mutate(chr_labels = gsub("chromosome_", "", chr_labels))

head(fst)

write.csv(fst, "satsuma_mean_fst_manhattan_plot.csv")

# --------------- now plot manhattan plots for each pairwise comparison -------------------------
# SANAKA x RUFINA
p.sanruf <- fst %>% 
  ggplot(aes(x=row, y=fst_san_ruf, col = chr_ordered)) + theme_bw()+
  geom_hline(yintercept = quantile(fst$fst_san_ruf,.999, na.rm = TRUE), color = "black", linetype = "dashed", size=.6) + 
  scale_color_manual(values=rep(c("black","grey50"), length(levels(factor(fst$chr_ordered)))/2+1))+
  labs(y="san ruf") +
  geom_point(size=0.9,shape=19,stroke=0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_x_continuous(expand = c(.01,0), breaks = chr_breaks$chr_breaks, 
                    labels = chr_breaks$chr_labels) +
  theme_classic() +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = .5),
        panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=7),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))
p.sanruf


#outl1<-filter(fst,fst_san_ruf>=quantile(fst_san_ruf,.999,na.rm=T))

#highlight outliers
#fst_out1 <- p.sanruf +
#  geom_point(data=outl1,shape=19,stroke=0.4,size=0.8,col="red")
#fst_out1

## SANAKA X MERRILLI
p.sanmer <- fst %>% 
  ggplot(aes(x=row, y=fst_san_mer, col = chr_ordered)) + theme_bw()+
  geom_hline(yintercept = quantile(fst$fst_san_mer,.999, na.rm = TRUE), color = "black", linetype = "dashed", size=.6) + 
  scale_color_manual(values=rep(c("black","grey50"), length(levels(factor(fst$chr_ordered)))/2+1))+
  geom_point(size=0.9,shape=19,stroke=0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(y="san mer") +
  scale_x_continuous(expand = c(.01,0), breaks = chr_breaks$chr_breaks, 
                     labels = chr_breaks$chr_labels) +
  theme_classic() +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = .5),
        panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=7),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))
p.sanmer

# outl2<-filter(fst,fst_san_mer>=quantile(fst_san_mer,.9995,na.rm=T))
# #highlight outliers
# fst_out2 <- p.sanmer +
#   geom_point(data=outl2,shape=19,stroke=0.4,size=0.8,col="red")
# fst_out2


## MAXIMA X RUFINA
p.maxruf <- fst %>% 
  ggplot(aes(x=row, y=fst_max_ruf, col = chr_ordered)) + theme_bw()+
  geom_hline(yintercept = quantile(fst$fst_max_ruf,.999, na.rm = TRUE), color = "black", linetype = "dashed", size=.6) + 
  scale_color_manual(values=rep(c("black","grey50"), length(levels(factor(fst$chr_ordered)))/2+1))+
  geom_point(size=0.9,shape=19,stroke=0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(y="max ruf") +
  scale_x_continuous(expand = c(.01,0),breaks = chr_breaks$chr_breaks, 
                     labels = chr_breaks$chr_labels) +
  theme_classic() +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = .5),
        panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=7),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))
p.maxruf

# outl3<-filter(fst,fst_max_ruf>=quantile(fst_max_ruf,.9995,na.rm=T))
# #highlight outliers
# fst_out3 <- p.maxruf +
#   geom_point(data=outl3,shape=19,stroke=0.4,size=0.8,col="red")
# fst_out3


## MAXIMA X MERRILLI
p.maxmer <- fst %>% 
  ggplot(aes(x=row, y=fst_max_mer, col = chr_ordered)) + theme_bw()+
  geom_hline(yintercept = quantile(fst$fst_max_mer,.999, na.rm = TRUE), color = "black", linetype = "dashed", size=.6) + 
  scale_color_manual(values=rep(c("black","grey50"), length(levels(factor(fst$chr_ordered)))/2+1))+
  geom_point(size=0.9,shape=19,stroke=0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  labs(y="max mer") +
  scale_x_continuous(expand = c(.01,0),breaks = chr_breaks$chr_breaks, 
                     labels = chr_breaks$chr_labels) +
  theme_classic() +
  theme(legend.position="none",
        panel.border=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = .5),
        panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text = element_text(size=7),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_line(size=0.2))
p.maxmer

# outl4<-filter(fst,fst_max_mer>=quantile(fst_max_mer,.9995,na.rm=T))
# #highlight outliers
# fst_out4 <- p.maxmer +
#   geom_point(data=outl3,shape=19,stroke=0.4,size=0.8,col="red")
# fst_out4




### stitch together in one plot
library(patchwork)
quartz(height=6,width=8)
p.sanruf + p.sanmer + p.maxruf + p.maxmer + plot_layout (ncol=1,heights=c(2,2,2,2))




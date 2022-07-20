# Zoom into high fst peaks plots
# 13 july 2022
# katherine carbeck

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)

stats_clean<-read.csv("AK/satsuma/stats_clean_for_plot.csv",header=T) # made in summary_stats_alaska_sosp_satsuma.R

# sosp gff
gff_all<-read.delim("AK/FST/Identify_candidate_genes/SongSparrow_draft2.homolog.gff", header=F, comment.char="#") #full gff
gff_all<-gff_all %>% filter (V3=="mRNA") #only keep the mRNA entries since these are the tracks we are interested in
gff_all$sosp_name<-substr(gff_all$V9,4,21) #make a new column with the RBFW gene number
gff_gene_forward<-subset(gff_all,V7=="+") #forward direction genes
gff_gene_reverse<-subset(gff_all,V7=="-") #reverse direction genes


#themes
theme_noxaxis<-  theme(axis.title.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x= element_blank(),
                       axis.line.x = element_blank())
#removes x axis 
theme_noyaxis<-  theme(axis.title.y = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y= element_blank(),
                       axis.line.y = element_blank())


data_cum <- stats_clean %>% 
  dplyr::group_by(chr_ordered) %>% 
  dplyr::summarise(max_bp = max(zefi_start)) %>%
  dplyr::mutate(bp_add = lag(cumsum(max_bp),na.rm=T, default = 0)) %>%
  dplyr::select(chr_ordered, bp_add)


my_join <- stats_clean %>% 
  inner_join(data_cum, by = "chr_ordered") %>% 
  mutate(bp_cum = zefi_start + bp_add)

axis_set <- my_join %>% 
  group_by(chr_ordered) %>% 
  summarize(center = mean(bp_cum))

mydf2 <- my_join


### for the fst pairwise plot use color for large subspecies and line type for small subspecies 
# #C29ED7 = merrilli (purple)
# #5AD199 = rufina
# #F57A7A = sanaka
# #FFD114 = maxima

# san ruf = solid, #5AD199
# san mer = solid, #C29ED7
# max ruf = dashed, #5AD199
# max mer = dashed, #C29ED7



##### peak on 391 is on ZEFI chr 17 ######
mydf3<-filter(mydf2,chr_ordered==17)
max(mydf3$bp_cum)

gg391_fst<- ggplot() +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_san_ruf), size= 1.4, color="grey10", linetype= "solid") + 
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_san_ruf), size= 1,   color="#5AD199", linetype= "solid") + 
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_san_ruf), size= .8,  color="#F57A7A", linetype= "dashed") + 
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_san_mer), size=1.4, color="grey10", linetype="solid") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_san_mer), size=1,   color="#F57A7A", linetype="solid") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_san_mer), size=.8,  color="#C29ED7", linetype="dashed") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_max_ruf),size=1.4, color="grey10", linetype="solid") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_max_ruf),size=1,   color="#FFD114", linetype="solid") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_max_ruf),size=.8,  color="#5AD199", linetype="dashed") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_max_mer),size=1.4, color="grey10", linetype="solid") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_max_mer),size=1,   color="#C29ED7", linetype="solid") +
  geom_line(mydf3,mapping=aes(x=bp_cum,y=fst_max_mer),size=.8,  color="#FFD114", linetype="dashed") +
 
  scale_x_continuous(expand = c(0,0), limits=c(968200000,970873220))+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),expand = c(0,0))+
  labs(y="Fst") +
  theme_classic() +
  theme_noxaxis
gg391_fst



######  gene plot #####
zefi.gff<-readGFF("AK/zefi/Taeniopygia_guttata.bTaeGut1_v1.p.107.chr.gff3")
head(zefi.gff)
zefi.genes.gff <- zefi.genes.gff %>% 
  dplyr::select(seqid, ID, start, end, strand, Name, description, gene_id)
head(zefi.genes.gff)
subdf <- filter(zefi.genes.gff,seqid=="17")
head(subdf)
subdf <- filter(subdf, start>11514430)
subdf$mid <- (subdf$start + subdf$end)/2



gg391_genev2<-ggplot() +
  geom_segment(data=filter(subdf,strand=="+"), aes(x=start,y=.05,xend=end,yend=.05),size=1,arrow=arrow(length = unit(0.2, "cm")), color="black") +
  geom_segment(data=filter(subdf,strand=="-"), aes(x=end,y=0,xend=start,yend=0),size=1,arrow=arrow(length = unit(0.2, "cm")), color="black") +
  geom_label_repel(data=filter(subdf,strand=="+"), 
                   aes(x=mid,y=.05,label = gene_id),
                   nudge_y      = -0.5,
                   direction    = "x",
                   angle        = 90,
                   vjust        = 0.2,
                   segment.size = 0.5) +
   geom_label_repel(data=filter(subdf,strand=="-"), 
                   aes(x=mid,y=0,label = gene_id),
                   nudge_y      = -0.5,
                   direction    = "x",
                   angle        = 90,
                   vjust        = 0.2,
                   segment.size = 0.5) +
  ylim(-0.5,.1) +
  theme_classic() +
  theme_noxaxis +
  theme_noyaxis

gg391_genev2



## Nuc diversity (pi) ##

# #C29ED7 = merrilli (purple)
# #5AD199 = rufina
# #F57A7A = sanaka
# #FFD114 = maxima

gg391_pi <-ggplot() +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=pi_maxima), size=1, color="#FFD114") +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=pi_sanaka), size=1, color="#F57A7A") +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=pi_rufina), size=1, color="#5AD199") +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=pi_merrilli), size=1, color="#C29ED7") +
  scale_y_continuous(limits=c(-0.00001,0.003),breaks=seq(0,0.003,.001),expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), limits=c(968200000,970873220))+
  # xlim(0,5000001) +
  labs(y="pi") +
  theme_classic() +
  theme_noxaxis

gg391_pi

## Tajimas D ##
gg391_taj <-ggplot() +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=taj_maxima), size=1, color="#FFD114") +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=taj_sanaka), size=1, color="#F57A7A") +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=taj_rufina), size=1, color="#5AD199" ) +
  geom_line(data=mydf3,mapping=aes(x=bp_cum,y=taj_merrilli), size=1, color="#C29ED7") +
  scale_y_continuous(limits=c(-3.1,3),breaks=seq(-3,3,1),expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), limits=c(968200000,970873220))+
  # xlim(0,2000000) +
  #scale_x_continuous(expand = c(0,0))+
  labs(y="Tajima's D",x="position") +
  theme_classic() 

gg391_taj



### put plots together ###


### put plots together ###
layout <- '
AAAAA
####B
CCCCC
DDDDD
'
wrap_plots(A=gg391_fst, B=gg391_genev2, C=gg391_pi, D=gg391_taj, design = layout)




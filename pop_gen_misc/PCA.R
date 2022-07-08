library(SNPRelate)
library(tidyverse)

###new PCA with vcf done with --very-sensitive-local aligned data.
#reformats the vcf file to gds file for use in futher analysis
#run on terminal:
snpgdsVCF2GDS(vcf.fn="/workdir/kcarbeck/SOSP_alaska_merged.vcf", out.fn="all_SNPs_rem_Q30.gds",
              method = c("biallelic.only"),compress.annotation="ZIP.max", 
              snpfirstdim=FALSE, verbose=TRUE)
       

# Run all below on local computer:
snpgdsSummary("all_SNPs_rem_Q30_4pops.gds")

genofile <- snpgdsOpen("all_SNPs_rem_Q30_4pops.gds")
miss <- snpgdsSampMissRate(genofile, sample.id=NULL, snp.id=NULL, with.id=TRUE)

pca <- snpgdsPCA(gdsobj = genofile, autosome.only=FALSE)
    # Principal Component Analysis (PCA) on genotypes:
    #   Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
    # # of samples: 40
    # # of SNPs: 11,223,039
    # using 1 thread
    # # of principal components: 32
    # PCA:    the sum of all selected genotypes (0,1,2) = 438082186
    
    
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
    # [1] 19.93  6.76  4.17  3.67  3.40  2.89
tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  PC3 = pca$eigenvect[,3],
                  stringsAsFactors = FALSE)
                  
# plot
tab %>%
  mutate(group = subspecies) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + 
  xlab("PC1 19.93%")+
  ylab("PC2 6.76%")+
  geom_point(aes(color=group, size=.1)) +
  scale_color_manual(values = c("#FFD114", "#C29ED7", "#5AD199", "#F57A7A")) +
  theme_classic()

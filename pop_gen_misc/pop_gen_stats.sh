# general pop gen stats
# author:katherine carbeck

### windowed FST ####

# make an environmental variable for the input vcf
VCF=/workdir/kcarbeck/filtered_alaska_merged.vcf

# run bcftools query to get sample names
bcftools query -l $VCF
bcftools query -l $VCF | grep "sanaka" > sanaka
bcftools query -l $VCF | grep "merrilli" > merrilli
bcftools query -l $VCF | grep "maxima" > maxima
bcftools query -l $VCF | grep "rufina" > rufina

# sanaka x merrilli
vcftools --gzvcf ${VCF} \
--weir-fst-pop sanaka \
--weir-fst-pop merrilli \
--fst-window-size 50000 --out FST_sanaka_merrilli &


# maxima x rufina
vcftools --gzvcf ${VCF} \
--weir-fst-pop maxima \
--weir-fst-pop rufina \
--fst-window-size 50000 --out FST_maxima_rufina &


# maxima x sanaka
vcftools --gzvcf ${VCF} \
--weir-fst-pop maxima \
--weir-fst-pop sanaka \
--fst-window-size 50000 --out FST_maxima_sanaka &


# maxima x merrilli
vcftools --gzvcf ${VCF} \
--weir-fst-pop maxima \
--weir-fst-pop merrilli \
--fst-window-size 50000 --out FST_maxima_merrilli


# sanaka x rufina
vcftools --gzvcf ${VCF} \
--weir-fst-pop sanaka \
--weir-fst-pop rufina \
--fst-window-size 50000 --out FST_sanaka_rufina &


# merrilli x rufina
vcftools --gzvcf ${VCF} \
--weir-fst-pop merrilli \
--weir-fst-pop rufina \
--fst-window-size 50000 --out FST_merrilli_rufina &




### per SNP FST ####
# sanaka x merrilli
vcftools --gzvcf ${VCF} \
--weir-fst-pop sanaka \
--weir-fst-pop merrilli \
--out FST_sanaka_merrilli &


# maxima x rufina
vcftools --gzvcf ${VCF} \
--weir-fst-pop maxima \
--weir-fst-pop rufina \
--out FST_maxima_rufina


# maxima x sanaka
vcftools --gzvcf ${VCF} \
--weir-fst-pop maxima \
--weir-fst-pop sanaka \
--out FST_maxima_sanaka &


# maxima x merrilli
vcftools --gzvcf ${VCF} \
--weir-fst-pop maxima \
--weir-fst-pop merrilli \
--out FST_maxima_merrilli &


# sanaka x rufina
vcftools --gzvcf ${VCF} \
--weir-fst-pop sanaka \
--weir-fst-pop rufina \
--out FST_sanaka_rufina &


# merrilli x rufina
vcftools --gzvcf ${VCF} \
--weir-fst-pop merrilli \
--weir-fst-pop rufina \
--out FST_merrilli_rufina



#### calculate FIS (inbreeding coefficient)
vcftools --vcf filtered_alaska_merged.vcf --het --out AK_FIS_output



### nucleotide diversity ###
vcftools --vcf filtered_alaska_merged.vcf --indv AK_maxima1_sorted.bam --indv AK_maxima2_sorted.bam --indv AK_maxima3_sorted.bam --indv AK_maxima4_sorted.bam --indv AK_maxima5_sorted.bam --indv AK_maxima6_sorted.bam --indv AK_maxima7_sorted.bam --indv AK_maxima8_sorted.bam --indv AK_maxima9_sorted.bam --indv AK_maxima10_sorted.bam --indv AK_maxima11_sorted.bam --indv AK_maxima12_sorted.bam --window-pi 25000 --out maxima_nucleotide_diversity &

vcftools --vcf filtered_alaska_merged.vcf --indv AK_rufina1_sorted.bam --indv AK_rufina2_sorted.bam --indv AK_rufina3_sorted.bam --indv AK_rufina4_sorted.bam --indv AK_rufina5_sorted.bam --indv AK_rufina6_sorted.bam --indv AK_rufina7_sorted.bam --indv AK_rufina8_sorted.bam --indv AK_rufina9_sorted.bam --indv AK_rufina10_sorted.bam --indv AK_rufina11_sorted.bam --indv AK_rufina12_sorted.bam --window-pi 25000 --out rufina_nucleotide_diversity &

vcftools --vcf filtered_alaska_merged.vcf --indv AK_merrilli1_sorted.bam --indv AK_merrilli2_sorted.bam --indv AK_merrilli3_sorted.bam --indv AK_merrilli4_sorted.bam --indv AK_merrilli5_sorted.bam --indv AK_merrilli6_sorted.bam --indv AK_merrilli7_sorted.bam --indv AK_merrilli8_sorted.bam --window-pi 25000 --out merrilli_nucleotide_diversity &

vcftools --vcf filtered_alaska_merged.vcf --indv AK_sanaka1_sorted.bam --indv AK_sanaka2_sorted.bam --indv AK_sanaka3_sorted.bam --indv AK_sanaka4_sorted.bam --indv AK_sanaka5_sorted.bam --indv AK_sanaka6_sorted.bam --indv AK_sanaka7_sorted.bam --indv AK_sanaka8_sorted.bam --window-pi 25000 --out sanaka_nucleotide_diversity &


### TajimaD ###
vcftools --vcf filtered_alaska_merged.vcf --indv AK_maxima1_sorted.bam --indv AK_maxima2_sorted.bam --indv AK_maxima3_sorted.bam --indv AK_maxima4_sorted.bam --indv AK_maxima5_sorted.bam --indv AK_maxima6_sorted.bam --indv AK_maxima7_sorted.bam --indv AK_maxima8_sorted.bam --indv AK_maxima9_sorted.bam --indv AK_maxima10_sorted.bam --indv AK_maxima11_sorted.bam --indv AK_maxima12_sorted.bam --TajimaD 25000 --out maxima_TD &

vcftools --vcf filtered_alaska_merged.vcf --indv AK_rufina1_sorted.bam --indv AK_rufina2_sorted.bam --indv AK_rufina3_sorted.bam --indv AK_rufina4_sorted.bam --indv AK_rufina5_sorted.bam --indv AK_rufina6_sorted.bam --indv AK_rufina7_sorted.bam --indv AK_rufina8_sorted.bam --indv AK_rufina9_sorted.bam --indv AK_rufina10_sorted.bam --indv AK_rufina11_sorted.bam --indv AK_rufina12_sorted.bam --TajimaD 25000 --out rufina_TD &

vcftools --vcf filtered_alaska_merged.vcf --indv AK_merrilli1_sorted.bam --indv AK_merrilli2_sorted.bam --indv AK_merrilli3_sorted.bam --indv AK_merrilli4_sorted.bam --indv AK_merrilli5_sorted.bam --indv AK_merrilli6_sorted.bam --indv AK_merrilli7_sorted.bam --indv AK_merrilli8_sorted.bam --TajimaD 25000 --out merrilli_TD &

vcftools --vcf filtered_alaska_merged.vcf --indv AK_sanaka1_sorted.bam --indv AK_sanaka2_sorted.bam --indv AK_sanaka3_sorted.bam --indv AK_sanaka4_sorted.bam --indv AK_sanaka5_sorted.bam --indv AK_sanaka6_sorted.bam --indv AK_sanaka7_sorted.bam --indv AK_sanaka8_sorted.bam --TajimaD 25000 --out sanaka_TD &


#### relatedness ####
vcftools --vcf raw_alaska_merged.vcf --relatedness


#katherine carbeck

# first filter for missing data
vcftools --vcf filtered_alaska_merged.vcf --max-missing 1 --recode --stdout | gzip > filtered_alaska_merged_no_missing.vcf.gz

# then filter for linkage using ldPruning script from Joana Meier (ldPruning.sh)
bash ldPruning.sh filtered_alaska_merged_no_missing.vcf.gz


vcftools --gzvcf filtered_alaska_merged_no_missing.vcf.gz --plink --out filtered_alaska_merged_no_missing.LDprunedPlinkformat
       #After filtering, kept 1989848 out of a possible 1989848 Sites
/programs/plink-1.9-x86_64-beta3.30/plink --file filtered_alaska_merged_no_missing.LDprunedPlinkformat --make-bed
       #1985604 variants and 40 people pass filters and QC.


/programs/admixture/admixture --cv=5 plink.bed 1 > log1.out &  
/programs/admixture/admixture --cv=5 plink.bed 2 > log2.out & 
/programs/admixture/admixture --cv=5 plink.bed 3 > log3.out & 
/programs/admixture/admixture --cv=5 plink.bed 4 > log4.out & 
/programs/admixture/admixture --cv=5 plink.bed 5 > log5.out & 


## identify the best value of k clusters (lowest cross-validation error)
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > admixture.cv.error

# 1 0.69853
# 2 0.57874
# 3 0.60478
# 4 0.63781
# 5 0.82334


## To make plotting easier, we can make a file with the individual names in one column and the species names in the second column. As the species name is in the individual name, it is very easy to extract the species name from the individual name:
awk '{split($1,name,"_"); print $1,name[1]}' plink.nosex > plink.list


paste plink.list alaska_merged_plink.3.Q


#Then copy the "admixture.cv.error" file and all of the .Q files to your personal computer to plot in R...

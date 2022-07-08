# followed https://github.com/simonhmartin/genomics_general
# July 14, 2021

### to run you need these files in your directory:
    # SongSparrow_reference.fasta
    # filtered_alaska_merged.vcf.gz 
    # genomics.py
    # parseVCF.py
    # parseVCFs.py
    # popfile.txt
    # popgenWindows.py


## first pre-process VCF file into .geno format
## need to download parseVCFs.py

bgzip filtered_alaska_merged.vcf
tabix filtered_alaska_merged.vcf.gz

python parseVCFs.py -i filtered_alaska_merged.vcf.gz \
--skipIndels --gtf flag=DP --threads 12 |
bgzip > output.geno.gz


## next calculate diversity and divergence stats with popgenWindows.py (requires genomics.py)
## The script popgenWindows.py computes some standard population genomic statistics in sliding windows: pi, FST and DXY. It requires the script genomics.py to be present in the same directory, or in your Python path.

python popgenWindows.py -w 50000 -m 25 -g output.geno.gz -o summarystats.csv -T 5 -f phased --writeFailedWindows -p POP_1 AK_rufina1_sorted.bam,AK_rufina2_sorted.bam,AK_rufina3_sorted.bam,AK_rufina4_sorted.bam,AK_rufina5_sorted.bam,AK_rufina6_sorted.bam,AK_rufina7_sorted.bam,AK_rufina8_sorted.bam,AK_rufina9_sorted.bam,AK_rufina10_sorted.bam,AK_rufina11_sorted.bam,AK_rufina12_sorted.bam -p POP_2 AK_sanaka1_sorted.bam,AK_sanaka2_sorted.bam,AK_sanaka3_sorted.bam,AK_sanaka4_sorted.bam,AK_sanaka5_sorted.bam,AK_sanaka6_sorted.bam,AK_sanaka7_sorted.bam,AK_sanaka8_sorted.bam -p POP_3 AK_merrilli1_sorted.bam,AK_merrilli2_sorted.bam,AK_merrilli3_sorted.bam,AK_merrilli4_sorted.bam,AK_merrilli5_sorted.bam,AK_merrilli6_sorted.bam,AK_merrilli7_sorted.bam,AK_merrilli8_sorted.bam -p POP_4 AK_maxima1_sorted.bam,AK_maxima2_sorted.bam,AK_maxima3_sorted.bam,AK_maxima4_sorted.bam,AK_maxima5_sorted.bam,AK_maxima6_sorted.bam,AK_maxima7_sorted.bam,AK_maxima8_sorted.bam,AK_maxima9_sorted.bam,AK_maxima10_sorted.bam,AK_maxima11_sorted.bam,AK_maxima12_sorted.bam --popsFile popfile.txt

# Script for prepping mapped reads for variant caller
# author: katherine carbeck

# A new Bam file will be produced for each step. In the end, you only need to keep samplename_sortedRGmark.bam

# Reserved a medium machine (24 cores)
# Bam files from mapping script above are indexed and sorted already
# copy into workdir:
\ls *.bam | parallel -j 20 cp -a {} /workdir/kcarbeck &



#### Prepare the genome: index it (fai and dict files)
# R:Refernece
# O:Output
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar CreateSequenceDictionary R=SongSparrow_reference.fasta O=SongSparrow_reference.dict
samtools faidx SongSparrow_reference.fasta


#### Add or Replace Read Groups #####
# https://gatk.broadinstitute.org/hc/en-us/articles/360037226472-AddOrReplaceReadGroups-Picard-
# add read group information: Assigns all the reads in a file to a single new read-group.
# use samtools view -H sample.bam | grep '@RG'
# samtools command above will print required fields for the AddOrReplaceReadGroups below
# RGLB: Read-Group library - can be the same for all samples
# RGPL: Read-Group platform
# RGPU: run barcode (will be different for each individual; i.e., barcode sequence)
# RGSM: Read-Group sample name (use read group info from previous steps/ same as input file)


mkdir /workdir/kcarbeck/ReadGroups
[ -f /workdir/kcarbeck/ReadGroups/AddOrReplaceReadGroupsCommands.txt ] && rm /workdir/kcarbeck/ReadGroups/AddOrReplaceReadGroupsCommands.txt
cd /workdir/kcarbeck/

# awk script ran on personal computer (change directory to data folder to access metadata file)
# used to assemble script for add or replace read groups and put in a text file to run in parallel
awk -v FS="\t" '
{
    if (NR > 1){
        samp=$3
        split($5,taxon," ")
        subsp=taxon[3]
        if (!subsp) {
            subsp="x"
        }
        barcode=$15
        print "java -Xmx5g -jar /programs/picard-tools-2.8.2/picard.jar AddOrReplaceReadGroups INPUT=AK_" subsp "_" samp "_sorted.bam OUTPUT=AK_" subsp "_" samp "_sortedRG.bam RGID=" samp " RGLB=" samp " RGPL=illumina RGPU=" barcode " RGSM=AK_" subsp "_" samp
    }
}
' Mme_sequencing_metadata.txt > AddOrReplaceReadGroupsCommands.txt

parallel -j 22 < /workdir/kcarbeck/ReadGroups/AddOrReplaceReadGroupsCommands.txt


###### mark duplicates - identifies duplicate reads: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
# Metrics File: File name to write duplicate metrics to
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of file handles to keep open when spilling read ends to a desk. keep this set at 1000
# uses a lot of memory (could probably run more than 8 at a time if I decreased -Xmx to ~10g?)

for RGBAM in *_sortedRG.bam; do

  #string replacement command
  OUT=${RGBAM/%_sortedRG.bam/_sorted_RGmark.bam}
  METRICS=${RGBAM/%_sortedRG.bam/.metrics.txt}

  #create index file
  echo "java -Xmx15g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=$RGBAM OUTPUT=$OUT METRICS_FILE=$METRICS MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000" >> /workdir/kcarbeck/markDuplicatesCommands.txt

done

parallel -j 8 < /workdir/kcarbeck/markDuplicatesCommands.txt



#### Validate files
# since running HaplotypeCaller, don't need to realign or fix mates unless there is an error in ValidateSamFile

[ -f /workdir/kcarbeck/validateCommands.txt ] && rm /workdir/kcarbeck/validateCommands.txt

cd /workdir/kcarbeck/

for MARKBAM in *RGmark.bam; do

  OUTFILE=${MARKBAM/%.bam/_validate}

  echo "java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar ValidateSamFile I=$MARKBAM OUTPUT=$OUTFILE MODE=SUMMARY" >> /workdir/kcarbeck/ValidateCommands.txt

done

parallel -j 22 < /workdir/kcarbeck/ValidateCommands.txt

# concatenate all validate output into one text file to see if there are any errors
cat *validate > summary_validate.txt



#### index .bam files for HaplotypeCaller
#can run multiple at a time
[ -f /workdir/kcarbeck/IndexCommands.txt ] && rm /workdir/kcarbeck/IndexCommands.txt

cd /workdir/kcarbeck/

for MARKBAM in *RGmark.bam; do

  #create index file
  echo "java -Xmx10g -jar /programs/picard-tools-2.8.2/picard.jar BuildBamIndex I=$MARKBAM" >> /workdir/kcarbeck/IndexCommands.txt

done

parallel -j 22 < /workdir/kcarbeck/IndexCommands.txt


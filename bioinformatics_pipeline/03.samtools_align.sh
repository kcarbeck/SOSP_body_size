#author: katherine carbeck
#samtools align to ref genome

#load your sample names into an array.
#in this example I will get the sample names from the adapter removal files I generated after a previous step.

bowtie2-build -f /workdir/kcarbeck/SongSparrow_reference.fasta SOSPindex

## this was done on medium 40 core machine -- and bowtie2 step completed about 19 individuals per day ##
## samtools sort was the limiting factor -- increasing RAM makes it faster
## there were issues with some files not finishing sorting -- would have to kill them and restart...still not sure why this happens

mkdir /workdir/kcarbeck/align

INDS=($(for i in /workdir/kcarbeck/*.settings; do echo $(basename -s .settings $i); done))
#basename - will remove the directory path and returns the file name. -s tell it which suffix to remove from the end of the file name (in this case .settings)
#note: the variable INDS will now contain an array of the sample names extracted from the files names.

#If you want to see what is stored in this variable you can type:
#echo ${INDS[@]}
#@=all the elements in the array

REFERENCE=SOSPindex

[ -f /workdir/kcarbeck/align/bowtie2Commands.txt ] && rm /workdir/kcarbeck/align/bowtie2Commands.txt

for SAMPLEID in ${INDS[@]}; do
  #declare variables. This makes it easier and neater to write your command line and you just have to change these for future projects.
  ONESEQ=/workdir/kcarbeck/${SAMPLEID}.pair1.truncated
  TWOSEQ=/workdir/kcarbeck/${SAMPLEID}.pair2.truncated
  USEQ=/workdir/kcarbeck/${SAMPLEID}.collapsed,/workdir/kcarbeck/${SAMPLEID}.collapsed.truncated,/workdir/kcarbeck/${SAMPLEID}.singleton.truncated
  OUTPUTSORTED=/workdir/kcarbeck/align/${SAMPLEID}_sorted.bam

  # align with bowtie - the output is piped directly into samtools to avoid having the intermediate .sam file.
  echo "Aligning $SAMPLEID with bowtie"
  #this just writes a line telling you which sample is being worked on.
  echo "bowtie2 -p 8 --phred33 --very-sensitive-local -x $REFERENCE -I 149 -X 900 --rg-id $SAMPLEID --rg SM:$SAMPLEID -1 $ONESEQ -2 $TWOSEQ -U $USEQ|samtools sort -@ 4 -m 12G -o $OUTPUTSORTED" >> /workdir/kcarbeck/align/bowtie2Commands.txt

done

parallel -j 4 < /workdir/kcarbeck/align/bowtie2Commands.txt



# index for loop
[ -f /workdir/kcarbeck/align/samtoolsIndexCommands.txt ] && rm /workdir/kcarbeck/align/samtoolsIndexCommands.txt

cd /workdir/kcarbeck/align/

for BAM in *.bam; do

  #create index file
  echo "samtools index $BAM" >> /workdir/kcarbeck/align/samtoolsIndexCommands.txt

done

parallel -j 32 < /workdir/kcarbeck/align/samtoolsIndexCommands.txt



# qualimap for loop
# may have to increase mem limit:
#/programs/qualimap_v2.2.1/qualimap bamqc -bam SAMPLE.bam --java-mem-size=30G -outfile SAMPLE.sorted.pdf

[ -f /workdir/kcarbeck/align/qualimapCommands.txt ] && rm /workdir/kcarbeck/align/qualimapCommands.txt

cd /workdir/kcarbeck/align/

for BAM in *.bam; do

  #string replacement command
  PDF=${BAM/%.bam/.pdf}
  # qualimap
  echo "/programs/qualimap_v2.2.1/qualimap bamqc -bam $BAM -outfile $PDF" >> /workdir/kcarbeck/align/qualimapCommands.txt

done

parallel -j 2 < /workdir/kcarbeck/align/qualimapCommands.txt


# to run and get log output:
# bash qualimap.sh 2>&1 | tee qualimap_may1_$(date +%Y%m%d-%Hh%Mm%Ss).log

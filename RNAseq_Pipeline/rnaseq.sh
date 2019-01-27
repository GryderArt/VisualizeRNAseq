#!/bin/bash

###########################################

#Declaring variables

HOME='/data/khanlab/projects/ChIP_seq/RNA_DATA/'
DATA='/data/khanlab/projects/DATA/'
REF='/data/Clinomics/Ref/khanlab/Index/ctat_genome_lib_STAR_STAR-Fusion/ref_genome.fa.star.idx'
RSEM='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/hg19'
CHR='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/final-coordinates.txt'
GENOME='/data/Clinomics/Ref/khanlab/ucsc.hg19.fasta'
SAMPLE="$HOME/$1"
BAM="$SAMPLE/$1.bam"
TBAM="$SAMPLE/$1.UCSC.bam"
############################################

#loading modules

module load STAR/2.5.3a
module load samtools
module load rsem
module load R
module load igvtools

#############################################

#cd $HOME
#mkdir $1
cd $HOME/$1


if [ -f "$HOME/$1/$1_R1.fastq.gz" -a -f "$HOME/$1/$1_R2.fastq.gz" ]
then
	STAR --genomeDir $REF --readFilesIn  $HOME/$1/$1_R1.fastq.gz $HOME/$1/$1_R2.fastq.gz --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimSegmentMin 12  --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10  --alignMatesGapMax 100000  --alignIntronMax 100000  --chimSegmentReadGapMax 3 --outFileNamePrefix $1 --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 2  --outSAMunmapped Within --quantMode TranscriptomeSAM
else
    	STAR --genomeDir $REF --readFilesIn  $HOME/$1/$1_R1.fastq.gz  --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimSegmentMin 12  --chimJunctionOverhangMin 12 --alignSJDBOverhangMin 10 --alignMatesGapMax 100000  --alignIntronMax 100000  --chimSegmentReadGapMax 3 --outFileNamePrefix $1 --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 2  --outSAMunmapped Within --quantMode TranscriptomeSAM
fi


if [ -f "$SAMPLE/$1Log.final.out" ]
then
	echo "mapping completed"
else
	echo "mapping incomplete, please make sure the input files are in the $DATA directory "
fi

mv $1Aligned.toTranscriptome.out.bam $1.UCSC.bam
mv $1Aligned.sortedByCoord.out.bam $1.bam

if [ -f "$HOME/$1/$1_R1.fastq.gz" -a -f "$HOME/$1/$1_R2.fastq.gz" ]

then
	rsem-calculate-expression --no-bam-output --paired-end -p $SLURM_CPUS_PER_TASK --estimate-rspd  --bam  $TBAM  $RSEM $1 
else 
	rsem-calculate-expression --no-bam-output --single-end -p $SLURM_CPUS_PER_TASK --estimate-rspd  --bam  $TBAM  $RSEM $1
fi


awk 'NR == 1; NR > 1 {print $0 |"sort -k1"}' $1.isoforms.results > $1sorted-isoform.results

join  -t $'\t' <(sort $CHR) <(sort $1sorted-isoform.results) | tac > $1_isoformswithcoordinates.txt

awk 'BEGIN {FS="\t"; OFS="\t"} {print $3,$4,$5,$2,$1,$10}' $1_isoformswithcoordinates.txt  > $1_temp.txt 

{ printf 'Chr\tStart\tStop\tGeneID\tTranscriptID\tTPM\n';cat $1_temp.txt ; } > $1.transcript.TPM.txt

rm $1_temp.txt $1_isoformswithcoordinates.txt $1sorted-isoform.results

if [ -f "$SAMPLE/$1.genes.results" ]
then
	echo "Counts generated"
else
	echo "Please check the input bam file"
fi

samtools index $BAM

java -jar /usr/local/apps/igvtools/2.3.98/igvtools.jar count $BAM  $1.tdf $GENOME

echo "tdf file generated"

chmod -R 775 $HOME/$1

echo "permissions changed"

echo "pipeline completed"



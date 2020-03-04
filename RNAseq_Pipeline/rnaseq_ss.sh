#!/bin/bash

############################################
### RNA-seq pipeline for ChIP-seq - 2020   #
### by Berkley Gryder gryderart@gmail.com  #
### with vineela.gangalapudi@nih.gov       #
### and sivasish.sindiri@nih.gov           #
############################################

source '/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/trap.sh'
############################################

#Declaring variables

####################################
LOCAL="/lscratch/$SLURM_JOBID/"
THREADS=$SLURM_CPUS_ON_NODE
####################################

HOME='/data/khanlab/projects/ChIP_seq/RNA_DATA'
DATA='/data/khanlab/projects/DATA'
REF='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/star_ucsc_fusion_index/ref_genome.fa.star.idx'
RSEM='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/hg19'
CHR='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/final-coordinates.txt'
CHRGENELEVEL='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/final_gene_coordinates.txt'
GENOME='/data/khanlab/projects/ChIP_seq/RNA_DATA/RSEM/ucsc.hg19.fasta'
SAMPLE="$HOME/$1"
BAM="$1.bam"
TBAM="$1.UCSC.bam"

############################################

#loading modules

module load STAR/2.5.3a
module load samtools
module load rsem
module load R
module load igvtools

#############################################
echo "Sample: $HOME/$1/$1_R1.fastq.gz"

mkdir -p $LOCAL/$1
cd $LOCAL/$1
#########################

#map BAMS if not done previously
if [ -f "$HOME/$1/$1.UCSC.bam" ]
then 
	echo "mapping completed previously"
else
	if [ -f "$HOME/$1/$1_R1.fastq.gz" -a -f "$HOME/$1/$1_R2.fastq.gz" ]
	then
		echo "running STAR in paired-end mode"
		echo "with $HOME/$1/$1_R1.fastq.gz and $HOME/$1/$1_R2.fastq.gz"
		STAR --genomeDir $REF --readFilesIn  $HOME/$1/$1_R1.fastq.gz $HOME/$1/$1_R2.fastq.gz --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimSegmentMin 12  --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10  --alignMatesGapMax 100000  --alignIntronMax 100000  --chimSegmentReadGapMax 3 --outFileNamePrefix $1 --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 2  --outSAMunmapped Within --quantMode TranscriptomeSAM
	else
		echo "running STAR in single-end mode"
		STAR --genomeDir $REF --readFilesIn  $HOME/$1/$1_R1.fastq.gz  --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --chimSegmentMin 12  --chimJunctionOverhangMin 12  --alignIntronMax 100000  --chimSegmentReadGapMax 3 --outFileNamePrefix $1 --runThreadN $SLURM_CPUS_PER_TASK --outFilterMismatchNmax 2  --outSAMunmapped Within --quantMode TranscriptomeSAM
	fi

	if [ -f "$SAMPLE/$1Log.final.out" ]
	then
		echo "mapping completed"
	else
		echo "mapping incomplete, please make sure the input files are in the $DATA directory "
	fi

	mv $1Aligned.toTranscriptome.out.bam $1.UCSC.bam
	mv $1Aligned.sortedByCoord.out.bam $1.bam
fi

#count transcripts and genes if not done previously
if [ -f "$HOME/$1/$1.genes.results" ]
then 
	echo "RSEM counting completed previously"
else
	if [ -f "$HOME/$1/$1_R1.fastq.gz" -a -f "$HOME/$1/$1_R2.fastq.gz" ]

	then
		rsem-calculate-expression --no-bam-output --paired-end -p $SLURM_CPUS_PER_TASK --estimate-rspd  --bam  $TBAM  $RSEM $1 
	else 
		rsem-calculate-expression --no-bam-output  -p $SLURM_CPUS_PER_TASK --estimate-rspd  --bam  $TBAM  $RSEM $1
	fi
	#check if complete
	if [ -f "$SAMPLE/$1.genes.results" ]
	then
		echo "Counts generated"
	else
		echo "Please check the input bam file"
	fi
fi

#gene level to EDEN format
join  -t $'\t' <(sort -k1,1  $CHRGENELEVEL) <(sort -k1,1 $1.genes.results) | tac > $1_geneswithcoordinates.txt
	echo "sort and joined chromosomal locations to gene level results"
awk 'BEGIN {FS="\t"; OFS="\t"} {print $2,$3,$4,$1,$9}' $1_geneswithcoordinates.txt  > $1_genetemp.txt
	echo "rearranged column order to temp gene file"
{ printf 'Chr\tStart\tStop\tGeneID\tTPM\n';cat $1_genetemp.txt ; } > $1.gene.TPM.txt

rm $1_genetemp.txt $1_geneswithcoordinates.txt

#isoforms to COLTRON format
awk 'NR == 1; NR > 1 {print $0 |"sort -k1"}' $1.isoforms.results > $1sorted.isoform.results
	echo "sorted isoform level results"
join  -t $'\t' <(sort $CHR) <(sort $1sorted.isoform.results) | tac > $1_isoformswithcoordinates.txt
	echo "joined chromosomal locations to isoform level results"
awk 'BEGIN {FS="\t"; OFS="\t"} {print $3,$4,$5,$2,$1,$10}' $1_isoformswithcoordinates.txt  > $1_temp.txt 
	echo "rearranged column order to temp isoform file"
{ printf 'Chr\tStart\tStop\tGeneID\tTranscriptID\tTPM\n';cat $1_temp.txt ; } > $1.transcript.TPM.txt

rm $1_temp.txt $1_isoformswithcoordinates.txt $1sorted.isoform.results

#index BAM and make TDF if not done previously
if [ -f "$HOME/$1/$1.tdf" ]
then 
	echo "TDF made previously"
else
	samtools index $BAM

	java -jar /usr/local/apps/igvtools/2.3.98/igvtools.jar count $BAM  $1.tdf $GENOME

	echo "tdf file generated"
fi

#make RPM tdf
if [ -f "$HOME/$1/$1.RPM.tdf" ]
then 
	echo "RPM scaled TDF made previously"
else

	total_mapped_reads=`samtools view -c -F 260 $BAM`
	echo "total mapped reads: $total_mapped_reads"
	scale_factor=`echo "scale=5; 1000000/$total_mapped_reads" | bc `
	echo "scale TDF factor: $scale_factor"

	/data/khanlab/projects/ChIP_seq/scripts/scaleTDF.pl -i $HOME/$1/$1.tdf -o $HOME/$1/$1.RPM.tdf -f $scale_factor
	echo "RPM scaled TDF complete"
fi

####################################	
mv $1.gene.TPM.txt $HOME/$1
mv $1.transcript.TPM.txt $HOME/$1
mv $1.genes.results $HOME/$1
mv $1.isoforms.results $HOME/$1
mv $1.tdf $HOME/$1
mv $1.RPM.tdf $HOME/$1
####################################

chmod -R 775 $HOME/$1
chgrp khanlab -R $HOME/$1
echo "permissions changed"
echo "pipeline completed!"

#!/bin/sh

################################################################################
# Sample sh script for ATAC-seq parsing for Biotagging Project                 #
# © Daria Gavriouchkina                                                        #
# Sauka-Spenlger lab                                                           #
################################################################################

# Requirements:
# Programs:
# * bowtie v.1.1.2
# * bedtools v.2.15.0
# * bedGraphToBigWig from UCSC genome browser website
# * python v.2.7.5
# Files :
# bowtie index for danRer7 zebrafish index
# gtf

# Sample overview of strategy for paired end sequenicng
module load bowtie/1.1.2
module load bedtools
module load ucsctools
module load python/2.7.5

A1=read1.fastq
A2=read2.fastq
AN=name

GENOME=danRer10
CHROM=Zv9.chrom.size
GTF=Zv9_ensGene.gtf

bowtie -S -p 4 -X 2000 -m 2 $GENOME -1 $A1 -2 $A2  --chunkmb 500 $AN\.sam

samtools view -bS $AN\.sam > $AN\.bam
samtools sort $AN\.bam $AN\.sort
genomeCoverageBed -bg -split -ibam $AN\.sort.bam -g $CHROM > $AN\.bg
bedGraphToBigWig $AN\.bg $CHROM $AN\.bw

sam2bwPE_zf.pl -sam $AN\.sam -build danRer10 -name smo_$AN

samtools sort -n $AN\.bam $AN\.nsort
bedtools bamtobed -bedpe -i $AN\.nsort.bam > $AN\.nsort.bam.bed
atac_bedpe_parse2.py $CHROM $AN\.nsort.bam.bed
macs2 callpeak -t $AN\.nsort.bam.pebed -f BED --name macs2_$AN --shiftsize=100 --nomodel --slocal 1000

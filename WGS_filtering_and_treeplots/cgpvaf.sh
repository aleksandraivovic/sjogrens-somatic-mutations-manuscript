#!/bin/bash

#bsub -G team154-grp -o $patient/log_cgpvaf.%J -e $patient/err_cgpvaf.%J -q long -R'select[mem>3000] rusage[mem=3000]' -M3000 -J ${patient}.'[1-24]'
patient=$1
opt=$2

REF=reference/human/GRCh37d5 # provide your own ref
SAMPLES=$(cat $patient/samples.txt | grep $patient | tr '\n' ' ')

declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

module load cgpVAFcommand

cgpVaf.pl \
	-d Bams \
	-o $patient \
	-a $opt \
	-mq 30 \
	-bq 25 \
	-g $REF/genome.fa \
	-be .sample.dupmarked.bam \
	-bo 1 \
	-b $patient/all_muts_filtered.bed \
	-nn PDv37is \
	-tn $SAMPLES \
	-chr ${arr[$LSB_JOBINDEX-1]}

#!/bin/bash

#bsub -G team154-grp -o $patient/log.%J -e $patient/err.%J -q normal -R'select[mem>3000] rusage[mem=3000]' -M3000 -J cgpvaf'[1-24]' ./cgpvaf_concatfarm5.sh
patient=$1
opt=$2
WD=$(pwd)

REF=reference/human/GRCh37d5 # provide your own
OUT=$WD/$patient
DATA=Bams

SAMPLES=$(cat $patient/samples.txt | grep $patient | tr '\n' ' ')

mkdir -p $OUT

declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")

module load cgpVAFcommand

cgpVaf.pl \
	-d $DATA \
	-o $OUT \
	-a $opt \
	-mq 30 \
	-bq 25 \
	-g $REF/genome.fa \
	-be .sample.dupmarked.bam \
	-bo 1 \
	-b $OUT/all_muts_filtered.bed \
	-nn PDv37is \
	-tn $SAMPLES \
	-ct 1

sed '/^##/d' $OUT/PDv37is*.tsv > $OUT/${patient}.${opt}.tsv && rm $OUT/*_${opt}_vaf*

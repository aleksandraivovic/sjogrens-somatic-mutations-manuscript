sample=$1
project=$2

mkdir $sample

bsub -G team268-grp -e $sample/log.%J -o $sample/log.%J -q normal -n4 -R'select[mem>36000] rusage[mem=36000]' -M36000 -J ${sample}-sort ./sort_bam.sh $sample $project

sample=$1

bsub -G team268-grp -e $sample/log.%J -o $sample/log.%J -q normal -n4 -R'select[mem>12000] rusage[mem=12000]' -M12000 -J ${sample}-mixcr ./mixcr_analyze.sh

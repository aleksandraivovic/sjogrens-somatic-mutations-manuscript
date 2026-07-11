sample=$1

ulimit -c unlimited

echo "Starting analysis for: $sample"

#cd $sample/

mixcr analyze -Xmx40g milab-human-dna-xcr-7genes-multiplex \
 --species hsa \
 ~/scratch126/$sample/${sample}_R1.fastq \
 ~/scratch126/$sample/${sample}_R2.fastq \
 ~/scratch126/$sample/${sample}.mixcr

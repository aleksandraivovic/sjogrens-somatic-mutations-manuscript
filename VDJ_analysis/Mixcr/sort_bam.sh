sample=$1

samtools view -u -F 1280 $sample/${sample}.sample.dupmarked.bam | samtools view -u -f 1 - | samtools sort -n -o - >> $sample/${sample}.rmdup.sorted.bam

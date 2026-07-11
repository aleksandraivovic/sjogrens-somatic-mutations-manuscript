sample=$1

bedtools bamtofastq -i $sample/${sample}.rmdup.sorted.bam -fq $sample/${sample}.rmdup.r1.fastq -fq2 $sample/${sample}.rmdup.r2.fastq


For VDJ extraction from whole genome and targeted DNA sequencing, the Mixcr algorithm was used. 
https://mixcr.com/mixcr/about/

Starting from BAM files, alignments were first sorted using the bash script sort_bam.sh and its wrapper run_sortbam.sh for running as an LFS job. 

Then, the sorted BAM files were reverted for FASTQ files before realignment with Mixcr using the bash script mkfq.sh and its wrapper run_mkfq.sh. 

Finally, Mixcr was run using a preset for DNA sequences with the script mixcr_analyze.sh and its wrapper run_mixcr.sh

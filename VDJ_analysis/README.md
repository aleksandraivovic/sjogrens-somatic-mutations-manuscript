For VDJ extraction from whole genome and targeted DNA sequencing, the Mixcr algorithm was used. 

For more detail on Mixcr, please refer to https://mixcr.com/mixcr/about/.

In the Mixcr/ subdirectory are the bash scripts used to run Mixcr on bam files, including sorting and reverting to fastq files before re-alignment with Mixcr. 

Starting from BAM files, alignments were first sorted using the bash script sort_bam.sh. Tts wrapper run_sortbam.sh for running as an LFS job is provided.

Then, the sorted BAM files were reverted for FASTQ files before realignment with Mixcr using the bash script mkfq.sh and its wrapper run_mkfq.sh. 

Finally, Mixcr was run using a preset for DNA sequences with the script mixcr_analyze.sh and its wrapper run_mixcr.sh

The TGS_VDJ/ subdirectory contains Mixcr output and R scripts used for downstream analysis and visualization of clones. 

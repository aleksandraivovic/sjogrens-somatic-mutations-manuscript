The CNVkit package was used to analyze copy number events across the targeted regions outlined in the SJOGRENS_BAITS.bed file. 

https://cnvkit.readthedocs.io/en/stable/quickstart.html

The analysis was run using a “flat” reference, in the absence of normal samples. 

# build flat reference
cnvkit.py batch *Tumor.bam -n -t SJOGRENS_BAITS.bed -f hg19.fasta \
    --access data/access-5kb-mappable.hg19.bed \
    --output-reference my_flat_reference.cnn -d example2/

# analyze against reference
cnvkit.py batch *Tumor.bam -r my_flat_reference.cnn -p 0 --scatter --diagram -d example4/

The .cnr files created as output were used for downstream analysis and plotting in the “cnv_plot.R” script.

Single nucleotide variants called by CaveMan and indels called by Pindel underwent several filtering steps, as previously described in publications from our group.


A unique set of filters was applied to CaveMan calls as described in detail here, to remove artefacts specific to lasercapture microdissection and low-input library preparation, as described by Mathijs Sanders:

https://github.com/MathijsSanders/SangerLCMFiltering


Subsequently, exact binomial and beta binomial models were applied to filter out germline artefacts and noise, respectively, as previously described by Tim Coorens here:

https://github.com/TimCoorens/Unmatched_NormSeq


Further scripts for phylogenetic tree construction are also provided. The order of implementation of this pipeline is:


1. LCM-specific Caveman filters as described here: //github.com/MathijsSanders/SangerLCMFiltering

2. Re-genotyping and VAF estimation using https://github.com/cancerit/vafCorrect/tree/dev:
	cgpvaf.sh and cgpvaf_concat.sh

3. Exact binomial and beta-binomial filtering (same as in https://github.com/TimCoorens/Unmatched_NormSeq):
	filtering.R

4. MPBoot for phylogenetic tree construction from filtered variants, using the mpboot package https://github.com/diepthihoang/mpboot:
	mpboot.sh

5. Visualizing trees from MPBoot output:
	treeplots.R

6. The trinucleotide context of each sample can be plotted with the following:
	trinuc_context.R


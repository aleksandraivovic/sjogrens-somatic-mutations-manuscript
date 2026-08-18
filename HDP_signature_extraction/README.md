Signatures were extracted using the HDP package as described by Nicola Roberts: https://github.com/nicolaroberts/hdp.

For the Sjögren's analysis, we used no priors for the signature extraction, as outlined in the "HDP_bySamle_noprior.R" script and as previously described in similar work. 

Signature components were run across multiple sampling chains and were extracted using the script "hdp_extract_chains_noprior_bySamp.R"

Signature deconvolution, reattribution, and postprocessing were done using the "hdp_noprior_sigfit_postprocess.R" script as previously described, using the SigFit R package: https://github.com/kgori/sigfit.

Input files are provided to reproduce the analysis. 

# Inigo Martincorena - 2016
# Wrapper for shearwaterML - Many samples simultaneously - No matched normal
# This script runs shearwaterML on a number of regions from a bedfile
#

## Environment - Modify genome_file!

args = commandArgs(TRUE)

# Bedfile with the list of regions captured
bedfile = args[1]
entry_start = as.numeric(args[2])
entry_end  = as.numeric(args[3])
tumourbams = read.table(args[4], header=0, stringsAsFactors=F)[,1]
normalbams = read.table(args[5], header=0, stringsAsFactors=F)[,1]
outdir = args[6]

regions = read.table(bedfile, header=0, stringsAsFactors=F)
regions = regions[regions[,1] %in% c(1:22,"X","Y"),] # Selecting only valid chrs (this also protects against a bedfile header)
colnames(regions) = c("chr","start","stop")
regions = regions[entry_start:entry_end,]

# Location of genome fasta file
genome_file = "GRCh37/genome.fa";
min_totest = 2 # Minimum number of mutant reads in the tumour samples to test a site

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")


## Subfunction: Creating the count table

loadAllData_PhredQual = function (files, regions, ..., mc.cores=1, minPhred=30, samFlag=3844, mapq=25, keep=3) {

    nucleotides = c("A", "T", "C", "G", "-", "INS", "a", "t", "c", "g", "_", "ins")
    lengths = regions$stop - regions$start + 1
    rows = sum(lengths)
    beg = cumsum(c(1, lengths[-length(lengths)]))
    end = cumsum(lengths)
    coordinates = data.frame(chr = as.vector(unlist(sapply(1:nrow(regions), function(i) rep(regions$chr[i], lengths[i])))), 
                             pos = as.vector(unlist(sapply(1:nrow(regions), function(i) regions$start[i]:regions$stop[i]))))
    c = mclapply(files, function(f) {
        test.matrix <- matrix(0, ncol = length(nucleotides), nrow = rows)
        for (j in 1:nrow(regions)) {
            test.matrix[beg[j]:end[j], ] = bam2R(f, regions$chr[j], 
                regions$start[j], regions$stop[j], q=minPhred, mask=samFlag, mq=mapq, keepflag=keep, ...)[, nucleotides]
        }
        mode(test.matrix) = "integer"
        test.matrix
    }, mc.cores = mc.cores)
    counts = array(0, c(length(files), rows, length(nucleotides)))
    mode(counts) = "integer"
    for (i in 1:length(files)) counts[i, , ] = c[[i]]
    return(list(counts=counts,coordinates=coordinates,nucleotides=nucleotides))
}


## Subfunction: Estimating rho (I use grid ML search but we could use method of moments estimator as mg14 does)

estimateRho_gridml = function(x, mu) {

    rhovec = 10^seq(-6,-0.5,by=0.05) # rho will be bounded within 1e-6 and 0.32
    mm = c(x[,3:4])
    #cov = c(x[,3:4])+c(x[,1:2])
    cov = c(x[,1:2]) # << New one, Andrew found this small bug
    ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=mm, size=cov, rho=rhoj, prob=mu, log=T)))
    rhovec[ll==max(ll)][1]
}


## Loading the bam counts

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

tumcounts_obj_all = loadAllData_PhredQual(files=tumourbams, regions=regions)
normcounts_obj = loadAllData_PhredQual(files=normalbams, regions=regions)
norm_total = apply(normcounts_obj$counts[,,1:6]+normcounts_obj$counts[,,7:12], c(2,3), sum) # Global counts from the reference panel of samples
norm_total = norm_total/rowSums(norm_total) # Relative global counts (we will exclude from testing any site with >=40% mismatches in the normal panel, removing the reference bases)
tum_total = apply(tumcounts_obj_all$counts[,,1:6]+tumcounts_obj_all$counts[,,7:12], c(2,3), sum) # Global counts from the tumour samples
tum_total = tum_total/rowSums(tum_total)
sample_list = sapply(strsplit(tumourbams,"/"), function(x) substr(x[length(x)],1,nchar(x[length(x)])-4)) 

## Likelihood-Ratio Test

mutations_allsamples = NULL
for (tum in 1:length(tumourbams)) { # Iteratively running the code for the 1-sample wrapper_shearwaterML.R code

    # Recreating the 1-sample tumcounts_obj object
    tumcounts_obj = tumcounts_obj_all; tumcounts_obj$counts = array(tumcounts_obj$counts[tum,,], dim=c(1,dim(tumcounts_obj$counts[tum,,])))
    
    mutsites = which((tumcounts_obj$counts[1,,1:6] + tumcounts_obj$counts[1,,7:12])>=min_totest & norm_total<0.4, arr.ind=T) # Selecting promising sites

    if (nrow(mutsites)>0) {
    
        mutations = data.frame(sampleID=sample_list[tum], chr=tumcounts_obj$coordinates$chr[mutsites[,1]], pos=tumcounts_obj$coordinates$pos[mutsites[,1]], ref=NA, mut=tumcounts_obj$nucleotides[mutsites[,2]], xfw=NA, xbw=NA, nfw=NA, nbw=NA, mu=NA, pval=NA)
        mutations$tum_globalvaf = apply(mutsites, 1, function(x) tum_total[x[1],x[2]])
        l = length(tumcounts_obj$nucleotides)/2
        indsm = cbind(mutsites[,2], mutsites[,2]+l)
        sites_gr = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos))
        seqs = as.character(scanFa(genome_file, sites_gr))
        mutations$ref = as.character(seqs)
        
        for (j in 1:nrow(mutations)) {
     
            inds = indsm[j,]
            norcounts = normcounts_obj$counts[,mutsites[j],]
            tumcounts = tumcounts_obj$counts[,mutsites[j],]
    
            # Shearwater test
            pseudo = .Machine$double.eps
            x.fw = tumcounts[inds[1]]
            x.bw = tumcounts[inds[2]]
            n.fw = sum(tumcounts[1:l])
            n.bw = sum(tumcounts[(l+1):(2*l)])

            X.fw = sum(norcounts[,inds[1]])
            X.bw = sum(norcounts[,inds[2]])
            N.fw = sum(norcounts[,1:l])
            N.bw = sum(norcounts[,(l+1):(2*l)])
            mu = max(X.fw+X.bw,pseudo)/max(N.fw+N.bw,pseudo)
            counts = cbind(rowSums(norcounts[,1:l]),rowSums(norcounts[,(l+1):(2*l)]),norcounts[,inds])
            rho = estimateRho_gridml(counts,mu)
            
            ###############################
            #### Here we can save rho  ####
            ###############################
            
            
            rdisp = (1 - rho)/rho
    
            prob0.fw = (X.fw + x.fw)/(N.fw + n.fw); prob0.fw[prob0.fw==0] = pseudo
            prob1s.fw = x.fw/(n.fw+pseudo); prob1s.fw[prob1s.fw==0] = pseudo
            prob1c.fw = X.fw/(N.fw+pseudo); prob1c.fw[prob1c.fw==0] = pseudo
            prob1s.fw = pmax(prob1s.fw,prob1c.fw) # Min error rate is that of the population (one-sided test)
            nu0.fw = prob0.fw * rdisp; nu1s.fw = prob1s.fw * rdisp; nu1c.fw = prob1c.fw * rdisp; 
    
            prob0.bw = (X.bw + x.bw)/(N.bw + n.bw); prob0.bw[prob0.bw==0] = pseudo
            prob1s.bw = x.bw/(n.bw+pseudo); prob1s.bw[prob1s.bw==0] = pseudo
            prob1c.bw = X.bw/(N.bw+pseudo); prob1c.bw[prob1c.bw==0] = pseudo
            prob1s.bw = pmax(prob1s.bw,prob1c.bw) # Min error rate is that of the population (one-sided test)
            nu0.bw = prob0.bw * rdisp; nu1s.bw = prob1s.bw * rdisp; nu1c.bw = prob1c.bw * rdisp; 
    
            # Likelihood-Ratio Tests
            LL.fw = logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(X.fw, N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, nu1s.fw, rdisp) - logbb(X.fw, N.fw, nu1c.fw, rdisp)
            pvals_fw = pchisq(-2*LL.fw, df=1, lower.tail=F)/2 # We divide by 2 as we are performing a 1-sided test
            LL.bw = logbb(x.bw, n.bw, nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, rdisp) - logbb(x.bw, n.bw, nu1s.bw, rdisp) - logbb(X.bw, N.bw, nu1c.bw, rdisp)
            pvals_bw = pchisq(-2*LL.bw, df=1, lower.tail=F)/2 # We divide by 2 as we are performing a 1-sided test
            pvals_both = pchisq(-2*(log(pvals_fw)+log(pvals_bw)),4,low=F) # Fisher's combined p-value
        
            # Saving the result
            mutations[j,6:11] = c(x.fw, x.bw, n.fw, n.bw, mu, pvals_both)
        }
        mutations_allsamples = rbind(mutations_allsamples, mutations)
    }
}

mutations_allsamples = mutations_allsamples[which(mutations_allsamples$pval<1e-3),] # We don't need to save highly non-significant p-values

cat("outdir=",outdir, "; entry_start=", entry_start,"; entry_end=",entry_end,"\n");
system(sprintf("mkdir -p %s", outdir))
cat("Going to write to: ", sprintf("%s/mismatches_%0.0f_%0.0f.txt", outdir, entry_start, entry_end), "\n");
write.table(mutations_allsamples, file=sprintf("%s/mismatches_%0.0f_%0.0f.txt", outdir, entry_start, entry_end), col.names=T, row.names=F, sep="\t", quote=F)
cat("Done!!!\n");


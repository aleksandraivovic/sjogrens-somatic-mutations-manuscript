patient=$1
opt=$2 # "snp" or "indel"

# need to install mpboot for tree building
# filtered variants are input for mpboot

mpboot-sse-1.1.0-Linux/bin/mpboot -s ${patient}/bb_filt/${opt}_for_MPBoot.fa -bb 1000

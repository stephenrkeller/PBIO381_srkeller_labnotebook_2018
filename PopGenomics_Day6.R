setwd("~/github/PBIO381_srkeller_labnotebook_2018/myresults/")

library(OutFLANK)
library(vcfR)
library(adegenet)

geno_in <- read.fwf("OTAU_2018_reads2snps_DP10GP95_biallelic_MAF01_Miss0.8.vcf.recode.vcf.geno", width=rep(1,72))
geno <- t(geno_in)

meta <- read.table("cols_data.txt", header=T) # Bring in the metadata
meta <- meta[order(meta$sample_ID),] # Reorder the metadata by Individual Number

# Provide the transposed genotype file, locus names, and population names.
OF_SNPs <- MakeDiploidFSTMat(geno, locusNames=seq(1, 133291, by=1), popNames=meta$population) 

OF_out <- OutFLANK(FstDataFrame=OF_SNPs, LeftTrimFraction=0.05, RightTrimFraction=0.1, Hmin=0.1, NumberOfSamples=3, qthreshold=0.2)

OutFLANKResultsPlotter(OF_out, withOutliers=T, NoCorr=T, Hmin=0.1, binwidth=0.01, titletext="Scan for local selective sweeps")


outliers <- which(OF_out$results$OutlierFlag=="TRUE")
print(outliers)

# figure out which genes (transcript IDs) these SNPs actually belong to
# we need our vcf file. 
# Make sure your vcf file is in the working directory of your R session, then read it in:

vcf1 <- read.vcfR("OTAU_2018_reads2snps_DP10GP95_biallelic_MAF01_Miss0.8.vcf.recode.vcf")

# extract the info on transcript ID and SNP position with the 'getFIX' command

vcfann <- as.data.frame(getFIX(vcf1))
vcfann[outliers,]


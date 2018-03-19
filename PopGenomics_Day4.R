################
# First, plot Q values from ADMIXTURE

install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)
devtools::install_github('royfrancis/pophelper')

library(pophelper)

setwd("/Users/srkeller/github/PBIO381_srkeller_labnotebook_2018/myresults/")

admixfiles=list.files(path=("ADMIXTURE/"),full.names=T)
admixlist=readQ(files=admixfiles,filetype="basic")

# metadata: sample id and pop from beetle.pop file
metadata=read.table("beetle.pop",header=F)

# format metadata to a data frame and ind variables as chars. for plotting
metadata2=data.frame(sampleid=metadata[,1], population=metadata[,2])

metadata2$sampleid=as.character(metadata2$sampleid)
metadata2$population=as.character(metadata2$population)

# add in the sample id to the different Q files for plotting
if(length(unique(sapply(admixlist,nrow)))==1)
  admixlist <- lapply(admixlist,"rownames<-",metadata2$sampleid)

head(admixlist[[3]])

p <- plotQ(admixlist[c(1,2,3,4)],
           returnplot=T,exportplot=T,quiet=T,basesize=11, imgoutput="join", 
           showyaxis=T, showticks=T, panelspacer=0.4, useindlab=F, showindlab=F, 
           grplab=metadata2[2], grplabsize=4, linesize=1, barsize=1, pointsize=3, 
           panelratio=c(4,1), divsize = 0.75, pointcol="white", showtitle=T, 
           titlelab="ADMIXTURE analysis on O. tauri, SNPs thinned to 1 per kb", 
           splab=c("K=2","K=3","K=4","K=5"), outputfilename="ADMIXTURE_Otauri",
           imgtype="pdf",height=3,width=25)

plot(p$plot[[1]])

#########################
# Second, run PCA on SNP data

install.packages("vcfR")
install.packages("adegenet")

library(vcfR)
library(adegenet)

vcf1 <- read.vcfR("OTAU_2018_reads2snps_DP10GP95_biallelic_MAF01_Miss0.8.vcf.recode.vcf")

gl1 <- vcfR2genlight(vcf1)
print(gl1)
gl1$ind.names
gl1$loc.names[1:10]
gl1$chromosome[1:3]
gl1$pop

#Nothing in the field termed "pop"...need to fix that:

meta <- read.table("cols_data.txt", header=T)

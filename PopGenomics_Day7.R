###################################################
### R script to plot and parse Bayescan results ###
###       PBIO/BIO381 -- Spring 2018            ###
###################################################

# First, scp your results files from the server to your local directory
# Then set the working directory in R to your local path where you saved your results

setwd ("~/github/PBIO381_srkeller_labnotebook_2018/myresults/Bayescan/")

# Import file with SNP ID's so you can re-associate them with the Bayescan output

sites <- read.table("sites_87K.txt", header=T)

# Import likelihood scores from the mcmc chain to assess sufficient burn-in and lack of autocorrelation

LL <- read.table ("OTAU_filtered.recode.baye.sel", header=T)

# Let's plot the log-likelihood from the 10k model to assess convergence after the burn-in:
plot(LL_$logL, type="o", col="gray", ylab="Log-Likehood")


# *** NOW, PLOT THE LIKELIHOOD FROM THE OTHER RUNS AND COMPARE: *** #
# *** PLOT ALL 3 RUNS ON THE SAME MULTI-PANEL FIGURE USING par(mfrow=c(,)) *** #
# *** TRY ZOOMING IN ON THE X-AXIS TOO LOOK FOR AUTOCORRELATION USING xlim=c(,) *** #


# If you're satisfied that the runs converged after burn-in, then proceed to look at the alpha and fst levels for each locus

# First, import Bayescan results for different levels of prior odds
Fst <- read.table ("OTAU_filtered.recode.baye_fst.txt", header=T)

# ONLY If you ran multiple models #
#Compare alpha from the 3 different models:
pairs(cbind(pr100$alpha, pr_1K$alpha, pr_10K$alpha), 
      labels=c("alpha_pr100", "alpha_pr1K", "alpha_pr10K"))


# Now, let's put results together with locus names from the sites file
Fst_sites <- cbind(sites, Fst)

# And take the -log10 of the q-values, which is the typical way they're reported for ease of visualizing
Fst_sites$log10q <- -log10(Fst_sites$qval)

# Use the False Discover Rate (FDR) to control multiple tests
# An FDR of 0.01 means that 1% of our positive results will be false positives 

FDR <- 0.01

# Plotting a 2 panel figure of alpha and fst
par(mfrow=c(1,2))
plot(Fst_sites$fst, Fst_sites$alpha, xlab="Fst", ylab="Alpha")
plot(Fst_sites$log10q, Fst_sites$alpha, xlab="-log10(q-value)", ylab="Alpha")
alphacand <- Fst_sites[which(Fst_sites$log10q>(-log10(FDR))),]
points(alphacand$log10q, alphacand$alpha, col="red", pch=16)
title("Bayescan: FDR = 0.01", outer=T)

# Now, plot the results in a 4 panel figure showing the distribution of alpha and how the signficiance of alpha vs. the associated q-values
par(mfrow=c(2,2))
hist(Fst_sites$alpha, breaks=500, col="blue", xlab="Alpha", ylab="Frequency", main="", xlim=c(quantile(Fst_sites$alpha, 0.025),quantile(Fst_sites$alpha, 0.975)))
hist(Fst_sites$fst, breaks=500, col="red", xlab="Fst", ylab="Frequency", main="", xlim=c(quantile(Fst_sites$fst, 0.025),quantile(Fst_sites$fst, 0.975)))
plot(Fst_sites$fst, Fst_sites$alpha, xlab="Fst", ylab="Alpha")
plot(Fst_sites$log10q, Fst_sites$alpha, xlab="-log10(q-value)", ylab="Alpha")
alphacand <- Fst_sites[which(Fst_sites$log10q>(-log10(FDR))),]
points(alphacand$log10q, alphacand$alpha, col="red", pch=16)
title("Bayescan: FDR = 0.01", outer=T)


# Get a list of the SNPs that are significant at a given prior odds and FDR:
candsnps_pr_Fst <- Fst_sites[which(Fst_sites$log10q>-log10(FDR)),]


# *** HOW MANY ARE THERE? WHICH TRANSCRIPTS? *** #
dim()
head()

# Now, export these candidates for annotation

write.table(candsnps_pr_Fst, "bayescan_candsnps.txt", quote=F, row.names=F, col.names=T)

# This script is for simulating expression data give a list of genes, their mean,
# and their standard deviation. To test is gene expression is different between 
# 2 groups, say control and mutant, one must know each group's mean, standard deviation,
# and sample size to conduct a t-test. If we only have summary information per gene,
# as we have in the "expr" data frame below, it is impossible to directly calculate
# p-value. However, as we do know the mean and standard deviation for each one, we can
# sample from a simulated distribution to obtain multiple technical replicates.
# Once we have multiple simulated repicates, we can then perform statistical tests
# for differences such as a t-test, or differential expression analysis.

#clear all data and detach any packages
pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list = ls(all=TRUE))

expr = data.frame(Gene = c("EZH2", "FOXA2", "HMX1", "HOXA5", "HOXA6", "PAX5", "S16"),
                  Mean_Cont = c(30.22627, 30.178621, 23.289217, 22.996525, 29.920652, 12.64944, 24.40246),
                  SD_Cont = c(1.1, 1.0, 1.5, 0.9, 1.2, 2.1, 1.1), 
                  Mean_Mutant = c(15.11313, 18.08931, 10.64461, 18.49826, 18.96033, 33.29888, 24.40246),
                  SD_Mutant = c(1.5, 0.8, 1.1, 1.0, 1.3, 1.8, 0.9))
print(expr)

get_counts_from_mean = function(x, n, features, mu, dv,
                                fname=deparse(substitute(features)),
                                mname=deparse(substitute(mu)),
                                sname=deparse(substitute(dv))
                                ){
  # print(expr[,sname])
  mat_cont = matrix(0, length(expr[,mname]), n)
  for (i in 1:length(expr[,mname])) {
    mat_cont[i,] = rnorm(n = n, mean = expr[,mname][i], sd = expr[,sname][i])
  }

  sim_cols <- paste0(mname, 1:n)
  colnames(mat_cont) <- sim_cols
  rownames(mat_cont) <- expr[,fname]

  return(mat_cont)
}

#Simulate the data
#
# Typically, with at least 20-30 replicates (either biological or simulated),
# a distribution tends toward a normal distribution, based on Central Limit Theorem.
# So, we will simulate 25 samlples, then test some of the distributions
n = 25
# simulation of control group
sim_control = get_counts_from_mean(expr, n, Gene, Mean_Cont, SD_Cont)
#simulation of "mutant" or "treatment" group
sim_mutant = get_counts_from_mean(expr, n, Gene, Mean_Mutant, SD_Mutant)
counts = cbind(sim_control, sim_mutant)

#now lets aggregate everything into a SummarizedExperiment object for easy analysis
require(SummarizedExperiment)
rowData = DataFrame(Gene = rownames(counts))

colData <- DataFrame(Treatment=c(
                     rep("control", length(colnames(sim_control))),
                     rep("mutant", length(colnames(sim_mutant)))
                     ),
                     row.names=colnames(counts))

se = SummarizedExperiment(assays=list(counts=counts),
                     rowData=rowData, colData=colData)

##univatriate visualization of FOXA2 gene in control group
test = assay(se)[rowData(se)$Gene == "FOXA2", colData(se)$Treatment == "control"]
hist(test)

#standard qqplot
qqplot(rnorm(n), test)

#qqplot with confidence interval
library(ggpubr)
ggqqplot(test)

#if p-vlue is significant, it means the test distribution is different that normal dist
# shapiro.test function test can only handle up to 5000 numbers in input.
# if you need to simulate more replicates (say 250,000), use the Kolmogorov-Smirnov test
# listed below
shapiro.test(test)

# similar to shapiro.test, if p-value is different, it means the 2 distributions
# are different (i.e. test dist and normal dist are different).
# The parameter "pnorm" is the cummulative normal distribution.
# You don't want to use something like "dnorm", which is the density function 
# of the normal distribution 
ks.test(test,"pnorm") #notice this is incorrectly showing a difference
ks.test(test,"pnorm", mean(test), sd(test)) #correct way

# we now have a simmulated data, which follows a normal distribution in the form of a 
# SummarizedExperiment (SE) object.
# We can use this SE object for statistics or as input for  differential expression 
# programs like DESeq2 or EdgeR
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
# biocLite("rhdf5")
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("scater")

# source("https://bioconductor.org/biocLite.R")
# biocLite("SRAdb")

# source("https://bioconductor.org/biocLite.R")
# biocLite("EnsDb.Hsapiens.v86")

library(EnsDb.Hsapiens.v86)
txdb <- EnsDb.Hsapiens.v86
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene

library(SRAdb)
require(readr)
library(tximport)
require(scater)

dir = "C:/Users/Chris/Documents/working/kallisto"
dir
files = Sys.glob(file.path(dir, "*.tsv"))
files
prefix = gsub(".*(SRR\\d*)_abundance\\.tsv", "\\1", files)
prefix
length(files)
#names(files) <- paste0(prefix, "_", 1:length(files))
names(files) <- paste0(prefix)
names(files)

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
#txi.kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene)
txi.kallisto

prj <- as.data.frame(read_csv("working/prjna236018_runinfo.csv"))

prj = subset(prj, Run!="Run")
rownames(prj) = prj$Run

prj = prj[ rownames(prj) %in% names(files), ]
prj


colData <- DataFrame(prj)
head(colData)

sce = SingleCellExperiment(assays = list(counts = txi.kallisto$counts),
                           colData = colData)
sce

exprs(sce) <- log2(
  calculateCPM(sce, use.size.factors = FALSE) + 1)
sce

keep_feature <- rowSums(exprs(sce) > 0) > 0
sce <- sce[keep_feature,]

colData(sce)$SampleGroup = gsub("_\\d+$", "", colData(sce)$SampleName)


sce = calculateQCMetrics(sce, feature_controls = list(eg = 1:40))

# scater_gui(sce)

plotQC(sce, type = "highest-expression", exprs_values = "counts")
plotQC(sce, type = "exprs-freq-vs-mean")
#plotQC(sce, type = "expl", method = "pairs", theme_size = 6)
plotQC(sce, type = "find-pcs", variable = "total_features",
       plot_type = "pcs-vs-vars")

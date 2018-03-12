require(scater)
data("sc_example_counts")
data("sc_example_cell_info")

head(sc_example_counts)

example_sce <- SingleCellExperiment(
  assays = list(counts = sc_example_counts), colData = sc_example_cell_info)

example_sce
counts(example_sce)
colnames(example_sce)

exprs(example_sce) <- log2(
  calculateCPM(example_sce, use.size.factors = FALSE) + 1)

keep_feature <- rowSums(exprs(example_sce) > 0) > 0
example_sce <- example_sce[keep_feature,]
example_sce

require(readr)
library(tximport)
dir = "/home/rhooads/local_kal"
dir
runs = list.files(dir)
runs
files = Sys.glob(file.path(dir, runs, "*.tsv"))
files
length(files)
names(files) <- paste0("sample", 1:length(files))
names(files)

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
txi.kallisto

txi.kallisto$counts

sce_example = readTxResults(samples = names(files),
                            files = files, type = "kallisto", txOut = TRUE)
sce_example


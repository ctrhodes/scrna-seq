###
# Download sra data

pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list = ls(all=TRUE))

# source("https://bioconductor.org/biocLite.R")
# biocLite("SRAdb") #, suppressUpdates=TRUE
# install.packages("RCurl")

library(RCurl)
#options(download.file.method="libcurl")
library(SRAdb)
mainDir = setwd("~/working")
mainDir = setwd("~/working")
mainDir

subDir <- "sra"
subDir
if (!grepl(subDir, mainDir)) {
  sra_path = file.path(mainDir, subDir)
  sra_path
}


if( !dir.exists(sra_path) ) {
  dir.create(sra_path)
} else {
  cat("sra directory already exists!")
}


if( ! file.exists('SRAmetadb.sqlite') ) {
  sqlfile <- getSRAdbFile()
} else {
  sqlfile = 'SRAmetadb.sqlite'
}

sqlfile <-'SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile)
sra_tables <- dbListTables(sra_con)
sra_tables
dbListFields(sra_con,"study")

setwd(sra_path)
rs = listSRAfile( c("SRP041736"), sra_con, fileType ='sra')
head(rs)
downloads = gsub("\\.sra", "", list.files(sra_path))
downloads
remaining = setdiff(rs$run, downloads)
remaining
# ri = getSRAinfo( c("SRP041736"), sra_con, sraType = "sra" )
# head(ri)
getSRAfile( c(remaining), sra_con, fileType ='sra')

##
# system list files and get file size
# move small files to a dir
# big files to diff dir

# # fix later
# while (length(remaining) > length(downloads) {
#   getSRAinfo( c("remaining"), sra_con, sraType = "sra" )
#   remove files with size 0
# rs = listSRAfile( c("SRP041736"), sra_con, fileType ='sra')
# head(rs)
# downloads = gsub("\\.sra", "", list.files(sra_path))
# downloads
# remaining = setdiff(rs$run, downloads)
# remaining
# }

# sra download:
# add non-zero exit status try as well as while loop 
# while (length(remaining) > length(downloaded))


################
# Move small files to "small file" dir
# Move large files to "large file" dir



sra_size
max(sra_size)
file.size(sras)

require(gdata)
library(filesstrings)
list.files()

sra_size = file.size(Sys.glob(file.path(sra_path, "*.sra")))
names(sra_size) = basename(Sys.glob(file.path(sra_path, "*.sra")))
sra_size
humanReadable(sra_size, standard = "IEC")
size_index = which(humanReadable(small_files, standard = "IEC") < "40 MiB")
small_files = names(small_files[size_index])
file.move( small_files, paste0(sra_dir, "/", "newdir") )


#########
# Download sra toolkit

library(RCurl)

getwd()
list.files()
toolkit_dir = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-"
if( Sys.info()['sysname'] == "Windows" ) {
  URL = paste0(toolkit_dir, "win64.zip")
} else if( Sys.info()['sysname'] == "Linux" ) {
  URL = paste0(toolkit_dir, "ubuntu64.tar.gz")
} else {
  URL = paste0(toolkit_dir, "mac64.tar.gz")
}

URL

if ( grepl("\\.zip", basename(URL)) ) {
  download.file(URL, destfile=basename(URL), method="libcurl")
  unzip(basename(URL))
} else {
  download.file(URL, destfile=basename(URL), method="libcurl")
  untar(basename(URL))
}

file.move( "C:\\Users\\LinLab\\Documents\\working\\sra\\sratoolkit.2.9.0-win64\\bin\\fastq-dump.exe", 
           paste0(sra_path, "/", "low_reads") )

setwd(file.path(sra_path, "low_reads"))
list.files()
dumps = Sys.glob(file.path(sra_path, "low_reads", "*.sra"))
dump_prefix = gsub("\\.sra", "", basename(dumps))
dump_prefix
for (i in 1:length(dump_prefix)) {
  print(dump_prefix[i])
  system( paste0("fastq-dump ", dump_prefix[i]) )
}


if( !dir.exists(file.path(mainDir, "fastq")) ) {
  dir.create(file.path(mainDir, "fastq"))
} else {
  cat("fastq directory already exists!")
}

library(filesstrings)

fastqs = Sys.glob(file.path(sra_path, "low_reads", "*fastq"))

file.move( fastqs, file.path(mainDir, "fastq") )

#########
# Download transcriptome and kallisto aligner


list.files()
kall_dir = "https://github.com/pachterlab/kallisto/releases/download/v0.44.0/"
if( Sys.info()['sysname'] == "Windows" ) {
  URL = paste0(kall_dir, "kallisto_windows-v0.44.0.zip")
} else if( Sys.info()['sysname'] == "Linux" ) {
  URL = paste0(kall_dir, "kallisto_linux-v0.44.0.tar.gz")
} else {
  URL = paste0(kall_dir, "kallisto_mac-v0.44.0.tar.gz")
}

URL

if ( grepl("\\.zip", basename(URL)) ) {
  download.file(URL, destfile=basename(URL), method="libcurl")
  unzip(basename(URL))
} else {
  download.file(URL, destfile=basename(URL), method="libcurl")
  untar(basename(URL))
}

# fix
# file.delete(basename(URL))
setwd(mainDir)
list.files()

#we are using the refSeq mRNA fasta file to reduce memory requirements. if you have enough memory, try mRNA.fa.gz
mrna_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz"
mrna_URL
download.file(mrna_URL, destfile=basename(mrna_URL), method="libcurl")
require(R.utils)
gunzip(basename(mrna_URL))

mrna = gsub("\\.gz", "", basename(mrna_URL))
mrna

ercc_URL = "https://raw.githubusercontent.com/roryk/tiny-test-data/master/genomes/ERCC/ERCC92/seq/ERCC92.fa"
ercc_URL
download.file(ercc_URL, destfile=basename(ercc_URL), method="libcurl")
basename(ercc_URL)

list.files()
getwd()
file.copy("refMrna.fa", "refMrna_ercc92.fa")
file.append("refMrna_ercc92.fa", "ERCC92.fa")
file.move("refMrna_ercc92.fa", file.path(mainDir, "kallisto"))

setwd(file.path(mainDir, "kallisto"))
list.files()
system( "kallisto index -i refMrna_ercc92.idx refMrna_ercc92.fa" )

#alternatively, use runKalliston in scater package
for (i in 1:length(list.files()) {
  system( "kallisto quant -i $INDEX_DIR/refMrna_ercc92.idx -o $OUT_DIR/${describer} -b 10 ${describer}$RD1_SUFFIX.fastq.gz ${describer}$RD2_SUFFIX.fastq.gz" )
}

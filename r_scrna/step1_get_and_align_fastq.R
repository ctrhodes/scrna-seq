
# set system to not sleep, especially if using windows

# source("https://bioconductor.org/biocLite.R")
# biocLite("scater") #, suppressUpdates=TRUE
# biocLite("SRAdb") #, suppressUpdates=TRUE
# install.packages("RCurl")

pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list = ls(all=TRUE))

library(RCurl)
library(SRAdb)
require(XML)
require(xml2)
library(rvest)
require(dplyr)
require(scater)

setRelWd <- function(rel_path){
  curr_dir <- getwd()
  abs_path <- file.path(curr_dir,rel_path)
  if(dir.exists(abs_path)){
    setwd(abs_path)
  }
  else
  {
    warning('Directory does not exist. Please create it first.')
  }
}

createRelWd <- function(rel_path){
  curr_dir <- getwd()
  abs_path <- file.path(curr_dir,rel_path)
  if(!dir.exists(abs_path)){
    dir.create(abs_path)
  }
  else
  {
    warning('Directory already exists.')
  }
}

goHomeDir <- function(){
  while (path.expand("~") != getwd() ) {
    message( paste0("leaving: ", getwd()) )
    setwd("..")
  }
  message( paste0("current directory is: ", getwd()) )
}


createRelWd("project")

setRelWd("project")
getwd()

if( ! file.exists('SRAmetadb.sqlite') ) {
  sqlfile <- getSRAdbFile()
} else {
  sqlfile = 'SRAmetadb.sqlite'
}

sra_con <- dbConnect(SQLite(),sqlfile)
sra_tables <- dbListTables(sra_con)
sra_tables
dbListFields(sra_con,"study")

study_accession = "SRP041736"

psi = listSRAfile( c(study_accession), sra_con, fileType ='sra')
psi
myGetSRAinfo = function (in_acc, sra_con, sraType = "sra") 
{
  sraFile <- listSRAfile(in_acc, sra_con = sra_con, fileType = sraType, 
                         srcType = "ftp")
  sraFileDir <- paste(na.omit(unique(dirname(sraFile$ftp))), 
                      "/", sep = "")
  sraFileBase <- na.omit(unique(basename(sraFile$ftp)))
  file_name = NULL
  file_size = NULL
  file_date = NULL
  require(RCurl)
  opts = curlOptions(header = TRUE)
  for (sraFileDir_1 in sraFileDir) {
    x <- getURL(sraFileDir_1, .opts = opts)
    x1 <- strsplit(x[1], "\n")[[1]]
    x2 <- sub("^.*\\s+anonymous\\s+", "", x1, perl = TRUE)
    file_name <- c(file_name, paste(sraFileDir_1, sub("^.*\\s+", 
                                                      "", x2, perl = TRUE), sep = ""))
    file_date1 = sub("^\\d+\\s{1}", "", x2, perl = TRUE)
    file_date <- c(file_date, substr(file_date1, 1, gregexpr("\\s", 
                                                             file_date1, perl = TRUE)[[1]][4] - 1))
    file_size <- c(file_size, ceiling(as.integer(sub("\\s+.*$", 
                                                     "", x2, perl = TRUE))/1024))
    Sys.sleep(0.5)
  }
  file_loc = paste0(file_name, sraFileBase)
  print(file_loc)
  file_info <- as.data.frame(cbind(file_name = file_loc, `size(KB)` = file_size,
                                   date = file_date))
  sraFileInfo <- merge(sraFile, file_info, by.x = "ftp", by.y = "file_name",
                       all.x = TRUE)
  return(sraFileInfo)
}

ps = myGetSRAinfo( c(study_accession), sra_con = sra_con)
head(ps)


# get metadata
#############

# dbListTables(sra_con)
tblnames = as.vector(dbGetQuery(
  sra_con,"select tbl_name from sqlite_master where type='table' and tbl_name not like '%ft%';")[,1])
tblnames

dbListFields(sra_con,"sra")

study_info = dbGetQuery(sra_con, 
                        paste("SELECT run_accession, library_name, sample_attribute", 
                              "FROM sra WHERE",
                              paste0("study_accession='", study_accession, "'"))
)
head(study_info)
require(tidyr)
study_info$cell_type = gsub("_\\w*", "", study_info$library_name)

tissue <- study_info$cell_type
tissue[grepl("K562", tissue) | grepl("HL60", tissue) | 
             grepl("2339", tissue)] <- "blood"
tissue[grepl("BJ", tissue) | grepl("Kera", tissue) | 
             grepl("2338", tissue)] <- "dermal"
tissue[grepl("iPS", tissue)] <- "pluripotent"
tissue[grepl("GW21+3", tissue) | grepl("GW21", tissue) | 
             grepl("NPC", tissue) | grepl("GW16", tissue)] <- "neural"

study_info = cbind(study_info, tissue)

study_info= inner_join(study_info, ps, by = c("run_accession" = "run") )
head(study_info)

# # metadata directly from Pollen

URL = "https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/pollen.rds"
if( ! file.exists(basename(URL)) ) {
  download.file(URL, destfile=basename(URL), method="libcurl", mode="wb")
}
pollen = readRDS("pollen.rds")

colData(pollen)$library_name = rownames(colData(pollen))
colData(pollen)$library_name = gsub("Hi_", "", colData(pollen)$library_name)
colData(pollen)$library_name = gsub("GW21.2", "GW21+3", colData(pollen)$library_name)

study_info$library_name
colData(pollen)$library_name

study_info = study_info[study_info$library_name %in% colData(pollen)$library_name, ]
head(study_info)
dim(study_info)

for (i in 1:length(study_info$sample_attribute)) {
  study_info$sample_attribute[i] = paste(unlist(strsplit(study_info$sample_attribute[i], " \\|\\| ")), collapse = ",")
}
study_info$`size(KB)` = as.numeric(as.character(study_info$`size(KB)`))
head(study_info)

study_info = study_info[with(study_info, order(library_name, -`size(KB)`)), ]
head(study_info)

write.table(study_info, paste0(study_accession, "_study_info.txt"), 
            quote = TRUE, sep = "\t",
            row.names = FALSE, col.names = TRUE)


# Finished metadata
#####################


study_info = read.table("SRP041736_study_info.txt", header = TRUE)

# remove low coverage files
use_cells = study_info[!duplicated(study_info$sample), ]
head(use_cells)
dim(use_cells)


createRelWd("fastq")
setRelWd("fastq")

myGetFASTQfile = function (in_acc, sra_con, destDir = getwd(), srcType = "ftp", 
          makeDirectory = FALSE, method = "curl", ascpCMD = NULL) 
{
  sraFiles = getFASTQinfo(in_acc, sra_con, srcType)
  if (makeDirectory == TRUE && !file.exists(destDir)) {
    tryCatch(dir.create(destDir), error = function(err) {
      stop("failed to create '", destDir, "': ", conditionMessage(err))
    })
  }
  message("Files are saved to: \n'", destDir, "'\n")
  if (srcType == "ftp") {
    if (missing(method)) 
      method <- ifelse(!is.null(getOption("download.file.method")), 
                       getOption("download.file.method"), "auto")
    fnames <- sraFiles$ftp
    fileinfo = NULL
    for (i in fnames) {
      try(download.file(i, destfile = file.path(destDir, basename(i)), 
                    method = method))
    }
  }
  else if (srcType == "fasp" & !is.null(ascpCMD)) {
    ascpR(ascpCMD, sraFiles$fasp, destDir)
  }
  return(sraFiles)
}

# myGetFASTQfile( c(low_reads$run), sra_con, srcType = 'ftp')

# accs = low_reads$run
accs = use_cells$run_accession

#check if we missed any
dnload_files = regmatches( list.files(), regexpr("SRR[0-9]{7}.*fastq.gz", list.files()))
downloads = unique( gsub("(SRR[0-9]{7}).*", "\\1", dnload_files) )
downloads

remaining = setdiff(accs, downloads)
remaining

while (length(remaining) > 0) {
  dnload_files = regmatches( list.files(), regexpr("SRR[0-9]{7}.*fastq.gz", list.files()))
  downloads = unique( gsub("(SRR[0-9]{7}).*", "\\1", dnload_files) )
  
  remaining = setdiff(accs, downloads)
  print(remaining)
  myGetFASTQfile( c(remaining), sra_con, srcType = 'ftp' )
}

# move non-paired end samples to different folder
single_end = regmatches( list.files(), regexpr("SRR[0-9]{7}.fastq.gz", list.files()) )
single_end

require("filesstrings")
file.move( single_end, file.path(getwd(), "single_end") )


#get file size of _1.fastq.gz and _2.fastq.gz and compare. Make sure size of one is ~ .8-1.2 times size of other

setRelWd("..")

# for windows only, check if perl is installed and in PATH
# linux and mac already have perl in path
is_perl_installed = function(){
tryCatch({shell("perl -v", intern = TRUE)}, 
         warning=function(e){cat("WARNING : perl not installed or not in PATH \n
                               Either add perl to PATH with add_perl_path() \n
                               Or install perl with install_perl() \n",conditionMessage(e), "\n")})
}
is_perl_installed()

install_perl <- function() { 
  n <- readline(prompt="Enter install type; i for interactive, s for silent: ")
  if(!grepl("^[iIsS]$",n))
  {
    return(install_type())
  }
  
  if( Sys.info()['sysname'] == "Windows" ) {
    URL = "http://strawberryperl.com/download/5.26.1.1/strawberry-perl-5.26.1.1-64bit.msi"
    if( ! file.exists(basename(URL)) ) {
      download.file(URL, destfile=basename(URL), method="libcurl", mode="wb")
      }
    }
  
  if ((n == "i") | (n == "I")) {
    CMD = paste0("msiexec.exe /i ", basename(URL))
    print(CMD)
    shell( CMD )
  } else {
    CMD = paste0("msiexec.exe /a ", basename(URL), " /qn /norestart /log strawberry.log")
    print(CMD)
    shell( CMD )
  }
}
install_perl()

add_perl_path = function(){
  if ( grepl("perl", shell("PATH", intern = TRUE)) ) {
    print("perl already in PATH!")
  } else {
    loc = shell("where perl", intern = TRUE)
    CMD = paste0( "set PATH=", loc,";%PATH%" )
    shell(CMD)
  }
}
add_perl_path()

#
# install and run fastqc on windows
# if you are using linux or mac, just use fastqcr package

fqc_dir = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"
if( Sys.info()['sysname'] == "macOS" ) {
  URL = paste0(fqc_dir, "fastqc_v0.11.7.dmg")
} else {
  URL = paste0(fqc_dir, "fastqc_v0.11.7.zip")
}

if ( grepl("\\.zip", basename(URL)) ) {
  download.file(URL, destfile=basename(URL), method="libcurl")
  unzip(basename(URL))
} else {
  download.file(URL, destfile=basename(URL), method="libcurl")
#  untar(basename(URL))
}

# run from directory containing fastq files


# either make symbolic link to fastqc dir (below), add .pl to fastqc extension, or add fastq dir as fcq() argument
# add input function in fqc() to ask user which to do

setRelWd("fastq")

# this block is only needed for windows 7 or windows 8 OS
# unix/ linux can call fastqc wrappers directly, or call fastqc directly in bash
span = 30
brks = seq(1, length(dnload_files), span)
require(filesstrings)
for (b in 1:length(brks)){
  print(brks[b])
  if( !dir.exists( paste0("sub_", brks[b])) ) {
    dir.create( paste0("sub_", brks[b]) )
  } else {
    cat("fastq subdirectory already exists!")
  }
  toMove = regmatches( list.files(), regexpr("SRR[0-9]{7}.*fastq.gz", list.files()))
  print( file.path(getwd(), toMove[1:span]) )
  try(
  file.move( file.path(getwd(), toMove[1:span]), file.path(getwd(), paste0("sub_", brks[b])) )
  )
}


subdirs = regmatches( list.dirs(), regexpr("sub_.*", list.dirs()))
subdirs
fq_dir = file.path(getwd(), subdirs)
fq_dir
fastqc_dir = file.path(path.expand("~/project"), "FastQC")
fastqc_dir
qc_dir = file.path(path.expand("~/project/fastq/FASTQC"))
qc_dir

# FastQC runs ALL samples in parallel. Only use if you have small file sizes or small number of large files

fqc = function ( fq_dir = getwd(), qc_dir = NULL, fastqc_dir = getwd() ) {
  if (is.null(qc_dir)) 
    qc_dir <- file.path(fq_dir, "FASTQC")
  if( !dir.exists(qc_dir)) {
    dir.create(qc_dir)
  } else {
    cat("FASTQC directory already exists!")
  }
  # interval = 1
  files = Sys.glob(file.path(fq_dir[interval], "*.fastq.gz"))
  for (i in 1:length(files)){
    CMD = paste0( "perl ", paste0(fastqc_dir, "/"), "fastqc ",  files[i], " --outdir ", qc_dir)
    print(CMD)
    shell(CMD, intern = TRUE)
    # tasks = shell("tasklist", intern = TRUE)
    # task_num = regmatches( tasks, regexpr("java.exe", tasks))
    # while ( length(task_num) > 1 ) {
    #   Sys.sleep(60)
    # }
    # interval += 1
  }
}

interval = 9:10
fqc(fq_dir = fq_dir, qc_dir = qc_dir, fastqc_dir = fastqc_dir)



require(fastqcr)
qc <- qc_aggregate(file.path(getwd(), "fastq", "FASTQC"))
qc

# Inspecting QC Problems
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# See which modules failed in the most samples
qc_fails(qc, "module")
fail_seq_qual = qc_fails(qc, "module", compact = FALSE) %>% 
  filter(module == "Per base sequence quality") %>% 
  select(sample)

fail_seq_qual$sample


#https://jermdemo.blogspot.com/2010/03/soft-trimming-in-r-using-shortread-and.html
#softTrim
#trim first position lower than minQuality and all subsequent positions
#omit sequences that after trimming are shorter than minLength
#left trim to firstBase, (1 implies no left trim)
#input: ShortReadQ reads
#       integer minQuality
#       integer firstBase
#       integer minLength
#output: ShortReadQ trimmed reads
library("ShortRead")
softTrim<-function(reads,minQuality,firstBase=1,minLength=5){
  qualMat<-as(FastqQuality(quality(quality(reads))),'matrix')
  qualList<-split(qualMat,row(qualMat))
  ends<-as.integer(lapply(qualList,
                          function(x){which(x < minQuality)[1]-1}))
  #length=end-start+1, so set start to no more than length+1 to avoid negative-length
  starts<-as.integer(lapply(ends,function(x){min(x+1,firstBase)}))
  #use whatever QualityScore subclass is sent
  newQ<-ShortReadQ(sread=subseq(sread(reads),start=starts,end=ends),
                   quality=new(Class=class(quality(reads)),
                               quality=subseq(quality(quality(reads)),
                                              start=starts,end=ends)),
                   id=id(reads))
  
  #apply minLength using srFilter
  lengthCutoff <- srFilter(function(x) {
    width(x)>=minLength
  },name="length cutoff")
  newQ[lengthCutoff(newQ)]
} 

setRelWd("fastq")
list.files()
for (i in 1:length(fail_seq_qual$sample)) {
  in_name = paste0(fail_seq_qual$sample[i], ".fastq.gz")
  out_name = paste0(fail_seq_qual$sample[i], ".trimmed.fastq.gz")
  print(i)
  print(in_name)
  print(out_name)
  reads<-readFastq(in_name)
  trimmedReads<-softTrim(reads=reads,minQuality=10,firstBase=1,minLength=20)
  writeFastq(trimmedReads,file=out_name)
}

qc_report(file.path(getwd(), "FASTQC"), result.file = "multi-qc-report" )

# qc report indicates smallest fastq reads are 25 bp
# well use this to build our kmer in kallisto index
qc_df = qc_stats(qc)

createRelWd("failed_fastq")

fail_qc = paste0(fail_seq_qual$sample, ".fastq.gz")

require("filesstrings")
file.move( fail_qc, file.path(getwd(), "failed_fastq") )

setwd("..")
getwd()



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

getwd()
createRelWd("genome")
setRelWd("genome/")

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
file.copy("Homo_sapiens.GRCh38.91.cdna.all.fa", "GRCh38.91_ercc92.fa")
file.append("GRCh38.91_ercc92.fa", "ERCC92.fa")
setwd("..")
require(filesstrings)
file.move(file.path(getwd(), "genome", "refMrna_ercc92.fa"), file.path(getwd(), "kallisto"))

setRelWd("kallisto")
list.files()

shell( "kallisto index -i refMrna_ercc92.idx -k 25 refMrna_ercc92.fa" )

setwd("..")
getwd()
setRelWd("fastq")
getwd()
list.files()

srr_samples = regmatches(list.files(), regexpr(".*\\.fastq.*", list.files()))
srr_samples
sample_names = unique( gsub(".*(SRR[0-9]{7}).*", "\\1", srr_samples) )
srr_base
srr_read1 = regmatches(srr_samples, regexpr(".*SRR[0-9]{7}_1.*", srr_samples))
srr_read1
srr_read2 = regmatches(srr_samples, regexpr(".*SRR[0-9]{7}_2.*", srr_samples))
srr_read2

targets = data.frame(sample_names, srr_read1, srr_read2)
head(targets)
write.table(targets, "targets_file.txt", row.names = FALSE, quote = FALSE, sep = "\t")
list.files()
setwd("..")

# Run kallisto
# is planning to use kallisto-sleuth, set n_bootstrap_samples > 0, otherwise leave 0
setRelWd("kallisto")
createRelWd("quant")
getwd()
list.files()

targets_file = file.path(path.expand("~/project/fastq"), "targets_file.txt")
targets_file

# format ensebl headers: GRCh38.91_ercc92.idx = gsub(" cdna.*", "", GRCh38.91_ercc92.idx)
transcript_index = "GRCh38.91_ercc92.idx"
transcript_index

# can run without generating bam files, only transcript abundance
runKallisto(targets_file = targets_file, transcript_index = transcript_index, single_end = FALSE,
              output_prefix = "quant/output", n_bootstrap_samples = 5)


# can also run kallisto and simultaneously generate bam files for RSeQC
# This option is also good if you run into memory problems with standard runKallisto() function
createRelWd("quant_pseudobam")
targets_info = read.table(targets_file, header = TRUE)
head(targets_info)

for (i in 1:length(targets_info$sample_names)) {
  cmd = paste0("kallisto quant -i ", transcript_index, " -o ",
               file.path(getwd(), "quant_pseudobam", paste0("output_", targets_info$sample_names[i])),
               " -b 5 --bias --pseudobam ",
               file.path(path.expand("~/project/fastq"), targets_info$srr_read1[i]), " ",
               file.path(path.expand("~/project/fastq"), targets_info$srr_read2[i])
  )
  print(cmd)
  shell(cmd)
}


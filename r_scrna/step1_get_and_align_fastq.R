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

# Move small files to "small file" dir
if( !dir.exists(file.path(sra_path, "low_reads")) ) {
  dir.create(file.path(sra_path, "low_reads"))
} else {
  cat("low_reads directory already exists!")
}

# Move large files to "large file" dir
if( !dir.exists(file.path(sra_path, "high_reads")) ) {
  dir.create(file.path(sra_path, "high_reads"))
} else {
  cat("high_reads directory already exists!")
}

require(gdata)
library(filesstrings)
list.files()

sra_size = file.size(Sys.glob(file.path(sra_path, "*.sra")))
names(sra_size) = basename(Sys.glob(file.path(sra_path, "*.sra")))
sra_size
humanReadable(sra_size, standard = "IEC")
size_index = which(humanReadable(small_files, standard = "IEC") < "40 MiB")
small_files = names(small_files[size_index])
file.move( small_files, file.path(sra_path, "low_reads" )


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

file.move( file.path(sra_path, "sratoolkit.2.9.0-win64/bin", "fastq-dump.exe"), 
          file.path(sra_path, "low_reads")

setwd(file.path(sra_path, "low_reads"))
list.files()
dumps = Sys.glob(file.path(sra_path, "low_reads", "*.sra"))
dumps
dump_prefix = gsub("\\.sra", "", basename(dumps))
dump_prefix
for (i in 1:length(dump_prefix)) {
#  print(dump_prefix[i])
  cmd = paste0("fastq-dump -I --split-3 ", dump_prefix[i])
  print(cmd)
  system( cmd )
}


if( !dir.exists(file.path(mainDir, "fastq")) ) {
  dir.create(file.path(mainDir, "fastq"))
} else {
  cat("fastq directory already exists!")
}

library(filesstrings)

# we only want to move files ending in either _1.fastq or _2.fastq. All other fastq files contain only orphaned reads that
# are not part of a matched pair
fastq1 = Sys.glob(file.path(sra_path, "low_reads", "*_1.fastq"))
fastq2 = Sys.glob(file.path(sra_path, "low_reads", "*_2.fastq"))
fastqs = c(fastq1, fastq2)
file.move( fastqs, file.path(mainDir, "fastq") )

# for windows only, check if perl is installed and in PATH
# linux and mac already have perl in path
is_perl_installed = function(){
tryCatch({system("perl -v", intern = TRUE)}, 
         error=function(e){cat("ERROR : perl not installed or not in PATH \n
                               Either add perl to PATH with add_perl_path() \n
                               Or install perl with install_perl()",conditionMessage(e), "\n")})
}
is_perl_installed()

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
    CMD = paste0("msiexec.exe /a ", basename(URL), " /qb /norestart /log strawberry.log")
    print(CMD)
    shell( CMD )
  }
}
install_perl()


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

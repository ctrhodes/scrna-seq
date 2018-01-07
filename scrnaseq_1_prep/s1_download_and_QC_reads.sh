startup up google cloud compute vm by cloning the following and following readme instructions:

https://github.com/ctrhodes/gce-startup.git

for step 7, run following script. The -f/--full option installs a gnome desktop on the vm allowing you to view HTML from fastQC reports:

gce-startup.sh -f

install fastQC:

ubuntu 16.04 LTS typically ships with java 1.8. but the google compute image does not seem to include this. So install java:

sudo apt-get install default-jre

and check version:

java -version

install fastQC in your "run" directory that is in your PATH (This should not be needed for java, but I like to be consistent where I install all my program binaries)

download fastQC:

curl -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.6.zip

install unzip to unzip fastqc:

sudo apt install unzip

unzip fastqc_v0.11.6.zip

make fastqc executable:

chmod 755 fastqc

install multiqc

multiQC is a python module. for anything python, best to use anaconda and its package manager conda:

cd /tmp

curl -O https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh

bash Anaconda2-5.0.1-Linux-x86_64.sh

source ~/.bashrc

conda list
conda search "^python$"
python --version
conda install -c bioconda multiqc

last we need some fastq files, as it is much cheaper to use previously published data on NCBIs SRA, we will download a project of low depth sequencing in the mammalian brain. I has 327 samples, but I stopped the downloaded (control-c) part way through to save time:

install programs

install edirect:

cd ~
/bin/bash
perl -MNet::FTP -e \
  '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
   $ftp->login; $ftp->binary;
   $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
gunzip -c edirect.tar.gz | tar xf -
rm edirect.tar.gz
builtin exit
export PATH=$PATH:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"
./edirect/setup.sh

add edirect to PATH:
echo "export PATH=\$PATH:\$HOME/edirect" >> $HOME/.bash_profile
source $HOME/.bash_profile


install sra-toolkit:
if not done already, create "run" directory:
mkdir -p $HOME/run
add "run" directory to path:
echo 'export PATH=$PATH:$HOME/run' >> ~/.bashrc
move into "run" for sra toolkit download:
cd $HOME/run

download "sratoolkit":
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

extract:
tar -vxzf sratoolkit.tar.gz

add sratoolkit to $PATH

check correctly installed:
which fastq-dump

install gnu parralel:
sudo apt-get install parallel


###
#download project/ multi fastq-dump
###

choose a project on Gene Expression Omnibus (GEO), usually individual samples found as GSE numbers in publications, then find bioproject. In this example, bioproject is PRJNA236018. First I make a folder to store sra project and partition at least 500GB of space on my disk.

mkdir ~/sra
cd ~/sra

QUERY="PRJNA236018"
esearch -db sra -query $QUERY | efetch --format runinfo |cut -d "," -f 1 > SRR.numbers

the SRR.numbers file has multiple blank lines and lines with the value "Run". Just ignore these, gnu parallel will skip them if there is an error

wc -l SRR.numbers
head SRR.numbers
SRR_LIST=$(cat SRR.numbers)
echo $SRR_LIST
parallel --jobs 3 "fastq-dump --split-files --origfmt --gzip {}" ::: $SRR_LIST

inspect for erroneous file downloads and re-download fastq files if needed.

run fastqc on all files:
cd ~/sra
fastqc *fastq

make a QC folder and move all fastqc output files there:
mkdir ~/sra/QC
mv *fastqc* ./QC/

cd ~/sra/QC
then run multiqc

multiqc .

this will produce a multiqc. html file where you can inspect the graphs of all fastq files at once interactively.

for my example, it generates a list of files not passing read quality, which I then remove in the ~/sra folder using grep *file | rm However, I do this to save time for this tutorial, and in the real world, I would likely trim the low quality bases off the read in question with Trimmomatic

We now have a single cell experiment that has read QC completed.

Now we go on to alignment, which I will demonstrate with STAR, but could use kallisto or salmon.


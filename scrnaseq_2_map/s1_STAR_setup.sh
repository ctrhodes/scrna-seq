#map high quality reads to genome with STAR

can either use the binaries (best) or build from scratch
# Get latest STAR bin and source from releases
wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar -xzf 2.5.3a.tar.gz
cd STAR-2.5.3a

cd STAR-2.5.3a/bin/Linux_x86_64

make STAR program executable:
chmod 755 STAR
create a sybomlic link (short-cut) to STAR

ln -S ./STAR $HOME/run/

test if installed:
STAR

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git
cd STAR/source

# Build STAR
make STAR

Another option, if you don't mind using an older version of STAR is to use apt install:
sudo apt-get install rna-star


#####
#build genome index
mkdir ~/genomes
cd ~/genomes

#download genome sequence (fasta) and annotation (gtf). GTF is technically a GFF extension. So when programs require GFF file (as with STAR genome assemly), often they mean GTF. I personally prefer GTF over older GFF formats
#doesn't matter which source you use for genome/annotation (i.e. ensembl, ucsc, nih, etc) just make sure you use the same source for fastq and gtf. I am using ensembl for everything:

wget http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/ENSEMBL/homo_sapiens/ENSEMBL.homo_sapiens.release-83/Homo_sapiens.GRCh38.dna.primary_assembly.fa

wget http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/ENSEMBL/homo_sapiens/ENSEMBL.homo_sapiens.release-83/Homo_sapiens.GRCh38.83.gtf

#Also, the demo project uses 92 ERCC spike-ins, so we have to add the sequence and annotation entries to the files:

git clone https://github.com/NCBI-Hackathons/HASSL_Homogeneous_Analysis_of_SRA_rnaSequencing_Libraries.git

cd HASSL_Homogeneous_Analysis_of_SRA_rnaSequencing_Libraries/Spike-Ins/ERCC92

mv ERCC92.fa ~/genomes
mv ERCC92.gtf ~/genomes

merge files together
cd ~/genomes/
cat *.fa > Homo_sapiens.GRCh38.83_ERCC.fa

cat *.gtf > Homo_sapiens.GRCh38.83_ERCC.gtf

now we need to add a dir for STAR to add is genome indexes

mkdir STAR_genomes

at this point it might be good in crease the number of cores and memory avaiable, requires at least 30 GB memory and you should use multiple cores. In this case, I am using a standard N1 vm with 16 cores and 60GB memory. After the index is built (about an hour), I will reduce resources to save money.

If you analyse multiple species (i.e. mouse and human), build each genome index in a seperate directory. For example, ~/STAR_genomes/mouse/ and ~/STAR_genomes/human. This is because the index files are all named SA. Also, for sjdOverhang, look for "sequence read length" in multiqc file, and take the largest one minus one (in my case, reads of 100 and 101 bp, so 101-1=100):

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/STAR_genomes --genomeFastaFiles ~/genomes/Homo_sapiens.GRCh38.83_ERCC.fa --sjdbGTFfile ~/genomes/Homo_sapiens.GRCh38.83_ERCC.gtf --sjdbOverhang 100

Now we align:
make output dir for aligned files:
mkdir ~/sra/aligned

run STAR, I am reducing my vm resources to use 8 cores and 30GB memory:




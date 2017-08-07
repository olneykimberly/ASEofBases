#!/bin/bash

#############
# What: This set of bash commands downloads initial data and does some processing
#
#               1. Convert Mapability file to bed format and create file to filter on mapability
#               2. Get protein coding gene annotations
#               3. Download & parse sample information from 1000 genomes individuals
#               4. Download & parse sample information from Geuvadis RNAseq data
#
#############



# Define directories
#Location of logfile
LOG_FILE=$PWD/get_sh.log
touch $LOG_FILE
#!! Changed absolute paths to relative paths
progdir="$(dirname "$PWD")"/2_prog
rawdir="$(dirname "$PWD")"/3_raw
mkdir -p "$(dirname "$PWD")"/3_raw #!! Make raw directory if it doesn't exist already, -p means mkdir if it doesn't exist

#Define locations of files, can be URLs or absolute paths to local files
gencodeURL=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
kgenomeURL=ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
kgenotypeURL=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz
gevodisURL=http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt

#Define functions
#function that logs everything it is given and prints it to the screen as well
function logIt {
    echo "[$(date -u)]: $*"
    echo "[$(date -u)]: $*" >> $LOG_FILE
}

function checkIt {
  if[[ -s $1]]; then
    logIt $1 looks good! Moving on
  else
    logIt $1 does not look right, exiting
    exit 1
  fi
}

logIt Program directory is $progdir
logIt Raw data directory is $rawdir


#############
## 2. Get protein coding gene annotations
#############
# This data is to develop a set of filters so that we can analyze only variants
# in protein coding regions.
#
# gencode.v19.annotation

cd $rawdir

# Download or get gencode gene annotation file
logIt Downlading gencode gene annotation file from $gencodeURL
wget $gencodeURL
checkIt gencode.v19.annotation.gtf.gz

#Keep only protein coding gene entries:
#!! zcat on ubuntu
logIt Keep only protein coding gene entries
gzcat $rawdir/gencode.v19.annotation.gtf.gz | awk '{if($3=="gene" && $20=="\"protein_coding\";"){print $0}}' > $rawdir/gencode.protein_coding.genes.v19.gtf

checkIt $rawdir/gencode.protein_coding.genes.v19.gtf

### TAKING CHR FROM GTF FILE
logIt Taking CHR from GTF File
for chr in {1..22}
do
pattern="chr$chr"
egrep -e $pattern $rawdir/gencode.protein_coding.genes.v19.gtf | sort -n -k4 -k5 | cut -f1,4,5,9 > $rawdir/gencode.chr$chr.gore
checkIt $rawdir/gencode.chr$chr.gore
done

pattern="chrX"
egrep -e $pattern $rawdir/gencode.protein_coding.genes.v19.gtf | sort -n -k4 -k5 | cut -f1,4,5,9 > $rawdir/gencode.chrX.gore
checkIt $rawdir/gencode.chrX.gore

pattern="chrY"
egrep -e $pattern $rawdir/gencode.protein_coding.genes.v19.gtf | sort -n -k4 -k5 | cut -f1,4,5,9 > $rawdir/gencode.chrY.gore
checkIt $rawdir/gencode.chrY.gore

#############
## 3. Download & parse sample information from 1000 genomes individuals
#############

logIt Download and parse sampel information from 1000 genomes individuals

# Download full 1000 genomes individuals list
logIt Downloading 10000 gnomes list from $kgenomeURL
wget $kgenomeURL
checkIt integrated_call_samples_v3.20130502.ALL.panel

# Make sure 1000 genomes sample list is in the raw data directory
logIt Moving 1000 gnomes sample list to raw directory
mv integrated_call_samples_v3.20130502.ALL.panel $rawdir/kg.integrated_call_samples_v3.20130502.ALL.panel

# Cut the first column of data
logIt Cutting the first column of data
cut -f1 $rawdir/kg.integrated_call_samples_v3.20130502.ALL.panel > $rawdir/kg.inds
checkIt $rawdir/kg.inds

# For each chromosome, download the 1000 genomes genotype file
#!! Only downloads one file, removing extra lines
logIt Downloading 1000 gnomes genotype file from $kgenotypeURL
wget $kgenotypeURL
checkIt ALL.chr13.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz

#############
## 4. Download & parse sample information from Geuvadis RNAseq data
#############

logIt Download and parse sample information from Geuvadis RNAseq data

# download Geuvadis sample info and move to raw data directory
logIt Downloading sample geuvadis data from $gevodisURL
wget $gevodisURL
mv E-GEUV-1.sdrf.txt $rawdir/geuvadis.E-GEUV-1.sdrf.txt

# parse to get list of bamfile names
logIt Parse to get list of bamfile names
cut -f1,29,31 $rawdir/geuvadis.E-GEUV-1.sdrf.txt | grep _1.fastq.gz | cut -f3 > $rawdir/geuvadis.bamlist
checkIt $rawdir/geuvadis.bamlist

# list of individual ids
logIt Get list of invidual ids
cut -f1,29,31 $rawdir/geuvadis.E-GEUV-1.sdrf.txt | grep _1.fastq.gz | cut -f1 > $rawdir/geuvadis.ind
checkIt $rawdir/geuvadis.ind

## get the list of individual ids in each population
## all Geuvadis samples are present in 1000 genomes
logIt Get the list of individual ids in each population
for pop in CEU FIN GBR TSI YRI
do
cut -f1,29,31,33 $rawdir/geuvadis.E-GEUV-1.sdrf.txt | grep _1.fastq.gz | grep $pop | cut -f1 > $rawdir/geuvadis.$pop.ind
checkIt $rawdir/geuvadis.$pop.ind
grep -f $rawdir/geuvadis.$pop.ind $rawdir/kg.inds > $rawdir/geuvadis.$pop.kG.ind
checkIt $rawdir/geuvadis.$pop.kG.ind
done

exit 0

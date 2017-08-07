#!/bin/bash

#Define paths
progdir="$(dirname "$PWD")"/2_prog
LOG_FILE=$PWD/preq.log
touch $LOG_FILE

#Set urls
bigWigURL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
htslibURL=https://github.com/samtools/htslib.git
angsdURL=http://popgen.dk/software/download/angsd/angsd0.563.tar.gz
vcftoolsURL=https://github.com/vcftools/vcftools.git
bigWigURL=http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign50mer.bigWig

#function that logs everything it is given and prints it to the screen as well
function logIt {
    echo "[$(date -u)]: $*"
    echo "[$(date -u)]: $*" >> $LOG_FILE
}

#function that checks if given file is 0 bytes or doesn't exist
function checkIt {
  if [ -s "$1" ]
  then
    logIt $1 looks good! Moving on
  else
    logIt $1 does not look right, exiting
    exit 1
  fi
}

#Program starts from here
logIt Starting
logIt Going into program directory
cd $progdir

#!! Makes BigWig directory if it doesn't exist.
logIt Making BigWig directory if it doesnt exist
mkdir -p BigWig

#!! Open the BigWig directory
logIt Heading into BigWig directory
cd BigWig

#!! Download the BigWig software in the new directory
logIt Downloading bigWigToBedGraph from $bigWigURL
wget $bigWigURL
checkIt bigWigToBedGraph

logIt Making bigWigToBedGraph executable
chmod +x bigWigToBedGraph

#!! Back into original directory
logIt Heading back to original directory
cd ..

#!! Install autoconf and zlibgen, its a prerequisite for other softwares
logIt Installing autoconf, zlib, git and tabix
sudo apt install -y autoconf zlib1g-dev libbz2-dev liblzma-dev git tabix


#!! install htslib
logIt Git cloning htslib from $htslibURL
git clone $htslibURL
cd htslib

logIt Compiling htslib
autoheader
autoconf
./configure
make
sudo make install
logIt htslib compiled and installed

cd ..

#!! install angsd
logIt Downloading angsd from $angsdURL
wget $angsdURL
logIt Extracting angsd
checkIt angsd0.563.tar.gz
tar -xf angsd0.563.tar.gz
cd angsd0.563
sudo make HTSSRC=../htslib
export angsd=$PWD/angsd
cd ..
logIt angsd has been installed

#!! Compile the 3 softwares
echo $PWD

logIt Checking if all the cpp files are here
checkIt VCFmergeGTF3.cpp
checkIt ieatgorV2.cpp
checkIt getliners.cpp

logIt Compiling all of the programs
g++ -O3 -o VCFmergeGTF3 VCFmergeGTF3.cpp -lz
g++ -O3 -o ieatgor ieatgorV2.cpp -lz
g++ -O3 -o getliners getliners.cpp -lz

#!! install vcftools
logIt Downloading vcftools from $vcftoolsURL
git clone $vcftoolsURL
cd vcftools
export PERL5LIB=$PWD/src/perl/
./autogen.sh
./configure
logIt Compiling vcftools
make
sudo make install

#############
## 1. Convert Mapability file to bed format and create file to filter on mapability
#############
# This data is to develop a set of filters so that we can exclude regions of the
# genome that are difficult to map to.
#
# check to get the correct version of bigWig for your system
#
# Note: the mapability file needs to already have been downloaded an put in the BigWig directory

#!! Go in the new directory
cd $progdir/BigWig/

# convert bigwig to bed format:
#!! Download the bigWig file
logIt Downloading bigwig file from $bigWigURL
wget $bigWigURL

logIt bigWig file downloaded, converting mapability file to bed format and create file to filter on mapability
$progdir/BigWig/bigWigToBedGraph wgEncodeCrgMapabilityAlign50mer.bigWig wgEncodeCrgMapabilityAlign50mer.bed
checkIt wgEncodeCrgMapabilityAlign50mer.bed

# parse bed file to target file:
# Depending on the the system used (linux, mac, PC) you may need to change sed 's/\t/:/' to sed 's/	/:/'
#egrep "\s1$" wgEncodeCrgMapabilityAlign50mer.bed | sed 's/\t/:/' | sed 's/\t/-/' | cut -f1 > wgEncodeCrgMapabilityAlign50mer.target
logIt Parse bed file to target file
egrep "\s1$" wgEncodeCrgMapabilityAlign50mer.bed | sed 's/	/:/' | sed 's/	/-/' | cut -f1 > wgEncodeCrgMapabilityAlign50mer.target
checkIt wgEncodeCrgMapabilityAlign50mer.target

mv wgEncodeCrgMapabilityAlign50mer.target $rawdir

# parse target file to individual chromosome files
logIt Parsing target files to individual chromosomes
TRGTfile=$rawdir/wgEncodeCrgMapabilityAlign50mer.target
for chr in {1..22}
do
pattern="^chr$chr:"
egrep -e $pattern $TRGTfile > $rawdir/chr$chr.50mer.target
checkIt $rawdir/chr$chr.50mer.target
wc -l $rawdir/chr$chr.50mer.target
done

pattern="^chrX"
egrep -e $pattern $TRGTfile > $rawdir/chrX.50mer.target
checkIt $rawdir/chrX.50mer.target
wc -l $rawdir/chrX.50mer.target

pattern="^chrY"
egrep -e $pattern $TRGTfile > $rawdir/chrY.50mer.target
checkIt $rawdir/chrY.50mer.target
wc -l $rawdir/chrY.50mer.target

logIt Installed all pre requesites!

exit 0

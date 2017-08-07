#!/bin/bash

#############
# What: This bash script will call another script and together will conduct filtering and
#		parsing of the genotype and RNAseq data.
#
# 		For all Individuals, chromosomes, and population calls script 3_makeData.sh
#############

#this function checks if a location is a url or local file. If local file it makes sure it exists and if link it downloads it. It returns aboslute path to file in both cases.
function locationCheck {
  regex='(https?|ftp|file)://[-A-Za-z0-9\+&@#/%?=~_|!:,.;]*[-A-Za-z0-9\+&@#/%=~_|]'
  if [[ $1 =~ $regex ]]
  then
      echo "Downloading link"
      wget $1
      result="$PWD${string##*/}"
  else
      checkIt $1
      result=$1
  fi
}

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

#online location of Geuvadis data
location=http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/

#define directories
progdir="$(dirname "$PWD")"/2_prog
rawdir="$(dirname "$PWD")"/3_raw
outdir="$(dirname "$PWD")"/4_data

#define input filtering files
INFOfileP=$rawdir/geuvadis.E-GEUV-1.sdrf.txt
KGINDSP=$rawdir/kg.inds
GTFfileP=$rawdir/gencode.protein_coding.genes.v19.gtf
TRGTfileP=$rawdir/wgEncodeCrgMapabilityAlign50mer.target
pop=CEU #pop= population, Central European Utah
chr=13

#location of BAM files
BAMSfileP=$rawdir/geuvadis.bamlist

#if file path is local then check file else download it and then check it
locationCheck INFOfileP
INFOfile=$result
checkIt INFOfile

locationCheck KGINDSP
KGINDS=$result
checkIt KGINDS

locationCheck GTFfileP
GTFfile=$result
checkIt GTFfile

locationCheck TRGTfileP
TRGTfile=$result
checkIt TRGTfile

locationCheck BAMSfileP
BAMSfile=$result
checkIt BAMSfile

### PROGRAMS
VCFmergeGTF="$(dirname "$PWD")"/2_prog/VCFmergeGTF3
GETLINERS="$(dirname "$PWD")"/2_prog/getliners
IEATGOR="$(dirname "$PWD")"/2_prog/ieatgor
angsd=$PWD/temp/angsdnew/angsd

#load modules needed for this analysis
#!! Not used here angsd=$progdir/angsd/0.563
#!! Not used with new installation of vcftools vcftools=$progdir/vcftools_0.1.12b/bin/vcftools

#!! remove first line of ind file
tail -n +2 "$KGINDS" > "$KGINDS.tmp" && mv "$KGINDS.tmp" "$KGINDS"

#for each autosome
#!! Only done for one autosome so not needed
#for chr in {1..22}
#do

vcftools --gzvcf $rawdir/ALL.chr$chr.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz --keep $rawdir/geuvadis.$pop.kG.ind --remove-filtered-all --remove-indels --mac 2 --max-alleles 2 --recode --recode-INFO-all --out $rawdir/kg.chr$chr.$pop.mac2

# vcftools						calls the vcftools modules
# gzvcf							tells the package that the vcf file is zipped
# --keep						keep only these individuals
# geuvadis.$pop.kG.in			list of individual sample IDs that I wish to look at, in my case the CEU population (Central European Utah)
# --remove-filtered-all 		Removes all sites with a FILTER flag other than PASS.
# --remove-indels 				Include or exclude sites that contain an indel. For these options "indel" means any variant that alters the length of the REF allele.
# --mac 2 						Include only sites with Minor Allele Count greater than or equal to the "--mac" value and less than or equal to the "--max-mac" value.
# --max-alleles 2 				One of these options may be used without the other. Allele count is simply the number of times that allele appears over all individuals at that site.
# --recode 						output
# --recode-INFO-all 			These options can be used with the above recode options to define an INFO key name to keep in the output file. This option can be used multiple times to keep more of the INFO fields.
# --out 						define output file
# kg.chr22.CEU.mac2				name of my output file

#	done

    #for each individual
while read ind
do
  echo "Results in directory:" $outdir
  echo "Population: " $pop
  echo "Chromosome: " $chr
  echo "Individual: " $ind

  # merge VCF and GTF files for each individual
  # note that VCFmergeGTF will only work on a phased VCF file, as it is.
  logIt merge VCF and GTF files for each individual
  echo $ind > $outdir/tmp.ind
  $VCFmergeGTF $rawdir/gencode.chr$chr.gore $rawdir/kg.chr$chr.$pop.mac2.recode.vcf $outdir/tmp.ind $outdir/tmp.merged.gz

  # keep only heterozygous sites
  logIt keep only heterozygous sites
  zcat $outdir/tmp.merged.gz | awk '$8 != 0' |  awk '$8 != 3' | sed '1d' > $outdir/tmp.het
  echo "Number of heterozygous protein coding SNPs for ind" $ind
  wc -l $outdir/tmp.het

  # list heterozygous positions
  logIt list heterozygous positions
  cut -f1,2 $outdir/tmp.het | sed 's/^/chr/' > $outdir/tmp.positions

  ### READING GEUVADIS DATA
  logIt READING GEUVADIS DATA
  grep $ind $BAMSfile > $outdir/tmp.bam
  echo -n $location | cat - $outdir/tmp.bam > $outdir/tmp.bamlist
  #!! fixed path
  wc -l $outdir/tmp.bamlist

  ### GET COUNTS, MAPQ FILTER, REMOVES DUPLICATES SEE LAPPALAINEN SUPPLEMENT P 47
  logIt GET COUNTS, MAPQ FILTER, REMOVES DUPLICATES SEE LAPPALAINEN SUPPLEMENT P 47
  $angsd -bam $outdir/tmp.bamlist -r chr$chr: -doCounts 4 -out $outdir/tmp.counts.protein_coding.snps -dumpCounts 4 -doQsdist 1 -nThreads 1 -minMapQ 40

  ### COMBINING COUNTS AND POSITIONS FROM ANGSD
  logIt COMBINING COUNTS AND POSITIONS FROM ANGSD
  gunzip $outdir/tmp.counts.protein_coding.snps.pos.gz
  gunzip $outdir/tmp.counts.protein_coding.snps.counts.gz
  paste $outdir/tmp.counts.protein_coding.snps.pos $outdir/tmp.counts.protein_coding.snps.counts | awk '$3 != 0' | sed '1d' > $outdir/tmp.counts.protein_coding.snps.chr$chr

  ### FILTERING GtfVcf DATA AND MERGING WITH COUNTS
  logIt FILTERING GtfVcf DATA AND MERGING WITH COUNTS
  cut -f2 $outdir/tmp.counts.protein_coding.snps.chr$chr > $outdir/tmp.keys
  # tmp.counts.protein_coding.snps.chr$chr is the chr, position, and counts
  # tmp.keys is a list of the positions

  logIt get liners running
  $GETLINERS -c 2 -k $outdir/tmp.keys -f $outdir/tmp.het > $outdir/tmp.covered.snps.info

  logIt paste
  paste $outdir/tmp.covered.snps.info $outdir/tmp.counts.protein_coding.snps.chr$chr | cut -f1-8,11-15 | sed 's/^/chr/' >  $outdir/tmp.data

  ### FILTER FOR ALIGNABILITY

  logIt FILTER FOR ALIGNABILITY
  $IEATGOR $rawdir/chr$chr.50mer.target $outdir/tmp.data > $outdir/chr$chr.ind$ind.data

  ### REMOVE INTERMEDIATE FILES
#rm $outdir/tmp.*

done <$KGINDS

#!! Wrong filename?
#rm $rawdir/kg.chr13.$pop.mac2.recode.vcf

#done

#compile stuff
#g++ -O3 -o VCFmergeGTF3 VCFmergeGTF3.cpp -lz
#g++ -O3 -o ieatgor ieatgorV2.cpp -lz
#g++ -O3 -o getliners getliners.cpp -lz

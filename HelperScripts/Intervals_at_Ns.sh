#PBS -S /bin/bash
#PBS -q batch
#PBS -N Intervals_at_Ns
#PBS -l nodes=1:ppn=6
#PBS -l walltime=4:00:00
#PBS -l mem=22gb
#PBS -M 
#PBS -m abe

### A script to break up WGS data into Intervals at N-strings
### E. Dittmar February 2020

### The following script takes in as input:
### 1.) Reference FASTA
### 2.) List of genomic regions to break into smaller intervals (e.g chromosome names)

### And Outputs:
### 1.) A new FASTA file with only the sequence specified by the input list + corresponding index and dictionary files
### 2.) Both interval and BED files with the interval coordinates of N strings above the minimum number specified
###	3.) A genome file for use with Bedtools functions
### 4.) A final .bed file with intervals that are of the approximate size specified and broken up at strings of N's. 
      # This file can be put directly into the CUSTOM_INTERVALS variable of your config to run Genomics_DB_Import and Genotype_GVCFs

####################################
########### DEPENDENCIES ###########
####################################

### Requires GATK (including jar filepath), Picard (including jar filepath), and Bedtools
### Replace with correct installation paths

#   What is the full file path for the GATK jar file?
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
GATK_JAR=/usr/local/apps/eb/GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8/gatk

#   What is the full file path for the Picard jar file?
module load picard/2.16.0-Java-1.8.0_144
PICARD_JAR=/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar

module load BEDTools/2.29.2-GCC-8.2.0-2.31.1
#BEDTOOLS=
#export PATH=${BEDTOOLS}:${PATH}

####################################
##### USER-SUPPLIED VARIABLES ######
####################################

### Include full filepaths for all files

# Reference FASTA file
REF_FASTA=

# List of Genomic Regions to make intervals from
# Names must match names in REF_FASTA
INPUT_INTERVALS=

# Where do you want output files?
OUTPUTDIR=

#	Minimum number of "N's" which will be used for scattering intervals
NUM_Ns_TOL=500

#	Approximately how large do you want intervals to be (in bp)?
INT_SIZE=10000000

#	Specify a temporary directory
TEMP_DIR=

# Name for final output .bed file (no .bed extension here)
OUT_NAME=


####################################
## CODE FOR INTERVALS FILE AT N's ##
####################################

set -o pipefail

#1. Use GATK to find regions with N's

# First, make new FASTA file with the subset of genomic regions specified in INPUT_INTERVALS
# This automatically creates index and dictionary files which are also required by "ScatterIntervalsByNs"
echo "Writing FASTA files for input intervals"
 gatk FastaReferenceMaker \
   -R "$REF_FASTA" \
   -O "$OUTPUTDIR/Input_Intervals.fasta" \
   -L "${INPUT_INTERVALS}"

# Then make intervals file showing N-strings + remaining "ACGT" intervals
echo "Finding genomic coordinates for N-strings greater than ${NUM_Ns_TOL}"
java -jar ${PICARD_JAR} ScatterIntervalsByNs \
    R="${OUTPUTDIR}/Input_Intervals.fasta" \
    OT=BOTH \
    OUTPUT="${OUTPUTDIR}/N_${NUM_Ns_TOL}_sep.interval_list" \
    N="${NUM_Ns_TOL}"

# Then convert to BED file format
echo "Converting N-sep list to BED format"
java -jar ${PICARD_JAR} IntervalListToBed \
	I="${OUTPUTDIR}/N_${NUM_Ns_TOL}_sep.interval_list" \
	O="${OUTPUTDIR}/N_${NUM_Ns_TOL}_sep_intervals.bed" \
	TMP_DIR="${TEMP_DIR}"


#2. Make a Genome file

echo "Making Genome File"
awk -v OFS='\t' {'print $1,$2'} "${OUTPUTDIR}/Input_Intervals.fasta.fai" > "${OUTPUTDIR}/Input_Intervals_genomeFile.txt"


#3. Make Intervals file with smaller-sized intervals for each chromosome/region broken up at N's

echo "Breaking up genome at N-strings into regions that are approximately ${INT_SIZE} bp in size"

bedtools makewindows \
-g "${OUTPUTDIR}/Input_Intervals_genomeFile.txt" \
-w $INT_SIZE | \

bedtools intersect -wo -sorted -F 0.5 \
-a stdin \
-b "${OUTPUTDIR}/N_${NUM_Ns_TOL}_sep_intervals.bed" | \

bedtools groupby -i stdin \
-g 1,2 -c 5,6 -o min,max | \

awk -v OFS='\t' {'print $1,$3,$4'} > "${OUTPUTDIR}/${OUT_NAME}.bed"

# Replace chromosome names with their original names from the "Input_Intervals" file
# For some reason, GATK's FASTA reference maker changes the chromosome names to numbers based on their order
# Here the names are changed back

for i in $(cat $INPUT_INTERVALS); do
  oldvar=$(sed -n "/${i}/=" $INPUT_INTERVALS)
  echo "Replacing chromosome name $oldvar to $i"
  sed -i 's/'^"${oldvar}"'\b/'"${i}"'/g' "${OUTPUTDIR}/${OUT_NAME}.bed"
done


### A script to break up WGS data into Intervals at N-strings
### E. Dittmar February 2020

### The following script takes in as input:
### 1.) Reference FASTA
### 2.) List of genomic regions to use (e.g chromosome names)

### And Outputs:
### 1.) A new FASTA file with only the sequence specified by the input list + corresponding index and dictionary files
### 2.) Both interval and BED files with the interval coordinates of N strings above the minimum number specified
###	3.) A genome file for use with Bedtools functions
### 4.) A final .bed file with intervals that are of the approximate size specified and broken up at strings of N's

####################################
########### DEPENDENCIES ###########
####################################

### Requires GATK (including jar filepath), Picard (including jar filepath), and Bedtools

#   What is the full file path for the GATK jar file?
module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
GATK_JAR=/usr/local/apps/eb/GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8/gatk

#   What is the full file path for the Picard jar file?
module load picard/2.16.0-Java-1.8.0_144
PICARD_JAR=/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar

module load BEDTools/2.29.2-GCC-8.2.0-2.31.1


####################################
##### USER-SUPPLIED VARIABLES ######
####################################

# Reference FASTA sequence
REF_FASTA=/scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta

# List of Genomic Regions to make intervals from
# Names must match names in REF_FASTA
INPUT_INTERVALS=/scratch/eld72413/SAM_seq/results2/VCF_results/intervals2.list

# Where do you want output files?
OUTPUTDIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/ScriptTest/New

#	Minimum number of "N's" which will be used for scattering intervals
NUM_Ns_TOL=3000

#	Approximately how large do you want intervals to be?
INT_SIZE=20000000

#	Specify a temporary directory
TEMP_DIR=/scratch/eld72413/SAM_seq/results2/VCF_results_new/N_Intervals/Int_Test

# Name for final output .bed file (no .bed extension here)
OUT_NAME=N_List_test


####################################
## CODE FOR INTERVALS FILE AT N's ##
####################################

set -o pipefail

#1. Use GATK to find regions with N's

# First, make new FASTA file with interval subset
# This automatically creates an index and dictionary which is required by "ScatterIntervalsByNs"
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


#3. Make Intervals file with regularly sized intervals for each chromosome/region

echo "Breaking up genome over N-strings that are approximately ${INT_SIZE} in size"

bedtools makewindows \
-g "${OUTPUTDIR}/Input_Intervals_genomeFile.txt" \
-w $INT_SIZE | \

bedtools intersect -wo -sorted -F 0.5 \
-a stdin \
-b "${OUTPUTDIR}/N_sep_intervals.bed" | \

bedtools groupby -i stdin \
-g 1,2 -c 5,6 -o min,max | \

awk -v OFS='\t' {'print $1,$3,$4'} > "${OUTPUTDIR}/${OUT_NAME}.bed"

# Replace chromosome names with their original names from the "Input_Intervals" file
# For some reason, FASTA reference maker changes the chromosome names to numbers based on the order of the intervals
# Here the names are changed back

for i in $(cat $INPUT_INTERVALS); do
  oldvar=$(sed -n "/${i}/=" $INPUT_INTERVALS)
  echo "Replacing chromosome name $oldvar to $i"
  sed -i 's/'^"${oldvar}"'\b/'"${i}"'/g' "${OUTPUTDIR}/${OUT_NAME}.bed"
done


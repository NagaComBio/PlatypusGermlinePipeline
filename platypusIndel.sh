#!/bin/bash
# Reading in the source file
# Command line inputs
source $CONFIG_FILE
umask 007

###############################################################################
# Preparing the Input files

# Reading for the Somatic projects, control.BAM and tumor.BA
#
if [[ $CANCERorGERMLINE == 1 ]]
  then
  CONTROL_BAM=$ANALYSIS_FOLDER/${PID}/$ALIGNMENT_FOLDER/${BAM_PREFIX_CONTROL}_${PID}_${BAM_SUFIX}
  TUMOR_BAM=$ANALYSIS_FOLDER/${PID}/$ALIGNMENT_FOLDER/${BAM_PREFIX_TUMOR}_${PID}_${BAM_SUFIX}
  BAM=$CONTROL_BAM,$TUMOR_BAM

# List of samples, first pid will be the target for the output
#
elif [[ $CANCERorGERMLINE == 3 ]] 
  then
  allBAM=()
  arrayPIDS=($(echo $PIDs | tr "#" " "))
  echo $PIDS
  for pid in ${arrayPIDS[@]}
  do
    echo $pid
    PID_BAM=$BAM_DIR/$pid/$ALIGNMENT_FOLDER/${BAM_PREFIX_GERMLINE}_${pid}_${BAM_SUFIX}
    allBAM+=(${PID_BAM})
  done
  BAM=$(printf "%s," "${allBAM[@]}")
  BAM=${BAM/%,/}
fi

# Output file

OUTPUT_FILE=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${RAWFILE_PREFIX}${PID}.vcf  

###############################################################################
# Library file for recent Platypus version
export C_INCLUDE_PATH=$DATA_SHARE_PATH/analysis/htslib/lib/include/
export LIBRARY_PATH=$DATA_SHARE_PATH/analysis/htslib/lib/lib/
export LD_LIBRARY_PATH=$DATA_SHARE_PATH/analysis/htslib/lib/lib/

module load htslib/1.3.2
module load bcftools/1.9
module load Platypus/0.8.1.1 
###############################################################################
# The actual Playpus call
if [[ $DEFAULT == 0 ]] 
then
Platypus.py callVariants \
	--refFile=$REF \
	--output=$OUTPUT_FILE \
	--bamFiles=$BAM \
	--nCPU=10 \
	--genIndels=1 \
	--logFileName=${OUTPUT_FILE}.log.txt \
	--genSNPs=1 \
	--minFlank=0 \
	--verbosity=0
fi

###############################################################################
## Separting the SNVs and Indels

INPUT_FILE=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${RAWFILE_PREFIX}${PID}.vcf
OUTPUT_SNV=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${SNVFILE_PREFIX}${PID}_raw.vcf.gz
OUTPUT_INDEL=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${INDELFILE_PREFIX}${PID}_raw.vcf.gz

# Separate multi-alleic sites - Have to check again the proper funcationality
bcftools norm -m-both -o ${INPUT_FILE}.temp1 $INPUT_FILE
bcftools norm -f $REF -o ${INPUT_FILE}.temp2 $INPUT_FILE.temp1

# Remove old files if present
# SNVs 
if [ -e $OUTPUT_SNV ] 
then 
  rm $OUTPUT_SNV
  rm $OUTPUT_SNV.tbi
fi
# INDELs
if [ -e $OUTPUT_INDEL ]
then
  rm $OUTPUT_INDEL
  rm $OUTPUT_INDEL.tbi
fi

# Separate snvs and indels based in the length of the REF or ALT
(grep "#" ${INPUT_FILE}.temp2 ; cat ${INPUT_FILE}.temp2 | grep -v "#" | awk -F '\t' '{if(length($4) == length($5) && length($4)==1 ){print $0}}') | ${BGZIP_BIN} > $OUTPUT_SNV && ${TABIX_BIN} -p vcf $OUTPUT_SNV
(grep "#" ${INPUT_FILE}.temp2 ; cat ${INPUT_FILE}.temp2 | grep -v "#" | awk -F '\t' '{if(length($4) > 1 || length($5) > 1 ){print $0}}') | ${BGZIP_BIN} > $OUTPUT_INDEL && ${TABIX_BIN} -p vcf $OUTPUT_INDEL

# bgzip raw file 
bgzip -f $OUTPUT_FILE
tabix -p vcf -f ${OUTPUT_FILE}.gz

# Remove temp files
rm $INPUT_FILE.temp1 $INPUT_FILE.temp2 

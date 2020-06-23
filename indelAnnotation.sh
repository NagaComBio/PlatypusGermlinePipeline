#!/bin/bash

source $CONFIG_FILE
module load htslib/1.3.2
module load perl/5.20.2

### first check for existence of BAM files and their indexes!
#testing=FILES_TO_EVALUATE
#type=eval_job
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nTest for files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1
#[[ ${ok} == 0 ]] && echo -e "\nEvaluation of bam files for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 1

### Check and create output file names
#testing=FILE_NAMES
#type=create
#source ${TOOL_CHECK_FILES}
#[[ $? != 0 ]] && echo -e "\nCreation of file names for PID: ${PID} had non zero exit status, exiting pipeline\n\n" && exit 2

### Test availability of previous files
#[[ ! -f ${VCF_RAW} ]] && echo -e "The raw file: ${VCF_RAW} was not found, exiting..." && exit 2

### Remove previous downstream files
#rm ${VCF_FINAL} ${VCF_SOMATIC} ${VCF_SOMATIC_FUNCTIONAL}

### Create temporary filenames used in this script:
#TMP_FILE_PREFIX=${RESULTS_PER_PIDS_DIR}/${PID}/${RESULTS_SUB_DIR}/${INDELFILE_PREFIX}${PID}
TMP_FILE_PREFIX=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${INDELFILE_PREFIX}${PID}
VCF_TMP=${VCF_FINAL%.gz}.tmp
FOR_ANNOVAR=${TMP_FILE_PREFIX}.ForAnnovar.bed
ANNOVAR_SEGDUP=${TMP_FILE_PREFIX}.Annovar.segdup
ANNOVAR_CYTOBAND=${TMP_FILE_PREFIX}.Annovar.cytoband
ANNOVAR_VARIANT=${TMP_FILE_PREFIX}.ForAnnovar.bed.variant_function
ANNOVAR_VARIANT_EXON=${TMP_FILE_PREFIX}.ForAnnovar.bed.exonic_variant_function


set -x
set -o pipefail

eval VCF_RAW=$VCF_RAW
eval VCF_TMP=$VCF_TMP
eval VCF_FINAL=$VCF_FINAL

### Annotate with polymorphisms (dbSNP and 1K genomes) and prepare annovar input file
zcat ${VCF_RAW} | 
${TOOL_ANNOTATE_VCF} -a - -b "${DBSNP_INDEL}" --columnName=${DBSNP_COL} --reportMatchType  --bAdditionalColumn=2 --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${KGENOME_INDEL}" --columnName=${KGENOMES_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "$PHASE3_INDEL" --columnName=${PHASE3_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${ExAC_INDEL}" --columnName=${ExAC_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${gnomAD_INDEL}" --columnName=${gnomAD_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${LC2_INDEL}" --columnName=${LC2_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${LC3_INDEL}" --columnName=${LC3_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${LC4_INDEL}" --columnName=${LC4_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
${TOOL_ANNOTATE_VCF} -a - -b "${CLINVAR}" --columnName=${CLINVAR_COL}  --reportMatchType --bAdditionalColumn=2  --reportBFeatCoord --padding=${INDEL_ANNOTATION_PADDING} --minOverlapFraction=${INDEL_ANNOTATION_MINOVERLAPFRACTION} --maxBorderDistanceSum=${INDEL_ANNOTATION_MAXBORDERDISTANCESUM} --maxNrOfMatches=${INDEL_ANNOTATION_MAXNROFMATCHES} |
tee ${VCF_TMP} |
${TOOL_VCF_TO_ANNOVAR} ${CHR_PREFIX} ${CHR_SUFFIX} > ${FOR_ANNOVAR}.tmp


if [[ "$?" != 0 ]]
then
    echo "There was a non-zero exit code in the polymorphism annotation pipe; exiting..." 
    exit 3
else
    mv ${VCF_TMP} ${VCF_FINAL%.gz}
    mv ${FOR_ANNOVAR}.tmp ${FOR_ANNOVAR}
fi


###### Basic annotation with Annovar
# Gene annotation with annovar
${ANNOVAR} --buildver=${ANNOVAR_BUILDVER} ${ANNOVAR_DBTYPE} ${FOR_ANNOVAR} ${ANNOVAR_DBPATH}

# segdup annotation with annovar
${ANNOVAR} --buildver=${ANNOVAR_BUILDVER} -regionanno -dbtype segdup --outfile=${ANNOVAR_SEGDUP} ${FOR_ANNOVAR} ${ANNOVAR_DBPATH}
av_segdup=`ls ${ANNOVAR_SEGDUP}*genomicSuperDups`

# cytoband annotation with annovar
${ANNOVAR} --buildver=${ANNOVAR_BUILDVER} -regionanno -dbtype band --outfile=${ANNOVAR_CYTOBAND} ${FOR_ANNOVAR} ${ANNOVAR_DBPATH}
av_cytoband=`ls ${ANNOVAR_CYTOBAND}*cytoBand`

${TOOL_NEW_COLS_TO_VCF} --vcfFile=${VCF_FINAL%.gz} \
    --newColFile="${TOOL_PROCESS_ANNOVAR} ${ANNOVAR_VARIANT} ${ANNOVAR_VARIANT_EXON} |" \
    --newColHeader=${ANNOVAR_GENEANNO_COLS} --chrPrefix=${CHR_PREFIX} --chrSuffix=${CHR_SUFFIX} --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
${TOOL_NEW_COLS_TO_VCF} --vcfFile="-" --newColFile=${av_segdup} --newColHeader=${ANNOVAR_SEGDUP_COL} --chrPrefix=${CHR_PREFIX} --chrSuffix=${CHR_SUFFIX} \
    --reportColumns="1" --bChrPosEnd="2,7,8"  |
${TOOL_NEW_COLS_TO_VCF} --vcfFile="-" --newColFile=${av_cytoband} --newColHeader=${ANNOVAR_CYTOBAND_COL} --chrPrefix=${CHR_PREFIX} --chrSuffix=${CHR_SUFFIX} \
    --reportColumns="1" --bChrPosEnd="2,7,8" > ${VCF_TMP}


if [[ "$?" == 0 ]]
then
    mv ${VCF_TMP} ${VCF_FINAL%.gz}
else
    echo "There was a non-zero exit code in the Annovar annotation pipe; temp file ${VCF_TMP} not moved back" 
    exit 4
fi

#indel_reliability_pipe=`${PERL_BIN} ${TOOL_CREATE_PIPES} ${VCF_FINAL%.gz} ${CONFIG_FILE} ${TOOL_ANNOTATE_VCF} INDEL_RELIABILITY`

#if [[ "$?" != 0 ]] || [[ -z "${indel_reliability_pipe}" ]]; then echo "problem when generating INDEL_RELIABILITY pipe. Exiting..."; exit 5; fi

#eval ${indel_reliability_pipe} > ${VCF_TMP}

#if [[ "$?" == 0 ]]
#then
#    mv ${VCF_TMP} ${VCF_FINAL%.gz}
#else
#    echo "There was a non-zero exit code in the INDEL_RELIABILITY pipe; temp file ${VCF_TMP} not moved back" 
#    exit 5
#fi

################################## Confidence annotation ############################################

### Get the real names of the columns created by Platypus
#control_column=`${SAMTOOLS_BIN} view -H ${FILE_CONTROL_BAM} | grep -m 1 SM: | ${PERL_BIN} -ne 'chomp;$_=~m/SM:(.+)\w/;print "$1\n";'`
#tumor_column=`${SAMTOOLS_BIN} view -H ${FILE_TUMOR_BAM} | grep -m 1 SM: | ${PERL_BIN} -ne 'chomp;$_=~m/SM:(.+)\w/;print "$1\n";'`

#$PERL_BIN $TOOL_PLATYPUS_CONF_ANNO --fileName=${VCF_FINAL%.gz} --controlColName=${control_column} --tumorColName=${tumor_column} ${CONFIDENCE_OPTS} > ${VCF_TMP}

#if [[ "$?" == 0 ]]
#then
#    mv ${VCF_TMP} ${VCF_FINAL%.gz}
#else
#    echo "There was a non-zero exit code in the confidence annotation; temp file ${VCF_TMP} not moved back" 
#    exit 6
#fi

### Remove temp files:
rm $ANALYSIS_FOLDER/${PID}/$SUBDIR/${INDELFILE_PREFIX}${PID}*Annovar*

${BGZIP_BIN} -f ${VCF_FINAL%.gz} && ${TABIX_BIN} -f -p vcf ${VCF_FINAL}

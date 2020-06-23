#!/bin/bash

source $CONFIG_FILE 
module load htslib/1.3.2 
module load perl/5.20.2

set -o pipefail
set -x

cd $ANALYSIS_FOLDER/${PID}/$SUBDIR

### Check if there are downstream files and delete them
if [[ -f ${SNVFILE_PREFIX}${PID}_allSNVdiagnosticsPlots.pdf ]]; then rm ${SNVFILE_PREFIX}${PID}_allSNVdiagnosticsPlots.pdf; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_germline_coding_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf ]]; then rm ${SNVFILE_PREFIX}${PID}_germline_coding_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_intermutation_distance_conf_${MIN_CONFIDENCE_SCORE}_to_10.png ]]; then rm ${SNVFILE_PREFIX}${PID}_intermutation_distance_conf_${MIN_CONFIDENCE_SCORE}_to_10.png; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_intermutation_distance_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt ]]; then rm ${SNVFILE_PREFIX}${PID}_intermutation_distance_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_MAF_conf_${MAFMINSCORE}_to_10.png ]]; then rm ${SNVFILE_PREFIX}${PID}_MAF_conf_${MAFMINSCORE}_to_10.png; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_perChromFreq_conf_${MIN_CONFIDENCE_SCORE}_to_10.png ]]; then rm ${SNVFILE_PREFIX}${PID}_perChromFreq_conf_${MIN_CONFIDENCE_SCORE}_to_10.png; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_purityEST.txt ]]; then rm ${SNVFILE_PREFIX}${PID}_purityEST.txt; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_snvs_with_context_conf_${MAFMINSCORE}_to_10.png ]]; then rm ${SNVFILE_PREFIX}${PID}_snvs_with_context_conf_${MAFMINSCORE}_to_10.png; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_somatic_coding_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf ]]; then rm ${SNVFILE_PREFIX}${PID}_somatic_coding_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_somatic_in_dbSNP_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt ]]; then rm ${SNVFILE_PREFIX}${PID}_somatic_in_dbSNP_conf_${MIN_CONFIDENCE_SCORE}_to_10.txt; fi
if [[ -f ${SNVFILE_PREFIX}${PID}_somatic_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf ]]; then rm ${SNVFILE_PREFIX}${PID}_somatic_snvs_conf_${MIN_CONFIDENCE_SCORE}_to_10.vcf; fi
if [[ -f ${SNVFILE_PREFIX}${PID}.vcf.gz ]]; then rm ${SNVFILE_PREFIX}${PID}.vcf.gz; fi
if [[ -f ${SNVFILE_PREFIX}${PID}.vcf.gz.tbi ]]; then rm ${SNVFILE_PREFIX}${PID}.vcf.gz.tbi; fi
annofile=$(ls ${SNVFILE_PREFIX}${PID}*Annovar* 2> /dev/null | wc -l)
if [[ "$annofile" != 0 ]]; then rm ${SNVFILE_PREFIX}${PID}*Annovar*; fi

if [[ "${REMOVE_INDEL_CALLS}" == 1 ]]
then
	rm ${INDELFILE_PREFIX}${PID}*.vcf
fi

### Check if we have all per chromosome snv files
### Or if the already concatenated file exists in which case there will be a warning that it was imposible to test if all chromosmes have been called !
zippedfile=${SNVFILE_PREFIX}${PID}_raw.vcf.gz
if [ ! -s $zippedfile ]
then
	okay=1
	allfiles=''

	for chridx in ${CHROMOSOME_INDICES[@]}
	do
	    chr=${CHR_PREFIX}${chridx}${CHR_SUFFIX}
	    perchrfile=${SNVFILE_PREFIX}${PID}.${chr}.vcf
	    if [ ! -f $perchrfile ]
	    then
		echo "Error: File $perchrfile does not exist."
		okay=0
	    else
		allfiles=${allfiles}' '$perchrfile
	    fi
	done
	if [[ "$okay" != 0 ]]
	then
		perl $TOOLS_DIR/headeredFileConcatenator.pl ${allfiles} | ${BGZIP_BIN} > ${zippedfile} && ${TABIX_BIN} -p vcf ${zippedfile} && rm ${allfiles}
	fi
else
	okay=1
	
	for chridx in ${CHROMOSOME_INDICES[@]}
	do
		chr=${CHR_PREFIX}${chridx}${CHR_SUFFIX}
		if ${TABIX_BIN} -l ${zippedfile} | grep -Fxq "${chr}"
		then
			echo "Chromosome ${chr} found in zipped file"
		else 
			echo "No entries found for chromosome ${chr} in the zipped file. Chromosome could be missing or there are no variants on this chromosme"
			okay=2
		fi
	done
fi
if [[ "$okay" == 0 ]]
then
    echo "There were missing per chromosome snv files. Exiting..."
    exit 2
elif [[ "$okay" == 2 ]]
then
	echo "There was at least one chomosome without variants in the zipped concatenated file. Review your config file and consider to rerun calling!"
fi

###### Annotate with polymorphisms (dbSNP and 1K genomes) and prepare annovar input file
### Do different things depending on if the file is already zipped or not...
zcat ${zippedfile} | 
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$DBSNP" --columnName=$DBSNP_COL --reportMatchType  --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$KGENOMESNP" --columnName=$KGENOMES_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$PHASE3_SNP" --columnName=$PHASE3_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$ExAC_SNP" --columnName=$ExAC_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$gnomAD_SNP" --columnName=$gnomAD_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$LC2_SNV" --columnName=$LC2_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$LC3_SNV" --columnName=$LC3_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$LC4_SNV" --columnName=$LC4_COL  --reportMatchType --bAdditionalColumn=2 |
perl $TOOLS_DIR/annotate_vcf.pl -a - -b "$CLINVAR" --columnName=$CLINVAR_COL --reportMatchType --bAdditionalColumn=2 |
tee ${SNVFILE_PREFIX}${PID}.vcf.tmp |
perl $TOOLS_DIR/vcf_to_annovar.pl $CHR_PREFIX $CHR_SUFFIX > ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed.tmp


if [[ "$?" != 0 ]]
then
    echo "There was a non-zero exit code in the polymorphism annotation pipe; exiting..." 
    exit 2
else
    mv ${SNVFILE_PREFIX}${PID}.vcf.tmp ${SNVFILE_PREFIX}${PID}.vcf
    mv ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed.tmp ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed
fi


###### Basic annotation with Annovar
# Gene annotation with annovar
$ANNOVAR --buildver=$ANNOVAR_BUILDVER ${ANNOVAR_DBTYPE} ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed $ANNOVAR_DBPATH

# segdup annotation with annovar
$ANNOVAR --buildver=$ANNOVAR_BUILDVER -regionanno -dbtype segdup --outfile=${SNVFILE_PREFIX}${PID}.Annovar.seqdup ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed $ANNOVAR_DBPATH
av_segdup=`ls ${SNVFILE_PREFIX}${PID}.Annovar.seqdup*genomicSuperDups`

# self-chain annotation with annovar
#$ANNOVAR --buildver=$ANNOVAR_BUILDVER -regionanno -dbtype bed --bedfile SelfChain.bed --outfile=${SNVFILE_PREFIX}${#PID}.Annovar.selfchain ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed $ANNOVAR_DBPATH
#av_selfchain=`ls ${SNVFILE_PREFIX}${PID}.Annovar.selfchain*bed`

# cytoband annotation with annovar
$ANNOVAR --buildver=$ANNOVAR_BUILDVER -regionanno -dbtype band --outfile=${SNVFILE_PREFIX}${PID}.Annovar.cytoband ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed $ANNOVAR_DBPATH
av_cytoband=`ls ${SNVFILE_PREFIX}${PID}.Annovar.cytoband*cytoBand`

perl $TOOLS_DIR/newCols2vcf.pl --vcfFile=${SNVFILE_PREFIX}${PID}.vcf \
    --newColFile="perl $TOOLS_DIR/processAnnovarOutput.pl ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed.variant_function ${SNVFILE_PREFIX}${PID}.ForAnnovar.bed.exonic_variant_function |" \
    --newColHeader=$ANNOVAR_GENEANNO_COLS --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX --reportColumns="3,4,5,6" --bChrPosEnd="0,1,2" |
perl $TOOLS_DIR/newCols2vcf.pl --vcfFile="-" --newColFile=$av_segdup --newColHeader=$ANNOVAR_SEGDUP_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
    --reportColumns="1" --bChrPosEnd="2,7,8"  |
# perl $TOOLS_DIR/newCols2vcf.pl --vcfFile="-" --newColFile=$av_selfchain --newColHeader=$ANNOVAR_SELFCHAIN_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
#    --reportColumns="1" --bChrPosEnd="2,7,8" | 
perl $TOOLS_DIR/newCols2vcf.pl --vcfFile="-" --newColFile=$av_cytoband --newColHeader=$ANNOVAR_CYTOBAND_COL --chrPrefix=$CHR_PREFIX --chrSuffix=$CHR_SUFFIX \
    --reportColumns="1" --bChrPosEnd="2,7,8" > ${SNVFILE_PREFIX}${PID}.vcf.tmp

if [[ "$?" == 0 ]]
then
    mv ${SNVFILE_PREFIX}${PID}.vcf.tmp ${SNVFILE_PREFIX}${PID}.vcf
else
    echo "There was a non-zero exit code in the Annovar annotation pipe; temp file ${SNVFILE_PREFIX}${PID}.vcf.tmp not moved back" 
    exit 2
fi

## commenting out this part
#snv_reliability_pipe=`perl $TOOLS_DIR/createpipes.pl ${SNVFILE_PREFIX}${PID}.vcf $CONFIG_FILE $TOOLS_DIR/annotate_vcf.pl SNV_RELIABILITY`

#if [[ "$?" != 0 ]] || [[ -z "$snv_reliability_pipe" ]]; then echo "problem when generating SNV_RELIABILITY pipe. Exiting..."; exit 2; fi

#eval $snv_reliability_pipe > ${SNVFILE_PREFIX}${PID}.vcf.tmp

#if [[ "$?" == 0 ]]
#then
#     mv ${SNVFILE_PREFIX}${PID}.vcf.tmp ${SNVFILE_PREFIX}${PID}.vcf
#else
#    echo "There was a non-zero exit code in the SNV_RELIABILITY pipe; temp file ${SNVFILE_PREFIX}${PID}.vcf.tmp not moved back" 
#    exit 2
#fi

# confidence annotation with/without germline
# For platypus Idel we dont run confidence for now
if [[ $RUN_CONFIDENCE == 1 ]]
then
  if [[ "$GERMLINE_AVAILABLE" == 1 ]]
  then
	# for somatic SNVs only: if there are overlapping reads (same name) at the position, count the base only once;
	# also remove counts for bases that are neither reference nor variant, and change DP4 field accordingly
	# this needs the BAM file to perform a pileup
	#sample=$TUMOR_LABEL
	pid=$PID
	#eval "tumorbam=$BAMFILE_BP"
	eval "tumorbamfullpath=$TUMOR_BAMFILE_FULLPATH_BP"           # $RESULTS_PER_PIDS_DIR/${pid}/$ALIGNMENT_SUBDIR/$tumorbam
	# and the basequal and mapqual used for mpileup-bcftools.
	# if not specified, set to defaults 30 and 13, respectively
	mapqual=`echo $MPILEUP_OPTS | perl -ne '($qual) = $_ =~ /\-q\s*(\d+)/;print $qual'`
	mapqual=${mapqual:-30}
	basequal=`echo $MPILEUP_OPTS | perl -ne '($qual) = $_ =~ /\-Q\s*(\d+)/;print $qual'`
	basequal=${basequal:-13}

	cat ${SNVFILE_PREFIX}${PID}.vcf | python $PIPELINE_DIR/filter_PEoverlap.py --alignmentFile=${tumorbamfullpath} --mapq=$mapqual --baseq=$basequal --qualityScore=phred | perl $PIPELINE_DIR/confidenceREAnnotation_SNVs.pl -i - ${CONFIDENCE_OPTS} > ${SNVFILE_PREFIX}${PID}.vcf.tmp
	if [[ "$?" == 0 ]]
	then
		if [ -f ${SNVFILE_PREFIX}${PID}.vcf.tmp ]; then mv ${SNVFILE_PREFIX}${PID}.vcf.tmp ${SNVFILE_PREFIX}${PID}.vcf; fi
		else
		echo "SNV confidenceAnnotation with germline pipe returned non-zero exit code; temp file ${SNVFILE_PREFIX}${PID}.vcf.tmp not moved back"
		exit 21
	fi
  else	# no germline information available
	perl $PIPELINE_DIR/confidenceAnnotation_SNVsNoGermline.pl ${SNVFILE_PREFIX}${PID}.vcf $CONFIG_FILE > ${SNVFILE_PREFIX}${PID}.vcf.tmp
	if [[ "$?" == 0 ]]
	then
		if [ -f ${SNVFILE_PREFIX}${PID}.vcf.tmp ]; then mv ${SNVFILE_PREFIX}${PID}.vcf.tmp ${SNVFILE_PREFIX}${PID}.vcf; fi
		else
		echo "confidenceAnnotation returned non-zero exit code; temp file ${SNVFILE_PREFIX}${PID}.vcf.tmp not moved back"
		exit 21
	fi
  fi
fi

${BGZIP_BIN} -f ${SNVFILE_PREFIX}${PID}.vcf && ${TABIX_BIN} -f -p vcf ${SNVFILE_PREFIX}${PID}.vcf.gz

annofile=$(ls ${SNVFILE_PREFIX}${PID}*Annovar* 2> /dev/null | wc -l)
if [[ "$annofile" != 0 ]]; then rm ${SNVFILE_PREFIX}${PID}*Annovar*; fi

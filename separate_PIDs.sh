#!/bin/bash

source $CONFIG_FILE

arrayPIDS=($(echo $PIDs | tr "#" " "))
#echo $PIDS
#DEFA_HEADER="#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT"
#DEEP_HEADER="DBSNP,1K_GENOMES,ANNOVAR_FUNCTION,GENE,EXONIC_CLASSIFICATION,ANNOVAR_TRANSCRIPTS,SEGDUP,CYTOBAND,MAPABILITY,HISEQDEPTH,SIMPLE_TANDEMREPEATS,REPEAT_MASKER,DUKE_EXCLUDED,DAC_BLACKLIST,SELFCHAIN,CpGislands,CgiMountains,Enhancers,miRNAs_snoRNAs,miRBase18,miRNAtargets,phastConsElem20bp,TFBScons,DNMT3_BS,COSMIC,COSMIC_indels"

ForSNVDeep_file=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${SNVFILE_PREFIX}${PID}.vcf.gz
ForIndelDeep_file=$ANALYSIS_FOLDER/${PID}/$SUBDIR/${INDELFILE_PREFIX}${PID}.vcf.gz
#AllVariants_file=$ANALYSIS_FOLDER/${PID}/$SUBDIR/AllVar_${PID}.vcf.gz

#(zcat $ForSNVDeep_file | head -n2000 | grep '#' | cut -f1-18 ; (zcat $ForSNVDeep_file | grep -v '#' | cut -f1-18 ; zcat $ForIndelDeep_file | grep -v '#' | cut -f1-18) | sort -V ) | bgzip -f > $AllVariants_file && tabix -f -p vcf $AllVariants_file

for pid in ${arrayPIDS[@]}
  do
    echo $pid
    if [[ ! -d $ANALYSIS_FOLDER/$pid/${SUBDIR}_Separate ]]
    then
      mkdir -p -m 770 $ANALYSIS_FOLDER/$pid/${SUBDIR}_Separate
    fi   
    snv_sep=$ANALYSIS_FOLDER/${pid}/${SUBDIR}_Separate/${SNVFILE_PREFIX}${pid}.vcf.gz
    indel_sep=$ANALYSIS_FOLDER/${pid}/${SUBDIR}_Separate/${INDELFILE_PREFIX}${pid}.vcf.gz
    #AllVariants_sep=$ANALYSIS_FOLDER/${pid}/${SUBDIR}_Separate/AllVar_${pid}.vcf.gz

    #perl $PIPELINE_DIR/pullColumnwithHeader.pl $DEFA_HEADER,${VCF_PREFIX_GERMLINE}${pid},$DEEP_HEADER $ForSNVDeep_file $snv_sep
    #perl $PIPELINE_DIR/pullColumnwithHeader.pl $DEFA_HEADER,${VCF_PREFIX_GERMLINE}${pid},$DEEP_HEADER $ForIndelDeep_file $indel_sep
    ln -s -f $ForSNVDeep_file $snv_sep
    ln -s -f $ForSNVDeep_file.tbi $snv_sep.tbi

    ln -s -f $ForIndelDeep_file $indel_sep
    ln -s -f $ForIndelDeep_file $indel_sep.tbi

    #ln -s -f $AllVariants_file $AllVariants_sep
    #ln -s -f $AllVariants_file.tbi $AllVariants_sep.tbi
  done


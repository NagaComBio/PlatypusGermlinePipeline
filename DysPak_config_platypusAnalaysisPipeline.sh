#!/bin/bash

# Config file for Platypus Pipeline
DATA_ROOT_PATH=ROOT_PATH
DATA_SHARE_PATH=ROOT_SHARE_PATH
PIPELINE_DIR=/home/$USER/projects/gitRepos/platypusgermlinepipeline
TOOLS_DIR=$PIPELINE_DIR

PROJECT_NAME=Dystonia_Pakistan
ANALYSIS_FOLDER=$DATA_ROOT_PATH/analysis/UniExome/UnknownGeneticDisorders/sequencing/exon_sequencing/results_per_pid
CLUSTER_EO=$ANALYSIS_FOLDER/../log_folder
#BAM_DIR=$DATA_ROOT_PATH/project/UniExome/$PROJECT_NAME/sequencing/exon_sequencing/view-by-pid
BAM_DIR=$ANALYSIS_FOLDER

# Sample type to analysis
# Cancer == 1
# Germline == 2
# Family == 3
CANCERorGERMLINE=3

# Only applies to GermLine sample
MULTICALL=0

# Input files
#ALIGNMENT_FOLDER=/blood/paired/merged-alignment
ALIGNMENT_FOLDER=alignment

RAWFILE_PREFIX=rawCalls_
INDELFILE_PREFIX=indel_
SNVFILE_PREFIX=snvs_

BAM_SUFIX=merged.bam.rmdup.bam

#For Germline sample
BAM_PREFIX_GERMLINE=blood
VCF_PREFIX_GERMLINE=sample_blood_

# Reference Genome 
REF=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/indexes/bwa/bwa06_1KGRef/hs37d5.fa

# Number of individual called together !! 
NUMBER_INDIVIDUALS=1

# Platypus running parameters.
DEFAULT=0

# subdir for platypus temp and output files
SUBDIR=platypus_familial

# Pipeline parts to run
RUN_PLATYPUS=1
RUN_INDELANNOTATION=1
RUN_SNVANNOTATION=1
RUN_SEPARATE_PIDs=1

############################################################################
### Used binaries
### Change here to use software not in $PATH or to explicitly use certain versions
SAMTOOLS_BIN=samtools
BCFTOOLS_BIN=bcftools
BGZIP_BIN=bgzip
TABIX_BIN=tabix

### Scripts used by pipeline with paths:
#<SCRIPTS
TOOL_PLATYPUS_INDEL_ANNOTATION=${PIPELINE_DIR}/indelAnnotation.sh
TOOL_ANNOTATE_VCF=${PIPELINE_DIR}/annotate_vcf.pl
TOOL_NEW_COLS_TO_VCF=${PIPELINE_DIR}/newCols2vcf.pl
TOOL_PROCESS_ANNOVAR=${PIPELINE_DIR}/processAnnovarOutput.pl
TOOL_VCF_TO_ANNOVAR=${PIPELINE_DIR}/vcf_to_annovar.pl
#>SCRIPTS

### File names
INDELFILE_PREFIX="indel_"
#<FILE_NAMES
VCF_RAW=${ANALYSIS_FOLDER}/'${PID}'/${SUBDIR}/${INDELFILE_PREFIX}'${PID}'_raw.vcf.gz
VCF_FINAL=${ANALYSIS_FOLDER}/'${PID}'/${SUBDIR}/${INDELFILE_PREFIX}'${PID}'.vcf.gz
#>FILE_NAMES

### reference genome and chromosomes to run
# Platypus runs on all chromosome,this is just a hack.
CHROMOSOME_INDICES=(1)

CHR_PREFIX='' # set to '' for samples with chr-less coordinates
CHR_SUFFIX='' # set to '.fa' for ill-named EMBL samples

### annovar
ANNOVAR=$DATA_SHARE_PATH/annovar/annovar_Sept2013/annotate_variation.pl
ANNOVAR_BUILDVER='hg19'
ANNOVAR_DBPATH='$DATA_SHARE_PATH/annovar/annovar_Nov2014/humandb/'
ANNOVAR_DBTYPE='-dbtype wgEncodeGencodeCompV19'

### Annotation (SNV and Indel)
KGENOMES_COL=1K_GENOMES
KGENOMES_ADDITIONAL_COLUMNS=5 # comma-separated list (no spaces!) of zero-based column numbers of columns from kgenomes vcf to be reported
DBSNP_COL=DBSNP
DBSNP_ADDITIONAL_COLUMNS=2 # comma-separated list (no spaces!) of zero-based column numbers of columns from dbsnp vcf to be reported
ANNOVAR_GENEANNO_COLS="ANNOVAR_FUNCTION,GENE,EXONIC_CLASSIFICATION,ANNOVAR_TRANSCRIPTS"
ANNOVAR_SEGDUP_COL=SEGDUP
#ANNOVAR_SELFCHAIN_COL=SELFCHAIN
ANNOVAR_CYTOBAND_COL=CYTOBAND
LC2_COL=LOCAL_CONTROL_2
LC3_COL=LOCAL_CONTROL_3
LC4_COL=LOCAL_CONTROL_4
PHASE3_COL=PHASE3_1K_GENOMES
ExAC_COL=gnomAD_EXOMES_v2.1
gnomAD_COL=gnomAD_GENOMES_v2.1

#<PIPE_CONFIG:SNV_RELIABILITY
MAPABILITY=$DATA_SHARE_PATH/annotation/hg19/wgEncodeCrgMapabilityAlign100mer.bedGraph.gz
HISEQDEPTH=$DATA_SHARE_PATH/annotation/hg19/HiSeqDepthTop10Pct.bed.gz
SIMPLE_TANDEMREPEATS=$DATA_SHARE_PATH/annotation/hg19/SimpleTandemRepeats.bed.gz:4
REPEAT_MASKER=$DATA_SHARE_PATH/annotation/hg19/RepeatMasker_UCSCdownload_20120514.txt_sorted.bed.gz
DUKE_EXCLUDED=$DATA_SHARE_PATH/annotation/hg19/DukeExcluded.bed.gz
DAC_BLACKLIST=$DATA_SHARE_PATH/annotation/hg19/DACBlacklist.bed.gz
SELFCHAIN=$DATA_SHARE_PATH/annotation/hg19/selfChain.bed.gz:4::--maxNrOfMatches=5
#>PIPE_CONFIG

#<PIPE_CONFIG:INDEL_RELIABILITY
MAPABILITY="$DATA_SHARE_PATH/annotation/hg19/wgEncodeCrgMapabilityAlign100mer.bedGraph.gz:::--breakPointMode --aEndOffset=1"
HISEQDEPTH=$DATA_SHARE_PATH/annotation/hg19/HiSeqDepthTop10Pct.bed.gz
SIMPLE_TANDEMREPEATS=$DATA_SHARE_PATH/annotation/hg19/SimpleTandemRepeats.bed.gz:4
REPEAT_MASKER=$DATA_SHARE_PATH/annotation/hg19/RepeatMasker_UCSCdownload_20120514.txt_sorted.bed.gz
DUKE_EXCLUDED=$DATA_SHARE_PATH/annotation/hg19/DukeExcluded.bed.gz
DAC_BLACKLIST=$DATA_SHARE_PATH/annotation/hg19/DACBlacklist.bed.gz
SELFCHAIN=$DATA_SHARE_PATH/annotation/hg19/selfChain.bed.gz:4::--maxNrOfMatches=5
#>PIPE_CONFIG

### confidence assignment
MAXDEPTH=150    # control coverage above which artefact region will be assumed. Default 150. Set accordingly high for targetted sequencing

### SNV annotation
DBSNP=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/databases/dbSNP/dbSNP_147/00-All.SNV.vcf.gz
KGENOMESNP=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/databases/1000genomes/ALL.wgs.phase1_integrated_calls.20101123.snps_chr.vcf.gz

### SNVs and Indel annotation
LC2_SNV=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/LocalControls/UniExome_playtpusMaster.SNVs.vcf.gz
LC2_INDEL=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/LocalControls/UniExome_playtpusMaster.Indels.vcf.gz

LC3_SNV=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/LocalControls/MMML_WGS_Platypus_minFlank_RAW.SNVs.vcf.gz
LC3_INDEL=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/LocalControls/MMML_WGS_Platypus_minFlank_RAW.Indels.vcf.gz

LC4_SNV=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/databases/LocalControls/pcawg_AF.All.CommonArtifacts.vcf.gz
LC4_INDEL=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/databases/LocalControls/pcawg_AF.All.CommonArtifacts.vcf.gz

PHASE3_SNP=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/1000Genome/Phase3/forNGSshare/ALL.chrALL.SNP.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
PHASE3_INDEL=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/1000Genome/Phase3/forNGSshare/ALL.chrALL.Indel.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

ExAC_SNP=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/v2.1/exomes/gnomad.exomes.r2.1.sites.SNV.vcf.gz
ExAC_INDEL=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/v2.1/exomes/gnomad.exomes.r2.1.sites.INDEL.vcf.gz

gnomAD_SNP=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/v2.1/genomes/gnomad.genomes.r2.1.sites.SNV.vcf.gz
gnomAD_INDEL=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/v2.1/genomes/gnomad.genomes.r2.1.sites.INDEL.vcf.gz

CLINVAR=$DATA_ROOT_PATH/analysis/UniExome/GermlineAnnotation/annotationInfo/clinVar/clinvar_20170530.vcf.gz
CLINVAR_COL=CLINVAR

### Indel annotation
DBSNP_INDEL=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/databases/dbSNP/dbSNP_147/00-All.INDEL.vcf.gz
KGENOME_INDEL=$DATA_SHARE_PATH/assemblies/hg19_GRCh37_1000genomes/databases/1000genomes/ALL.wgs.phase1_integrated_calls.20101123.indels_plain.vcf.gz
INDEL_ANNOTATION_PADDING=10
INDEL_ANNOTATION_MINOVERLAPFRACTION=0.7
INDEL_ANNOTATION_MAXBORDERDISTANCESUM=20

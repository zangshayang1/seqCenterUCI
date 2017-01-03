#!/bin/bash
#$ -N SeanTest
#$ -q cbcl-a64,cbcl-i24,bio,pub64,free64
#$ -pe openmp 8-64
#$ -m beas


<<COMMENT
========================================================================== CONFIG
COMMENT

<<COMMENT
=========================================== set dir
COMMENT

PARENT_DIR=test_run_latest
RAW_DATA_DIR=raw_data
PRE_QC_DIR=preQC
POST_QC_DIR=postQC
TRIMMED_READS_DIR=trimmed_reads
ALIGNMENT_DIR=star_alignment_2pass
FEATURE_COUNTS_DIR=featureCounts_output


<<COMMENT
=========================================== data sourcing
COMMENT

WEBPAGE=dataSourcePage.htm
WEB_PYPARSER=GenomeDataSourcedTo.py
RAW_DATA_URLs=fastq_urls.txt
HASH_CODE_LIST=md5sum_codes.txt
DATA_EXTENSION=fastq


<<COMMENT
=========================================== fastqc settings
COMMENT

PREQC_MODE=
POSTQC_MODE=PEtrimed


<<COMMENT
=========================================== trimmomatic settings
COMMENT

MODE=PE


<<COMMENT
=========================================== STAR settings
COMMENT

MAXIMUM_READS=30000000
INDEXED_GENOME_PATH=/cbcl/szang/seqCenter/refGenome/Mus_musculus/Ensembl/NCBIM37/StarIndex



<<COMMENT
=========================================== featureCounts settings
COMMENT

ANNOTATION_PATH=/cbcl/szang/seqCenter/refGenome/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf
FEATURECOUNTS_LIST=./featureCountsList
FEATURECOUNTS_COMBINED=featureCountsCombined.txt


<<COMMENT
=========================================== environment settings
COMMENT

TERM=xterm-256color



<<COMMENT
========================================================================== EXE
COMMENT


<<COMMENT
=========================================== modules
COMMENT

module load anaconda/2-2.3.0

source ./RNAseq_module.sh



<<COMMENT
=========================================== directories
COMMENT

mkdir $PARENT_DIR
mkdir $PARENT_DIR/$RAW_DATA_DIR
mkdir $PARENT_DIR/$PRE_QC_DIR
mkdir $PARENT_DIR/$POST_QC_DIR
mkdir $PARENT_DIR/$TRIMMED_READS_DIR
mkdir $PARENT_DIR/$ALIGNMENT_DIR
mkdir $PARENT_DIR/$FEATURE_COUNTS_DIR


<<COMMENT
=========================================== download data
COMMENT

# $WEB_PYPARSER: $1 resource webpage file; $2 name the output list of urls; $3 name the output list of hashcode;
python $WEB_PYPARSER $WEBPAGE $RAW_DATA_URLs $HASH_CODE_LIST

# Download_Fastq: $1 input list of downloadable urls; $2 output dir;
Download_Fastq $RAW_DATA_URLs $PARENT_DIR/$RAW_DATA_DIR

# Check_md5sum: $1 input rawData.txt.gz dir; $2 original hashcode;
Check_md5sum $PARENT_DIR/$RAW_DATA_DIR $HASH_CODE_LIST

# Decompress_gzFiles: $1 input rawData.gz dir; 
Decompress_gzFiles $PARENT_DIR/$RAW_DATA_DIR

# Modify_Extensions: $1 input rawData.txt dir; $2 specify extension (fastq);
Modify_Extensions $PARENT_DIR/$RAW_DATA_DIR $DATA_EXTENSION

# Remove_PrNotRecog: $1 input rawData.fastq dir;
Remove_PrNotRecog $PARENT_DIR/$RAW_DATA_DIR

# FastQC: $1 input rawData.fastq dir; $2 output QC dir;
FastQC $PARENT_DIR/$RAW_DATA_DIR $PARENT_DIR/$PRE_QC_DIR $PREQC_MODE

# TrimReads: $1 SE || PE; $2 input rawData.fastq dir; $3 output trimedData.fastq dir
TrimReads $MODE $PARENT_DIR/$RAW_DATA_DIR $PARENT_DIR/$TRIMMED_READS_DIR

# Post-trimming QC: $1 input trimedData.fastq dir; $2 output dir; $3 specify searching pattern
FastQC $PARENT_DIR/$TRIMMED_READS_DIR $PARENT_DIR/$POST_QC_DIR $POSTQC_MODE

# RunStar: $1 input trimedData.fastq dir; $2 output aligned.bam dir; $3 indexed genome; $4 Maximum number of input reads
RunStar $PARENT_DIR/$TRIMMED_READS_DIR $PARENT_DIR/$ALIGNMENT_DIR $INDEXED_GENOME_PATH $MAXIMUM_READS

# FeatureCounts: $1 input Aligned_data dir; $2 output featureCounts dir; $3 annotation file
FeatureCounts $PARENT_DIR/$ALIGNMENT_DIR $PARENT_DIR/$FEATURE_COUNTS_DIR $ANNOTATION_PATH

# Combine separated featureCounts outcomes
python Combinator.py $FEATURECOUNTS_LIST $PARENT_DIR/$FEATURE_COUNTS_DIR/$FEATURECOUNTS_COMBINED

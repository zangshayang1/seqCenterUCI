#!/bin/bash
#$ -N SeanTest
#$ -q cbcl-a64,cbcl-i24,bio,pub64,free64
#$ -pe openmp 8-64
#$ -m beas


# hahahaha

##################
# 0. SET CONSTANT
##################

# env var

TERM=xterm-256color



# job root/sub sir

PARENT_DIR=test_run
RAW_DATA_DIR=raw_data
PRE_QC_DIR=preQC
POST_QC_DIR=postQC
TRIMMED_READS_DIR=trimmed_reads
ALIGNMENT_DIR=star_alignment_2pass
FEATURE_COUNTS_DIR=featureCounts_output



# crawled urls storage

RAW_DATA_URLs=raw_urls



# input pattern of downstream process [decompress, rename, qc, featureCounts]

COMPRESSED_RAW_DATA_PATTERN=*.txt.gz
DECOMPRESSED_RAW_DATA_PATTERN=*.txt
RAW_READS_PATTERN=*.fastq
STAR_MAPPED_SORTED_RESULTS_PATTERN=*.sortedByCoord.out.bam



# name featureCounts output file

FEATURE_COUNTS_OUTPUT=test_run_featureCounts



# additional labels 

TRIM_PAIR_LABEL=_paired
TRIM_UNPAIR_LABEL=_unpaired



# pre-running settings

PAIR_END=True
STAR_INDEXED=True
STAR_TWO_PASS=True
MAXIMUM_READS=15640877



# reference and package path

PYTHON_CRAWLER_PATH=./genome_spider.py
UNINDEXED_GENOME_PATH=refGenome/Mus_musculus/Ensembl/NCBIM37/Sequence/WholeGenomeFasta/genome.fa
INDEXED_GENOME_PATH=refGenome/Mus_musculus/Ensembl/NCBIM37/StarIndex
TRIMMOMATIC_PATH=/data/apps/trimmomatic/0.35/trimmomatic-0.35.jar
ANNOTATION_PATH=refGenome/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf


#################
# 1. load modules 
#################

module load anaconda/2-2.3.0
module load fastqc/0.11.2
module load java/1.8.0.51
module load STAR/2.5.1b-static
module load subread/1.5.0-p3


############
# 2. make dir
############

mkdir $PARENT_DIR
mkdir $PARENT_DIR/$RAW_DATA_DIR
mkdir $PARENT_DIR/$PRE_QC_DIR
mkdir $PARENT_DIR/$POST_QC_DIR
mkdir $PARENT_DIR/$TRIMMED_READS_DIR
mkdir $PARENT_DIR/$ALIGNMENT_DIR
mkdir $PARENT_DIR/$FEATURE_COUNTS_DIR


######################################
# 3. download dataset
# RAW_DATA_URLs is a List of data urls
######################################

python $PYTHON_CRAWLER_PATH $PARENT_DIR/$RAW_DATA_URLs


while read line
do
    wget -P $PARENT_DIR/$RAW_DATA_DIR $line   	
       # -P output_dir_prefix
       # the output name will follow whatever it was originally named
done < $PARENT_DIR/$RAW_DATA_URLs


#///////////////////////////
# need to add md5sum
#///////////////////////////

##########################
# 4. decompress and rename
##########################


for data in $PARENT_DIR/$RAW_DATA_DIR/$COMPRESSED_RAW_DATA_PATTERN
do
    gunzip $data
    # automatically replace .gz files with output .txt
done

for data in $PARENT_DIR/$RAW_DATA_DIR/$DECOMPRESSED_RAW_DATA_PATTERN
do 
    prefix_length=$((${#data}-4))   
    prefix=${data:0:$prefix_length}
    mv $data $prefix.fastq
done


#///////////////////////////
# separate output into independent dir
#///////////////////////////



################
# 5. pre trim QC
################

fastqc -o $PARENT_DIR/$PRE_QC_DIR $PARENT_DIR/$RAW_DATA_DIR/$RAW_READS_PATTERN 


#///////////////////////////
# sometimes trimming params can only be decided after you check out the QC results.
#///////////////////////////


#///////////////////////////
# separate output into independent dir
#///////////////////////////



#################
# 5. trim   # separate output into different dirs
#################

reads1_list=($PARENT_DIR/$RAW_DATA_DIR/*READ1$RAW_READS_PATTERN)
reads2_list=($PARENT_DIR/$RAW_DATA_DIR/*READ2$RAW_READS_PATTERN)

len1=${#reads1_list[@]}
len2=${#reads2_list[@]}

if [ $len1==$len2 ]; then
    echo "input size match! ready to loop into Trimmomatic..."
else
    echo "input sizes do not match! exiting..."
    exit
fi


for ((i = 0; i < $len1 && i < $len2; i++))
do
    # ALL The Following $readsin/out Stand For Both Path/File
    # Establish Unified Output Dir

    output_dir=$PARENT_DIR/$TRIMMED_READS_DIR
    
    echo "Trimmomatic New Round Started... Output Dir = $output_dir"
    
    readsin1=${reads1_list[i]}
    readsin2=${reads2_list[i]}
    
    # String Manipulation Fails When / Is Found In The Path
    # So Basename Taken
    base1=$(basename $readsin1)
    base2=$(basename $readsin2)
    
    prefix1_length=$((${#base1}-6))  #remove .fastq
    prefix2_length=$((${#base2}-6))
    prefix1=${base1:0:$prefix1_length}
    prefix2=${base2:0:$prefix2_length}

    # String Manu: $strVar1$strVar2 concatenates but $strVar1"_paired.fastq" will not concatenate without ""
    # load TRIM_PAIR/UNPAIR_LABEL in output files
    readsout1=$output_dir/$prefix1$TRIM_PAIR_LABEL".fastq"
    readsout2=$output_dir/$prefix1$TRIM_UNPAIR_LABEL".fastq"
    readsout3=$output_dir/$prefix2$TRIM_PAIR_LABEL".fastq"
    readsout4=$output_dir/$prefix2$TRIM_UNPAIR_LABEL".fastq"
    
    if [ PAIR_END ]; then
        java -jar $TRIMMOMATIC_PATH PE -threads 24 \
             -phred33 $readsin1 $readsin2 $readsout1 $readsout2 $readsout3 $readsout4 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    else
        echo "under dev..."
    fi
done

echo "Done trimming!"


#///////////////////////////
# separate output into independent dir
#///////////////////////////


####################
# 6. post-trim QC
####################

echo Ready to do post-QC...

fastqc -o $PARENT_DIR/$POST_QC_DIR $PARENT_DIR/$TRIMMED_READS_DIR/$RAW_READS_PATTERN




#///////////////////////////
# NOTE: consider technical duplicates before mapping
#///////////////////////////
#///////////////////////////
# separate output into independent dir
#///////////////////////////



#########################################
# 7. alignment - star, local var - label_length - needs to be set according to the name of input 
###################################################################################

# first exam if the genome.fa is STAR-indexed or not
# if not, the following builds the STAR-indexed genome

if [ STAR_INDEXED ]; then
    echo "STAR indexed genome already exists. Starting Mapping..."
else
    echo "Generating STAR genome indexing..."

    # --sjdbOverhang followed by (reads length - 1), reads length = 100 in our case
    STAR --runThreadN 24 \
         --runMode genomeGenerate \
         --genomeDir $INDEXED_GENOME_PATH \
         --genomeFastaFiles $UNINDEXED_GENOME_PATH \
         --sjdbGTFfile $ANNOTATION_PATH \
         --sjdbOverhang 99 
fi


# preparing paired lists of inputs - read1 & read2

reads1_list=($PARENT_DIR/$TRIMMED_READS_DIR/*READ1*$TRIM_PAIR_LABEL$RAW_READS_PATTERN)
reads2_list=($PARENT_DIR/$TRIMMED_READS_DIR/*READ2*$TRIM_PAIR_LABEL$RAW_READS_PATTERN)

len1=${#reads1_list[@]}
len2=${#reads2_list[@]}

if [ $len1==$len2 ]; then
    echo input size match! ready to loop into STAR...
else
    echo input sizes do not match! exiting...
    exit
fi


for ((i = 0; i < $len1 && i < $len2; i++))
do
    echo "Mapping New Round Started... Output Dir =$PARENT_DIR/$ALIGNMENT_DIR"

    label_length=18   # take the first 18 char of every set of inputs' name to label the output set

    readsin1=${reads1_list[i]}
    readsin2=${reads2_list[i]}

    base=$(basename $readsin1)
    label=${base:0:$label_length}_  # underscore _ added as separator

    if [ $STAR_TWO_PASS ]; then     
        STAR --runThreadN 24 \
             --genomeDir $INDEXED_GENOME_PATH \
             --readFilesIn $readsin1 $readsin2 \
             --outFileNamePrefix $PARENT_DIR/$ALIGNMENT_DIR/$label \
             --outSAMtype BAM SortedByCoordinate \
             --twopass1readsN $MAXIMUM_READS
             # sample-p01 mapped 12640877 reads here use a number that exceed the 1pass mapped reads to make sure they go through 2pass
    else
        STAR --runThreadN 24 \
             --genomeDir $INDEXED_GENOME_PATH \
             --readFilesIn $readsin1 $readsin2 \
             --outFileNamePrefix $PARENT_DIR/$ALIGNMENT_DIR/$label \
             --outSAMtype BAM SortedByCoordinate
    fi

done


# compressed alignment summaries into a .tar file

if [ STAR_TWO_PASS ]; then
    star_alignment_summary_compressed_label=star_alignment_2pass.summary.tar
else
    star_alignment_summary_compressed_label=star_alignment_1pass.summary.tar
fi

tar -cf $PARENT_DIR/$ALqIGNMENT_DIR/$star_alignment_summary_compressed_label $PARENT_DIR/$ALIGNMENT_DIR/*.final.out



#///////////////////////////
# separate output into independent dir
#///////////////////////////



##################
# 8. featureCounts  # separate output
##################

echo "starting feature counting..."

echo $PARENT_DIR/$ALIGNMENT_DIR/$STAR_MAPPED_SORTED_RESULTS_PATTERN

featureCounts -T 24 \
	      -a $ANNOTATION_PATH \
	      -o $PARENT_DIR/$FEATURE_COUNTS_DIR/$FEATURE_COUNTS_OUTPUT \
	      $PARENT_DIR/$ALIGNMENT_DIR/$STAR_MAPPED_SORTED_RESULTS_PATTERN

tar -cf $PARENT_DIR/$FEATURE_COUNTS_DIR/$FEATURE_COUNTS_OUTPUT.tar $PARENT_DIR/$FEATURE_COUNTS_DIR/*

echo "done featureCounts."





//////////////// Jenny's //////////////
fastqc -o $OUTFOLDER -t 16 -q --noextract ${file}

java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE -threads 64 ../../${file} ../../${BASE}-READ2-Sequences.txt.gz ${BASE}_f_p.fq.gz ${BASE}_f_up.fq.gz ${BASE}_r_p.fq.gz ${BASE}_r_up.fq.gz ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-PE.fa:2:30:10 TRAILING:15 MINLEN:20

tophat2 -p $CORES -G $GTFfile -r 100 --mate-std-dev 125 --no-coverage-search  -o $OUTFOLDER --library-type=fr-firststrand --transcriptome-index=$KNOWN $GENOME $READ1 $READ2
/////////////////////////////////////////

#!/bin/bash
#$ -N SeanTest
#$ -q cbcl-a64,cbcl-i24,bio,pub64,free64
#$ -pe openmp 8-64
#$ -m beas



<<COMMENT
function: Download_Fastq
params: $1 a list of fastq urls; $2 output dir prefix
returns: None
COMMENT

Download_Fastq() {
while read line
do
    wget -nv -P $2 $line
       # -P output_dir_prefix
       # -nv non-verbose
       # the output name will follow whatever it was originally named
done < $1
}



<<COMMENT
function: Decompress_gzFiles
params: $1 a list of files.gz - input dir ending w/o '/';
returns: None
COMMENT

Decompress_gzFiles() {
for data in $1/*
do
    gunzip $data
    # automatically replace .gz files with output .txt
done
}



<<COMMENT
function: Modify_Extensions
params: $1 a list of file.txt - input dir ending w/o '/'; $2 specify the extension you want for your data to follow (e.g: .fastq)
COMMENT

Modify_Extensions() {
for data in $1/*
do
    prefix=$(echo $data | cut -f1 -d '.')
    mv $data $prefix.$2
    # prefix includes path + data basename but excludes (the extension) everything after the first "." found
done
}


<<COMMENT
function: Remove_PrNotRecog, self-explanatory
params: $1 fastq dir
returns: None
COMMENT
Remove_PrNotRecog() {

rm $1/*PrNotRecog*

}



<<COMMENT
function: Check_md5sum
params: $1 input dir ending without '/'; $2 original_code
returns bool
COMMENT

Check_md5sum() {

for data in $1/*
do
    code=$(md5sum $data | cut -f1 -d ' ')
    echo $code >> md5sumListGenerated 
done

hashtemp=$(md5sum ./md5sumListGenerated | cut -f1 -d ' ')
hashoriginal=$(md5sum $2 | cut -f1 -d ' ')

if [[ $hashtemp == $hashoriginal ]]; then
    echo "md5sum check passed. Continue..."
else
    echo "Error occurred during data downloading..."
    exit
fi

}





<<COMMENT
function: just a wrapper around fastqc since it is a one-line thing and takes a list of inputs already
params: $1 input dir ending w/o '/'; $2 output dir; $3 help to specify input (*.fastq for raw data; *_Paired*.fastq for PE trimed data.)
return: None
COMMENT

FastQC() {
module load fastqc/0.11.2

if [[ $3 == "PEtrimed" ]]; then
    search_range=$1/*_Paired*.fastq
elif [[ $3 == "SEtrimed" ]]; then
    echo "under dev..."
else
    search_range=$1/*
fi

for data in $search_range
do
    name=$(basename $data)
    prefix=$(echo $name | cut -f1 -d '.')
    mkdir $2/$prefix
    fastqc -o $2/$prefix -t 16 -q $data
done
}




<<COMMENT
function: TrimReads just a wrapper around PE_Trim and SE_Trim
params: $1 specify SE:single-end || PE:pair-end; $2 input root dir ending without "/"; $3 output dir ending without "/"
returns: None
COMMENT

TrimReads() {

if [[ $1 == "SE" ]]; then
    # TRY NOT TO DO PUT DOUBLE QUOTE AROUND $1
    echo "single-end mode is still under dev."
elif [[ $1 == "PE" ]]; then
    PE_Trim $2 $3   
else
    echo "You must specify running mode as SE or PE."
    exit
fi
} 


<<COMMENT
<IMPORTANT>
function: the main func that exe PE trimmomatic
The default settings include 
1.path to call trimmomatic
2.the path for IlluminaClip:TrueSeq-3
3.the separator "READ1"&&"READ2" to distinguish two sets of input reads
4.set label "_Paired" && "_Unpaired", which will go as reference into mapping.
<IMPORTANT>
params: $1 the first argument passed from wrapper TrimReads; $2 the third argument "output dir" passed from wrapper TrimReads.
returns: None
COMMENT

PE_Trim() {

module load java/1.8.0.51

reads1_list=($1/*READ1*)
reads2_list=($1/*READ2*)

len1=${#reads1_list[@]}
len2=${#reads2_list[@]}

if [[ $len1==$len2 ]]; then
    echo "input size match! ready to loop into Trimmomatic..."
else
    echo "input sizes do not match! exiting..."
    exit
fi

for ((i = 0; i < $len1 && i < $len2; i++))
do
    echo "Trimmomatic New Round Started... "

    readsin1=${reads1_list[i]}
    readsin2=${reads2_list[i]}
    # label output originated from input plus "_Paired" with the third argument specifying if the name keeps the input path or not
    readsout1=$2/$(Labeled_With $readsin1 _Paired ExcludePath)
    readsout2=$2/$(Labeled_With $readsin1 _Unpaired ExcludePath)
    readsout3=$2/$(Labeled_With $readsin2 _Paired ExcludePath)
    readsout4=$2/$(Labeled_With $readsin2 _Unpaired ExcludePath)   
#    echo $readsin1
#    echo $readsin2
#    echo $readsout1
#    echo $readsout4

    java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE -threads 24 \
         -phred33 $readsin1 $readsin2 $readsout1 $readsout2 $readsout3 $readsout4 \
         ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-PE.fa:2:30:10 TRAILING:15 MINLEN:20
done

echo "Done trimming!"
}



<<COMMENT
function: Labeled_With 
params: $1 a filename ending with .extension; $2 a label string to be added right before the extension; $3 Flag default set to be IncludePath
returns: the filename (including path) + label + .extension
COMMENT

Labeled_With() {

PATH_SETTING=${3:-IncludePath}

if [[ $PATH_SETTING == "ExcludePath" ]]; then
    base=$(basename $1)
    prefix=$(echo $base | cut -f1 -d '.')
else
    prefix=$(echo $1 | cut -f1 -d '.')
fi

extension=$(echo $1 | cut -f2 -d '.')

echo $prefix$2.$extension
}




<<COMMENT
function: StarIndex, just a wrapper around the original cmd tool
params: $1 unindexed genome.fa; $2 annotation file; $3 (read length - 1); $4 output dir
returns: None
COMMENT

StarIndex() {

mkdir $4

STAR --runThreadN 24 \
     --runMode genomeGenerate \
     --genomeDir $4 \
     --genomeFastaFiles $1 \
     --sjdbGTFfile $2 \
     --sjdbOverhang $3
}




<<COMMENT
function: RunStar
params: $1 input dir ending w/o '/'; $2 output dir; $3 indexed genome; $4 specify maximum reads for two pass run
<IMPORTANT>
1. One default setting that inherits from TrimReads is the "_Paired" label.
2. Besides "READ1" && "READ2", "L1" && "L2" will be used to separate technical dups
<IMPORTANT>
returns: None
COMMENT

RunStar() {

module load STAR/2.5.1b-static

reads1_list1=($1/*L1*READ1*_Paired*)
reads1_list2=($1/*L2*READ1*_Paired*)
reads2_list1=($1/*L1*READ2*_Paired*)
reads2_list2=($1/*L2*READ2*_Paired*)

len1=${#reads1_list1[@]}
len2=${#reads2_list1[@]}

if [[ $len1==$len2 ]]; then
    echo input size match! ready to loop into STAR...
else
    echo input sizes do not match! exiting...
    exit
fi

for ((i = 0; i < $len1 && i < $len2; i++))
do
    echo "Mapping New Round Started..."    
    read1=${reads1_list1[i]}
    read1_dup=${reads1_list2[i]}
    read2=${reads2_list1[i]}
    read2_dup=${reads2_list2[i]}

    base1=$(basename $read1)
    base2=$(basename $read2)
    label=$(Common_Char $base1 $base2)
    mkdir $2/$label
    # --outFileNamePrefix -> underscore _ serves as separator in multiple output names
    STAR --runThreadN 24 \
         --genomeDir $3 \
         --readFilesIn $read1,$read1_dup $read2,$read2_dup \
         --outFileNamePrefix $2/$label/$label"_" \
         --outSAMtype BAM SortedByCoordinate \
	 --twopassMode Basic \
         --twopass1readsN $4
done
}



<<COMMENT
function Common_char 
params: $1 $2 are two input strings
return: the first N common characters in their names
COMMENT

Common_Char() {
i=0
while [[ ${1:i:1} == ${2:i:1} ]]
do
    i=$((i+1))
done
echo ${1:0:i}
}




<<COMMENT
function: a basic wrapper around featureCounts
params: $1 input dir ending without '/'; $2 output root dir ending without '/'; $3 annotation file;
returns: None
COMMENT

FeatureCounts() {

module load subread/1.5.0-p3

for dir in $1/*
do
    data=$dir/*Aligned.sortedByCoord.out.bam
    prefix=$(basename $dir)
    mkdir $2/$prefix
    featureCounts -T 24 \
                  -a $3 \
                  -o $2/$prefix/featureCounts_output \
                  $data
    echo $2/$prefix/featureCounts_output >> ./featureCountsList
done

echo "feature counting done!"
}


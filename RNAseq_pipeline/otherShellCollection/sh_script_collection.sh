Shell Script Collection

#=====================================

# How to loop through 2 arguments at once

#=====================================

im1_dir=(/cbcl/szang/seqCenter/test_run/raw_data/*READ1*.fastq)
im2_dir=(/cbcl/szang/seqCenter/test_run/raw_data/*READ2*.fastq)

for ((i = 0; i < ${#im1_dir[@]} && i < ${#im2_dir[@]}; i++))
do
    echo "${im1_dir[i]}" "${im2_dir[i]}"
done




#=====================================

# How to loop through all the arguments

#=====================================

#!/bin/bash

echo "You start with $# positional parameters"

# Loop until all parameters are used up

while [ "$1" != "" ]; do

    echo "Parameter 1 equals $1"

    echo "You now have $# positional parameters"

    # Shift all the parameters down by one

    shift

done

#=====================================

# How to read in str arguments

#=====================================

#!/bin/bash

echo -n "Enter some text > "

read text

echo "You entered: $text”



#=====================================

# How to get optional arguments for

# e.g. 

# command -a /cbcl/szang/yvonne/gencode.v19.annotation.gtf -i /cbcl/szang/yvonne/indexex/genome -# p /cbcl/szang/yvonne/add_qsub/ 

#=====================================


while getopts "a:i:p:" opt; do
    case $opt in
        a) ANNOTATION_FILE=$OPTARG ;;
        i) INDEX_PATH_PREFIX=$OPTARG ;;
        p) INPUT_PATH=$OPTARG ;;
    esac
done

echo $#

if [ $(( $# - $OPTIND )) -lt 2 ]; then
    echo "Usage: `basename` [options] ARG1 ARG2 ARG3"
    exit 1
fi


ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}
ARG3=${@:$OPTIND+2:1}

echo $ANNOTATION_FILE
echo $INDEX_PATH_PREFIX
echo $INPUT_PATH
echo $ARG1 
echo $ARG2 
echo $ARG3


#######################
# function build
#######################

say() 
{
CONTENT=$1
echo $CONTENT
}


say


#######################
# system_page
#######################

#!/bin/bash

# system_page - A script to produce a system information HTML file

##### Constants

TITLE="System Information for $HOSTNAME"
RIGHT_NOW=$(date +"%x %r %Z")
TIME_STAMP="Updated on $RIGHT_NOW by $USER"

##### Functions

function system_info
{
    echo "<h2>System release info</h2>"
    echo "<p>Function not yet implemented</p>"

}   # end of system_info


function show_uptime
{
    echo "<h2>System uptime</h2>"
    echo "<pre>"
    uptime
    echo "</pre>"

}   # end of show_uptime


function drive_space
{
    echo "<h2>Filesystem space</h2>"
    echo "<pre>"
    df
    echo "</pre>"

}   # end of drive_space


function home_space
{
    # Only the superuser can get this information

    if [ "$(id -u)" = "0" ]; then
        echo "<h2>Home directory space by user</h2>"
        echo "<pre>"
        echo "Bytes Directory"
        du -s /home/* | sort -nr
        echo "</pre>"
    fi

}   # end of home_space


function write_page
{
    cat <<- _EOF_
    <html>
        <head>
        <title>$TITLE</title>
        </head>
        <body>
        <h1>$TITLE</h1>
        <p>$TIME_STAMP</p>
        $(system_info)
        $(show_uptime)
        $(drive_space)
        $(home_space)
        </body>
    </html>
_EOF_

}

function usage
{
    echo "usage: system_page [[[-f file ] [-i]] | [-h]]"
}


##### Main

interactive=
filename=~/system_page.html

while [ "$1" != "" ]; do
    case $1 in
        -f | --file )           shift
                                filename=$1
                                ;;
        -i | --interactive )    interactive=1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done


# Test code to verify command line processing

if [ "$interactive" = "1" ]; then
	echo "interactive is on"
else
	echo "interactive is off"
fi
echo "output file = $filename"


if [ "$interactive" = "1" ]; then

    response=

    echo -n "Enter name of output file [$filename] > "
    read response
    if [ -n "$response" ]; then
        filename=$response
    fi

    if [ -f $filename ]; then
        echo -n "Output file exists. Overwrite? (y/n) > "
        read response
        if [ "$response" != "y" ]; then
            echo "Exiting program."
            exit 1
        fi
    fi
fi

# Write page (comment out until testing is complete)

write_page > $filename
o



##############################
# Tophat process
##############################

#!/bin/bash
#$ -N tophat
#$ -q cbcl-a64,cbcl-i24,bio,pub64,free64
#$ -pe openmp 8-64
#$ -m beas


module load bowtie2/2.2.7
module load tophat/2.1.0
module load samtools/1.1


ANNOTATION_FILE="/cbcl/szang/Chip-seq/tophat_IDX_ANTT_package/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"

INDEX_PATH_PREFIX="/cbcl/szang/Chip-seq/tophat_IDX_ANTT_package/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome"

SEARCH_PATTERN="ppd*.fastq"

for file in ./$SEARCH_PATTERN
do
  in_filename=$(basename $file)
  echo "in_filename is $in_filename"
  pref_len=$((${#in_filename}-6))
  label=${in_filename:0:$pref_len}
  echo "label is $label"
  out_filename=$label"_TophatOut"
  echo "out_filename is $out_filename"
  echo "ready to run tophat ..."
  tophat -p 32 -o ./$out_filename -g 1 --library-type fr-secondstrand -G $ANNOTATION_FILE $INDEX_PATH_PREFIX $file
  echo "$label mapping completed."
  cd $out_filename/
  sorted_bam_name=$label"_sorted.bam"
  sorted_sam_name=$label"_sorted.sam"
  echo "sorted bam name is $sorted_bam_name, sorted sam name is $sorted_sam_name."
  samtools sort -o $sorted_bam_name -@ 32 accepted_hits.bam
  echo "$label sorting completed."
  samtools view -h $sorted_bam_name > $sorted_sam_name
  echo "$label SAM is ready."
  cd ../
  echo "------------ preparing next sample ------------"
done



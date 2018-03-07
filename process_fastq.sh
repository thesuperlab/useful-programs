#!/bin/bash
#
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 12 Oct. 2017
#
# Process raw fastq files for further analysis. Includes g-unzipping, quality trimming, quality filtering,
# clipping short sequences, and converting final result to fasta format.
# INDIVIDUAL IDS MUST NOT HAVE DASH - OR UNDERSCORE _ IN THE ID NAME
#
# Usage:
#	./process_fastq.sh <input_dir> <output_dir>

input_dir=$1  
output_dir=$2

# make sure directory names end in a forward slash
if [ ""${input_dir:$(expr ${#string} - 1):1} != "/" ]; then
	input_dir=$input_dir/
fi
if [ ""${output_dir:$(expr ${#string} - 1):1} != "/" ]; then
	output_dir=$output_dir/
fi
if [ ! -e $output_dir ]; then
	mkdir $output_dir
fi

fastx_path=/usr/local/Cellar/fastx_toolkit/0.0.14_1/bin/
trimmer_cmd="${fastx_path}fastq_quality_trimmer -Q 33 -t 25" 
filter_cmd="${fastx_path}fastq_quality_filter -Q 33 -q 25 -p 80"
fasta_cmd="${fastx_path}fastq_to_fasta -Q 33 -r"
fastx_clipper="${fastx_path}fastx_clipper -l 5"

for in_filename in $( ls ${input_dir}); do
    id=$( echo $in_filename | sed s/[-_].*$// | tr ' ' '\n' | sort -u | tr '\n' ' ' | sed 's/ //g' )
	out_filename="debug.fna"
    if [[ $in_filename == *"R1"* ]]; then
        out_filename="${output_dir}${id}_R1.fna"
    elif [[ $in_filename == *"R2"* ]]; then
        out_filename="${output_dir}${id}_R2.fna"
    else
        echo "${id} not paired-end"
        out_filename="${output_dir}${id}.fna"
    fi
    gunzip -kc ${input_dir}$in_filename  | $trimmer_cmd | $filter_cmd | $fastx_clipper | $fasta_cmd | sed 's/>/>'${id}'_/g' > $out_filename
done

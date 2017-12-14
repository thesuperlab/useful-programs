#!/bin/bash
#
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# Updated: 11 Dec. 2017
#
# Combine same samples from multiple sequencing runs into single file per sample..
#
# Usage:
#	./combine_filtered_fastas.sh sorted_unique_fasta_filenames.txt input_dir/ output_dir/

input_file=$1
input_dir=$2 
output_dir=$3
dir_name $input_dir
dir_name $output_dir
if [ ! -e $output_dir ]; then
	mkdir $output_dir
fi
for filename in $(cat $input_file); do
	echo "processing $filename"
	for inter_dir in $(ls $input_dir); do
		input_filename=${input_dir}${inter_dir}/${filename}
		output_filename=${output_dir}${filename}
		if [ -e $input_filename ]; then
			cat $input_filename >> $output_filename
		fi
	done
done

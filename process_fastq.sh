#!/bin/bash
#
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 12 Oct. 2017
#
# Process raw fastq files for further analysis. Includes g-unzipping, quality trimming, quality filtering, 
# clipping short sequences, combining forward and reverse reads, and converting final result to fasta format.
# INDIVIDUAL IDS MUST NOT HAVE DASH - OR UNDERSCORE _ IN THE ID NAME
#
# Usage:
#	./process_fastq.sh <input_dir> <output_dir> <batch_size>

input_dir=$1  
output_dir=$2
batch_size=$3  # number of individuals per output file

# make sure directory names end in a forward slash
if [ ""${input_dir:$(expr ${#string} - 1):1} != "/" ]; then
	$input_dir=$input_dir/
fi
if [ ""${output_dir:$(expr ${#string} - 1):1} != "/" ]; then
	$output_dir=$output_dir/
if [ ! -e $output_dir ]; then
	mkdir $output_dir
echo "made dir $dir"

fastx_path=/usr/local/Cellar/fastx_toolkit/0.0.14_1/bin/
trimmer_cmd="${fastx_path}fastq_quality_trimmer -Q 33 -t 25" 
filter_cmd="${fastx_path}fastq_quality_filter -Q 33 -q 25 -p 80"
fasta_cmd="${fastx_path}fastq_to_fasta -Q 33 -r"
fastx_clipper="${fastx_path}fastx_clipper -l 5"

indiv_ids=$(ls $input_dir | sed s/[-_].*$// | tr ' ' '\n' | sort -u | tr '\n' ' ')
batch_indiv_count=0
batch_number=1

for id in ${indiv_ids[@]}; do
		let "batch_indiv_count+=1"
		if [ $batch_indiv_count -eq 1 ]; then
			first_indiv_id=$id
		fi
	
		batch_filename=${output_dir}/${batch_number}
		gunzip -kc ${input_dir}${id}*  | $trimmer_cmd | $filter_cmd | $fastx_clipper | $fasta_cmd | sed 's/>/>'${id}'_/g' >> $batch_filename
		
		# start new batch for next individual
		if [ $batch_indiv_count -eq $batch_size ]; then
			let "batch_number+=1"
			batch_indiv_count=0  
			# clean up intermediate file
			if [ $batch_size -eq 1 ]; then
				output_filename=${output_dir}${id}.fasta
			else
				output_filename=${output_dir}${first_indiv_id}-${id}.fasta
			fi
			mv $batch_filename $output_filename
		fi
done

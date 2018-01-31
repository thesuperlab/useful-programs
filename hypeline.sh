#!/bin/bash
# hypeline: a pipeline for haplotype calling using whatshap
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 30 Nov. 2017
#
# Usage:
#	./hypeline.sh <reference_name> <input_dir> <output_dir> 

input_dir=$2
output_dir=$3
snps_dir=snp_sites/
# make sure directory names end in a forward slash
if [ ""${input_dir:$(expr ${#string} - 1):1} != "/" ]; then
	$input_dir=$input_dir/
fi
if [ ""${output_dir:$(expr ${#string} - 1):1} != "/" ]; then
	$output_dir=$output_dir/
fi
reference_name=$1  # do not include .fna extension
reference_file=${reference_name}.fna
if [ ! -e ${reference_name}.1.bt2 ]; then
	bowtie2-build $reference_file $reference_name
fi
if [ ! -e ${reference_file}.fai ]; then
	samtools faidx $reference_file
fi
inter_dir=intermediates/
indiv_haps_dir=${inter_dir}individual_haps/
combined_haps_dir=${inter_dir}combined_haps/
directories=($output_dir $inter_dir $indiv_haps_dir $combined_haps_dir $snps_dir)
for dir in ${directories[@]}; do
	if [ ! -e $dir ]; then
		mkdir $dir
		echo "made dir $dir"
	fi
done

for filename in $(ls $input_dir); do
	id=$(echo $filename | sed -e 's/\..*//g')

	sorted_bam=${inter_dir}${id}.sorted.bam	
	bowtie2 -x $reference_name -f ${input_dir}${filename} | samtools view -S -b | samtools sort - -o $sorted_bam

	filtered_vcf=${inter_dir}${id}.filtered.vcf
	freebayes --min-base-quality 3 --min-mapping-quality 1 -f $reference_file $sorted_bam | bcftools filter -e 'QUAL < 20' | bcftools filter -e 'DP < 5' -o $filtered_vcf

	phased_vcf=${inter_dir}${id}.phased.vcf
	whatshap phase --ignore-read-groups --reference $reference_file -o $phased_vcf $filtered_vcf $sorted_bam

	bgzip $phased_vcf
	zipped_vcf=${phased_vcf}.gz
	tabix $zipped_vcf

	hap1=${indiv_haps_dir}${id}.hap1.fasta
	hap2=${indiv_haps_dir}${id}.hap2.fasta
	bcftools consensus -H 1 -f $reference_file $zipped_vcf > $hap1
	bcftools consensus -H 2 -f $reference_file $zipped_vcf > $hap2
done

haps_to_single_loci.py $indiv_haps_dir $combined_haps_dir $reference_file

for filename in $(ls $combined_haps_dir); do
	id=$(echo $filename | sed -e 's/\..*//g')
	aligned_file=${combined_haps_dir}${id}.aligned.haps
	mafft --retree 1 --maxiterate 0 --quiet ${combined_haps_dir}${filename} > $aligned_file
	haps_filename=${output_dir}${id}.haps
	remove_ref_seqs.py $aligned_file $haps_filename $reference_file
	snp-sites -o ${snps_dir}${id}.snps $haps_filename
done

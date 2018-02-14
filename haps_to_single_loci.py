#!/usr/local/bin/python3
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 30 Nov. 2017
#
# Usage:
#   haps_to_single_loci.py <input_dir> <output_dir> <reference_filename>
import Bio.SeqIO
import collections
import os
import sys


def main():

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    input_dir = input_dir if input_dir.endswith("/") else input_dir + "/"
    output_dir = output_dir if output_dir.endswith("/") else output_dir + "/"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    reference_filename = sys.argv[3]
    hap1_suffix = '.hap1.fna'
    hap2_suffix = '.hap2.fna'
    loci = collections.defaultdict(dict)
    input_filenames = os.listdir(input_dir)
    individual_ids = {filename.split('.')[0] for filename in input_filenames}
    for individual in individual_ids:
        for filename, hap_number in ((input_dir + individual + hap1_suffix, '_1'), (input_dir + individual + hap2_suffix, '_2')):
            for seq_record in Bio.SeqIO.parse(filename, 'fasta'):
                locus_id = seq_record.id
                seq_record.id = individual + hap_number
                loci[locus_id][seq_record.id] = seq_record
    reference_sequences = {seq_record.id: seq_record for seq_record in Bio.SeqIO.parse(reference_filename, 'fasta')}
    for locus_id in loci:
        sequences = [loci[locus_id][seq_id] for seq_id in sorted(loci[locus_id].keys())]
        sequences.insert(0, reference_sequences[locus_id])
        with open(output_dir + locus_id + '.haps.ref.fna', 'w') as file:
            Bio.SeqIO.write(sequences, file, 'fasta')

main()

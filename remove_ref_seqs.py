#!/usr/local/bin/python3
#
# Author: Kelly Sovacool
# Email: kellysovacool@uky.edu
# 11 Dec. 2017
#
# Remove reference sequences from a fasta file
#
# Usage:
#   remove_ref_seqs.py input_filename output_filename reference_filename

import Bio.SeqIO
import sys

def main():

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    reference_filename = sys.argv[3]
    ref_seq_ids = set(seq_record.id for seq_record in Bio.SeqIO.parse(reference_filename, 'fasta'))
    seqs_to_keep = [seq_record for seq_record in Bio.SeqIO.parse(input_filename, 'fasta') if seq_record.id not in ref_seq_ids]
    Bio.SeqIO.write(seqs_to_keep, output_filename, "fasta")

main()

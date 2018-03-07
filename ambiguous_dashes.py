#!/usr/local/bin/python3
"""
Convert non-nucleotides to dashes. If an entire sequence is all dashes, throw it out.
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
Date: 14 Feb. 2018

Usage:
    ./ambiguous_dashes.py <input_filename> <output_filename>

Options:
    -h --help   display this message

"""
import Bio.Seq
import Bio.SeqIO
import docopt
import os


def main(args):
    nucleotides = "agtcAGTC"
    out_records = list()
    with open(args['<input_filename>'], 'r') as infile:
        in_records = Bio.SeqIO.parse(infile, 'fasta')
        for record in in_records:
            new_seq = "".join(nuc if nuc in nucleotides else '-' for nuc in record.seq)
            if ''.join(set(new_seq)) != '-': # throw out sequences that are entirely dashes
                record.seq = Bio.Seq.Seq(new_seq)
                out_records.append(record)
            else:
                print("{} from {} has no nucleotides".format(new_seq.id, input_filename))
    Bio.SeqIO.write(out_records, args['<output_filename>'], 'fasta')


if __name__ == "__main__":
    main(docopt.docopt(__doc__))

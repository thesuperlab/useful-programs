#!/usr/local/bin/python3
""" Convert files in nexus format to fasta format 
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
06 Feb. 2018

Usage:
    ./nexus2fasta.py input_dir_or_file output_dir_or_file
"""
import os
import sys


def main():
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    if os.path.isdir(input_path):
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        input_path = input_path if input_path.endswith("/") else input_path + "/"
        output_path = output_path if output_path.endswith("/") else output_path + "/"
        for input_filename in os.listdir(input_path):
            output_filename = input_filename.split(".")[0] + '.fna'
            convert_nexus_to_fasta(input_path + input_filename, output_path + output_filename)
    elif os.path.isfile(input_path):
        convert_nexus_to_fasta(input_path, output_path)
    else:
        print("{} is neither file nor directory".format(input_path))


def convert_nexus_to_fasta(input_filename, output_filename):
    if os.path.exists(output_filename):
        print("output filename {} already exists.".format(output_filename))
        answer = input("overwrite? (y or n) ")
        answer = answer.lower()
    else:
        answer = "y"
    if "y" in answer or "yes" in answer:
        with open(output_filename, 'w') as outfile:
            with open(input_filename, 'r') as infile:
                beginning = True
                end = False
                for line in infile:
                    split_line = line.split()
                    if split_line:
                        if not beginning and not end and len(split_line) >= 2:
                            sequence = ">" + split_line[0].strip("'") + '\n' + max(split_line, key=lambda word: len(word)).strip("'") + "\n"
                            outfile.write(sequence)
                        elif split_line[0] == 'MATRIX':
                            beginning = False
                        elif split_line[0] == ';':
                            end = True
main()

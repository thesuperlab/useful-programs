#!/usr/bin/python
"""
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
Updated: 06 Oct. 2017

Takes a subsample file and creates bash scripts for running structure on the DLX.

Usage:
    ./make_strx_scripts.py <subsample.strx.txt>
"""
import sys


def main():
    structure_path = '/home/klso224/tools/structure-2.3.4/structure'  # change this to the location of structure on your account
    subsample_filename = sys.argv[1]
    email="kellysovacool@uky.edu"
    k_min = 2
    k_max = 8  # change this to the maximum k value you want to test
    runs = 16  # should be no greater than the number of threads available
    N = 0  # individuals
    L = 0  # loci
    with open(subsample_filename, 'r') as file:
        is_first = True
        for line in file:
            N += 1
            if is_first:
                is_first = False
                L = len(line.split()) - 1
    N = N / 2
    for k in range(k_min, k_max + 1):
        script = '#!/bin/bash\n#SBATCH --mail-type=ALL\n#SBATCH --mail-user=' + email + "\nSBATCH_NODELIST: $SBATCH_NODELIST\n\n"
        for run in range(1, runs + 1):
            script += structure_path + " -i " + subsample_filename + " -K " + str(k) + " -L " + str(L) + " -N " + str(N) + " -o " + subsample_filename.split('.')[0] + '_k' + str(k) + '_run' + str(run) + '.out &\nsleep 10\n\n'
        script += 'wait'
        with open('strx_' + subsample_filename.split('.')[0] + '_k' + str(k) + '.sh', 'w') as script_file:
            script_file.write(script)

main()

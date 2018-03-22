#!/usr/local/bin/python3
"""
Author: Kelly Sovacool
Email: kellysovacool@uky.edu
28 Sep 2017

Usage: 
    snp_subsample.py <snp-sites-dir> <output-filename-base> [--abundance-filter=<percent-cutoff> --output-filtered-fasta-dir=<filtered_dir> --skip-filter --size=<subsample-size> --all-snps-all-loci --missing-cutoff=<percent-cutoff> --output-format=<fasta-or-strx> --ignore-indels]
    snp_subsample.py --help

Options:
    -h --help                                   display this incredibly helpful message.
    --abundance-filter=<percent-cutoff>         filter out SNPs below a percent abundance cutoff [default: 0.05].
    --missing-cutoff=<percent-cutoff>           filter out SNP sites with more than <percent-cutoff> of individuals missing data [default: 0.5].
    --skip-filter                               skip the filtering step (e.g. if using already-filtered data).
    --size=<subsample-size>                     number of subsamples [default: 1].
    --all-snps-all-loci                         output a file containing all snps from all loci.
    --output-filtered-fasta-dir=<filtered_dir>  output filtered snp-sites fastas to a directory.
    --output-format=<fasta-or-strx>             output format for subsamples [default: strx].
    --ignore-indels                             if an indel is present at a snp site, ignore that site. default behavior is to report indels as missing data.
"""
import Bio.Seq
import Bio.SeqIO
import docopt
import os
import random

# Redesign: clean up main fcn, better OOP design


def main(args):
    args['--ignore-indels'] = '--ignore-indels' in args
    args['--skip-filter'] = '--skip-filter' in args
    all_individual_ids = set()
    for fasta_filename in os.listdir(args['<snp-sites-dir>']):  # build set of individual ids
        with open(args['<snp-sites-dir>'] + fasta_filename, 'r') as infile:
            all_individual_ids.update(record.id for record in Bio.SeqIO.parse(infile, 'fasta'))
    bad_individual_ids = set()
    for id in all_individual_ids:  # throw out individuals that aren't haplotyped (usually they're not real individuals)
        id_no_number = id.split('_')[0]
        if id_no_number + '_1' not in all_individual_ids or id_no_number + '_2' not in all_individual_ids:
            bad_individual_ids.add(id)
    all_individual_ids.difference_update(bad_individual_ids)

    loci = Loci()  # build loci
    for fasta_filename in sorted(os.listdir(args['<snp-sites-dir>'])):
        with open(args['<snp-sites-dir>'] + fasta_filename, 'r') as infile:
            locus_id = fasta_filename.split('.')[0]
            locus = Locus(locus_id, Bio.SeqIO.parse(infile, 'fasta'), all_individual_ids, bad_individual_ids)
            if not args['--skip-filter']:
                locus.filter(float(args['--abundance-filter']), float(args['--missing-cutoff']), ignore_indels=args['--ignore-indels'])
                if args['--output-filtered-fasta-dir']:
                    with open(args['--output-filtered-fasta-dir'] + fasta_filename, 'w') as filtered_fasta:
                        for seq_id, sequence in locus.individuals.items():
                            filtered_fasta.write('>' + seq_id + '\n')
                            filtered_fasta.write(sequence + '\n')
        loci.append(locus)
    output_file_extension = output_format_extension(args['--output-format'])

    if args['--size']:  # take subsample(s)
        for number in range(1, int(args['--size']) + 1):
            loci.new_random_poly_sites()
            with open(args['<output-filename-base>'] + str(number) + output_file_extension, 'w') as outfile:
                for id in sorted(all_individual_ids):
                    line = id + '\t' if output_file_extension != '.strx.txt' else id.split('_')[0] + '\t'
                    for locus in loci:
                        nuc = locus.individuals[id][locus.random_poly_site] if id in locus.individuals else '-'
                        if output_file_extension == '.strx.txt':
                            nuc = nucleotide_to_structure(nuc)
                        line += nuc + '\t'
                    line += '\n'
                    outfile.write(line)

    if args['--all-snps-all-loci']:
        with open('all-snps-all-loci' + output_file_extension, 'w') as outfile:
            for id in sorted(all_individual_ids):
                line = id + '\t' if output_file_extension != '.strx.txt' else id.split('_')[0] + '\t'
                for locus in loci:
                    seq = locus.individuals[id] if id in locus.individuals else '-'*len(locus.polymorphic_sites)
                    if output_file_extension == '.strx.txt':
                        strx_seq = ''
                        for nuc in seq:
                            strx_seq += nucleotide_to_structure(nuc) + '\t'
                        line += strx_seq
                    else:
                        line += seq
                line += '\n'
                outfile.write(line)


class Loci(list):
    def __init__(self, iterable=None):
        if iterable:
            super().__init__(iterable)
        else:
            super().__init__()

    def filter(self, abundance_cutoff=0.05, missing_cutoff=0.5, ignore_indels=False):
        for locus in self:
            locus.filter(abundance_cutoff, missing_cutoff, ignore_indels=ignore_indels)

    def new_random_poly_sites(self):
        for locus in self:
            locus.new_random_poly_site()


class Locus:
    def __init__(self, id, seq_records, all_individual_ids, bad_individual_ids):
        self.id = id
        self.individuals = {seq_rec.id: str(seq_rec.seq) for seq_rec in seq_records}  # id: sequence

        length = max(len(seq) for seq in self.individuals.values())
        for missing_id in all_individual_ids - set(self.individuals.keys()):
            #print('individual', missing_id, 'not in locus', self.id, '-- filling in nucleotides as dashes')
            self.individuals[missing_id] = '-' * length
        for bad_id in bad_individual_ids:
            self.individuals.pop(bad_id, False)

        self.polymorphic_sites = []
        first_individual = True
        for indiv_id, sequence in self.individuals.items():
            i = 0
            for nucleotide in sequence:
                if first_individual:
                    self.polymorphic_sites.append(PolymorphicSite())
                if nucleotide not in self.polymorphic_sites[i].snps and nucleotide in SNP.nucleotides:
                    self.polymorphic_sites[i].snps[nucleotide] = SNP(nucleotide)
                    self.polymorphic_sites[i].snps[nucleotide].individual_ids.add(indiv_id)
                i += 1
            first_individual = False
        self.random_poly_site = None
        self.new_random_poly_site()  # use same poly site for all subsamples in this run

    def __str__(self):
        return self.id

    def filter(self, abundance_cutoff, missing_cutoff, ignore_indels=False):
        original_count_sites = len(self.polymorphic_sites)
        removed_sites = 0
        i = 0
        while i < (original_count_sites - removed_sites):
            pmsite = self.polymorphic_sites[i]
            if (ignore_indels and pmsite.indel_or_missing_exists) or pmsite.has_low_abundance(len(self.individuals), abundance_cutoff) or pmsite.number_individuals_missing_data >= missing_cutoff:
                self.polymorphic_sites.pop(i)
                for id, sequence in self.individuals.items():
                    self.individuals[id] = sequence[:i] + sequence[i+1:]  # update sequence for individuals
                removed_sites += 1
            i += 1
        print(removed_sites, ' SNP sites thrown out for contig ', self.id, '. ', len(self.polymorphic_sites), ' SNP sites remain.', sep='')

    def new_random_poly_site(self):
        self.random_poly_site = random.randrange(len(self.polymorphic_sites))


class PolymorphicSite:

    def __init__(self):
        self.snps = {}  # nucleotide: SNP object

    def __len__(self):
        return len(self.snps)

    def has_low_abundance(self, total_individuals, abundance_cutoff):
        most_popular_snp = max(self.snps.values(), key=lambda snp: len(snp.individual_ids))
        others = set(self.snps.values())
        others.remove(most_popular_snp)
        num_individuals = len({indiv for snp in others for indiv in snp.individual_ids})
        return num_individuals / total_individuals <= abundance_cutoff

    @property
    def number_individuals_missing_data(self):
        individuals = set()
        for nuc, snp in self.snps.items():
            if nuc not in SNP.nucleotides:
                individuals.update(snp.individual_ids)
        return len(individuals)

    @property
    def indel_or_missing_exists(self):
        return bool(len({nuc for nuc in self.snps if nuc == '-'}))


class SNP:
    nucleotides = {'A': 1, 'C': 2, 'G': 3, 'T': 4}

    def __init__(self, nucleotide, individual_ids=None):
        super().__init__()
        self.nucleotide = nucleotide
        self.individual_ids = {id for id in individual_ids} if individual_ids else set()

    def __hash__(self):
        return hash((self.nucleotide, tuple(sorted(self.individual_ids))))


def output_format_extension(format):
    if format == 'strx':
        extension = '.strx.txt'
    elif format == 'fasta':
        extension = '.fasta'
    else:
        raise ValueError('unrecognized output format {}'.format(format))
    return extension


def nucleotide_to_structure(nucleotide):
    if nucleotide in SNP.nucleotides:
        structure_int = SNP.nucleotides[nucleotide]
    else:
        structure_int = -9
    return str(structure_int)


def check_directory(dir_name):
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        if not dir_name.endswith('/'):
            dir_name = dir_name + '/'
        return dir_name


if __name__ == "__main__":
    arguments = docopt.docopt(__doc__)
    for arg in ("<snp-sites-dir>", "--output-filtered-fasta-dir"):
        if arguments[arg]:
            arguments[arg] = check_directory(arguments[arg])
    main(arguments)


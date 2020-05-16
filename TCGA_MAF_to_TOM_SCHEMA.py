import sys
from typing import List
from os import listdir, path
import re
import shutil
from pathlib import Path
from distutils.dir_util import copy_tree

INPUT_FILE_PATTERN = r'.+\.bed$'
OUTPUT_FIELDS_TO_BE_UNIQUE = [0, 1, 8, 9]  # chrom, start, ref, alt
try :
    input_dir_location = sys.argv[1]
    output_dir_location = sys.argv[2]
except Exception:
    print('this program requires\n'
          '1. the path to the directory containing the files to process.\n'
          '2. the output directory path.\n'
          f'This module will convert files matching the regex {INPUT_FILE_PATTERN}'
          ' by removing duplicated variants, splitting REF AL1 AL1 into REF ALT AL1 AL2, and removing from REF and ALT ' 
          'equal sequences of nucleotides preceding the variant.\n'
          f'The files that don\'t match {INPUT_FILE_PATTERN} are copied as is into the destination folder.')
    sys.exit(1)
if input_dir_location[-1] != path.sep:
    input_dir_location += path.sep
if output_dir_location[-1] != path.sep:
    output_dir_location += path.sep
# create clean output dir if it doesn't exists
odp = Path(output_dir_location)
if odp.exists() and odp.is_dir():
    shutil.rmtree(odp)
odp.mkdir(parents=True, exist_ok=True)
log_file = open('./log.log', mode='w')

valid_chromosomes = {'chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY', 'chrMT', 'chr23', 'chr24', 'chr25',
                     '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25'}


def minimal_representation_of_variant(fields_of_variant: List[str]) -> str:
    fields_to_be_unique = [fields_of_variant[i] for i in OUTPUT_FIELDS_TO_BE_UNIQUE]
    return '\t'.join(fields_to_be_unique)


def split_variants(input_line: str, skip_weird_chromosomes:bool=True) -> List[List[str]]:
    """
    In general this method returns two arrays formed like below:
    original variants [chrom: reference allele] + [alt1] + [al1, al2] + original_variant[id:]
    original variants [chrom: reference allele] + [alt2] + [al1, al2] + original_variant[id:]
    However, if alt1 and alt2 have the same value, then only one array is returned, with both al1 and al2 = to 1.
    """
    fields = input_line.split('\t')
    output = list()
    if skip_weird_chromosomes and fields[0] not in valid_chromosomes:
        return output
    else:
        trailing_fields = fields[11:]
        if fields[8] != fields[9]:  # ref != alt1
            output.append(
                fields[0:10] + ['1', '0'] + trailing_fields
            )
        if fields[8] != fields[10]: # ref != alt2
            if fields[9] == fields[10]: # alt1 == alt2
                # i.e. if there's a second variant, but it's equal to the first one
                output[0][11] = '1'
            else:
                output.append(
                    fields[0:9] + [fields[10], '0', '1'] + trailing_fields
                )
        return output
OUTPUT_IDX = {
    'chrom': 0,
    'start': 1,
    'stop': 2,
    'strand': 3,
    'mut_type': 7,
    'ref': 8,
    'alt': 9,
    'al1': 10,
    'al2': 11,
    'id': 12
}


def remove_common_allele_prefix(ref: str, alt: str, start:str, stop:str):
    idx = 0
    while idx < len(ref) and idx < len(alt) and ref[idx] == alt[idx]:
        idx += 1
    # update start
    start = int(start)+idx
    # safety measure... you never know
    if int(stop) < start:
        stop = start
    return ref[idx:], alt[idx:], str(start), str(stop)


def remove_slash_for_empty_alleles(allele_string) -> str:
    if allele_string == '-':
        return ''
    else:
        return allele_string


def remove_novel_variant_ids(id_string:str) -> str:
    return '' if id_string == 'novel' else id_string


def convert_1_to_0_based_coordinates(start:str, stop:str, mut_type:str):
    start, stop = int(start), int(stop)
    start -= 1
    if mut_type == 'INS':
        stop = start
    return str(start), str(stop)


def find_nth_occurrence_in_string(what: str, in_str: str, n_th):
    start = in_str.find(what)
    while start >= 0 and n_th > 1:
        start = in_str.find(what, start + len(what))
        n_th -= 1
    return start

# example input
pos ='0     1           2           3   4       5       6                   7   8   9   10  11      12                              13                              14      15      16'
s1 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	G	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-10A-01D-A19M-08	null	null	055f269a-df3a-4063-a414-59e6a33cbba2'
s2 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	G	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
s4 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	T	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
s5 = '7	140453139	140453163	1	CNTRL	11064	Missense_Mutation	INS	TAGCTAGACCAAAATCACCTATTT	TAGCTAGACCAAAATCACCTATTT	TAGCTAGACCAAAATCACCTATTTTAGCTAGACCAAAATCACCTATTT	rs121913368'
s52 = '7	140453139	140453193	1	CNTRL	11064	Missense_Mutation	INS	TAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATAT	TAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATAT	TAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATATTAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATAT	rs121913370'
s6 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	A	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
s7 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	G	G	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
# the following line is not a real case
s8 = 'chrGL000193.1	8	8	+	CNTRL	11064	Missense_Mutation	SNP	G	GAA	G	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'


def transform_line(input_l: str, already_transformed_outputs: set):
    output_lines = list()
    # separate variants
    replacement_lines: List[List[str]] = split_variants(input_l, skip_weird_chromosomes=True)
    for var in replacement_lines:
        # remove prefix nucleotides
        var[OUTPUT_IDX['ref']], var[OUTPUT_IDX['alt']] , var[OUTPUT_IDX['start']], var[OUTPUT_IDX['stop']]= remove_common_allele_prefix(
            var[OUTPUT_IDX['ref']],
            var[OUTPUT_IDX['alt']],
            var[OUTPUT_IDX['start']],
            var[OUTPUT_IDX['stop']])
        var[OUTPUT_IDX['ref']] = remove_slash_for_empty_alleles(var[OUTPUT_IDX['ref']])
        var[OUTPUT_IDX['alt']] = remove_slash_for_empty_alleles(var[OUTPUT_IDX['alt']])
        var[OUTPUT_IDX['start']], var[OUTPUT_IDX['stop']] = convert_1_to_0_based_coordinates(var[OUTPUT_IDX['start']], var[OUTPUT_IDX['stop']], var[OUTPUT_IDX['mut_type']])
        var[OUTPUT_IDX['id']] = remove_novel_variant_ids(var[OUTPUT_IDX['id']])
        # return the transformed variant only if it new in already_transformed_output
        var_rep = minimal_representation_of_variant(var)
        if var_rep not in already_transformed_outputs:
            already_transformed_outputs.add(var_rep)
            # recompose line with all fields
            output_lines.append(
                '\t'.join(var)
            )
    return output_lines


def transform_files():
    print(f'reading content of dir {input_dir_location}')
    print(f'transformed files will be written into {output_dir_location}')
    print('if needed, you can see the list of files transformed in the file log.log in the same directory of the script.')
    # find files to_transform
    num_files = 0
    names_files_to_transform = list()
    for file in listdir(input_dir_location):
        num_files += 1
        if re.match(INPUT_FILE_PATTERN, file):
            names_files_to_transform.append(file)
        elif path.isfile(input_dir_location+file):
            shutil.copy(input_dir_location+file, output_dir_location)
        else:
            copy_tree(input_dir_location+file, output_dir_location+file)

    print(f'{num_files} total files of {input_dir_location} to write into {output_dir_location}')
    print(f'{len(names_files_to_transform)} files to transform')

    for file_name in names_files_to_transform:
        print(f'transforming {file_name}', end='\t...\t', file=log_file)
        with open(input_dir_location+file_name, 'r', encoding='utf-8') as input_file:
            with open(output_dir_location+file_name, 'w', encoding='utf-8') as output_file:
                transformed_variants = set()
                # transform input file line by line
                num_lines = 0
                for line in input_file:
                    num_lines += 1
                    line = line.rstrip('\n')
                    for output_lines in transform_line(line, transformed_variants):
                        output_file.write(output_lines+'\n')
                if len(transformed_variants) == 0 and num_lines != 0:
                    print(f'was reading file {input_file.name}')
                    raise Exception(f'was reading file {input_file.name} but hey, no region was transformed!')
        print('done', file=log_file)


try :
    transform_files()
finally:
    log_file.close()

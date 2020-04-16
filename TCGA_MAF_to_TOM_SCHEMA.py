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


def minimal_representation_of_variant(fields_of_variant: List[str]) -> str:
    fields_to_be_unique = [fields_of_variant[i] for i in OUTPUT_FIELDS_TO_BE_UNIQUE]
    return '\t'.join(fields_to_be_unique)


def split_variants(input_line: str) -> List[List[str]]:
    fields = input_line.split('\t')
    output = list()
    trailing_fields = fields[11:]
    if fields[8] != fields[9]:
        output.append(
            fields[0:10] + ['1', '0'] + trailing_fields
        )
    if fields[8] != fields[10]:
        if fields[9] == fields[10]:
            # i.e. if there's a second variant, but it's equal to the first one
            output[0][11] = '1'
        else:
            output.append(
                fields[0:9] + [fields[10], '0', '1'] + trailing_fields
            )
    return output


def remove_prefix(ref: str, alt: str):
    idx = 0
    while idx < len(ref) and idx < len(alt) and ref[idx] == alt[idx]:
        idx += 1
    return ref[idx:], alt[idx:]


def find_nth_occurrence_in_string(what: str, in_str: str, n_th):
    start = in_str.find(what)
    while start >= 0 and n_th > 1:
        start = in_str.find(what, start + len(what))
        n_th -= 1
    return start

# example input
s1 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	G	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-10A-01D-A19M-08	null	null	055f269a-df3a-4063-a414-59e6a33cbba2'
s2 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	G	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
s3='0       1           2           3   4       5       6                   7   8   9   10  11      12                              13                              14      15      16'
s4 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	T	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
s5 = '7	140453139	140453163	1	CNTRL	11064	Missense_Mutation	INS	TAGCTAGACCAAAATCACCTATTT	TAGCTAGACCAAAATCACCTATTT	TAGCTAGACCAAAATCACCTATTTTAGCTAGACCAAAATCACCTATTT	rs121913368'
s52 = '7	140453139	140453193	1	CNTRL	11064	Missense_Mutation	INS	TAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATAT	TAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATAT	TAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATATTAGCTAGACCAAAATCACCTATTTTTACTGTGAGGTCTTCATGAAGAAATATAT	rs121913370'
s6 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	A	A	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'
s7 = 'chr9	123936008	123936008	+	CNTRL	11064	Missense_Mutation	SNP	G	G	G	null	TCGA-BJ-A2NA-01A-12D-A19J-08	TCGA-BJ-A2NA-11A-11D-A19J-08	null	null	c88a3e4d-4316-4bc3-b2fa-0ac3dd76e558'


def transform_line(input_l: str, already_transformed_outputs: set):
    output_lines = list()
    # separate variants
    replacement_lines: List[List[str]] = split_variants(input_l)
    for var in replacement_lines:
        # remove prefix nucleotides
        var[8], var[9] = remove_prefix(var[8], var[9])
        # check if such a variant has already been transformed
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
    # find files to_transform
    num_files = 0
    names_files_to_transform = list()
    for file in listdir(input_dir_location):
        num_files +=1
        if re.match(INPUT_FILE_PATTERN, file):
            names_files_to_transform.append(file)
        elif path.isfile(input_dir_location+file):
            shutil.copy(input_dir_location+file, output_dir_location)
        else:
            copy_tree(input_dir_location+file, output_dir_location+file)

    print(f'{num_files} children to copy into the {output_dir_location}')
    print(f'{len(names_files_to_transform)} files to transform')

    for file_name in names_files_to_transform:
        print(f'transforming {file_name}', end='\t...\t')
        with open(input_dir_location+file_name, 'r', encoding='utf-8') as input_file:
            with open(output_dir_location+file_name, 'w', encoding='utf-8') as output_file:
                transformed_variants = set()
                # transform input file line by line
                for line in input_file:
                    line = line.rstrip('\n')
                    for output_lines in transform_line(line, transformed_variants):
                        output_file.write(output_lines+'\n')
        print('done')


transform_files()
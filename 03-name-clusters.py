#!/bin/env python3

import os
import re

from os.path import join as opj
from shutil import copyfile

from collections import Counter

from misc import list_of_files_at_path


organelles = ['cp', 'mt']

for organelle in organelles:

    input_dir = opj('.', 'data', 'output', 'Embryophyta', f'{organelle}')
    output_dir_1 = opj('.', 'data', 'output', 'Embryophyta', f'{organelle}-1-renamed')
    output_dir_2 = opj('.', 'data', 'output', 'Embryophyta', f'{organelle}-2-combined')

    os.makedirs(name=output_dir_1, exist_ok=True)
    os.makedirs(name=output_dir_2, exist_ok=True)

    file_paths, e = list_of_files_at_path(input_dir)
    used_names = list()
    for fp in file_paths:
        seq_count = 0
        with open(fp, 'r') as f:
            names = list()
            for line in f:
                if line.startswith('>'):
                    seq_count += 1
        # cluster_name = os.path.basename(f.name)
        # if seq_count < 2:
        #     print(f'{cluster_name} -> {cluster_name}')
        #     copyfile(fp, opj(output_dir_1, f'{cluster_name}.fasta'))
        #     continue
        with open(fp, 'r') as f:
            names = list()
            for line in f:
                if not line.startswith('>'):
                    continue
                # r = re.findall('^>(\\w+)\\|(\\w+)\\|\\|(.+)\\|\\|(\\w+)\\|(.+)\\|\\|(\\d+:\\d+)\\|\\|(.+)\\|(\\d+)\\|(\\d+)$', line)
                r = re.findall('^>(\\w+)\\|(\\w+)\\|\\|(.+)\\|\\|(\\w+)\\|(.+)\\|\\|(.+)\\|(\\d+)\\|(\\d+)$', line)
                names.append(r[0][5].replace('/', '__').replace(')', '').replace('(', ''))
                names_unique = set(names)
            counts = [names.count(x) for x in names_unique]
            choice = [item for items, c in Counter(names).most_common() for item in [items] * c][0]

            if choice in ('ndhl',):
                choice = 'ndhi'

            if choice in ('ycf10',):
                choice = 'cema'

            if choice in ('tic214',):
                choice = 'ycf1'

            if choice in ('pafii',):
                choice = 'ycf4'

            if choice in ('psb30',):
                choice = 'ycf12'

            if choice in ('ribosomal_protein_s15',):
                choice = 'rps15'

            if choice in ('ribosomal_protein_l32',):
                choice = 'rpl32'

            if choice in ('photosystem_ii_protein_k',):
                choice = 'psbk'

            # print(counts, names_unique, choice)
            used_names.append(choice)
            c = used_names.count(choice)
            choice += '___' + str(c)
            copyfile(fp, opj(output_dir_1, f'{choice}.fasta'))
            # print(counts, names_unique, choice)
            print(f'{cluster_name} -> {choice}')

    file_paths, e = list_of_files_at_path(output_dir_1)

    for fp in file_paths:
        cluster_name = os.path.basename(fp.split('___')[0])
        if '_cluster_' in cluster_name:
            continue
        cluster_name = cluster_name.split('.fasta')[0]
        cluster_name = cluster_name.split('_')[0]
        cluster_name = cluster_name.split('-')[0]
        fp_combined = opj(output_dir_2, f'{cluster_name}.fasta')
        with open(fp_combined, 'w') as f_combined:
            print(fp_combined)

    for fp in file_paths:
        cluster_name = os.path.basename(fp.split('___')[0])
        if '_cluster_' in cluster_name:
            continue
        cluster_name = cluster_name.split('.fasta')[0]
        cluster_name = cluster_name.split('_')[0]
        cluster_name = cluster_name.split('-')[0]
        fp_combined = opj(output_dir_2, f'{cluster_name}.fasta')
        with open(fp, 'r') as f:
            with open(fp_combined, 'a') as f_combined:
                for line in f:
                    f_combined.write(line)
                    # if line.startswith('>'):
                    #     # r = re.findall('^>(\\w+)\\|(\\w+)\\|\\|(.+)\\|\\|(\\w+)\\|(.+)\\|\\|(\\d+:\\d+)\\|\\|(.+)\\|(\\d+)\\|(\\d+)$', line)
                    #     r = re.findall('^>(\\w+)\\|(\\w+)\\|\\|(.+)\\|\\|(\\w+)\\|(.+)\\|\\|(.+)\\|(\\d+)\\|(\\d+)$', line)
                    #     organism = r[0][2]
                    #     f_combined.write(f'>{organism}\n')
                    # else:
                    #     f_combined.write(line)

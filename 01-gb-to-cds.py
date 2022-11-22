#!/bin/env python3

import os
from pprint import pp
from Bio import SeqIO
from ncbi_taxonomy_local import Taxonomy

input_dir = os.path.join('.', 'data', 'input')
output_dir = os.path.join('.', 'data', 'output')

acc_to_exclude = {'NC_059016.1', 'NC_059015.1', 'NC_059012.1', 'NC_023209.1', 'NC_034982.1', 'NC_050882.1', 'NC_062584.1', 'NC_039660.1'}
# acc_to_exclude = []

# ----------------------------------------------------------------------------
taxonomy = Taxonomy()
taxonomy.init('~/NCBI_TAXONOMY_DB')  # a path to a directory where the database will be stored.
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
tax_def_Viridiplantae = {
    'name': 'Viridiplantae',
    'taxa_to_include': {'Viridiplantae'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Chlorophyta = {
    'name': 'Chlorophyta',
    'taxa_to_include': {'Chlorophyta'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Streptophyta = {
    'name': 'Streptophyta',
    'taxa_to_include': {'Streptophyta'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Embryophyta = {
    'name': 'Embryophyta',
    'taxa_to_include': {'Embryophyta'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Tracheophyta = {
    'name': 'Tracheophyta',
    'taxa_to_include': {'Tracheophyta'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Spermatophyta = {
    'name': 'Spermatophyta',
    'taxa_to_include': {'Spermatophyta'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Liliopsida = {
    'name': 'Liliopsida',
    'taxa_to_include': {'Liliopsida'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

tax_def_Pentapetalae = {
    'name': 'Pentapetalae',
    'taxa_to_include': {'Pentapetalae'},
    'taxa_to_exclude': set(),
    'taxids_to_include': set(),
    'taxids_to_exclude': set()
}

# ----------------------------------------------------------------------------
tax_defs = [
    # tax_def_Viridiplantae,
    # tax_def_Chlorophyta,
    # tax_def_Streptophyta,
    tax_def_Embryophyta,
    # tax_def_Tracheophyta,
    # tax_def_Spermatophyta,
    # tax_def_Liliopsida,
    # tax_def_Pentapetalae,
]
# ----------------------------------------------------------------------------
for tax_def in tax_defs:

    tax_def_name = tax_def['name']
    tax_dir = os.path.join(output_dir, tax_def_name)
    if os.path.exists(tax_dir):
        continue

    records = SeqIO.parse(os.path.join(input_dir, 'short_list.gb'), 'gb')

    taxa_to_include = tax_def['taxa_to_include']
    taxa_to_exclude = tax_def['taxa_to_exclude']
    taxids_to_include = tax_def['taxids_to_include']
    taxids_to_exclude = tax_def['taxids_to_exclude']

    for _ in taxa_to_include:
        _ = taxonomy.taxids_for_name(_)
        _ = int(_['tax_ids'][0]['tax_id'])
        taxids_to_include.add(_)

    for _ in taxa_to_exclude:
        _ = taxonomy.taxids_for_name(_)
        _ = int(_['tax_ids'][0]['tax_id'])
        taxids_to_exclude.add(_)

    for _ in taxids_to_include.copy():
        _ = taxonomy.all_descending_taxids(_)
        if _ is None:
            continue
        _ = set([int(x) for x in _])
        taxids_to_include = taxids_to_include.union(_)

    for _ in taxids_to_exclude.copy():
        _ = taxonomy.all_descending_taxids(_)
        if _ is None:
            continue
        _ = set([int(x) for x in _])
        taxids_to_exclude = taxids_to_exclude.union(_)

    taxids_to_include = taxids_to_include - taxids_to_exclude
    # ----------------------------------------------------------------------------

    os.makedirs(name=tax_dir, exist_ok=True)

    # out_file_base = os.path.join(tax_dir, f'{tax_def_name}')
    out_file_base = os.path.join(tax_dir, '')

    with open(f'{out_file_base}cp.fasta', 'w') as _: print(f'Creating {_.name}.')
    with open(f'{out_file_base}mt.fasta', 'w') as _: print(f'Creating {_.name}.')

    cds_descriptions = list()
    cds_descriptions_temp = list()
    cds_seqs_temp = list()
    seq_ids = list()

    for r in records:

        if r.id in acc_to_exclude:
            continue

        r_source = r.features[0]
        assert r_source.type == 'source'
        _ = r_source.qualifiers['db_xref']
        assert len(_) == 1
        r_taxid = int(_[0].replace('taxon:', ''))

        if r_taxid not in taxids_to_include:
            continue

        r_organism = r_source.qualifiers['organism'][0]
        r_phylum = taxonomy.higher_rank_for_taxid(r_taxid, 'phylum')
        r_order = taxonomy.higher_rank_for_taxid(r_taxid, 'order')
        r_family = taxonomy.higher_rank_for_taxid(r_taxid, 'family')
        r_organelle = r_source.qualifiers['organelle'][0]

        if 'plastid' in r_organelle:
            r_organelle = 'cp'
        elif 'mitochondrion' in r_organelle:
            r_organelle = 'mt'
        else:
            print(f'Unexpected organelle: {r_organelle}')
            break

        taxid_organism = taxonomy.scientific_name_for_taxid(r_taxid)
        assert r_organism == taxid_organism

        # r_description = f'{r_phylum}|{r_order}|{r_family}||{r_organism}||{r_taxid}|{r_organelle}|{r.id}'
        r_description = f'{r_order}|{r_family}||{r_organism}||{r_organelle}|{r.id}'
        r_description = r_description.replace(' ', '_')

        r_cds_count = 0
        r_features = sorted(r.features, key=lambda x: min(x.location.nofuzzy_start, x.location.nofuzzy_end))
        for feature in r_features:
            if feature.type == 'CDS':
                r_cds_count += 1
                ft = feature

                cds_nt = ''
                for part in ft.location.parts:
                    start = part.start
                    end = part.end
                    if part.strand == -1:
                        cds_nt += str(r.seq[start:end].reverse_complement()).lower()
                    elif part.strand == 1:
                        cds_nt += str(r.seq[start:end]).lower()

                cds_codon_start = abs(int(ft.qualifiers['codon_start'][0])) - 1
                cds_nt = cds_nt[cds_codon_start:]
                cds_nt_len = len(cds_nt)
                cds_codon_len = (cds_nt_len // 3) * 3
                cds_nt = cds_nt[:cds_codon_len]

                cds_diff_nt_aa = 0
                if 'translation' not in ft.qualifiers:
                    assert 'pseudo' in ft.qualifiers
                    cds_aa = None
                else:
                    assert 'translation' in ft.qualifiers
                    cds_aa = ft.qualifiers['translation'][0]

                assert 'gene' or 'product' or 'protein_id' in ft.qualifiers

                if 'gene' in ft.qualifiers:
                    cds_name_key = 'gene'
                elif 'product' in ft.qualifiers:
                    cds_name_key = 'product'
                else:
                    cds_name_key = 'protein_id'

                cds_name = ft.qualifiers[cds_name_key][0]
                if cds_name in ('hypothetical protein', 'orf'):
                    # cds_name = ft.qualifiers['protein_id'][0]
                    continue
                if cds_aa is None:
                    # cds_name = f'{cds_name}|PSEUDOGENE'
                    continue
                if cds_name.lower().startswith('orf'):
                    # cds_name = ft.qualifiers['protein_id'][0]
                    continue

                cds_name = cds_name.lower()

                ft_start = min(ft.location.nofuzzy_start, ft.location.nofuzzy_end)
                ft_end = max(ft.location.nofuzzy_start, ft.location.nofuzzy_end)

                # cds_description = f'{r_description}||{ft_start}:{ft_end}||{cds_name}'
                cds_description = f'{r_description}||{cds_name}'
                cds_description = cds_description.replace(' ', '_')

                if len(cds_nt) > len(cds_aa) * 3 + 3:
                    print(cds_description)
                    continue

                cds_description_temp = cds_description

                cds_descriptions_temp.append(cds_description)
                cds_seqs_temp.append(cds_nt)
                i = cds_descriptions_temp.count(cds_description)
                cds_description += '|' + str(i)

                seq_id = f'{r_description}||{cds_nt}'
                seq_ids.append(seq_id)
                j = seq_ids.count(seq_id)
                cds_description += '|' + str(j)

                if j > 1:
                    continue

                # if i > 1:
                #     idx = cds_descriptions_temp.index(cds_description_temp)
                #     cds_nt_old = cds_seqs_temp[idx]
                #     if len(cds_nt_old) > len(cds_nt):
                #         if cds_nt in cds_nt_old:
                #             cds_nt = cds_nt_old

                cds_descriptions.append(cds_description)

                with open(f'{out_file_base}{r_organelle}.fasta', 'a') as _:
                    for ln in ('>' + cds_description, cds_nt):
                        _.write(f'{ln}\n')

    cds_descriptions_unique = set(cds_descriptions)
    cds_description_count = len(cds_descriptions)
    cds_description_count_unique = len(cds_descriptions_unique)
    # print(cds_description_count - cds_description_count_unique)

    cds_description_duplicates = cds_descriptions.copy()
    for x in cds_descriptions_unique:
        cds_description_duplicates.remove(x)

    assert cds_description_count == cds_description_count_unique

    print(f'Done {tax_def_name}\n')

from Bio import SeqIO

records = SeqIO.parse('GenBank/short_list.gb', 'gb')

for record in SeqIO.parse('GenBank/short_list.gb', 'gb'):
    if 'mitochond' in record.description:
        gene_source = 'mt'
    else:
        gene_source = 'cl'
    rec_name = record.description.split('complete')
    rec_name_short = rec_name[0] + gene_source
    for feature in record.features:
        if feature.type == 'CDS':
            gene_seq = ''
            for part in feature.location.parts:
                start = part.start
                end = part.end
                if part.strand == -1:
                    gene_seq += str(record.seq[start:end].reverse_complement()).lower()
                elif part.strand == 1:
                    gene_seq += str(record.seq[start:end]).lower()

            if feature.qualifiers.get('gene'):
                gene_name_type = 'gene'
                gene_name = feature.qualifiers.get('gene')[0]
            elif not feature.qualifiers.get('gene'):
                if feature.qualifiers.get('protein_id'):
                    gene_name_type = 'protein_id'
                    gene_name = feature.qualifiers.get('protein_id')[0]
            gene_desc = f'{rec_name_short} {gene_name} {gene_name_type} {feature.location.start}-{feature.location.end}'.replace(' ','_').replace(',','').replace('>','GT|').replace('<','LT|').lower()
            # Lt = Less Than <  GT = Greater Than >
            gene_final = f'>{gene_desc}'
            gene_list = [gene_final, gene_seq, '']

            with open(f'genes/{gene_source}.fa', 'a') as gene_file:
                for i in gene_list:
                    gene_file.write(i + '\n')
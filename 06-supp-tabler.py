import pandas as pd

oldcsv = pd.read_csv('overlap.csv', header=None)
taxcsv = pd.read_csv('Mitochondria_Dataset1_MatrixOccupanceclades.tsv', sep='\t')
f = pd.DataFrame([['Name', ' TaxID', 'Mito Accession', ' Seq Length', 'Chloro Accession', ' Seq Length']])

for i in range(1,len(oldcsv[1])):
    a = str(oldcsv[1][i]).strip("[]").replace("'",'')
    b = a.split(',')
    c = str(oldcsv[2][i]).strip("[]").replace("'",'')
    d = c.split(',')
    e = pd.DataFrame([[b[0],b[3],b[1],b[4],d[1],d[4]]])
    f = pd.concat([f,e], ignore_index=True)

# f.to_csv('refseq.csv', index=False, header=False)

g = []
h = []

for x in range(len(taxcsv)):
    g.append(taxcsv['taxid'][x].tolist())

for y in range(1,len(f)):
    if all(z == int(f[1].iloc[y]) for z in g):
        h.append(y)

f.drop(f.index[h], inplace=True)
i = f.reset_index(drop=True)

i.to_csv('table.csv', index=False, header=False)
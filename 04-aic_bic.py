'''
Inputs: Maximum Likelihood from MT_Partioned_sp.iqtree, MT_Partioned_spp.iqtree,
        CP_Partioned_sp.iqtree, CPandMT_Partioned_sp.iqtree

Outputs: Excel file (aic_bic.xlsx) containing AIC, AICc, BIC scores

Maximum Likelihood Data:
    Taxa: 226
    Genes = 79

    MT_Partioned_sp.iqtree (MT_SP)
        (ll) Log-likelihood:    -624371.5275
        (k) parameters:         14847
        (n) total sites:        58295
        (AIC) score:            1278437.0550
        (AICc) score:           1288584.9711
        (BIC) score:            1411663.2185
    (ll, k, n)
    (-655263.1180, 877, 58295)

    MT_Partioned_spp.iqtree (MT_SPP)
        (ll) Log-likelihood:    -655263.1180
        (k) parameters:         877
        (n) total sites:        58295
        (AIC) score:            1312280.2361
        (AICc) score:           1312307.0576
        (BIC) score:            1320149.7953
    (ll, k, n)
    (-655263.1180, 877, 58295)

    CP_Partioned_sp.iqtree (CP_SPP)
        (ll) Log-likelihood:    -2378871.5762
        (k) parameters:         1317
        (n) total sites:        103806
        (AIC) score:            4760377.1525
        (AICc) score:           4760411.0258
        (BIC) score:            4772954.8700
    (ll, k, n)
    (-2378871.5762, 1317, 103806)

    CPandMT_Partioned_sp.iqtree (CPandMT_SPP)
        (ll) Log-likelihood:    -3059725.7470
        (k) parameters:         1628
        (n) total sites:        162101
        (AIC) score:            6122707.4940
        (AICc) score:           6122740.5466
        (BIC) score:            6138980.9411
    (ll, k, n)
    (-3059725.7470, 1628, 162101)
'''

import math
import pandas as pd

                                            #(ll, k, n)
MT_SPP = (-655263.1180, 877, 58295)         #Mitochondrian SPP
MT_SP = (-624371.5275, 14847, 58295)        #Mitochondrian SP
CP_SPP = (-2378871.5762, 1317, 103806)      #Chloroplast_SPP
MTCP_SPP = (-3059725.7470, 1628, 162101)    #Mito + Chloro
MUCP = (MT_SP[0] + CP_SPP[0], MT_SP[1] + CP_SPP[1], MT_SP[2] + CP_SPP[2]) #MT_SP + CP_SPP
MPCP = (MT_SPP[0] + CP_SPP[0], MT_SPP[1] + CP_SPP[1], MT_SPP[2] + CP_SPP[2]) #MT_SPP + CP_SPP


ORG_DATA = [MT_SP, MT_SPP, CP_SPP, MTCP_SPP, MUCP, MPCP]    #list containing data
COLUMN_NAMES = ['MT_SP', 'MT_SPP', 'CP_SPP', 'MTCP_SPP', 'MT_SP + CP_SPP', 'MT_SPP + CP_SPP']  #column names
INDEX_NAMES = ['AIC', 'AICc', 'BIC']    #row names

df = pd.DataFrame(columns=COLUMN_NAMES, index=INDEX_NAMES)  #create empty data frame, labeled columns/rows

def func(name, org):    #calculations of AIC, AICc, BIC from ORG_DATA, fill empty dataframe
    aic = (2 * org[1]) - (2 * org[0])						# 2k - 2ll
    aicc = aic + (((2 * org[1]) * (org[1] + 1)) / (org[2] - org[1] - 1))	# aic + ((2k(k + 1)) / (n - k - 1))
    bic = (org[1] * math.log(org[2])) - (2 * org[0])				# (k * ln(n)) - 2ll
    df[name] = [aic,aicc,bic]
    return

for i in range(len(ORG_DATA)):      #iterate through ORG_DATA list
    func(COLUMN_NAMES[i], ORG_DATA[i])

df.to_excel('aic_bic.xlsx')     #save dataframe as excel file

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# BL microRNA barplot, percent EBV (Fig 1a?)
s = pd.read_csv('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/\
    Papers/The_EBV_microRNA_targetome/data_for_figures/figure1/A/blood871418-suppl2.csv', index_col=0)  # metadata from original BL paper (blood 2019)

m = pd.read_table('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/\
    Papers/The_EBV_microRNA_targetome/data_for_figures/figure1/A/bl_mir_cpm_primary_tumors_only.tsv', index_col=0)  # BL c.p.m. microRNA

d = {'pos':[],'neg':[]}

for i,j in zip(s['Patient barcode'], s['EBV status']):
    if 'pos' in j:
        d['pos'].append(i)
    else:
        d['neg'].append(i)

m2 = m[list(set(d['pos']) & set(m.columns))]

m2 = m2[np.sum(m2.loc[[i for i in m2.index if 'ebv' in i]]).sort_values().index]

sums = np.sum(m2.loc[[i for i in m2.index if 'ebv' in i]]) / 10000  # convert c.p.m. to c.p.hundred for percentage plotting

fig = plt.figure(figsize=(24,5), dpi=300)
ax = plt.subplot()
ax2 = ax.twinx()
ax2.bar([i-.2 for i in range(len(m2.columns))], sums, facecolor='#88D498', ec='k', linewidth=.5)
ax2.set_yticks(range(0, 101, 25))
ax.set_xlim([-1, 68.5])
ax.set_yticks([])
ax.set_xticks([])
ax2.grid(ls='--',alpha=.5)
ax2.set_axisbelow(True)
plt.savefig('bl_sum_ebv_mirs.svg',dpi=300)
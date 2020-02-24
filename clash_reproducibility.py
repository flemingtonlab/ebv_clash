import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations 
from scipy import stats



def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i or 'MTRNR' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    return df

def get_counts(path):
    df = process_df(path)
    df = pd.DataFrame(df.groupby(['mrna', 'mir'])['count'].sum()).reset_index()
    df.index = [f'{i}_{j}' for i,j in zip(df['mir'], df['mrna'])]
    return df

snu_mrna_clash_paths = [

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/ebv/SNU719_CLASH_1/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_1_1_comp_human_ebv.clash'
    '.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/ebv/SNU719_CLASH_2/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_2_1_comp_human_ebv.clash'
    '.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/ebv/SNU719_CLASH_3/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_3_1_comp_human_ebv.clash'
    '.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

]

akata_mrna_clash_paths = [

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/ebv/Akata_Latency1/Latency1_R1_comp_human_ebv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/ebv/Akata_Latency2/Latency2_R1_comp_human_ebv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/ebv/Akata_Latency3/Latency3_R1_comp_human_ebv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated')

]

hybrids = set()
for path in akata_mrna_clash_paths:
    counts = get_counts(path)
    hybrids = hybrids | set(counts.index)

akata = pd.DataFrame(index=hybrids)
for path in akata_mrna_clash_paths:
    counts = get_counts(path)
    akata[path.split('/')[-2]] = counts['count']


akata = akata.fillna(0)
akata['mean'] = np.mean(akata, 1)
akata = akata.sort_values('mean')

hybrids = set()
for path in snu_mrna_clash_paths:
    counts = get_counts(path)
    hybrids = hybrids | set(counts.index)

snu = pd.DataFrame(index=hybrids)
for path in snu_mrna_clash_paths:
    counts = get_counts(path)
    snu[path.split('/')[-2]] = counts['count']


snu = snu.fillna(0)
snu['mean'] = np.mean(snu, 1)
snu = snu.sort_values('mean')

for x, y in combinations(snu.columns[:3], 2):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(snu[x], snu[y], s=1, c='k', alpha=.3)
    ax.set_ylim([1, 100000])
    ax.set_xlim([1, 100000])
    ax.set_xscale('log')
    ax.set_yscale('log')
    pearson, p = stats.pearsonr(snu[x], snu[y])
    plt.savefig(f'x-{x}_y-{y}_pearson-{pearson}_p-{p}.png', dpi=500)

for x, y in combinations(akata.columns[:3], 2):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(akata[x], akata[y], s=1, c='k', alpha=.3)
    ax.set_ylim([10, 100000])
    ax.set_xlim([10, 100000])
    ax.set_xscale('log')
    ax.set_yscale('log')
    pearson, p = stats.pearsonr(akata[x], akata[y])
    plt.savefig(f'x-{x}_y-{y}_pearson-{pearson}_p-{p}.png', dpi=500)

snu.to_csv('snu719.mir_mrna.counts.tsv', sep='\t')
akata.to_csv('akata.mir_mrna.counts.tsv', sep='\t')



import pandas as pd
import numpy as np
from scipy.stats import ttest_rel

path = '/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/alignments/ebv/SNU719_CLASH_2/180306_I283_FCHNYYGBBXX_L1_SNU719_CLASH_2_1_comp_human_ebv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'

def process_df(path, vir='ebv'):

    df = pd.read_table(path, index_col=0)
    blacklist = set([i for i in df['mrna'] if 'MT-' in i]) | set(['SLMAP','HASPIN','VMP1','MTRNR2L12'])
    df = df[~df['mrna'].isin(blacklist)]
    df['species'] = df['mir'].map(lambda x: vir if vir in x else 'host')
    df['count'] = 1000000 * df['count'] / np.sum(df['count'])
    return df


def count_unique_interactions(df, min_counts=10):

    summed = df.groupby(['mir', 'mrna'])['count'].sum().reset_index()
    uniques = summed[summed['count'] > min_counts].groupby(['mir'])[['mrna']].nunique().sort_values('mrna')
    return uniques


def count_major_interactions_by_species(df, vir='ebv', min_counts=10):

    virus = df[df['species'] == vir].groupby(['mir', 'mrna'])['count'].sum().reset_index()
    virus = set(virus[virus['count'] > min_counts]['mrna'])
    host = df[df['species'] == 'host'].groupby(['mir', 'mrna'])['count'].sum().reset_index()
    host = set(host[host['count'] > min_counts]['mrna'])
    
    return virus, host


def calculate_average_binding_energy_by_species_weighted(path, vir='ebv'):

    df = process_df(path, vir=vir)
    virus = df[df['species'] == vir].groupby(['mir','mrna','mrna_feature_start']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
    host = df[df['species'] == 'host'].groupby(['mir','mrna','mrna_feature_start']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()

    vir_binding_list = []
    for count, binding_energy in zip(virus['count'], virus['binding_energy']):
        for n in range(round(count)):
            vir_binding_list.append(binding_energy)

    host_binding_list = []
    for count, binding_energy in zip(host['count'], host['binding_energy']):
        for n in range(round(count)):
            host_binding_list.append(binding_energy)

    return np.array(vir_binding_list), np.array(host_binding_list)



def calculate_average_binding_energy_by_species_unweighted(path, vir='ebv', min_counts=10):

    df = process_df(path, vir=vir)
    virus = df[df['species'] == vir].groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
    virus_filtered = virus[virus['count'] >= min_counts]
    host = df[df['species'] == 'host'].groupby(['mir','mrna']).agg({'count':sum, 'binding_energy':min}).sort_values('count').reset_index()
    host_filtered = host[host['count'] >= min_counts]

    return np.array(virus_filtered['binding_energy']), np.array(host_filtered['binding_energy'])

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


mghv_mrna_clash_paths = [

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/mhv68/SRR8395242/SRR8395242_1_comp_mouse_mhv68.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/mhv68/SRR8395243/SRR8395243_1_comp_mouse_mhv68.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/mhv68/SRR8395244/SRR8395244_1_comp_mouse_mhv68.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated')

]


kshv_mrna_clash_paths = [
    
    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/kshv/SRR5876950/SRR5876950_1_comp_human_kshv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),

    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/kshv/SRR5876951/SRR5876951_1_comp_human_kshv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated'),
 
    ('/Volumes/Flemington_Lab_Documents/20_lab_member_directories/3_Nate/Projects/The_EBV_microRNA_targetome/Clash_data/'
    'alignments/kshv/SRR5876952/SRR5876952_1_comp_human_kshv.clash.blast.natescript.hyb.mir92hg_converted.tsv.annotated')
    
    ]





earr, harr = [], []
for path in snu_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='ebv', min_counts=50)
    print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))
print(ttest_rel(earr, harr))

earr, harr = [], []
for path in snu_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_weighted(path, vir='ebv')
    print('ebv:', np.mean(e), '    hsa: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))
print(ttest_rel(earr, harr))



from scipy.stats import ttest_rel
earr, harr = [], []
for path in mghv_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='mghv', min_counts=500)
    print('mghv:', np.mean(e), '    mmu: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))
print(ttest_rel(earr, harr))


earr, harr = [], []
for path in mghv_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_weighted(path, vir='mghv')
    print('mghv:', np.mean(e), '    mmu: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))
print(ttest_rel(earr, harr))


earr, harr = [], []
for path in kshv_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_unweighted(path, vir='kshv', min_counts=100)
    print('kshv:', np.mean(e), '    hsa: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))
print(ttest_rel(earr, harr))


earr, harr = [], []
for path in kshv_mrna_clash_paths:

    e,h = calculate_average_binding_energy_by_species_weighted(path, vir='kshv')
    print('kshv:', np.mean(e), '    hsa: ', np.mean(h))
    earr.append(np.mean(e))
    harr.append(np.mean(h))
print(ttest_rel(earr, harr))




def count_seed_match(path, vir='ebv', min_counts=10):
    
    
    df = process_df(path, vir=vir)
    df = df[df['count'] >= min_counts]
    df['a1'] = df['a1'].map(lambda x: 0 if x=='unknown' or x=='0' else 1)
    df = df.groupby(['mir','mrna','a1','m2_5', 'm6_7', 'm8'])['count'].sum().sort_values().reset_index()
    df['seed'] = np.sum(df[['a1', 'm2_5', 'm6_7', 'm8']], 1)
    df = df.groupby(['mir', 'mrna']).agg({'count': sum, 'seed': max}).reset_index()
    df['species']= df['mir'].map(lambda x: vir if vir in x else 'host')
    df = df.groupby(['species', 'seed'])['count'].count()
    df = df.reset_index().set_index('species')
    virus, host = {}, {}
    vir_total = np.sum(df.loc[vir, 'count'])
    host_total = np.sum(df.loc['host', 'count'])
    virus['no_seed'] = df[df['seed']<6].loc[vir]['count'].sum() / vir_total
    host['no_seed'] = df[df['seed']<6].loc['host']['count'].sum() / host_total
    for i in range(6, 9):
        try:
            virus[f'{i}mer'] = df[df['seed'] == i].loc[vir]['count'] / vir_total
        except KeyError:
            virus[f'{i}mer'] = 0
        try:
            host[f'{i}mer'] = df[df['seed'] == i].loc['host']['count'] / host_total
        except KeyError:
            host[f'{i}mer'] = 0

    return virus, host

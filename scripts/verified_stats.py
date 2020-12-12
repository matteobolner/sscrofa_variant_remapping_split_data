import pandas as pd
import numpy as np
import requests, sys
import json
#def get_gene_coords(input_df, output_file):
df = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/verified_vars/ENSSSCG00000000038.verified')
complete_df = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/complete_coords/ENSSSCG00000000038.complete')
updated_ids = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/remapped_ids.csv', sep = ',')

complete_df
remapped_list = []
for i,j in zip(df['seq_v10'], df['seq_v11']):
    if i == j :
        remapped_list.append('True')
    else:
        seq_id = 0
        for b,c in zip(i,j):
            if b == c:
                seq_id += 1
            else:
                seq_id += 0
        if seq_id >= 5 and i[5]==j[5]:
            remapped_list.append('True')
        elif seq_id >= 6 and i[5]!=j[5]:
            remapped_list.append('diff_central_base')
        elif seq_id <=5 and i[5] !=j[5]:
            remapped_list.append('False')
        else:
            remapped_list.append(seq_id)
df['var_correctly_remapped'] = remapped_list
df = df.drop(columns=['chromosome', 'are_segments_identical', 'are_the_central_bases_identical', 'strand'])
df

merged_df = df.merge(complete_df, left_on='var_pos_10', right_on='final_var_pos')
merged_df = merged_df.drop(columns=['chromosome', 'start', 'end', 'final_var_pos'])

merged_df

temp_df = updated_ids.loc[updated_ids['id_10_2'] == merged_df['ensembl_id'][0]]
temp_df
updated_id = temp_df['updated_id'].iloc[0]
updated_id
merged_df['ensembl_id_11'] = updated_id
if updated_id != 'NOT_FOUND':
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/"+ updated_id + "?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json()
    start_11 = decoded['start']
    end_11 = decoded['end']
    name_11 = decoded['display_name']
    desc_11 = decoded['description']
merged_df['gene_start'] = start_11
merged_df['gene_end'] = end_11
merged_df['gene_name_11'] = name_11

info_data = {'gene_name': [name_11], 'ensembl_id' : [updated_id], 'chromosome' : [merged_df['chromosome_11'].iloc[0]], 'gene_start': [start_11], 'gene_end': [end_11], 'strand': [merged_df['strand'].iloc[0]], 'description': [desc_11]}
info_df = pd.DataFrame(info_data, columns =['gene_name', 'ensembl_id', 'chromosome', 'gene_start', 'gene_end', 'strand', 'description'])
info_df
df = pd.DataFrame (data, columns = ['First Column Name','Second Column Name',...])
merged_df
get_stats(snakemake.input[0], snakemake.output[0])

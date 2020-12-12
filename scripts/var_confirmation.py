import pandas as pd
import numpy as np
import requests
import sys

def update_ids(remapped_file):
    #remapped_df = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/remapped_vars/ENSSSCG00000016940.remapped', sep = ',')
    remapped_df = pd.read_csv(remapped_file, sep = ',')
    updated_ids = pd.read_csv('/home/pelmo/sscrofa_variant_remapping/data/remapped_ids.csv', sep = ',')
    temp_df = updated_ids.loc[updated_ids['id_10_2'] == remapped_df['ensembl_id_10'][0]]
    updated_id = temp_df['updated_id'].iloc[0]
    remapped_df['ensembl_id_11'] = updated_id
    if updated_id != 'NOT_FOUND':
        server = "https://rest.ensembl.org"
        ext = "/lookup/id/"+ updated_id + "?"
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
          r.raise_for_status()
          sys.exit()
        decoded = r.json()
        chr_11 = decoded['seq_region_name']
        remapped_df['chromosome_11'] = chr_11
        remapped_df
    else:
        remapped_df['chromosome_11'] = 'NOT_FOUND'
    if remapped_df['chromosome_11'][0] == 'AEMK02000679.1':
        remapped_df['chromosome_11'] = '3'
    return(remapped_df)



def verify_vars(updated_df, output_file):
    #input_file = "/home/pelmo/sscrofa_variant_remapping/data/remapped_vars/ENSSSCG00000031069.remapped"
    df = updated_df
    df = df.dropna(subset = ['var_pos_11'])
    df['var_pos_11'] = df['var_pos_11'].astype(int)
    df = df[['chromosome', 'chromosome_11', 'var_pos_10', 'var_pos_11', 'strand']]
    df['coordinates_v10'] = df['chromosome'].astype(str) + ":" + df['var_pos_10'].astype(str) + ".." + df['var_pos_10'].astype(str) + ":" + df['strand'].astype(str)
    if df['chromosome_11'][0] != "NOT_FOUND":
        df['coordinates_v11'] = df['chromosome_11'].astype(str) + ":" + df['var_pos_11'].astype(str) + ".." + df['var_pos_11'].astype(str) + ":" + df['strand'].astype(str)

        server_v10 = "https://may2017.rest.ensembl.org"
        server_v11 = "https://rest.ensembl.org"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}

        decoded_dict_v10 = {}
        decoded_dict_v11 = {}

        segment_comparison_list = []
        base_comparison_list = []

        for i,j in zip(df['coordinates_v10'], df['coordinates_v11']):
            ext_v10 = "/sequence/region/sscrofa/" + i + "?expand_3prime=5;expand_5prime=5"
            ext_v11 = "/sequence/region/sscrofa/" + j + "?expand_3prime=5;expand_5prime=5"

            r_v10 = requests.get(server_v10 + ext_v10, headers={ "Content-Type" : "application/json"})
            r_v11 = requests.get(server_v11 + ext_v11, headers={ "Content-Type" : "application/json"})
            if not r_v10.ok:
              r_v10.raise_for_status()
              sys.exit()

            if not r_v11.ok:
              r_v11.raise_for_status()
              sys.exit()

            decoded_v10 = r_v10.json()
            decoded_v11 = r_v11.json()

            decoded_dict_v10[i] = decoded_v10['seq']
            decoded_dict_v11[j] = decoded_v11['seq']
            #temp_seq_v10 = Seq(decoded_v10['seq'])
            #temp_seq_v11 = Seq(decoded_v11['seq'])
            #print(temp_seq_v11.find(decoded_v10['seq']))
            if decoded_v10['seq'] == decoded_v11['seq']:
                segment_comparison_list.append("True")
            else:
                segment_comparison_list.append("False")

            if decoded_v10['seq'][5] == decoded_v11['seq'][5]:
                base_comparison_list.append("True")
            else:
                base_comparison_list.append("False")
        seq_df_v10 = pd.DataFrame.from_dict(decoded_dict_v10, orient='index', columns = ['seq_v10'])
        seq_df_v11 = pd.DataFrame.from_dict(decoded_dict_v11, orient='index', columns = ['seq_v11'])

        merged_dfs_first_part = df.merge(seq_df_v10, left_on = 'coordinates_v10', right_index =True )
        merged_dfs = merged_dfs_first_part.merge(seq_df_v11, left_on = 'coordinates_v11', right_index =True )

        merged_dfs.drop(['coordinates_v10', 'coordinates_v11'], axis =1)
        merged_dfs = merged_dfs[['chromosome', 'chromosome_11', 'strand', 'var_pos_10', 'var_pos_11', 'seq_v10', 'seq_v11']]

        merged_dfs['are_segments_identical'] = segment_comparison_list
        merged_dfs['are_the_central_bases_identical'] = base_comparison_list
        merged_dfs.to_csv(output_file, index = False)
    else:
        df.to_csv(output_file, index = False)
    return()


verify_vars(update_ids(snakemake.input[0]), snakemake.output[0])

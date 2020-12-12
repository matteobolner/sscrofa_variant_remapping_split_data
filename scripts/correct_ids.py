import pandas as pd
remap_df = pd.read_csv("/home/pelmo/sscrofa_variant_remapping/data/id_remap.txt")
remap_df['updated_id'] = "NOID"
for i,row in remap_df.iterrows():
    if row['id_11_1'] != row['id_11_1']:
        remap_df.at[i, 'updated_id'] = row['id_10_2']
    elif row['id_11_1'] == 'NOT FOUND':
        remap_df.at[i, 'updated_id'] = str(row['id_10_2']) + "_NOT_UPDATED"
    else:
        remap_df.at[i, 'updated_id'] = row['id_11_1']


remap_df.to_csv("/home/pelmo/sscrofa_variant_remapping/data/remapped_ids.csv")

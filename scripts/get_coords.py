import pandas as pd
import os
import requests
import json

server = "https://may2017.rest.ensembl.org"

def get_coords(starting_file, gene_file, vars_file, complete_coords_file):
    starting_df = pd.read_csv(starting_file, sep='|', names=['gene_name','ensembl_id','chromosome','start','end','var_pos'])
    starting_df = starting_df.drop_duplicates()
    gene_coords = starting_df.head(n=1)[['chromosome','start','end']]
    ensembl_id = starting_df['ensembl_id'].iloc[0]
    ext = "/lookup/id/" + ensembl_id + "?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    strand = decoded['strand']
    if strand == +1:
        starting_df['final_var_pos'] = (starting_df['start'] + starting_df['var_pos']) - 1
    elif strand == -1:
        starting_df['final_var_pos'] = (starting_df['end'] - starting_df['var_pos']) + 1

    starting_df['strand'] = strand

    gene_coords = gene_coords['chromosome'].astype(str) + ':' + gene_coords['start'].astype(str) + '-' + gene_coords['end'].astype(str)
    vars_coords = starting_df['chromosome'].astype(str) + ":" + starting_df['final_var_pos'].astype(str)+ "-" + starting_df['final_var_pos'].astype(str)

    gene_coords.to_csv(gene_file, header = False, index = False)
    vars_coords.to_csv(vars_file, header = False, index = False)
    starting_df.to_csv(complete_coords_file, index=False, sep = ',')
    return()

get_coords(snakemake.input[0], snakemake.output.gene_files, snakemake.output.var_files, snakemake.output.complete_coords)

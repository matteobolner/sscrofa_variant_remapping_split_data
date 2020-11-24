import pandas as pd
import numpy as np
def merge_results(report_file, initial_coords_file, remapped_vars_file, unmapped_vars_file):
    report_df = pd.read_csv(report_file, sep='\t')
    complete_df = pd.read_csv(initial_coords_file)
    report_df_trimmed = report_df[['source_start', 'mapped_start']]
    report_df_trimmed.columns = ['var_pos_10', 'var_pos_11']
    final_df = complete_df.merge(report_df_trimmed, left_on='final_var_pos', right_on='var_pos_10')
    final_df = final_df.drop(columns=['final_var_pos'])
    final_df.columns = ['gene_name_10','ensembl_id_10', 'chromosome', 'gene_start_10', 'gene_end_10', 'relative_var_pos_10', 'strand', 'var_pos_10', 'var_pos_11']
    final_df = final_df.replace(r'^\s*$', np.nan, regex=True)
    unmapped_df = final_df[final_df['var_pos_11'].isna()]
    unmapped_df.to_csv(unmapped_vars_file, index=False, na_rep = 'NO_MAPPING')
    final_df.to_csv(remapped_vars_file, index = False, na_rep = 'NaN')

    return()

merge_results(snakemake.input[0], snakemake.input[1], snakemake.output[0], snakemake.output[1])

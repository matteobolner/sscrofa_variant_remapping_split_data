IDS, = glob_wildcards("data/starting_files/{id}")

rule all:
    input:
        #expand("data/gene_coords_remapped/{id}.gene.report", id=IDS),
        #expand("data/gene_coords_remapped/{id}.gene.annot", id=IDS)
        #expand("data/gene_coords/{id}.gene", id=IDS),
        #expand("data/var_coords/{id}.vars", id=IDS)
        #expand("data/complete_coords/{id}.complete", id=IDS)
        expand("data/remapped_vars/{id}.remapped", id=IDS),


rule get_coords:
    input:
        "data/starting_files/{id}"

    output:
        gene_files = "data/gene_coords/{id}.gene",
        var_files = "data/var_coords/{id}.vars",
        complete_coords = "data/complete_coords/{id}.complete"
    script:
        "scripts/get_coords.py"

rule merge_results:
    input:
        "data/var_coords_remapped_report/{id}.vars.report",
        "data/complete_coords/{id}.complete"
    output:
        remapped_vars_file = "data/remapped_vars/{id}.remapped",
        unmapped_vars_file = "data/unmapped_vars/{id}.unmapped"
    script:
        "scripts/merge_results.py"



'''rule remap_genes:
    input:
        "data/var_coords/{id}.vars"
    output:
        annot = "data/var_coords_remapped/{id}.vars.annot",
        report = "data/var_coords_remapped/{id}.vars.report"
    shell:
        "bash scripts/remap_vars.sh {input} {output.annot} {output.report}"'''

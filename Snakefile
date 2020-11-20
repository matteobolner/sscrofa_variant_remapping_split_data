IDS, = glob_wildcards("data/starting_files/{id}")

rule all:
    input:
        expand("data/gene_coords_remapped/{id}.gene.report", id=IDS),
        expand("data/gene_coords_remapped/{id}.gene.annot", id=IDS)
        #expand("data/gene_coords/{id}.gene", id=IDS),
        #expand("data/var_coords/{id}.vars", id=IDS)

rule get_coords:
    input:
        "data/starting_files/{id}"

    output:
        gene_files = "data/gene_coords/{id}.gene",
        var_files = "data/var_coords/{id}.vars",
        complete_coords = "data/complete_coords/{id}.complete"
    script:
        "scripts/get_coords.py"



rule remap_genes:
    input:
        "data/gene_coords/{id}.gene"
    output:
        annot = "data/gene_coords_remapped/{id}.gene.annot",
        report = "data/gene_coords_remapped/{id}.gene.report"
    shell:
        "bash scripts/remap_gene.sh {input} {output.annot} {output.report}"

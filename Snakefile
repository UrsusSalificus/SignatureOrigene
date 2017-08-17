
'''
rule setup:
    shell:
        'bash scripts/setup'
'''

SPECIES, = glob_wildcards("config/species/{any}")

rule all:
    input:
        expand("data/genomes/{species}_genomes.fna", species = SPECIES),
        expand("data/genomes/{species}_feature_table.fna", species = SPECIES)


rule download:
    input:
        "config/species/{species}"
    output:
        "data/genomes/{species}_genomes.fna",
        "data/genomes/{species}_feature_table.fna"
    shell:
        """
        no_path=$( basename {input} )
        bash download_genomes.sh $no_path {output}
        """

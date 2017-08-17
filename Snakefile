SPECIES, = glob_wildcards("config/species/{any}")
WINDOWS, = glob_wildcards("config/windows/{any}")
GS, = glob_wildcards("config/gs/{any}")
KMER, = glob_wildcards("config/kmer/{any}")

if "FCGRs" in GS:
    rule all:
        input:
            expand("data/following/FCGRs/{windows}/{species}_done.txt",
                windows = WINDOWS, species = SPECIES),
            expand("data/genomes/{species}_feature_table.fna", species = SPECIES)
else:
    rule all:
        input:
            expand("data/following/CGRs/{windows}/{species}_done.txt",
                windows = WINDOWS, species = SPECIES),
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
        bash scripts/download_genomes.sh $no_path {output}
        """

rule CGR:
    input:
        expand("data/genomes/{species}_genomes.fna", species = SPECIES),
        expand("config/windows/{windows}", windows = WINDOWS)
    output:
        "data/following/CGRs/{windows}/{species}_done.txt"
    threads:
        3
    shell:
        "python3 scripts/windowed_CGR.py {input} {threads} {output}"


if "FCGRs" in GS:
    rule FCGR:
        input:
            expand("data/following/CGRs/{windows}/{species}_done.txt",
                windows = WINDOWS, species = SPECIES),
            expand("data/genomes/{species}_genomes.fna", species = SPECIES),
            expand("config/windows/{windows}", windows = WINDOWS),
            expand("config/kmer/{kmer}", kmer = KMER)
        output:
            "data/following/FCGRs/{windows}/{species}_done.txt"
        threads:
            3
        shell:
            "python3 scripts/windowed_FCGR.py {input} {threads} {output}"




######

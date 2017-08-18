SPECIES, = glob_wildcards("config/species/{any}")
WINDOWS, = glob_wildcards("config/windows/{any}")
GS, = glob_wildcards("config/gs/{any}")
KMER, = glob_wildcards("config/kmer/{any}")
FEATURES, = glob_wildcards("config/features/{any}")

if "FCGRs" in GS:
    if "CDS" in FEATURES:
        rule FCGR_CDS:
            input:
                expand("data/following/FCGRs/{windows}/{species}_done.txt",
                    windows = WINDOWS, species = SPECIES),
                expand("files/features/{windows}/{species}_CDS.txt",
                    windows = WINDOWS, species = SPECIES)
    else:
        rule FCGR_LCR:
            input:
                expand("data/following/FCGRs/{windows}/{species}_done.txt",
                    windows = WINDOWS, species = SPECIES),
if "DFTs" in GS:
    if "CDS" in FEATURES:
        rule DFT_CDS:
            input:
                expand("data/following/DFTs/{windows}/{species}_done.txt",
                    windows = WINDOWS, species = SPECIES),
                expand("files/features/{windows}/{species}_CDS.txt",
                    windows = WINDOWS, species = SPECIES)
    else:
        rule DFT_LCR:
            input:
                expand("data/following/CGRs/{windows}/{species}_done.txt",
                    windows = WINDOWS, species = SPECIES)

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

if "CDS" in FEATURES:
    rule CDS_extract:
        input:
            expand("data/genomes/{species}_genomes.fna", species = SPECIES),
            expand("data/genomes/{species}_feature_table.fna", species = SPECIES),
            expand("config/windows/{windows}", windows = WINDOWS)
        output:
            "files/features/{windows}/{species}_CDS.txt"
        shell:
            "python3 scripts/extract_CDS.py {input} {output}"

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

if "DFTs" in GS:
    rule DFT:
        input:
            expand("data/following/CGRs/{windows}/{species}_done.txt",
                windows = WINDOWS, species = SPECIES),
            expand("data/genomes/{species}_genomes.fna", species = SPECIES),
            expand("config/windows/{windows}", windows = WINDOWS)
        output:
            "files/DFTs/{windows}/{species}/{species}_DFTs.txt"
        threads:
            3
        shell:
            "python3 scripts/windowed_DFT.py {input} {threads} {output}"

    rule DFT_dist:
        input:
            expand("files/DFTs/{windows}/{species}_DFTs.txt",
                windows = WINDOWS, species = SPECIES),
            expand("data/genomes/{species}_genomes.fna", species = SPECIES),
            expand("config/windows/{windows}", windows = WINDOWS)
        output:
            "files/distances/euclidean/{windows}/{species}_dist_matrix.txt"
        shell:
            "python3 scripts/euclidean_dist.py {input} {output}"








##############

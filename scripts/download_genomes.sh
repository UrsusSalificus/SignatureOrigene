#!/usr/bin/env bash
# Will download whole concatenated genomes

species='h_sapiens d_rerio c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli'
accession='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/035/GCA_000002035.4_GRCz11/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/051/215/GCA_001051215.1_ASM105121v1/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/023/665/GCF_000023665.1_ASM2366v1/'

# Must know the number of species/genomes, and where we are (script folder)
number_genomes=$( echo $species | wc -w )
way_back=$( pwd )

for each_genome in $( seq 1 1 $number_genomes ) ; do
    each_species=$( echo $species | cut -d" " -f$each_genome )
    each_accession=$( echo $accession | cut -d" " -f$each_genome )
    out_folder=$'../data/genomes/'$each_species
    if [ ! -d "$out_folder" ]; then
        mkdir -p $out_folder
    fi
    cd $out_folder

    wget -r -np -nd -A "genomic.fna.gz" $each_accession

    all_tar=$( find *.gz )

    # Will then decompress everything
    for each_tar in $all_tar ; do
        IFS='.' read -r -a array <<< "$each_tar"
        outfile_name=$( echo $each_species'_'${array[1]}'.'${array[2]} )
        gunzip < $each_tar > $outfile_name
    done


    cd $way_back
done


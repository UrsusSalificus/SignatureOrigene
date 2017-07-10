#!/usr/bin/env bash
# Will download whole concatenated genomes

species='h_sapiens m_musculus c_elegans d_melanogaster a_thaliana s_cerevisiae e_coli'
accession='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.36_GRCh38.p10/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.3_TAIR10/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/'

# Must know the number of species/genomes, and where we are (script folder)
number_genomes=$( echo $species | wc -w )
out_folder=$'../data/genomes/'
if [ ! -d "$out_folder" ]; then
    mkdir -p $out_folder
fi
cd $out_folder

for each_genome in $( seq 1 1 $number_genomes ) ; do
    each_species=$( echo $species | cut -d" " -f$each_genome )
    each_accession=$( echo $accession | cut -d" " -f$each_genome )

    # Download with download time stamp (--no...), recursively (-r), but only on the first directory (-l 1),
    # without ascending to parent directories (-np), and without recreating the whole parent directory hierarchy (-nd)
    # and finally looking only for a certain pattern (-A #list_of_pattern)
    wget --no-use-server-timestamps -r -l 1 -np -nd -A "genomic.fna.gz" $each_accession
    wget --no-use-server-timestamps -r -l 1 -np -nd -A "feature_table.txt.gz" $each_accession

    all_tar=$( find *.gz )

    # Will then decompress everything
    for each_tar in $all_tar ; do
        IFS='.' read -r -a array <<< "$each_tar"
        outfile_name=$( echo $each_species'_'${array[1]}'.'${array[2]} )
        gunzip < $each_tar > $outfile_name
        rm $each_tar
    done
done

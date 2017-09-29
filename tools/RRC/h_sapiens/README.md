Genetic map of Homo sapiens GRCh37 obtained from the 
[HapMap](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz).

`
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
tar -xzf genetic_map_HapMapII_GRCh37.tar.gz
rm *tar.gz
cat genetic_map* >> all_genetic_maps.txt
rm genetic_map_GRCh37_chr*
mv README.txt MAP_README.txt
`

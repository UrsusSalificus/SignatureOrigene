import CGR_functions as fn
import re
import math

# To makes up with any updates, we use the list of species from the downloading script:
species=[]
for_species='download_genomes.sh'
with open(for_species, 'r') as species_file:
    for each_line in species_file:
        if each_line.split('=')[0] == 'species':
            for each_species in each_line.split('=')[1].strip().strip('\'').split(' '):
                species.append(each_species)

# Due to non-consistent pattern of file name, whole genome ('genomic') is used in multiple
# fasta files (CDS or RNA only, and any nucleotides) names. We must thus reconstruct the exact path.
species_genomes=[]
for each_species in species:
    all_files = fn.extract_path('../data/genomes/', str(each_species + '*'))
    index = all_files[0].split('/')[3].split('_')[2] + '_' + all_files[0].split('/')[3].split('_')[3]
    pattern_genome = str(each_species + '_' + index + '_genomic*')
    species_genomes.append(fn.extract_path('../data/genomes/', pattern_genome)[0])

window_size = 150000

when_removed = []
for each_path in range(len(species_genomes)):
    when_removed.append([])
    records = fn.fetch_fasta(species_genomes[each_path])
    if len(records) > 1:
        for each_chromosomes in records:
            cleaned_seq = re.sub('[Nn]', '', str(each_chromosomes.seq))
            n_windows = math.floor(len(cleaned_seq) / window_size)
            when_removed[each_path].append(n_windows)
    else:
        cleaned_seq = re.sub('[Nn]', '', str(records[0].seq))
        n_windows = math.floor(len(cleaned_seq) / window_size)
        when_removed[each_path].append(n_windows)
    print(str('Ended' + species_genomes[each_path] + '!'))

when_avoided= []
for each_path in range(len(species_genomes)):
    when_avoided.append([])
    records = fn.fetch_fasta(species_genomes[each_path])
    if len(records) > 1:
        for each_chromosomes in records:
            if len(str(each_chromosomes.seq)) > window_size:
                uncleaned_seq = str(each_chromosomes.seq)
                n_windows = math.floor(len(uncleaned_seq) / window_size)
                for each_window_start in range (0, n_windows*window_size, window_size):
                    window = uncleaned_seq[each_window_start:each_window_start+window_size]
                    if 'N' in window or 'n' in window:
                        when_avoided[each_path].append(0)
                    else:
                        when_avoided[each_path].append(1)
            else:
                when_avoided[each_path].append(0)

    else:
        uncleaned_seq = str(records[0].seq)
        n_windows = math.floor(len(uncleaned_seq) / window_size)
        for each_window_start in range(0, n_windows * window_size, window_size):
            window = uncleaned_seq[each_window_start:each_window_start + window_size]
            if 'N' in window or 'n' in window:
                when_avoided[each_path].append(0)
            else:
                when_avoided[each_path].append(1)
    print(str('Ended' + species_genomes[each_path] + '!'))

for each_species in range(len(when_avoided)):
    print(species[each_species])
    print(sum(when_removed[each_species]))
    print(sum(when_avoided[each_species]))

len(when_removed[5])
species[5]


import CGR_functions
import math
import numpy

# Beta globin cds of human is used to perform this test
#fasta_file = '../data/h_sapiens_globin.fasta'
fasta_file = '../../test2'

records = CGR_functions.fetch_fasta(fasta_file)

records_seq = str(records[0].seq)

CGR = CGR_functions.CGR_coordinates(records_seq, '')   # '' -> don't want an outfile
k_size = 8


#####################################################################################################################3
# THE FUNCTION:


# We will now search fo the boundaries of the grids by binary search
# If CGR is a string = a path
if isinstance(CGR, str):
    x_coord = []
    y_coord = []
    with open(CGR, 'r') as CGR_file:
        for each_coord in CGR_file:
            x_coord.append(each_coord.split()[0])
            y_coord.append(each_coord.split()[1])
    coordinates = [x_coord, y_coord]
# Else it can be used as it is
else:
    coordinates = CGR

# We take out the k_size-1 first coordinates, are these are the coordinates of words smaller than our
# wanted k-mer size, and would add small errors later on when counting the frequencies
for each_xy in [0,1]:
    coordinates[each_xy] = coordinates[each_xy][k_size-1::]

# We will then sort our x coordinates
sort_x = sorted(coordinates[0])


# Make sure the decimals are precise enough
decimals = int(math.pow(10, k_size + 2))
# Calculate the number of different k-mer we will compute (= number of grid we divide the CGR with)
grid_size = int(math.pow(4, k_size))

# Compute all the possible starting/ending coordinates of all these grids
start_coord_ranges = range(0, int(decimals), int(decimals / math.sqrt(grid_size)))
end_coord_ranges = range(int(decimals / math.sqrt(grid_size)), int(decimals + decimals / math.sqrt(grid_size)),
                         int(decimals / math.sqrt(grid_size)))


#########################################
# Binary search won't work -> don't search for exact match !
# Instead, can walk throuhg the coordinates (lentgh = n), and keep a mark everytime you are a boundary!
# You only have to this both x and y !
# Final O-notation = O(n) + O(n) + O(m)   where m = number of grids
#########################################


# For x coordinates:
x_boundaries = []
which_boundary = 0

for each_coordinates in range(len(sort_x)):
    # If we are a the last coordinate, we must mark any remaining boundaries as empty
    if each_coordinates == len(sort_x) - 1:
        x_boundaries.append(each_coordinates)
        while len(x_boundaries) != len(start_coord_ranges) + 1:
            x_boundaries.append('empty')
        break
    # First, we see if the i coordinates is a "boundary" of a grid (if bigger/outside of the i (which_boundary) grid)
    if sort_x[each_coordinates] >= start_coord_ranges[which_boundary] / decimals:
        # If we are just before the last grid, we cannot have problem of empty grid left:
        if which_boundary == len(start_coord_ranges)-1:
            x_boundaries.append(each_coordinates)
            which_boundary += 1
        # In some cases, the i coordinates is not only bigger than the i boundary, but also of i+1
        # (and even the next, so on and so on). In this case, we know that the i boundary
        # does not contain any coordinates ('empty'), and that we can go to the next boundary (which_boundary +1)
        elif sort_x[each_coordinates] >= start_coord_ranges[which_boundary+1] / decimals:
            while sort_x[each_coordinates] >= start_coord_ranges[which_boundary+1] / decimals:
                x_boundaries.append('empty')
                which_boundary += 1
        # If not at the before-last grid, or if only bigger to i grid, it is the boundary of the i grid,
        # we then note where we are in the index, and go for the next starting boundary (which_boundary +1)
        else:
            x_boundaries.append(each_coordinates)
            which_boundary += 1
        # If we ended up being at the last grid, the last boundary is simply the last element
        # and we can stop here (break)
        if which_boundary == len(start_coord_ranges):
            x_boundaries.append(len(sort_x))
            break



# We will need the order of indices of x coordinates to do the same operations on y:
numpy_x = numpy.array(coordinates[0])
numpy_y = numpy.array(coordinates[1])
sort_index_x = numpy.argsort(numpy_x)
y_boundaries = []

for each_column in range(len(x_boundaries)-1):
    y_boundaries.append([])

    # We will compute the y boundaries differently depending on whether the x column is empty or not:
    if not x_boundaries[each_column] == 'empty':
        # For each column of grids, we can sort the y of the sorted x, and do the same steps as before
        # First, we must find all the y corresponding to the x in the right grid:

        # In case we are at the last column
        if each_column == len(x_boundaries)-1:
            corresponding_x_indexes = [x_boundaries[each_column], len(sort_x)]

        # Now in cases next column is an empty column, need the first non-empty one to know the corresponding x
        next_column = each_column + 1
        if x_boundaries[next_column] == 'empty':
            while True:
                # In cases we have only empty columns left,
                if next_column == len(x_boundaries)-1:
                    corresponding_x_indexes = [x_boundaries[each_column], len(sort_x)]
                    break
                # Else we must find which is the next non-empty
                if x_boundaries[next_column] == 'empty':
                    next_column += 1
                else:
                    corresponding_x_indexes = [x_boundaries[each_column], x_boundaries[next_column]]
                    break
        # Otherwise, it's simply the x between the boundaries
        else:
            corresponding_x_indexes = [x_boundaries[each_column], x_boundaries[next_column]]

        # Now that we have the corresponding x indexes, we can find the y of these corresponding x, sort them
        # and do the same operation that was performed on x before
        each_column_y = numpy_y[sort_index_x[corresponding_x_indexes[0]:corresponding_x_indexes[1]]]
        sort_y = numpy.sort(each_column_y)
        which_boundary = 0
        # Now do the same operation (see comments for the x coordinates operation)
        for each_coordinates in range(len(sort_y)):
            if sort_y[each_coordinates] >= start_coord_ranges[which_boundary] / decimals:
                # Penultimate grid
                if which_boundary == len(start_coord_ranges) - 1:
                    y_boundaries[each_column].append(each_coordinates)
                    which_boundary += 1
                # Empty grid later on:
                elif sort_y[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                    while not which_boundary == len(start_coord_ranges) - 1 \
                            and sort_y[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                        y_boundaries[each_column].append('empty')
                        which_boundary += 1
                    y_boundaries[each_column].append(each_coordinates)
                    which_boundary += 1
                # Any not special grid
                else:
                    y_boundaries[each_column].append(each_coordinates)
                    which_boundary += 1
                # Ultimate grid:
                if which_boundary == len(start_coord_ranges):
                    y_boundaries[each_column].append(len(sort_y))
                    break
            # If we are at the last coordinate
            if each_coordinates == len(sort_y) - 1:
                y_boundaries[each_column].append(len(sort_y))
                while len(y_boundaries[each_column]) != len(start_coord_ranges) + 1:
                    y_boundaries[each_column].append('empty')
                break

    # Else, the column is empty, and must be marked as such:
    else:
        while len(y_boundaries[each_column]) != len(start_coord_ranges)+1:
            y_boundaries[each_column].append('empty')

FCGR = []
# We know have to use the y boundaries to count the frequencies:

for each_column in range(len(y_boundaries)):
    for each_kmer in range(len(y_boundaries[each_column])-1):
        next_kmer = each_kmer + 1
        # If empty, this mean there is no counts for this kmer.
        if y_boundaries[each_column][each_kmer] == 'empty':
            FCGR.append(0)
        # As we subtract, we must have a non-empty index to compare with:
        elif y_boundaries[each_column][each_kmer + 1] == 'empty':
            find_non_empty = True
            while y_boundaries[each_column][next_kmer] == 'empty':
                # If the last kmer is also empty it simply means that each_kmer was the last boundary
                if next_kmer == len(start_coord_ranges):
                    FCGR.append(0)
                    find_non_empty = False
                    break
                else:
                    next_kmer += 1
            # If we find non-empty boundaries in the next kmers, can use it to count the number of kmer in grid i
            # Note that we add +1 to the index, as we are not as a starting position
            if find_non_empty:
                # Note: Python start indexing at 0.
                # To easily count using subtraction we must add +1 to the next index to get the real index.
                FCGR.append( (y_boundaries[each_column][next_kmer]) - y_boundaries[each_column][each_kmer])
        # Else can simply use the next kmer index to count the number of coordinates in the grid
        else:
            # Note: Python start indexing at 0.
            # To easily count using subtraction we must add +1 to the next index to get the real index.
            FCGR.append( (y_boundaries[each_column][next_kmer]) - y_boundaries[each_column][each_kmer])

for each in range(len(FCGR)):
    if FCGR[each] != 0:
        print(each)


sum(FCGR)

for each in range(len(FCGR)):
    print(str(FCGR_labels[each] + ' : T = ' + str(lol[each]) + ' F = ' + str(FCGR[each]) ))



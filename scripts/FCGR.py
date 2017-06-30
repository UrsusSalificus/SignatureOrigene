import CGR_functions
import math
import numpy

# Beta globin cds of human is used to perform this test
# fasta_file = '../data/h_sapiens_globin.fasta'
fasta_file = '../../test2'

records = CGR_functions.fetch_fasta(fasta_file)

records_seq = str(records[0].seq)

CGR = CGR_functions.CGR_coordinates(records_seq, '')   # '' -> don't want an outfile
k_size = 2


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
    # We search the starting boundary of next grid, if we find it, we note where we are in the index,
    # and go for the next starting boundary (each_column +1)
    if sort_x[each_coordinates] >= start_coord_ranges[which_boundary] / decimals:
        # If we get at the last grid, the last boundary is simply the last element
        # and we can stop here (break)
        if which_boundary == len(start_coord_ranges):
            x_boundaries.append(len(sort_x))
            break
        if sort_x[each_coordinates] >= start_coord_ranges[which_boundary+1] / decimals:
            while sort_x[each_coordinates] >= start_coord_ranges[which_boundary+1] / decimals:
                x_boundaries.append('empty')
                which_boundary += 1
        else:
            x_boundaries.append(each_coordinates)
            which_boundary += 1
    if each_coordinates == len(sort_x) - 1:
        x_boundaries.append(each_coordinates)
        while len(x_boundaries) != len(start_coord_ranges):
            x_boundaries.append('empty')


# We will need the order of indices of x coordinates to do the same operations on y:
numpy_x = numpy.array(coordinates[0])
numpy_y = numpy.array(coordinates[1])
sort_index_x = numpy.argsort(numpy_x)
y_boundaries = []

for each_column in range(len(x_boundaries) -1):
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
        # Now do the same operation
        # TODO: finish this part
        for each_coordinates in range(len(sort_y)):
            # We search the starting boundary of next grid, if we find it, we note where we are in the index,
            # and go for the next starting boundary (each_column +1)
            if sort_y[each_coordinates] >= start_coord_ranges[which_boundary] / decimals:
                # If we get at the last grid, the last boundary is simply the last element
                # and we can stop here (break)
                if which_boundary == len(start_coord_ranges):
                    y_boundaries[each_column].append(len(sort_y))
                    break
                if sort_y[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                    while sort_x[each_coordinates] >= start_coord_ranges[which_boundary + 1] / decimals:
                        y_boundaries[each_column].append('empty')
                        which_boundary += 1
                else:
                    y_boundaries[each_column].append(each_coordinates)
                    which_boundary += 1
            if each_coordinates == len(sort_x) - 1:
                y_boundaries[each_column].append(each_coordinates)
                while len(y_boundaries[each_column]) != len(start_coord_ranges):
                    y_boundaries[each_column].append('empty')

    # Else, the column is empty, and must be marked as such:
    else:
        while len(y_boundaries[each_column]) != len(start_coord_ranges):
            y_boundaries[each_column].append(0)


    if each_column_y.any():

    else:
        while len(y_boundaries[each_column]) != len(start_coord_ranges):
            y_boundaries[each_column].append(0)

FCGR = []
# We know have to use the y boundaries to count the frequencies:
for each_column in range(len(y_boundaries)):
    for each_kmer in range(len(y_boundaries[each_column])-1):
        FCGR.append(y_boundaries[each_column][each_kmer + 1] - y_boundaries[each_column][each_kmer])

sum(FCGR)

for each in range(len(FCGR)):
    print(str(FCGR_labels[each] + ' : T = ' + str(lol[each]) + ' F = ' + str(FCGR[each]) ))



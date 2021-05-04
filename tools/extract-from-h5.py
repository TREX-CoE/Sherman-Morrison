#!/usr/bin/env python
import h5py
import sys

# Helper script to extract a few cycles from a large dataset. This will be
# especially useful for the Verificarlo CI, since vfc_ci_cycles.txt can be
# used to both extract (with this script), and read the small dataset (in a CI
# run).


    # Parse arguments

if len(sys.argv) != 4:
    sys.stderr.write(
        "Error: Wrong number of arguments. Usage : extract_h5.py "\
        "<source_dataset.hdf5> <cycles_list.txt> <destination_dataset.hdf5>\n"
    )
    exit(1)

source_dataset_path = sys.argv[1]
cycles_list_path = sys.argv[2]
destination_dataset_path = sys.argv[3]


    # Read the cycles list

cycles_list = []

try:
    f = open(cycles_list_path)
    for line in f:
        cycles_list.extend([cycle for cycle in line.split()])
except IOError:
    sys.stderr.write("Error: Could not read " + cycles_list_path + "\n")
    exit(1)
finally:
    f.close()


    # Read the source dataset, and extract the cycles to the destination dataset

try:
    fs = h5py.File(source_dataset_path, "r")

except IOError:
    sys.stderr.write("Error: Could not read " + source_dataset_path + "\n")
    exit(1)

fd = h5py.File(destination_dataset_path, "w")

# Copy cycles groups

groups = [
    "slater_matrix_dim",
    "nupdates",
    "slater_matrix",
    "slater_inverse",
    "col_update_index",
    "updates"
]

for cycle in cycles_list:
    cycle_name = "cycle_" + cycle

    new_cycle = fd.create_group(cycle_name)

    # Copy all datasets
    for group_name in groups:
        fs.copy(cycle_name + "/" + group_name, new_cycle)


print("Dataset successfully exported to %s" % source_dataset_path)

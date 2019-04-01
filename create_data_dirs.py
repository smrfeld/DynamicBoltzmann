import os
import sys

if len(sys.argv) != 2:
	print("Error: must specify single argument as name of directory (from current directory).")
	sys.exit(1)

dir_name = sys.argv[1]
print("Making directory structure in: " + str(dir_name))

# Make
if not os.path.exists(dir_name):
	os.makedirs(dir_name)

# Make subdirs
subdir_name = dir_name + "/diff_eq_rhs"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/moments"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)

subdir_name = dir_name + "/ixn_params"
if not os.path.exists(subdir_name):
	os.makedirs(subdir_name)
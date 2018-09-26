import os

# No samples
n_sample = 10

for i in range(0,n_sample):
	dir_name = "lattice_v%03d" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%03d/counts/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%03d/lattice/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
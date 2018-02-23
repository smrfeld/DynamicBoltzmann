import os

# No samples
n_sample = 100

for i in range(0,n_sample):
	dir_name = "lattice_v%02d" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%02d/counts/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%02d/lattice/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	dir_name = "lattice_v%02d/nns/" % i
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
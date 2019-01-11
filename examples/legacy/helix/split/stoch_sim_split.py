import numpy as np
import os

# Number splits
n_splits = 8;

# Size of each
split_size = 25;

# Remove existing dirs
for i in range(1,n_splits+1):
	os.system("rm -rf stoch_sim_"+str(i))
	os.system("mkdir stoch_sim_"+str(i))
	os.system("mkdir stoch_sim_"+str(i)+"/lattice_v001")
	os.system("mkdir stoch_sim_"+str(i)+"/lattice_v001/lattice")

# Move
for i in range(0,n_splits):
	for j in range(0,split_size):
		name_from = "stoch_sim/lattice_v001/lattice/%04d.txt"%(i*split_size+j)
		name_to = "stoch_sim_"+str(i+1)+"/lattice_v001/lattice/%04d.txt"%j
		print("copying " + name_from + " to: " + name_to)
		os.system("cp " + name_from + " " + name_to)

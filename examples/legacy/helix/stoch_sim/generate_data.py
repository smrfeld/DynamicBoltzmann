import numpy as np

# Length
box_length = 1000

# MinNo particles
minNo = 150

# amplitudes
amp = 100

# Period
p = 60.0

def fa(idx):
	return int(minNo + amp * np.cos(2 * np.pi * idx / p))

def fb(idx):
	return int(minNo + amp * np.sin(2 * np.pi * idx / p))

# No time
n = 200

# Trajs
a = []
b = []
c = []
for i in range (0,n):
	a.append(fa(i))
	b.append(fb(i))
	if i==0:
		c.append(0)
	else:
		c.append(c[-1]+1)

# Dir to write
dwrite = "lattice_v001/lattice/"


# Write
for i in range(0,n):
	# lattice reset
	latt = {}
	# fill
	for x in range(1,1+a[i]):
		latt[x] = "A"
	for x in range(300,300+b[i]):
		latt[x] = "B"
	for x in range(600,600+c[i]):
		latt[x] = "C"

	# write
	f = open(dwrite+"%04d.txt"%i,'w')
	for key,val in latt.items():
		f.write(str(key) + " " + str(val) + "\n")
	f.close()


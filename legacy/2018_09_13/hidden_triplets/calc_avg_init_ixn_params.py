import os

have = 0.0
jave = 0.0

for i in range(1,21):

	f = open("bimol_annihilation/lattice_v%02d/lattice/init.txt"%i)
	for line in f:
		if line == "":
			continue
		line_s = line.split(" ")
		if line_s[0] == 'h':
			have += float(line_s[1][:-2])
		else:
			jave += float(line_s[1][:-2])

	f.close()

print ("h ave: " + str(have / 20))
print ("J ave: " + str(jave / 20))

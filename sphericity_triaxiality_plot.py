inertiaPath = '/Users/brandyn/astro/Triaxiality/dat/'
plotPath = "/Users/brandyn/astro/Triaxiality/plots/"

inertiaFile = inertiaPath+"inertia_tensors_DM.dat"

fof = []
a_dm = []
b_dm = []
lambda_a_dm = []
lambda_b_dm = []
lambda_c_dm = []
with open(inertiaFile,'r') as input:
	for line in input.readlines():
		fields = line.split()
		fof.append(int(fields[0]))
		a_dm.append((float(fields[3])/float(fields[1]))**0.5)
		b_dm.append((float(fields[2])/float(fields[1]))**0.5)
		lambda_a_dm.append(float(fields[3])**0.5)
		lambda_b_dm.append(float(fields[2])**0.5)
		lambda_c_dm.append(float(fields[1])**0.5)

inertiaFile = inertiaPath+"inertia_tensors_AGN80.dat"

a_agn80 = []
b_agn80 = []
lambda_a_agn80 = []
lambda_b_agn80 = []
lambda_c_agn80 = []
with open(inertiaFile,'r') as input:
	for line in input.readlines():
		fields = line.split()
		fof.append(int(fields[0]))
		a_agn80.append((float(fields[3])/float(fields[1]))**0.5)
		b_agn80.append((float(fields[2])/float(fields[1]))**0.5)
		lambda_a_agn80.append(float(fields[3])**0.5)
		lambda_b_agn80.append(float(fields[2])**0.5)
		lambda_c_agn80.append(float(fields[1])**0.5)

inertiaFile = inertiaPath+"inertia_tensors_AGN87.dat"

a_agn87 = []
b_agn87 = []
lambda_a_agn87 = []
lambda_b_agn87 = []
lambda_c_agn87 = []
with open(inertiaFile,'r') as input:
	for line in input.readlines():
		fields = line.split()
		fof.append(int(fields[0]))
		a_agn87.append((float(fields[3])/float(fields[1]))**0.5)
		b_agn87.append((float(fields[2])/float(fields[1]))**0.5)
		lambda_a_agn87.append(float(fields[3])**0.5)
		lambda_b_agn87.append(float(fields[2])**0.5)
		lambda_c_agn87.append(float(fields[1])**0.5)
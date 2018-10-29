import matplotlib.pyplot as plt
import numpy as np

#read data in
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

#process data
a = np.array([np.array(a_dm), np.array(a_agn80), np.array(a_agn87)])
b = np.array([np.array(b_dm), np.array(b_agn80), np.array(b_agn87)])
lambda_a = np.array([np.array(lambda_a_dm), np.array(lambda_a_agn80), np.array(lambda_a_agn87)])
lambda_b = np.array([np.array(lambda_b_dm), np.array(lambda_b_agn80), np.array(lambda_b_agn87)])
lambda_c = np.array([np.array(lambda_c_dm), np.array(lambda_c_agn80), np.array(lambda_c_agn87)])

T = (1-b**2.)/(1-a**2.)

#plot data

plt.figure()
colors = ['green','blue','#ff00ff']
labels = ["DMO", "AGN 8.0", "AGN 8.7"]
plt.hist(np.transpose(a), bins=25,color=colors,histtype="step", fill=False, label=labels)
plt.legend()
plt.title("Minor:major axis ratio",fontsize=18)
plt.xlabel(r"$a$", fontsize=20)
plt.ylabel("Number of clusters",fontsize=18)
plt.tick_params(axis='both',labelsize=14)
plt.tight_layout()
plt.savefig(plotPath+"minor_axis.png")

plt.figure()
colors = ['green','blue','#ff00ff']
labels = ["DMO", "AGN 8.0", "AGN 8.7"]
plt.hist(np.transpose(b), bins=25,color=colors,histtype="step", fill=False, label=labels)
plt.legend()
plt.title("Intermediate:major axis ratio",fontsize=18)
plt.xlabel(r"$b$", fontsize=20)
plt.ylabel("Number of clusters",fontsize=18)
plt.tick_params(axis='both',labelsize=14)
plt.tight_layout()
plt.savefig(plotPath+"intermediate_axis.png")

plt.figure()
plt.hist(np.transpose(lambda_a), bins=25,color=colors,histtype="step", fill=False, label=labels)
plt.legend()
plt.title("Minor axis",fontsize=18)
#plt.xlabel(r"$b$", fontsize=20)
plt.ylabel("Number of clusters",fontsize=18)
plt.tick_params(axis='both',labelsize=14)
plt.tight_layout()
plt.savefig(plotPath+"unscaled_minor_axis.png")

plt.figure()
plt.hist(np.transpose(lambda_b), bins=25,color=colors,histtype="step", fill=False, label=labels)
plt.legend()
plt.title("Intermediate axis",fontsize=18)
#plt.xlabel(r"$b$", fontsize=20)
plt.ylabel("Number of clusters",fontsize=18)
plt.tick_params(axis='both',labelsize=14)
plt.tight_layout()
plt.savefig(plotPath+"unscaled_intermediate_axis.png")

plt.figure()
plt.hist(np.transpose(lambda_c), bins=25,color=colors,histtype="step", fill=False, label=labels)
plt.legend()
plt.title("Major axis",fontsize=18)
#plt.xlabel(r"$b$", fontsize=20)
plt.ylabel("Number of clusters",fontsize=18)
plt.tick_params(axis='both',labelsize=14)
plt.tight_layout()
plt.savefig(plotPath+"unscaled_major_axis.png")

plt.figure()
plt.hist(np.transpose(T), bins=25, color=colors, histtype="step", fill=False, label=labels)
plt.legend()
plt.xlabel(r"$T$", fontsize=20)
plt.ylabel("Number of clusters", fontsize=18)
plt.tick_params(axis='both', labelsize=14)
plt.tight_layout()
plt.savefig(plotPath+"triaxiality.png")

#plt.grid(True)
#plt.show()


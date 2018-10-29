import numpy as np 
import matplotlib.pyplot as plt 

catPath = '/Users/brandyn/astro/Triaxiality/dat/'
plotPath = "/Users/brandyn/astro/Triaxiality/plots/"

m200 = {}
bins = np.logspace(12, 15,6)
widths = (bins[1:] - bins[:-1])

catFile = catPath+"cluster_catalog_DM.dat"
with open(catFile,'r') as input:
	m200temp = []
	for line in input.readlines():
		fields = line.split()
		m200temp.append(10.**float(fields[1]))

	m200['DM'] = np.array(m200temp)

hist_dm = np.histogram(m200['DM'],bins=bins)
#hist_norm = hist[0]/widths

catFile = catPath+"cluster_catalog_AGN80.dat"
with open(catFile,'r') as input:
	m200temp = []
	for line in input.readlines():
		fields = line.split()
		m200temp.append(10.**float(fields[1]))

	m200['AGN80'] = np.array(m200temp)

hist_agn80 = np.histogram(m200['AGN80'],bins=bins)

catFile = catPath+"cluster_catalog_AGN87.dat"
with open(catFile,'r') as input:
	m200temp = []
	for line in input.readlines():
		fields = line.split()
		m200temp.append(10.**float(fields[1]))

	m200['AGN87'] = np.array(m200temp)

hist_agn87 = np.histogram(m200['AGN87'],bins=bins)

catFile = catPath+"cluster_catalog_REF.dat"
with open(catFile,'r') as input:
	m200temp = []
	for line in input.readlines():
		fields = line.split()
		m200temp.append(10.**float(fields[1]))

	m200['REF'] = np.array(m200temp)

hist_ref = np.histogram(m200['REF'],bins=bins)



cndens_dm = np.array(hist_dm[0])/(400.**3)
lcndens_dm = np.log10(cndens_dm)
hmf_dm = np.zeros(cndens_dm.shape, np.float)
hmf_dm[0:-1] = np.diff(cndens_dm)/np.diff(np.log10(bins[:-1]))
hmf_dm[-1] = (cndens_dm[-1] - cndens_dm[-2])/(np.log10(bins[-2]) - np.log10(bins[-3]))

cndens_agn80 = np.array(hist_agn80[0])/(400.**3)
lcndens_agn80 = np.log10(cndens_agn80)
hmf_agn80 = np.zeros(cndens_agn80.shape, np.float)
hmf_agn80[0:-1] = np.diff(cndens_agn80)/np.diff(np.log10(bins[:-1]))
hmf_agn80[-1] = (cndens_agn80[-1] - cndens_agn80[-2])/(np.log10(bins[-2]) - np.log10(bins[-3]))

cndens_agn87 = np.array(hist_agn87[0])/(400.**3)
lcndens_agn87 = np.log10(cndens_agn87)
hmf_agn87 = np.zeros(cndens_agn87.shape, np.float)
hmf_agn87[0:-1] = np.diff(cndens_agn87)/np.diff(np.log10(bins[:-1]))
hmf_agn87[-1] = (cndens_agn87[-1] - cndens_agn87[-2])/(np.log10(bins[-2]) - np.log10(bins[-3]))

cndens_ref = np.array(hist_ref[0])/(400.**3)
lcndens_ref = np.log10(cndens_ref)
hmf_ref = np.zeros(cndens_ref.shape, np.float)
hmf_ref[0:-1] = np.diff(cndens_ref)/np.diff(np.log10(bins[:-1]))
hmf_ref[-1] = (cndens_ref[-1] - cndens_ref[-2])/(np.log10(bins[-2]) - np.log10(bins[-3]))

ratio_agn80 = cndens_agn80/cndens_dm
ratio_agn87 = cndens_agn87/cndens_dm
ratio_ref =   cndens_ref/cndens_dm

#print hmf_dm
#print hist_norm
#print widths
#print hist
'''
plt.figure()
plt.scatter(bins[:-1], -hmf_dm, color='green')
plt.scatter(bins[:-1], -hmf_agn80, color='blue')
plt.scatter(bins[:-1], -hmf_agn87,color='#ff00ff')
plt.scatter(bins[:-1], -hmf_ref, color='red')
plt.xscale('log')
plt.yscale('log')
plt.ylim([10**-8, 10**-2])
plt.show()
'''
print hist_dm[0]
#plt.figure()
fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
ax1.plot(hist_dm[1][1:],lcndens_dm,color='green',label="DMO")
ax1.plot(hist_agn80[1][1:],lcndens_agn80,color='blue', label="AGN 8.0")
ax1.plot(hist_agn87[1][1:],lcndens_agn87,color='#ff00ff', label="AGN 8.7")
ax1.plot(hist_ref[1][1:],lcndens_ref,color='red', label="REF")
ax1.set_ylabel("number density")
ax1.legend()

ax2.plot(hist_agn80[1][1:], ratio_agn80, color='blue')
ax2.plot(hist_agn87[1][1:], ratio_agn87, color='#ff00ff')
ax2.plot(hist_ref[1][1:], ratio_ref, color='red')
ax2.set_xlabel(r"$M_{200}$ ($M_{\odot}$)")
ax2.set_ylabel("Ratio")
plt.xscale('log')
#plt.yscale('log')
plt.savefig(plotPath+"hmf.png")
#plt.show()
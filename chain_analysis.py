import numpy as np
import corner
import matplotlib.pyplot as plt
import cPickle
import colossus.cosmology.cosmology as Cosmology

fitsPath = "/Users/brandyn/astro/Triaxiality/MCMC/chains/"
fitsFile = fitsPath+"triax_spherical_sige05_fullgrid_sphericalfit_convergence.pkl"
cosmo = Cosmology.setCosmology('WMAP7-ML')

with open(fitsFile,'rb') as input:
	chains = cPickle.load(input)
	m200 = 200.*(4./3.)*np.pi*(chains["r200"])**3.*cosmo.rho_c(0.25)*1000**3*cosmo.h**2

	#samples = np.transpose([m200, chains["c"], chains["a"], chains["b"], np.sin(chains["theta"]), chains["phi"]])
	samples = np.transpose([m200, chains['c']])

#print chains["r200"]
print samples.shape
plt.figure()
corner.hist2d(chains["r200"], chains['c'])
plt.title("spherical spherical fit")
plt.xlabel("r200")
plt.ylabel("c")
plt.show()
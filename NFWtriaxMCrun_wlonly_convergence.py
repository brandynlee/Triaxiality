import numpy as np
import colossus.cosmology.cosmology as Cosmology
from radial_model_2d import RadialModel2D
import nfw
import NFWtriaxMCmodel_wlonly_spherical_convergence as NFWmodel
import cPickle
import dataHandle

cosmo = Cosmology.setCosmology('WMAP7-ML')
r200true = 2.
ctrue = 4.
atrue = 1.
btrue = 1.
thetatrue = 0.
phitrue = 0.
m200true = 200.*(4./3.)*atrue*btrue*np.pi*(r200true)**3.*cosmo.rho_c(0.25)*1000**3*cosmo.h**2

sige = 0.05
#sige=.05
#outPath = "/Users/brandyn/astro/CLASHtxtFiles/"
outPath = '/work/04144/blee/lonestar/Triaxiality/dat/'

x = np.linspace(-3.,3.,100)
y = np.linspace(-3.,3.,100)
xx,yy = np.meshgrid(x,y)

cluster = nfw.nfwProfile(parameters={'M':m200true, 'c':ctrue, 'a':atrue, 'b':btrue, 'theta':thetatrue, 'phi':phitrue}, z_lens=0.25, mass_definition='200c', cosmology='WMAP7-ML',z_source=1.)
kappa = cluster.convergence(xx,yy)

noise = np.random.normal(0.,sige,size=[100,100])/2**0.5

kappa = kappa + noise

zl = 0.25
zs=1.

nwalkers = 48
lenburn = 200
lenchain = 500
npix = 100

fitsPath = "/work/04144/blee/lonestar/Triaxiality/MCMC/chains/"

#r200chain, cchain, achain, bchain, thetachain, phichain = NFWtriaxMCmodel_wlonly.MCfit(nwalkers, lenburn, lenchain, xgal=xx, ygal=yy, eps1=g1, eps2=g2, zl=zl,zs=1.)
#r200chain, cchain = NFWtriaxMCmodel_wlonly_spherical.MCfit(nwalkers, lenburn, lenchain, xgal=xx, ygal=yy, eps1=g1, eps2=g2, zl=zl,zs=1.)
r200chain, cchain = NFWmodel.MCfit(nwalkers, lenburn, lenchain, xgal=xx, ygal=yy, k=kappa, zl=zl, zs=1.)
'''
print r200chain[-1]
print cchain[-1]

print achain[-1]
print bchain[-1]
print thetachain[-1]
print phichain[-1]

print "avg r200: "+str(np.average(r200chain))
print "avg c: "+str(np.average(cchain))
print "avg a: "+str(np.average(achain))
print "avg b: "+str(np.average(bchain))
print "avg theta:"+str(np.average(thetachain))
print "avg phi: "+str(np.average(phichain))
'''
#cluster_chains = {"r200": r200chain, "c": cchain, "a":achain, "b": bchain, "theta": thetachain, "phi":phichain}
cluster_chains = {"r200": r200chain, "c": cchain}
outFile = fitsPath+"triax_spherical_sige05_fullgrid_sphericalfit_convergence.pkl"
with open(outFile,'wb') as output:
	cPickle.dump(cluster_chains,output)
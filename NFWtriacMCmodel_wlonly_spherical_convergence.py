import numpy as np
import sys
#sys.path.append('/Users/brandyn/astro/baryonimpact')
#sys.path.append('/Users/brandyn/astro/')
#import dataHandle
#from Clusters import *
#import triaxNFW
#from astropy.io import fits
#from astropy.coordinates import SkyCoord
import emcee
import scipy.optimize as opt
#import corner
import matplotlib.pyplot as plt
import colossus.cosmology.cosmology as Cosmology
from radial_model_2d import RadialModel2D
import nfw

cosmo = Cosmology.setCosmology('WMAP7-ML')
#zl = 0.250 #for sim and ideal clusters
#zs = 1.
fov = 30.
#ndens = 30.
ndens = 30
#sige = 0.25
sige = .05
#sige=.2
r200min=.1
r200max = 5.
m200min = 200.*(4./3.)*np.pi*(r200min)**3.*cosmo.rho_c(0.25)*1000**3*cosmo.h**2
m200max = 200.*(4./3.)*np.pi*(r200max)**3.*cosmo.rho_c(0.25)*1000**3*cosmo.h**2
cmin=1.
#cmax = 30.
cmax = 30.
amin = .2
amax = 1.
bmin = .2
bmax = 1.
thetamin = 0.
thetamax = np.pi/2.
sinthetamin = np.sin(thetamin)
sinthetamax = np.sin(thetamax)
phimin = 0.
phimax = np.pi/2.
npix = 100

def nfwmodel(xgal,ygal,r200,c,zl, zs):
	#clusterTest = triaxNFW.triaxNFW(c = c, r200 = r200, a=a, b = b, theta=theta, phi=phi, z=zl, zs=zs, cosmoName='WMAP7-ML')
	m200 = 200.*(4./3.)*np.pi*(r200)**3.*cosmo.rho_c(zl)*1000**3*cosmo.h**2
	clusterTest = nfw.nfwProfile(parameters={'M':m200, 'c':c, 'a':1., 'b':1., 'theta':0., 'phi':0.}, z_lens=0.25, mass_definition='200c', cosmology='WMAP7-ML',z_source=1.)

	kappa = clusterTest.convergence(xgal,ygal)
	kappa = kappa.flatten()
	
	return kappa

def lnprior(p,zl):
	r200, c = p
	m200 = 200.*(4./3.)*np.pi*(r200)**3.*cosmo.rho_c(zl)*1000**3*cosmo.h**2
	#sintheta = np.sin(theta)
	if not (m200min < m200 < m200max and cmin < c < cmax):
		return -np.inf
	else:
		return 0.

	#qx,qy,q = returnq(a,b,theta,phi)
	#psi = returnpsi(a,b,theta,phi)
	#gauss1 = np.log(1.0/(np.sqrt(2*np.pi)*orient_err))-0.5*(psi-orient)**2/orient_err**2
	#gauss2 = np.log(1.0/(np.sqrt(2*np.pi)*axis_ratio_error))-0.5*(q-axis_ratio)**2/axis_ratio_err**2
	#return gauss1 + gauss2

def lnlike(p, xgal,ygal,k,zl, zs):
	print p
	r200, c = p
	kappa = nfwmodel(xgal,ygal,r200,c,zl, zs)
	loglike = (-.5)*np.sum(((kappa-convergence)/sige)**2.)
	
	print loglike
	#sig = (1-model**2.)*sige
	#print sig
	#sige_comp = sige/2**0.5
	#return -0.5*(np.sum(((redg-model)/(sige))**2.))
	#lnlikefunc = np.sum(((redg-model)/sig)**2.+2*np.log(sig))
	#lnlikefunc = np.nan_to_num(lnlikefunc)
	#return lnlikefunc
	return loglike

def lnprob(p, xgal,ygal,k, zl, zs):
	lp = lnprior(p,zl)
	if not np.isfinite(lp):
		return -np.inf
	return lp + lnlike(p,xgal,ygal,k,zl,zs)

def maxlike(xgal,ygal,k, r200_start, c_start, a_start, b_start, theta_start, phi_start, zl, zs):
	nll = lambda *args: -lnlike(*args)
	#return opt.minimize(nll, [r200_start, c_start], args = (rgal, redg),method='Powell')
	return opt.minimize(nll,[r200_start,c_start, a_start, b_start, theta_start, phi_start],args = (xgal,ygal,k,zl,zs),bounds=[(r200min,r200max),(cmin,cmax), (amin, amax), (bmin, bmax), (thetamin, thetamax), (phimin, phimax)],method='SLSQP')

#xgal and ygal need to be on a meshgrid
def MCfit(nwalkers, lenburn, lenchain,xgal, ygal, k, zl,zs, maxlikeReturn = False,):
	rin = 0.5
	rout = 3.

	#rgal = np.sqrt(xgal**2.+ygal**2.)
	#keepVals = np.logical_and(rgal<rout, rgal>rin)
	#eps1 = eps1[keepVals]
	#eps2 = eps2[keepVals]
	#rgal = rgal[keepVals]
	#xgal = xgal[keepVals]
	#ygal = ygal[keepVals]

	r200_start = 1.5
	c_start = 4.
	#a_start = .8
	#b_start = .81
	#phi_start = 0.
	#theta_start = 0.

	r200_start = np.random.uniform(r200min, r200max, size=nwalkers)
	c_start = np.random.uniform(cmin, 10., size=nwalkers)
	#a_start = np.random.uniform(amin, amax, size=nwalkers)
	#b_start = np.random.uniform(bmin, bmax, size=nwalkers)
	#theta_start = np.random.uniform(thetamin, thetamax, size=nwalkers)
	#phi_start = np.random.uniform(phimin, phimax, size=nwalkers)

	
	#result = maxlike(xgal,ygal,eps1,eps2,r200_start,c_start,a_start, b_start, theta_start, phi_start, zl, zs)
	ndim = 2
	#pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
	
	pos = np.transpose([r200_start, c_start])
	sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(xgal,ygal,k,zl,zs), threads=48)

	pos,prob,state = sampler.run_mcmc(pos,lenburn)
	sampler.reset()
	pos,prob,state = sampler.run_mcmc(pos,lenchain)

	r200chain = sampler.flatchain[:,0]
	cchain = sampler.flatchain[:,1]
	'''
	achain = sampler.flatchain[:,2]
	bchain = sampler.flatchain[:,3]
	thetachain = sampler.flatchain[:,4]
	phichain = sampler.flatchain[:,5]
	'''

	return (r200chain, cchain)

def returnfABC(a,b,theta,phi):
	sintheta = np.sin(theta)
	costheta = np.cos(theta)
	sintheta2 = sintheta**2
	costheta2 = costheta**2
	sinphi = np.sin(phi)
	cosphi = np.cos(phi)
	cosphi2 = cosphi**2
	sinphi2 = sinphi**2
	sin2phi = np.sin(2*phi)

	a2 = a**2
	b2 = b**2

	f = sintheta2*(cosphi2/a2+sinphi2/b2) + costheta2
	A = costheta2*(sinphi2/a2+cosphi2/b2) + sintheta2/(a2*b2)
	B = costheta*sin2phi*(1/a2 - 1/b2)
	C = sinphi2/b2 + cosphi2/a2

	return (f,A,B,C)

def returnq(a,b,theta,phi):
	f,A,B,C = returnfABC(a,b,theta,phi)
	qX = (2*f/(A+C - ((A-C)**2+B**2)**.5))**.5
	qY = (2*f/(A+C + ((A-C)**2+B**2)**.5))**.5

	q = qX/qY

	return (qX, qY, q)

def returnpsi(a,b,theta,phi):
	f,A,B,C = returnfABC(a,b,theta,phi)
	psi = 0.5*np.arctan2(B,A - C)
	return psi
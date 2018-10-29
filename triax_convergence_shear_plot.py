import matplotlib.pyplot as plt 
import numpy as np 
import colossus.cosmology.cosmology as Cosmology
from radial_model_2d import RadialModel2D
import nfw

def rebin2(inArray, shape):
	#rebin of 2D array into shape
	sh = shape[0],inArray.shape[0]//shape[0],shape[1],inArray.shape[1]//shape[1]
	return inArray.reshape(sh).mean(-1).mean(1)

cosmo = Cosmology.setCosmology('WMAP7-ML')
r200true = 1.8
ctrue = 4.
atrue = .6
btrue = .6
phitrue = 0.
thetatrue = np.pi/2.
m200true = 200.*(4./3.)*atrue*btrue*np.pi*(r200true)**3.*cosmo.rho_c(0.25)*1000**3*cosmo.h**2

nbin=50
shape = np.array([nbin,nbin])
x = np.linspace(-3.,3.,100)
y = np.linspace(-3.,3.,100)
xx,yy = np.meshgrid(x,y)

cluster = nfw.nfwProfile(parameters={'M':m200true, 'c':ctrue, 'a':atrue, 'b':btrue, 'theta':thetatrue, 'phi':phitrue}, z_lens=0.25, mass_definition='200c', cosmology='WMAP7-ML',z_source=1.)
kappa,gamma1,gamma2 = cluster.conv_shear_map(xx,yy)

gammaMag = (gamma1**2. + gamma2**2.)**0.5
gammaPhi = .5*np.arctan2(gamma2,gamma1)
gammaPhi = np.pi/2. - gammaPhi
gammax = gammaMag*np.cos(gammaPhi)
gammay = gammaMag*np.sin(gammaPhi)


gammaPhi = rebin2(gammaPhi, shape)
gammaMag = rebin2(gammaMag, shape)
gammax = rebin2(gammax, shape)
gammay = rebin2(gammay, shape)


gammaCap = .022
gammaXCapped = np.asarray(np.where(gammaMag<=gammaCap,gammax,np.zeros((nbin,nbin),dtype=float)))
gammaYCapped = np.asarray(np.where(gammaMag<=gammaCap,gammay,np.zeros((nbin,nbin),dtype=float)))

levels=[ .01, .02, .04, .08, .16, .32, .64]
plt.figure()
plt.contourf(yy,xx,kappa,levels, cmap=plt.cm.coolwarm)
plt.axes().set_aspect('equal')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xlabel('Mpc')
plt.ylabel('Mpc')


xx = rebin2(xx,shape)
yy = rebin2(yy,shape)

plt.figure()
plt.quiver(xx,yy,gammaXCapped, gammaYCapped, pivot='mid', headlength=1, headwidth=1, scale=.1)
plt.axes().set_aspect('equal')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xlabel('Mpc')
plt.ylabel('Mpc')
plt.show()
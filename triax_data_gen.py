import numpy as np
import colossus.cosmology.cosmology as Cosmology
from radial_model_2d import RadialModel2D
import nfw

cosmo = Cosmology.setCosmology('WMAP7-ML')
r200true = 2.
ctrue = 4.
atrue = .4
btrue = .4
thetatrue = 0.
phitrue = 0.
m200true = 200.*(4./3.)*atrue*btrue*np.pi*(r200true)**3.*cosmo.rho_c(0.25)*1000**3*cosmo.h**2
print m200true

sige = 0.05
#sige=.05
outPath = "/Users/brandyn/astro/CLASHtxtFiles/"
#outPath = '/work/04144/blee/lonestar/Triaxiality/dat/'

x = np.linspace(-3.,3.,100)
y = np.linspace(-3.,3.,100)
xx,yy = np.meshgrid(x,y)

cluster = nfw.nfwProfile(parameters={'M':m200true, 'c':ctrue, 'a':atrue, 'b':btrue, 'theta':thetatrue, 'phi':phitrue}, z_lens=0.25, mass_definition='200c', cosmology='WMAP7-ML',z_source=1.)
g1,g2 = cluster.reduced_shear(xx,yy)

noise1 = np.random.normal(0.,sige,size=[100,100])/2**0.5
noise2 = np.random.normal(0.,sige,size=[100,100])/2**0.5

eps1 = noise1+g1
eps2 = noise2+g2

np.save(outPath+"x.npy", xx.flatten())
np.save(outPath+"y.npy", yy.flatten())

np.savetxt(outPath+"triax_prolate_los_sige05_fullgrid.txt", np.transpose([eps1.flatten(), eps2.flatten()]), delimiter=',', fmt='%s')

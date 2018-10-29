import numpy as np
from radial_model_2d import RadialModel2D
import nfw
import matplotlib.pyplot as plt

'''
def charOverdensity(Delta=200, c):
        delta_c = (Delta/3.)*(c**3.)/(np.log(1. + c) - c/(1. + c))
        return delta_c #unitless

def nfwProfile(x,y,):
'''

nfw = nfw.nfwProfile(parameters={'M':1E15, 'c':4, 'a':0.4, 'b':1., 'theta':0., 'phi':np.pi}, z_lens=0.25, mass_definition='200c', cosmology='WMAP7-ML',z_source=1.)


zlos = np.linspace(-2.,2., 10)
xlos = np.linspace(-2.,2.,10)
ylos = np.linspace(-2.,2.,10)

xxlos, yylos = np.meshgrid(xlos,ylos)

#print nfw.density_3d(zlos,xlos,ylos)

#print nfw.surface_density(xxlos,yylos)

kappa,gamma1,gamma2 = nfw.conv_shear_map(xxlos,yylos)

#print kappa
#print gamma1
#print gamma2

#print gamma1/(1-kappa)
#print gamma2/(1-kappa)

#print nfw.reduced_shear(xxlos,yylos)

fig = plt.figure()

ax = fig.add_subplot(111)
ax.contour(xxlos,yylos,kappa)
ax.set_aspect('equal')
#ax.set_xlim([-5.,5.])
#ax.set_ylim([-5.,5.])
#ax.set_xticks(np.linspace(-4,4,8))
#ax.set_yticks(np.linspace(-4,4,8))
plt.grid()


gammaMag = (gamma1**2.+gamma2**2.)**0.5
gammaPhi = (0.5*np.arctan2(gamma2,gamma1))

#gammaMag = (gamma1_fft**2.+gamma2_fft**2.)**0.5
#gammaPhi = (0.5*np.arctan2(gamma2_fft,gamma1_fft))
gammaX = gammaMag*np.cos(gammaPhi)
gammaY = gammaMag*np.sin(gammaPhi)

plt.quiver(xxlos, yylos, gammaY, gammaX, pivot='mid', headlength=1, headwidth=1)
plt.show()
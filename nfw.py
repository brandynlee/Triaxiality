 #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 22:32:04 2018

@author: MattFong
"""

from radial_model_2d import RadialModel2D


import colossus.cosmology.cosmology as Cosmology

from colossus.halo.mass_so import M_to_R
from colossus.halo.mass_so import densityThreshold

import numpy as np


'''
############################################################################
                                  NFW
############################################################################
'''
#NFWProfile
class nfwProfile(RadialModel2D):
    def __init__(self, parameters, z_lens, mass_definition, cosmology=None, esp = None,z_source=None):
        #super().__init__(z_lens, func=self.density_3d, params=None, z_source=None)
        #Profile1D.__init__(self, zL, mdef, chooseCosmology)
        RadialModel2D.__init__(self, z_lens=z_lens, func=self.density_3d, params=parameters, z_source=z_source, grid=True,
         cosmology=cosmology,mass_definition=mass_definition)

        if isinstance(mass_definition, str) :
            self.mass_definition = mass_definition
        else :
            # Need to implement test of valid mass_definitions
            raise TypeError('mass_definition should be a string')

        if isinstance(cosmology, str) :
            self.cosmology_str = cosmology
        else :
            raise TypeError('cosmology should be a string')
            
        self.mass_definition = mass_definition
        self.chooseCosmology = cosmology
        listCosmologies = ['planck15-only', 'planck15', 'planck13-only', \
    'planck13', 'WMAP9-only', 'WMAP9-ML', 'WMAP9', 'WMAP7-only', 'WMAP7-ML', \
    'WMAP7', 'WMAP5-only', 'WMAP5-ML', 'WMAP5', 'WMAP3-ML', 'WMAP3', \
    'WMAP1-ML', 'WMAP1', 'bolshoi', 'millennium', 'powerlaw']        
        if cosmology is None:
            raise Exception('A name for the cosmology must be set.')
        if cosmology not in listCosmologies:
            msg = 'Cosmology must be one of ' + str(listCosmologies)    
            raise Exception(msg)
        if cosmology in listCosmologies:
            self.cosmo = Cosmology.setCosmology(cosmology)
        #input [M] = M_dot/h
        self.r_mass_definition = M_to_R(self.M_mass_definition*self.cosmo.h, self.zL, self.mass_definition)/1E3/self.cosmo.h/(self.a*self.b)**(1/3.) #Mpc from kpc/h
        self.Delta = int(mass_definition[:-1])
        
        
        #[rho_mass_definition] = M_dot/Mpc^3 from M_{\odot}h^2/kpc^3
        self.rho_mass_definition = (densityThreshold(self.zL, self.mass_definition) * ((1E3)**3.) *(self.cosmo.h)**2.)/self.Delta

        self.rs = self.r_mass_definition/self.c #Mpc
        
        if esp == None:
            self.esp = 1E-5
        else:
            self.esp = esp
        return
        
    '''
    ############################################################################
                               Analytic NFW profile
    ############################################################################
    '''
    def charOverdensity(self):
        sigma_c = (self.Delta/3.)*(self.c**3.)/(np.log(1. + self.c) - self.c/(1. + self.c))
        return sigma_c #unitless
    
    def density_3d(self, zlos, xlos, ylos):
        #R in Mpc
        #[sigma_c] = unitless
        #[rho_mdef] = M_dot / Mpc^3
        def rotation(
            zlos, xlos, ylos):
            sinphi = np.sin(self.phi)
            cosphi = np.cos(self.phi)
            sintheta = np.sin(self.theta)
            costheta = np.cos(self.theta)

            xpa = -sinphi*xlos - cosphi*costheta*ylos + cosphi*sintheta*zlos
            ypa = cosphi*xlos - sinphi*costheta*ylos + sinphi*sintheta*zlos
            zpa = sintheta*ylos + costheta*zlos

            return(xpa, ypa, zpa)

        xpa,ypa,zpa = rotation(zlos,xlos,ylos)
        R = ((xpa/self.a)**2.+(ypa/self.b)**2.+zpa**2.)**0.5
        const =  self.rho_mass_definition * self.charOverdensity() 
        rhoForm = 1./( (R/self.rs) * (1. + R/self.rs)**2.)
        return (const * rhoForm)
    ''' 
    def surface_density(self, r):
        #[r] = Mpc
        
        #[rs] = Mpc        
        rs = self.rs
        expSig = np.empty(len(r))
        for i in range(len(r)):
            if r[i]/self.rs < 1.0 - self.esp:
                expSig[i] = (1./((r[i]/rs)**2. - 1.))*(1. - 2.*np.arctanh(np.sqrt(\
                (1. - (r[i]/rs))/(1. + (r[i]/rs))))/np.sqrt(1. - (r[i]/rs)**2.))
            #if r[i]/rs == 1.0:
            #    expSig[i] = 1./3.
            if r[i]/rs >= 1.0 - self.esp and r[i]/rs <= 1.0 + self.esp:
                expSig[i] = 1./3.
            if r[i]/self.rs > 1.0 + self.esp:
                expSig[i] = (1./((r[i]/rs)**2. - 1.))*(1. - 2.*np.arctan(np.sqrt(\
                ((r[i]/rs) - 1.)/(1. + (r[i]/rs))))/np.sqrt((r[i]/rs)**2. - 1.))
        const = 2.*(self.rs)*self.charOverdensity()*(self.rho_mdef)
        #[Sigma] = M_dot / Mpc^2 
        return (expSig * const)
    def mean_surface_density(self, r):
        #[r] = Mpc
        x = r/self.rs
        const = 4.*self.rs*self.charOverdensity()*self.rho_mdef
        if type(x) is np.ndarray:
            #print x
            fnfw = np.piecewise(x,[x<1., x==1., x>1.], \
                        [lambda x:(2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x), \
                         lambda x:(1.+np.log(0.5)), \
                         lambda x:(2./np.sqrt(x*x-1.)*np.arctan2(np.sqrt(x-1.),np.sqrt(1.+x))+np.log(x/2.))/(x*x)])
            return const*fnfw
    
        else:
            if x<1:
                return const*(2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x)
            elif x==1:
                return const*(1.+np.log(0.5))
            else:
                return const*(2./np.sqrt(x*x-1.)*np.arctan2(np.sqrt(x-1.),np.sqrt(1.+x))+np.log(x/2.))/(x*x)
        
        
    
    def delta_sigma(self, r):
        #[r] = Mpc
        rs = self.rs
        
        x = r/rs
        expG = np.empty(len(x))
        for i in range(len(x)):
            if x[i] < 1.0 - self.esp:
                expG[i] = 8.0*np.arctanh(np.sqrt((1.0-(x[i]))/(1.0+(x[i]))))/(((x[i])**2.0)*np.sqrt(1.0-(x[i])**2.0))\
                + (4.0/((x[i])**2.0))*np.log((x[i])/2.0) \
                - 2.0/((x[i])**2.0 - 1.0)\
                + 4.0*np.arctanh(np.sqrt((1.0-(x[i]))/(1.0+(x[i]))))/(((x[i])**2.0 - 1.0)*np.sqrt(1.0-(x[i])**2.0))
            if x[i] > 1.0 - self.esp and x[i] < 1.0 + self.esp:
                expG[i] = (10.0 / 3.0) + (4.0 * np.log(1.0/2.0))
            if x[i] > 1.0 + self.esp:
                expG[i] = 8.0*np.arctan(np.sqrt(((x[i])-1.0)/(1.0+(x[i]))))/(((x[i])**2.0)*np.sqrt((x[i])**2.0 - 1.0))\
                + (4.0/((x[i])**2.0))*np.log((x[i])/2.0) \
                - 2.0/((x[i])**2.0 - 1.0)\
                + 4.0*np.arctan(np.sqrt(((x[i])-1.0)/(1.0+(x[i]))))/(((x[i])**2.0 - 1.0)**(3.0/2.0))
        
        
        #[rho_mdef] = M_dot/Mpc^3 
        #[charOverdensity] = unitless
        charOverdensity = self.charOverdensity()
        Const = (rs) * (self.rho_mdef) * (charOverdensity)
        # [Const] = 1 
        Const = Const
        return Const*expG
    
    '''

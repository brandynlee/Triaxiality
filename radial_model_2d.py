'''General profile class that inherits from Models'''
'''
import Triaxiality
from Triaxiality.models import Model
import scipy
from scipy import integrate
from Triaxiality.core.weak_lensing_functions import sigma_crit
RELATIVE_ERROR = 1E-6
import numpy as np
from Triaxiality.models.cosmo_dependent_models.cosmology_model import angular_diameter_distance
from Triaxiality.models.cosmo_dependent_models.cosmology_model import angular_diameter_distance_two_objects
'''
from model import Model
import scipy
from scipy import integrate
RELATIVE_ERROR=1E-6
from cosmology_model import angular_diameter_distance
from cosmology_model import angular_diameter_distance_two_objects
import numpy as np
import colossus.cosmology.cosmology as Cosmology
import math

#some constants
G = 4.5177E-48 #units of Mpc^3/solar mass/s^2
clight = 9.716E-15 #units of Mpc/s

class RadialModel2D(Model) :
    """
    Generalized superclass for 1D model profiles. It inherits from Models.

    Attributes
    ----------
    z_lens: float
        Redshift of the lensing object

    mass_definition: string
        Definition of mass, e.g. 'm200c'

    cosmology: string
        Label of cosmology, e.g. 'WMAP7'

    func: callable
        Functional form of the model, should be wrapped by the class

    params: list
        List of Parameter objects (default to None)

    z_source : float or iterable of floats (default to None)
        List of source redshifts or single source redshift
    
    """
    def __init__(self, z_lens, func, params=None, z_source=None, grid=False, cosmology=None,mass_definition=None) :
        
        if isinstance(z_lens, float) :
            self.z_lens = z_lens
        else :
            raise TypeError('z_lens should be a float')


        if isinstance(z_source, float) or \
           (np.iterable(z_source) and all(isinstance(z_s,float) for z_s in z_source) ) or\
           (z_source is None):
            self.z_source = z_source
        else :
            raise TypeError('z_source should be a float or list of floats')

        self.params = params
        try:
            assert type(self.params == dict)
        except:
            print(self.params)
        self.M_mass_definition = params['M'] #[M] = M_dot
        self.c = params['c']
        self.a = params['a']
        self.b = params['b']
        self.theta = params['theta']
        self.phi = params['phi']
        self.zL = z_lens
        self.grid = grid

        if isinstance(cosmology, str) :
            self.cosmology_str = cosmology
        else :
            raise TypeError('cosmology should be a string')
            
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
        
        
        # Need to implement r as the independent_var for all Profile1D models
        #super().__init__(func, params=params)
        #Model.__init__(self, )

        if z_source is not None :
            self.sigma_crit = self.calculate_sigma_crit_with_cosmology()

        f,A,B,C = self.returnfABC()
        self.f = f
        self.A = A
        self.B = B
        self.C = C

        q,qX,qY = self.returnq()
        self.q = q
        self.qX = qX
        self.qY = qY

        self.psi = self.returnpsi()
                
    def density_3d(self,z,x,y):
        """
        3D volume density profile, :math:'\rho'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -------
        rho: ndarray
            :math:'\rho', the 3D density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^3'.
            It has the same dimensions as r. 

        """
        return self.func(x,y, self.params)
    
    def surface_density(self, x,y):
        """
        Projected surface density profile, :math:'\Sigma'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -------
        sigma: ndarray
            :math:'\Sigma', the surface density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'. 
            It has the same dimensions as r.
        
        """

        if self.grid:
            xlos = list(x.flatten())
            ylos = list(y.flatten())
            npix = int(len(xlos)**0.5)

        int_lim = 15.
        surfacemass = []
        for i in range(len(xlos)):
            surfacemasstemp, surfacemass_err = scipy.integrate.quad(self.density_3d, -int_lim, int_lim, args=(xlos[i],ylos[i]))
            surfacemass.append(surfacemasstemp)

        if self.grid:
            surfacemass = np.array(surfacemass).reshape((npix,npix))

        return surfacemass

    
    def convergence(self, x,y):
        """
        Convergence, or dimensionless surface mass density, :math:'\kappa=\Sigma/\Sigma_{crit}'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -------
        kappa: ndarray
            :math:'\kappa', the convergence, which is unitless. It has the same dimensions of r.
        
        """

        sigma = self.surface_density(x,y)
        kappa = sigma/self.sigma_crit
        return kappa
    
    def shear(self, x,y):
        """
        Tangential shear, :math:'\gamma_{t}=\Delta\Sigma / \Sigma_{crit}'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -------
        gamma: ndarray
            :math:'\gamma_t', the tangential shear, which is unitless. It has the same dimensions as r.
        
        """

        kappa = self.convergence(x,y)
        fov = np.maximum(x) - np.minimum(x)

        kappaF = np.fft.fft2(kappa)
        pix = fov/kappaF.shape[0]

        #get the frequencies using pixel size as spacing
        freqx = np.fft.fftfreq(kappaF.shape[0],d=pix)
        freqy = np.fft.fftfreq(kappaF.shape[1],d=pix)
        freqX, freqY= np.meshgrid(freqx, freqy)

        #initialize and then calculate shear in fourier space
        gamma1F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)
        gamma2F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)

        gamma1F = kappaF*(freqX**2-freqY**2)/(freqX**2+freqY**2)
        gamma2F = 2.*kappaF*(freqX*freqY)/(freqX**2+freqY**2)

        #replace bad elements
        gamma1F = np.nan_to_num(gamma1F)
        gamma2F = np.nan_to_num(gamma2F)

        #fft back to position space
        gamma1 = np.fft.ifft2(gamma1F).real
        gamma2 = np.fft.ifft2(gamma2F).real

        #gammaMag = (gamma1**2+gamma2**2)**(1./2)
        #gammaPhi = (1./2)*np.arctan2(gamma2,gamma1)
        #gammaPhi = np.pi/2. - gammaPhi
        #gamma1 = np.transpose(gammaMag*np.cos(2.*gammaPhi))
        #gamma2 = np.transpose(gammaMag*np.sin(2.*gammaPhi))
        #print gammaMag
        #print gamma2
        #print gammaPhi

        gamma1 = np.transpose(gamma1)
        gamma2 = np.transpose(gamma2)
        return gamma1, gamma2

    def conv_shear_map(self,x,y):
        kappa = self.convergence(x,y)
        fov = np.amax(x) - np.amin(x)

        kappaF = np.fft.fft2(kappa)
        pix = fov/kappaF.shape[0]

        #get the frequencies using pixel size as spacing
        freqx = np.fft.fftfreq(kappaF.shape[0],d=pix)
        freqy = np.fft.fftfreq(kappaF.shape[1],d=pix)
        freqX, freqY= np.meshgrid(freqx, freqy)

        #initialize and then calculate shear in fourier space
        gamma1F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)
        gamma2F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)

        gamma1F = kappaF*(freqX**2-freqY**2)/(freqX**2+freqY**2)
        gamma2F = 2.*kappaF*(freqX*freqY)/(freqX**2+freqY**2)

        #replace bad elements
        gamma1F = np.nan_to_num(gamma1F)
        gamma2F = np.nan_to_num(gamma2F)

        #fft back to position space
        gamma1 = np.fft.ifft2(gamma1F).real
        gamma2 = np.fft.ifft2(gamma2F).real

        #gammaMag = (gamma1**2+gamma2**2)**(1./2)
        #gammaPhi = (1./2)*np.arctan2(gamma2,gamma1)
        #gammaPhi = np.pi/2. - gammaPhi
        #gamma1 = np.transpose(gammaMag*np.cos(2.*gammaPhi))
        #gamma2 = np.transpose(gammaMag*np.sin(2.*gammaPhi))
        #print gammaMag
        #print gamma2
        #print gammaPhi

        gamma1 = np.transpose(gamma1)
        gamma2 = np.transpose(gamma2)

        kappa = np.array(kappa)
        gamma1 = np.array(gamma1)
        gamma2 = np.array(gamma2)

        #gammaMag = (gamma1**2.+gamma2**2.)**0.5
        #gammaPhi = (0.5*np.arctan2(gamma2,gamma1))
        #gammaPhi = np.pi/2. - gammaPhi

        #gamma1 = gammaMag*np.cos(2*gammaPhi)
        #gamma2 = gammaMag*np.sin(2*gammaPhi)

        return kappa,gamma1, gamma2
    
    def reduced_shear(self, x,y):
        """
        Tangential reduced shear, :math:'g_t'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -------
        redg: ndarray
            :math:'g_t', the tangential reduced shear, which is unitless. It has the same dimensions as r.
        
        """
        # Make sure that self.sigma_crit is an attribute that is
        # filled in (otherwise, we didn't specify z_sources in the
        # constructor
        kappa,gamma1,gamma2 = self.conv_shear_map(x,y)

        redg1 = gamma1/(1-kappa)
        redg2 = gamma2/(1-kappa)
        return redg1,redg2

    def wl_and_sl(self,x,y):
        kappa, gamma1, gamma2 = self.conv_shear_map(x,y)

        g1 = gamma1/(1-kappa)
        g2 = gamma2/(1-kappa)

        sigma = kappa*self.sigma_crit

        arcsec = 3.908 / 1000 #number of Mpc per arcsecond for just the one cluster
        area = .03**2 #area of a pixel in Mpc^2
        area_arc = area / arcsec**2 #area of pixel in arcsecond^2
        x_arc = x / arcsec #coordinates in arcseconds
        y_arc = y / arcsec 

        r_arc = (x_arc**2.+y_arc**2.)**0.5

        M10 = np.sum(sigma[r_arc<10.])*area_arc
        M20 = np.sum(sigma[r_arc<20.])*area_arc
        M30 = np.sum(sigma[r_arc<30.])*area_arc
        M40 = np.sum(sigma[r_arc<40.])*area_arc

        return g1,g2,M10,M20,M30,M40

    def calculate_sigma_crit_with_cosmology(self) :
        # calculate da lens
        # calculate da source
        # calculate da lens source
        def sigma_crit(d_a_lens, d_a_source, d_a_lens_source):
            sigmaC = clight**2*d_a_source/(4*math.pi*G*d_a_lens_source*d_a_lens)

            return sigmaC

        d_a_lens = angular_diameter_distance(self.z_lens, self.cosmology_str)
        d_a_source = angular_diameter_distance(self.z_source, self.cosmology_str)
        d_a_lens_source = angular_diameter_distance_two_objects(self.z_lens, self.z_source, self.cosmology_str)

        return sigma_crit(d_a_lens, d_a_source, d_a_lens_source)

    def returnfABC(self):
        sintheta = np.sin(self.theta)
        costheta = np.cos(self.theta)
        sintheta2 = sintheta**2
        costheta2 = costheta**2
        sinphi = np.sin(self.phi)
        cosphi = np.cos(self.phi)
        cosphi2 = cosphi**2
        sinphi2 = sinphi**2
        sin2phi = np.sin(2*self.phi)

        a2 = self.a**2
        b2 = self.b**2

        f = sintheta2*(cosphi2/a2+sinphi2/b2) + costheta2
        A = costheta2*(sinphi2/a2+cosphi2/b2) + sintheta2/(a2*b2)
        B = costheta*sin2phi*(1/a2 - 1/b2)
        C = sinphi2/b2 + cosphi2/a2

        return (f,A,B,C)

    def returnq(self):
        qX = (2*self.f/(self.A+self.C - ((self.A-self.C)**2+self.B**2)**.5))**.5
        qY = (2*self.f/(self.A+self.C + ((self.A-self.C)**2+self.B**2)**.5))**.5

        q = qX/qY

        return (qX, qY, q)

    def returnpsi(self):
        psi = 0.5*np.arctan2(self.B,self.A - self.C)
        return psi

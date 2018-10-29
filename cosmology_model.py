
import astropy.cosmology as apycosmo
from astropy import constants
import numpy as np

def get_cosmology( cosmology_str ) :
    """
    Currently needs to be implemented to parse for astropy cosmology options
    
    cosmology_str : str
        Label of cosmology


    """
    #print("getting cosmology %s, need to implement options from astropy"%cosmology_str)

    return apycosmo.FlatLambdaCDM(H0=70.4, Om0=0.272)



def angular_diameter_distance(redshift, cosmology_str):
    cosmology = get_cosmology(cosmology_str)
    return cosmology.angular_diameter_distance( redshift )

def angular_diameter_distance_two_objects(z_lens, z_source, cosmology_str) :
    """
    Angular diameter distance between two objects
    Hogg+00 Eqn (19)
    
    Parameters
    ----------
    
    z_lens : float
        Lens redshift or redshifts
    
    z_source : float
        Source redshift or redshifts
    
    cosmology : str
        Name of cosmology
        
    """
    
    
    cosmology = get_cosmology(cosmology_str)
    #Need cosmology to be defined else the code breaks when using cosmo.H(0).to...
    Omega_k = 0 #for flat LCDM
    hubble_radius= constants.c.to('Mpc/s')/cosmology.H0.to('/s')                  
    ang_diameter_lens = cosmology.angular_diameter_distance(z_lens)
    ang_diameter_source = cosmology.angular_diameter_distance(z_source)
    #return (ang_diameter_source*np.sqrt(1.+ Omega_k*ang_diameter_lens**2/hubble_radius**2) \
    #        - ((1. + z_lens)/(1.+z_source))*ang_diameter_lens*np.sqrt(1.+ Omega_k*ang_diameter_source**2/hubble_radius**2))/(1.+z_source)
    return ((1.+z_source)*ang_diameter_source \
            - (1. + z_lens)*ang_diameter_lens)/(1.+z_source)



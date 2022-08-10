from neutronpy.config import Config
import numpy as np

class probability_distributions:

    def oned_distribution(self,x):
        """
        Function which gives the initial distribution of neutrons
        
        Args:
            x (float) : This parameter gives the neutron distribution at a location in the admitted 
            
        Returns:
            distribution (float) : Initial distribution at location x
        """
        _lambda = 100.0 # Half the inverse variance squared
        L = Config.critical_length_1d
        return np.exp(- _lambda *(( x - L /2.0) /( L /2.0) )**2.0) 

    def twod_distribution_special(self,x,y):

        L = Config.critical_length_2d
        return (1- x/L)*(1-y/L)*x*y/((L/4)**2)

    def twod_distribution_gaussian(self,x,y):

        L = Config.critical_length_2d
        return np.exp( ( ( x - L /2.0)**2.0 + ( y - L /2.0)**2.0 ) /( L)**2.0)
    
    def threed_distribution_gaussian(self,x,y,z):

        L = Config.critical_length_3d
        return np.exp( ( ( x - L /2.0)**2.0 + ( y - L /2.0)**2.0 + ( z - L /2.0)**2.0) /( 1.5*L)**2.0)
    
    def threed_distribution_special(self,x,y,z):

        L = Config.critical_length_3d
        return (1 -x/L)*(1 - y/L)*(1 - z/L)*x*y*z/((L/2)**3)
    
    def cylindrical_distribution_special(self,r,z):

        r1 = Config.critical_radius_cylindrical
        L = Config.critical_length_cylindrical

        return (1 -( r/r1)**2)*np.sin(np.pi*z/L)

    def cylindrical_distribution_gaussian(self,r,z):

        r1 = Config.critical_radius_cylindrical
        L = Config.critical_length_cylindrical

        return ((1 -( r / r1 ) ^2)**3)*((1 -(1 -2*z/L)**2)**3)

    def spherical_distribution_special(self,r):

        r1 = Config.critical_radius_spherical
        return 1 -( r/r1)**2

    def spherical_distribution_gaussian(self,r):
        r1 = Config.critical_radius_spherical
        return np.exp(-50*(r/r1)**2)


distributions = probability_distributions()

def choose_distribution(dimension : str, symmetry :str, distribution_type : str):

    if dimension == '1d':
        return distributions.oned_distribution

    elif dimension == '2d':
        if distribution_type == 'gaussian':
            return distributions.twod_distribution_gaussian
        else:
            return distributions.twod_distribution_special
        
    if dimension == '3d':

        if symmetry == 'cartesian':
                
            if distribution_type == 'gaussian':
                return distributions.threed_distribution_gaussian
            else:
                return distributions.threed_distribution_special
        
        elif symmetry == 'cylindrical':

            if distribution_type == 'gaussian':
                return distributions.cylindrical_distribution_gaussian
            else:
                return distributions.cylindrical_distribution_special
            
        elif symmetry == 'spherical':
            if distribution_type == 'gaussian':
                return distributions.spherical_distribution_gaussian
            else:
                return distributions.spherical_distribution_special


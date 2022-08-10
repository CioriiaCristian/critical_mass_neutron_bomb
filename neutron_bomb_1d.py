from neutronpy import core
from neutronpy.config import Config
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import pandas as pd

L = Config.critical_length_1d
x_domain = np.linspace(0,L,Config.number_of_data_points)
x_bar_domain = x_domain*(np.pi/L) # The conversion is needed for the simpson algorithm
time_domain = np.linspace(0,Config.time_limit_1d,Config.number_of_data_points)


def integrand(p,x_bar):
    """
    Integrand to compute fourier coefficients at time t = 0

    Args:
        x(float) : This parameter gives the neutron distribution at a location in the admitted domain
        p (integer) : the label of the fourier coefficient
    Returns:
        integrand : integrand evaluated for p coefficient at location x 
    """
    
    initial_distribution = core.choose_distribution(dimension='1d',distribution_type='gaussian', symmetry='')
    
    return (2/np.pi)*initial_distribution((L/np.pi)*x_bar)*np.sin(p*x_bar)

def neutron_distribution(ap, x, t):
    """
    Function which gives the neutron distribution at location x and time t
    
    Args:
        x (float) : location in the admitted domain
        t (float) : time
    Returns:
        neutron_distribution (float) : neutron distribution at time t and location x
    """
    _sum = 0
    for i in list(range(1,Config.coefficient_size_1d)):
        _sum+= ap[i]*np.exp((Config.eta - Config.mu *( i * np.pi / L )**2) *t)*np.sin(i*np.pi*x/L)
    return _sum



def main():

    ap = np.zeros(Config.coefficient_size_1d) # the list containing the fourier coefficients
    for i in list(range(1,Config.coefficient_size_1d)):
        ap[i] = simpson(integrand(i,x_bar_domain),x_bar_domain)

    dataframe = pd.DataFrame(data=ap)
    dataframe.to_csv('output\\1d_neutron_bomb\\coefficients.csv')

    values = np.zeros((Config.number_of_data_points,Config.number_of_data_points))
    for idx1,x_loc in enumerate(x_domain):
        for idx2,time_loc in enumerate(time_domain):
            values[idx1,idx2] = neutron_distribution(ap,x_loc,time_loc)

    axis_1,axis_2 = np.meshgrid(x_domain,time_domain)
    fig = plt.figure()
    ax = plt.axes(projection ="3d")
    ax.plot_surface(axis_1, axis_2, values, cstride=1, cmap='hot')
    plt.savefig('output\\1d_neutron_bomb\\1d_figure.png')
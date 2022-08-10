import numpy as np
from neutronpy import core
from neutronpy.config import Config
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.special import jn_zeros,jv
import pandas as pd
from tqdm import tqdm


L = Config.critical_length_cylindrical 
r1 = Config.critical_radius_cylindrical
z_domain = np.linspace(0,L,Config.number_of_data_points)
z_bar_domain = z_domain*(np.pi/L) # The conversion is needed for the simpson algorithm
r_domain = np.linspace(0,r1,Config.number_of_data_points)
r_bar_domain = r_domain*(np.pi/r1) # The conversion is needed for the simpson algorithm
time_domain = np.linspace(0,Config.time_limit_3d_cylindrical,Config.number_of_data_points)


def integrand(q,r0):
    
    _alpha_q = jn_zeros(0,q)[-1]
    return (2/(np.pi**2))/((jv(1, _alpha_q))**2)*jv(0, _alpha_q * r0/np.pi)*r0*(1-(r0/np.pi)**2)

def neutron_distribution(aq, r, z, t):
    
    n = Config.coefficient_size_cylindrical
    _sum = 0
    for q in list(range(1,n)):
        _alpha_q = jn_zeros(0,q)[-1]
        _sum+= aq[q-1]*np.exp((Config.eta*(r1*L)**2 - Config.mu*( (_alpha_q*L)**2 + (np.pi*r1)**2))*t/((r1*L)**2))*np.sin(np.pi*z/L)*jv(0,_alpha_q*r/r1)
    
    return _sum

def main():

    n = Config.coefficient_size_cylindrical
    aq = []
    for q in range(1,n):
        aq.append(simpson(integrand(q,r_bar_domain),r_bar_domain))
    
    dataframe = pd.DataFrame(data=aq)
    dataframe.to_csv('output\\3d_neutron_bomb_cylindrical\\coefficients.csv')

    print('Computing the neutron density for the cylindrical symmetry scenario')
    values = np.zeros((r_domain.size,z_domain.size,time_domain.size))
    for t_idx,t in enumerate(time_domain[-2:-1]):
        for r_idx, r_loc in enumerate(tqdm(r_domain)):
            for z_idx,z_loc in enumerate(z_domain):
                values[r_idx,z_idx,t_idx] = neutron_distribution(aq,r_loc,z_loc,t)
    
    axis_1,axis_2 = np.meshgrid(r_domain,time_domain)
    fig = plt.figure()
    ax = plt.axes(projection ="3d")
    ax.plot_surface(axis_1, axis_2, values[:,:,0], cstride=1, cmap='hot')
    plt.savefig('output\\3d_neutron_bomb_cylindrical\\cylindrical_symmetry_figure.png')

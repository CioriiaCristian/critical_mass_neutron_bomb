import numpy as np
from neutronpy import core
from neutronpy.config import Config
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from mpl_toolkits.mplot3d.axes3d import Axes3D
import pandas as pd
from scipy.optimize import root_scalar


R = Config.critical_radius_neumann
r_domain = np.linspace(1e-6,R,Config.number_of_data_points)
time_domain = np.linspace(0,Config.time_limit_neumann,Config.number_of_data_points)

def _critical_radius(r):
    _factor = np.sqrt(Config.eta/Config.mu)
    return -1 + r*_factor*(1/np.tan(r*_factor)) + (3/2)*(r/Config._lambda_t)

def k_equation(k):
    _beta = (1/R)*(1-(3*R)/(2*Config._lambda_t))
    return (1/np.tan(k*R))- (1/k)*_beta


def neutron_distribution(r,t):

    K = root_scalar(k_equation, bracket = [0.5,np.pi/R]).root
    alpha = Config.mu*(K**2) - Config.eta
    # Impose the n(R,t=0) = 1
    A = R/np.sin(K*R)

    return A * np.exp(-alpha*t)*(np.sin(K*r)/r)

def main ():

    statements = []

    """
    The domain of the cotan function is (0,np.pi), thus we have to look for values of r between 0 and np.pi * np.sqrt(Config.mu/Config.eta)
    """
    
    R_crit = root_scalar(_critical_radius, bracket = [1e-6,np.pi*np.sqrt(Config.mu/Config.eta)]).root

    statements.append('The critical radius is ' + str(R_crit) + ' meters')

    """
    The domain of the cotan function is (0,np.pi), thus we have to look for values of r between 0 and np.pi/R
    """
    K = root_scalar(k_equation, bracket = [0.5,np.pi/R]).root

    statements.append('The k factor for operational radius = ' + str(R) + ' m is ' + str(K) + 'm^-1')

    alpha = Config.mu*(K**2) - Config.eta

    statements.append('The exponential factor is equal to ' + str(alpha) + ' s^-1')

    dataframe = pd.DataFrame(data=statements)
    dataframe.to_csv('output\\3d_neutron_bomb_neuman\\statements.csv')


    values = np.zeros((Config.number_of_data_points,Config.number_of_data_points))
    for idx1,i in enumerate(r_domain):
        for idx2,j in enumerate(time_domain):
            values[idx1,idx2] = neutron_distribution(i,j)

    axis_1,axis_2 = np.meshgrid(r_domain,time_domain)
    fig = plt.figure()
    ax = plt.axes(projection ="3d")
    ax.plot_surface(axis_1, axis_2, values, cstride=1, cmap='binary')
    plt.savefig('output\\3d_neutron_bomb_neuman\\neumann_figure.png')

main()
import numpy as np
from neutronpy import core
from neutronpy.config import Config
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from mpl_toolkits.mplot3d.axes3d import Axes3D
import pandas as pd


r1 = Config.critical_radius_spherical
r_domain = np.linspace(0, r1, Config.number_of_data_points)
r_bar_domain = r_domain * (
    np.pi / r1
)  # The conversion is needed for the simpson algorithm
time_domain = np.linspace(0, Config.time_limit_spherical, Config.number_of_data_points)


def integrand(p, r_bar):
    """
    Integrand to compute fourier coefficients at time t = 0

    Args:
        r_bar (float) : This parameter gives the neutron distribution at a location in the admitted domain <- (0,np.pi)
        p (integer) : the label of the fourier coefficient
    Returns:
        integrand : integrand evaluated for p coefficient at location r_bar
    """

    initial_distribution = core.choose_distribution(
        dimension="3d", symmetry="spherical", distribution_type="special"
    )

    return (
        (2 * r1 / (np.pi**2))
        * r_bar
        * (initial_distribution(r_bar * (r1 / np.pi)))
        * np.sin(p * r_bar)
    )


def neutron_distribution(ap, r, t):
    """
    Function which gives the neutron distribution

    Args:
        r (float) : radial distance
        t (float) : time
    Returns:
        neutron_distribution (float) : neutron distribution at time t and radial distance r
    """
    n = Config.coefficient_size_spherical
    _sum = 0
    if r < 1e-8:
        # 0 radius edge case
        for i in list(range(1, n)):
            _sum += (
                ap[i - 1]
                * np.exp((Config.eta - Config.mu * (i * np.pi / r1) ** 2) * t)
                * (i * np.pi / r1)
            )
        return _sum
    else:
        for i in list(range(1, n)):
            _sum += (
                ap[i - 1]
                * (1 / r)
                * np.exp((Config.eta - Config.mu * (i * np.pi / r1) ** 2) * t)
                * np.sin(i * np.pi * r / r1)
            )
        return _sum


def main():

    n = Config.coefficient_size_spherical
    ap = []  # the list containing the fourier coefficients
    for i in list(range(1, n)):
        ap.append(simpson(integrand(i, r_bar_domain), r_bar_domain))

    dataframe = pd.DataFrame(data=ap)
    dataframe.to_csv("output\\3d_neutron_bomb_spherical\\coefficients.csv")

    values = np.zeros((Config.number_of_data_points, Config.number_of_data_points))
    for idx1, i in enumerate(r_domain):
        for idx2, j in enumerate(time_domain):
            values[idx1, idx2] = neutron_distribution(ap, i, j)

    axis_1, axis_2 = np.meshgrid(r_domain, time_domain)
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(axis_1, axis_2, values, cstride=1, cmap="hot")
    plt.savefig("output\\3d_neutron_bomb_spherical\\spherical_symmetry_figure.png")

from neutronpy import core
from neutronpy.config import Config
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import pandas as pd
from tqdm import tqdm

L = Config.critical_length_2d
x_domain = np.linspace(0, L, Config.number_of_data_points)
x_bar_domain = x_domain * (
    np.pi / L
)  # The conversion is needed for the simpson algorithm
y_domain = np.linspace(0, L, Config.number_of_data_points)
y_bar_domain = y_domain * (
    np.pi / L
)  # The conversion is needed for the simpson algorithm
time_domain = np.linspace(0, Config.time_limit_2d, Config.number_of_data_points)


def integrand(p, q, x_bar, y_bar):

    initial_distribution = core.choose_distribution(
        dimension="2d", symmetry="cartesian", distribution_type="special"
    )

    return (
        (4 / (np.pi**2))
        * initial_distribution((L / np.pi) * x_bar, (L / np.pi) * y_bar)
        * np.sin(p * x_bar)
        * np.sin(q * y_bar)
    )


def neutron_distribution(apq, x, y, t):
    _sum = 0
    for i in range(0, Config.coefficient_size_2d):
        for j in range(Config.coefficient_size_2d):
            _sum += (
                apq[i, j]
                * np.exp(
                    Config.eta * t
                    - Config.mu * ((i * np.pi / L) ** 2 + (j * np.pi / L) ** 2) * t
                )
                * np.sin(i * np.pi * x / L)
                * np.sin(j * np.pi * y / L)
            )
    return _sum


def main():

    n = Config.coefficient_size_2d
    apq = np.zeros((n, n))  # the list containing the fourier coefficients
    for i in list(range(n)):
        for j in list(range(n)):
            # Integration over y axis
            _integrated_y_strips = []
            for x_bar in x_bar_domain:
                _integrated_y_strips.append(
                    simpson(integrand(i, j, x_bar, y_bar_domain), y_bar_domain)
                )
            # Integration over x axis
            apq[i, j] = simpson(np.array(_integrated_y_strips), x_bar_domain)

    dataframe = pd.DataFrame(data=apq)
    dataframe.to_csv("output\\2d_neutron_bomb\\coefficients.csv")

    values = np.zeros((x_domain.size, y_domain.size, time_domain.size))
    print("Performing calculations for the 2d neutron bomb in cartesian coordinates")
    for t_idx, t in enumerate(time_domain[-2:-1]):
        for x_idx, x_loc in enumerate(tqdm(x_domain)):
            for y_idx, y_loc in enumerate(y_domain):
                values[x_idx, y_idx, t_idx] = neutron_distribution(apq, x_loc, y_loc, t)

    axis_1, axis_2 = np.meshgrid(x_domain, y_domain)
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(axis_1, axis_2, values[:, :, 0], cstride=1, cmap="hot")
    ax.set_xlabel('X axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel(r'Neutron density  $( m^{-3})$')
    plt.savefig("output\\2d_neutron_bomb\\2d_figure.png")

from neutronpy import core
from neutronpy.config import Config
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import pandas as pd
from tqdm import tqdm

N = int(Config.number_of_data_points / 2)  # This is to save time for the computations
L = Config.critical_length_3d
x_domain = np.linspace(0, L, N)
x_bar_domain = x_domain * (
    np.pi / L
)  # The conversion is needed for the simpson algorithm
y_domain = np.linspace(0, L, N)
y_bar_domain = y_domain * (
    np.pi / L
)  # The conversion is needed for the simpson algorithm
z_domain = np.linspace(0, L, N)
z_bar_domain = y_domain * (
    np.pi / L
)  # The conversion is needed for the simpson algorithm
time_domain = np.linspace(0, Config.time_limit_3d_cartesian, N)


def integrand(p, q, r, x_bar, y_bar, z_bar):

    initial_distribution = core.choose_distribution(
        dimension="3d", symmetry="cartesian", distribution_type="special"
    )

    return (
        ((2 / np.pi) ** 3)
        * initial_distribution(
            (L / np.pi) * x_bar, (L / np.pi) * y_bar, (L / np.pi) * z_bar
        )
        * np.sin(p * x_bar)
        * np.sin(q * y_bar)
        * np.sin(r * z_bar)
    )


def neutron_distribution(apqr, x, y, z, t):

    n = Config.coefficient_size_3d
    _sum = 0
    for i in range(1, n):
        for j in range(1, n):
            for k in range(1, n):
                _sum += (
                    apqr[i, j, k]
                    * np.exp(
                        Config.eta * t
                        - Config.mu
                        * (
                            (i * np.pi / L) ** 2
                            + (j * np.pi / L) ** 2
                            + (k * np.pi / L) ** 2
                        )
                        * t
                    )
                    * np.sin(i * np.pi * x / L)
                    * np.sin(j * np.pi * y / L)
                    * np.sin(i * np.pi * z / L)
                )
    return _sum


def main():

    n = Config.coefficient_size_3d
    apqr = np.zeros((n, n, n))  # the list containing the fourier coefficients
    flat_apqr = []
    print(
        "Computing coefficients for the 3d cartesian neutron diffusion solution (4 progress bars)"
    )
    for i in list(range(1, n)):
        for j in tqdm(list(range(1, n))):
            for k in list(range(1, n)):
                _integrated_y_strips = []
                # Integration over z axis
                for x_bar in x_bar_domain:
                    # Fix the x
                    _integrated_z_slice = []
                    for y_bar in y_bar_domain:
                        # Fix the y
                        _integrated_z_slice.append(
                            simpson(
                                integrand(i, j, k, x_bar, y_bar, z_bar_domain),
                                z_bar_domain,
                            )
                        )

                    # Integration over y axis
                    _integrated_y_strips.append(
                        simpson(np.array(_integrated_z_slice), y_bar_domain)
                    )
                # Integration over x
                apqr[i, j, k] = simpson(np.array(_integrated_y_strips), x_bar_domain)
                flat_apqr.append([i, j, k, apqr[i, j, k]]) 
                # The flat apqr is necessary for comparison with the paper results

    dataframe = pd.DataFrame(data=flat_apqr)
    dataframe.to_csv("output\\3d_neutron_bomb_cartesian\\coefficients.csv")

    values = np.zeros((x_domain.size, y_domain.size, z_domain.size, time_domain.size))
    # We only plot the neutron distribution at the last step of the time domain 
    for t_idx, t in enumerate(time_domain[-2:-1]):
        print(
            "Computing the neutron density for the 3d cartesian neutrom bomb problem"
        )
        for x_idx, x_loc in enumerate(tqdm(x_domain)):
            for y_idx, y_loc in enumerate(y_domain):
                for z_idx, z_loc in enumerate(z_domain):
                    values[x_idx, y_idx, z_idx, t_idx] = neutron_distribution(
                        apqr, x_loc, y_loc, z_loc, t
                    )

    axis_1, axis_2 = np.meshgrid(x_domain, y_domain)
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    # The neutron distribution is plotted only for z = L/2
    ax.plot_surface(axis_1, axis_2, values[:, :, int(N / 2), 0], cstride=1, cmap="hot")
    ax.set_title('Neutron distribution at z = L/2')
    ax.set_xlabel('X axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel(r'Neutron density  $( m^{-3})$')
    plt.savefig("output\\3d_neutron_bomb_cartesian\\3d_figure.png")

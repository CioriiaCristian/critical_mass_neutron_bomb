import os

# Path
path = "output/"

paths = [
    "output\\1d_neutron_bomb\\",
    "output\\2d_neutron_bomb\\",
    "output\\3d_neutron_bomb_cartesian\\",
    "output\\3d_neutron_bomb_cylindrical\\",
    "output\\3d_neutron_bomb_spherical\\",
    "output\\3d_neutron_bomb_neuman\\",
]

isdir = os.path.isdir(path)

if isdir == False:
    os.mkdir(path)
    for _path in paths:
        os.mkdir(_path)

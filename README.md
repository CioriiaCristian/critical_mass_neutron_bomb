# neutronpy

neutronpy is a python library developed so as to provide a way to compute neutron diffusion in a variety of circumstances pertaining to the creation of a runaway neutron generation (a.k.a atomic bomb).

When in the critical_mass_neutron_bomb project directory, open the terminal and run:

## Installation
First, set up the Python virtual environment and install all dependencies:

```shell
$   pip install -r requirements
$   python3 setup.py
```
Then run the program via

```shell
$   python3 run.py 
```

This will run all scenarios of a runaway reaction, as described in the article [Neutron Diffusion](https://www.researchgate.net/publication/323035158_Neutron_diffusion) by _Graham W. Griffiths_ . 

The scenarios include neutron diffusion:

1. Cartesian symmetry
- 1d
- 2d
- 3d
2. Cylindrical symmetry
3. Spherical symmetry (Dirichlet boundary conditions)
4. Spherical symmetry (Von Neumann boundary conditions)

The _run.py_ file will generate all the results of the paper for the different diffusion scenarios and save them in the **output** folder. 


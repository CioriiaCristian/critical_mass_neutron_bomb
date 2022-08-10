
class Config:

    mu =2.3446e+5 # diffusion constant
    eta=1.8958e+8 # neutron rate of formation
    _lambda_t = 0.0360 
    
    # All lenghts are in meters
    critical_length_1d = 0.111 
    critical_length_2d = 15.7/100
    critical_length_3d = 19.2/100 
    critical_length_cylindrical = 19.2/100 
    critical_radius_cylindrical = 0.104
    critical_radius_spherical = 0.1204 
    critical_radius_neumann = 8.5/100
    
    coefficient_size_1d = 30
    coefficient_size_2d = 6
    coefficient_size_3d = 5
    coefficient_size_cylindrical = 10
    coefficient_size_spherical = 30

    number_of_data_points = 100

    # All time limits are in seconds
    time_limit_1d = 1e-5 
    time_limit_2d = 1e-7
    time_limit_3d_cartesian = 2e-7
    time_limit_3d_cylindrical = 1e-5
    time_limit_spherical = 2e-6
    time_limit_neumann = 3e-6
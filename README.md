# Flow modeling in a flat channel with electrodes
![Modeling example](info/Studio_Project(1).gif?raw=true)

In this project, you can find scripts for calculating the Navier-Stokes equation and the Convective diffusion equation to understand speed, pressure, and concentration distribution. You can change parameters such as electrode width and height, distance between electrodes, etc.

## Repository structure
    .                   
    ├── info/                            # dir with usefull information about modeling process
    │   ├── Borders.png                  # name of borders
    │   └── channel_size.png             # channel params info                 
    ├── script/                          # dir with script
    │   ├── channel_shape.py             # create array for channel 
    │   ├── laplace.py                   # functions to calc Laplace equation 
    │   ├── navier_stokes.py             # functions to calc Navier-Stokes equation  
    │   ├── pressure_poisson_periodic.py # functions to calc Poisson equation 
    │   └── plot_gif.py                  # create gif from vector of calculated data 
    ├── calc_diffusion.py                # main file for init Diffusion equation calculations
    ├── params_channel_shape.py          # file to set channel params
    ├── params_conc.py                   # file to set concentration params
    ├── params_to_calc.py                # file to set initial params for Navier-Stoke
    └── README.md


## Equations
All the equations that were used in this project are here - > ![Equations](equations.md)

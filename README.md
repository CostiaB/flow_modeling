# Flow modeling in a flat channel with electrodes
![Modeling example](info/Studio_Project(1).gif?raw=true)

In this project, you can find scripts for calculating the Navier-Stokes equation and the Convective diffusion equation to understand speed, pressure, and concentration distribution. You can change parameters such as electrode width and height, distance between electrodes, etc. To run this scripts you need torch and CUDA.

## How to use
### To calc Navier-Stokes equation
1. If need change channel and electrodes params in params_channel_shape.py
2. If need change initial values for Navier-Stokes equation in params_to_calc.py
3. Start calc_navier.py
Arguments for calc_navier.py:

| **Arg** | **Description** | **Type** | **Default value**|
| ------ | ------ | ------ | ------ |
|--nt| Number of iterations when calc Navier-Stokes| int | 20 |
|--p_it| Number of iterations when calc Poisson | int | 1500 |
|--freq| Frequency | float | 0.1 |
|--freq_points| Number of points for each freq | int | 1000 |
|--sine| Flag to calc Navier-Stokes with periodical conditions | boolean | False |

Example:
```md
python3 calc_navier.py --nt 200000 --p_it 50 --freq 1 --sine --freq_points 1000
```
### Convective diffusion
1. If need change channel and electrodes params in params_channel_shape.py
2. If need change initial values for Convective diffusion in params_conc.py
3. Start calc_diffusion.py 
Arguments for calc_navier.py:

| **Arg** | **Description** | **Type** | **Default value**|
| ------ | ------ | ------ | ------ |
|--N_it| Number of diffusion steps| int | 5000 |
|--stop_it| Max number of "laplace2d_parralel_C0" iterations | int | 3000 |
|--l1norm_target| Target L1 norm value | float | 1e-7 |
|--sine| Flag to calc with Navier-Stokes data | boolean | False |
|--on_bot| Flag to allocate all electrodes on bottom line | boolean | False |

Example:
```md
calc_diffusion.py --stop_it 100000 --N_it 100 --l1norm_target 1e-07 --on_bot
```




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

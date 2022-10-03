'''
LiI density 133,85 g/mole
electrolyte density 1000 + 133.85*4 = 1535.4 g/l = 1535.4 kg/m^3 
dynamic viscosity 0.66 mkPa*s
kinematic viscosity      1.25e-6   mu/rho  =   0.00042985541    m^2/s
'''


dx = 1
dy = 1

dt = 1e-9

rho = 1.5354e3 / 1e18
nu = 1.25e-6 * 1e12
D = 1.5e-9

F = 0.000015
P0 = 0.0015 * 1e-9

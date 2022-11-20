import numpy as np
import matplotlib.pyplot as plt

u = np.load('/home/lapa/flow_modeling_1/results/navier_2022_11_08_11_10/u_gr.npy')
ch_it = np.load('/home/lapa/flow_modeling_1/results/navier_2022_11_08_11_10/ch_it.npy')
ts = np.load('/home/lapa/flow_modeling_1/results/navier_2022_11_08_11_10/ts.npy')

dt = 1e-9
nt = 1000
w = 1/(2*np.pi)
freq = 1 
w = freq * (2*np.pi)
t = nt / 200
nu = 1.25e-6 * 1e12

rho = 1.5354e3 / 1e18
eta = nu * rho 
dx = 0.05
F = 0.0015 


a = -(F/1500)/dx * 1/rho * np.exp(1j*w*t)
a = F/6.28e6

u_max = 1j*a/w * np.exp(-1j*w*t) * dt


tmp = [1j*a/w * np.exp(-1j*w*t*6.28) for t in ts]

a = F/1e6
tmp = [1j*a/w * np.exp(-1j*w*t) for t in ts]

plt.plot(tmp)
plt.ylim(min(tmp), max(tmp))
print([(1j*a/w * np.exp(-1j*w*t) ) for t in ch_it])
u_m = a * np.exp(-1j*w*t) *(450**2) / (12  * nu)

print(np.abs(u_max))
f_17_5 = (450**2/(12*eta) *  (F/1500)/dx) / 1e18

print(f_17_5)

print(f'u real {u[-2][650//2, 1400//2]}, u_max {np.abs(u_max)}, u_m {np.abs(u_m)}')

u_vec = [uu[100:550, 600].mean() for uu in u]

formula = r'$\displaystyle \bar{v} = -\frac{h^{2}}{12\eta} * \frac{dp}{dx}$'

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 22})
plt.figure(figsize=(15, 10), dpi=250)

plt.hlines(f_17_5, 0, len(u), color='r', linestyles='dashed')
plt.text(len(u)*0.95, f_17_5*0.8, formula , fontsize=20, 
         bbox=dict(facecolor='red', alpha=0.25))
plt.plot(u_vec)
plt.plot(tmp)
plt.ylim(min(u_vec), max(u_vec))
plt.xlabel('Iterations', fontsize=20)
plt.ylabel(r'Velocity ($\displaystyle \frac{m}{s}$)', fontsize=20)
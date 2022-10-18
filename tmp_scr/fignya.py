import numpy as np
from params_channel_shape import *
from params_to_calc import *
import matplotlib.pyplot as plt
from script.channel_shape import channel_shape
import matplotlib

step = 10
width = 100

p_gr = np.load('p_gr.npy')
u_gr = np.load('u_gr.npy')
v_gr = np.load('v_gr.npy')


u = u_gr[9]
v = v_gr[9]
p = p_gr[9]
x = np.linspace(0, nx, nx)
y = np.linspace(0, ny, ny)
X, Y = np.meshgrid(x, y)


XX = X[:int(X.shape[0]):step, 0:int(X.shape[1]):step]
YY = Y[:int(Y.shape[0]):step, 0:int(Y.shape[1]):step]
uu = u[:int(X.shape[0]):step, 0:int(X.shape[1]):step]
vv = v[:int(Y.shape[0]):step, 0:int(Y.shape[1]):step]
pp = p[:int(Y.shape[0]):step, 0:int(Y.shape[1]):step]

p_gif = channel_shape(field = field,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  anode_value=anode_value+10000, cathode_value=cathode_value+10000,
                  border_value=border_value+10000, nx=nx,
                  ny=ny)

fig = plt.figure(figsize = (50,20),dpi=100)
plt.contourf(XX, YY, pp, alpha=0.5,levels=50)
#plt.quiver(XX, YY, np.abs(pp))
plt.quiver(XX, YY, uu, vv)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","black"])
plt.imshow(p_gif, alpha = 0.1, cmap=cmap)

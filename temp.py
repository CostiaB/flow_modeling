import numpy as np
import gc
import os
import time
import shutil
import matplotlib.pyplot as plt
from script.navier_stokes import navier_stokes
from script.channel_shape import channel_shape
from script.plot_gif import plot_gif
from params_channel_shape import *
from params_to_calc import *
    
const = 10000

field = np.ones((ny, nx))
p = channel_shape(field = field,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  anode_value=anode_value+10000, cathode_value=cathode_value+10000,
                  border_value=border_value+10000, nx=nx,
                  ny=ny)



p[bottom_right[0]:, (nx-bottom_left[1])-10] = const
p[(ny-top_left[0]-50):, top_left[1]] = const
p[0:bottom_right[0]+50, bottom_right[1]] = const # C

p[0:(ny-top_left[0]), bottom_left[1]+50] = const #D

p[(ny-top_left[0])-30:(ny-top_left[0])-29, bottom_left[1]:top_left[1]] = const #!!! E

p[bottom_right[0]+29:bottom_right[0]+30, bottom_right[1]+1:(nx-bottom_left[1])] = const #!!! F



'''

     
        # dp/dy = 0 at bottom channel border
        
        p[bottom_right[0]-1:bottom_right[0], bottom_right[1]:(nx-bottom_left[1])] = \
                p[bottom_right[0]:bottom_right[0]+1, bottom_right[1]:(nx-bottom_left[1])]
        
        # dp/dx at left bottom border
        p[0:(ny-top_left[0]), bottom_left[1]-1] = p[0:(ny-top_left[0]), bottom_left[1]]'''

fig = plt.figure(figsize=(50, 20), dpi=500)
cmap = plt.get_cmap('RdBu')
plt.imshow(p[::-1,:])
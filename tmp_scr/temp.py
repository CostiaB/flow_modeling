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




        
        
#top E DONE!!!
#p[(ny-top_left[0])-1:(ny-top_left[0]), bottom_left[1]+1:top_left[1]+2] = const*2
#p[(ny-top_left[0])-2:(ny-top_left[0])-1, bottom_left[1]+1:top_left[1]+2] = -const*2 

p[(ny-top_left[0]):(ny-top_left[0])+1, bottom_left[1]+1:top_left[1]+2] = const*2
p[(ny-top_left[0])-1:(ny-top_left[0]), bottom_left[1]+1:top_left[1]+2] = -const*2 



#botton F DONE!!!
#p[bottom_right[0]:bottom_right[0]+1, bottom_right[1]-1:(nx-bottom_left[1])-1] = const*2 
#p[bottom_right[0]+1:bottom_right[0]+2, bottom_right[1]-1:(nx-bottom_left[1])-1] = -const*2 

p[bottom_right[0]-1:bottom_right[0], bottom_right[1]-1:(nx-bottom_left[1])-1] = const*2 
p[bottom_right[0]:bottom_right[0]+1, bottom_right[1]-1:(nx-bottom_left[1])-1] = -const*2 



#D DONE!!!
p[0:(ny-top_left[0]), bottom_left[1]] = const*2
p[0:(ny-top_left[0]), bottom_left[1]+1] = -const*2



#A DONE!!!
p[bottom_right[0]:, (nx-bottom_left[1])-1] = const*2
p[bottom_right[0]:, (nx-bottom_left[1])-2] = -const*2 



#B DONE!!!
#p[(ny-top_left[0]):, top_left[1]+1] = const*2
p[(ny-top_left[0]):, top_left[1]+2] = -const*2



# C DONE!!!
p[0:bottom_right[0], bottom_right[1]-1] = const*2
p[0:bottom_right[0], bottom_right[1]-2] = -const*2

# P
p[-1, top_left[1]+1:(nx-bottom_left[1])] = const*2
p[0, bottom_left[1]:bottom_right[1]] = const*2

t = 0
g = 1
if h==0:
    h = 1
elif h<0:
    t = 1
    g = 0
#bottom cathode        
bc = p[ny-top_left[0]-g*h-1:ny-top_left[0]-t*h+1, top_left[1]-s-w-d-w-1:top_left[1]-s-w-d+1]
#top cathode
tc = p[bottom_right[0]+t*h-1:bottom_right[0]+g*h+1, bottom_right[1]+s+w+d-1:bottom_right[1]+s+w+d+w+1] 
    
#p[bottom_right[0]+t*h-1:bottom_right[0]+g*h+1, bottom_right[1]+s+w+d-1:bottom_right[1]+s+w+d+w+1]  = -const*2


p[ny-top_left[0]-g*h-1:ny-top_left[0]-t*h, top_left[1]-s-w-d-w-1:top_left[1]-s-w-d+1] = -const*2


fig = plt.figure(figsize=(50, 20), dpi=200)
cmap = plt.get_cmap('RdBu')
plt.imshow(p[::-1,:])
plt.savefig('/home/lapa/hhjhj.png')



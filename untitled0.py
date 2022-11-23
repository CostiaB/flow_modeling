import numpy as np
from script.channel_shape import channel_shape
from params_channel_shape import *
import matplotlib.pyplot as plt


p = channel_shape(field=field,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,top_right=top_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, nx=nx,
                  ny=ny)

center_x, center_y = p.shape[0]//2, p.shape[1]//2 
height = 10
length = 30
        


'''
p[(center_x-height//2):(center_x+height//2), (center_y-length//2)-1] = 1000
p[(center_x-height//2):(center_x+height//2), (center_y-length//2)] = 10

p[(center_x-height//2):(center_x+height//2), (center_y+length//2)] = 1000
p[(center_x-height//2):(center_x+height//2), (center_y+length//2)-1] = 10

p[(center_x-height//2),(center_y-length//2):(center_y+length//2)] = 1000
p[(center_x-height//2)+1,(center_y-length//2):(center_y+length//2)] = 10

p[(center_x+height//2),(center_y-length//2):(center_y+length//2)] = 1000
p[(center_x+height//2)-1,(center_y-length//2):(center_y+length//2)] = 10'''


fig = plt.figure(figsize = (30,12), dpi=300)
cmap = plt.get_cmap('RdBu')
plt.imshow(p[::-1,:])
#plt.imshow(p[::-1,:], cmap=cmap, vmin=1, vmax=100)
#plt.savefig('/home/lapa/Figure 2022-10-31 112958.png')
        
tmp = [i for i,x in enumerate(p) if sum(x)>0]
y_start, y_stop = tmp[0], tmp[-1]
tmp = np.where(p[y_start]>0)[0]
x_start, x_stop = tmp[0], tmp[-1]


print(y_start, y_stop, x_start)
print(y_start, y_stop, x_stop)



'''
tmp = [[i,sum(x)] for i,x in enumerate(p_gif) if sum(x)==1880]
x_start, x_stop = tmp[0], tmp[-1]
tmp = np.where(p_gif[x_start]>0)[0]
y_start, y_stop = tmp[0], tmp[-1]
'''
#p_gif[y_start:y_stop, 10] = 1000
'''channel_start_x = np.where(p_gif==1)[0][0]
channel_start_y = np.where(p_gif[channel_start_x]==1)[0][0]
channel_end_x = np.where(p_gif==1)[-1][-1]
p_gif[0:100, channel_start_x] = 1'''




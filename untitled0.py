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

# dp/dy = 0 at bottom channel border E
   
p[(ny-top_left[0]-1), bottom_left[1]+2:top_left[1]] = -100
p[(ny-top_left[0]-2), bottom_left[1]+2:top_left[1]] = 100

# dp/dy = 0 at bottom channel border F

p[bottom_right[0], bottom_right[1]:(nx-bottom_left[1]-2)] = -100
p[bottom_right[0]+1, bottom_right[1]:(nx-bottom_left[1]-2)] = 100

# dp/dx at left bottom border D


p[0:(ny-top_left[0]-2), bottom_left[1]] = -100
p[0:(ny-top_left[0]-2), bottom_left[1]+1] = 100

# dp/dx at right top border A


p[(bottom_right[0]+2):, (nx-bottom_left[1]-1)] = -100
p[(bottom_right[0]+2):, (nx-bottom_left[1]-2)] = 100

# dp/dx at left top border B


p[(ny-top_left[0]):, top_left[1]] = -100
p[(ny-top_left[0]):, top_left[1]+1] = 100

# dp/dx at right bottom border C

p[0:bottom_right[0], bottom_right[1]-1] = -100
p[0:bottom_right[0], bottom_right[1]-2] = 100


# dp/dx fill emty space in borders
# D/E
# center
p[(ny-top_left[0]-1), bottom_left[1]] = -100
p[(ny-top_left[0]-2), bottom_left[1]+1] = 100
# side
p[(ny-top_left[0]-1), bottom_left[1]+1] = -100
p[(ny-top_left[0]-2), bottom_left[1]] = -100

#B/E
# center
p[(ny-top_left[0]-1), top_left[1]] = -100
p[(ny-top_left[0]-2), top_left[1]+1] = 100
# side
p[(ny-top_left[0]-1), top_left[1]+1] = 100
p[(ny-top_left[0]-2), top_left[1]] = 100

# C/F
# center
p[bottom_right[0], bottom_right[1]-1] = -100
p[bottom_right[0]+1, bottom_right[1]-2] = 100
# side
p[bottom_right[0], bottom_right[1]-2] = 100
p[bottom_right[0]+1, bottom_right[1]-1] = 100

#A/F
# center
p[bottom_right[0], (nx-bottom_left[1]-1)] = -100
p[bottom_right[0]+1, (nx-bottom_left[1]-2)] = 100
# side
p[bottom_right[0], (nx-bottom_left[1]-2)] = -100
p[bottom_right[0]+1, (nx-bottom_left[1]-1)] = -100


#electrode
g = 1
t = 0

#Top anode
#right
p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w)] = -100
p[(ny-top_left[0]-g*h)-2: (ny-top_left[0]-t*h)-2, (top_left[1]-s-w)-1] = 100
#left
p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s)-1] = -100
p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s)] = 100
#top
p[(ny-top_left[0]-g*h)-1, (top_left[1]-s-w) : (top_left[1]-s)] = -100
p[(ny-top_left[0]-g*h)-2, (top_left[1]-s-w) : (top_left[1]-s)] = 100

#field[(ny-top_left[0]-g*h) : (ny-top_left[0]-t*h), (top_left[1]-s-w) : (top_left[1]-s)] = anode_value

#Top cathode
#right
p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w-d-w)] = -100
p[(ny-top_left[0]-g*h)-2: (ny-top_left[0]-t*h)-2, ((top_left[1]-s-w-d-w))-1] = 100
#left
p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w-d)-1] = -100
p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w-d)] = 100
#top
p[(ny-top_left[0]-g*h)-1, (top_left[1]-s-w-d-w) : (top_left[1]-s-w-d)] = -100
p[(ny-top_left[0]-g*h)-2, (top_left[1]-s-w-d-w) : (top_left[1]-s-w-d)] = 100
#field[ny-top_left[0]-g*h:ny-top_left[0]-t*h, top_left[1]-s-w-d-w:top_left[1]-s-w-d] = cathode_value


#bottom anode
#right
p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s)] = -100
p[(bottom_right[0]+t*h)+2: (bottom_right[0]+g*h)+2, (bottom_right[1]+s)-1] = 100
#left
p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w)-1] = -100
p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w)] = 100
#top
p[(bottom_right[0]+g*h), (bottom_right[1]+s) : (bottom_right[1]+s+w)] = -100
p[(bottom_right[0]+g*h)+1, (bottom_right[1]+s) : (bottom_right[1]+s+w)]= 100

#field[bottom_right[0]+t*h:bottom_right[0]+g*h, bottom_right[1]+s:bottom_right[1]+s+w] = anode_value

#bottom cathode
#right
p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d)] = -100
p[(bottom_right[0]+t*h)+2: (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d)-1] = 100
#left
p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d+w)-1] = -100
p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d+w)] = 100
#top
p[(bottom_right[0]+g*h), (bottom_right[1]+s+w+d) : (bottom_right[1]+s+w+d+w)] = -100
p[(bottom_right[0]+g*h)+1, (bottom_right[1]+s+w+d) : (bottom_right[1]+s+w+d+w)]= 100

#field[bottom_right[0]+t*h:bottom_right[0]+g*h, bottom_right[1]+s+w+d:bottom_right[1]+s+w+d+w] = cathode_value

fig = plt.figure(figsize = (30,12), dpi=300)
cmap = plt.get_cmap('RdBu')
#plt.imshow(p[::-1,:])
plt.imshow(p[::-1,:], cmap=cmap, vmin=-10, vmax=10)
plt.savefig('/home/lapa/Figure 2022-10-31 112958.png')
        
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




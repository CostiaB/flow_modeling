import numpy as np
import itertools
from script.channel_shape import channel_shape


def laplace2d_parralel_C0(c, 
                          freq,
                          w, h, d, s,
                          bottom_left, top_left, bottom_right,
                          anode_value, cathode_value,
                          border_value,
                          ny, nx,
                          b,
                          rho, nu, 
                          dx, dy,
                          P0, stepcount,
                          nit, l1norm_target,
                          stop_it=2500):
    
    l1norm = 1
    cn = np.empty_like(c)

    dx2 = dx**2
    dy2 = dy**2

    

    for i in itertools.count(0, 1):
        cn = c.copy()
        c[1:-1, 1:-1] = ((dy2 * (cn[1:-1, 2:] + cn[1:-1, 0:-2]) +
                         dx2 * (cn[2:, 1:-1] + cn[0:-2, 1:-1])) /
                        (2 * (dx2 + dy2)))
            
                 
        
        c = channel_shape(field=c,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value,
                  nx=nx,
                  ny=ny)
        
        # dp/dy = 0 at top channel border E
        
        #p[(ny-top_left[0])-1:(ny-top_left[0]), bottom_left[1]+1:top_left[1]+2] =\
            #p[(ny-top_left[0])-2:(ny-top_left[0])-1, bottom_left[1]+1:top_left[1]+2]   
            
        c[(ny-top_left[0]):(ny-top_left[0])+1, bottom_left[1]+1:top_left[1]+2] =\
            c[(ny-top_left[0])-1:(ny-top_left[0]), bottom_left[1]+1:top_left[1]+2]    

        # dp/dy = 0 at bottom channel border F
        
        #p[bottom_right[0]:bottom_right[0]+1, bottom_right[1]-1:(nx-bottom_left[1])-1] =\
        #    p[bottom_right[0]+1:bottom_right[0]+2, bottom_right[1]-1:(nx-bottom_left[1])-1]
        
        c[bottom_right[0]-1:bottom_right[0], bottom_right[1]-1:(nx-bottom_left[1])-1] =\
            c[bottom_right[0]:bottom_right[0]+1, bottom_right[1]-1:(nx-bottom_left[1])-1]
        
        
        # dp/dx at left bottom border D
        
        #p[0:(ny-top_left[0]), bottom_left[1]] = p[0:(ny-top_left[0]), bottom_left[1]+1]
        c[0:(ny-top_left[0]), bottom_left[1]-1] = c[0:(ny-top_left[0]), bottom_left[1]]
        
        # dp/dx at right top border A
        
        #p[bottom_right[0]:, (nx-bottom_left[1])-1] = p[bottom_right[0]:, (nx-bottom_left[1])-2]
        c[bottom_right[0]:, (nx-bottom_left[1])] = c[bottom_right[0]:, (nx-bottom_left[1])-1]
        
        # dp/dx at left top border B
        
        #p[(ny-top_left[0]):, top_left[1]+1] = p[(ny-top_left[0]):, top_left[1]+2]
        c[(ny-top_left[0]):, top_left[1]] = c[(ny-top_left[0]):, top_left[1]+1]
        
        # dp/dx at right bottom border C
        #p[0:bottom_right[0], bottom_right[1]-1] = p[0:bottom_right[0], bottom_right[1]-2]
        c[0:bottom_right[0], bottom_right[1]] = c[0:bottom_right[0], bottom_right[1]-1]
        
        #check whether we can implement this condicions via channel shape or 
        #directly thought shape of channel
        
        # border (gradC,n) == 0 
        #x,y = border(w, 0, d,l)
        #c[0, x] = c[1, x]
        
        #c[-1, :] = c[-2, :]  #?????
        
        
        #c[np.where(mask_cathodes==1)] = 0
        #c[np.where(mask_anodes==1)] = ca

           
       
        if i%10 == 0 and i > 1:
            l1norm = np.abs((np.sum(np.abs(c[:]) - np.abs(cn[:])) /
            np.sum(np.abs(cn[:]))))
            if  l1norm < l1norm_target or i > stop_it:
                #c[np.where(mask_cathodes==1)] = 0
                #c[np.where(mask_anodes==1)] = 0
                break
         

    return c




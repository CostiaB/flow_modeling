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




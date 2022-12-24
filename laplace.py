import torch
import itertools
from script.channel_shape import channel_shape

def laplace2d_parralel_C0(c,
                          w, h, d, s,
                          bottom_left, top_left, 
                          bottom_right, top_right,
                          anode_value, cathode_value,
                          border_value,
                          ny, nx,
                          b,
                          rho, nu, 
                          dx, dy,
                          P0, stepcount,
                          l1norm_target,
                          stop_it=2500,
                          silent=False, on_bot=False):
    
    l1norm = 1
    cn = torch.empty_like(c, device='cuda', requires_grad=False, dtype=torch.double)

    dx2 = dx**2
    dy2 = dy**2

    

    for i in itertools.count(0, 1):
        cn = c.clone()
        c[1:-1, 1:-1] = ((dy2 * (cn[1:-1, 2:] + cn[1:-1, 0:-2]) +
                         dx2 * (cn[2:, 1:-1] + cn[0:-2, 1:-1])) /
                        (2 * (dx2 + dy2)))
            
                 
        
        c = channel_shape(field=c,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right, top_right=top_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, nx=nx,
                  ny=ny, on_bot=on_bot)
        
        if on_bot:
            
            #For electrodes when all allocated at bottom line  
            # dp/dy = 0 at bottom channel border E
            c[(ny-top_left[0]-1), bottom_left[1]+2:top_left[1]] =\
                c[(ny-top_left[0]-2), bottom_left[1]+2:top_left[1]] 
                
            # dp/dy = 0 at bottom channel border F
            c[bottom_right[0], bottom_right[1]:(bottom_right[1]+s)] =\
                c[bottom_right[0]+1, bottom_right[1]:(bottom_right[1]+s)] 
    
            c[bottom_right[0], (bottom_right[1]+s+w):(bottom_right[1]+s+w+d)] =\
                c[bottom_right[0]+1, (bottom_right[1]+s+w):(bottom_right[1]+s+w+d)] 
    
            c[bottom_right[0], (bottom_right[1]+s+w+d+w):(top_left[1]-s-w-d-w)] =\
                c[bottom_right[0]+1, (bottom_right[1]+s+w+d+w):(top_left[1]-s-w-d-w)] 
                
            c[bottom_right[0], (top_left[1]-s-w-d):(top_left[1]-s-w) ] =\
                c[bottom_right[0]+1, (top_left[1]-s-w-d):(top_left[1]-s-w)] 
                
            c[bottom_right[0], (top_left[1]-s):(nx-bottom_left[1]-2)] =\
                c[bottom_right[0]+1, (top_left[1]-s):(nx-bottom_left[1]-2)] 
                
                
        else:            
                
            # dp/dy = 0 at bottom channel border E
    
            c[(ny-top_left[0]-1), bottom_left[1]+2:(top_left[1]-s-w-d-w)] =\
                c[(ny-top_left[0]-2), bottom_left[1]+2:(top_left[1]-s-w-d-w)] 
    
            c[(ny-top_left[0]-1), (top_left[1]-s-w-d):(top_left[1]-s-w)] =\
                c[(ny-top_left[0]-2), (top_left[1]-s-w-d):(top_left[1]-s-w)] 
    
            c[(ny-top_left[0]-1), (top_left[1]-s):top_left[1]] =\
                c[(ny-top_left[0]-2), (top_left[1]-s):top_left[1]] 
    
            
            # dp/dy = 0 at bottom channel border F
            
    
            c[bottom_right[0], bottom_right[1]:(bottom_right[1]+s)] =\
                c[bottom_right[0]+1, bottom_right[1]:(bottom_right[1]+s)] 
    
            c[bottom_right[0], (bottom_right[1]+s+w):(bottom_right[1]+s+w+d)] =\
                c[bottom_right[0]+1, (bottom_right[1]+s+w):(bottom_right[1]+s+w+d)] 
    
            c[bottom_right[0], (bottom_right[1]+s+w+d+w):(nx-bottom_left[1]-2)] =\
                c[bottom_right[0]+1, (bottom_right[1]+s+w+d+w):(nx-bottom_left[1]-2)] 

        # dp/dx at left bottom border D


        c[0:(ny-top_left[0]-2), bottom_left[1]] =\
            c[0:(ny-top_left[0]-2), bottom_left[1]+1] 

        # dp/dx at right top border A


        c[(bottom_right[0]+2):, (nx-bottom_left[1]-1)] =\
            c[(bottom_right[0]+2):, (nx-bottom_left[1]-2)] 

        # dp/dx at left top border B


        c[(ny-top_left[0]):, top_left[1]] =\
            c[(ny-top_left[0]):, top_left[1]+1] 

        # dp/dx at right bottom border C

        c[0:bottom_right[0], bottom_right[1]-1] =\
            c[0:bottom_right[0], bottom_right[1]-2] 


        # dp/dx fill emty space in borders
        # D/E
        # center
        c[(ny-top_left[0]-1), bottom_left[1]] =\
            c[(ny-top_left[0]-2), bottom_left[1]+1] 
        # side
        c[(ny-top_left[0]-1), bottom_left[1]+1] =\
            c[(ny-top_left[0]-2), bottom_left[1]]

        #B/E
        # center
        c[(ny-top_left[0]-1), top_left[1]] =\
            c[(ny-top_left[0]-2), top_left[1]+1] 
        # side
        c[(ny-top_left[0]-1), top_left[1]+1] =\
            c[(ny-top_left[0]-2), top_left[1]] 

        # C/F
        # center
        c[bottom_right[0], bottom_right[1]-1] =\
            c[bottom_right[0]+1, bottom_right[1]-2] 
        # side
        c[bottom_right[0], bottom_right[1]-2] =\
            c[bottom_right[0]+1, bottom_right[1]-1] 

        #A/F
        # center
        c[bottom_right[0], (nx-bottom_left[1]-1)] =\
            c[bottom_right[0]+1, (nx-bottom_left[1]-2)] 
        # side
        c[bottom_right[0], (nx-bottom_left[1]-2)] =\
            c[bottom_right[0]+1, (nx-bottom_left[1]-1)]
            
        
        if not silent:    
            if i % 1000 == 0:
                print(i)
        
        if i % 10 == 0 and i > 1:
           
            l1norm = torch.abs((torch.sum(torch.abs(c[:]) - torch.abs(cn[:])) /
                                torch.sum(torch.abs(cn[:]))))
            if  l1norm < l1norm_target or i > stop_it:
                break
         

    return c




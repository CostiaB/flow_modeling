import numpy as np
import torch
from script.channel_shape import channel_shape

def build_up_b(rho, dt, dx, dy, u, v, w, h, d, s, 
               bottom_left, bottom_right, top_left, top_right,
               anode_value, cathode_value, border_value,
               nx, ny):
    
    b = torch.zeros_like(u, requires_grad=False, dtype=torch.double)
    b = channel_shape(field=b,
                      w=w, h=h, d=d, s=s,
                      bottom_left=bottom_left, top_left=top_left,
                      bottom_right=bottom_right, top_right=top_right,
                      anode_value=anode_value, cathode_value=cathode_value,
                      border_value=border_value, 
                      nx=nx, ny=ny)
    
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    
    # Periodic BC Pressure @ x = 2
    b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1,-2]) / (2 * dx) +
                                    (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
                          ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 -
                          2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
                               (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) -
                          ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2))
    
    # Periodic BC Pressure @ x = 0
    b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) +
                                   (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) -
                         ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 -
                         2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
                              (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx))-
                         ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2))
    
    #b[315:326, 315:326] = 0 # test
    return b

def pressure_poisson_periodic(p_z, 
                              freq,
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
                              nit,
                              dt):
    pn = torch.empty_like(p_z, requires_grad=False, dtype=torch.double)
    p = p_z.clone()

    dx2 = dx**2
    dy2 = dy**2
    

    
    for q in range(nit):
        pn = p.clone()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx2) /
                         (2 * (dx2 + dy2)) -
                         dx2 * dy2 / (2 * (dx2 + dy2)) * b[1:-1, 1:-1])
        
        # Periodic BC Pressure @ x = 2
        p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 +
                        (pn[2:, -1] + pn[0:-2, -1]) * dx**2) /
                       (2 * (dx**2 + dy**2)) -
                       dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, -1])

        # Periodic BC Pressure @ x = 0
        p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 +
                       (pn[2:, 0] + pn[0:-2, 0]) * dx**2) /
                      (2 * (dx**2 + dy**2)) -
                      dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 0])     
       
        
        
        p = channel_shape(field=p,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right, top_right=top_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, 
                  nx=nx, ny=ny)
       
        
        # One electrode in the middle borders

        
        
        
        
        
        # dp/dy = 0 at top channel border E
 
        '''
        For nonlinear channel with non frequency 
        '''    
        #p[(ny-top_left[0]):(ny-top_left[0])+1, bottom_left[1]+1:top_left[1]+2] =\
        #    p[(ny-top_left[0])-1:(ny-top_left[0]), bottom_left[1]+1:top_left[1]+2]  
        
        
        p[(ny-top_left[0]):(ny-top_left[0])+1, :] = p[(ny-top_left[0])-1:(ny-top_left[0]), :] 

        # dp/dy = 0 at bottom channel border F
        
        '''
        For nonlinear channel with non frequency 
        '''    
        #p[bottom_right[0]-1:bottom_right[0], bottom_right[1]-1:(nx-bottom_left[1])-1] =\
        #    p[bottom_right[0]:bottom_right[0]+1, bottom_right[1]-1:(nx-bottom_left[1])-1]
        
        p[bottom_right[0]-1:bottom_right[0], :] = p[bottom_right[0]:bottom_right[0]+1, :]
        
        # dp/dx at left bottom border D
        

        #p[0:(ny-top_left[0]), bottom_left[1]-1] = p[0:(ny-top_left[0]), bottom_left[1]]
        
        # dp/dx at right top border A
        
       
        #p[bottom_right[0]:, (nx-bottom_left[1])] = p[bottom_right[0]:, (nx-bottom_left[1])-1]
        
        # dp/dx at left top border B
        
        
        #p[(ny-top_left[0]):, top_left[1]] = p[(ny-top_left[0]):, top_left[1]+1]
        
        # dp/dx at right bottom border C

        #p[0:bottom_right[0], bottom_right[1]] = p[0:bottom_right[0], bottom_right[1]-1]
        
        #check whether we can implement this condicions via channel shape or 
        #directly thought shape of channel
        
        
        #100 points for ech frequency
        '''
        For nonlinear channel with non frequency 
        '''
        '''
        time_mult = 1/(100 * freq)
        
        freq_cos = np.cos(2 * np.pi * (stepcount * time_mult) * freq)  
        
        
        p[-1, top_left[1]+1:(nx-bottom_left[1])] = P0 * freq_cos #check id second part is correct
        p[0, bottom_left[1]:bottom_right[1]] = - P0 * freq_cos #check id second part is correct
        '''
        '''
        For expirement with flat channel and with the constant frequency
        '''
        
        freq_cos = 1
        #p[315:326, 315:326] = 0 #test
        #p[315, 315:326] = p[316, 315:326] #test
        #p[325, 315:326] = p[324, 315:326] #test
        
        #p[315:326, 325] = p[315:326, 324] #test
        #p[315:326, 315] = p[315:326, 316] #test    
        
        
        #p[98:551, 0] = P0 * freq_cos
        #p[98:551, -1] = P0 * freq_cos
        #p[250:399, 1139] = P0 * freq_cos
        
        #p[-1, top_left[1]:] = 0
        #p[0, bottom_left[1]:bottom_right[1]] = +P0      
              
    

    return p
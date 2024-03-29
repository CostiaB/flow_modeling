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
    '''
    # Periodic BC Pressure @  y = top value
    b[-1, 1:-1] = (rho * (1 / dt * ((u[0, 1:-1] - u[-2, 1:-1]) / (2 * dx) +
                                    (v[-1, 2:] - v[-1, 0:-2]) / (2 * dy)) -
                          ((u[0, 1:-1] - u[-2, 1:-1]) / (2 * dx))**2 -
                          2 * ((u[-1, 2:] - u[-1, 0:-2]) / (2 * dy) *
                               (v[0, 1:-1] - v[-2, 1:-1]) / (2 * dx)) -
                          ((v[-1, 2:] - v[-1, 0:-2]) / (2 * dy))**2))
    
    # Periodic BC Pressure @ y = 0
    b[0, 1:-1] = (rho * (1 / dt * ((u[1, 1:-1] - u[-1, 1:-1]) / (2 * dx) +
                                   (v[0, 2:] - v[0, 0:-2]) / (2 * dy)) -
                         ((u[1, 1:-1] - u[1, 1:-1]) / (2 * dx))**2 -
                         2 * ((u[0, 2:] - u[0, 0:-2]) / (2 * dy) *
                              (v[1, 1:-1] - v[-1, 1:-1]) / (2 * dx))-
                         ((v[0, 2:] - v[0, 0:-2]) / (2 * dy))**2))
    
    '''
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
        '''
        # Periodic BC Pressure @ y = top
        p[-1, 1:-1] = (((pn[0, 1:-1] + pn[-2, 1:-1])* dy**2 +
                        (pn[-1, 2:] + pn[-1, 0:-2]) * dx**2) /
                       (2 * (dx**2 + dy**2)) -
                       dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[-1, 1:-1])

        # Periodic BC Pressure @ y = 0
        p[0, 1:-1] = (((pn[1, 1:-1] + pn[-1, 1:-1])* dy**2 +
                       (pn[0, 2:] + pn[0, 0:-2]) * dx**2) /
                      (2 * (dx**2 + dy**2)) -
                      dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[0, 1:-1])     
       
        '''
        
        p = channel_shape(field=p,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right, top_right=top_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, 
                  nx=nx, ny=ny)
       
        
        
        # dp/dy = 0 at bottom channel border E
           
        p[(ny-top_left[0]-1), bottom_left[1]+2:top_left[1]] =\
            p[(ny-top_left[0]-2), bottom_left[1]+2:top_left[1]] 

        # dp/dy = 0 at bottom channel border F

        p[bottom_right[0], bottom_right[1]:(nx-bottom_left[1]-2)] =\
            p[bottom_right[0]+1, bottom_right[1]:(nx-bottom_left[1]-2)] 

        # dp/dx at left bottom border D


        p[0:(ny-top_left[0]-2), bottom_left[1]] =\
            p[0:(ny-top_left[0]-2), bottom_left[1]+1] 

        # dp/dx at right top border A


        p[(bottom_right[0]+2):, (nx-bottom_left[1]-1)] =\
            p[(bottom_right[0]+2):, (nx-bottom_left[1]-2)] 

        # dp/dx at left top border B


        p[(ny-top_left[0]):, top_left[1]] =\
            p[(ny-top_left[0]):, top_left[1]+1] 

        # dp/dx at right bottom border C

        p[0:bottom_right[0], bottom_right[1]-1] =\
            p[0:bottom_right[0], bottom_right[1]-2] 


        # dp/dx fill emty space in borders
        # D/E
        # center
        p[(ny-top_left[0]-1), bottom_left[1]] =\
            p[(ny-top_left[0]-2), bottom_left[1]+1] 
        # side
        p[(ny-top_left[0]-1), bottom_left[1]+1] =\
            p[(ny-top_left[0]-2), bottom_left[1]]

        #B/E
        # center
        p[(ny-top_left[0]-1), top_left[1]] =\
            p[(ny-top_left[0]-2), top_left[1]+1] 
        # side
        p[(ny-top_left[0]-1), top_left[1]+1] =\
            p[(ny-top_left[0]-2), top_left[1]] 

        # C/F
        # center
        p[bottom_right[0], bottom_right[1]-1] =\
            p[bottom_right[0]+1, bottom_right[1]-2] 
        # side
        p[bottom_right[0], bottom_right[1]-2] =\
            p[bottom_right[0]+1, bottom_right[1]-1] 

        #A/F
        # center
        p[bottom_right[0], (nx-bottom_left[1]-1)] =\
            p[bottom_right[0]+1, (nx-bottom_left[1]-2)] 
        # side
        p[bottom_right[0], (nx-bottom_left[1]-2)] =\
            p[bottom_right[0]+1, (nx-bottom_left[1]-1)]
            
            
        #electrode
        g = 1
        t = 0
        
        #Top anode
        #right
        p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w)] =\
            p[(ny-top_left[0]-g*h)-2: (ny-top_left[0]-t*h)-2, (top_left[1]-s-w)-1]
        #left
        p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s)-1] =\
            p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s)]
        #top
        p[(ny-top_left[0]-g*h)-1, (top_left[1]-s-w) : (top_left[1]-s)] =\
            p[(ny-top_left[0]-g*h)-2, (top_left[1]-s-w) : (top_left[1]-s)]
        
        #Top cathode
        #right
        p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w-d-w)] =\
            p[(ny-top_left[0]-g*h)-2: (ny-top_left[0]-t*h)-2, ((top_left[1]-s-w-d-w))-1]
        #left
        p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w-d)-1] =\
            p[(ny-top_left[0]-g*h)-2 : (ny-top_left[0]-t*h)-2, (top_left[1]-s-w-d)]
        #top
        p[(ny-top_left[0]-g*h)-1, (top_left[1]-s-w-d-w) : (top_left[1]-s-w-d)] =\
            p[(ny-top_left[0]-g*h)-2, (top_left[1]-s-w-d-w) : (top_left[1]-s-w-d)]
        
        #bottom anode
        #right
        p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s)] =\
            p[(bottom_right[0]+t*h)+2: (bottom_right[0]+g*h)+2, (bottom_right[1]+s)-1]
        #left
        p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w)-1] =\
            p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w)]
        #top
        p[(bottom_right[0]+g*h), (bottom_right[1]+s) : (bottom_right[1]+s+w)] =\
            p[(bottom_right[0]+g*h)+1, (bottom_right[1]+s) : (bottom_right[1]+s+w)]
        
        
        #bottom cathode
        #right
        p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d)] =\
            p[(bottom_right[0]+t*h)+2: (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d)-1]
        #left
        p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d+w)-1] =\
            p[(bottom_right[0]+t*h)+2 : (bottom_right[0]+g*h)+2, (bottom_right[1]+s+w+d+w)]
        #top
        p[(bottom_right[0]+g*h), (bottom_right[1]+s+w+d) : (bottom_right[1]+s+w+d+w)] =\
            p[(bottom_right[0]+g*h)+1, (bottom_right[1]+s+w+d) : (bottom_right[1]+s+w+d+w)]
              
    

    return p

from tqdm import tqdm
import numpy as np
import os
from script.pressure_poisson_periodic import pressure_poisson_periodic, build_up_b
from script.channel_shape import channel_shape

def save_file(file, name, save_folder):
    np.save(os.path.join(save_folder, name), file)

def navier_stokes(u, v, freq,
                  w, h, d, s,
                  bottom_left, top_left, bottom_right,
                  anode_value, cathode_value,
                  border_value,
                  ny, nx,
                  p, b,
                  rho, nu, 
                  dt, dx, dy,
                  F, P0, nt, n_save=10, p_it=1500, save_folder='', stepcount=1,
                  save_vectors=True,
                  silent=False):
    
    
    if save_vectors:
        u_gr = []
        v_gr = []
        p_gr = []
        ch_it = []

   
    
    #check if p initialization is working
    '''    
        b = build_up_b(rho, dt, dx, dy, u, v, w, h, d, s, 
                       bottom_left, bottom_right, top_left, 
                       anode_value, cathode_value, border_value,
                       nx, ny)
        
        p = pressure_poisson_periodic(p,
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
                                  1, dt)'''
    

    for t in tqdm(range(nt)):
        un = u.copy()
        vn = v.copy()
        b = build_up_b(rho, dt, dx, dy, u, v, w, h, d, s, 
                       bottom_left, bottom_right, top_left, 
                       anode_value, cathode_value, border_value,
                       nx, ny)
        p = pressure_poisson_periodic(p,
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
                              p_it, dt)
        
        

        u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx * 
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy * 
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                        dt / (2 * rho * dx) * 
                        (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                        nu * (dt / dx**2 * 
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                        dt / dy**2 * 
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))


        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx * 
                        (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy * 
                        (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) *(p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        nu * (dt / dx**2 *
                        (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 * 
                        (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])) +
                        0*F * dt)


      # Wall  and electrodes BC: u,v = 0 
        v = channel_shape(field=v,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, nx=nx,
                  ny=ny)
        u = channel_shape(field=u,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, nx=nx,
                  ny=ny)
        #v[-1, top_left[1]:] = 0 #F*np.cos(2*np.pi*(stepcount*dt)*freq)  
        #v[0, bottom_left[1]:bottom_right[1]] = 0 # -F*np.cos(2*np.pi*(stepcount*dt)*freq) 
        
        if p.max() > 5:
            if not silent:
                print(f'\nit {t}, b_max : {b.max():.2E}, p_max : {p.max():.2E}, u_max : {u.max():.2E}, v_max : {v.max():.2E}')
            print('Value is too big! Calculation was stopped! ')
            break
        
        
            
        #v[0, bottom_left[1]:bottom_right[1]] = F*np.cos(2*np.pi*(stepcount*dt)*freq)
        #v[-1, top_left[1]:-1] = -F*np.cos(2*np.pi*(stepcount*dt)*freq) 
        if t % n_save == 0:
            if not silent:
                print(f'\nit {t}, b_max : {b.max():.2E}, p_max : {p.max():.2E}, u_max : {u.max():.2E}, v_max : {v.max():.2E}')
            if save_vectors:    
                u_gr.append(u.copy())
                v_gr.append(v.copy())
                p_gr.append(p.copy())
                ch_it.append(t)
            
        stepcount += 1
      
        #udiff = (numpy.sum(u) - numpy.sum(un)) / numpy.sum(u)
        
        
    '''    
    v = channel_shape(field=v,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  andode_value=0, cathode_value=0,
                  border_value=0,
                  ny=ny)
    
    u = channel_shape(field=u,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right,
                  andode_value=0, cathode_value=0,
                  border_value=0,
                  ny=ny)'''
    if not silent:
        print('Calculation complete')
        
    if save_vectors:
        files = [u_gr, v_gr, p_gr, ch_it]
        names = ['u_gr', 'v_gr', 'p_gr', 'ch_it']
        
        for file, name in zip(files, names):
            save_file(file, name, save_folder)
    
        print('Calculated data was saved')
    return u, v, p
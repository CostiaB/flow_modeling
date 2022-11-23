from tqdm import tqdm
import numpy as np
import os
from script.pressure_poisson_periodic import pressure_poisson_periodic, build_up_b
from script.channel_shape import channel_shape
from copy import copy

def save_file(file, name, save_folder):
    np.save(os.path.join(save_folder, name), file)

def navier_stokes(u, v, freq,
                  w, h, d, s,
                  bottom_left, top_left,
                  bottom_right, top_right,
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
        ts = [0]

   
    
    #check if p initialization is working
        
    b = build_up_b(rho, dt, dx, dy, u, v, w, h, d, s, 
                       bottom_left, bottom_right,
                       top_left,  top_right,
                       anode_value, cathode_value, border_value,
                       nx, ny)
    p = pressure_poisson_periodic(p,
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
                              p_it, dt)
    
    
    freq_step = 0
    J = copy(F)
    
    for t in tqdm(range(nt)):
        
        
        #v[315:326, 315:326] = 0 #test
        #u[315:326, 315:326] = 0 #test
        
 
        if t % 10 == 0:
            time_mult = 1/(100 * freq)
            freq_cos = np.cos(2 * np.pi * (freq_step * time_mult) * freq)  
            F = J * freq_cos
            freq_step += 1
            
        
        
        un = u.copy()
        vn = v.copy()
        
        b = build_up_b(rho, dt, dx, dy, u, v, w, h, d, s, 
                       bottom_left, bottom_right,
                       top_left,  top_right,
                       anode_value, cathode_value, border_value,
                       nx, ny)
        p = pressure_poisson_periodic(p,
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
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + F * dt)


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
        
                # Periodic BC u @ x = 2     
        u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx * 
                      (un[1:-1, -1] - un[1:-1, -2]) -
                       vn[1:-1, -1] * dt / dy * 
                      (un[1:-1, -1] - un[0:-2, -1]) -
                       dt / (2 * rho * dx) *
                      (p[1:-1, 0] - p[1:-1, -2]) + 
                       nu * (dt / dx**2 * 
                      (un[1:-1, 0] - 2 * un[1:-1,-1] + un[1:-1, -2]) +
                       dt / dy**2 * 
                      (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)
    
        # Periodic BC u @ x = 0
        u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *
                     (un[1:-1, 0] - un[1:-1, -1]) -
                      vn[1:-1, 0] * dt / dy * 
                     (un[1:-1, 0] - un[0:-2, 0]) - 
                      dt / (2 * rho * dx) * 
                     (p[1:-1, 1] - p[1:-1, -1]) + 
                      nu * (dt / dx**2 * 
                     (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
                      dt / dy**2 *
                     (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)
    
        # Periodic BC v @ x = 2
        v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx *
                      (vn[1:-1, -1] - vn[1:-1, -2]) - 
                       vn[1:-1, -1] * dt / dy *
                      (vn[1:-1, -1] - vn[0:-2, -1]) -
                       dt / (2 * rho * dy) * 
                      (p[2:, -1] - p[0:-2, -1]) +
                       nu * (dt / dx**2 *
                      (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
                       dt / dy**2 *
                      (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))
    
        # Periodic BC v @ x = 0
        v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx *
                     (vn[1:-1, 0] - vn[1:-1, -1]) -
                      vn[1:-1, 0] * dt / dy *
                     (vn[1:-1, 0] - vn[0:-2, 0]) -
                      dt / (2 * rho * dy) * 
                     (p[2:, 0] - p[0:-2, 0]) +
                      nu * (dt / dx**2 * 
                     (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
                      dt / dy**2 * 
                     (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))


      # Wall  and electrodes BC: u,v = 0 
        v = channel_shape(field=v,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right, top_right=top_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, nx=nx,
                  ny=ny)
        u = channel_shape(field=u,
                  w=w, h=h, d=d, s=s,
                  bottom_left=bottom_left, top_left=top_left,
                  bottom_right=bottom_right, top_right=top_right,
                  anode_value=anode_value, cathode_value=cathode_value,
                  border_value=border_value, nx=nx,
                  ny=ny)
        
        
        
        v[-1, top_left[1]:] = 0 #F*np.cos(2*np.pi*(stepcount*dt)*freq)  
        v[0, bottom_left[1]:bottom_right[1]] = 0 # -F*np.cos(2*np.pi*(stepcount*dt)*freq) 
        
        u[-1, top_left[1]:] = 0 #F*np.cos(2*np.pi*(stepcount*dt)*freq)  
        u[0, bottom_left[1]:bottom_right[1]] = 0 # -F*np.cos(2*np.pi*(stepcount*d
        '''
        v[315:326, 315:326] = 0 #test
        u[315:326, 315:326] = 0 #test
        '''
        if any([p.max(), u.max(), v.max()]) > 1:
            if not silent:
                print(f'\nit {t}, b_max : {b.max():.2E}, p_max : {p.max():.2E}, u_max : {u.max():.2E}, v_max : {v.max():.2E}')
            print('Value is too big! Calculation was stopped! ')
            break
        
        
            
        #v[0, bottom_left[1]:bottom_right[1]] = F*np.cos(2*np.pi*(stepcount*dt)*freq)
        #v[-1, top_left[1]:-1] = -F*np.cos(2*np.pi*(stepcount*dt)*freq) 
        
        

    
        
        if t % n_save == 0:
            if not silent:
                print(f'\nit {t}, b_max : {b.max():.2E}, p_max : {p.max():.2E}, u_max : {u.max():.2E}, v_max : {v.max():.2E}')
                
                w = freq/(2*np.pi)
                a = (w/(2*nu))**(1/2)*(1 + 1j)
                eta = nu * rho 
                
                f_17_5 = (450**2/(12*eta) *  (J/1500)/dx)/1e18
                print(f'u_max: {f_17_5 :.2E}, u: {u[100:550, 600].mean():.2E}')
            if save_vectors:    
                u_gr.append(u.copy())
                v_gr.append(v.copy())
                p_gr.append(p.copy())
                ch_it.append(t)
                ts.append(freq_step * time_mult)
            
        stepcount += 1
        if stepcount % 500 == 0:
            udiff = (np.sum(u) - np.sum(un)) / np.sum(u)
            print(udiff)
        
        
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
        files = [u_gr, v_gr, p_gr, ch_it, ts]
        names = ['u_gr', 'v_gr', 'p_gr', 'ch_it', 'ts']
        
        for file, name in zip(files, names):
            save_file(file, name, save_folder)
    
        print('Calculated data was saved')
    return u, v, p
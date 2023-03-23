import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from tqdm import tqdm
import numpy as np
import torch
import os
from script.pressure_poisson_periodic import pressure_poisson_periodic, build_up_b
from script.channel_shape import channel_shape
from copy import copy
from collections import deque


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
                  F, P0, nt, max_vec_len,
                  n_save=10, p_it=1500,  save_folder='', stepcount=1,
                  save_vectors=True,
                  silent=False,
                  sine=False,
                  freq_points=1000
                  ):
    
    
    if save_vectors:
        u_gr = deque(maxlen=max_vec_len)
        v_gr = deque(maxlen=max_vec_len)
        p_gr = deque(maxlen=max_vec_len)
        ch_it = deque(maxlen=max_vec_len)
        ts = deque([0], maxlen=max_vec_len)
        
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
    time_mult = freq / 100
    J = copy(F)
    
    for t in tqdm(range(nt)):
        
        if sine:
            if t % freq_points == 0:
                freq_cos = np.cos(2 * np.pi * (freq_step * time_mult))
                F = J * freq_cos
                freq_step += 1
           
        
        un = u.clone()
        vn = v.clone()
        
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
        
        

       #center
       #y
        c_y = (bottom_right[0], ny-top_left[0])
        '''
              c_y[0]+1:c_y[1]-1 #1:-1
              c_y[0]:c_y[1]-2 # 0:-2
              c_y[0]+2:c_y[1] # 2:
        '''

        c_x = (bottom_right[1], top_left[1])
        '''
       c_x[0]:c_x[1] #1:-1 need to check; now it was made one point wider
       c_x[0]-1:c_x[1]-1 #0:-2
       c_x[0]+1:c_x[1]+1 #2:
        '''       
                          
     

        u[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] = (un[c_y[0]+1:c_y[1]-1 , c_x[0]-1:c_x[1]+1 ] -
                          un[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] * dt / dx * 
                          (un[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] - un[c_y[0]+1:c_y[1]-1, c_x[0]-2:c_x[1]]) -
                          vn[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] * dt / dy * 
                          (un[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] - un[c_y[0]:c_y[1]-2, c_x[0]-1:c_x[1]+1 ]) -
                          dt / (2 * rho * dx) * 
                          (p[c_y[0]+1:c_y[1]-1, c_x[0]:c_x[1]+2] - p[c_y[0]+1:c_y[1]-1, c_x[0]-2:c_x[1]]) +
                          nu * (dt / dx**2 * 
                          (un[c_y[0]+1:c_y[1]-1, c_x[0]:c_x[1]+2] - 2 * un[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] + un[c_y[0]+1:c_y[1]-1, c_x[0]-2:c_x[1]]) +
                          dt / dy**2 * 
                          (un[c_y[0]+2:c_y[1], c_x[0]-1:c_x[1]+1 ] - 2 * un[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1 ] + un[c_y[0]:c_y[1]-2, c_x[0]-1:c_x[1]+1 ])) + 
                           F * dt)

        v[c_y[0]+1:c_y[1]-1 , c_x[0]-1:c_x[1]+1] = (vn[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1] -
                        un[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1] * dt / dx * 
                        (vn[c_y[0]+1:c_y[1]-1 , c_x[0]-1:c_x[1]+1] - vn[c_y[0]+1:c_y[1]-1 , c_x[0]-2:c_x[1]]) -
                        vn[c_y[0]+1:c_y[1]-1 , c_x[0]-1:c_x[1]+1] * dt / dy * 
                        (vn[c_y[0]+1:c_y[1]-1, c_x[0]-1:c_x[1]+1] - vn[c_y[0]:c_y[1]-2, c_x[0]-1:c_x[1]+1]) -
                        dt / (2 * rho * dy) *(p[c_y[0]+2:c_y[1], c_x[0]-1:c_x[1]+1] - p[c_y[0]:c_y[1]-2, c_x[0]-1:c_x[1]+1]) +
                        nu * (dt / dx**2 *
                        (vn[c_y[0]+1:c_y[1]-1 , c_x[0]:c_x[1]+2] - 2 * vn[c_y[0]+1:c_y[1]-1 , c_x[0]-1:c_x[1]+1] + vn[c_y[0]+1:c_y[1]-1 , c_x[0]-2:c_x[1] ]) +
                        dt / dy**2 * 
                        (vn[c_y[0]+2:c_y[1], c_x[0]-1:c_x[1]+1] - 2 * vn[c_y[0]+1:c_y[1]-1 , c_x[0]-1:c_x[1]+1] + vn[c_y[0]:c_y[1]-2, c_x[0]-1:c_x[1]+1])) +
                        0 * F * dt)



       #left part
        l_y = (0, ny-top_left[0])
        l_x = (bottom_left[1], bottom_right[1])
        '''
       #y
       l_y[0]+1:l_y[1]-1 #1:-1
       l_y[0]:l_y[1]-2 # 0:-2
       l_y[0]+2:l_y[1] # 2:
       #x

       l_x[0]+1:l_x[1]-1 #1:-1 need to check; now it was made one point wider
       l_x[0]:l_x[1]-2 #0:-2
       l_x[0]+2:l_x[1] #2:
        '''
        u[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] = (un[l_y[0]+1:l_y[1]-1 , l_x[0]+1:l_x[1]-1 ] -
                          un[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] * dt / dx * 
                          (un[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] - un[l_y[0]+1:l_y[1]-1, l_x[0]:l_x[1]-2]) -
                          vn[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] * dt / dy * 
                          (un[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] - un[l_y[0]:l_y[1]-2, l_x[0]+1:l_x[1]-1 ]) -
                          dt / (2 * rho * dx) * 
                          (p[l_y[0]+1:l_y[1]-1, l_x[0]+2:l_x[1]] - p[l_y[0]+1:l_y[1]-1, l_x[0]:l_x[1]-2]) +
                          nu * (dt / dx**2 * 
                          (un[l_y[0]+1:l_y[1]-1, l_x[0]+2:l_x[1]] - 2 * un[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] + un[l_y[0]+1:l_y[1]-1, l_x[0]:l_x[1]-2]) +
                          dt / dy**2 * 
                          (un[l_y[0]+2:l_y[1], l_x[0]+1:l_x[1]-1 ] - 2 * un[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1 ] + un[l_y[0]:l_y[1]-2, l_x[0]+1:l_x[1]-1 ])) + 
                          0* F * dt)

        v[l_y[0]+1:l_y[1]-1 , l_x[0]+1:l_x[1]-1] = (vn[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1] -
                        un[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1] * dt / dx * 
                        (vn[l_y[0]+1:l_y[1]-1 , l_x[0]+1:l_x[1]-1] - vn[l_y[0]+1:l_y[1]-1 , l_x[0]:l_x[1]-2 ]) -
                        vn[l_y[0]+1:l_y[1]-1 , l_x[0]+1:l_x[1]-1] * dt / dy * 
                        (vn[l_y[0]+1:l_y[1]-1, l_x[0]+1:l_x[1]-1] - vn[l_y[0]:l_y[1]-2, l_x[0]+1:l_x[1]-1]) -
                        dt / (2 * rho * dy) *(p[l_y[0]+2:l_y[1], l_x[0]+1:l_x[1]-1] - p[l_y[0]:l_y[1]-2, l_x[0]+1:l_x[1]-1]) +
                        nu * (dt / dx**2 *
                        (vn[l_y[0]+1:l_y[1]-1 , l_x[0]+2:l_x[1]] - 2 * vn[l_y[0]+1:l_y[1]-1 , l_x[0]+1:l_x[1]-1] + vn[l_y[0]+1:l_y[1]-1 , l_x[0]:l_x[1]-2 ]) +
                        dt / dy**2 * 
                        (vn[l_y[0]+2:l_y[1], l_x[0]+1:l_x[1]-1] - 2 * vn[l_y[0]+1:l_y[1]-1 , l_x[0]+1:l_x[1]-1] + vn[l_y[0]:l_y[1]-2, l_x[0]+1:l_x[1]-1])) +
                        F * dt)




       #right part
        r_y = (bottom_right[0], ny)  # found mistake x and y were swapped
        r_x = (top_left[1], top_right[1]) 
        '''
       #y
       r_y[0]+1:r_y[1]-1 #1:-1
       r_y[0]:r_y[1]-2 # 0:-2
       r_y[0]+2:r_y[1] # 2:
       #x

       r_x[0]+1:r_x[1]-1 #1:-1 need to check; now it was made one point wider
       r_x[0]:r_x[1]-2 #0:-2
       r_x[0]+2:r_x[1] #2:
        '''
        u[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] = (un[r_y[0]+1:r_y[1]-1 , r_x[0]+1:r_x[1]-1 ] -
                          un[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] * dt / dx * 
                          (un[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] - un[r_y[0]+1:r_y[1]-1, r_x[0]:r_x[1]-2]) -
                          vn[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] * dt / dy * 
                          (un[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] - un[r_y[0]:r_y[1]-2, r_x[0]+1:r_x[1]-1 ]) -
                          dt / (2 * rho * dx) * 
                          (p[r_y[0]+1:r_y[1]-1, r_x[0]+2:r_x[1]] - p[r_y[0]+1:r_y[1]-1, r_x[0]:r_x[1]-2]) +
                          nu * (dt / dx**2 * 
                          (un[r_y[0]+1:r_y[1]-1, r_x[0]+2:r_x[1]] - 2 * un[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] + un[r_y[0]+1:r_y[1]-1, r_x[0]:r_x[1]-2]) +
                          dt / dy**2 * 
                          (un[r_y[0]+2:r_y[1], r_x[0]+1:r_x[1]-1 ] - 2 * un[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1 ] + un[r_y[0]:r_y[1]-2, r_x[0]+1:r_x[1]-1 ])) + 
                          0* F * dt)

        v[r_y[0]+1:r_y[1]-1 , r_x[0]+1:r_x[1]-1] = (vn[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1] -
                        un[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1] * dt / dx * 
                        (vn[r_y[0]+1:r_y[1]-1 , r_x[0]+1:r_x[1]-1] - vn[r_y[0]+1:r_y[1]-1 , r_x[0]:r_x[1]-2 ]) -
                        vn[r_y[0]+1:r_y[1]-1 , r_x[0]+1:r_x[1]-1] * dt / dy * 
                        (vn[r_y[0]+1:r_y[1]-1, r_x[0]+1:r_x[1]-1] - vn[r_y[0]:r_y[1]-2, r_x[0]+1:r_x[1]-1]) -
                        dt / (2 * rho * dy) *(p[r_y[0]+2:r_y[1], r_x[0]+1:r_x[1]-1] - p[r_y[0]:r_y[1]-2, r_x[0]+1:r_x[1]-1]) +
                        nu * (dt / dx**2 *
                        (vn[r_y[0]+1:r_y[1]-1 , r_x[0]+2:r_x[1]] - 2 * vn[r_y[0]+1:r_y[1]-1 , r_x[0]+1:r_x[1]-1] + vn[r_y[0]+1:r_y[1]-1 , r_x[0]:r_x[1]-2 ]) +
                        dt / dy**2 * 
                        (vn[r_y[0]+2:r_y[1], r_x[0]+1:r_x[1]-1] - 2 * vn[r_y[0]+1:r_y[1]-1 , r_x[0]+1:r_x[1]-1] + vn[r_y[0]:r_y[1]-2, r_x[0]+1:r_x[1]-1])) +
                        F * dt)
        
                # Periodic BC u @ y = top  
                #r_x[1]
                #l_x[0]+1:l_x[1]-1 #1:-1 need to check; now it was made one point wider
                #l_x[0]:l_x[1]-2 #0:-2
                #l_x[0]+2:l_x[1] #2:
                # r_x[0]+1:r_x[1]-1 #1:-1 need to check; now it was made one point wider
                # r_x[0]:r_x[1]-2 #0:-2
                # r_x[0]+2:r_x[1] #2:
                
        u[-1, r_x[0]+1:r_x[1]-1] = (un[-1, r_x[0]+1:r_x[1]-1] - un[-1, r_x[0]+1:r_x[1]-1] * dt / dx * 
                      (un[-1, r_x[0]+1:r_x[1]-1] - un[-2, r_x[0]+1:r_x[1]-1]) -
                       vn[-1, r_x[0]+1:r_x[1]-1] * dt / dy * 
                      (un[-1, r_x[0]+1:r_x[1]-1] - un[-1, r_x[0]:r_x[1]-2 ]) -
                       dt / (2 * rho * dx) *
                      (p[0, l_x[0]+1:l_x[1]-1] - p[-2, r_x[0]+1:r_x[1]-1]) + 
                       nu * (dt / dx**2 * 
                      (un[0, l_x[0]+1:l_x[1]-1] - 2 * un[-1, r_x[0]+1:r_x[1]-1] + un[-2, r_x[0]+1:r_x[1]-1]) +
                       dt / dy**2 * 
                      (un[-1, r_x[0]+2:r_x[1]] - 2 * un[-1, r_x[0]+1:r_x[1]-1] + un[-1, r_x[0]:r_x[1]-2 ])))
    
        # Periodic BC u @ y = 0
        u[0, l_x[0]+1:l_x[1]-1] = (un[0, l_x[0]+1:l_x[1]-1] - un[0, l_x[0]+1:l_x[1]-1] * dt / dx *
                    (un[0, l_x[0]+1:l_x[1]-1] - un[-1, r_x[0]+1:r_x[1]-1]) -
                     vn[0, l_x[0]+1:l_x[1]-1] * dt / dy * 
                    (un[0, l_x[0]+1:l_x[1]-1] - un[0, r_x[0]:r_x[1]-2 ]) - 
                     dt / (2 * rho * dx) * 
                    (p[1, l_x[0]+1:l_x[1]-1] - p[-1, r_x[0]+1:r_x[1]-1]) + 
                     nu * (dt / dx**2 * 
                    (un[1, l_x[0]+1:l_x[1]-1] - 2 * un[0, l_x[0]+1:l_x[1]-1] + un[-1, r_x[0]+1:r_x[1]-1]) +
                     dt / dy**2 *
                    (un[0, r_x[0]+2:r_x[1]] - 2 * un[0, l_x[0]+1:l_x[1]-1] + un[0, r_x[0]:r_x[1]-2 ])))
    
        # Periodic BC v @ y = top
        v[-1, r_x[0]+1:r_x[1]-1 ] = (vn[-1, r_x[0]+1:r_x[1]-1 ] - un[-1, r_x[0]+1:r_x[1]-1 ] * dt / dx *
                      (vn[-1, r_x[0]+1:r_x[1]-1 ] - vn[-2, r_x[0]+1:r_x[1]-1 ]) - 
                       vn[-1, r_x[0]+1:r_x[1]-1 ] * dt / dy *
                      (vn[-1, r_x[0]+1:r_x[1]-1 ] - vn[-1, r_x[0]:r_x[1]-2 ]) -
                       dt / (2 * rho * dy) * 
                      (p[-1, r_x[0]+2:r_x[1] ] - p[-1, r_x[0]:r_x[1]-2 ]) +
                       nu * (dt / dx**2 *
                      (vn[0, l_x[0]+1:l_x[1]-1] - 2 * vn[-1, r_x[0]+1:r_x[1]-1 ] + vn[-2, r_x[0]+1:r_x[1]-1 ]) +
                       dt / dy**2 *
                      (vn[-1, r_x[0]+2:r_x[1] ] - 2 * vn[-1, r_x[0]+1:r_x[1]-1 ] + vn[-1, r_x[0]:r_x[1]-2 ]))+ F * dt)
    
        # Periodic BC v @ y = 0
        v[0, l_x[0]+1:l_x[1]-1] = (vn[0, l_x[0]+1:l_x[1]-1] - un[0, l_x[0]+1:l_x[1]-1] * dt / dx *
                     (vn[0, l_x[0]+1:l_x[1]-1] - vn[-1, r_x[0]+1:r_x[1]-1 ]) -
                      vn[0, l_x[0]+1:l_x[1]-1] * dt / dy *
                     (vn[0, l_x[0]+1:l_x[1]-1] - vn[0, l_x[0]:l_x[1]-2]) -
                      dt / (2 * rho * dy) * 
                     (p[0, l_x[0]+2:l_x[1]] - p[0, l_x[0]:l_x[1]-2]) +
                      nu * (dt / dx**2 * 
                     (vn[1, l_x[0]+1:l_x[1]-1] - 2 * vn[0, l_x[0]+1:l_x[1]-1] + vn[-1, r_x[0]+1:r_x[1]-1 ]) +
                      dt / dy**2 * 
                     (vn[0, l_x[0]+2:l_x[1]] - 2 * vn[0, l_x[0]+1:l_x[1]-1] + vn[0, l_x[0]:l_x[1]-2])) + F * dt)

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
        
        if any([p.max(), u.max(), v.max()]) > 1:
            if not silent:
                print(f'\nit {t}, b_max : {b.max():.2E}, p_max : {p.max():.2E}, u_max : {u.max():.2E}, v_max : {v.max():.2E}')
            print('Value is too big! Calculation was stopped! ')
            break
        
        if t % n_save == 0:
            if not silent:
                print(f'\nit {t}, b_max : {b.max():.2E}, p_max : {p.max():.2E}, u_max : {u.max():.2E}, v_max : {v.max():.2E}')
            if save_vectors:    
                u_gr.append(u.clone().cpu())
                v_gr.append(v.clone().cpu())
                p_gr.append(p.clone().cpu())
                ch_it.append(t)
                ts.append(freq_step * time_mult)
        stepcount += 1
        if stepcount % 500 == 0:
            udiff = (torch.sum(u) - torch.sum(un)) / torch.sum(u)
            print(udiff)

    if not silent:
        print('Calculation complete')
        
    if save_vectors:
        files = [u_gr, v_gr, p_gr, ch_it, ts]
        names = ['u_gr', 'v_gr', 'p_gr', 'ch_it', 'ts']
        
        for file, name in zip(files, names):
            save_file(file, name, save_folder)
    
        print('Calculated data was saved')
    return u, v, p

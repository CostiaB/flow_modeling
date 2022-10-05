def main(N_it, p_it, freq, stop_it):
    import numpy as np
    import gc
    import os
    import time
    import shutil
    import matplotlib.pyplot as plt
    from tqdm import tqdm

    from script.navier_stokes import navier_stokes
    from script.channel_shape import channel_shape
    from script.laplace import laplace2d_parralel_C0
    from script.plot_gif import plot_gif
    import params_channel_shape
    import params_to_calc
    import params_conc
    
    dx = params_to_calc.dx
    dy = params_to_calc.dy

    dt = params_to_calc.dt

    rho = params_to_calc.rho
    nu = params_to_calc.nu
    D = params_to_calc.D

    F = params_to_calc.F
    P0 = params_to_calc.P0
    
    nx = params_channel_shape.nx
    ny = params_channel_shape.ny

    field = params_channel_shape.field
    anode_value = params_channel_shape.anode_value
    cathode_value = params_channel_shape.cathode_value
    border_value = params_channel_shape.border_value
    w = params_channel_shape.w
    h = params_channel_shape.h
    d = params_channel_shape.d
    s = params_channel_shape.s
    bottom_left = params_channel_shape.bottom_left
    top_left = params_channel_shape.top_left
    bottom_right = params_channel_shape.bottom_right

    v = params_channel_shape.v
    u = params_channel_shape.u
    p = params_channel_shape.p
    b = params_channel_shape.b

    c0 = params_conc.c0
    ca = params_conc.ca
    c = params_conc.c
    c_anode_value = params_conc.c_anode_value
    c_cathode_value = params_conc.c_cathode_value
    resp = []
    ch_it = []
    c_arr = []

    try:
        
        n_save = N_it // 200 if N_it > 200 else 1
        save_folder = './results/'
        save_folder += time.strftime('%Y_%m_%d_%H_%M', time.localtime())
        if not os.path.exists(save_folder):
            os.mkdir(save_folder)
            
         
        for a in tqdm(range(N_it)):
            u, v, p = navier_stokes(u, v, freq,
                              w, h, d, s,
                              bottom_left, top_left, bottom_right,
                              anode_value, cathode_value,
                              border_value,
                              ny, nx,
                              p, b,
                              rho, nu, 
                              dt, dx, dy,
                              F, P0, 1, n_save, p_it, save_folder=save_folder,
                              stepcount=a, save_vectors=False)


            ''' 
            conc =  laplace2d_parralel_C0(c, 
                                  freq,
                                  w, h, d, s,
                                  bottom_left, top_left, bottom_right,
                                  c_anode_value, c_cathode_value,
                                  border_value,
                                  ny, nx,
                                  b,
                                  rho, nu, 
                                  dx, dy,
                                  P0=P0, stepcount=a,
                                  nit=77, l1norm_target=1,
                                  stop_it=stop_it)
    
             c = conc.copy()
            '''
    
         ##loop across number of time steps
            for it in range(1000):
                cn = c.copy()
                
                #print('A: ', D * dt / dx**2 * (cn[1:-1, 2:] - 2 * cn[1:-1, 1:-1] + cn[1:-1, 0:-2]).mean())
                #print('B: ', D * dt / dy**2 * (cn[2:, 1:-1] - 2 * cn[1:-1, 1:-1] + cn[0:-2, 1:-1]).mean())
                #print('C: ', ( (u[1:-1,1:-1]* dt / dx * (cn[1:-1, 1:-1] - cn[1:-1, 0:-2])) ).mean())
                #print('D: ', (v[1:-1,1:-1] * dt / dy * (cn[1:-1, 1:-1] - cn[0:-2, 1:-1])).mean())
                
                c[1:-1, 1:-1] = cn[1:-1,1:-1] - \
                                    D * dt / dx**2 * \
                                    (cn[1:-1, 2:] - 2 * cn[1:-1, 1:-1] + cn[1:-1, 0:-2]) - \
                                    D * dt / dy**2 * \
                                    (cn[2:, 1:-1] - 2 * cn[1:-1, 1:-1] + cn[0:-2, 1:-1]) + \
                                    ( (u[1:-1,1:-1]* dt / dx * (cn[1:-1, 1:-1] - cn[1:-1, 0:-2])) )+ \
                                        (v[1:-1,1:-1] * dt / dy * (cn[1:-1, 1:-1] - cn[0:-2, 1:-1]))
            
                    
    
                c = channel_shape(field=c,
                            w=w, h=h, d=d, s=s,
                            bottom_left=bottom_left, top_left=top_left,
                            bottom_right=bottom_right,
                            anode_value=ca, cathode_value=0,
                            border_value=0,
                            nx=nx,
                            ny=ny)
            
            
                
            if a % n_save == 0:
                c_arr.append(c)
                ch_it.append(a)
                
            t = 0
            g = 1
            if h==0:
                h = 1
            elif h<0:
                t = 1
                g = 0
            #bottom cathode        
            bc = c[ny-top_left[0]-g*h-1:ny-top_left[0]-t*h, top_left[1]-s-w-d-w-1:top_left[1]-s-w-d+1]
            #top cathode
            tc = c[bottom_right[0]+t*h:bottom_right[0]+g*h+1, bottom_right[1]+s+w+d-1:bottom_right[1]+s+w+d+w+1]
            
    
    
            integr = tc.sum() - bc.sum()
            print(f'integration {freq} : {integr:.2E}')
            resp.append(integr)
            
            
        fig = plt.figure(figsize=(50, 20), dpi=200)
        plt.plot(resp)
        plt.savefig(save_folder + '/conc.png')
        print('Consentration plot was saved')
        
        p_gif = channel_shape(field = field,
                          w=w, h=h, d=d, s=s,
                          bottom_left=bottom_left, top_left=top_left,
                          bottom_right=bottom_right,
                          anode_value=anode_value+10000, 
                          cathode_value=cathode_value+10000,
                          border_value=border_value+10000, nx=nx,
                          ny=ny)
        
        name = 'concentrarion'
        np.save(save_folder +'/c_arr', c_arr)
        print('Concentration data was saved')
        gc.collect()
        imgs = np.load(save_folder +'/c_arr.npy')
        vmin = 0 if len(imgs[imgs<0]) == 0 else np.mean(imgs[imgs<0])
        vmax = 0 if len(imgs[imgs>0]) == 0 else np.mean(imgs[imgs>0])
        if 0 not in [vmin, vmax]:    
            mean = np.mean(np.abs([vmin, vmax]))
            vmin, vmax = -mean, mean
        if vmax > 0.2:
            vmin, vmax = -0.02, 0.02
        plot_gif(imgs, p_gif, name, vmin, vmax, ch_it, save_folder)
        gc.collect()
        print(f'Gif {name} was saved')
        
        shutil.copy('params_to_calc.py', save_folder)
        shutil.copy('params_channel_shape.py', save_folder)
        shutil.copy('params_conc.py', save_folder)
        print('Params was saved')
        print(f'Directory name: {save_folder}')
        
        
        
        gc.collect()
    except KeyboardInterrupt:
        os.rmdir(save_folder)
        print(f'\nNothing to save! Folder {save_folder} was removed')
        
if __name__ == '__main__':
    
    import argparse
    
    N_it = 200
    p_it = 1500
    freq = .1
    stop_it = 3000
    
    parser = argparse.ArgumentParser(description='Create params for calculation')
    parser.add_argument('--N_it', default=N_it, type=int)
    parser.add_argument('--p_it', default=p_it, type=int)
    parser.add_argument('--freq', default=freq, type=float)
    parser.add_argument('--stop_it', default=stop_it, type=float)
    
    args = parser.parse_args()
    main(N_it=args.N_it, p_it=args.p_it, freq=args.freq, stop_it=stop_it)
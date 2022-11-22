def main(nt, p_it, freq):
    import numpy as np
    import torch
    import gc
    import os
    import time
    import shutil

    from script.navier_stokes import navier_stokes
    from script.channel_shape import channel_shape
    from script.plot_gif import plot_gif
    from script.save_params import write_params 
    
    import params_channel_shape
    import params_to_calc
    
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
    top_right = params_channel_shape.top_right

    v = params_channel_shape.v
    u = params_channel_shape.u
    p = params_channel_shape.p
    b = params_channel_shape.b

    
    

    try:
        max_vec_len = 260
        n_save = nt // max_vec_len if nt > max_vec_len else 1
        save_folder = './results/'
        save_folder += 'navier_'
        save_folder += time.strftime('%Y_%m_%d_%H_%M', time.localtime())
        if not os.path.exists(save_folder):
            os.mkdir(save_folder)
            
        u, v, p = navier_stokes(u, v, freq,
                          w, h, d, s,
                          bottom_left, top_left,
                          bottom_right, top_right,
                          anode_value, cathode_value,
                          border_value,
                          ny, nx,
                          p, b,
                          rho, nu, 
                          dt, dx, dy,
                          F, P0, nt, max_vec_len, n_save,
                          p_it, save_folder=save_folder, save_vectors=True)
        
        p_gif = channel_shape(field = field,
                          w=w, h=h, d=d, s=s,
                          bottom_left=bottom_left, top_left=top_left,
                          bottom_right=bottom_right, top_right=top_right,
                          anode_value=anode_value+10000, cathode_value=cathode_value+10000,
                          border_value=border_value+10000, nx=nx,
                          ny=ny).to('cpu').numpy()
        np.save(save_folder + '/' + 'channel_for_gif', p.to('cpu'))
        
        names = ['p_init_0', 'u_init_0', 'v_init_0']
        files = ['p_gr.npy', 'u_gr.npy', 'v_gr.npy']

        names = ['u_init_0']
        files = ['u_gr.npy']

        ch_it = np.load(save_folder + '/ch_it.npy').tolist()
        
        
        for file, name in zip(files, names):
            gc.collect()
            imgs = np.load(save_folder + '/' + file, allow_pickle=True)
            imgs = np.array([i.numpy() for i in imgs])
            vmin = 0 if len(imgs[imgs<0]) == 0 else np.mean(imgs[imgs<0])
            vmax = 0 if len(imgs[imgs>0]) == 0 else np.mean(imgs[imgs>0])
        
            if 0 not in [vmin, vmax]:    
                mean = np.mean(np.abs([vmin, vmax]))
                vmin, vmax = -mean, mean
                
            if vmax > 0.2:
                vmin, vmax = -0.02, 0.02
            #vmin, vmax = 1e-19, -1e-19
            plot_gif(imgs, p_gif, name, vmin, vmax, ch_it, save_folder)
            gc.collect()
            print(f'Gif {name} was saved')
        
        shutil.copy('params_to_calc.py', save_folder)
        shutil.copy('params_channel_shape.py', save_folder)
        write_params(save_folder, p_it=p_it, nt=nt, freq=freq)
        print('Params was saved')
        print(f'Directory name: {save_folder}')
        
        
        import matplotlib.pyplot as plt
        x = np.linspace(0, 2, nx)
        y = np.linspace(0, 2, ny)
        X, Y = np.meshgrid(x, y)
        
        u = u.to('cpu')
        v = v.to('cpu')
        

        fig = plt.figure(figsize = (20,15), dpi=400)
        st = 15
        plt.quiver(X[::st, ::st], Y[::st, ::st], u[::st, ::st],
                   v[::st, ::st], width=1e-3)
        plt.savefig(save_folder + '/final_speed.png')
        
        fig = plt.figure(figsize = (20,15), dpi=400)
        st = 1
        plt.quiver(X[200:400:st, 200:400:st], Y[200:400:st, 200:400:st],
                   u[200:400:st, 200:400:st], v[200:400:st, 200:400:st], width=1e-3)
        plt.savefig(save_folder + '/final_speed_center.png')
        
        
        fig = plt.figure(figsize = (20,15))
        plt.plot(u[:, 500])
        plt.savefig(save_folder + '/u_slice.png')
        
        
        
        gc.collect()
    except KeyboardInterrupt:
        os.rmdir(save_folder)
        print(f'\nNothing to save! Folder {save_folder} was removed')
        
if __name__ == '__main__':
    
    import argparse
    
    nt = 20000
    p_it = 50
    freq = .1
    
    parser = argparse.ArgumentParser(description='Create params for calculation')
    parser.add_argument('--nt', default=nt, type=int)
    parser.add_argument('--p_it', default=p_it, type=int)
    parser.add_argument('--freq', default=freq, type=float)
    
    args = parser.parse_args()
    main(nt=args.nt, p_it=args.p_it, freq=args.freq)
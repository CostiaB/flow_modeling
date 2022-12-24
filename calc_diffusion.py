def plot_conc(c, ca,
              w, h, d, s,
              bottom_left, top_left,
              bottom_right, top_right,
              nx, ny, save_folder):

    import matplotlib.pyplot as plt
    import numpy as np
    
    fig, ax = plt.subplots(figsize=(16, 16))
    plt.imshow(c[250:400, 200:1200])
    plt.colorbar()
    plt.savefig(save_folder + '/conc_ch.png')

    print(ca)	
    scale = 0.0001

    x, y = bottom_left[1]+2, top_left[1]

    tmp = c[250:400, x:y]

    a = np.zeros(c.shape[1])
    a[x:y] = ca
    a[(top_left[1]-s-w) : (top_left[1]-s)] = ca - (h+ 1)* scale
    a[(top_left[1]-s-w-d-w) : (top_left[1]-s-w-d)] = ca - (h+ 1)* scale
    a = a[x:y]

    fig, ax = plt.subplots(figsize=(16,16))
    plt.plot(tmp[:75, :].mean(axis=0))
    plt.plot(a)
    plt.ylim(0, ca*1.05)
    plt.savefig(save_folder + '/conc_top.png')


    fig, ax = plt.subplots(figsize=(16,16))
    plt.plot(c[bottom_right[0], x:y])
    plt.plot(a)
    plt.ylim(0, ca*1.05)
    plt.savefig(save_folder + '/conc_top_border.png')


    x, y = bottom_right[1], (nx-bottom_left[1]-2)

    tmp = c[250:400, x:y]


    a = np.zeros(c.shape[1])
    a[x:y] = 1 * scale 
    a[(bottom_right[1]+s) : (bottom_right[1]+s+w)] = (1+h) * scale
    a[(bottom_right[1]+s+w+d) : (bottom_right[1]+s+w+d+w)] = (1+h) * scale
    a = a[x:y]

    fig, ax = plt.subplots(figsize=(16,16))
    plt.plot(tmp[75:, :].mean(axis=0))
    plt.plot(a)
    plt.ylim(0, ca*1.05)
    plt.savefig(save_folder + '/conc_bot.png')

    fig, ax = plt.subplots(figsize=(16,16))
    plt.plot(c[(ny-top_left[0]-1), x:y])
    plt.plot(a)
    plt.ylim(0, ca*1.05)
    plt.savefig(save_folder + '/conc_bot_border.png')



def main(N_it, stop_it, l1norm_target, sine, on_bot):
    import numpy as np
    import gc
    import os
    import time
    import shutil
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    import torch

    from script.channel_shape import channel_shape
    from script.laplace import laplace2d_parralel_C0
    from script.plot_gif import plot_gif
    from script.save_params import write_params

    import params_channel_shape
    import params_to_calc
    import params_conc

    dx = params_to_calc.dx
    dy = params_to_calc.dy

    dt = params_to_calc.dt

    rho = params_to_calc.rho
    nu = params_to_calc.nu
    D = params_to_calc.D

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
    b = params_channel_shape.b

    ca = params_conc.ca
    c = params_conc.c
    c_anode_value = params_conc.c_anode_value
    c_cathode_value = params_conc.c_cathode_value
    resp = []
    ch_it = []
    c_arr = []
    freq = 0

    try:
        print('Calculation start with current params:')
        print(
            f'N_it: {N_it}, stop_it: {stop_it}, l1norm_target: {l1norm_target}')

        save_folder = './results/'
        if not os.path.exists(save_folder):
            os.mkdir(save_folder)
        save_folder += 'diffusion_'
        save_folder += time.strftime('%Y_%m_%d_%H_%M', time.localtime())
        if not os.path.exists(save_folder):
            os.mkdir(save_folder)

        if sine == False:
            conc = laplace2d_parralel_C0(c,
                                         w, h, d, s,
                                         bottom_left, top_left,
                                         bottom_right, top_right,
                                         c_anode_value, c_cathode_value,
                                         border_value,
                                         ny, nx,
                                         b,
                                         rho, nu,
                                         dx, dy,
                                         P0=P0, stepcount=1,
                                         l1norm_target=l1norm_target,
                                         stop_it=stop_it, on_bot=on_bot)

            p_gif = channel_shape(field=field,
                                  w=w, h=h, d=d, s=s,
                                  bottom_left=bottom_left, top_left=top_left,
                                  bottom_right=bottom_right, top_right=top_right,
                                  anode_value=anode_value+10000,
                                  cathode_value=cathode_value+10000,
                                  border_value=border_value+10000, nx=nx,
                                  ny=ny, on_bot=on_bot).to('cpu').numpy()

            conc = conc.to('cpu').numpy()
            conc = conc[::-1, :]
            
            np.save(save_folder + '/conc', conc)
            print('Concentration data was saved')
            
            vmin = conc.min()
            vmax = conc.max()
            plt.imshow(conc)
            plt.colorbar()
            plt.savefig(save_folder + '/conc.png')
            plot_conc(conc, ca,
                      w, h, d, s,
                      bottom_left, top_left,
                      bottom_right, top_right,
                      nx, ny, save_folder)

            print('Consentration plots was saved')
            gc.collect()

        else:
            gc.collect()

            inp = input("Type path to navier data: ")
            if len(inp) == 0:
                print('Empty input!')
                path = './results/navier_2022_11_24_23_12/'
            else:
                path = inp
            print(f'Data source: {path}')
            u_vec = np.load(os.path.join(path, 'u_gr.npy'), allow_pickle=True)
            v_vec = np.load(os.path.join(path, 'v_gr.npy'), allow_pickle=True)
            
            import sys
            sys.path.append(path)
            import start_params
            
            nt = start_params.nt
            freq = start_params.freq
            freq_points = start_params.freq_points
            f_steps = 100
            t_calc = (nt / (freq_points * f_steps))/freq
            vec_len = len(u_vec)
            dt = (t_calc/vec_len)/N_it
            print(f'freq = {freq}; dt = {dt}')

            for a in tqdm(range(len(u_vec))):
                
                u = u_vec[a].cuda()
                v = v_vec[a].cuda()

                conc = laplace2d_parralel_C0(c,
                                             w, h, d, s,
                                             bottom_left, top_left,
                                             bottom_right, top_right,
                                             c_anode_value, c_cathode_value,
                                             border_value,
                                             ny, nx,
                                             b,
                                             rho, nu,
                                             dx, dy,
                                             P0=P0, stepcount=a,
                                             l1norm_target=l1norm_target,
                                             stop_it=stop_it,
                                             silent=True, on_bot=on_bot)

                c = conc.clone()
                gc.collect()
                # loop across number of time steps
                for it in range(N_it):

                    #cn = c.copy()
                    cn = c.clone()
                    c[1:-1, 1:-1] = cn[1:-1, 1:-1] - \
                        D * dt / dx**2 * \
                        (cn[1:-1, 2:] - 2 * cn[1:-1, 1:-1] + cn[1:-1, 0:-2]) - \
                        D * dt / dy**2 * \
                        (cn[2:, 1:-1] - 2 * cn[1:-1, 1:-1] + cn[0:-2, 1:-1]) + \
                        ((u[1:-1, 1:-1] * dt / dx * (cn[1:-1, 1:-1] - cn[1:-1, 0:-2]))) + \
                        (v[1:-1, 1:-1] * dt / dy *
                         (cn[1:-1, 1:-1] - cn[0:-2, 1:-1]))

                    c = channel_shape(field=c,
                                      w=w, h=h, d=d, s=s,
                                      bottom_left=bottom_left, top_left=top_left,
                                      bottom_right=bottom_right, top_right=top_right,
                                      anode_value=ca, cathode_value=0,
                                      border_value=0,
                                      nx=nx,
                                      ny=ny, on_bot=on_bot)

                c_arr.append(c.clone().cpu())
                ch_it.append(a)

                t = 0
                g = 1
                if h == 0:
                    h = 1
                elif h < 0:
                    t = 1
                    g = 0

                # bottom cathode
                bc = c[ny-top_left[0]-g*h-1:ny-top_left[0]-t*h,
                       top_left[1]-s-w-d-w-1:top_left[1]-s-w-d+1]
                # top cathode
                tc = c[bottom_right[0]+t*h:bottom_right[0]+g*h+1,
                       bottom_right[1]+s+w+d-1:bottom_right[1]+s+w+d+w+1]

                integr = tc.sum() - bc.sum()
                resp.append(integr.cpu())
                gc.collect()

            gc.collect()
            np.save(save_folder + '/resp', resp)
            fig = plt.figure(figsize=(50, 20), dpi=200)
            plt.plot(resp)
            plt.savefig(save_folder + '/conc.png')
            print('Consentration plot was saved')

            p_gif = channel_shape(field=field,
                                  w=w, h=h, d=d, s=s,
                                  bottom_left=bottom_left, top_left=top_left,
                                  bottom_right=bottom_right, top_right=top_right,
                                  anode_value=anode_value+10000,
                                  cathode_value=cathode_value+10000,
                                  border_value=border_value+10000, nx=nx,
                                  ny=ny, on_bot=on_bot).to('cpu').numpy()

            name = 'concentrarion'
            np.save(save_folder + '/c_arr', c_arr)

            print('Concentration data was saved')
            gc.collect()
            imgs = np.load(save_folder + '/c_arr.npy', allow_pickle=True)

            np.save('imgs', imgs)
            imgs = np.array([i.numpy() for i in imgs])
            vmin = 0
            vmax = 0 if len(imgs[imgs > 0]) == 0 else np.mean(imgs[imgs > 0])
            
            plot_gif(imgs, p_gif, name, vmin, vmax,
                     ch_it, save_folder, is_conc=True)
            
            gc.collect()
            print(f'Gif {name} was saved')
            
            c = c.to('cpu').numpy()
            c = c[::-1, :]

            plot_conc(c, ca,
                      w, h, d, s,
                      bottom_left, top_left,
                      bottom_right, top_right,
                      nx, ny, save_folder)

        shutil.copy('params_to_calc.py', save_folder)
        shutil.copy('params_channel_shape.py', save_folder)
        shutil.copy('params_conc.py', save_folder)
        write_params(save_folder,
                     N_it=N_it, stop_it=stop_it, l1norm_target=l1norm_target,
                     sine=sine, freq=freq)
        print('Params was saved')
        print(f'Directory name: {save_folder}')

        gc.collect()
    except:
        shutil.rmtree(save_folder)
        print(f'\nNothing to save! Folder {save_folder} was removed')


if __name__ == '__main__':

    import argparse

    N_it = 5000
    stop_it = 3000
    l1norm_target = 1e-7

    parser = argparse.ArgumentParser(
        description='Create params for calculation')
    parser.add_argument('--N_it', default=N_it, type=int,
                        help=f'Number of diffusion steps, default: {N_it}')
    parser.add_argument('--stop_it', default=stop_it, type=int,
                        help=f'Max number of "laplace2d_parralel_C0" iterations, default: {stop_it}')
    parser.add_argument('--l1norm_target', default=l1norm_target, type=float,
                        help=f'Target L1 norm value, default: {l1norm_target}')
    parser.add_argument('--sine', action='store_true',
                        help='To calc for navier data')
    parser.add_argument('--on_bot', action='store_true',
                        help='To allocate all electrodes on bottom line')
    args = parser.parse_args()

    main(N_it=args.N_it, stop_it=args.stop_it,
         l1norm_target=args.l1norm_target,
         sine=args.sine, on_bot=args.on_bot)

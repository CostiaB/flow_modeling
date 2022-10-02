from script.plot_gif import plot_gif
import numpy as np
import gc

names = ['p_init_0', 'u_init_0', 'v_init_0']
files = ['p_gr.npy', 'u_gr.npy', 'v_gr.npy']
p_gif = np.load('channel_for_gif.npy')
ch_it = np.load('ch_it.npy').tolist()

for file, name in zip(files, names):
    gc.collect()
    imgs = np.load(file)
    vmin, vmax = np.mean(imgs[imgs<0]), np.mean(imgs[imgs>0])
    mean = np.mean(np.abs([vmin, vmax]))
    #vmin, vmax = -mean, mean
    vmin, vmax = -0.02, 0.02
    plot_gif(imgs, p_gif, name, vmin, vmax, ch_it)
    gc.collect()
    print(name, 'saved')

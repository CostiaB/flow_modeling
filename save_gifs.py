from script.plot_gif import plot_gif
import numpy as np
import gc

path = './results/navier_2022_11_10_20_15/'
names = ['p_init_0', 'u_init_0', 'v_init_0']
files = ['p_gr.npy', 'u_gr.npy', 'v_gr.npy']
p_gif = np.load(path+'channel_for_gif.npy')
ch_it = np.load(path+'ch_it.npy').tolist()

for file, name in zip(files, names):
    gc.collect()
    imgs = np.load(path+file)
    #vmin, vmax = np.mean(imgs[imgs<0]), np.mean(imgs[imgs>0])
    vmin, vmax = np.mean(imgs[:3][imgs[:3]<0]), np.mean(imgs[:3][imgs[:3]>0])
    mean = np.mean(np.abs([vmin, vmax]))
    vmin, vmax = -mean, mean
    #vmin, vmax = -0.02, 0.02
    plot_gif(imgs, p_gif, 'home/lapa/' + name, vmin, vmax, ch_it)
    gc.collect()
    print(name, 'saved')

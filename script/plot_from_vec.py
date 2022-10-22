import numpy as np
from plot_gif import plot_gif
import matplotlib.pyplot as plt
import matplotlib

path = '/home/lapa/Downloads/2022_10_15_21_12-20221020T054807Z-001/2022_10_15_21_12/'


p = np.load(path + 'p_gr.npy')
u = np.load(path + 'u_gr.npy')
v = np.load(path + 'v_gr.npy')
ch_it = np.load(path + 'ch_it.npy').tolist()
channel = np.load(path + 'channel_for_gif.npy')

for imgs, file_name in zip([p, u, v], ['p', 'u', 'v']):
    vmin = 0 if len(imgs[imgs<0]) == 0 else np.mean(imgs[imgs<0])
    vmax = 0 if len(imgs[imgs>0]) == 0 else np.mean(imgs[imgs>0])
    if 0 not in [vmin, vmax]:    
        mean = np.mean(np.abs([vmin, vmax]))
        vmin, vmax = -mean, mean
        if vmax > 0.2:
            vmin, vmax = -0.02, 0.02
    vmin /= 100
    vmax /= 100
    plot_gif(imgs, channel, file_name, vmin, vmax, ch_it, path, 2)
    
'''    
vmin = 0 if len(p[p<0]) == 0 else np.mean([p<0])
vmax = 0 if len(p[p>0]) == 0 else np.mean(p[p>0])
if 0 not in [vmin, vmax]:    
    mean = np.mean(np.abs([vmin, vmax]))
    vmin, vmax = -mean, mean
    if vmax > 0.02:
        vmin, vmax = -0.02, 0.02  
vmin, vmax = 4e-15, -4e-15
fig, ax = plt.subplots(figsize=(16,10))
cmap = plt.get_cmap('RdBu')
plt.imshow(p[::-1,:], vmin=vmin, vmax=vmax, cmap=cmap)
plt.colorbar()
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","black"])
#plt.imshow(channel[::-1,:], alpha = 1, cmap=cmap)
'''
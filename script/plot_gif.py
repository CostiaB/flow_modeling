import gif 
import matplotlib
import matplotlib.pyplot as plt
import gc


def plot_gif(vector, p, file_name, vmin=-0.01, vmax=0.01, ch_it=None,
             save_folder=''):

    @gif.frame
    def helper_plot_1(frame, i, p):
        fig, ax = plt.subplots(figsize=(16,10))
        cmap = plt.get_cmap('RdBu')
        plt.imshow(frame[::-1,:], vmin=vmin, vmax=vmax, cmap=cmap)
        plt.colorbar()
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","black"])
        plt.imshow(p[::-1,:], alpha = 0.1, cmap=cmap)
        if ch_it == None:
            ax.set_title(i)
        else:
            ax.set_title(ch_it[i])
        

    frames = []
    for idx, frame in enumerate(vector):
        frames.append(helper_plot_1(frame, idx, p))
    
    duration = len(frames)*0.1
    file_name = save_folder + '/' + file_name + '.gif'
    
    gif.save(frames, file_name, 
            duration=duration)
    gc.collect()

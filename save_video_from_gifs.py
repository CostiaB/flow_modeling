import moviepy.editor as mp
import glob
import re

path = '/home/lapa/Downloads/ggg/'
vector_types = ['u_gr', 'v_gr', 'p_gr']

for vec_type in vector_types:
    clips = []
    name = path + vec_type + '[0-9]*.gif'
    print(name)
    gif_list = glob.glob(name)
    print(gif_list)
    
    regex = r'{}([0-9]*).*.pkl'.format(vec_type)
    gif_list.sort(key=lambda x: re.findall(regex, x))
    for gif in gif_list:
        clips.append(mp.VideoFileClip(gif))
    final = mp.concatenate_videoclips(clips)
    final.write_videofile(path + vec_type + '.mp4')
    print(f'Film {vec_type} was saved')

def main(gif_path):
    import moviepy.editor as mp
    clip = mp.VideoFileClip(gif_path)
    name = gif_path.split('.')[0]
    clip.write_videofile(name+'.mp4')
    print('Video was saved!')

    
if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Params for video saving')
    parser.add_argument('--gif_path',  type=str)
    args = parser.parse_args()
    main(gif_path=args.gif_path)
    
    
    

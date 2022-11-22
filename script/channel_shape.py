def channel_shape(field,
                  w, h, d, s,
                  bottom_left, top_left, 
                  bottom_right, top_right,
                  anode_value, cathode_value,
                  border_value,
                  nx, ny):
    '''
    field - target field
    w - electodes width,
    h - height, 
    d - distance between electrodes,
    s - distance from corner
     
    sizes of correspoding rectangulars in tuples
    (y,x)
    
    bottom_left 
    top_left 
    bottom_right
    top_right
    
    Values on electrodes
    
    andode_value
    cathode_value
    
    border_value - value on inner channel border
    ny  -  length of the channel
    
    '''
#(y,x)
    
    #rectanguars
    #left
    field[0:bottom_left[0], 0:bottom_left[1]] = border_value
    #bottom
    field[0:bottom_right[0], bottom_right[1]:] = border_value
    
    #top
    field[ny-top_left[0]:, 0:top_left[1]] = border_value
    #right
    field[ny-top_right[0]:, top_right[1]:] = border_value
 
    #h>0 or <0
    t = 0
    g = 1
    if h == 0:
        h = 1
    elif h < 0:
        t = 1
        g = 0
    
    return field
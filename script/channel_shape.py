def channel_shape(field,
                  w, h, d, s,
                  bottom_left, top_left, bottom_right,
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
    field[:, (nx-bottom_left[1]):] = border_value

    #right border
    field[bottom_right[0]:, (nx-bottom_left[1])] = border_value

    #right side of lower part
    field[0:bottom_right[0], bottom_right[1]] = border_value

    #bottom of channel
    field[bottom_right[0]-1, bottom_right[1]:-1] = border_value

    #left border
    field[0:(ny-top_left[0]), bottom_left[1]-1] = border_value
    
    #channel top
    field[ny-top_left[0], bottom_left[1]:top_left[1]] = border_value

    #channel top left side
    field[(ny-top_left[0]):, top_left[1]] = border_value

    
    #TO DO add each side to make border
 
    #h>0 or <0
    t = 0
    g = 1
    if h == 0:
        h = 1
    elif h < 0:
        t = 1
        g = 0
    #top anode
    field[ny-top_left[0]-g*h:ny-top_left[0]-t*h, top_left[1]-s-w:top_left[1]-s] = anode_value
    
    #top cathode
    field[ny-top_left[0]-g*h:ny-top_left[0]-t*h, top_left[1]-s-w-d-w:top_left[1]-s-w-d] = cathode_value

    #bottom anode
    field[bottom_right[0]+t*h:bottom_right[0]+g*h, bottom_right[1]+s:bottom_right[1]+s+w] = anode_value

    #bottom cathode
    field[bottom_right[0]+t*h:bottom_right[0]+g*h, bottom_right[1]+s+w+d:bottom_right[1]+s+w+d+w] = cathode_value
    
    
    return field
def channel_shape(field,
                  w, h, d, s,
                  bottom_left, top_left, 
                  bottom_right, top_right,
                  anode_value, cathode_value,
                  border_value,
                  nx, ny, on_bot=False):
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
    
    if on_bot:  
        #For electrodes when all allocated at bottom line    
        #top anode
        field[bottom_right[0]+t*h : bottom_right[0]+g*h,
              (top_left[1]-s-w) : (top_left[1]-s)] = anode_value
        
        #top cathode
        field[bottom_right[0]+t*h : bottom_right[0]+g*h,
              top_left[1]-s-w-d-w : top_left[1]-s-w-d] = cathode_value       
    else:    
        #top anode
        field[(ny-top_left[0]-g*h) : (ny-top_left[0]-t*h),
              (top_left[1]-s-w) : (top_left[1]-s)] = anode_value
        
        #top cathode
        field[ny-top_left[0]-g*h : ny-top_left[0]-t*h,
              top_left[1]-s-w-d-w : top_left[1]-s-w-d] = cathode_value
    

    #bottom anode
    field[bottom_right[0]+t*h : bottom_right[0]+g*h,
          bottom_right[1]+s : bottom_right[1]+s+w] = anode_value

    #bottom cathode
    field[bottom_right[0]+t*h : bottom_right[0]+g*h,
          bottom_right[1]+s+w+d : bottom_right[1]+s+w+d+w] = cathode_value
    
    
    
    
    
    return field

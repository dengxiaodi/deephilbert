import math

def hilbert_index_to_xy(z, r):    
    if r <= 0 :
        return -1, -1
    
    # positions of resolution r = 1 hilbert curve
    
    positions = (
        [0, 0], # quadrant 0
        [0, 1], # quadrant 1
        [1, 1], # quadrant 2
        [1, 0]  # quadrant 3
    )
    
    # the last two bits of hibert index z represents 
    # the quadrant of the index
    
    quadrant = z & 0x03
    z >>= 2
    rmin = 1
    w = 2
    x, y = positions[quadrant]
    while z > 0:
        quadrant = z & 0x03
        if quadrant == 0 : # left bottom 
                x, y = y, x
        elif quadrant == 1 : # left top
                x, y = x, y + w
        elif quadrant == 2 : # right top
                x, y = x + w, y + w
        elif quadrant == 3 : # right bottom
                x, y = w * 2 - y - 1, w - x - 1
        
        z >>= 2
        rmin += 1
        w <<= 1
        
    if((rmin ^ r) & 0x01) :
        x, y = y, x
    
    return x, y

def xy_to_hilbert_index(x, y, r) :
    # quadrants
    
    quadrants = [
        0, # (0, 0) 
        1, # (0, 1)
        3, # (1, 0)
        2  # (1, 1)
    ]
    
    # determin the minimum resolution rmin, 
    # the index of rmin will be the same with r
    
    rmin = int(math.floor(math.log(max(x, y), 2))) + 1
    if r < rmin :
        return -1
    
    # check the parities of rmin with r, 
    # exchange x, y if they are different
    
    if((rmin ^ r) & 0x01) :
        x, y = y, x
    
    w = 1 << (rmin - 1)
    z = 0
    while rmin > 0 :
        quadrant = quadrants[(int(x / w) << 1) | int(y / w)]
        if quadrant == 0 : # left bottom 
            x, y = y, x
        elif quadrant == 1 : # left top
            x, y = x, y - w
        elif quadrant == 2 : # right top
            x, y = x - w, y - w
        elif quadrant == 3 : # right bottom
            x, y = w - y - 1, w * 2 - x -1
        w >>= 1
        z <<= 2
        z |= quadrant
        rmin -= 1
    
    return z

# generate a hilbert point list of resolution r

def hilbert_indexes(r) :
    n = 4 ** r
    coordinates = []
    for i in range(n) :
        coordinates += [(hilbert_index_to_xy(i, r))]
        
    return coordinates

# generate a hilbert map list of length n vector with resolution r
# data points between adjacent hilbert points are linear interpolated

def stair_map(coordinates, index_start, index_end, x) :
    if index_start == index_end :
        return coordinates[index_start]
    
    x1, y1 = coordinates[index_start]
    x2, y2 = coordinates[index_end]
    
    offset = x - index_start
    if x1 == x2 :
        if y2 > y1 :
            return x1, y1 + offset
        else :
            return x1, y1 - offset
    elif y1 == y2 :
        if x2 > x1 :
            return x1 + offset, y1
        else :
            return x1 - offset, y1 
    return -1, -1
    
def hilbert_map(n, r) :
    index_count = 4 ** r
    
    # rescale n vector positions
    
    pos_list = [i * 1.0 * (index_count - 1) / (n - 1) for i in range(n)]
    
    # assign hilbert segment start indexe for each position
    
    pos_list = [(int(math.floor(x)), int(math.ceil(x)), x) for x in pos_list]
    
    # get hilbert index coordinates
    
    index_coords = hilbert_indexes(r)
    
    # get the mapped position coordinates
    
    pos_coords = [(stair_map(index_coords, pos[0], pos[1], pos[2])) for pos in pos_list]
    
    return pos_coords

def hilbert_plot(v, r) :
    data_len = len(v)
    data_coords = hilbert_map(data_len, r)
    data = [(coord[0][0], coord[0][1], coord[1]) for coord in zip(data_coords, v)]
    
    return data
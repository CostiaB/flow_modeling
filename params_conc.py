import numpy as np
from params_channel_shape import *

c0 = 1
ca = 0.03*1e-12
concentration = 0.03*1e-12
c = np.ones((ny, nx))*concentration
c_anode_value = 0.03*1e-12
c_cathode_value = 0

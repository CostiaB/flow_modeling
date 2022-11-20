import numpy as np


scale = 1
nx = int((400+300+400+300)/scale)
ny = int((250+150+250)/scale)

field = np.ones((ny, nx))
anode_value = 0
cathode_value = 0
border_value = 0
w = int(40/scale)
h = -1
d = int(150/scale)
s = int(50/scale)
bottom_left = (int(0/scale), int(1/scale))
top_left = (int(250/scale), int((nx-2)/scale))
bottom_right = (int(250/scale), int(0/scale))

v = np.zeros((ny, nx)) / 100
u = np.zeros((ny, nx)) / 100
p = np.zeros((ny, nx)) / 100000
b = np.zeros((ny, nx)) / 100000

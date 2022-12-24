import torch


scale = 1
nx = int((400+300+400+300)/scale)
ny = int((250+150+250)/scale)

field = torch.ones((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
anode_value = 0
cathode_value = 0
border_value = 0
w = int(40/scale)
h = int(20/scale)
d = int(150/scale)
s = int(50/scale)
bottom_left = (int(650/scale), int(10/scale))
top_left = (int(250/scale), int((nx-260)/scale))
bottom_right = (int(250/scale), int(260/scale))
top_right = (int(650/scale), int(-10/scale))

v = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
u = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
p = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
b = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)

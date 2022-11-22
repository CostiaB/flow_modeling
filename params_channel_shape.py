import torch


scale = 1
nx = int((400+300+400+300)/scale)
ny = int((250+150+250)/scale)

field = torch.ones((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
anode_value = 0
cathode_value = 0
border_value = 0
w = int(40/scale)
h = 20
d = int(150/scale)
s = int(50/scale)
bottom_left = (int(100/scale), int(400/scale))
bottom_right = (int(100/scale), int(300/scale)) #diff from point (0;ny)
top_left = (int(100/scale), int(400/scale)) #form point (0;0)
top_right = (int(100/scale), int(300/scale)) #diff from point (0;ny)




v = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
u = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
p = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)
b = torch.zeros((ny, nx), device='cuda', requires_grad=False, dtype=torch.double)

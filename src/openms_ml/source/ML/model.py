import torch
import torch.nn.functional as F


class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.layer_1 = torch.nn.Linear(100, 256)
        self.layer_2 = torch.nn.Linear(256, 1)

    def forward(self, x):
        x = F.relu(self.layer_1(x))
        x = F.relu(self.layer_2(x))
        return x


traced_net = torch.jit.trace(Net(), torch.randn(1, 100))
torch.jit.save(traced_net, 'net.pt')
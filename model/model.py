import torch
from torchvision.ops import MLP

class LinearModel(torch.nn.Module):

    def __init__(self, sizes):

        super(LinearModel, self).__init__()

        layers = []

        for idx, (input, output) in enumerate(sizes):
            layers.append(torch.nn.Linear(input, output))

            if idx < len(sizes) - 1:
                layers.append(torch.nn.ReLU())
                layers.append(torch.nn.Dropout(p=0.1))

        print(layers)

        self.layers = torch.nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)



class MLP(torch.nn.Module):

    def __init__(self, sizes):

        super(MLP, self).__init__()

        layers = []

        for idx, (input, output) in enumerate(sizes[:-1]):
            layers.append(MLP(input, output, dropout=0.1))

        #self.layers.append(nn.Linear(sizes[-1][0], sizes[-1][1]))
        self.layers = torch.nn.Sequential(*layers)


    def forward(self, x):
        return self.layers(x)




import torch

class LinearModel(torch.nn.Module):

    def __init__(self, sizes):

        super(LinearModel, self).__init__()

        layers = []

        for idx, (input, output) in enumerate(sizes):
            layers.append(torch.nn.Linear(input, output))

            if idx < len(sizes) - 1:
                layers.append(torch.nn.LeakyReLU())
                layers.append(torch.nn.Dropout(p=0.1))

        self.layers = torch.nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)





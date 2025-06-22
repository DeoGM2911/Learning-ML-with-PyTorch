import torch
from torch import nn
from torch.nn import functional as F
import numpy as np
from .functional import res_body, save_hyperparameters


class ResNetBlock(nn.Module):
    """The ResNet block"""
    def __init__(self, in_channels, channels, depth, groups, down_sample=False):
        super(ResNetBlock, self).__init__()
        self.net = nn.Sequential()
        save_hyperparameters(self)
        
        # Add the layers in
        prev_in_channels = in_channels
        for i in range(depth):
            self.net.add_module(f"block_body_{i}", res_body(self.channels, self.groups[i]))
            self.net.add_module(f"batch_norm_{i}", nn.LazyBatchNorm1d())
            self.net.add_module(f"activation_{i}", nn.ReLU())
        if channels[-1] != in_channels:
            self.net.add_module(f"conv1x1", nn.LazyConv1d(self.channels[-1], kernel_size=1))
    
    def forward(self, X):
        out = X
        for i in range(self.depth):
            out = self.net.get_submodule(f"block_body_{i}")(out)
            out = self.net.get_submodule(f"batch_norm_{i}")(out)
            out = self.net.get_submodule(f"activation_{i}")(out)
        
        iden = F.relu(X)
        if self.down_sample:
            out = F.avg_pool1d(out, 2)
            iden = F.avg_pool1d(iden, 2)
        
        # If the numbers of channels are different, do a 1 - 1D conv
        if out.size(1) != iden.size(1):
            iden = self.net.conv1x1(iden)
        return out + iden


class PolymerCNN(nn.Module):
    def __init__(self, out_size, channels, groups, depth, depths, hidden):
        super(PolymerCNN, self).__init__()
        save_hyperparameters(self)

        # the network
        self.stem = nn.Sequential()
        self.body = nn.Sequential()
        self.head = nn.Sequential()
        
        # the stem
        self.stem.append(nn.Sequential(
            nn.LazyConv1d(4, kernel_size=5, stride=2, padding=1),
            nn.LazyBatchNorm1d(),
            nn.ReLU())
        )
        
        # the body
        for i in range(depth):
            down_sample = False
            if i % 2 == 1:
                down_sample = True
            
            self.body.append(
                ResNetBlock(4, channels[3*i: 3*i + 3], depths[i], groups[i], down_sample),
            )
        
        # the head
        self.head.append(nn.Sequential(
            nn.Flatten(),
            nn.LazyLinear(hidden),
            nn.LazyBatchNorm1d(),
            nn.ReLU(),
            nn.LazyLinear(out_size)
        ))
        

    def forward(self, X):
        out = self.stem(X)
        out = self.body(out)
        out = self.head(out)
        return out


class MLP(nn.Module):
    def __init__(self, out_size, hiddens):
        super(MLP, self).__init__()
        save_hyperparameters(self)
        self.net = nn.Sequential()
        for hidden_size in hiddens:
            self.net.append(nn.Sequential(
                nn.LazyLinear(hidden_size),
                nn.Dropout(),
                nn.LazyBatchNorm1d(),
                nn.ReLU()
            ))
        self.final_layer = nn.LazyLinear(out_size)
    
    def forward(self, X):
        out = self.net(X)
        return out


class PolymerRNN(nn.Module):
    def __init__(self, input_size, out_size, num_layer, hidden_size, hiddens, device):
        super(PolymerRNN, self).__init__()
        save_hyperparameters(self)
        
        self.rnn = nn.LSTM(
            input_size=input_size,
            hidden_size=hidden_size,
            num_layers=num_layer,
            batch_first=True
        )
        
        self.fc = nn.Sequential()
        for hidden in hiddens:
            self.fc.append(nn.Sequential(
                nn.LazyLinear(hidden),
                nn.LazyBatchNorm1d(),
                nn.ReLU()
            ))
        
        self.final_fc = nn.LazyLinear(out_size)
    
    def forward(self, X):
        shape = (self.num_layer, X.size(0), self.hidden)
        h0, c0 = torch.zeros(shape).to(self.device), torch.zeros(shape).to(self.device)
        out, _ = self.rnn(X, (h0, c0))
        if len(self.hiddens) != 0:
            out = self.fc(out)
        out = self.final_fc(out)
        return out

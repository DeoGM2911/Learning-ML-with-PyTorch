import torch
from torch import nn
from model.functional import save_hyperparameters


class MADLoss(nn.Module):
    def __init__(self):
        super(MADLoss, self).__init__()        
    
    def mad_loss(self, y_hat, y):
        return torch.mean(torch.abs(y_hat - y))
    
    def __call__(self, y_hat, y):
        return self.mad_loss(y_hat, y)


class weightedMADLoss(nn.Module):
    def __init__(self, weights, threshold):
        super(weightedMADLoss, self).__init__()
        assert len(weights) == 2
        save_hyperparameters(self)
        
    def weighted_mad(self, y_hat, y):
        # Compute the weights for this batch
        _weights = torch.zeros(y.size(0))
        _weights[:] = self.weights[1]
        _weights.masked_fill_(y < self.threshold, self.weights[0])
        
        # Compute the loss
        loss = _weights * torch.abs(y_hat - y)
        return torch.mean(loss)
    
    def __call__(self, y_hat, y):
        return self.weighted_mad(y_hat, y, self.weights, self.threshold)
    
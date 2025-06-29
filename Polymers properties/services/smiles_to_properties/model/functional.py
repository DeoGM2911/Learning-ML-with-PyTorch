import torch
from torch import nn
import inspect


def save_hyperparameters(obj, ignore=[]):
    """Register all attributes when initialize the module"""
    frame = inspect.currentframe().f_back
    args, _, _, values = inspect.getargvalues(frame)
    for arg in args:
        if arg != "self" and arg not in ignore:
            setattr(obj, arg, values[arg])


def res_body(channels, groups):
    """
    Return the residual 1D convolution block's body inspired by the ResNet architecture
    """
    assert len(channels) == 3
    blk = nn.Sequential(
        nn.LazyConv1d(channels[0], kernel_size=1),
        nn.LazyBatchNorm1d(),
        nn.ReLU(),
        nn.LazyConv1d(channels[1], kernel_size=3, padding=1, groups=groups),
        nn.LazyBatchNorm1d(),
        nn.ReLU(),
        nn.LazyConv1d(channels[2], kernel_size=1),
        nn.LazyBatchNorm1d(),
        nn.ReLU()
    )
    return blk

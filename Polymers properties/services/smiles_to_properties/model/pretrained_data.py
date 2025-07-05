from ..model.data import PolymerDataset
import torch

TG_MAX = 472.25
RG_MAX = 34.672905605
RG_MIN = 9.7283551


def convert_tg(outs):
    return outs.cpu() * TG_MAX - 273.15

def convert_rg(outs):
    return outs.cpu() * (RG_MAX - RG_MIN) + RG_MIN
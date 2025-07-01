from ..model.data import PolymerDataset
import torch

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
DATA = PolymerDataset("./services/smiles_to_properties/data/train.csv", device)

TG_MAX = DATA.tg_max 
RG_MAX = DATA.rg_max 
RG_MIN = DATA.rg_min


def convert_tg(outs):
    return outs.cpu() * TG_MAX - 273.15

def convert_rg(outs):
    return outs.cpu() * (RG_MAX - RG_MIN) + RG_MIN
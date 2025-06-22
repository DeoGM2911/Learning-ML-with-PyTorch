import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle 
from utils.vectorize import EMBED_DIM, mol_vectorize
import torch
import numpy as np


class PolymerDataset():
    LOWER_PROP = {"tg": "Tg", "tc": "Tc", "ffv": "FFV", 
                    "density": "Density", "rg": "Rg"}
    
    def __init__(self, file_path, device):
        self.original_df = pd.read_csv(file_path)
        df = shuffle(self.original_df.copy(True))
        self.device = device
        
        # Split the data based on the properties. Enforce the property to be non-NULL
        df_tg = df[df["Tg"].notnull()][["SMILES", "Tg"]]
        df_ffv = df[df["FFV"].notnull()][["SMILES", "FFV"]]
        df_tc = df[df["Tc"].notnull()][["SMILES", "Tc"]]
        df_density = df[df["Density"].notnull()][["SMILES", "Density"]]
        df_rg = df[df["Rg"].notnull()][["SMILES", "Rg"]]
        self.train_test_split(["tg", "ffv", "tc", "density", "rg"], 
                            [df_tg, df_ffv, df_tc, df_density, df_rg])
        
        # Convert from degree celcius to kelvin and scale to be between 0 and 1
        self.tg_max = self.get_data("tg")["Tg"].max()
        self.get_data("tg")["Tg"] = (self.get_data("tg")["Tg"] + 273.15) / self.tg_max
        self.get_data("tg", False)["Tg"] = (self.get_data("tg", False)["Tg"] + 273.15) / self.tg_max

        # Rescale the radius of gyration
        self.rg_max = self.get_data("rg")["Rg"].max()
        self.rg_min = self.get_data("rg")["Rg"].min()
        self.get_data("rg")["Rg"] = (self.get_data("rg")["Rg"] - self.rg_min) / (self.rg_max - self.rg_min)
        self.get_data("rg", False)["Rg"] = (self.get_data("rg", False)["Rg"] - self.rg_min) / (self.rg_max - self.rg_min)
        
        # Vectorize the SMILES
        self.vectorize(["tg", "tc", "ffv", "density", "rg"])
        
    def get_data(self, property: str, train=True):
        assert property in ["tg", "tc", "ffv", "density", "rg"]
        idx = 1
        if train:
            idx = 0
        
        return getattr(self, f"{property}_dataset")[idx]
    
    def train_test_split(self, names, dfs, test_size=0.2):
        for name, df in zip(names, dfs):
            setattr(self, f"{name}_dataset", train_test_split(df, test_size=test_size))
    
    def create_vocab(self):
        self.vocab = {}
        self.rev_vocab = {}
        for i, s in enumerate(self.original_df["SMILES"]):
            self.vocab[s] = i
            self.rev_vocab[i] = s
    
    def get_targets(self, property, train=True):
        assert property in ["tg", "tc", "ffv", "density", "rg"]
        idx = 1
        if train:
            idx = 0
        
        y = getattr(self, f"{property}_dataset")[idx][self.LOWER_PROP[property]]
        
        # Converrt to Torch tensors
        return torch.from_numpy(y.to_numpy().astype(np.float32)).view(-1, 1).to(self.device)
    
    def get_features(self, property, train=True, cnn_model=False):
        assert property in ["tg", "tc", "ffv", "density", "rg"]
        idx = 1
        if train:
            idx = 0
        
        smiles = getattr(self, f"{property}_embed_{idx}")
        
        X = torch.from_numpy(smiles.astype(np.float32))
        if cnn_model:
            X = X.view(len(X), 1, -1)
        
        return X.to(self.device)
    
    def vectorize(self, properties):
        for property in properties:
            df = getattr(self, f'{property}_dataset')
            
            for i in range(2):
                setattr(self, f"{property}_embed_{i}", mol_vectorize(
                    getattr(self, f'{property}_dataset')[i]["SMILES"], EMBED_DIM
                ))
    
    def convert_tg(self, outs):
        return outs.cpu() * self.tg_max - 273.15

    def convert_rg(self, outs):
        return outs.cpu() * (self.rg_max - self.rg_min) + self.rg_min

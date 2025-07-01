import sys
import os
from time import perf_counter

# Append the smiles_to_properties dir to sys.path
# sys.path.append(f"{os.getcwd()}/services/smiles_to_properties/")

from ..utils.pretrained import load_pretrained, predict
from ..utils.vectorize import mol_vectorize, EMBED_DIM
import pandas as pd
import torch
import numpy as np
from pprint import pprint


def test_models(test_file):
    # Read in the test files
    test_data = pd.read_csv(test_file)
    
    # Embedding
    X = mol_vectorize(test_data["SMILES"], EMBED_DIM)
    
    # Convert to Torch tensor
    X = torch.from_numpy(X.astype(np.float32))
    
    # Load in the models
    models = load_pretrained()
    
    # get predictions
    start = perf_counter()
    y = predict(models, X)
    end = perf_counter()
    return y, end - start


if __name__ == "__main__":
    y, time = test_models(f"{os.getcwd()}/services/smiles_to_properties/data/test_data.csv")
    print(f"Time needed: {time:.2f}")
    pprint(y, indent=2)
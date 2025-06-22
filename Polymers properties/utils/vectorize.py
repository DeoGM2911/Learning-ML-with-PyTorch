from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.DataStructs import ConvertToNumpyArray
import pandas as pd
import numpy as np

# Default word model for vectorizing
EMBED_DIM = 2048


def mol_vectorize(documents:pd.DataFrame, embed_dim: int):
    """
    Vectorize a SMILES documenet of molecues. Note that the provided word model 
    **MUST** have the embedding dimension equal to the passed `embed_dim` param.
    
    # Parameters:
    
    - word_model: the word2vec model to convert sentences representing molecues into 
    vectors
    - documents: the list of SMILES string representation of the molecues. Expected to be
    cleaned pandas data frame
    - embed_dim: the embedding dimension of the word_model
    
    # Return:
    
    - a numpy array of embedded molecues
    """
    # Vectorize each molecue
    embeddings = np.zeros((len(documents), embed_dim))
    for i, molecue in enumerate(documents):
        mol = Chem.MolFromSmiles(molecue)

        # Convert to "sentence" of atom environments
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=embed_dim)

        # Generate embedding
        fp_array = np.zeros((embed_dim,), dtype=int)
        ConvertToNumpyArray(fp, fp_array)
        embeddings[i] = fp_array
    return embeddings

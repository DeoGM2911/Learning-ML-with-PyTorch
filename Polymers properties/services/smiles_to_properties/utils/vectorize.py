from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import ConvertToNumpyArray
import pandas as pd
import numpy as np

# Default word model for vectorizing
EMBED_DIM = 2048


def mol_vectorize_text(smiles, embed_dim):
    """
    given a SMILES string, return the embedding for the it.
    """
    embedding = np.zeros((1, embed_dim))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=embed_dim)
    mol = Chem.MolFromSmiles(smiles)

    # Convert to "sentence" of atom environments
    fp = gen.GetFingerprint(mol)

    # Generate embedding
    fp_array = np.zeros((embed_dim,), dtype=int)
    ConvertToNumpyArray(fp, fp_array)
    embedding[0] = fp_array
    return embedding


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
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=embed_dim)
    for i, molecue in enumerate(documents):
        mol = Chem.MolFromSmiles(molecue)

        # Convert to "sentence" of atom environments
        fp = gen.GetFingerprint(mol)

        # Generate embedding
        fp_array = np.zeros((embed_dim,), dtype=int)
        ConvertToNumpyArray(fp, fp_array)
        embeddings[i] = fp_array
    return embeddings

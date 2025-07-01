from rdkit import Chem
from .validate_txt import clean_smiles


def inchi_to_smiles(inchi):
    # Convert to mol
    if not (mol := Chem.MolFromInchi(inchi)):
        return "<NOT FOUND>"
    
    # Convert to smiles
    smiles = Chem.MolToSmiles(mol)
    return clean_smiles(smiles)
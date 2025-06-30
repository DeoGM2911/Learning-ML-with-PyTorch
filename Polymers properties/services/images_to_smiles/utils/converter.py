from rdkit import Chem


def inchi_to_smiles(inchi):
    # Convert to mol
    if not (mol := Chem.MolFromInchi(inchi)):
        return "<NOT FOUND>"
    
    # Convert to smiles
    smiles = Chem.MolToSmiles(mol)
    return smiles
from rdkit import Chem


def inchi_to_smiles(inchi):
    # Convert to mol
    mol = Chem.MolFromInchi()
    
    # Convert to smiles
    smiles = Chem.MolToSmiles()
    return smiles
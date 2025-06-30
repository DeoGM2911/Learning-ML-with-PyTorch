from rdkit import Chem
from rdkit import RDLogger

# Suppress all RDKit warnings and info messages
RDLogger.DisableLog('rdApp.*')


def is_inchi(s):
    try:
        Chem.MolFromInchi(s, treatWarningAsError=True)
        return True
    except Exception:
        return False


def is_smiles(s):
    if Chem.MolFromSmiles(s) is None:
        return False
    return True


def is_smiles_polymer(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False
    # Check for wildcard atoms (possible polymer ends)
    return any(atom.GetSymbol() == '*' for atom in mol.GetAtoms())
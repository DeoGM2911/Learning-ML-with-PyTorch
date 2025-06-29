from bigsmiles import BigSMILES
from rdkit import Chem


def is_inchi(s):
    try:
        Chem.MolFromInchi(s)
        return True
    except Exception:
        return False


def is_smiles(s):
    try:
        Chem.MolFromSmiles(s)
        return True
    except Exception:
        return False


def is_smiles_polymer(smiles):
    try:
        parsed = BigSMILES(smiles)
        return True
    except Exception:
        return False
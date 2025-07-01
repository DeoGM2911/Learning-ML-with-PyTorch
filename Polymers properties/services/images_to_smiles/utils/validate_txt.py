from rdkit import Chem
from rdkit import RDLogger
import re

# Suppress all RDKit warnings and info messages
RDLogger.DisableLog('rdApp.*')


def clean_smiles(smiles):
    # Replace [X] or X with wildcard [*]
    # add * to the start if not already there
    if not smiles.startswith('*'):
        smiles = '*' + smiles
    return re.sub(r'\[?X\]?', '[*]', smiles)


def is_inchi(s):
    try:
        Chem.MolFromInchi(s, treatWarningAsError=True)
        return True
    except Exception:
        return False


def is_smiles(s):
    s = clean_smiles(s)
    if Chem.MolFromSmiles(s) is None:
        return False
    return True


def is_smiles_polymer(smiles):
    smiles = clean_smiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False
    # Check for wildcard atoms (possible polymer ends)
    return any(atom.GetSymbol() == '*' for atom in mol.GetAtoms())
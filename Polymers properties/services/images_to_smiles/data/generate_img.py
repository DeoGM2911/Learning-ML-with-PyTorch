from rdkit import Chem
from rdkit.Chem import Draw
from pathlib import Path

l = ["[*]CC[*]",
            "[*]C(C)C[*]",
            "[*]C(Cl)C[*]",
            "[*]OC(=O)c1ccc(C(=O)O[*])cc1",
            "[*]C(C)(C(=O)OC)C[*]",
            "[*]OC(=O)CCCCCC(=O)NCCCCC[*]",
            "[*]C(F)(F)C(F)(F)[*]",
            "[*]OC(=O)c1ccc(C(=O)O[*])cc1",
            "[*]c1ccc(N)c(/c1)[*]",
            "CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(c(=O)n2C)C", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
            "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"]

current = Path(__file__).resolve()
img_path = current.parent.parent / "images"

for i, chem in enumerate(l):
    # Create a molecule from SMILES
    mol = Chem.MolFromSmiles(chem)

    # Generate image
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(img_path / f"/{i}.png")
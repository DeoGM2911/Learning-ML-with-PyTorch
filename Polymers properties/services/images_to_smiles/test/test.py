from pathlib import Path
import sys
import os
from pprint import pprint

# Append the smiles_to_properties dir to sys.path
# sys.path.append(f"{os.getcwd()}/services/images_to_smiles/")

current = Path(__file__).resolve()
utils = current.parent.parent
data = current.parent.parent / "data"
img = current.parent.parent / "images"

import pandas as pd
from ..utils.converter import inchi_to_smiles
from ..utils.validate_txt import *
from ..models.decimer import predict_smiles, is_chemical_image


def test_converter():
    inchi = pd.read_csv(data / "inchi.csv", delimiter="!")
    its = []
    # Convert it to smiles
    for inchi_str in inchi["InChI"]:
        its.append(inchi_to_smiles(inchi_str))
    # Convert to a DF
    inchi["SMILES"] = pd.Series(its, name="smiles")
    inchi.to_csv("inchi_smiles.csv", index=False)


def test_validate_text():
    smiles_samples = ["[*]CC[*]",
            "[*]C(C)C[*]",
            "[*]C(Cl)C[*]",
            "[*]OC(=O)c1ccc(C(=O)O[*])cc1",
            "[*]C(C)(C(=O)OC)C[*]",
            "[*]OC(=O)CCCCCC(=O)NCCCCC[*]",
            "[*]C(F)(F)C(F)(F)[*]",
            "[*]OC(=O)c1ccc(C(=O)O[*])cc1",
            "[*]c1ccc(N)c(/c1)[*]"]
    fake_samples = ["dkljasf", "hello world!", "I'm noob", "jklahsdkjbc", "carbon",
                    "methanol", "este"]
    non_polymer_smiles = ["CC(=O)Oc1ccccc1C(=O)O", "Cn1cnc2c1c(=O)n(c(=O)n2C)C", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
                          "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"]

    for sample in smiles_samples:
        assert is_smiles(sample) == True
        assert is_smiles_polymer(sample) == True
    for sample in fake_samples:
        # assert is_inchi(sample) == False
        assert is_smiles(sample) == False
        assert is_smiles_polymer(sample) == False
    for sample in non_polymer_smiles:
        assert is_smiles(sample) == True
        assert is_smiles_polymer(sample) == False


def test_image_to_smiles():
    img_paths = [img / f"{i}.png" for i in range(13)]
    
    with open(data / "img_to_smiles.txt", "w") as out_file:
        for mol_img in img_paths:
            out_file.write(f"{predict_smiles(mol_img)}\n")


def test_is_chemical():
    result = {}
    for file in img.iterdir():
        result[file] = is_chemical_image(img / file)
    
    pprint(result, indent=4)


if __name__ == "__main__":
    test_converter()
    test_validate_text()
    test_image_to_smiles()
    test_is_chemical()
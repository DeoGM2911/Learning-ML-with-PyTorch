import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from smiles_to_properties.utils.pretrained import load_pretrained, predict
from smiles_to_properties.utils.vectorize import mol_vectorize, EMBED_DIM
from images_to_smiles.models.decimer import predict_smiles, is_chemical_image
from images_to_smiles.utils.validate_txt import is_smiles_polymer
from pathlib import Path


if __name__ == "__main__":
    # Load in the model
    models = load_pretrained()
    
    # Convert the image into smiles if possible, else ignore
    current = Path(__file__).resolve()
    dest = Path(__file__).resolve().parent / "test_outputs"
    imgs = current.parent.parent / "images_to_smiles" / "images"
    for img in imgs.iterdir():
        if is_chemical_image(img):
            smiles = predict_smiles(img)
            if is_smiles_polymer(smiles):
                # Generate prediction
                X = mol_vectorize([smiles], EMBED_DIM)
                with open(str(dest / f"{str(img).split('.')[0]}.txt"), "w") as f:
                    f.write(str(predict(models, X)))
            else:
                print("Not a polymer!")
        else:
            print("Not a chemical!")

from ..smiles_to_properties.utils.pretrained import load_pretrained, predict
from ..smiles_to_properties.utils.vectorize import mol_vectorize, EMBED_DIM
from ..images_to_smiles.models.decimer import predict_smiles, is_chemical_image
from ..images_to_smiles.utils.validate_txt import is_smiles_polymer
from pathlib import Path


if __name__ == "__main__":
    # Load in the model
    models = load_pretrained()
    
    # Convert the image into smiles if possible, else ignore
    imgs = Path(__file__).resolve().parent.parent / "images_to_smiles" / "images"
    dest = Path(__file__).resolve().parent / "test_outputs" 
    for img in imgs.iterdir():
        if is_chemical_image(img)[0]:
            smiles = predict_smiles(img)
            if is_smiles_polymer(smiles):
                # Generate prediction
                X = mol_vectorize([smiles], EMBED_DIM)
                with open(f"{str(dest)}/{img.stem}.txt", "w") as f:
                    f.write(str(predict(models, X)))
            else:
                print("Not a polymer!")
        else:
            print("Not a chemical!")

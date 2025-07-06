from DECIMER import predict_SMILES
from decimer_image_classifier import DecimerImageClassifier
from PIL import Image
from ..utils.validate_txt import clean_smiles
import numpy as np

# Threshold for detecting chemical image
THRESHOLD = 0.005


def predict_smiles(image_path):
    """
    Using the DECIMER transformer, generate the SMILES representation of the 
    chemical in the image.
    """
    return clean_smiles(predict_SMILES(np.array(image_path) if isinstance(image_path, Image.Image) else image_path))


def is_chemical_image(img_file):
    """
    Using the DECIMER classifier to check if an image is a chemical molecue.
    Return the prediction (true/false) and the confidence.
    """
    clf = DecimerImageClassifier()
    img = Image.open(img_file) if isinstance(img_file, str) else img_file
    is_struct = clf.is_chemical_structure(img, THRESHOLD)
    score = clf.get_classifier_score(img)
    return is_struct, score

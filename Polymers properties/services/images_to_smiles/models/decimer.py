from DECIMER import predict_SMILES
from decimer_image_classifier import DecimerImageClassifier


def predict_smiles(image_path):
    """
    Using the DECIMER transformer, generate the SMILES representation of the 
    chemical in the image.
    """
    return predict_SMILES(image_path, True)


def is_chemical_image(img_file):
    """
    Using the DECIMER classifier to check if an image is a chemical molecue.
    Return the prediction (true/false) and the confidence.
    """
    clf = DecimerImageClassifier()
    is_struct = clf.is_chemical_structure(img_file)
    score = clf.get_classifier_score(img_file)
    return is_struct, score
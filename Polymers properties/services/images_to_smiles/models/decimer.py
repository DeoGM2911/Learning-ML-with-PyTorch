from decimer_segmentation import segment_image
from transformers import TrOCRProcessor, VisionEncoderDecoderModel
from decimer_image_classifier import DecimerImageClassifier
from PIL import Image
import torch


# Load model and processor (can take time on first run)
print("Loading DECIMER model...")
model_name = "Kohulan/DECIMER-Image_Transformer"
processor = TrOCRProcessor.from_pretrained(model_name)
model = VisionEncoderDecoderModel.from_pretrained(model_name)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
model.eval()
print("Model loaded.")
clf = DecimerImageClassifier()


def predict_smiles(image_path):
    """
    Using the DECIMER transformer, generate the SMILES representation of the 
    chemical in the image.
    """
    # Load image
    image = Image.open(image_path).convert("RGB")

    # Preprocess
    pixel_values = processor(images=image, return_tensors="pt").pixel_values
    pixel_values = pixel_values.to(device)

    # Predict
    with torch.no_grad():
        generated_ids = model.generate(pixel_values)
        smiles = processor.batch_decode(generated_ids, skip_special_tokens=True)[0]
    return smiles.strip()


def is_chemical_image(img_file):
    """
    Using the DECIMER classifier to check if an image is a chemical molecue.
    Return the prediction (true/false) and the confidence.
    """
    is_struct = clf.is_chemical_structure(img_file)
    score = clf.get_classifier_score(img_file)
    return is_struct, score
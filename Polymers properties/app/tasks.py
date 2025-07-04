# tasks.py
from app import celery
from PIL import Image
from io import BytesIO


@celery.task
def add(x, y):
    return x + y


@celery.task
def process_prediction(text_input, file_bytes=None):
    """Do your heavy inference here."""
    pil_img = None
    if file_bytes:
        # convert the raw bytes back to PIL
        pil_img = Image.open(BytesIO(file_bytes)).convert('RGB')

        # now you can run your preprocess pipeline on `pil_img`, e.g.:
        # pil_img = pil_img.resize((224, 224))
        # arr = np.array(pil_img) / 255.0
        # …
    else:
        # Expect the text input to be a string
        ...
    
    # TODO: your model inference using text_input and/or pil_img
    reply = f"Processed: {text_input or 'image data'}"
    return reply
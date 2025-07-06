# tasks.py
from app import celery
from PIL import Image
from io import BytesIO
from services.smiles_to_properties.utils.pretrained import MODELS, predict
from services.smiles_to_properties.utils.vectorize import mol_vectorize_text, EMBED_DIM
from services.images_to_smiles.models.decimer import predict_smiles, is_chemical_image
from services.images_to_smiles.utils.validate_txt import clean_smiles, is_smiles, is_smiles_polymer


def generate_reply(props, error):
    """Generate a reply based on the predicted properties."""
    if props:
        return f"""<!-- Prediction Results Card -->
<section class="prediction-results">
  <h2>Predicted Properties</h2>
  <table class="properties-table">
    <thead>
      <tr>
        <th>Property</th>
        <th>Value</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>Glass Transition Temperature</td>
        <td>{props['tg'].item():.2f} °C</td>
      </tr>
      <tr>
        <td>Fractional Free Volume (FFV)</td>
        <td>{props['ffv'].item():.4f}</td>
      </tr>
      <tr>
        <td>Thermal Conductivity</td>
        <td>{props['tc'].item():.2f} W/m·K</td>
      </tr>
      <tr>
        <td>Density</td>
        <td>{props['density'].item():.3f} g/cm<sup>3</sup></td>
      </tr>
      <tr>
        <td>Radius of Gyration</td>
        <td>{props['rg'].item():.2f} Å</td>
      </tr>
    </tbody>
  </table>
  <p class="note">
    <em>Note:</em> Sometimes the model may fail to detect non-polymer structures.
  </p>
</section>
        """
    else:
        return error


@celery.task
def process_prediction(text_input, file_bytes=None):
    global MODELS
    """Do your heavy inference here."""
    pil_img = None
    if file_bytes:
        # convert the raw bytes to PIL Image
        try:
            pil_img = Image.open(BytesIO(file_bytes)).convert('RGB')
        except Exception as e:
            return generate_reply(None, 'Invalid image file.')
        else:
            # Check if the image is a chemical structure
            if not is_chemical_image(pil_img)[0]:
                return generate_reply(None, 'The provided image is not a chemical structure.')
    else:
        # validate the text input
        if not is_smiles(text_input):
            return generate_reply(None, 'Invalid SMILES input.')
        text_input = clean_smiles(text_input.strip())
        if not is_smiles_polymer(text_input):
            return generate_reply(None, 'Input is not a valid polymer.')

    smiles = text_input.upper() if text_input else None
    # If the user provide an image
    if pil_img:
        # Convert it to a SMILES before predicting
        smiles = predict_smiles(pil_img)
        if not is_smiles_polymer(smiles):
            return generate_reply(None, 'The image does not contain a valid polymer structure.')
    
    # Now predict the properties
    props = predict(MODELS, mol_vectorize_text(smiles, EMBED_DIM))
    
    return generate_reply(props, None)

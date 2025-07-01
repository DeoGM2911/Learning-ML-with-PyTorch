import sys
from pathlib import Path

# sys.path.append(str(Path(__file__).resolve().parent.parent))

import torch
import numpy as np
from ..model.models import PolymerCNN
from pathlib import Path
from ..model.pretrained_data import convert_rg, convert_tg

def load_pretrained():
    # CNN models
    tg_cnn = PolymerCNN(1, (16, 32, 32, 64, 128, 256), ((4, 8), (2, 4, 4, 8)), 2, (2, 4), 256)
    ffv_cnn = PolymerCNN(1, (32, 64, 64, 64, 128, 128), ((4, 8), (8, 16, 16, 32)), 2, (2, 4), 512)
    tc_cnn = PolymerCNN(1, (16, 32, 32, 64, 128, 256), ((4, 8), (2, 4, 4, 8)), 2, (2, 4), 256)
    density_cnn = PolymerCNN(1, (16, 32, 32, 64, 128, 256), ((4, 8), (2, 4, 4, 8)), 2, (2, 4), 256)
    rg_cnn = PolymerCNN(1, (16, 32, 32, 64, 128, 256), ((4, 8), (2, 4, 4, 8)), 2, (2, 4), 256)
    
    models = {'tg': tg_cnn, 'ffv': ffv_cnn, 'tc': tc_cnn, 'density': density_cnn, 'rg': rg_cnn}
    device = torch.device("cuda:0" if torch.cuda.is_available() else 'cpu')
    
    # Resolve the path
    current = Path(__file__).resolve()
    root = current.parent.parent
    path_pretrain = root / 'pretrained'
    
    # Load in the models
    for property in models:
        model = models[property]
        model.load_state_dict(torch.load(path_pretrain / f'{property}_cnn.pth', map_location=device, weights_only=True))
        model.eval()
    return models


def predict(models, X):
    if type(X) == torch.Tensor:
        X = X.numpy()
    X = torch.from_numpy(X.astype(np.float32)).view(len(X), 1, -1)
    tg_out = convert_tg(models["tg"](X).detach())
    ffv_out = models["ffv"](X).detach()
    tc_out = models["tc"](X).detach()
    density_out = models["density"](X).detach()
    rg_out = convert_rg(models["rg"](X).detach())
    return {'tg': tg_out.numpy(), 'ffv': ffv_out.numpy(), 'tc': tc_out.numpy(), 
            'density': density_out.numpy(), 'rg': rg_out.numpy()}

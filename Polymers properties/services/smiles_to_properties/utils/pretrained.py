import torch
from model.models import PolymerCNN


def load_pretrained():
    # CNN models
    tg_cnn = PolymerCNN(1, (32, 32, 64, 64, 128, 256), ((4, 8), (2, 4, 8, 8)), 2, (2, 4), 256)
    ffv_cnn = PolymerCNN(1, (32, 64, 64, 128, 128, 256), ((4, 8), (8, 16, 16, 32)), 2, (2, 4), 512)
    tc_cnn = PolymerCNN(1, (16, 32, 32, 64, 128, 256), ((4, 8), (2, 4, 8, 8)), 2, (2, 4), 256)
    density_cnn = PolymerCNN(1, (32, 32, 64, 64, 128, 256), ((4, 8), (2, 4, 8, 8)), 2, (2, 4), 256)
    rg_cnn = PolymerCNN(1, (16, 32, 32, 64, 128, 256), ((4, 8), (2, 4, 8, 8)), 2, (2, 4), 256)
    
    models = {'tg': tg_cnn, 'ffv': ffv_cnn, 'tc': tc_cnn, 'density': density_cnn, 'rg': rg_cnn}
    device = torch.device("cuda:0" if torch.cuda.is_available() else 'cpu')
    # Load in the models
    for property in models:
        model = models[property]
        model.load_state_dict(torch.load(f'{property}_cnn.pt', map_location=device))
    
    return models
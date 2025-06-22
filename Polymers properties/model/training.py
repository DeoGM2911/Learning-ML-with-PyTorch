import numpy as np
from torch import optim
from model.data import PolymerDataset


def configure(model, device, loss_func, lr):
    model.to(device)
    criterion = loss_func
    optimizer = optim.Adam(model.parameters(), lr=lr)
    return criterion, optimizer


def train(models, optimizers, criterions, dataset: PolymerDataset, cnn_model, num_epochs):
    properties = ["tg", "tc", "ffv", "density", "rg"]
    # To record the loss history
    train_history = {p: np.zeros(num_epochs) for p in properties}
    test_history = {p: np.zeros(num_epochs) for p in properties}
    
    for epoch in range(num_epochs):
        for property in properties:
            # Training
            model = models[property]
            model.train()
            optimizers[property].zero_grad()
            
            # Forward pass
            y_hat = model(dataset.get_features(property, True, cnn_model))
            train_loss = criterions[property](y_hat, dataset.get_targets(property))
            
            # Backward step
            train_loss.backward()
            optimizers[property].step()
            
            # Record the train loss
            train_history[property][epoch] = train_loss.item()

            # Validate
            model.eval()
            y_hat = model(dataset.get_features(property, False, cnn_model))
            test_loss = criterions[property](y_hat, dataset.get_targets(property, False))
            test_history[property][epoch] = test_loss.item()
    
    return train_history, test_history
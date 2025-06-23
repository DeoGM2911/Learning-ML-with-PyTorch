import numpy as np
from torch import optim
from model.data import PolymerDataset
import matplotlib.pyplot as plt
from utils.visualize import visualize_loss
from sklearn.utils import shuffle


def configure(model, device, loss_func, lr):
    model.to(device)
    criterion = loss_func()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    return criterion, optimizer


def train(models, optimizers, criterions, dataset: PolymerDataset, cnn_model, rnn_model, num_epochs, batch_size):
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
            
            X_train = dataset.get_features(property, True, cnn_model, rnn_model)
            y_train = dataset.get_targets(property)
            # Shuffle the data
            X_train, y_train = shuffle(X_train, y_train)
            N_train = int(np.ceil(len(y_train) // batch_size))
            train_loss_epoch = []
            for i in range(N_train):
                end = min((i+1) * batch_size, len(y_train))
                X_batch, y_batch = X_train[i * batch_size: end], y_train[i * batch_size: end]
                
                # Forward pass
                y_hat = model(X_batch)
                train_loss = criterions[property](y_hat, y_batch)
                
                # Backward step
                train_loss.backward()
                optimizers[property].step()
                train_loss_epoch.append(train_loss.item())
            
            # Record the train loss
            train_history[property][epoch] = np.mean(train_loss_epoch).item()

            # Validate
            model.eval()
            X_test = dataset.get_features(property, False, cnn_model, rnn_model)
            y_test = dataset.get_targets(property, False)
            N_test = int(np.ceil(len(y_test) // batch_size))
            test_loss_epoch = []
            for i in range(N_test):
                end = min((i+1) * batch_size, len(y_test))
                X_batch, y_batch = X_test[i * batch_size: end], y_test[i * batch_size: end]
                y_hat = model(X_batch)
                test_loss = criterions[property](y_hat, y_batch)
                test_loss_epoch.append(test_loss.item())
            
            test_history[property][epoch] = np.mean(test_loss_epoch).item()
        
        # For keeping track of time
        print(f"Epoch: {epoch}")
        
    return train_history, test_history

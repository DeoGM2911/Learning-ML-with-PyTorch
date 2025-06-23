import matplotlib.pyplot as plt
import pandas.plotting as pd_plt


def visualize_loss(train_hist, test_hist, name):
    fig, axes = plt.subplots(2, 3, figsize=(12, 6))
    axes = axes.flatten()
    
    for i, key in enumerate(train_hist.keys()):
        axes[i].set_title(f"Training with {key}")
        axes[i].plot(train_hist[key], label="train")
        axes[i].plot(test_hist[key], label="test")
        axes[i].set_xlabel("epoch")
        axes[i].set_ylabel("loss")
        axes[i].legend()
    
    axes[-1].axis('off')
    
    plt.tight_layout()
    plt.show()
    fig.savefig(f"../images/{name}.png")


def scatter_matrix(data):
    pd_plt.scatter_matrix(data, alpha=0.5)

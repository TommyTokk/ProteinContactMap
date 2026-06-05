#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np


def main():
    csv_file = sys.argv[1]
    fn = os.path.splitext(os.path.basename(csv_file))[0]
    try:
        df = np.loadtxt(csv_file, delimiter=',', dtype=np.float32)
    except:

        df = pd.read_csv(csv_file, header=None, dtype=np.float32).values

    if df.ndim != 2:
        raise ValueError(f"Expected a 2D contact map, got shape {df.shape}")

    # Contact maps are square; crop to the largest common square if needed.
    if df.shape[0] != df.shape[1]:
        min_dim = min(df.shape[0], df.shape[1])
        df = df[:min_dim, :min_dim]


    os.makedirs("plots", exist_ok=True)


    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    

    im = ax.imshow(df, cmap="YlOrRd", aspect='equal', interpolation='nearest')
    ax.set_box_aspect(1)
    

    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    ax.set_title("Contact map", fontweight="bold")
    

    plt.savefig(f"plots/{fn}.png", dpi=300,
                pil_kwargs={'optimize': True})
    plt.close(fig)


if __name__ == "__main__":
    main()
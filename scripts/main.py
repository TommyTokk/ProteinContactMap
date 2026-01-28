#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np


def main():
    csv_file = sys.argv[1]
    fn = csv_file.split(".")[0].split("/")[2]
    try:
        df = np.loadtxt(csv_file, delimiter=',', dtype=np.float32)
    except:

        df = pd.read_csv(csv_file, header=None, dtype=np.float32).values


    os.makedirs("plots", exist_ok=True)


    fig, ax = plt.subplots(figsize=(10, 10))
    

    im = ax.imshow(df, cmap="YlOrRd", aspect='auto', interpolation='nearest')
    

    plt.colorbar(im, ax=ax)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    ax.set_title("Contact map", fontweight="bold")
    

    plt.savefig(f"plots/{fn}.png", dpi=300, bbox_inches='tight', 
                pil_kwargs={'optimize': True})
    plt.close(fig)


if __name__ == "__main__":
    main()
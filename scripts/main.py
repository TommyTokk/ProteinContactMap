#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os


def main():
    csv_file = sys.argv[1]
    fn = csv_file.split(".")[0].split("/")[2]

    df = pd.read_csv(csv_file, header=None)

    # Create plots directory if it doesn't exist
    os.makedirs("plots", exist_ok=True)

    sns.heatmap(df, cmap="YlOrRd", cbar=True)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel("")
    plt.ylabel("")

    plt.title("Contact map", fontweight="bold")
    plt.savefig(f"plots/{fn}.png", dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()

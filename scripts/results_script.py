#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sys import argv



def analyse_model(model_path: str) -> dict:
    # Get the model name from the path
    model_path = Path(model_path)
    model_name = model_path.name

    model_results = {}

    file_list = sorted(list(model_path.glob("*.csv")))
    
    seq_res_df = None
    omp_res_df = None
    mpi_res_df = None

    for file in file_list:
        file_name = file.stem

        if "sequential" in file_name:
            seq_res_df = handle_sequential(file)
        elif "openmp" in file_name:
            omp_res_df = handle_omp(file)
        elif "mpi" in file_name:
            mpi_res_df = handle_mpi(file)
        else:
            raise ValueError(f"Unknown file name: {file_name}")


    
    
    return (seq_res_df, omp_res_df, mpi_res_df)
        

        
def handle_sequential(file: Path) -> pd.DataFrame:

    model_size = int(file.stem.split("_")[-1][:-3])
    print(f"{file.stem} | Sequential | Model size: {model_size}")
    df = pd.read_csv(file)

    # Throw away the first 5 runs
    df = df.iloc[5:]
    
    min_time = df["alg_time"].min()
    max_time = df["alg_time"].max()
    seq_avg_time = df["alg_time"].mean()
    seq_std_time = df["alg_time"].std()

    # Create a new DataFrame to store the results with index as key
    result_df = pd.DataFrame({
        "min_time": [min_time],
        "max_time": [max_time],
        "avg_time": [seq_avg_time],
        "std_time": [seq_std_time],
        "model_size": [model_size],
    })

    return result_df

def handle_omp(file: Path) -> pd.DataFrame:
    model_size = int(file.stem.split("_")[-1][:-3])
    print(f"{file.stem} | OpenMP | Model size: {model_size}")
    df = pd.read_csv(file)

    df = df.iloc[5:]
    
    # Group by the number of threads and calculate the average and std time for each group
    omp_results = df.groupby("num_threads")[["alg_time"]].agg(["mean", "std", "min", "max"]).reset_index()
    omp_results.columns = ["num_threads", "avg_time", "std_time", "min_time", "max_time"]
    omp_results["model_size"] = model_size

    return omp_results

def handle_mpi(file: Path) -> pd.DataFrame:
    model_size = int(file.stem.split("_")[-1][:-3])
    print(f"{file.stem} | MPI | Model size: {model_size}")
    df = pd.read_csv(file)

    df = df.iloc[5:]
    
    # Group by the number of processes and number of threads, then calculate the average and std time for each group
    mpi_results = df.groupby(["num_processors", "num_threads"])[["alg_time"]].agg(["mean", "std", "min", "max"]).reset_index()
    mpi_results.columns = ["num_processors", "num_threads", "avg_time", "std_time", "min_time", "max_time"]
    mpi_results["model_size"] = model_size
    return mpi_results

def omp_speedup(omp_res_df: pd.DataFrame, seq_res_df: pd.DataFrame) -> pd.DataFrame:
    # Group by the number of threads and calculate the speedup for each group using minimum times
    omp_speedup_df = omp_res_df.copy()
    omp_speedup_df["speedup"] = seq_res_df["min_time"].iloc[0] / omp_speedup_df["min_time"]
    return omp_speedup_df

def omp_efficiency(omp_res_df: pd.DataFrame, seq_res_df: pd.DataFrame) -> pd.DataFrame:
    # Group by the number of threads and calculate the efficiency for each group using minimum times
    omp_efficiency_df = omp_res_df.copy()
    omp_efficiency_df["efficiency"] = seq_res_df["min_time"].iloc[0] / (omp_efficiency_df["min_time"] * omp_efficiency_df["num_threads"])
    return omp_efficiency_df

def omp_cost(omp_res_df: pd.DataFrame, seq_res_df: pd.DataFrame) -> pd.DataFrame:
    # Group by the number of threads and calculate the cost for each group using minimum times
    omp_cost_df = omp_res_df.copy()
    omp_cost_df["cost"] = omp_cost_df["min_time"] * omp_cost_df["num_threads"]
    return omp_cost_df

def mpi_speedup(mpi_res_df: pd.DataFrame, seq_res_df: pd.DataFrame) -> pd.DataFrame:
    # Group by the number of processes and number of threads, then calculate the speedup for each group using minimum times
    mpi_speedup_df = mpi_res_df.copy()
    mpi_speedup_df["speedup"] = seq_res_df["min_time"].iloc[0] / mpi_speedup_df["min_time"]
    return mpi_speedup_df

def mpi_efficiency(mpi_res_df: pd.DataFrame, seq_res_df: pd.DataFrame) -> pd.DataFrame:
    # Group by the number of processes and number of threads, then calculate the efficiency for each group using minimum times
    mpi_efficiency_df = mpi_res_df.copy()
    mpi_efficiency_df["efficiency"] = (seq_res_df["min_time"].iloc[0] /((mpi_efficiency_df["min_time"]) * (mpi_efficiency_df["num_processors"] * mpi_efficiency_df["num_threads"])))
    return mpi_efficiency_df

def mpi_cost(mpi_res_df: pd.DataFrame, seq_res_df: pd.DataFrame) -> pd.DataFrame:
    # Group by the number of processes and number of threads, then calculate the cost for each group using minimum times
    mpi_cost_df = mpi_res_df.copy()
    mpi_cost_df["cost"] = mpi_cost_df["min_time"] * (mpi_cost_df["num_processors"] * mpi_cost_df["num_threads"])
    return mpi_cost_df

def plot_mpi_metric(mpi_all_df: pd.DataFrame, metric: str, metric_label: str, model_name: str):
    # Set seaborn style
    sns.set_theme(style="whitegrid")
    
    # Build a global color palette so all subplots share the same thread->color mapping
    all_threads = sorted(mpi_all_df["num_threads"].unique(), key=int)
    all_threads_str = [str(t) for t in all_threads]
    bright_colors = sns.color_palette("bright", n_colors=len(all_threads_str))
    thread_palette = dict(zip(all_threads_str, bright_colors))
    
    # Plot the MPI metric for each number of processors
    num_procs = sorted(mpi_all_df["num_processors"].unique())
    n_plots = len(num_procs)
    
    # Calculate number of columns for 2-row layout
    n_cols = (n_plots + 1) // 2
    n_rows = 2 if n_plots > 1 else 1
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8*n_cols, 6*n_rows))
    
    # Flatten axes array for easier iteration
    if n_plots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, n_proc in enumerate(num_procs):
        ax = axes[idx]
        proc_data = mpi_all_df[mpi_all_df["num_processors"] == n_proc].copy()
        
        if not proc_data.empty:
            proc_data["num_threads"] = proc_data["num_threads"].astype(str)
            
            sns.lineplot(
                data=proc_data,
                x="model_size",
                y=metric,
                hue="num_threads",
                hue_order=all_threads_str,
                palette=thread_palette,
                marker='o',
                markersize=4,
                linewidth=2,
                ax=ax
            )

        ax.set_xlabel("Model Size (residues)", fontsize=11)
        ax.set_ylabel(metric_label, fontsize=11)
        ax.set_title(f"{n_proc} Processors", fontsize=12, fontweight='bold')
        ax.set_xscale('log')
        ax.set_xticks([1e2, 1e3, 1e4])
        ax.set_xticklabels([r'$10^2$', r'$10^3$', r'$10^4$'])
        # Add darker vertical lines at placeholder (power-of-10) sizes
        for pw in [1e2, 1e3, 1e4]:
            ax.axvline(x=pw, color='dimgrey', linewidth=0.8, alpha=0.8, linestyle='-')
        # Add lighter vertical lines at each model size
        for ms in sorted(mpi_all_df["model_size"].unique()):
            ax.axvline(x=ms, color='lightgrey', linewidth=0.5, alpha=0.6, linestyle='--')
        # Remove per-subplot legend
        if ax.get_legend():
            ax.get_legend().remove()

    # Hide unused subplots
    for i in range(len(num_procs), len(axes)):
        axes[i].axis('off')

    # Single shared legend for the whole figure
    legend_handles = [plt.Line2D([0], [0], color=thread_palette[t], marker='o',
                      markersize=4, linewidth=2, label=t) for t in all_threads_str]
    fig.legend(handles=legend_handles, title="Threads", fontsize=9,
              loc='center right', bbox_to_anchor=(1.0, 0.5), borderaxespad=1)

    fig.suptitle(f"MPI {metric_label} - {model_name}", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.92, 1])
    
    # Save the plot
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    plt.savefig(results_dir / f"{model_name}_mpi_{metric}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved MPI {metric} plot to results/{model_name}_mpi_{metric}.png")

def plot_omp_metric(omp_all_df: pd.DataFrame, metric: str, metric_label: str, model_name: str):
    # Set seaborn style (if not set globally, or to ensure consistency)
    sns.set_theme(style="whitegrid")
    
    # Plot the OMP metric with model sizes on X-axis
    plt.figure(figsize=(10, 6))
    
    if not omp_all_df.empty:
        # Convert num_threads to string and sort properly
        plot_data = omp_all_df.copy()
        plot_data["num_threads"] = plot_data["num_threads"].astype(str)
        hue_order = sorted(plot_data["num_threads"].unique(), key=lambda x: int(x))

        sns.lineplot(
            data=plot_data,
            x="model_size",
            y=metric,
            hue="num_threads",
            hue_order=hue_order,
            palette="bright",
            marker='o',
            markersize=4,
            linewidth=2
        )
    
    plt.xlabel("Model Size (residues)", fontsize=11)
    plt.ylabel(metric_label, fontsize=11)
    plt.title(f"OpenMP {metric_label} - {model_name}", fontsize=14, fontweight='bold')
    plt.xscale('log')
    ax = plt.gca()
    ax.set_xticks([1e2, 1e3, 1e4])
    ax.set_xticklabels([r'$10^2$', r'$10^3$', r'$10^4$'])
    # Add darker vertical lines at placeholder (power-of-10) sizes
    for pw in [1e2, 1e3, 1e4]:
        ax.axvline(x=pw, color='dimgrey', linewidth=0.8, alpha=0.8, linestyle='-')
    # Add lighter vertical lines at each model size
    for ms in sorted(omp_all_df["model_size"].unique()):
        ax.axvline(x=ms, color='lightgrey', linewidth=0.5, alpha=0.6, linestyle='--')
    # Rebuild legend with clean title
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, title="Threads", fontsize=9, loc='best')

    # Save the plot
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    plt.savefig(results_dir / f"{model_name}_omp_{metric}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved OMP {metric} plot to results/{model_name}_omp_{metric}.png")


def plot_omp_avg_time(omp_all_df: pd.DataFrame, seq_all_df: pd.DataFrame, model_name: str):
    """Plot average execution time vs model size for OpenMP, with one line per thread count
    and the sequential baseline."""
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))

    # Plot sequential baseline (lighter style)
    seq_sorted = seq_all_df.sort_values("model_size")
    plt.plot(seq_sorted["model_size"], seq_sorted["avg_time"],
             marker='s', markersize=5, linewidth=1.5, color='grey',
             linestyle='--', alpha=0.7, label='Sequential', zorder=5)

    if not omp_all_df.empty:
        plot_data = omp_all_df.copy()
        plot_data["num_threads"] = plot_data["num_threads"].astype(str)
        hue_order = sorted(plot_data["num_threads"].unique(), key=lambda x: int(x))

        sns.lineplot(
            data=plot_data,
            x="model_size",
            y="avg_time",
            hue="num_threads",
            hue_order=hue_order,
            palette="bright",
            marker='o',
            markersize=4,
            linewidth=2
        )

    plt.xlabel("Model Size (residues)", fontsize=11)
    plt.ylabel("Average Time (s)", fontsize=11)
    plt.title(f"OpenMP Average Execution Time - {model_name}", fontsize=14, fontweight='bold')
    plt.xscale('log')
    plt.yscale('log')
    ax = plt.gca()
    ax.set_xticks([1e2, 1e3, 1e4])
    ax.set_xticklabels([r'$10^2$', r'$10^3$', r'$10^4$'])
    for pw in [1e2, 1e3, 1e4]:
        ax.axvline(x=pw, color='dimgrey', linewidth=0.8, alpha=0.8, linestyle='-')
    for ms in sorted(omp_all_df["model_size"].unique()):
        ax.axvline(x=ms, color='lightgrey', linewidth=0.5, alpha=0.6, linestyle='--')
    # Rebuild legend: combine sequential + thread lines, place outside plot
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, title="Threads", fontsize=9,
                  loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)

    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    plt.savefig(results_dir / f"{model_name}_omp_avg_time.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved OMP avg_time plot to results/{model_name}_omp_avg_time.png")


def plot_mpi_avg_time(mpi_all_df: pd.DataFrame, seq_all_df: pd.DataFrame, model_name: str):
    """Plot average execution time vs model size for MPI, with one subplot per processor count
    and one line per thread count, plus the sequential baseline."""
    sns.set_theme(style="whitegrid")

    # Build a global color palette so all subplots share the same thread->color mapping
    all_threads = sorted(mpi_all_df["num_threads"].unique(), key=int)
    all_threads_str = [str(t) for t in all_threads]
    bright_colors = sns.color_palette("bright", n_colors=len(all_threads_str))
    thread_palette = dict(zip(all_threads_str, bright_colors))

    num_procs = sorted(mpi_all_df["num_processors"].unique())
    n_plots = len(num_procs)

    n_cols = (n_plots + 1) // 2
    n_rows = 2 if n_plots > 1 else 1

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8*n_cols, 6*n_rows))

    if n_plots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    seq_sorted = seq_all_df.sort_values("model_size")

    for idx, n_proc in enumerate(num_procs):
        ax = axes[idx]

        # Plot sequential baseline (lighter style)
        ax.plot(seq_sorted["model_size"], seq_sorted["avg_time"],
                marker='s', markersize=5, linewidth=1.5, color='grey',
                linestyle='--', alpha=0.7, label='Sequential', zorder=5)

        proc_data = mpi_all_df[mpi_all_df["num_processors"] == n_proc].copy()

        if not proc_data.empty:
            proc_data["num_threads"] = proc_data["num_threads"].astype(str)

            sns.lineplot(
                data=proc_data,
                x="model_size",
                y="avg_time",
                hue="num_threads",
                hue_order=all_threads_str,
                palette=thread_palette,
                marker='o',
                markersize=4,
                linewidth=2,
                ax=ax
            )

        ax.set_xlabel("Model Size (residues)", fontsize=11)
        ax.set_ylabel("Average Time (s)", fontsize=11)
        ax.set_title(f"{n_proc} Processors", fontsize=12, fontweight='bold')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xticks([1e2, 1e3, 1e4])
        ax.set_xticklabels([r'$10^2$', r'$10^3$', r'$10^4$'])
        for pw in [1e2, 1e3, 1e4]:
            ax.axvline(x=pw, color='dimgrey', linewidth=0.8, alpha=0.8, linestyle='-')
        for ms in sorted(mpi_all_df["model_size"].unique()):
            ax.axvline(x=ms, color='lightgrey', linewidth=0.5, alpha=0.6, linestyle='--')
        # Remove per-subplot legend
        if ax.get_legend():
            ax.get_legend().remove()

    for i in range(len(num_procs), len(axes)):
        axes[i].axis('off')

    # Single shared legend for the whole figure
    seq_handle = plt.Line2D([0], [0], color='grey', marker='s', markersize=5,
                            linewidth=1.5, linestyle='--', alpha=0.7, label='Sequential')
    thread_handles = [plt.Line2D([0], [0], color=thread_palette[t], marker='o',
                      markersize=4, linewidth=2, label=t) for t in all_threads_str]
    fig.legend(handles=[seq_handle] + thread_handles, title="Threads", fontsize=9,
              loc='center right', bbox_to_anchor=(1.0, 0.5), borderaxespad=1)

    fig.suptitle(f"MPI Average Execution Time - {model_name}", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.92, 1])

    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    plt.savefig(results_dir / f"{model_name}_mpi_avg_time.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved MPI avg_time plot to results/{model_name}_mpi_avg_time.png")


def main():

    benchmarks_path = Path(argv[1]) if len(argv) > 1 else Path("benchmark_results")

    # Lists to collect all results
    all_omp_stats = []
    all_mpi_stats = []
    all_seq_stats = []

    for model_path in sorted(benchmarks_path.iterdir()):
        if model_path.is_dir():
            print(f"\n{'='*60}")
            print(f"Analysing model: {model_path.name}")
            print('='*60)
            seq_res_df, omp_res_df, mpi_res_df = analyse_model(model_path)

            if seq_res_df is None:
                print(f"Skipping {model_path.name}: No sequential results found.")
                continue
            
            # Store base stats
            seq_res_df["model_name"] = model_path.name
            all_seq_stats.append(seq_res_df)

            # Calculate metrics for this model
            if omp_res_df is not None:
                # Calculate metrics in place
                omp_speedup_df = omp_speedup(omp_res_df, seq_res_df)
                omp_efficiency_df = omp_efficiency(omp_res_df, seq_res_df)
                omp_cost_df = omp_cost(omp_res_df, seq_res_df)
                
                # Combine everything into one dataframe for this model
                omp_res_df["model_name"] = model_path.name
                omp_res_df["speedup"] = omp_speedup_df["speedup"]
                omp_res_df["efficiency"] = omp_efficiency_df["efficiency"]
                omp_res_df["cost"] = omp_cost_df["cost"]
                
                all_omp_stats.append(omp_res_df)

            if mpi_res_df is not None:
                # Calculate metrics in place
                mpi_speedup_df = mpi_speedup(mpi_res_df, seq_res_df)
                mpi_efficiency_df = mpi_efficiency(mpi_res_df, seq_res_df)
                mpi_cost_df = mpi_cost(mpi_res_df, seq_res_df)
                
                # Combine everything into one dataframe for this model
                mpi_res_df["model_name"] = model_path.name
                mpi_res_df["speedup"] = mpi_speedup_df["speedup"]
                mpi_res_df["efficiency"] = mpi_efficiency_df["efficiency"]
                mpi_res_df["cost"] = mpi_cost_df["cost"]
                
                all_mpi_stats.append(mpi_res_df)
        else:
            print(f"Skipping non-directory: {model_path.name}")

    # Process SEQ results
    seq_all_df = pd.concat(all_seq_stats, ignore_index=True)
    seq_all_df = seq_all_df[["model_name", "model_size", "avg_time", "std_time", "min_time", "max_time"]]

    # Process MPI results
    mpi_all_df = pd.concat(all_mpi_stats, ignore_index=True)
    mpi_all_df = mpi_all_df.sort_values(["model_size", "num_processors", "num_threads"])
    mpi_all_df = mpi_all_df[["model_name", "model_size", "num_processors", "num_threads", 
                              "avg_time", "std_time", "min_time", "max_time", "speedup", "efficiency", "cost"]]
    
    # Process OMP results
    omp_all_df = pd.concat(all_omp_stats, ignore_index=True)
    omp_all_df = omp_all_df.sort_values(["model_size", "num_threads"])
    omp_all_df = omp_all_df[["model_name", "model_size", "num_threads", 
                              "avg_time", "std_time", "min_time", "max_time", "speedup", "efficiency", "cost"]]
    
    print("\n" + "="*80)
    print("MPI STATS - ALL MODELS")
    print("="*80)
    print(mpi_all_df.to_string(index=False))
    
    print("\n" + "="*80)
    print("OMP STATS - ALL MODELS")
    print("="*80)
    print(omp_all_df.to_string(index=False))
    
    # Save to CSV files
    results_dir = Path("results")
    results_dir.mkdir(exist_ok=True)
    
    seq_all_df.to_csv(results_dir / "sequential_stats_all.csv", index=False)
    mpi_all_df.to_csv(results_dir / "mpi_stats_all.csv", index=False)
    omp_all_df.to_csv(results_dir / "omp_stats_all.csv", index=False)
    
    print(f"\nSaved statistics to {results_dir}")

    # Plot MPI Metrics
    plot_mpi_metric(mpi_all_df, "speedup", "Speedup", "Injection")
    plot_mpi_metric(mpi_all_df, "efficiency", "Efficiency", "Injection")
    plot_mpi_metric(mpi_all_df, "cost", "Cost", "Injection")

    # Plot OMP Metrics
    plot_omp_metric(omp_all_df, "speedup", "Speedup", "Injection")
    plot_omp_metric(omp_all_df, "efficiency", "Efficiency", "Injection")
    plot_omp_metric(omp_all_df, "cost", "Cost", "Injection")

    # Plot Average Time
    plot_omp_avg_time(omp_all_df, seq_all_df, "Injection")
    plot_mpi_avg_time(mpi_all_df, seq_all_df, "Injection")

    #Print the name of the model and the sizes in increasing order
    print("\n" + "="*80)
    print("MODEL SIZES (in increasing order)")
    print("="*80)
    model_sizes = mpi_all_df[["model_name", "model_size"]].drop_duplicates().sort_values("model_size")
    print(model_sizes.to_string(index=False))

if __name__ == "__main__":
    main()
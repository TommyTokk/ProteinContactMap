# Parallelisation analysis of Protein Contact Maps calculation

Author: Tommaso Tocchini<br>
Student ID: 618135<br>
A.Y.:2025-2026

---
# Introduction
The repository contains the project for the module of **Computing and Networking: Resources and Tools**, from the master degree of Informatics for Digital Health at Univeristy of Pisa.

The project address the challenge to analyse the performances of various kinds of execution mode, comparing:
- The **sequential version** of the code
- The **OpenMP implementation** of the code *varying the number of threads*
- An **hybrid implementation of the code using MPI + OpenMP** *varying both the number of threads and processors*

# Structure of the repository
```
.
├── README.md
├── makefile
├── mpi_main.cpp
├── omp_main.cpp
├── seq_main.cpp
├── bin/
├── build/
├── lib/
├── src/
├── data/
├── scripts/
├── results/
├── benchmark_results/
└── plots/
```

# Problem Description
A protein contact map is a simplified two-dimensional binary representation of a protein’s three-dimensional structure. By expressing tertiary structure as a 2D matrix, it preserves key spatial relationships and topology while discarding absolute Cartesian coordinates. This abstraction reduces dimensionality and removes rotational and translational variance, making contact maps well suited to computational workflows, especially modern deep learning methods for structure prediction, fold recognition, and design. For a protein consisting of $N$ amino acid residues, the contact map is represented as an $N \times N$ symmetric, binary matrix $C$. The individual matrix elements $C_{ij}$ are defined based on a spatial proximity threshold:
$$\begin{equation}
    C_{ij} = \begin{cases} 
        1 & \text{if } d_{ij} \le d_c \\ 
        0 & \text{otherwise} 
    \end{cases}
\end{equation}$$
where $d_{ij}$ denotes the spatial distance between designated reference points of residues $i$ and $j$, and $d_c$ represents a predefined threshold distance.The selection of $d_c$ typically ranges from $8\,\AA{}$ to $12\,\AA{}$, depending on the specific reference points selected and the physical interactions of interest (such as hydrogen bonding, hydrophobic packing, or salt bridges).

## Computation of the distances
In the project, the distances are computed using the $C_\alpha$ atoms. Given the coordinates $(x_i, y_i, z_i)$ and $(x_j, y_j, z_j)$ for the chosen $C_\alpha$ atoms of residues $i$ and $j$, the distance is defined as:
$$
\begin{equation}
    d_{ij} = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2}
\end{equation}
$$

For computational reason, the real computation coded is $$(d_{i,j})^2 = (x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2$$

In this way the individual matrix entry $C_{i,j}$ is calculated as
$$\begin{equation}
    C_{ij} = \begin{cases} 
        1 & \text{if } d_{ij}^2 \le d_c^2 \\ 
        0 & \text{otherwise} 
    \end{cases}
\end{equation}$$

# Project details
To compare the different execution policies, three approaches have been tested:
- Execution with all the optimisations applied
- Execution removing the optimisations from the sequential code
- Execution with extra workload injected

## Execution with optimisations
In this experiments multiple kinds of optimisations have been applied to the code. First of all, instead of using an Array of Structure approach, I used a Structure of Array approach; in this way the memory is better exploited since less access are made. To do so the *Model* structure has been used:
```c++
typedef struct Model{
    std::vector<float> X;
    std::vector<float> Y;
    std::vector<float> Z;

    void resize(int n){
        X.resize(n);
        Y.resize(n);
        Z.resize(n);
    }
}Model;
```

### Sequential code
In this experiment, the method used in the ```seq_main``` file is
```std::vector<uint8_t> get_residue_distances_soaV2(const Model& model, int size)```.
This method exploits a Structure of Arrays (SoA) layout and removes if/else
branches as a code optimization. Additional optimizations are enabled by
compiling with the -O3 flag, which improves performance by leveraging, for
example, vectorization.

### OMP code
For the OMP implementation, the method used in ```omp_main``` is ```std::vector<uint8_t> get_residue_distances_omp_soa(const Model& model, int size, int n_threads)```. As with the sequential code, this implementation exploits the SoA layout and the removal of if/else branches to improve performance. The parallel for is declared with **static** scheduling. *Since the number of chunks is not explicitly specified, the iteration space is divided into approximately equal-sized contiguous chunks, with at most one chunk assigned to each thread*. In the inner loop, the *SIMD* directive is used to instruct the compiler to generate SIMD instructions. This may be redundant, since the ```-O3``` compilation flag should already enable vectorization.

### Hybrid code
For the hybrid implementation, a combination of MPI and OpenMP (OMP) has been used. The idea is that the **master process** (*rank 0*) splits the workload among the worker processes and then participates in the computation itself. At the end, the master gathers the results using the ```MPI_Gatherv``` function.
To achieve an effective partitioning, the matrix has been divided into row chunks. Given an $N \times N$ matrix and $P$ workers, each worker processes an $(N/P) \times N$ submatrix. The function used to compute the distances in the submatrix is:
```c++
std::vector<uint8_t> get_residue_distances_mpi_soa(
    const Model& model,
    int size,
    size_t starting_row,
    size_t count,
    int n_threads)
```

Inside this function, OMP threads are used to parallelize the inner computation.

### Run the experiment
<div style="background: #fff3cd; padding: 12px; border-left: 4px solid #ffc107; border-radius: 4px; color: black;">
<strong style="color: black;">⚠️ Warning:</strong> 
Before running this experiment, ensure that each main file uses the correct methods and that you compile with the -O3 flag in the CXXFLAGS line of the makefile.
</div>
<br>

The repository includes a Makefile to simplify compilation and execution. After verifying that the correct methods are being used, run
```
make
```

to compile all versions of the code. This will create three files inside the bin folder.
To run a complete benchmark, the ```benchmark_all``` script is provided. This script can be used to execute multiple iterations of the various main programs, saving information about each run. The script accepts four arguments:
```
--all  : Run all benchmarks (sequential, OpenMP, and MPI)
--seq  : Run only the sequential benchmark
--omp  : Run only the OpenMP benchmark
--mpi  : Run only the MPI benchmark
```
To run the full benchmark suite, execute:
```./scripts/benchmark_all.sh --all```

Unless specified otherwise, this will create a *benchmark_results* folder containing all information about the executions.
Finally, to obtain statistics such as *speedup*, *efficiency*, and *cost*, and to generate plots of these statistics, run the ```results_script.py``` file.

## Execution without optimisation
Instead of adopting the Structure of Arrays (SoA) layout, I chose to use the *Atom* structure. 
```c++
typedef struct Atom{
    char record_name[7]; //Contains Atom or HETATM
    int serial_number; //Serial number
    char name[5]; //Atom name
    char alternate_locator_indicator[2];
    char residue_name[4]; // Name of the residue
    char chain_id[2]; //Chain identifier
    int res_seq; // Residue serial number
    char code_residue_insert[2]; // Code for insertions of residue
    float x,y,z; // Coordinates of the atom
    float occupancy; //Occupancy
    float temp; // Temperature factor
    char segment_id[5]; // Segment identifier
    char element_sym[3]; // Element symbol
    char charge[3]; // Charge
    int model;
}Atom;
```

This choice means that atomic coordinates are no longer stored contiguously in memory, which increases the number of memory accesses and introduces additional execution overhead.
I also introduced an explicit if statement, which incurs a slight overhead due to branch evaluation. Furthermore, I removed the ```-O3``` optimisation flag, so optimisation techniques such as vectorisation and loop unrolling are no longer applied by the compiler.

### Sequential code
In this experiment, the only code that needs to be modified is the sequential code. The method used in the ```seq_main``` file is:
```c++
std::vector<uint8_t> dm = get_residue_distances(alphas_vec)
```
This function is implemented without any form of optimization.

### Run the experiment
<div style="background: #fff3cd; padding: 12px; border-left: 4px solid #ffc107; border-radius: 4px; color: black;">
<strong style="color: black;">⚠️ Warning:</strong> 
Before running this experiment, ensure that each main file uses the correct methods and that you compile the sequential code <b>WITHOUT</b> the O3 flag. <br>The OMP code and the hybrid code <b>MUST</b> be compiled using the O3 flag instead.
</div>
<br>

First of all, to compile the OMP code and the hybrid code using the optimisations run the command
```bash
make omp mpi
```

after that, remove the ```-O3``` flag from the compilation flags (*you can just comment line 5 and uncomment line 8 in the Makefile*).

After having modified the Makefile, run 
```
make seq
```
to compile the sequential code without the compiler optimisation.

<div style="background: #de97ff; padding: 12px; border-left: 4px solid #bd07ff; border-radius: 4px; color: black;">
<strong style="color: black;">Note:</strong> 
Before compile the sequential code, be sure to use the not optimised version of the PCM computation in the main file.
</div>
<br>

Finally, to run the benchmark and toget the stats of the experiment, run
```
./scripts/benchmark_all.sh
./scripts/results_script.py
```

## Execution with workload injection
For the final experiment, I injected additional code into the distance calculation to simulate a heavier computational workload. Specifically, I added 100 floating-point additions and subtractions:
```c++
// Extra operations to increase execution time
for (int z = 0; z < N_ITERATIONS; z++) {
    dist_sq += 0.000001f;
    dist_sq -= 0.000001f;
}
```

This modification was applied to all versions of the code. The rationale was that increasing the computational cost per iteration might improve the relative performance of the parallelised implementations.

### Sequential code
In this experiment, the method used in the seq_main file is:
```c++
std::vector<uint8_t> get_residue_distances_seq_inj(const Model& model, int size)
```
This function follows the same logic of the functions that uses SoA optimisation, but adds, after the computation of the distance, the extra workload.

### OMP code
For the OMP implementation, the method used in ```omp_main``` is 
```c++
std::vector<uint8_t> get_residue_distances_omp_inj(const Model& model, int size, int n_threads)
```
Also in this case, the method has the same logic of the one used in the fully optimised experiment, but with the addition of the extra workload.

### Hybrid implementation
For the hybrid implementation, the function invoked in ```mpi_main``` is
```c++
std::vector<uint8_t> get_residue_distances_mpi_inj(const Model& model, int size, size_t starting_row, size_t count, int n_threads)
```
In this scenario as well, the function follows the same logic as the one used in the fully optimized experiment, with the only difference being the inclusion of the additional workload.

### Run the experiment
<div style="background: #fff3cd; padding: 12px; border-left: 4px solid #ffc107; border-radius: 4px; color: black;">
<strong style="color: black;">⚠️ Warning:</strong> 
Before running this experiment, ensure that each main file uses the correct methods and that you compile with the -O3 flag in the CXXFLAGS line of the makefile.
</div>
<br>

As in the previous experiments, after verifying that you are using the correct methods, run:
```
make
```
to compile all versions of the code. Then execute:
```
./scripts/benchmark_all.sh
./scripts/results_script.py
```
to run the benchmark and process the results.
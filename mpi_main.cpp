#include "lib/utils.hpp"
#include <chrono>
#include <iostream>
#include "mpi.h"
#include <numeric>


int main(int argc, char const **argv){

    // Initialize MPI
    MPI_Init(&argc, (char***)&argv);

    // Get the number of processes and the rank of the current process
    int n_procs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    
    // Check if the input and output paths and number of threads are provided
    if(argc <= 3){
        if(rank == 0) printf("[ERROR] USAGE: ./main <input_path> <output_dir_path> n_threads\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Vector to hold the alpha carbon atoms (only filled on the root process)
    std::vector<Atom> alphas_vec;

    // Get the input file path, output directory, and number of threads from command line arguments
    const char *file_path = argv[1];
    const char *output_dir = argv[2];
    const int n_threads = atoi(argv[3]);

    // Variables for measuring execution time
    double parallel_start, parallel_end, parallel_time;

    // Vectors to hold the counts and starting indices for each process
    std::vector<int> counts(n_procs), starts(n_procs);
    
    // Variables to hold the number of alpha carbon atoms and the filename (only filled on the root process)
    int alphas_size;
    std::string pdb_filename;

    if(rank == 0){ // Only the root process loads the data

        // Open the input file and check if it was opened successfully
        FILE *fptr = fopen(file_path, "r");
        if (!fptr) {
            printf("[ERROR] Cannot open file\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Get filename without path and extension
        pdb_filename = get_filename(file_path);

        // Load atoms from the input file
        std::unordered_map<int, std::vector<Atom>> atoms = load_atoms_from_file(fptr);
        
        fclose(fptr);// Close the file after loading

        // Get alpha carbon atoms from the loaded atoms and check if there are any
        std::unordered_map<int, std::vector<Atom>> alphas = get_alphas(atoms);
        if(alphas.empty()){
            printf("[ERROR] No alpha carbon atoms found in the input file.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Get the first entry of the alphas map to access the vector of alpha carbon atoms and its size
        auto iterator = alphas.begin();
        alphas_vec = iterator->second;
        alphas_size = alphas_vec.size();

        // Calculate the bounds for each process
        for(int i = 0; i < n_procs; i++){
            size_t start, count;
            get_processor_bounds(alphas_size, n_procs, i, start, count);
            starts[i] = start;
            counts[i] = count;
        }
    }

    //Broadcasting the number of elements in the array
    MPI_Bcast(&alphas_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Create Model with SOA layout
    Model m;
    m.resize(alphas_size);

    // Only the root process fills the Model struct with the coordinates of the alpha carbon atoms
    if (rank == 0) {
        for(int i = 0; i < alphas_size; i++){
            m.X[i] = alphas_vec[i].x;
            m.Y[i] = alphas_vec[i].y;
            m.Z[i] = alphas_vec[i].z;
        }
    }

    parallel_start = MPI_Wtime();

    // Broadcast the Coordinate vectors (SOA)
    MPI_Bcast(m.X.data(), alphas_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m.Y.data(), alphas_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(m.Z.data(), alphas_size, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Broadcasting the counts
    MPI_Bcast(counts.data(), n_procs, MPI_INT, 0, MPI_COMM_WORLD);

    // Scatter the starting indices to each process
    int my_start;
    int my_count = counts[rank];
    MPI_Scatter(starts.data(), 1, MPI_INT,
                &my_start, 1, MPI_INT,
                0, MPI_COMM_WORLD);

    // Calculate the disposals
    std::vector<int> displs(n_procs);
    std::exclusive_scan(counts.begin(), counts.end(), displs.begin(), 0);

    // Calculate the local distance matrix for each process
    // Comment the following line during benchmarking to test the MPI version of the code with injected workload
    std::vector<uint8_t> local_dm = get_residue_distances_mpi_soa(m, alphas_size, my_start, my_count, n_threads);

    // Uncomment the following line to test the MPI version of the code with injected workload
    //std::vector<uint8_t> local_dm = get_residue_distances_mpi_inj(m, alphas_size, my_start, my_count, n_threads);

    // Count non-zero elements in local_dm for debugging
    int nz = 0;
    for(auto el : local_dm){
        if(el > 0) nz++;
    }

    // Recalculate the counts
    std::vector<int> flat_counts(n_procs);
    for (int i = 0; i < flat_counts.size(); i++){
        flat_counts[i] = counts[i]*alphas_size;
    }

    // Recalculate the disposals
    std::vector<int> flat_displs(n_procs);
    std::exclusive_scan(flat_counts.begin(), flat_counts.end(), 
                   flat_displs.begin(), 0);

    // Verify local_dm size matches the expected count
    if (local_dm.size() != static_cast<size_t>(flat_counts[rank])) {
        printf(
            "[ERROR] Rank %d: local_dm size mismatch. Expected %d, got %lu\n",
            rank,
            flat_counts[rank],
            static_cast<unsigned long>(local_dm.size())
        );

        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Gather the local distance matrices into the all_results vector on the root process
    std::vector<uint8_t> all_results;
    if(rank == 0){
        int total = flat_displs[n_procs-1] + flat_counts[n_procs-1];
        all_results.resize(total);
        printf("Allocated all_results with size: %d\n", total);
    }

    // Gather the local distance matrices into the all_results vector on the root process
    MPI_Gatherv(
        flat_counts[rank] > 0 ? local_dm.data() : nullptr,// Handle the case where a process has no data to send
        flat_counts[rank], // send count
        MPI_UNSIGNED_CHAR, // send type
        rank == 0? all_results.data() : nullptr,// receive buffer (only valid on root) 
        flat_counts.data(), // receive counts
        flat_displs.data(), // displacements
        MPI_UNSIGNED_CHAR, // receive type
        0, // root
        MPI_COMM_WORLD // communicator 
    );

    parallel_end = MPI_Wtime();
    parallel_time = parallel_end - parallel_start;


    // Save the data
    if(rank == 0){
        int total_atoms = alphas_size;
        
        // Verify size
        if(all_results.size() != total_atoms * total_atoms){
            printf("ERROR: Size mismatch! Expected %d, got %lu\n", 
                total_atoms * total_atoms, all_results.size());
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Save the matrix directly (already in flat row-major format)

        // Comment this line during benchmarking to avoid usless wait
        save_distance_matrix(all_results, total_atoms, output_dir, pdb_filename);
    }

    if(rank == 0){
        printf("alg-time | %.10f | s\n", parallel_time);
    }

    MPI_Finalize();
    
    return EXIT_SUCCESS;
}

#include "lib/utils.hpp"
#include <chrono>
#include <iostream>
#include "mpi.h"
#include <numeric>


int main(int argc, char const **argv){

    

    //MPI init
    MPI_Init(&argc, (char***)&argv);

    // Create MPI datatype for Atom struct
    MPI_Datatype MPI_ATOM;
    MPI_Type_contiguous(sizeof(Atom), MPI_BYTE, &MPI_ATOM);
    MPI_Type_commit(&MPI_ATOM);

    //Getting the number of proces
    int n_procs, rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    std::vector<Atom> alphas_vec;

    const char *file_path = argv[1];
    const char *output_dir = argv[2];
    const int n_threads = atoi(argv[3]);

    if(argc <= 3){
        if(!rank){
            printf("[ERROR] USAGE: ./main <input_path> <output_dir_path> n_threads\n");
            return EXIT_FAILURE;
        }
        
    }

    double service_start, service_end, service_time;
    double parallel_start, parallel_end, parallel_time;
    
    service_start = MPI_Wtime();

    std::vector<int> counts(n_procs), starts(n_procs);
    
    int alphas_size;
    std::string pdb_filename;

    if(rank == 0){ 
        FILE *fptr = fopen(file_path, "r");
        if (!fptr) {
            printf("[ERROR] Cannot open file\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        pdb_filename = get_filename(file_path);

        std::unordered_map<int, std::vector<Atom>> atoms = load_atoms_from_file(fptr);
        
        fclose(fptr);

        std::unordered_map<int, std::vector<Atom>> alphas = get_alphas(atoms);

        if(alphas.empty()){
            printf("[ERROR] No alpha carbon atoms found in the input file.\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        auto iterator = alphas.begin();
        alphas_vec = iterator->second;
        alphas_size = alphas_vec.size();

        for(int i = 0; i < n_procs; i++){
            size_t start, count;
            get_processor_bounds(alphas_size, n_procs, i, start, count);
            starts[i] = start;
            counts[i] = count;
        }
    }

    //Broadcasting the number of elements in the array
    MPI_Bcast(&alphas_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    alphas_vec.resize(alphas_size);


    //Broadcast the elements
    MPI_Bcast(alphas_vec.data(), alphas_size, MPI_ATOM, 0, MPI_COMM_WORLD);

    // Create Model with SOA layout
    Model m;
    m.resize(alphas_size);
    for(int i = 0; i < alphas_size; i++){
        m.X[i] = alphas_vec[i].x;
        m.Y[i] = alphas_vec[i].y;
        m.Z[i] = alphas_vec[i].z;
    }

    //Broadcasting the counts
    MPI_Bcast(counts.data(), n_procs, MPI_INT, 0, MPI_COMM_WORLD);

    int my_start;
    int my_count = counts[rank];

    MPI_Scatter(starts.data(), 1, MPI_INT,
                &my_start, 1, MPI_INT,
                0, MPI_COMM_WORLD);

    //Calculate the disposals
    std::vector<int> displs(n_procs);
    std::exclusive_scan(counts.begin(), counts.end(), displs.begin(), 0);

    service_end = MPI_Wtime();
    parallel_start = MPI_Wtime();

    printf("Rank: %d, starting row: %d, count: %d, end row: %d\n", rank, my_start, my_count, (my_start+my_count));

    std::vector<uint8_t> local_dm = get_residue_distances_mpi_soa(m, alphas_size, my_start, my_count, n_threads);

    int nz = 0;
    for(auto el : local_dm){
        if(el > 0) nz++;
    }

    printf("Rank: %d, non zero el: %d\n", rank, nz);

    // Recalculate the counts
    std::vector<int> flat_counts(n_procs);
    for (int i = 0; i < flat_counts.size(); i++){
        flat_counts[i] = counts[i]*alphas_size;
    }

    // Recalculate the disposals
    std::vector<int> flat_displs(n_procs);
    std::exclusive_scan(flat_counts.begin(), flat_counts.end(), 
                   flat_displs.begin(), 0);

    // Gather the data
    std::vector<uint8_t> all_results;
    if(rank == 0){
        int total = flat_displs[n_procs-1] + flat_counts[n_procs-1];
        all_results.resize(total);
        printf("Allocated all_results with size: %d\n", total);
    }

    MPI_Gatherv(
        local_dm.data(), 
        flat_counts[rank], 
        MPI_UNSIGNED_CHAR, 
        rank == 0? all_results.data() : nullptr, 
        flat_counts.data(), 
        flat_displs.data(), 
        MPI_UNSIGNED_CHAR, 
        0, 
        MPI_COMM_WORLD
    );

    parallel_end = MPI_Wtime();

    service_time = service_end - service_start;
    parallel_time = parallel_end - parallel_start;

    if(rank == 0){
        printf("\n=== Timing Results ===\n");
        printf("Service time: %.6f(s)\n", service_time);
        printf("Parallel time: %.6f(s)\n", parallel_time);
        printf("Total time: %.6f(s)\n", service_time + parallel_time);
        
        // Amdahl's Law calculations
        double f_serial = service_time / (service_time + parallel_time);
        double f_parallel = parallel_time / (service_time + parallel_time);
        double theoretical_speedup = 1.0 / (f_serial + f_parallel / (n_procs * n_threads));
        double efficiency = (theoretical_speedup / (n_procs * n_threads)) * 100.0;
        
        printf("Serial fraction: %.4f\n", f_serial);
        printf("Parallel fraction: %.4f\n", f_parallel);
        printf("Total parallelism: %d (procs) x %d (threads) = %d\n", n_procs, n_threads, n_procs * n_threads);
        printf("Theoretical max speedup (Amdahl): %.2fx\n", theoretical_speedup);
        printf("Efficiency: %.2f%%\n", efficiency);
        printf("======================\n\n");
    }

    // Save the data
    if(rank == 0){
        int total_atoms = alphas_size;
        
        // Verify size
        if(all_results.size() != total_atoms * total_atoms){
            printf("ERROR: Size mismatch! Expected %d, got %zu\n", 
                total_atoms * total_atoms, all_results.size());
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Save the matrix directly (already in flat row-major format)
        //save_distance_matrix(all_results, total_atoms, output_dir, pdb_filename);
    }

    


    /*TODO:
    
        [1,2,3,4,5,6,7,8] -> 
        -> [
            [0, 1, 2, 3, 4, 5, 6, 7]
            [1, 0, 1, 2, 3, 4, 5, 6]
            [2, 1, 0, 1, 2, 3, 4 ,5]
            [3, 2, 1, 0, 1, 2, 3, 4]
            [4, 3, 2, 1, 0, 1, 2, 3]
            [5, 4, 3, 2, 1, 0, 1, 2]
            [6, 5, 4, 3, 2, 1, 0, 1]
            [7, 6 ,5 ,4, 3, 2, 1, 0]
            ]


        r0 -> [1,2] ->
                    -> [
                        [0, 1]
                        [1, 0]
                        ] X

        r1 -> [2,3] -> 
                    -> [[0,1] [10]]

        r0 -> [1,2] -> 
                    -> [
                        [0, 1, 2, 3, 4, 5, 6, 7]
                        [1, 0, 1, 2, 3, 4, 5, 6]
                        ]
    
    */



    MPI_Type_free(&MPI_ATOM);
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}

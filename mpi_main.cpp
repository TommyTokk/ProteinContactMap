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

    printf("Rank: %d, starting row: %d, count: %d, end row: %d\n", rank, my_start, my_count, (my_start+my_count));

    std::vector<std::vector<float>> local_dm = get_residue_distances_mpi(alphas_vec, my_start, my_count, n_threads);


    int nz = 0;
    for(const auto& arr : local_dm){
        for(const auto el: arr){
            if(el > 0) nz++;
        }
    }

    printf("Rank: %d, non zero el pre: %d\n", rank, nz);

    // Flatten the partial results
    std::vector<float> flattened;
    flattened.reserve(my_count * alphas_size);

    for(const auto& arr : local_dm){
        flattened.insert(flattened.end(), arr.begin(), arr.end());
    }

    nz = 0;
    for(auto el : flattened){
        if(el > 0) nz += 1;
    }

    printf("Rank: %d, non zero el: %d\n", rank, nz);

    // Recalculate the counts
    std::vector<int> flat_counts(n_procs);
    for (int i = 0; i < flat_counts.size(); i++){
        printf("%d, %d\n", counts[i], alphas_size);
        flat_counts[i] = counts[i]*alphas_size;
    }

    // Recalculate the disposals
    std::vector<int> flat_displs(n_procs);
    std::exclusive_scan(flat_counts.begin(), flat_counts.end(), 
                   flat_displs.begin(), 0);

    // Gather the data
    std::vector<float> all_results;
    if(rank == 0){
        int total = flat_displs[n_procs-1] + flat_counts[n_procs-1];
        all_results.resize(total);
        printf("Allocated all_results with size: %d\n", total);
    }

    MPI_Gatherv(
        flattened.data(), 
        flat_counts[rank], 
        MPI_FLOAT, 
        rank == 0? all_results.data() : nullptr, 
        flat_counts.data(), 
        flat_displs.data(), 
        MPI_FLOAT, 
        0, 
        MPI_COMM_WORLD
    );

    // Rebuild the data
    if(rank == 0){
        int total_atoms = alphas_size;  // 10510
        
        // Verify size
        if(all_results.size() != total_atoms * total_atoms){
            printf("ERROR: Size mismatch! Expected %d, got %zu\n", 
                total_atoms * total_atoms, all_results.size());
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        
        // Rebuild the matrix from flattened row-major data
        std::vector<std::vector<float>> distance_matrix(total_atoms, std::vector<float>(total_atoms));
        
        for(int row = 0; row < total_atoms; row++){
            for(int col = 0; col < total_atoms; col++){
                distance_matrix[row][col] = all_results[row * total_atoms + col];
            }
        }
        
        // Save the matrix
        save_distance_matrix(distance_matrix, output_dir, pdb_filename);
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

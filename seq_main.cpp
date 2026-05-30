#include "lib/utils.hpp"
#include <chrono>
#include <iostream>
#include <omp.h>


int main(int argc, char const *argv[]){
    
    // Variables for measuring execution time
    double service_start, service_end, service_time;

    // Check if the input and output paths are provided
    if(argc <= 2){
        printf("[ERROR] USAGE: ./main <input_path> <output_dir_path>\n");
        return EXIT_FAILURE;
    }

    // Get the input file path and output directory from command line arguments
    const char *file_path = argv[1];
    const char *output_dir = argv[2];

    // Open the input file
    FILE *fptr = fopen(file_path, "r");

    // Check if the file was opened successfully
    if(!fptr){
        printf("[ERROR] An error occurred during loading the file.\n");
        return EXIT_FAILURE;
    }

    // Get filename without path and extension
    std::string pdb_filename = get_filename(file_path);

    // Load atoms from the input file
    std::unordered_map<int, std::vector<Atom>>atoms = load_atoms_from_file(fptr);

    // Get alpha carbon atoms from the loaded atoms
    std::unordered_map<int, std::vector<Atom>> alphas = get_alphas(atoms);

    // Check if there are alpha carbon atoms and compute the distance matrix
    if(!alphas.empty()){

        // Get the first entry of the alphas map to access the vector of alpha carbon atoms
        auto iterator = alphas.begin();

        // Get the vector of alpha carbon atoms and its size
        std::vector<Atom> alphas_vec = iterator->second;
        int alphas_size = alphas_vec.size();

        // Create a Model struct and fill it with the coordinates of the alpha carbon atoms
        Model m;
        m.resize(alphas_size);
        for(int i = 0; i < alphas_size; i++){
            m.X[i] = alphas_vec[i].x;
            m.Y[i] = alphas_vec[i].y;
            m.Z[i] = alphas_vec[i].z;
        }

        // Measure the execution time of the distance matrix computation
        service_start = omp_get_wtime(); 
        
        // Compute the distance matrix using the SoA optimized version
        std::vector<uint8_t> dm = get_residue_distances_soaV2(m, alphas_size);
        
        /* 
        Uncomment this line if exectution without SoA optimisation is required
        WARNING: This is the version used in the report in the section "3.2 Execution without optimisations",
        if these line is uncommented, please remember to remove the "-O3" flag from the Makefile, otherwise the execution time will be too short and it will not be possible to see the difference between the two versions.
        std::vector<uint8_t> dm = get_residue_distances(alphas_vec);
        */

        /*
        Uncomment this line if the workload injection is required

        std::vector<uint8_t> dm = get_residue_distances_seq_inj(m, alphas_size);
        */

        service_end = omp_get_wtime();

        //Commment this line if you don't want to save the distance matrix to a file (e.g. during benchmarking) 
        save_distance_matrix(dm, alphas_size, output_dir, pdb_filename);

    }else{
        printf("[ERROR] No alpha carbon atoms found in the input file.\n");
        fclose(fptr);
        return EXIT_FAILURE;
    }

    // Calculate the execution time and print it
    service_time = service_end - service_start; 

    printf("alg-time | %.10f | s\n", service_time);

    

    fclose(fptr);
    return EXIT_SUCCESS;
}

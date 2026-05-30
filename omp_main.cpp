#include "lib/utils.hpp"
#include <chrono>
#include <iostream>


int main(int argc, char const *argv[]){

    // Variables for measuring execution time
    double parallel_start, parallel_end, parallel_time;

    // Check if the input and output paths and number of threads are provided
    if(argc <= 3){
        printf("[ERROR] USAGE: ./main <input_path> <output_dir_path> n_threads\n");
        return EXIT_FAILURE;
    }

    // Get the input file path, output directory, and number of threads from command line arguments
    const char *file_path = argv[1];
    const char *output_dir = argv[2];
    const int n_threads = atoi(argv[3]);

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
        std::vector<Atom> alphas_vec = iterator->second;
        int alphas_size = alphas_vec.size();

        // Create a Model struct and fill it with the coordinates of the alpha carbon atoms
        Model m;
        m.resize(alphas_size);

        // Fill the Model struct with the coordinates of the alpha carbon atoms
        for(int i = 0; i < alphas_size; i++){
            m.X[i] = alphas_vec[i].x;
            m.Y[i] = alphas_vec[i].y;
            m.Z[i] = alphas_vec[i].z;
        }
        parallel_start = omp_get_wtime();

        // Compute the distance matrix using the OMP SoA optimized version
        std::vector<uint8_t> dm = get_residue_distances_omp_soa(m, alphas_size, n_threads);


        /* Use the following lines to test the OMP version of the code with injected workload
        std::vector<uint8_t> dm = get_residue_distances_omp_inj(m, alphas_size, n_threads);
        */

        parallel_end = omp_get_wtime();

        // Comment this line during benchmarking to avoid usless wait
        save_distance_matrix(dm, alphas_size, output_dir, pdb_filename);
    }else{
        printf("[ERROR] No alpha carbon atoms found in the input file.\n");
        fclose(fptr);
        return EXIT_FAILURE;
    }

    parallel_time = parallel_end - parallel_start;

    printf("alg-time | %.10f | s\n", parallel_time);




    fclose(fptr);
    return EXIT_SUCCESS;
}

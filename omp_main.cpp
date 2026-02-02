#include "lib/utils.hpp"
#include <chrono>
#include <iostream>


int main(int argc, char const *argv[]){
    double service_start, service_end, service_time;
    double parallel_start, parallel_end, parallel_time;

    service_start = omp_get_wtime(); 

    if(argc <= 3){
        printf("[ERROR] USAGE: ./main <input_path> <output_dir_path> n_threads\n");
        return EXIT_FAILURE;
    }



    const char *file_path = argv[1];
    const char *output_dir = argv[2];
    const int n_threads = atoi(argv[3]);

    FILE *fptr = fopen(file_path, "r");

    if(!fptr){
        printf("[ERROR] An error occurred during loading the file.\n");
        return EXIT_FAILURE;
    }

    // Get filename without path and extension
    std::string pdb_filename = get_filename(file_path);

    std::unordered_map<int, std::vector<Atom>>atoms = load_atoms_from_file(fptr);

    std::unordered_map<int, std::vector<Atom>> alphas = get_alphas(atoms);

    if(!alphas.empty()){

        auto iterator = alphas.begin();
        std::vector<Atom> alphas_vec = iterator->second;

        int alphas_size = alphas_vec.size();

        Model m;

        m.resize(alphas_size);

        for(int i = 0; i < alphas_size; i++){
            m.X[i] = alphas_vec[i].x;
            m.Y[i] = alphas_vec[i].y;
            m.Z[i] = alphas_vec[i].z;
        }

        service_end = omp_get_wtime();
        parallel_start = omp_get_wtime();

        //std::vector<std::vector<float>> dm = get_residue_distances_omp(alphas_vec,0,alphas_vec.size(), n_threads);
        //std::vector<uint8_t> dm = get_residue_distances_omp_opt(alphas_vec, n_threads);
        std::vector<uint8_t> dm = get_residue_distances_omp_soa(m, alphas_size, n_threads);

        parallel_end = omp_get_wtime();

        //save_distance_matrix(dm, alphas_size, output_dir, pdb_filename);
    }else{
        printf("[ERROR] No alpha carbon atoms found in the input file.\n");
        fclose(fptr);
        return EXIT_FAILURE;
    }

    service_time = service_end - service_start;
    parallel_time = parallel_end - parallel_start;

    double total_time = service_time + parallel_time;

    double lb = service_time + (parallel_time/n_threads);
    double ub = total_time/lb;

    printf("Service time: %.6f(s)\n", service_time);
    printf("Parallel time: %.6f(s)\n", parallel_time);
    printf("Lower bound time: %.6f(s)\n", lb);
    printf("Total time: %.6f(s)\n", total_time);
    printf("Upper bound: %.6f\n", ub);

    



    fclose(fptr);
    return EXIT_SUCCESS;
}

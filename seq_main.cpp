#include "lib/utils.hpp"
#include <chrono>
#include <iostream>


int main(int argc, char const *argv[]){
    if(argc <= 2){
        printf("[ERROR] USAGE: ./main <input_path> <output_dir_path>\n");
        return EXIT_FAILURE;
    }

    const char *file_path = argv[1];
    const char *output_dir = argv[2];

    FILE *fptr = fopen(file_path, "r");

    if(!fptr){
        printf("[ERROR] An error occurred during loading the file.\n");
        return EXIT_FAILURE;
    }

    // Get filename without path and extension
    std::string pdb_filename = get_filename(file_path);

    std::unordered_map<int, std::vector<Atom>>atoms = load_atoms_from_file(fptr);

    std::unordered_map<int, std::vector<Atom>> alphas = get_alphas(atoms);

    const auto start{std::chrono::steady_clock::now()};

    std::unordered_map<int, std::vector<std::vector<float>>> dm = get_residue_distances(alphas);

    const auto finish{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{finish - start};
    
    std::cout<<std::to_string(elapsed_seconds.count()) + " (s)"<<std::endl;
    
    save_distance_matrix(dm, output_dir, pdb_filename);

    fclose(fptr);
    return EXIT_SUCCESS;
}

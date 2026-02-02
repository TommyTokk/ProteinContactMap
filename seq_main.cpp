#include "lib/utils.hpp"
#include <chrono>
#include <iostream>


int main(int argc, char const *argv[]){
    const auto start{std::chrono::steady_clock::now()};
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

        

        //std::vector<uint8_t> dm = get_residue_distances_opt(alphas_vec);
        std::vector<uint8_t> dm = get_residue_distances_soaV2(m, alphas_size);

        //save_distance_matrix(dm, alphas_size, output_dir, pdb_filename);

    }else{
        printf("[ERROR] No alpha carbon atoms found in the input file.\n");
        fclose(fptr);
        return EXIT_FAILURE;
    }

    const auto finish{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{finish - start};
        
    std::cout<<std::to_string(elapsed_seconds.count()) + "(s)"<<std::endl;

    

    fclose(fptr);
    return EXIT_SUCCESS;
}

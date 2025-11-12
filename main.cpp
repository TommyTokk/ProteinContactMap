#include "lib/utils.hpp"


int main(int argc, char const *argv[]){
    if(argc <= 1){
        printf("[ERROR] USAGE: ./main <input_path> <output_dir_path\n");
        return EXIT_FAILURE;
    }

    const char *file_path = argv[1];

    FILE *fptr = fopen(file_path, "r");

    if(!fptr){
        printf("[ERROR] An error occured during loading the file.\n");
        return EXIT_FAILURE;
    }

    std::unordered_map<int, std::vector<Atom>>atoms = load_atoms_from_file(fptr);

    std::unordered_map<int, std::vector<Atom>> alphas = get_alphas(atoms);

    std::unordered_map<int, std::vector<std::vector<float>>> dm = get_residue_distances(alphas);


    for(const std::pair<int, std::vector<std::vector<float>>> t : dm){
        int model = t.first;
        printf("Model: %d\n", model);
    }

    

    //Compute the distances
    

    //free_atoms_vector(atoms);
    fclose(fptr);
    return(EXIT_SUCCESS);
}

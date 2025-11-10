#include "lib/utils.hpp"

int main(int argc, char const *argv[]){
    if(argc == 1) return EXIT_FAILURE;

    const char *file_path = argv[1];

    FILE *fptr = fopen(file_path, "r");

    if(!fptr){
        printf("[ERROR] An error occured during loading the file.\n");
        return EXIT_FAILURE;
    }

    std::queue<Atom *> atoms = load_atom_queue_from_file(fptr);

    

    free_queue(atoms);
    fclose(fptr);
    return(EXIT_SUCCESS);
}

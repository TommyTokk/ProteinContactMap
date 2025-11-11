#include "lib/utils.hpp"

int main(int argc, char const *argv[]){
    if(argc == 1) return EXIT_FAILURE;

    const char *file_path = argv[1];

    FILE *fptr = fopen(file_path, "r");

    if(!fptr){
        printf("[ERROR] An error occured during loading the file.\n");
        return EXIT_FAILURE;
    }

    std::vector<Atom *> atoms = load_atom_queue_from_file(fptr);

    std::vector<Atom *> alphas = get_alphas(atoms);

    std::vector<Atom *> alphas_by_residue = get_alphas_by_residues(atoms);

    for(Atom *atom: alphas){
        printf("%s | %d | (%f, %f, %f)\n", atom->name, atom->res_seq, atom->x, atom->y, atom->z);
    }
    //Insert the computation part

    //Compute the distances
    

    free_atoms_vector(atoms);
    fclose(fptr);
    return(EXIT_SUCCESS);
}

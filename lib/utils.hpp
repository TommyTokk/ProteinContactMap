#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.hpp"
#include <queue>

std::vector<Atom *>load_atom_queue_from_file(FILE *fptr);
std::vector<Atom *> get_alphas(std::vector<Atom *> v);
std::vector<std::pair<int, Atom *>> get_alphas_by_residues(std::vector<Atom *> v);
int **get_contact_map(std::vector<Atom *> alphas);
void free_atoms_vector(std::vector<Atom *> queue);
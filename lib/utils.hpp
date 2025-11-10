#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.hpp"
#include <queue>

std::queue<Atom *>load_atom_queue_from_file(FILE *fptr);
void free_queue(std::queue<Atom *> queue);
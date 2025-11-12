#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.hpp"
#include <queue>
#include <unordered_map>
#include<math.h>
#include <fstream>
#include <iostream>

std::unordered_map<int, std::vector<Atom>> load_atoms_from_file(FILE *fptr);
std::unordered_map<int, std::vector<Atom>> get_alphas(std::unordered_map<int, std::vector<Atom>> map);
std::unordered_map<int, std::vector<std::vector<float>>> get_residue_distances(std::unordered_map<int, std::vector<Atom>> alphas);


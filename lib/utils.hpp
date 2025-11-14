#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.hpp"
#include <queue>
#include <unordered_map>
#include<math.h>
#include <fstream>
#include <iostream>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

std::unordered_map<int, std::vector<Atom>> load_atoms_from_file(FILE *fptr);
std::unordered_map<int, std::vector<Atom>> get_alphas(std::unordered_map<int, std::vector<Atom>> map);
std::unordered_map<int, std::vector<std::vector<float>>> get_residue_distances(std::unordered_map<int, std::vector<Atom>> alphas);
std::string get_filename(const char* path);
void save_csv(const std::vector<std::vector<float>> distance_matrix, const char *filepath);
void save_distance_matrix(const std::unordered_map<int, std::vector<std::vector<float>>> &dm, const char *output_dir, const std::string &pdb_filename);

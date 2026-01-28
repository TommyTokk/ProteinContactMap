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
#include <unistd.h>
#include <cstdio>

std::unordered_map<int, std::vector<Atom>> load_atoms_from_file(FILE *fptr);
std::unordered_map<int, std::vector<Atom>> get_alphas(std::unordered_map<int, std::vector<Atom>> map);
std::vector<std::vector<float>> get_residue_distances(std::vector<Atom> alphas);
std::vector<std::vector<float>> get_residue_distances_omp(std::vector<Atom> alphas, size_t start, size_t size, int n_threads);
std::vector<std::vector<float>> get_residue_distances_mpi(const std::vector<Atom>& alphas, size_t starting_row, size_t count, int n_threads);
void get_processor_bounds(size_t total_elements, int num_processors, int my_rank, size_t& start_index, size_t& count);
std::string get_filename(const char* path);
void save_csv(const std::vector<std::vector<float>> &distance_matrix, const char *filepath);
void save_distance_matrix(std::vector<std::vector<float>> &dm, const char *output_dir, const std::string &pdb_filename);

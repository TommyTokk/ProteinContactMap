#include "../lib/utils.hpp"


std::unordered_map<int, std::vector<Atom>> load_atoms_from_file(FILE *fptr){
	char line[81];
	std::vector<Atom> atoms;
	int model = 1;
	std::unordered_map<int, std::vector<Atom>> map;

	while (fgets(line, sizeof(line), fptr)){

		// Read the atom model
		if (strncmp(line, "MODEL", 5) == 0){
			sscanf(line + 10, "%4d", &model);
		}

		if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0){
			Atom atom;

			// Storing the record name
			strncpy(atom.record_name, line + 0, 6);
			atom.record_name[6] = '\0';

			// Storing the serial number
			sscanf(line + 6, "%5d", &atom.serial_number);

			// Storing the name
			strncpy(atom.name, line + 13, 4);
			atom.name[5] = '\0';

			// Storing the alternate locator
			strncpy(atom.alternate_locator_indicator, line + 16, 1);
			atom.alternate_locator_indicator[1] = '\0';

			// Storing the name of the residue
			strncpy(atom.residue_name, line + 17, 3);
			atom.residue_name[4] = '\0';

			// Storing the chain ID
			strncpy(atom.chain_id, line + 21, 1);
			atom.chain_id[2] = '\0';

			// Storing the residue sequence number
			sscanf(line + 22, "%4d", &atom.res_seq);

			// Storing the code for insertion of the residues
			strncpy(atom.code_residue_insert, line + 26, 1);
			atom.code_residue_insert[2] = '\0';

			// Storing the coordinates
			// Storing the x coordinate
			sscanf(line + 30, "%8f", &atom.x);

			// Storing the y coordinate
			sscanf(line + 38, "%8f", &atom.y);

			// Storing the z coordinate
			sscanf(line + 46, "%8f", &atom.z);

			// Storing the occupancy
			sscanf(line + 54, "%6f", &atom.occupancy);

			// Storing the temperature factor
			sscanf(line + 60, "%6f", &atom.temp);

			// Storing the segment id
			strncpy(atom.segment_id, line + 72, 4);
			atom.segment_id[5] = '\0';

			// Storing the element symbol
			strncpy(atom.element_sym, line + 76, 2);
			atom.element_sym[3] = '\0';

			// Storing the charge
			strncpy(atom.charge, line + 78, 2);
			atom.charge[3] = '\0';

			// Storing the model info
			atom.model = model;

			// Load the atom in the map
			map[model].push_back(atom);
			
		}
	}

	return map;
}

std::unordered_map<int, std::vector<Atom>> get_alphas(std::unordered_map<int, std::vector<Atom>> map){

	std::unordered_map<int, std::vector<Atom>> res;

	for (const std::pair<int, std::vector<Atom>>t : map){
		int k = t.first;
		std::vector<Atom> atoms = t.second;

		//Finding the alphas
		for(auto atom: atoms){
			if(strcasestr(atom.name, "CA")){
				res[k].push_back(atom);
			}
		}
	}
	return res;
}


std::unordered_map<int, std::vector<std::vector<float>>> get_residue_distances(std::unordered_map<int, std::vector<Atom>> alphas_map){
    std::unordered_map<int, std::vector<std::vector<float>>> res;
    
    for(const std::pair<int, std::vector<Atom>>t : alphas_map){
        int model = t.first;
        const std::vector<Atom>& alphas = t.second;
        size_t num_atoms = alphas.size();

        if(num_atoms > 0){
            std::vector<std::vector<float>> distance_matrix(num_atoms, std::vector<float>(num_atoms));

            for(size_t i = 0; i < num_atoms; i++){
                for(size_t j = 0; j < num_atoms; j++){
                    float x_diff = alphas[i].x - alphas[j].x;
                    float y_diff = alphas[i].y - alphas[j].y;
                    float z_diff = alphas[i].z - alphas[j].z;

                    distance_matrix[i][j] = sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2));
                }
            }
            res[model] = std::move(distance_matrix);
        }
    }

    return res;
}

//std::string get_filename(const std::string& path) {
//    size_t last_slash = path.find_last_of("/\\");
//    if (last_slash != std::string::npos) {
//        return path.substr(last_slash + 1);
//    }
//    return path;
//}

#include "../lib/utils.hpp"

/**
 * Loads the atoms from a PDB file and stores them in a map, where the key is the model number and the value is a vector of atoms belonging to that model. The function reads the PDB file line by line, parsing the relevant information for each atom and storing it in an Atom struct. The atoms are then grouped by their model number in the map.
 * 
 * @param fptr A file pointer to the opened PDB file
 * @return An unordered map where the key is the model number and the value is a vector of Atom structs representing the atoms in that model
 */
std::unordered_map<int, std::vector<Atom>> load_atoms_from_file(FILE *fptr){
	char line[81];
	int model = 1;
	std::unordered_map<int, std::vector<Atom>> map;

	while (fgets(line, sizeof(line), fptr)){

		// Read the atom model
		if (strncmp(line, "MODEL", 5) == 0){
			sscanf(line + 10, "%4d", &model);
		}

		if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0){
			Atom atom = {};

			// Storing the record name
			strncpy(atom.record_name, line + 0, 6);
			atom.record_name[6] = '\0';

			// Storing the serial number
			sscanf(line + 6, "%5d", &atom.serial_number);

			// Storing the name
			strncpy(atom.name, line + 13, 4);
			atom.name[4] = '\0';

			// Storing the alternate locator
			strncpy(atom.alternate_locator_indicator, line + 16, 1);
			atom.alternate_locator_indicator[1] = '\0';

			// Storing the name of the residue
			strncpy(atom.residue_name, line + 17, 3);
			atom.residue_name[3] = '\0';

			// Storing the chain ID
			strncpy(atom.chain_id, line + 21, 1);
			atom.chain_id[1] = '\0';

			// Storing the residue sequence number
			sscanf(line + 22, "%4d", &atom.res_seq);

			// Storing the code for insertion of the residues
			strncpy(atom.code_residue_insert, line + 26, 1);
			atom.code_residue_insert[1] = '\0';

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
			atom.segment_id[4] = '\0';

			// Storing the element symbol
			strncpy(atom.element_sym, line + 76, 2);
			atom.element_sym[2] = '\0';

			// Storing the charge
			strncpy(atom.charge, line + 78, 2);
			atom.charge[1] = '\0';

			// Storing the model info
			atom.model = model;

			// Load the atom in the map
            
			map[model].push_back(atom);
			
		}
	}

	return map;
}

/**
 * Filters the input map of atoms to extract only the alpha carbon atoms (CA). The function iterates through the input map, checking each atom's name for the presence of "CA" (case-insensitive). If an atom is identified as an alpha carbon, it is added to a new map that groups the alpha carbon atoms by their model number. The resulting map contains only the alpha carbon atoms from the original input.
 * 
 * @param map An unordered map where the key is the model number and the value is a vector of Atom structs representing the atoms in that model
 * @return An unordered map where the key is the model number and the value is a vector of Atom structs representing only the alpha carbon atoms in that model
 */
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

//================================================
// SEQUENTIAL IMPLEMENTATIONS
//================================================
/**
 * Naive implementation of the distance matrix calculation. It uses a vector unsigned integers to store the distance matrix.
 * The distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms.
 * 
 * The function is created not thinking about any kind of optimisation.
 * 
 * @param alphas The vector of alpha carbon atoms
 * @return A 2D vector representing the distance matrix, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise.
 */
std::vector<uint8_t> get_residue_distances(const std::vector<Atom>& alphas) {
    size_t num_atoms = alphas.size();
    std::vector<uint8_t> distance_vector(num_atoms * num_atoms);
    
    float cutoff_sq = 64.0f; 

    for(size_t i = 0; i < num_atoms; i++) {
        // Retrieve coordinates of the i-th atom
        float xi = alphas[i].x;
        float yi = alphas[i].y;
        float zi = alphas[i].z;
        
        // Pre-calculate row offset
        size_t row_offset = i * num_atoms;
        
        // Calculate distances for the current row
        for(size_t j = 0; j < num_atoms; j++) {
            // Retrieve coordinates of the j-th atom
            float dx = xi - alphas[j].x;
            float dy = yi - alphas[j].y;
            float dz = zi - alphas[j].z;
            
            // Compute squared distance and compare with cutoff
            float dist_sq = dx*dx + dy*dy + dz*dz;
            
            // Store result in the distance vector
            if (dist_sq <= cutoff_sq) {
                distance_vector[row_offset + j] = 1;
            } else {
                distance_vector[row_offset + j] = 0;
            }
        }
    }
    
    return distance_vector;
}



/**
 * Optimized version of the distance matrix calculation using Structure of Arrays (SoA) approach. It uses a vector unsigned integers to store the distance matrix.
 * 
 * The distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms.
 * 
 * The function is optimized for better cache performance by using a Structure of Arrays (SoA) approach, where the coordinates of the atoms are stored in separate vectors. This allows for better spatial locality when accessing the coordinates during distance calculations.
 * 
 * The function stores directly the boolean value of the contact, avoiding if/else statements taht could lead to branch evaluation penalties.
 * 
 * @param model The model containing the coordinates of the alpha carbon atoms in a Structure of Arrays (SoA) format
 * @param size The number of alpha carbon atoms
 * @return A vector representing the distance matrix, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and
 */
std::vector<uint8_t> get_residue_distances_soaV2(const Model& model, int size){
    std::vector<uint8_t> distance_vector(size * size);
    
    float cutoff_sq = 64.0f;
    
    // Calculate distances for each pair of atoms
    for(size_t i = 0; i < size; i++) {
        // Retrieve coordinates of the i-th atom
        float xi = model.X[i];
        float yi = model.Y[i];
        float zi = model.Z[i];
        
        // Set self-distance to 1 (contact with itself)
        distance_vector[i * size + i] = 1;
        
        for(size_t j = 0; j < size; j++) {
            // Retrieve coordinates of the j-th atom
            float dx = xi - model.X[j];
            float dy = yi - model.Y[j];
            float dz = zi - model.Z[j];
            
            // Compute squared distance and compare with cutoff
            float dist_sq = dx*dx + dy*dy + dz*dz;

            // Store result in the distance vector (1 if contact, 0 otherwise)
            uint8_t contact = (dist_sq <= cutoff_sq);
            
            // Store the contact result in the distance vector
            distance_vector[i * size + j] = contact;
        }
    }
    
    return distance_vector;
}

/**
 * Implementation of the distance matrix calculation with workload injection. It uses a vector unsigned integers to store the distance matrix.
 * 
 * The distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms.
 * 
 * The function is similar to the optimized version using Structure of Arrays (SoA) approach, but it includes a workload injection loop that performs redundant calculations to artificially increase the execution time. Used for banchmarking and testing the performance of the implementation under heavier workloads.
 * 
 * @param model The model containing the coordinates of the alpha carbon atoms in a Structure of Arrays (SoA) format
 * @param size The number of alpha carbon atoms
 * @return A vector representing the distance matrix, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise
 */
std::vector<uint8_t> get_residue_distances_seq_inj(const Model& model, int size){
    std::vector<uint8_t> distance_vector(size * size);
    
    float cutoff_sq = 64.0f;
    
    for(size_t i = 0; i < size; i++) {
        // Retrieve coordinates of the i-th atom
        float xi = model.X[i];
        float yi = model.Y[i];
        float zi = model.Z[i];
        
        // Set self-distance to 1 (contact with itself)
        distance_vector[i * size + i] = 1;
        
        for(size_t j = 0; j < size; j++) {
            // Retrieve coordinates of the j-th atom
            float dx = xi - model.X[j];
            float dy = yi - model.Y[j];
            float dz = zi - model.Z[j];
            
            // Compute squared distance and compare with cutoff
            float dist_sq = dx*dx + dy*dy + dz*dz;

            //Useless code to increase the execution time
            // N_ITERATIONS is defined in the utils.hpp file and can be adjusted to increase or decrease the amount of workload injected
            for(int z = 0; z < N_ITERATIONS; z++){
                dist_sq += 0.000001f;
                dist_sq -= 0.000001f;
            }

            // Store result in the distance vector (1 if contact, 0 otherwise)
            uint8_t contact = (dist_sq <= cutoff_sq);
            
            // Store the contact result in the distance vector
            distance_vector[i * size + j] = contact;
        }
    }
    
    return distance_vector;
}
//================================================
// SEQUENTIAL IMPLEMENTATIONS
//================================================

//================================================
// OMP IMPLEMENTATIONS
//================================================

/**
 * WARNING: This function is not used in the current implementation, and it is not optimized for performance. It is provided as a reference for a possible OpenMP implementation of the distance matrix calculation using the AoS approach
 * 
 * The distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms.
 * 
 * The function uses OpenMP to parallelize the outer loop over the atoms, allowing for concurrent computation of distance rows. Each thread computes rows of the distance matrix.
 * 
 * The function is optimized for better cache performance by using a Structure of Arrays (SoA) approach, where the coordinates of the atoms are stored in separate vectors. This allows for better spatial locality when accessing the coordinates during distance calculations.
 * 
 * @param alphas The vector of alpha carbon atoms
 * @param n_threads The number of threads to use for parallelization
 * @return A vector representing the distance matrix, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise
 */
[[unused]]
std::vector<uint8_t> get_residue_distances_omp_opt(const std::vector<Atom>& alphas, int n_threads) {
    // Get the number of atoms and initialize the distance vector
    size_t num_atoms = alphas.size();
    std::vector<uint8_t> distance_vector(num_atoms * num_atoms);
    
    float cutoff_sq = 64.0f; 

    // Parallelize the outer loop using OpenMP
    #pragma omp parallel for num_threads(n_threads)
    for(size_t i = 0; i < num_atoms; i++) {
        // Retrieve coordinates of the i-th atom
        float xi = alphas[i].x;
        float yi = alphas[i].y;
        float zi = alphas[i].z;
        
        // Pre-calculate row offset
        size_t row_offset = i * num_atoms;
        
        for(size_t j = 0; j < num_atoms; j++) {
            // Retrieve coordinates of the j-th atom
            float dx = xi - alphas[j].x;
            float dy = yi - alphas[j].y;
            float dz = zi - alphas[j].z;
            
            // Compute squared distance and compare with cutoff
            float dist_sq = dx*dx + dy*dy + dz*dz;
            
            // Store result in the distance vector (1 if contact, 0 otherwise)
            distance_vector[row_offset + j] = (dist_sq <= cutoff_sq);
        }
    }
    
    return distance_vector;
}


/**
 * Optimized version of the distance matrix calculation using OpenMP for parallelization and Structure of Arrays (SoA) approach. It uses a vector unsigned integers to store the distance matrix.
 * 
 * The distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms.
 * 
 * The function uses OpenMP to parallelize the outer loop over the atoms, allowing for concurrent computation of distance rows. Each thread computes rows of the distance matrix.
 * 
 * The function is optimized for better cache performance by using a Structure of Arrays (SoA) approach, where the coordinates of the atoms are stored in separate vectors. This allows for better spatial locality when accessing the coordinates during distance calculations.
 * 
 * @param model The model containing the coordinates of the alpha carbon atoms in a Structure of Arrays (SoA) format
 * @param size The number of alpha carbon atoms
 * @param n_threads The number of threads to use for parallelization
 * @return A vector representing the distance matrix, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise.
 */
std::vector<uint8_t> get_residue_distances_omp_soa(const Model& model, int size, int n_threads){
    std::vector<uint8_t> distance_vector(size * size);
    
    float cutoff_sq = 64.0f;
    
    // Parallelize the outer loop using OpenMP
    // schedule(static) is used because each iteration has uniform cost; it assigns contiguous blocks of rows to threads, reducing scheduling overhead and improving locality when writing the output matrix.
    #pragma omp parallel for num_threads(n_threads) schedule(static)
    for(size_t i = 0; i < size; i++) {

        // Retrieve coordinates of the i-th atom
        float xi = model.X[i];
        float yi = model.Y[i];
        float zi = model.Z[i];
        
        // Pre-calculate row offset
        size_t row_offset = i * size;

        #pragma omp simd
        // The inner loop is vectorized using OpenMP SIMD directive, which allows the compiler to generate SIMD instructions for the distance calculations, further improving performance by processing multiple distance calculations in parallel.
        for(size_t j = 0; j < size; j++) {
            // Retrieve coordinates of the j-th atom
            float dx = xi - model.X[j];
            float dy = yi - model.Y[j];
            float dz = zi - model.Z[j];
            
            // Compute squared distance and compare with cutoff
            float dist_sq = dx*dx + dy*dy + dz*dz;
            
            // Store result in the distance vector (1 if contact, 0 otherwise)
            distance_vector[row_offset + j] = (dist_sq <= cutoff_sq);
        }
    }
    
    return distance_vector;
}

/**
 * Implementation of the distance matrix calculation with workload injection using OpenMP for parallelization. It uses a vector unsigned integers to store the distance matrix.
 * 
 * The distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms.
 * 
 * The function is similar to the optimized version using OpenMP and Structure of Arrays (SoA) approach, but it includes a workload injection loop that performs redundant calculations to artificially increase the execution time. Used for banchmarking and testing the performance of the implementation under heavier workloads.
 * 
 * @param model The model containing the coordinates of the alpha carbon atoms in a Structure of Arrays (SoA) format
 * @param size The number of alpha carbon atoms
 * @param n_threads The number of threads to use for parallelization
 * @return A vector representing the distance matrix, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise
 */
std::vector<uint8_t> get_residue_distances_omp_inj(const Model& model, int size, int n_threads){
    std::vector<uint8_t> distance_vector(size * size);
    
    float cutoff_sq = 64.0f;
    
    #pragma omp parallel for num_threads(n_threads) schedule(static)
    for(size_t i = 0; i < size; i++) {
        float xi = model.X[i];
        float yi = model.Y[i];
        float zi = model.Z[i];
        

        size_t row_offset = i * size;

        #pragma omp simd
        for(size_t j = 0; j < size; j++) {
            float dx = xi - model.X[j];
            float dy = yi - model.Y[j];
            float dz = zi - model.Z[j];
            
            float dist_sq = dx*dx + dy*dy + dz*dz;

            //Useless code to increase the execution time
            for(int z = 0; z < N_ITERATIONS; z++){
                dist_sq += 0.000001f;
                dist_sq -= 0.000001f;
            }
            
            // Store result in the distance vector (1 if contact, 0 otherwise)
            distance_vector[row_offset + j] = (dist_sq <= cutoff_sq);
        }
    }
    
    return distance_vector;
}
//================================================
// OMP IMPLEMENTATIONS
//================================================

//================================================
// MPI IMPLEMENTATIONS
//================================================

/**
 * 
 * Implementation of the distance matrix calculation using MPI for distributed parallelization and OpenMP for shared-memory parallelization. It uses a vector of vectors of floats to store the distance matrix.
 * 
 * WARNING: This function is not used in the current implementation, and it is not optimized for performance. It is provided as a reference for a possible MPI implementation of the distance matrix calculation, but it is not recommended for use in production code due to its inefficiency and the fact that it uses a vector of vectors to store the distance matrix, which can lead to poor cache performance and increased memory usage.
 * 
 * NOTE: For the functions used in the main code, please refer to the get_residue_distances_mpi_soa and get_residue_distances_mpi_inj functions, which are optimized for performance and use a flat vector to store the distance matrix.
 */
[[unused]]
std::vector<std::vector<float>> get_residue_distances_mpi(const std::vector<Atom>& alphas, size_t starting_row, size_t count, int n_threads){
    size_t num_atoms = alphas.size();
    std::vector<std::vector<float>> distance_matrix(count, std::vector<float>(num_atoms));

    #pragma omp parallel for num_threads(n_threads)
    for(size_t i = 0; i < count; i++){ 
        size_t global_row = starting_row + i;
        
        for(size_t c = 0; c < num_atoms; c++){
            float x_diff = alphas[global_row].x - alphas[c].x;
            float y_diff = alphas[global_row].y - alphas[c].y;
            float z_diff = alphas[global_row].z - alphas[c].z;

            float distance_sq = (x_diff*x_diff) + (y_diff*y_diff) + (z_diff*z_diff);

            distance_matrix[i][c] = (distance_sq <= 64.0f) ? 1.0f : 0.0f;
        }
    }

    return distance_matrix;
}

/**
 * Optimized version of the distance matrix calculation using MPI for distributed parallelization and OpenMP for shared-memory parallelization, with Structure of Arrays (SoA) approach. It uses a vector unsigned integers to store the distance matrix.
 * 
 * The partial distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms, for the subset of rows assigned to the current MPI process.
 * 
 * The function uses OpenMP to parallelize the outer loop over the assigned rows of atoms, allowing for concurrent computation of distance rows. Each thread computes rows of the distance matrix.
 * 
 * The function is optimized for better cache performance by using a Structure of Arrays (SoA) approach, where the coordinates of the atoms are stored in separate vectors. This allows for better spatial locality when accessing the coordinates during distance calculations.
 * 
 * @param model The model containing the coordinates of the alpha carbon atoms in a Structure of Arrays (SoA) format
 * @param size The total number of alpha carbon atoms
 * @param starting_row The index of the first row of the distance matrix that this MPI process is responsible for computing
 * @param count The number of rows of the distance matrix that this MPI process is responsible for computing
 * @param n_threads The number of threads to use for parallelization
 * @return A vector representing the partial distance matrix computed by this MPI process, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise
 */
std::vector<uint8_t> get_residue_distances_mpi_soa(const Model& model, int size, size_t starting_row, size_t count, int n_threads){
    std::vector<uint8_t> distance_vector(count * size, 0); 
    
    float cutoff_sq = 64.0f;

    // Parallelize the outer loop using OpenMP, where each thread computes a subset of the assigned rows. The loop iterates over the rows assigned to this MPI process, and for each row, it computes the distances to all other atoms.
    #pragma omp parallel for num_threads(n_threads)
    for(size_t i = 0; i < count; i++){
        // Calculate the global row index in the distance matrix for the current iteration
        size_t global_row = starting_row + i;
        
        // Retrieve the coordinates of the atom corresponding to the current row
        float xi = model.X[global_row];
        float yi = model.Y[global_row];
        float zi = model.Z[global_row];
        
        // Pre-calculate the row offset for storing results in the distance vector
        size_t row_offset = i * size;
        
        // Vectorize the inner loop using OpenMP SIMD directive, which allows the compiler to generate SIMD instructions for the distance calculations, further improving performance by processing multiple distance calculations in parallel.
        #pragma omp simd
        for(size_t j = 0; j < size; j++){
            float dx = xi - model.X[j];
            float dy = yi - model.Y[j];
            float dz = zi - model.Z[j];
            
            float dist_sq = dx*dx + dy*dy + dz*dz;
            distance_vector[row_offset + j] = (dist_sq <= cutoff_sq);
        }
    }
    

    return distance_vector;
}

/**
 * Implementation of the distance matrix calculation using MPI for distributed parallelization and OpenMP for shared-memory parallelization, with Structure of Arrays (SoA) approach and workload injection. It uses a vector unsigned integers to store the distance matrix.
 * 
 * The partial distance matrix is stored in a 1D vector, where the element at index (i * num_atoms + j) corresponds to the distance between the i-th and j-th alpha carbon atoms, for the subset of rows assigned to the current MPI process.
 * 
 * The function is similar to the optimized version using MPI and OpenMP with Structure of Arrays (SoA) approach, but it includes a workload injection loop that performs redundant calculations to artificially increase the execution time. Used for banchmarking and testing the performance of the implementation under heavier workloads.
 * 
 * @param model The model containing the coordinates of the alpha carbon atoms in a Structure of Arrays (SoA) format
 * @param size The total number of alpha carbon atoms
 * @param starting_row The index of the first row of the distance matrix that this MPI process is responsible for computing
 * @param count The number of rows of the distance matrix that this MPI process is responsible for computing
 * @param n_threads The number of threads to use for parallelization
 * @return A vector representing the partial distance matrix computed by this MPI process, where each element is 1 if the distance between the corresponding atoms is less than or equal to 8 angstroms, and 0 otherwise
 */
std::vector<uint8_t> get_residue_distances_mpi_inj(const Model& model, int size, size_t starting_row, size_t count, int n_threads){
    std::vector<uint8_t> distance_vector(count * size, 0); 
    
    float cutoff_sq = 64.0f;

    #pragma omp parallel for num_threads(n_threads)
    for(size_t i = 0; i < count; i++){
        // Calculate the global row index in the distance matrix for the current iteration
        size_t global_row = starting_row + i;
        
        // Retrieve the coordinates of the atom corresponding to the current row
        float xi = model.X[global_row];
        float yi = model.Y[global_row];
        float zi = model.Z[global_row];
        
        // Pre-calculate the row offset for storing results in the distance vector
        size_t row_offset = i * size;
        
        #pragma omp simd
        for(size_t j = 0; j < size; j++){
            // Retrieve coordinates of the j-th atom
            float dx = xi - model.X[j];
            float dy = yi - model.Y[j];
            float dz = zi - model.Z[j];

            // Compute the injected workload to artificially increase the execution time. The loop performs redundant calculations that do not affect the final result, but they consume CPU cycles, allowing for testing the performance of the implementation under heavier workloads. 
            // N_ITERATIONS is defined in the utils.hpp file and can be adjusted to increase or decrease the amount of workload injected.
            for(int z = 0; z < N_ITERATIONS; z++){
                dx += 0.000001f;
                dx -= 0.000001f;
            }
            
            // Compute squared distance and compare with cutoff
            float dist_sq = dx*dx + dy*dy + dz*dz;
            distance_vector[row_offset + j] = (dist_sq <= cutoff_sq);
        }
    }
    

    return distance_vector;
}
//================================================
// MPI IMPLEMENTATIONS
//================================================


//================================================
// UTILS FUNCTIONS
//================================================

void get_processor_bounds(size_t total_elements, int num_processors, int my_rank, 
                          size_t& start_index, size_t& count) {
    
    size_t base_count = total_elements / num_processors;
    size_t remainder = total_elements % num_processors;

    //Calculate the number of elements for this processor
    if (my_rank < remainder) {
        count = base_count + 1;
    } else {
        count = base_count;
    }

    // The offset is: (rank * base) + (how many 'extra' items appeared before me)
    start_index = (my_rank * base_count) + std::min((size_t)my_rank, remainder);
}

std::string get_filename(const char* path) {
    std::string str_path(path);
    size_t last_slash = str_path.find_last_of("/\\");
    size_t last_dot = str_path.find_last_of(".");
    
    std::string filename = (last_slash != std::string::npos) 
        ? str_path.substr(last_slash + 1) 
        : str_path;
    
    if (last_dot != std::string::npos && last_dot > last_slash) {
        filename = filename.substr(0, last_dot - (last_slash + 1));
    }
    
    return filename;
}

void save_distance_matrix(
    const std::vector<uint8_t>& dm,
    size_t num_atoms,
    const char* output_dir,
    const std::string& pdb_filename
) {
    // Create subdirectory: output_dir/pdb_filename/
    const std::string subdir = std::string(output_dir) + "/" + pdb_filename;

    if (mkdir(subdir.c_str(), 0755) != 0 && errno != EEXIST) {
        std::fprintf(
            stderr,
            "Error: Could not create directory %s: %s\n",
            subdir.c_str(),
            std::strerror(errno)
        );
        return;
    }

    // Output path: output_dir/pdb_filename/pdb_filename.csv
    const std::string output_path = subdir + "/" + pdb_filename + ".csv";

    std::printf("Saving distance matrix to: %s\n", output_path.c_str());

    // Verify matrix size
    const size_t expected_size = num_atoms * num_atoms;

    if (dm.size() != expected_size) {
        std::fprintf(
            stderr,
            "Error: Matrix size mismatch. Expected %zu, got %zu\n",
            expected_size,
            dm.size()
        );
        return;
    }

    // Open file using POSIX I/O
    int fd = open(output_path.c_str(), O_CREAT | O_WRONLY | O_TRUNC, 0644);

    if (fd < 0) {
        std::fprintf(
            stderr,
            "Error: Could not open file %s: %s\n",
            output_path.c_str(),
            std::strerror(errno)
        );
        return;
    }

    // Each row contains:
    // num_atoms digits + (num_atoms - 1) commas + 1 newline
    std::string row_buffer;
    row_buffer.reserve(num_atoms * 2 + 1);

    for (size_t r = 0; r < num_atoms; ++r) {
        row_buffer.clear();

        const size_t row_offset = r * num_atoms;

        for (size_t c = 0; c < num_atoms; ++c) {
            row_buffer.push_back(dm[row_offset + c] ? '1' : '0');

            if (c + 1 < num_atoms) {
                row_buffer.push_back(',');
            }
        }

        row_buffer.push_back('\n');

        ssize_t written = write(fd, row_buffer.data(), row_buffer.size());

        if (written < 0) {
            std::fprintf(
                stderr,
                "Error: Failed while writing to file %s: %s\n",
                output_path.c_str(),
                std::strerror(errno)
            );

            close(fd);
            return;
        }

        if (static_cast<size_t>(written) != row_buffer.size()) {
            std::fprintf(
                stderr,
                "Error: Partial write while writing to file %s\n",
                output_path.c_str()
            );

            close(fd);
            return;
        }
    }

    if (close(fd) != 0) {
        std::fprintf(
            stderr,
            "Error: Failed to close file %s: %s\n",
            output_path.c_str(),
            std::strerror(errno)
        );

        return;
    }

    std::printf("Matrix saved successfully!\n");
}

//================================================
// UTILS FUNCTIONS
//================================================


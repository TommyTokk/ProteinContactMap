#include "../lib/utils.hpp"


std::queue<Atom *> load_atom_queue_from_file(FILE *fptr){
    char line[81];
    std::queue<Atom *> atom_queue;
    while (fgets(line, sizeof(line), fptr)) {
        if (strncmp(line, "ATOM", 4) == 0 || strncmp(line, "HETATM", 6) == 0){
            Atom *atom = (Atom *)calloc(1, sizeof(Atom));

            if(!atom){
                printf("Error during the allocation of the Atom structure!\n");
                return atom_queue;
            }

            // Storing the record name
            strncpy(atom->record_name, line + 0, 6);
            atom->record_name[6] = '\0';

            // Storing the serial number
            sscanf(line + 6, "%5d", &atom->serial_number);

            //Storing the name
            strncpy(atom->name, line + 13, 4);
            atom->name[5] = '\0';


	    //Storing the alternate locator
	    strncpy(atom->alternate_locator_indicator, line + 16, 1);
	    atom->alternate_locator_indicator[1] = '\0';

	    //Storing the name of the residue
	    strncpy(atom->residue_name, line + 17, 3);
	    atom->residue_name[4] = '\0';

	    //Storing the chain ID
	    strncpy(atom->chain_id, line + 21, 1);
	    atom->chain_id[2] = '\0';

	    //Storing the residue sequence number
	    sscanf(line + 22, "%4d",&atom->res_seq);

	    //Storing the code for insertion of the residues
	    strncpy(atom->code_residue_insert, line + 26, 1);
	    atom->code_residue_insert[2] = '\0';

	    //Storing the coordinates
	    //Storing the x coordinate
	    sscanf(line + 30, "%8f", &atom->x);

	    //Storing the y coordinate
	    sscanf(line + 38, "%8f", &atom->y);

	    //Storing the z coordinate
	    sscanf(line + 46, "%8f", &atom->z);

	    //Storing the occupancy
	    sscanf(line + 54, "%6f", &atom->occupancy);
	    
	    //Storing the temperature factor
	    sscanf(line + 60, "%6f", &atom->temp);

	    //Storing the segment id
	    
	    strncpy(atom->segment_id, line + 72, 4);
	    atom->segment_id[5] = '\0';

	    //Storing the element symbol
	    strncpy(atom->element_sym, line + 76, 2);
	    atom->element_sym[3] = '\0';

	    //Storing the charge 
	    
	    strncpy(atom->charge, line + 78, 2);
	    atom->charge[3] = '\0';
		
	    //Load the atom in the queue
	    atom_queue.push(atom);

	}
    }

    return atom_queue;
}


void free_queue(std::queue<Atom *> queue){
	while(!queue.empty()){
		Atom *atom = queue.front();
		queue.pop();

		free(atom);
	}
}

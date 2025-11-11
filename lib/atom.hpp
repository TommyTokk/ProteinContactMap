typedef struct Atom{
    char record_name[7]; //Contains Atom or HETATM
    int serial_number; //Serial number
    char name[5]; //Atom name
    char alternate_locator_indicator[2];
    char residue_name[4]; // Name of the residue
    char chain_id[2]; //Chain identifier
    int res_seq; // Residue serial number
    char code_residue_insert[2]; // Code for insertions of residue
    float x,y,z; // Coordinates of the atom
    float occupancy; //Occupancy
    float temp; // Temperature factor
    char segment_id[5]; // Segment identifier
    char element_sym[3]; // Element symbol
    char charge[3]; // Charge
}Atom;
#include "Hitran.ih"

void Hitran::readFixedWidthFile() {
    ifstream file(filepath);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filepath);
    }

    string line;
    while (getline(file, line)) {
        // Parse the fixed-width fields (based on _HITRANformat)
        vector<string> row;
        try {
            // The widths as per the Python version: [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15, 15, 15, 6, 12, 1, 7, 7]
            row.push_back(line.substr(0, 2));    // Molecule_ID
            row.push_back(line.substr(2, 1));    // Isotopologue_ID
            row.push_back(line.substr(3, 12));   // nu
            row.push_back(line.substr(15, 10));  // S
            row.push_back(line.substr(25, 10));  // A
            row.push_back(line.substr(35, 5));   // gamma_air
            row.push_back(line.substr(40, 5));   // gamma_self
            row.push_back(line.substr(45, 10));  // E_l
            row.push_back(line.substr(55, 4));   // n_air
            row.push_back(line.substr(59, 8));   // del_air
            row.push_back(line.substr(67, 15));  // global_u
            row.push_back(line.substr(82, 15));  // global_l
            row.push_back(line.substr(97, 15));  // local_u
            row.push_back(line.substr(112, 15)); // local_l
            row.push_back(line.substr(127, 6));  // err_ind
            row.push_back(line.substr(133, 12)); // References
            row.push_back(line.substr(145, 1));  // line_mixing
            row.push_back(line.substr(146, 7));  // g_u
            row.push_back(line.substr(153, 7));  // g_l

            // Add the parsed row to the dataset
            data.push_back(row);
        } catch (const exception& e) {
            cerr << "Error parsing line: " << line << endl;
        }
    }

    file.close();
}
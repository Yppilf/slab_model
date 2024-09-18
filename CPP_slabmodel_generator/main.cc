#include "main.ih"

void testMolecule(){
    try {
        Molecule water("H2O", 298.15); // Create a Molecule object for water
        cout << "Molar weight of H2O: " << water.getMolarWeight() << " g/mol" << endl;

        Molecule methane("CH4", 298.15); // Create a Molecule object for methane
        cout << "Molar weight of CH4: " << methane.getMolarWeight() << " g/mol" << endl;

        Molecule ethylene("C2H4", 298.15); // Create a Molecule object for ethylene
        cout << "Molar weight of C2H4: " << ethylene.getMolarWeight() << " g/mol" << endl;

        Molecule sulfurDioxide("SO2", 298.15); // Create a Molecule object for sulfur dioxide
        cout << "Molar weight of SO2: " << sulfurDioxide.getMolarWeight() << " g/mol" << endl;
    }
    catch (const invalid_argument& e) {
        cerr << e.what() << endl;
    }
}

void testHitran() {
    try {
        // Create a Hitran object with a test HITRAN file and isotopologues
        Hitran hitran("../hitran/H2O.par", "H2O", {1});

        // Load the HITRAN data within a specific wavelength range (1.0 to 30.0 microns)
        hitran.loadHitranData(1.0, 30.0);

        // Retrieve and print the filtered data
        const auto& filteredData = hitran.getFilteredData();

        // Print filtered data (example: printing the first few rows)
        for (size_t i = 0; i < min<size_t>(filteredData.size(), 5); ++i) {
            for (const string& val : filteredData[i]) {
                cout << val << " ";
            }
            cout << endl;
        }
    }
    catch (const exception& e) {
        cerr << e.what() << endl;
    }
}

int main () {
    testMolecule();
    testHitran();
}
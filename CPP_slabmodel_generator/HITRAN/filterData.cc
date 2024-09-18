#include "Hitran.ih"

void Hitran::filterData(double lowerLam, double upperLam) {
    filteredData.clear();
    
    for (const auto& row : data) {
        int isotopologueID = stoi(row[1]);  // Isotopologue_ID is the second column
        double nu = stod(row[2]);  // nu is in the third column (convert to double)

        // Convert nu (cm^-1) to wavelength (microns)
        double lambda = 1e4 / nu;

        // Check if the row is within the required wavelength range and isotopologue filter
        if (find(isotopologues.begin(), isotopologues.end(), isotopologueID) != isotopologues.end() &&
            lambda >= lowerLam && lambda <= upperLam) {
            filteredData.push_back(row);
        }
    }
}
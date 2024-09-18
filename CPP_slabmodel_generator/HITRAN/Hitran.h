#ifndef HITRAN_H
#define HITRAN_H

#include <string>
#include <vector>

class Hitran {
    std::string filepath;
    std::string molecule;
    std::vector<int> isotopologues;
    std::vector<std::vector<std::string>> data;  // Store the full dataset as strings for simplicity
    std::vector<std::vector<std::string>> filteredData;  // Store filtered data based on isotopologues and wavelength range

public:
    Hitran(const std::string& filepath, const std::string& molecule, const std::vector<int>& isotopologues);

    // Load and process the HITRAN file within wavelength bounds
    void loadHitranData(double lowerLam, double upperLam);

    // Get the filtered data (for now returning a vector of rows, later can add more structure)
    const std::vector<std::vector<std::string>>& getFilteredData() const;

private:
    void readFixedWidthFile();
    void filterData(double lowerLam, double upperLam);
};

#endif

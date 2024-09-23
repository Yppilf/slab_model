#include "SlabModel.ih"

const vector<vector<string>>& SlabModel::readHitranData() {
    // Use the Hitran class to load the data with wavelength filtering
    Hitran hitran(HITRANfile, molecule, isotopologues);
    
    // Load data within the specified wavelength range
    hitran.loadHitranData(lower_wavelength, upper_wavelength);

    const vector<vector<string>>& filteredData = hitran.getFilteredData();
    // Use the filtered data for spectrum generation

    return filteredData;
}
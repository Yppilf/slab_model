#include "SlabModel.ih"

void SlabModel::generateSlab(const string& output_filename, const string& mode) {
    // Step 1: Read HITRAN data, passing the wavelength range
    readHitranData();

    // Step 2: Compute the spectrum
    double R_grid = 1e5;
    computeSpectrum(R_grid);

    // Step 3: Optionally write the data to a file
    if (mode == "both" || mode == "file") {
        writeToFile(output_filename);
    }

    // Return the data if needed
    if (mode == "both" || mode == "return") {
        cout << "Slab model generation complete." << endl;
    }
}
#ifndef SLABMODEL_H
#define SLABMODEL_H

#include "../HITRAN/Hitran.h"
#include <tuple>

// Define a structure to hold line data
struct LineData {
    int i;
    double gu;
    double Eu;
    double A;
    double nu;
    double FNLTE;
    // Additional fields as necessary
};

// Define a structure to hold level data
struct LevelData {
    int i;
    double g;
    double E;
    double pop;
    double ltepop;
    // Additional fields as necessary
};

// Define a structure to hold the slab data
struct Slab {
    std::string directory;
    int model_number;
    double NH;
    double nColl;
    int ne;
    int nHe;
    int nHII;
    int nHI;
    int nH2;
    double dust_to_gas;
    double vturb;
    double Tg;
    double Td;
    int species_index;
    int species_number;
    std::string species_name;
    double abundance;
    double dv;
    std::vector<double> overlapLTE;
    std::vector<double> overlapTauLTE;
    std::vector<double> overlapNLTE;
    std::vector<double> overlapTauNLTE;
    std::vector<double> overlapFreq;
    double overlapR;
    int nlevels;
    std::vector<LevelData> leveldata;
    int nlines;
    std::vector<LineData> linedata;

    void add_model(const Slab& model) {
        // Implementation to add model to slab (if needed)
    }

    void write_to_file(const std::string& filename, const std::string& mode, bool verbose) const {
        // Implement file writing based on mode
    }
};


class SlabModel {
    double column_density;     // Gas column density in molecules per cm²
    double temperature;        // Gas temperature in Kelvin
    double vturb;              // Turbulent velocity (in km/s)
    std::string molecule;      // Molecule name (e.g., "CO2")
    double mol_mass;           // Molecular mass in atomic mass units (amu)
    std::string HITRANfile;    // Path to the HITRAN file
    std::string QTpath;        // Path to the partition sums (QT) directory
    std::vector<int> isotopologues;  // List of isotopologues to include
    double lower_wavelength;   // Lower wavelength limit (in microns)
    double upper_wavelength;   // Upper wavelength limit (in microns)

public:
    // Constructor with the wavelength range
    SlabModel(double column_density, double temperature, double vturb, const string& molecule,
        double mol_mass, const vector<int>& isotopologues, double lower_wavelength, double upper_wavelength);

    // Method to generate a slab model
    void generateSlab(const std::string& output_filename, const std::string& mode = "both");

    // Placeholder for convolving the data
    void convolveData();

private:
    const std::vector<std::vector<std::string>>& readHitranData();
    void computeSpectrum(double R_grid);
    void writeToFile(const std::string& filename);

    std::vector<double> boltzmann_distribution(const std::vector<double>& E, const std::vector<double>& g, double T, double Q);
    double fetch_QT(const std::string& molecule, int isotopologue, double T, const std::string& QTpath);
    std::vector<double> profile_function(const std::vector<double>& nu_grid, double nu, double gamma);
    std::vector<double> planck(double T, const std::vector<double>& nu_grid);
    std::vector<double> growth_function(const std::vector<double>& tau);
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> run_0D_point(double R_grid, 
        const std::vector<double>& nu_0, 
        const std::vector<double>& Aul, 
        const std::vector<double>& E_u, 
        const std::vector<double>& E_l, 
        const std::vector<double>& gu, 
        const std::vector<double>& gl, 
        const std::string& QTpath, 
        const std::vector<int>& isotopolog, 
        const std::vector<double>& nu_grid, 
        const std::string& mode);

    std::vector<double> fetch_QT(const std::string& mol, int iso, double T, const std::string& QTpath, const std::string& T_limit, bool verbose);
    Slab convert_to_slab_format(
        const std::vector<double>& nu_grid,
        const std::vector<double>& I_nu_line,
        const std::vector<double>& tau_ret,
        const std::vector<double>& I_nu,
        const std::vector<double>& tau_lbl,
        double R_grid,
        double Ng,
        double Tg,
        double vturb,
        const std::vector<double>& pop_l,
        const std::vector<double>& pop_u,
        const std::vector<LineData>& mol_data,
        const std::string& molecule,
        const std::string& mode
    );
};


#endif

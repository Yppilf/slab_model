#include "SlabModel.ih"

SlabModel::SlabModel(double column_density, double temperature, double vturb, const string& molecule,
                     double mol_mass, const vector<int>& isotopologues, double lower_wavelength, double upper_wavelength)
    : column_density(column_density), temperature(temperature), vturb(vturb), molecule(molecule),
      mol_mass(mol_mass), isotopologues(isotopologues),
      lower_wavelength(lower_wavelength), upper_wavelength(upper_wavelength) 
{
    HITRANfile = "../hitran/"+molecule+".par";
    QTpath = "./QTpy/";
}
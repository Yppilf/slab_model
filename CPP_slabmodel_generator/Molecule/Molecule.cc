#include "Molecule.ih"

Molecule::Molecule(const std::string& name, double temperature) 
    : name(name), temperature(temperature) {
    molarWeight = calculateMolarWeight();
}
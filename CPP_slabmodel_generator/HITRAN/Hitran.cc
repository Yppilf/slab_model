#include "Hitran.ih"

Hitran::Hitran(const string& filepath, const string& molecule, const vector<int>& isotopologues)
    : filepath(filepath), molecule(molecule), isotopologues(isotopologues) {}
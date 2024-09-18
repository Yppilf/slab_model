#include "Molecule.ih"

double Molecule::getAtomicWeight(const string& element) {
    // Periodic table with common atomic weights
    static const map<string, double> atomicWeights = {
        {"H", 1.008}, {"He", 4.0026}, {"Li", 6.94}, {"Be", 9.0122}, {"B", 10.81},
        {"C", 12.011}, {"N", 14.007}, {"O", 15.999}, {"F", 18.998}, {"Ne", 20.180},
        {"Na", 22.990}, {"Mg", 24.305}, {"Al", 26.982}, {"Si", 28.085}, {"P", 30.974},
        {"S", 32.06}, {"Cl", 35.45}, {"K", 39.098}, {"Ca", 40.078},
        // Add more elements as needed
    };

    auto it = atomicWeights.find(element);
    if (it != atomicWeights.end()) {
        return it->second;
    } else {
        throw invalid_argument("Unknown element: " + element);
    }
}
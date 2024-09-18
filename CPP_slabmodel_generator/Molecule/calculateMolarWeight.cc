#include "Molecule.ih"

// Private method to calculate molar weight
double Molecule::calculateMolarWeight() {
    map<std::string, int> elements = extractElements(name);
    double totalWeight = 0.0;

    for (const auto& element : elements) {
        double atomicWeight = getAtomicWeight(element.first);
        totalWeight += atomicWeight * element.second;
    }

    return totalWeight;
}
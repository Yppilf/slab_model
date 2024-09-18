#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>
#include <map>

class Molecule {
    std::string name;
    double temperature;
    double molarWeight;

public:
    Molecule(const std::string& name, double temperature);

    double getMolarWeight() const;      // Getter function

private:
    // Private method to calculate molar weight internally (replaces external library)
    double calculateMolarWeight();
    double getAtomicWeight(const std::string& element);
    std::map<std::string, int> extractElements(const std::string& molecule) const;
};

#endif

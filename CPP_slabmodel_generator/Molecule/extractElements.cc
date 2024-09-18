#include "Molecule.ih"

map<string, int> Molecule::extractElements(const string& molecule) const {
    map<string, int> elementCounts;
    regex pattern("([A-Z][a-z]*)(\\d*)");
    auto matches_begin = sregex_iterator(molecule.begin(), molecule.end(), pattern);
    auto matches_end = sregex_iterator();

    for (sregex_iterator i = matches_begin; i != matches_end; ++i) {
        smatch match = *i;
        string element = match.str(1);
        string countStr = match.str(2);
        int count = countStr.empty() ? 1 : stoi(countStr);

        elementCounts[element] += count;
    }

    return elementCounts;
}
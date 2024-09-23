#include "SlabModel.ih"

std::vector<double> SlabModel::boltzmann_distribution(const std::vector<double>& E, const std::vector<double>& g, double T, double Q) {
    const double kb = 1.380649e-23; // Boltzmann constant in J/K
    std::vector<double> distribution;
    distribution.reserve(E.size());
    for (size_t i = 0; i < E.size(); ++i) {
        double exp_term = std::exp(-E[i] / (kb * T));
        distribution.push_back(g[i] * exp_term / Q);
    }
    return distribution;
}
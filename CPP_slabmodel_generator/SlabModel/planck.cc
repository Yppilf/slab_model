#include "SlabModel.ih"

std::vector<double> SlabModel::planck(double T, const std::vector<double>& nu_grid) {
    const double kb = 1.380649e-23; // Boltzmann constant in J/K
    std::vector<double> intensity;
    intensity.reserve(nu_grid.size());
    for (const auto& nu : nu_grid) {
        double hnu = 6.62607015e-34 * nu; // Planck's constant * frequency
        double exponent = hnu / (kb * T);
        intensity.push_back(2.0 * hnu * nu * nu / (std::pow(3e10, 2)) * 1.0 / (std::exp(exponent) - 1.0));
    }
    return intensity;
}
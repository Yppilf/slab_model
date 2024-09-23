#include "SlabModel.ih"

std::vector<double> SlabModel::profile_function(const std::vector<double>& nu_grid, double nu, double gamma) {
    std::vector<double> profile;
    profile.reserve(nu_grid.size());
    for (const auto& nu_grid_point : nu_grid) {
        double x = (nu_grid_point - nu) / gamma;
        profile.push_back(std::exp(-x * x));
    }
    return profile;
}
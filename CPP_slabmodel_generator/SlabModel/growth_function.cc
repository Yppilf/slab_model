#include "SlabModel.ih"

std::vector<double> SlabModel::growth_function(const std::vector<double>& tau) {
    std::vector<double> growth;
    growth.reserve(tau.size());
    for (const auto& tau_val : tau) {
        growth.push_back(1.0 - std::exp(-tau_val));
    }
    return growth;
}
#include "SlabModel.ih"

const double kb = 1.380649e-23; // Boltzmann constant in J/K
const double amu = 1.66053906660e-27; // Atomic mass unit in kg
const double c = 2.99792458e10; // Speed of light in cm/s

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>>
SlabModel::run_0D_point(double R_grid, 
                        const std::vector<double>& nu_0, 
                        const std::vector<double>& Aul, 
                        const std::vector<double>& E_u, 
                        const std::vector<double>& E_l, 
                        const std::vector<double>& gu, 
                        const std::vector<double>& gl, 
                        const std::string& QTpath, 
                        const std::vector<int>& isotopolog, 
                        const std::vector<double>& nu_grid, 
                        const std::string& mode) {
    
    size_t N_lines = nu_0.size();
    size_t N_nu_lines = nu_grid.size();
    std::vector<double> vg_grid(N_lines);
    std::vector<double> pop_u(N_lines);
    std::vector<double> pop_l(N_lines);
    std::vector<double> Tauuu(N_lines);
    std::vector<double> sigma(N_lines);
    std::vector<double> tau_ret(N_nu_lines, 0.0);
    std::vector<double> I_nu_line(N_nu_lines, 0.0);
    std::vector<double> I_nu(N_lines, 0.0);

    // Calculate the turbulent velocity grid
    for (size_t i = 0; i < N_lines; ++i) {
        vg_grid[i] = std::sqrt(vturb * vturb + 2.0 * kb * temperature / (mol_mass * amu));
    }

    // Level population calculation
    double Q = fetch_QT(molecule, isotopolog[0], temperature, QTpath);
    auto pop_u_vec = boltzmann_distribution(E_u, gu, temperature, Q);
    auto pop_l_vec = boltzmann_distribution(E_l, gl, temperature, Q);

    for (size_t i = 0; i < N_lines; ++i) {
        pop_u[i] = pop_u_vec[i];
        pop_l[i] = pop_l_vec[i];
    }

    // Line optical depths
    for (size_t i = 0; i < N_lines; ++i) {
        Tauuu[i] = (Aul[i] * std::pow(c / nu_0[i], 3) / (8.0 * std::pow(M_PI, 1.5) * vg_grid[i]) * column_density * (gu[i] / gl[i] * pop_l[i] - pop_u[i]));
        sigma[i] = vg_grid[i] * nu_0[i] / c;
    }

    // RT calculation
    if (mode == "overlap" || mode == "both") {
        for (size_t inu = 0; inu < N_lines; ++inu) {
            std::vector<double> mask(N_nu_lines, 0.0);
            for (size_t j = 0; j < N_nu_lines; ++j) {
                mask[j] = std::abs(nu_grid[j] - nu_0[inu]) < (20.0 / 2.0 * sigma[inu]);
            }
            auto prof = profile_function(nu_grid, nu_0[inu], vg_grid[inu] * nu_0[inu] / c);
            for (size_t j = 0; j < N_nu_lines; ++j) {
                if (mask[j]) {
                    tau_ret[j] += Tauuu[inu] * prof[j];
                }
            }
        }
        auto planck_vals = planck(temperature, nu_grid);
        for (size_t i = 0; i < N_nu_lines; ++i) {
            I_nu_line[i] = planck_vals[i] * (1.0 - std::exp(-tau_ret[i]));
        }
    }

    if (mode == "line_by_line" || mode == "both") {
        auto planck_vals = planck(temperature, nu_0);
        auto growth_vals = growth_function(Tauuu);
        for (size_t i = 0; i < N_lines; ++i) {
            I_nu[i] = planck_vals[i] * growth_vals[i] * vg_grid[i] * nu_0[i] / c;
        }
    }

    return {I_nu_line, tau_ret, I_nu, Tauuu, pop_l, pop_u};
}
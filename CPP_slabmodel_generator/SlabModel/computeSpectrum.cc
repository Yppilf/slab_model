#include "SlabModel.ih"

// Constants
const double c = 2.99792458e10; // Speed of light in cm/s
const double um = 1e-6; // Micrometers to meters
const double percm2_to_perm2 = 1e4; // Conversion factor from cm² to m²
const double kmps = 1e5; // Conversion factor from km/s to cm/s

vector<double> stringToDoubles(const vector<std::string>& strVec) {
    std::vector<double> dblVec;
    for (const auto& str : strVec) {
        dblVec.push_back(stod(str));
    }
    return dblVec;
}

// Helper function to convert string to double
double stringToDouble(const string& s) {
    try {
        return stod(s);
    } catch (const invalid_argument&) {
        throw runtime_error("Invalid data in mol_data: unable to convert string to double.");
    }
}

// Conversion function
std::vector<LineData> convertMolData(const std::vector<std::vector<std::string>>& mol_data) {
    std::vector<LineData> line_data;

    for (const auto& row : mol_data) {
        if (row.size() < 5) {
            throw std::runtime_error("mol_data row does not have enough fields.");
        }

        LineData line;
        line.i = std::stoi(row[0]);  // Convert index
        line.gu = stringToDouble(row[1]);  // Convert g_u (statistical weight of upper level)
        line.Eu = stringToDouble(row[2]);  // Convert E_u (energy of upper state)
        line.A = stringToDouble(row[3]);  // Convert A (Einstein A coefficient)
        line.nu = stringToDouble(row[4]);  // Convert nu (frequency)

        line_data.push_back(line);  // Add this line to the vector
    }

    return line_data;
}


void SlabModel::computeSpectrum(double R_grid) {
    const std::vector<std::vector<std::string>>& mol_data = readHitranData();
    size_t N_lines = mol_data.size();

    if (N_lines == 0) {
        throw std::runtime_error("No HITRAN data available for processing.");
    }

    // Define wavelength grid
    double wave_spec_low = lower_wavelength;
    double wave_spec_high = upper_wavelength;

    size_t N_nu_lines = 1 + static_cast<size_t>(std::log10(wave_spec_high / wave_spec_low) / std::log10(1.0 + 1.0 / R_grid));
    std::vector<double> nu_grid(N_nu_lines);
    double dnu_ovlp = (std::log10(c / (wave_spec_low * um)) - std::log10(c / (wave_spec_high * um))) / N_nu_lines;

    for (size_t i = 0; i < N_nu_lines; ++i) {
        nu_grid[N_nu_lines - 1 - i] = std::pow(10, std::log10(c / (wave_spec_high * um)) + i * dnu_ovlp);
    }

    // Convert nu_grid to wavelength grid
    std::vector<double> W_grid(N_nu_lines);
    for (size_t i = 0; i < N_nu_lines; ++i) {
        W_grid[i] = c / nu_grid[i] / um;
    }

    // Prepare arrays from HITRAN data
    std::vector<double> lambda_0(N_lines);
    std::vector<double> nu_0(N_lines);
    std::vector<double> Aul(N_lines);
    std::vector<double> gu(N_lines);
    std::vector<double> gl(N_lines);
    std::vector<double> E_u(N_lines);
    std::vector<double> E_l(N_lines);

    for (size_t i = 0; i < N_lines; ++i) {
        const auto& row = mol_data[i];
        if (row.size() < 19) {
            throw std::runtime_error("HITRAN data row does not contain enough columns.");
        }

        lambda_0[i] = std::stod(row[2]); // Assuming 3rd column is wavelength
        nu_0[i] = c / (lambda_0[i] * um);
        Aul[i] = std::stod(row[3]); // Assuming 4th column is Aul
        gu[i] = std::stod(row[14]); // Assuming 15th column is gu
        gl[i] = std::stod(row[15]); // Assuming 16th column is gl
        E_u[i] = std::stod(row[8]); // Assuming 9th column is E_u
        E_l[i] = std::stod(row[7]); // Assuming 8th column is E_l
    }

    // Convert units
    column_density *= percm2_to_perm2;
    vturb *= kmps;

    // Step 4: Calculate the spectrum
    std::vector<double> I_nu_line, tau_ret, I_nu, tau_lbl, pop_l, pop_u;

    auto [I_nu_line, tau_ret, I_nu, Tauuu, pop_l, pop_u] = run_0D_point(R_grid, nu_0, Aul, E_u, E_l, gu, gl, QTpath, isotopologues, nu_grid, "both");

    // Step 5: Convert to slab format
    std::vector<LineData> mol_data = convertMolData(mol_data);
    Slab slab = convert_to_slab_format(
        nu_grid, I_nu_line, tau_ret, I_nu, tau_lbl, R_grid, column_density, temperature, vturb, pop_l, pop_u, mol_data, molecule, "both"
    );

    // For demonstration, output a placeholder message
    std::cout << "Spectrum computation and conversion complete." << std::endl;
}

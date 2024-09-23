#include "SlabModel.ih";

std::string trim(const std::string& str) {
    const auto strBegin = str.find_first_not_of(" \t");
    if (strBegin == std::string::npos) return "";
    const auto strEnd = str.find_last_not_of(" \t");
    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

std::unordered_map<std::string, double> read_QT_file(const std::string& file_path) {
    std::unordered_map<std::string, double> QTdict;
    std::ifstream file(file_path);

    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + file_path);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string json_content = buffer.str();
    file.close();

    size_t pos = 0;
    while ((pos = json_content.find('{')) != std::string::npos) {
        size_t end_pos = json_content.find('}', pos);
        if (end_pos == std::string::npos) break;
        std::string json_object = json_content.substr(pos, end_pos - pos + 1);
        json_content.erase(pos, end_pos - pos + 1);

        // Simple parsing assuming JSON is well-formed
        size_t key_start = json_object.find('\"') + 1;
        size_t key_end = json_object.find('\"', key_start);
        std::string key = json_object.substr(key_start, key_end - key_start);

        size_t value_start = json_object.find(':', key_end) + 1;
        size_t value_end = json_object.find(',', value_start);
        if (value_end == std::string::npos) value_end = json_object.find('}', value_start);
        std::string value_str = json_object.substr(value_start, value_end - value_start);
        double value = std::stod(trim(value_str));

        QTdict[key] = value;
    }

    return QTdict;
}

std::vector<double> SlabModel::fetch_QT(const std::string& mol, int iso, double T, const std::string& QTpath, const std::string& T_limit, bool verbose) {
    // Example QT_mol_id, QT_niso, and QT_Tmax as described
    std::unordered_map<std::string, int> QT_mol_id = {
        {"H2O", 1}, {"CO2", 2}, {"O3", 3}, {"N2O", 4}, {"CO", 5}, {"CH4", 6},
        {"O2", 7}, {"NO", 8}, {"SO2", 9}, {"NO2", 10}, {"NH3", 11}, {"HNO3", 12},
        {"OH", 13}, {"HF", 14}, {"HCl", 15}, {"HBr", 16}, {"HI", 17}, {"ClO", 18},
        {"OCS", 19}, {"H2CO", 20}, {"HOCl", 21}, {"N2", 22}, {"HCN", 23}, {"CH3Cl", 24},
        {"H2O2", 25}, {"C2H2", 26}, {"C2H6", 27}, {"PH3", 28}, {"COF2", 29}, {"SF6", 30},
        {"H2S", 31}, {"HCOOH", 32}, {"HO2", 33}, {"ClONO2", 35}, {"NO+", 36}, {"HOBr", 37},
        {"C2H4", 38}, {"CH3OH", 39}, {"CH3Br", 40}, {"CH3CN", 41}, {"CF4", 42}, {"C4H2", 43},
        {"HC3N", 44}, {"H2", 45}, {"CS", 46}, {"SO3", 47}, {"C2N2", 48}, {"COCl2", 49},
        {"SO", 50}, {"CH3F", 51}, {"GeH4", 52}, {"CS2", 53}, {"CH3I", 54}, {"NF3", 55},
        {"C3H4", 56}, {"CH3", 57}
    };

    std::vector<int> QT_niso = {0, 9, 13, 18, 5, 9, 4, 6, 3, 4, 2, 2, 2, 3, 2, 4, 4, 2, 2, 6, 3, 2, 3, 3, 2, 1, 3, 3, 1, 2, 1, 3, 1, 1, 1, 2, 1, 2, 3, 1, 2, 4, 1, 1, 6, 2, 4, 1, 2, 2, 3, 1, 5, 4, 2, 1, 1, 1};

    std::vector<double> QT_Tmax = {0, 5000., 5000., 5000., 5000., 5000., 5000., 6000., 6000., 6000.,
        5000., 5000., 3500., 3500., 3500., 3500., 5000., 3500., 5000., 5000., 3500., 5000., 5000.,
        1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.,
        5000., 5000., 5000., 5000., 5000., 9000., 9000., 9000., 9000., 9000., 9000., 9000., 9000., 9000.,
        2500., 2500., 2500., 2500., 7500., 7500., 7500., 7500., 7500., 7500., 5000., 5000., 5000., 5000., 5000., 5000., 5000.,
        1000., 1000., 6000., 6000., 3500., 3500., 9000., 5000., 5000., 6000., 6000., 6000., 6000., 6000., 6000., 6000.,
        6000., 6000., 5000., 5000., 5000., 5000., 5000., 5000., 3500., 5000., 5000., 5000., 5000., 5000., 6000., 6000.,
        5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000., 5000.};

    // Map molecule name to index
    auto mol_it = QT_mol_id.find(mol);
    if (mol_it == QT_mol_id.end()) {
        throw std::invalid_argument("Molecule not found in QT_mol_id.");
    }
    int mol_id = mol_it->second;

    if (iso <= 0 || iso > QT_niso[mol_id]) {
        throw std::out_of_range("Isotopologue number is out of range.");
    }

    std::string file_path = QTpath + std::to_string(mol_id) + "_" + std::to_string(iso) + ".json";
    auto QTdict = read_QT_file(file_path);

    std::vector<double> QT;

    if (T < 1.0 || T > QT_Tmax[mol_id]) {
        if (T_limit == "warn") {
            if (T < 1.0) {
                if (verbose) {
                    std::cerr << "WARNING: " << T << " falls below the temperature range 1 to " << QT_Tmax[mol_id] << "\n";
                    std::cerr << "Resetting T = 1\n";
                }
                T = 1.0;
            } else {
                if (verbose) {
                    std::cerr << "WARNING: " << T << " falls above the temperature range 1 to " << QT_Tmax[mol_id] << "\n";
                    std::cerr << "Resetting T = " << QT_Tmax[mol_id] << "\n";
                }
                T = QT_Tmax[mol_id];
            }
        } else {
            throw std::out_of_range("Temperature is out of range.");
        }
    }

    if (QTdict.find(std::to_string(static_cast<int>(T))) == QTdict.end() ||
        QTdict.find(std::to_string(static_cast<int>(T) + 1)) == QTdict.end()) {
        throw std::runtime_error("Temperature not found in QT file.");
    }

    double Q1 = QTdict[std::to_string(static_cast<int>(T))];
    double Q2 = QTdict[std::to_string(static_cast<int>(T) + 1)];
    double QT_value = Q1 + (Q2 - Q1) * (T - static_cast<int>(T));

    QT.push_back(QT_value);

    return QT;
}

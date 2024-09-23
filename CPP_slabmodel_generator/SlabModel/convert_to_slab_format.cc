#include "SlabModel.ih"

Slab SlabModel::convert_to_slab_format(
    const std::vector<double>& nu_grid,
    const std::vector<double>& I_nu_line,
    const std::vector<double>& tau_ret,
    const std::vector<double>& I_nu,
    const std::vector<double>& tau_lbl,
    double R_grid,
    double Ng,
    double Tg,
    double vturb,
    const std::vector<double>& pop_l,
    const std::vector<double>& pop_u,
    const std::vector<LineData>& mol_data,
    const std::string& molecule,
    const std::string& mode
) {
    Slab slab;
    slab.directory = "";
    slab.model_number = 1;
    slab.NH = Ng / 1e-4; // Adjust based on conversion factor
    slab.nColl = Ng / 1e-4; // Adjust based on conversion factor
    slab.ne = 1;
    slab.nHe = 1;
    slab.nHII = 1;
    slab.nHI = 1;
    slab.nH2 = 1;
    slab.dust_to_gas = 1;
    slab.vturb = vturb * 1e-1; // Adjust km/s to m/s
    slab.Tg = Tg;
    slab.Td = 0;
    slab.species_index = 1;
    slab.species_number = 1;
    slab.species_name = molecule;
    slab.abundance = 1;
    slab.dv = vturb * 1e3; // km/s to m/s

    if (mode == "overlap" || mode == "both") {
        slab.overlapLTE = I_nu_line;
        std::reverse(slab.overlapLTE.begin(), slab.overlapLTE.end());
        slab.overlapTauLTE = tau_ret;
        std::reverse(slab.overlapTauLTE.begin(), slab.overlapTauLTE.end());
        slab.overlapNLTE.assign(I_nu_line.size(), 0.0);
        slab.overlapTauNLTE.assign(I_nu_line.size(), 0.0);
        slab.overlapFreq.assign(nu_grid.begin(), nu_grid.end());
        std::reverse(slab.overlapFreq.begin(), slab.overlapFreq.end());
        slab.overlapR = R_grid;
    }

    if (mode == "line_by_line" || mode == "both") {
        slab.nlevels = 2 * mol_data.size();
        for (size_t i = 0; i < mol_data.size(); ++i) {
            LevelData level_data;
            level_data.i = 2 * i;
            level_data.g = mol_data[i].gu;
            level_data.E = mol_data[i].Eu;
            level_data.pop = pop_u[i];
            level_data.ltepop = 0.0; // Adjust if necessary
            slab.leveldata.push_back(level_data);

            level_data.i = 2 * i + 1;
            level_data.g = mol_data[i].gu;
            level_data.E = mol_data[i].Eu;
            level_data.pop = pop_l[i];
            level_data.ltepop = 0.0; // Adjust if necessary
            slab.leveldata.push_back(level_data);
        }

        slab.nlines = mol_data.size();
        for (size_t i = 0; i < mol_data.size(); ++i) {
            LineData line_data;
            line_data.i = static_cast<int>(i);
            line_data.gu = mol_data[i].gu;
            line_data.Eu = mol_data[i].Eu;
            line_data.A = mol_data[i].A;
            line_data.nu = mol_data[i].nu * 1e-9; // Convert to GHz
            line_data.FNLTE = I_nu[i] * 1e3; // Convert to appropriate units
            slab.linedata.push_back(line_data);
        }
    }

    return slab;
}

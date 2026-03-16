/*!
 * \file CFluidModel_tests.cpp
 * \brief Unit tests for the fluid model classes.
 * \author E.Bunschoten
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "catch.hpp"
#include <sstream>
#include "../../../SU2_CFD/include/fluid/CFluidModel.hpp"
#include "../../../SU2_CFD/include/fluid/CIdealGas.hpp"
#include "../../../SU2_CFD/include/fluid/CDataDrivenFluid.hpp"

std::array<su2double, 4> FluidModelChecks(CFluidModel* fluid_model, const su2double val_p, const su2double val_T, const su2double delta_rho = 1e-2, const su2double delta_e = 100) {
  /*--- Check consistency of reverse look-up ---*/
  {
    fluid_model->SetTDState_PT(val_p, val_T);

    const su2double val_rho_fluidmodel = fluid_model->GetDensity();
    const su2double val_e_fluidmodel = fluid_model->GetStaticEnergy();

    fluid_model->SetTDState_rhoe(val_rho_fluidmodel, val_e_fluidmodel);
    CHECK(Approx(fluid_model->GetPressure()) == val_p);
    CHECK(Approx(fluid_model->GetTemperature()) == val_T);
  }
  /*--- Check internal consistency between primary and derived fluid properties ---*/
  fluid_model->SetTDState_PT(val_p, val_T);
  const su2double val_rho = fluid_model->GetDensity();
  const su2double val_e = fluid_model->GetStaticEnergy();

  const su2double dTdrho_e = fluid_model->GetdTdrho_e();
  const su2double dPdrho_e = fluid_model->GetdPdrho_e();
  const su2double dTde_rho = fluid_model->GetdTde_rho();
  const su2double dPde_rho = fluid_model->GetdPde_rho();
  
  fluid_model->SetTDState_rhoe(val_rho + delta_rho, val_e);
  su2double T_plus = fluid_model->GetTemperature();
  su2double p_plus = fluid_model->GetPressure();

  fluid_model->SetTDState_rhoe(val_rho - delta_rho, val_e);
  su2double T_minus = fluid_model->GetTemperature();
  su2double p_minus = fluid_model->GetPressure();
  const su2double dTdrho_e_FD = (T_plus - T_minus) / (2 * delta_rho);
  const su2double dPdrho_e_FD = (p_plus - p_minus) / (2 * delta_rho);

  CHECK(dTdrho_e == Approx(dTdrho_e_FD));
  const su2double err_rel_dTdrho_e = std::abs((dTdrho_e - dTdrho_e_FD) / dTdrho_e);
  const su2double err_rel_dPdrho_e = std::abs((dPdrho_e - dPdrho_e_FD) / dPdrho_e);
  CHECK(err_rel_dTdrho_e < 1e-3);
  CHECK(err_rel_dPdrho_e < 1e-3);
  CHECK(dPdrho_e == Approx(dPdrho_e_FD));
  
  
  fluid_model->SetTDState_rhoe(val_rho, val_e + delta_e);
  T_plus = fluid_model->GetTemperature();
  p_plus = fluid_model->GetPressure();

  fluid_model->SetTDState_rhoe(val_rho, val_e - delta_e);
  T_minus = fluid_model->GetTemperature();
  p_minus = fluid_model->GetPressure();
  const su2double dTde_rho_FD = (T_plus - T_minus) / (2 * delta_e);
  const su2double dPde_rho_FD = (p_plus - p_minus) / (2 * delta_e);
  const su2double err_rel_dTde_rho = std::abs((dTde_rho - dTde_rho_FD) / dTde_rho); 
  const su2double err_rel_dPde_rho = std::abs((dPde_rho - dPde_rho_FD) / dPde_rho);
  CHECK(err_rel_dTde_rho < 1e-3);
  CHECK(err_rel_dPde_rho < 1e-3);
  CHECK(dTde_rho == Approx(dTde_rho_FD));
  CHECK(dPde_rho == Approx(dPde_rho_FD));
  

  return {err_rel_dTdrho_e, err_rel_dPdrho_e, err_rel_dTde_rho, err_rel_dPde_rho};
}

std::array<su2double, 4> FluidModelChecks_Prho(CFluidModel* fluid_model, const su2double val_p_ass, const su2double val_rho_ass, const su2double delta_rho = 1e-2, const su2double delta_e = 100) {
  /*--- Check consistency of reverse look-up ---*/
  {
    fluid_model->SetTDState_Prho(val_p_ass, val_rho_ass);

    const su2double val_rho_fluidmodel = fluid_model->GetDensity();
    const su2double val_e_fluidmodel = fluid_model->GetStaticEnergy();

    fluid_model->SetTDState_rhoe(val_rho_fluidmodel, val_e_fluidmodel);
    CHECK(Approx(fluid_model->GetPressure()) == val_p_ass);
    CHECK(Approx(fluid_model->GetDensity()) == val_rho_ass);
  }
  /*--- Check internal consistency between primary and derived fluid properties ---*/
  fluid_model->SetTDState_Prho(val_p_ass, val_rho_ass);
  const su2double val_rho = fluid_model->GetDensity();
  const su2double val_e = fluid_model->GetStaticEnergy();

  const su2double dTdrho_e = fluid_model->GetdTdrho_e();
  const su2double dPdrho_e = fluid_model->GetdPdrho_e();
  const su2double dTde_rho = fluid_model->GetdTde_rho();
  const su2double dPde_rho = fluid_model->GetdPde_rho();
  
  fluid_model->SetTDState_rhoe(val_rho + delta_rho, val_e);
  su2double T_plus = fluid_model->GetTemperature();
  su2double p_plus = fluid_model->GetPressure();

  fluid_model->SetTDState_rhoe(val_rho - delta_rho, val_e);
  su2double T_minus = fluid_model->GetTemperature();
  su2double p_minus = fluid_model->GetPressure();
  const su2double dTdrho_e_FD = (T_plus - T_minus) / (2 * delta_rho);
  const su2double dPdrho_e_FD = (p_plus - p_minus) / (2 * delta_rho);

  CHECK(dTdrho_e == Approx(dTdrho_e_FD));
  const su2double err_rel_dTdrho_e = std::abs((dTdrho_e - dTdrho_e_FD) / dTdrho_e);
  const su2double err_rel_dPdrho_e = std::abs((dPdrho_e - dPdrho_e_FD) / dPdrho_e);
  CHECK(err_rel_dTdrho_e < 1e-3);
  CHECK(err_rel_dPdrho_e < 1e-3);
  CHECK(dPdrho_e == Approx(dPdrho_e_FD));
  
  
  fluid_model->SetTDState_rhoe(val_rho, val_e + delta_e);
  T_plus = fluid_model->GetTemperature();
  p_plus = fluid_model->GetPressure();

  fluid_model->SetTDState_rhoe(val_rho, val_e - delta_e);
  T_minus = fluid_model->GetTemperature();
  p_minus = fluid_model->GetPressure();
  const su2double dTde_rho_FD = (T_plus - T_minus) / (2 * delta_e);
  const su2double dPde_rho_FD = (p_plus - p_minus) / (2 * delta_e);
  const su2double err_rel_dTde_rho = std::abs((dTde_rho - dTde_rho_FD) / dTde_rho); 
  const su2double err_rel_dPde_rho = std::abs((dPde_rho - dPde_rho_FD) / dPde_rho);
  CHECK(err_rel_dTde_rho < 1e-3);
  CHECK(err_rel_dPde_rho < 1e-3);
  CHECK(dTde_rho == Approx(dTde_rho_FD));
  CHECK(dPde_rho == Approx(dPde_rho_FD));
  

  return {err_rel_dTdrho_e, err_rel_dPdrho_e, err_rel_dTde_rho, err_rel_dPde_rho};
}

std::array<su2double, 4> FluidModelChecks_Der(CFluidModel* fluid_model, const su2double val_rho_ass, const su2double val_e_ass) {
  /*--- Check the derivatives extrapolation  ---*/

  fluid_model->SetTDState_rhoe(val_rho_ass, val_e_ass);
  const su2double dTdrho_e_val=fluid_model->GetdTdrho_e(); 
  const su2double dPdrho_e_val=fluid_model->GetdPdrho_e();
  const su2double dTde_rho_val=fluid_model->GetdTde_rho(); 
  const su2double dPde_rho_val=fluid_model->GetdPde_rho();

  return {dTdrho_e_val, dPdrho_e_val, dTde_rho_val, dPde_rho_val};
}


static inline std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

std::vector<std::array<double, 3>> load_csv_col3_skip_header(const std::string& path, char sep = ','){
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    in.imbue(std::locale::classic()); // forza '.' come separatore decimale

    std::vector<std::array<double, 3>> out;
    std::string line;

    // skip prima riga (header)
    if (!std::getline(in, line)) return out;

    size_t line_no = 1;
    while (std::getline(in, line)) {
        ++line_no;
        line = trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;

        std::array<double, 3> row{};
        std::stringstream ss(line);
        std::string cell;

        for (int k = 0; k < 3; ++k) {
            if (!std::getline(ss, cell, sep)) {
                throw std::runtime_error("Line " + std::to_string(line_no) +
                                         ": expected 3 columns, got fewer: " + line);
            }
            cell = trim(cell);
            try {
                row[k] = std::stod(cell);
            } catch (...) {
                throw std::runtime_error("Line " + std::to_string(line_no) +
                                         ": cannot parse double from '" + cell + "'");
            }
        }

        // opzionale: se vuoi essere severo e rifiutare colonne extra:
        // if (std::getline(ss, cell, sep)) throw std::runtime_error("Line ... too many columns");

        out.push_back(row);
    }

    return out;
}

TEST_CASE("Test case for ideal gas fluid model") {
  CIdealGas* fluid_model = new CIdealGas(1.4, 287.0);
  auto res1=FluidModelChecks(fluid_model, 101325, 300.0);
  auto res2=FluidModelChecks(fluid_model, 1e6, 600.0);

  delete fluid_model;
}

TEST_CASE("Test case for data-driven fluid model") {
  

  string NameFile="/home/ricca/src_SU2/SU2_dev_NS_CheckConv/UnitTests/SU2_CFD/fluid/Test13.csv";
  string LoadFile="/home/ricca/src_SU2/SU2_dev_NS_CheckConv/UnitTests/SU2_CFD/fluid/TROVA_TET4_Scaled_2PH_Nozzle_Expected_Exp.txt";

  std::stringstream config_options;
  config_options << "SOLVER=EULER" << std::endl;
  //config_options << "KIND_TURB_MODEL=SA" << std::endl;
  //config_options << "SA_OPTIONS= NONE" << std::endl;
  //config_options << "REYNOLDS_NUMBER=1e6" << std::endl;
  config_options << "FLUID_MODEL=DATADRIVEN_FLUID" << std::endl;
  config_options << "USE_PINN=NO" << std::endl;
  config_options << "INTERPOLATION_METHOD=LUT" << std::endl;
  config_options << "FILENAMES_INTERPOLATOR=(/home/ricca/src_SU2/SU2_dev_NS_CheckConv/UnitTests/SU2_CFD/fluid/LUT_2PH_TET4_EXP_Adapt_Ref_Add_Ref_V2.drg)" << std::endl;
  config_options << "DATADRIVEN_NEWTON_ITER=75" << std::endl;
  config_options << "DATADRIVEN_NEWTON_TOL=1E-10" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION=0.4" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_MEDIUM=0.15" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_LOW=0.2" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_ITER_MULT=0.333" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_HIGH_Y=3.25e5" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_MEDIUM_Y=2.75e5" << std::endl;
  config_options << "CONV_NUM_METHOD_FLOW=ROE" << std::endl;
  config_options << "MUSCL_FLOW= YES" << std::endl;

  /*--- Setup ---*/

  CConfig* config = new CConfig(config_options, SU2_COMPONENT::SU2_CFD, false);

  /*--- Define fluid model ---*/
  CDataDrivenFluid* fluid_model = new CDataDrivenFluid(config, false);

  su2double delta_rho = 1e-4, delta_e = 10;

  /*--- Check fluid model consistency for several combinations of pressure-temperature. ---*/

  /*
  std::vector<su2double> val_P_vec{19.31e5, 1.37e5};
  std::vector<su2double> val_T_vec{493.75, 404.62};
  
  std::ofstream out(NameFile, std::ios::out); 
  out << std::setprecision(8);
  out << "drho=" << delta_rho << "[kg/m3]" << "\n";
  out << "de=" << delta_e << "[J/kg]" << "\n";
  out << "P [Pa],T [K],err rel dTdrho_e [%],err rel dPdrho_e [%], err rel dTde_rho [%],err rel dPde_rho [%]\n";
  std::array<su2double, 4> res{};
  for (size_t i=0; i<val_P_vec.size(); i++) {
    res=FluidModelChecks(fluid_model, val_P_vec[i], val_T_vec[i],delta_rho,delta_e);
    out << val_P_vec[i] << ","<< val_T_vec[i] << ","<< res[0]*100.0 << ","<< res[1]*100.0 << ","<< res[2]*100.0 << ","<< res[3]*100.0 << ","<< "\n";

  }

  */

  auto expansion=load_csv_col3_skip_header(LoadFile,',');
  std::vector<double> rho_vec, P_vec;
  rho_vec.reserve(expansion.size());
  P_vec.reserve(expansion.size());

  for (const auto& r : expansion) {
    rho_vec.push_back(r[0]);
    P_vec.push_back(r[2]);
  }

  std::ofstream out(NameFile, std::ios::out); 
  out << std::setprecision(8);
  out << "drho=" << delta_rho << "[kg/m3]" << "\n";
  out << "de=" << delta_e << "[J/kg]" << "\n";
  out << "rho [kg/m3],P [Pa],err rel dTdrho_e [%],err rel dPdrho_e [%], err rel dTde_rho [%],err rel dPde_rho [%]\n";
  std::array<su2double, 4> res{};
  for (size_t i=0; i<P_vec.size(); i++) {
    res=FluidModelChecks_Prho(fluid_model, P_vec[i], rho_vec[i],delta_rho,delta_e);
    out << rho_vec[i] << ","<< P_vec[i] << ","<< res[0]*100.0 << ","<< res[1]*100.0 << ","<< res[2]*100.0 << ","<< res[3]*100.0 << ","<< "\n";
  }
  delete config;
  delete fluid_model;
}

TEST_CASE("Test case for first derivatives") {
  

  string NameFile="/home/ricca/src_SU2/SU2_dev_NS_CheckConv/UnitTests/SU2_CFD/fluid/Test10_Derivatives.csv";
  string LoadFile="/home/ricca/src_SU2/SU2_dev_NS_CheckConv/UnitTests/SU2_CFD/fluid/TROVA_TET4_Scaled_2PH_Nozzle_Expected_Exp.txt";

  std::stringstream config_options;
  config_options << "SOLVER=EULER" << std::endl;
  //config_options << "KIND_TURB_MODEL=SA" << std::endl;
  //config_options << "SA_OPTIONS= NONE" << std::endl;
  //config_options << "REYNOLDS_NUMBER=1e6" << std::endl;
  config_options << "FLUID_MODEL=DATADRIVEN_FLUID" << std::endl;
  config_options << "USE_PINN=NO" << std::endl;
  config_options << "INTERPOLATION_METHOD=LUT" << std::endl;
  config_options << "FILENAMES_INTERPOLATOR=(/home/ricca/src_SU2/SU2_dev_NS_CheckConv/UnitTests/SU2_CFD/fluid/LUT_2PH_TET4_EXP_Adapt_Ref_Add_Ref_V2.drg)" << std::endl;
  config_options << "DATADRIVEN_NEWTON_ITER=75" << std::endl;
  config_options << "DATADRIVEN_NEWTON_TOL=1E-10" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION=0.4" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_MEDIUM=0.15" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_LOW=0.2" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_ITER_MULT=0.333" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_HIGH_Y=3.25e5" << std::endl;
  config_options << "DATADRIVEN_NEWTON_EXTRA_RELAXATION_MEDIUM_Y=2.75e5" << std::endl;
  config_options << "CONV_NUM_METHOD_FLOW=ROE" << std::endl;
  config_options << "MUSCL_FLOW= YES" << std::endl;

  /*--- Setup ---*/

  CConfig* config = new CConfig(config_options, SU2_COMPONENT::SU2_CFD, false);

  /*--- Define fluid model ---*/
  CDataDrivenFluid* fluid_model = new CDataDrivenFluid(config, false);

  su2double delta_rho = 1e-4, delta_e = 10;

  /*--- Check fluid model consistency for several combinations of pressure-temperature. ---*/

  /*
  std::vector<su2double> val_P_vec{19.31e5, 1.37e5};
  std::vector<su2double> val_T_vec{493.75, 404.62};
  
  std::ofstream out(NameFile, std::ios::out); 
  out << std::setprecision(8);
  out << "drho=" << delta_rho << "[kg/m3]" << "\n";
  out << "de=" << delta_e << "[J/kg]" << "\n";
  out << "P [Pa],T [K],err rel dTdrho_e [%],err rel dPdrho_e [%], err rel dTde_rho [%],err rel dPde_rho [%]\n";
  std::array<su2double, 4> res{};
  for (size_t i=0; i<val_P_vec.size(); i++) {
    res=FluidModelChecks(fluid_model, val_P_vec[i], val_T_vec[i],delta_rho,delta_e);
    out << val_P_vec[i] << ","<< val_T_vec[i] << ","<< res[0]*100.0 << ","<< res[1]*100.0 << ","<< res[2]*100.0 << ","<< res[3]*100.0 << ","<< "\n";

  }

  */

  auto expansion=load_csv_col3_skip_header(LoadFile,',');
  std::vector<double> rho_vec, e_vec;
  rho_vec.reserve(expansion.size());
  e_vec.reserve(expansion.size());

  for (const auto& r : expansion) {
    rho_vec.push_back(r[0]);
    e_vec.push_back(r[1]);
  }

  std::ofstream out(NameFile, std::ios::out); 
  out << std::setprecision(8);
  out << "rho [kg/m3],e [J/kg],dTdrho_e [1/K],dPdrho_e [J/kg], dTde_rho [kgK/J],dPde_rho [kg/m3]\n";
  std::array<su2double, 4> res{};
  for (size_t i=0; i<rho_vec.size(); i++) {
    res=FluidModelChecks_Der(fluid_model, rho_vec[i], e_vec[i]);
    out << rho_vec[i] << ","<< e_vec[i] << ","<< res[0] << ","<< res[1] << ","<< res[2] << ","<< res[3] << ","<< "\n";
  }
  delete config;
  delete fluid_model;
}

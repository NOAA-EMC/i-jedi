/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/VariableChange/VaderCookbook.h"

#include <map>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "vader/vader.h"

namespace detail {
static std::map<std::string, std::vector<std::string>> cookbook() {
  return {
      // ps: from delp
      {"air_pressure_at_surface",    {"SurfaceAirPressure_A"}},
      // P: from delp, from ps (and ak/bk)
      {"air_pressure_levels",        {"AirPressureAtInterface_B",
                                      "AirPressureAtInterface_A"}},
      // p: from pe
      {"air_pressure",               {"AirPressure_A"}},
      // delp: from p
      {"air_pressure_thickness",     {"AirPressureThickness_A"}},
      // ln(p): from p
      {"ln_air_pressure",            {"LnAirPressure_A"}},
      // ln(pe): from pe
      {"ln_air_pressure_at_interface", {"LnAirPressureAtInterface_A"}},
      // p^kappa: from pe and ln(p)
      {"air_pressure_to_kappa",        {"AirPressureToKappa_A"}},
      // rh: from t q and qsat
      {"relative_humidity",            {"RelativeHumidity_A"}},
      // r: from q
      {"water_vapor_mixing_ratio_wrt_dry_air", {"WaterVaporMixingRatioWrtDryAir_C"}},
      // es: from t
      {"svp",                                  {"SaturationVaporPressure_A"}},
      // qsat: from t, and es
      {"saturation_water_vapor_mixing_ratio_wrt_moist_air", {"SaturationSpecificHumidity_A"}},
      // phi_e from phi tv ln(p) and ln(pe)
      {"geopotential_levels",                  {"GeopotentialLevels_A"}},
      // z: from phi
      {"geopotential_height",                  {"GeopotentialHeight_A", "GeopotentialHeight_B"}},
      // z_surf: from phi_surf
      {"geopotential_height_at_surface",       {"GeopotentialHeightAtSurface_A"}},
      // ze: from phi_e
      {"geopotential_height_levels",           {"GeopotentialHeightLevels_A",
                    "GeopotentialHeightLevels_B"}},
      // h_amsl: from z_surf
      {"height_above_mean_sea_level_at_surface",  {"HeightAboveMeanSeaLevelAtSurface_A"}},
      // pt: from t and pkz, from t and ps
      {"air_potential_temperature",               {"AirPotentialTemperature_B",
                                                   "AirPotentialTemperature_A"}},
      // tv: from t and q
      {"virtual_temperature",                     {"AirVirtualTemperature_A"}},
      {"eastward_wind_at_surface_adjacent_layer", {"EastwardWindAtSurfaceAdjacentLayer_A"}},
      {"northward_wind_at_surface_adjacent_layer", {"NorthwardWindAtSurfaceAdjacentLayer_A"}},
      {"wind_speed_at_surface_adjacent_layer",     {"WindSpeedAtSurfaceAdjacentLayer_A"}},
      {"eastward_wind_at_10m",                    {"EastwardWindAt10m_A"}},
      {"northward_wind_at_10m",                   {"NorthwardWindAt10m_A"}},
      {"wind_speed_at_10m",                       {"WindSpeedAt10m_A", "WindSpeedAt10m_B",
                       "WindSpeedAt10m_C"}},
      {"wind_to_direction_at_10m",                {"WindToDirectionAt10m_A"}},
      {"ratio_of_wind_at_10m_to_wind_at_surface_adjacent_layer",
       {"RatioOfWind10mToSurfaceAdjacent_A"}},
      // t: from tv and q; needed when tv is the control variable
      {"air_temperature",                         {"AirTemperature_B"}},
      // o3mr: from o3 mole fraction; needed when the control variable is ppmv
      {"ozone_mass_mixing_ratio",                 {"OzoneMassMixingRatio_A"}},
      // o3 mole fraction: from o3mr; lets the trajectory hold a variable the state does not
      {"mole_fraction_of_ozone_in_air",           {"MoleFractionOfOzoneInAir_A"}},
      // sst: from t
      {"sea_surface_temperature",                 {"SeaSurfaceTemperature_A"}},
      // t_insitu: NL from thetao (sea_water_potential_temperature), salinity, depth, lat/lon
      //           (SeaWaterTemperature_A); TL/AD via linearized TEOS-10 Jacobian
      //           (SeaWaterTemperature_B). Both are needed so the variable change
      //           works in the variational inner loop as well as the NL observer.
      // TODO(guillaumevernieres): re-enable once the vader insitu (GSW ocean) recipe PR is merged.
      // {"sea_water_temperature",                  {"SeaWaterTemperature_A",
      //                                             "SeaWaterTemperature_B"}},
  };
}
}  // namespace detail

namespace ijedi {

std::map<std::string, std::vector<std::string>> vaderDefaultCookbook() {
  return detail::cookbook();
}

oops::Variables requiredVariables(vader::Vader & vader, const oops::Variables & stateVars,
                                  const oops::Variables & incVars) {
  oops::Variables invalidIncVars = incVars;
  invalidIncVars -= stateVars;
  ASSERT_MSG(invalidIncVars.size() == 0,
             "Vader::requiredVariables: incVars contains variables not present in stateVars");

  oops::Variables requiredVars = stateVars;
  for (const auto & candidateVar : stateVars) {
    if (incVars.has(candidateVar)) continue;
    oops::Variables trialVars = requiredVars;
    trialVars -= candidateVar;
    oops::Variables ingredientVars = trialVars;
    oops::Variables neededVars = stateVars;
    vader::Vader::vaderPlanType plan;
    vader.planVariables(ingredientVars, neededVars, plan);
    if (neededVars.size() == 0) requiredVars = trialVars;
  }
  return requiredVars;
}

}  // namespace ijedi

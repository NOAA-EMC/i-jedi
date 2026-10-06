/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "vader/RecipeBase.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

class SeaIceThicknessParameters : public vader::RecipeParametersBase {
  OOPS_CONCRETE_PARAMETERS(SeaIceThicknessParameters, RecipeParametersBase)

 public:
  oops::RequiredParameter<std::string> name{"recipe name", this};
};

/*! \brief Thickness per ice area from a per-cell volume and the ice fraction
 *
 *  \details thickness = volume / sea_ice_area_fraction where the fraction is positive, 0 elsewhere
 *           (as soca's read of CICE history files: hi_h / aice_h, hs_h / aice_h). Nonlinear only:
 *           the variational analysis uses the thicknesses directly.
 */
class SeaIceThicknessBase : public vader::RecipeBase {
 public:
    SeaIceThicknessBase(const std::string & volume, const std::string & thickness)
      : volume_(volume), thickness_(thickness) {}

    oops::Variable product() const override { return oops::Variable{thickness_}; }
    oops::Variables ingredients() const override;
    size_t productLevels(const atlas::FieldSet &) const override { return 1; }
    atlas::FunctionSpace productFunctionSpace(const atlas::FieldSet &) const override;
    void executeNL(atlas::FieldSet &) override;

 private:
    const std::string volume_;
    const std::string thickness_;
};

class SeaIceThickness_A : public SeaIceThicknessBase {
 public:
    static const char Name[];
    typedef SeaIceThicknessParameters Parameters_;
    SeaIceThickness_A(const Parameters_ &, const vader::VaderConfigVars &)
      : SeaIceThicknessBase("sea_ice_volume", "sea_ice_thickness") {}
    std::string name() const override { return Name; }
};

class SeaIceSnowThickness_A : public SeaIceThicknessBase {
 public:
    static const char Name[];
    typedef SeaIceThicknessParameters Parameters_;
    SeaIceSnowThickness_A(const Parameters_ &, const vader::VaderConfigVars &)
      : SeaIceThicknessBase("sea_ice_snow_volume", "sea_ice_snow_thickness") {}
    std::string name() const override { return Name; }
};

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

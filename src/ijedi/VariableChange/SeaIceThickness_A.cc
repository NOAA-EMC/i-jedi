/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/VariableChange/SeaIceThickness.h"

#include <string>
#include <vector>

#include "oops/util/for_each.h"
#include "oops/util/Logger.h"

namespace ijedi
{
// -------------------------------------------------------------------------------------------------

const char SeaIceThickness_A::Name[] = "SeaIceThickness_A";
const char SeaIceSnowThickness_A::Name[] = "SeaIceSnowThickness_A";

static vader::RecipeMaker<SeaIceThickness_A> makerSeaIceThickness_(SeaIceThickness_A::Name);
static vader::RecipeMaker<SeaIceSnowThickness_A>
  makerSeaIceSnowThickness_(SeaIceSnowThickness_A::Name);

// -------------------------------------------------------------------------------------------------

oops::Variables SeaIceThicknessBase::ingredients() const
{
    return oops::Variables{std::vector<std::string>{volume_, "sea_ice_area_fraction"}};
}

// -------------------------------------------------------------------------------------------------

atlas::FunctionSpace SeaIceThicknessBase::productFunctionSpace(
    const atlas::FieldSet & afieldset) const
{
    return afieldset.field("sea_ice_area_fraction").functionspace();
}

// -------------------------------------------------------------------------------------------------

void SeaIceThicknessBase::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering " << name() << "::executeNL" << std::endl;

    util::for_each_column(
        [&](const auto volume, const auto fraction, auto thickness) {
            thickness(0) = fraction(0) > 0.0 ? volume(0) / fraction(0) : 0.0;
        },
        afieldset.field(volume_),
        afieldset.field("sea_ice_area_fraction"),
        afieldset.field(thickness_));

    oops::Log::trace() << "leaving " << name() << "::executeNL" << std::endl;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/VariableChange/SeaSurfaceTemperature.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "oops/util/for_each.h"
#include "oops/util/Logger.h"

namespace ijedi
{
// ------------------------------------------------------------------------------------------------

// Static attribute initialization
const char SeaSurfaceTemperature_A::Name[] = "SeaSurfaceTemperature_A";
const oops::Variables SeaSurfaceTemperature_A::Ingredients{std::vector<std::string>{
                                                           "sea_water_potential_temperature"}};

// Register the maker
static vader::RecipeMaker<SeaSurfaceTemperature_A> makerSSTempToSSTemp_(
                    SeaSurfaceTemperature_A::Name);

SeaSurfaceTemperature_A::SeaSurfaceTemperature_A(const Parameters_ & params,
                             const vader::VaderConfigVars & configVariables):
                                            configVariables_{configVariables}
{
    oops::Log::trace() << "SeaSurfaceTemperature_A::SeaSurfaceTemperature_A(params)"
        << std::endl;
}

std::string SeaSurfaceTemperature_A::name() const
{
    return SeaSurfaceTemperature_A::Name;
}

oops::Variable SeaSurfaceTemperature_A::product() const
{
    return oops::Variable{"sea_surface_temperature"};
}

oops::Variables SeaSurfaceTemperature_A::ingredients() const
{
    return SeaSurfaceTemperature_A::Ingredients;
}

size_t SeaSurfaceTemperature_A::productLevels(const atlas::FieldSet & afieldset) const
{
    return 1;
}

atlas::FunctionSpace SeaSurfaceTemperature_A::productFunctionSpace
                                                (const atlas::FieldSet & afieldset) const
{
    return afieldset.field("sea_water_potential_temperature").functionspace();
}

// -------------------------------------------------------------------------------------------------

void SeaSurfaceTemperature_A::executeNL(atlas::FieldSet & afieldset)
{
    oops::Log::trace() << "entering SeaSurfaceTemperature_A::executeNL function"
        << std::endl;

    // Get fields
    atlas::Field temperature = afieldset.field("sea_water_potential_temperature");
    atlas::Field surface_temperature = afieldset.field("sea_surface_temperature");

    util::for_each_column(
        [&](
            const auto t,
            auto sst) {
                // For now, just use the first level of sea water potential temperature
                // as sea surface temperature
                sst(0) = t(0);
            },
        temperature,
        surface_temperature);

    oops::Log::trace() << "leaving SeaSurfaceTemperature_A::executeNL function" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void SeaSurfaceTemperature_A::executeTL(atlas::FieldSet & afieldsetTL,
                                        const atlas::FieldSet & afieldsetTraj) {
    oops::Log::trace() << "entering SeaSurfaceTemperature_A::executeTL function"
        << std::endl;

    // Get fields
    atlas::Field temperature_tl = afieldsetTL.field("sea_water_potential_temperature");
    atlas::Field sst_tl = afieldsetTL.field("sea_surface_temperature");

    util::for_each_column(
        [&](
            const auto dt,
            auto dsst) {
                // For now, just use the first level of sea water potential temperature
                // as sea surface temperature
                dsst(0) = dt(0);
            },
        temperature_tl,
        sst_tl);

    oops::Log::trace() << "leaving SeaSurfaceTemperature_A::executeTL function" << std::endl;
}

// -------------------------------------------------------------------------------------------------


void SeaSurfaceTemperature_A::executeAD(atlas::FieldSet & afieldsetAD,
                                        const atlas::FieldSet & afieldsetTraj) {
    oops::Log::trace() << "entering SeaSurfaceTemperature_A::executeAD function"
        << std::endl;

    // Get fields
    atlas::Field temperature_ad = afieldsetAD.field("sea_water_potential_temperature");
    atlas::Field sst_ad = afieldsetAD.field("sea_surface_temperature");

    util::for_each_column(
        [&](
            auto dt,
            auto dsst) {
                // For now, just use the first level of sea water potential temperature
                // as sea surface temperature
                dt(0) += dsst(0);
                dsst(0) = 0.;
            },
        temperature_ad,
        sst_ad);

    oops::Log::trace() << "leaving SeaSurfaceTemperature_A::executeAD function" << std::endl;
}

}  // namespace ijedi

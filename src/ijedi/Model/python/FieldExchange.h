/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "ijedi/Model/python/PyModelBridge.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

/// \brief Move fields between i-jedi's MPI-partitioned atlas FieldSet and the global,
///        level-major buffers a Python model works with.
///
/// Both directions are collective and must be called on every task. The buffers are populated
/// on the root task only, which is the task whose bridge does the encoding and decoding.
///
/// These operate on a FieldSet rather than on a State so that Increment reuses them unchanged
/// when the tangent linear and adjoint arrive.

/// \brief Gather \p fields to the root task, renaming JEDI long names to the adapter's names
///        and multiplying by \p scaling on the way out.
void gatherToGlobal(const atlas::FunctionSpace & functionSpace,
                    const atlas::FieldSet & fields,
                    const std::map<std::string, std::string> & nameMap,
                    const std::map<std::string, double> & scaling,
                    PyModelBridge::Fields & global);

/// \brief Scatter \p global from the root task back into \p fields, undoing the renaming and
///        the scaling. Variables absent from \p global are left untouched, so a model that
///        returns fewer variables than it was given does not zero the rest of the state.
void scatterFromGlobal(const atlas::FunctionSpace & functionSpace,
                       const PyModelBridge::Fields & global,
                       const std::map<std::string, std::string> & nameMap,
                       const std::map<std::string, double> & scaling,
                       atlas::FieldSet & fields);

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

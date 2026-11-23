/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

namespace oops { class Variables; }
namespace vader { class Vader; }

namespace ijedi {

std::map<std::string, std::vector<std::string>> vaderDefaultCookbook();
oops::Variables requiredVariables(vader::Vader &, const oops::Variables &, const oops::Variables &);

}  // namespace ijedi

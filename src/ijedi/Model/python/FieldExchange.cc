/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/python/FieldExchange.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/option.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

namespace {

std::string adapterName(const std::map<std::string, std::string> & nameMap,
                        const std::string & jediName) {
  const auto it = nameMap.find(jediName);
  return it == nameMap.end() ? jediName : it->second;
}

double scaleFor(const std::map<std::string, double> & scaling, const std::string & jediName) {
  const auto it = scaling.find(jediName);
  return it == scaling.end() ? 1.0 : it->second;
}

}  // namespace

// -------------------------------------------------------------------------------------------------

void gatherToGlobal(const atlas::FunctionSpace & functionSpace,
                    const atlas::FieldSet & fields,
                    const std::map<std::string, std::string> & nameMap,
                    const std::map<std::string, double> & scaling,
                    PyModelBridge::Fields & global) {
  global.clear();

  for (const atlas::Field & field : fields) {
    const std::string & jediName = field.name();
    const int nLevels = field.shape(1);

    atlas::Field globalField = functionSpace.createField<double>(
        atlas::option::name(jediName) | atlas::option::levels(nLevels) |
        atlas::option::global());
    functionSpace.gather(field, globalField);

    // Empty on every task but the root, which is where the bridge runs.
    const size_t nNodes = static_cast<size_t>(globalField.shape(0));
    if (nNodes == 0) continue;

    const double scale = scaleFor(scaling, jediName);
    const auto view = atlas::array::make_view<double, 2>(globalField);

    std::vector<double> buffer(static_cast<size_t>(nLevels) * nNodes);
    for (int k = 0; k < nLevels; ++k) {
      for (size_t n = 0; n < nNodes; ++n) {
        buffer[static_cast<size_t>(k) * nNodes + n] = view(n, k) * scale;
      }
    }
    global[adapterName(nameMap, jediName)] = std::move(buffer);
  }
}

// -------------------------------------------------------------------------------------------------

void scatterFromGlobal(const atlas::FunctionSpace & functionSpace,
                       const PyModelBridge::Fields & global,
                       const std::map<std::string, std::string> & nameMap,
                       const std::map<std::string, double> & scaling,
                       atlas::FieldSet & fields) {
  for (atlas::Field & field : fields) {
    const std::string & jediName = field.name();
    const std::string & fromName = adapterName(nameMap, jediName);
    const int nLevels = field.shape(1);

    atlas::Field globalField = functionSpace.createField<double>(
        atlas::option::name(jediName) | atlas::option::levels(nLevels) |
        atlas::option::global());

    // Seed the global field with the current contents. This is collective and must happen
    // on every task and for every field, whatever the model returned - matching it to what
    // the root happens to find in the map would deadlock. It also gives the right answer for
    // a variable the model does not return: the field round-trips unchanged rather than
    // being zeroed.
    functionSpace.gather(field, globalField);

    // Only the root holds the gathered values, and only the root can see what came back.
    const size_t nNodes = static_cast<size_t>(globalField.shape(0));
    if (nNodes > 0) {
      const auto it = global.find(fromName);
      if (it != global.end()) {
        const std::vector<double> & buffer = it->second;
        if (buffer.size() != static_cast<size_t>(nLevels) * nNodes) {
          throw eckit::BadValue("scatterFromGlobal: the model returned '" + fromName +
                                "' with the wrong number of values for '" + jediName + "'",
                                Here());
        }
        const double scale = scaleFor(scaling, jediName);
        auto view = atlas::array::make_view<double, 2>(globalField);
        for (int k = 0; k < nLevels; ++k) {
          for (size_t n = 0; n < nNodes; ++n) {
            view(n, k) = buffer[static_cast<size_t>(k) * nNodes + n] / scale;
          }
        }
      }
    }

    functionSpace.scatter(globalField, field);
  }

  functionSpace.haloExchange(fields);
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

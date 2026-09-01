/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/python/StubPyModelBridge.h"

#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

namespace {

/// Read one half of the declaration's variable lists. Each entry is either a bare adapter
/// name or a mapping from the adapter's name to the JEDI long name it corresponds to, which
/// is the same shape the real adapters report.
std::map<std::string, std::string> readVariables(const eckit::Configuration & config,
                                                 const std::string & key) {
  std::map<std::string, std::string> variables;
  if (!config.has(key)) return variables;
  const eckit::LocalConfiguration sub(config, key);
  for (const std::string & name : sub.keys()) {
    variables[name] = sub.getString(name);
  }
  return variables;
}

}  // namespace

// -------------------------------------------------------------------------------------------------

StubPyModelBridge::StubPyModelBridge(const eckit::Configuration & config) {
  const eckit::LocalConfiguration declConf(config, "declaration");

  declaration_.description = declConf.getString("description", "stub python model");
  declaration_.timestep = util::Duration(declConf.getString("timestep"));
  if (declConf.has("levels")) {
    declaration_.levelsPa = declConf.getDoubleVector("levels");
  }
  declaration_.inputVariables = readVariables(declConf, "input variables");
  declaration_.forcingVariables = readVariables(declConf, "forcing variables");
  declaration_.supportsLinear = false;

  if (declaration_.inputVariables.empty()) {
    throw eckit::BadValue("StubPyModelBridge: the declaration lists no input variables",
                          Here());
  }

  oops::Log::info() << "StubPyModelBridge: " << declaration_.description
                    << ", internal timestep " << declaration_.timestep << std::endl;
}

// -------------------------------------------------------------------------------------------------

void StubPyModelBridge::encode(const Fields & inputs, const util::DateTime & validTime) {
  oops::Log::trace() << classname() << " encode at " << validTime << std::endl;
  for (const auto & entry : declaration_.inputVariables) {
    if (inputs.find(entry.first) == inputs.end()) {
      throw eckit::BadValue("StubPyModelBridge::encode: the declared input variable '" +
                            entry.first + "' was not supplied", Here());
    }
  }
  state_ = inputs;
  encoded_ = true;
}

// -------------------------------------------------------------------------------------------------

void StubPyModelBridge::advance(int steps) {
  if (!encoded_) {
    throw eckit::BadValue("StubPyModelBridge::advance called before encode", Here());
  }
  // Persistence: the state is unchanged by the passage of time.
  oops::Log::trace() << classname() << " advance " << steps << " internal steps" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void StubPyModelBridge::decode(Fields & outputs) const {
  if (!encoded_) {
    throw eckit::BadValue("StubPyModelBridge::decode called before encode", Here());
  }
  outputs = state_;
}

// -------------------------------------------------------------------------------------------------

void StubPyModelBridge::reset() {
  state_.clear();
  encoded_ = false;
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

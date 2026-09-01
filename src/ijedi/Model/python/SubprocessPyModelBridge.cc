/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Model/python/SubprocessPyModelBridge.h"

#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/parser/JSONParser.h"
#include "eckit/value/Value.h"

#include "oops/util/Logger.h"

namespace ijedi {

// -------------------------------------------------------------------------------------------------

namespace {

// The file descriptor the server expects its end of the socketpair on. stdin, stdout and
// stderr are left alone so the frameworks can write to them without corrupting the protocol.
constexpr int kControlFd = 3;

std::string jsonEscape(const std::string & value) {
  std::string out;
  out.reserve(value.size() + 2);
  for (char c : value) {
    switch (c) {
      case '"':  out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n";  break;
      case '\r': out += "\\r";  break;
      case '\t': out += "\\t";  break;
      default:   out += c;      break;
    }
  }
  return out;
}

/// Render the adapter's own configuration block as JSON. Values are passed through as
/// strings; an adapter coerces what it needs, which keeps this free of any assumption about
/// what a particular model's configuration looks like.
std::string adapterConfigJson(const eckit::Configuration & config) {
  std::string out = "{";
  if (config.has("adapter config")) {
    const eckit::LocalConfiguration sub(config, "adapter config");
    bool first = true;
    for (const std::string & key : sub.keys()) {
      if (!first) out += ",";
      first = false;
      out += "\"" + jsonEscape(key) + "\":\"" + jsonEscape(sub.getString(key)) + "\"";
    }
  }
  return out + "}";
}

std::string valueToString(const eckit::Value & v) {
  return v.as<std::string>();
}

}  // namespace

// -------------------------------------------------------------------------------------------------

SubprocessPyModelBridge::SubprocessPyModelBridge(const eckit::Configuration & config,
                                                 const eckit::mpi::Comm & comm) {
  isRoot_ = (comm.rank() == 0);

  std::string response;
  bool succeeded = true;
  std::string rootError;

  if (isRoot_) {
    // Bring the adapter up, then ask it what it is. Both failures - a missing adapter and a
    // model that will not load - surface here, before any forecast has started.
    try {
      start(config);
      const std::string init = "{\"op\":\"init\",\"adapter\":\"" + jsonEscape(adapter_) +
                               "\",\"config\":" + adapterConfigJson(config) + "}";
      call(init, {});
      response = call("{\"op\":\"declare\"}", {});
    } catch (const std::exception & e) {
      succeeded = false;
      rootError = e.what();
    }
  }

  response = shareDeclaration(comm, response, succeeded, rootError);

  const eckit::Value parsed = eckit::JSONParser::decodeString(response);
  const eckit::Value decl = parsed["declaration"];

  declaration_.description = valueToString(decl["description"]);
  declaration_.timestep = util::Duration(static_cast<int64_t>(decl["timestep_seconds"]));
  declaration_.supportsLinear = decl["supports_linear"];

  const eckit::Value levels = decl["levels_pa"];
  for (size_t i = 0; i < levels.size(); ++i) {
    declaration_.levelsPa.push_back(levels[i]);
  }

  for (const char * which : {"input_variables", "forcing_variables"}) {
    const eckit::Value vars = decl[which];
    const eckit::Value keys = vars.keys();
    auto & target = std::string(which) == "input_variables"
                    ? declaration_.inputVariables : declaration_.forcingVariables;
    for (size_t i = 0; i < keys.size(); ++i) {
      const std::string key = valueToString(keys[i]);
      target[key] = valueToString(vars[keys[i]]);
    }
  }

  if (isRoot_) {
    oops::Log::info() << "SubprocessPyModelBridge: " << declaration_.description
                      << " (pid " << childPid_ << "), internal timestep "
                      << declaration_.timestep << std::endl;
  }
}

// -------------------------------------------------------------------------------------------------

std::string SubprocessPyModelBridge::shareDeclaration(const eckit::mpi::Comm & comm,
                                                      const std::string & json,
                                                      bool rootSucceeded,
                                                      const std::string & rootError) {
  if (comm.size() == 1) {
    if (!rootSucceeded) throw eckit::Exception(rootError, Here());
    return json;
  }

  int ok = rootSucceeded ? 1 : 0;
  comm.broadcast(ok, 0);
  if (ok == 0) {
    throw eckit::Exception(isRoot_
        ? rootError
        : std::string("SubprocessPyModelBridge: the root task could not start the model "
                      "server; its error is reported on that task."), Here());
  }

  size_t length = json.size();
  comm.broadcast(length, 0);
  std::vector<char> buffer(length);
  if (isRoot_) std::copy(json.begin(), json.end(), buffer.begin());
  comm.broadcast(buffer, 0);
  return std::string(buffer.begin(), buffer.end());
}

// -------------------------------------------------------------------------------------------------

SubprocessPyModelBridge::~SubprocessPyModelBridge() {
  stop();
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::start(const eckit::Configuration & config) {
  adapter_ = config.getString("adapter");

  // The interpreter is the one in the model's own virtual environment, not the one JEDI was
  // built against - that is the whole point of running out of process.
  python_ = config.getString("python executable", "");
  if (python_.empty()) {
    const char * envPrefix = ::getenv("IJEDI_PYTHON_ENV");
    if (envPrefix != nullptr) python_ = std::string(envPrefix) + "/bin/python";
  }
  if (python_.empty()) {
    throw eckit::BadValue("SubprocessPyModelBridge: no python interpreter configured. Set "
                          "'python executable' in the model block, or IJEDI_PYTHON_ENV to a "
                          "virtual environment created by tools/setup_python_model_env.sh",
                          Here());
  }

  script_ = config.getString("server script", IJEDI_PYTHON_MODULE_DIR "/ijedi_model_server.py");

  int fds[2];
  if (::socketpair(AF_UNIX, SOCK_STREAM, 0, fds) != 0) {
    throw eckit::Exception(std::string("SubprocessPyModelBridge: socketpair failed: ") +
                           ::strerror(errno), Here());
  }

  const pid_t pid = ::fork();
  if (pid < 0) {
    ::close(fds[0]);
    ::close(fds[1]);
    throw eckit::Exception(std::string("SubprocessPyModelBridge: fork failed: ") +
                           ::strerror(errno), Here());
  }

  if (pid == 0) {
    // Child. Put our end of the socket on the agreed descriptor and hand over to Python.
    ::close(fds[0]);
    if (fds[1] != kControlFd) {
      ::dup2(fds[1], kControlFd);
      ::close(fds[1]);
    }
    // PYTHONPATH would drag the host's spack-stack packages into the model's environment and
    // defeat its isolation, which is exactly what this transport exists to provide.
    ::unsetenv("PYTHONPATH");
    ::execl(python_.c_str(), python_.c_str(), script_.c_str(), static_cast<char *>(nullptr));
    // Only reached if exec failed; the parent sees the socket close.
    ::_exit(127);
  }

  ::close(fds[1]);
  socket_ = fds[0];
  childPid_ = pid;
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::stop() {
  if (socket_ >= 0) {
    const std::string bye = "{\"op\":\"shutdown\"}";
    const uint64_t length = bye.size();
    // Best effort: if the child has already gone this fails harmlessly.
    ::send(socket_, &length, sizeof(length), MSG_NOSIGNAL);
    ::send(socket_, bye.data(), bye.size(), MSG_NOSIGNAL);
    ::close(socket_);
    socket_ = -1;
  }
  if (childPid_ > 0) {
    int status = 0;
    ::waitpid(childPid_, &status, 0);
    childPid_ = -1;
  }
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::writeAll(const void * data, size_t bytes) const {
  const char * p = static_cast<const char *>(data);
  while (bytes > 0) {
    const ssize_t n = ::send(socket_, p, bytes, MSG_NOSIGNAL);
    if (n <= 0) {
      if (n < 0 && errno == EINTR) continue;
      throw eckit::Exception("SubprocessPyModelBridge: the model server closed the "
                             "connection while being written to. Its stderr is in the log "
                             "above.", Here());
    }
    p += n;
    bytes -= static_cast<size_t>(n);
  }
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::readAll(void * data, size_t bytes) const {
  char * p = static_cast<char *>(data);
  while (bytes > 0) {
    const ssize_t n = ::recv(socket_, p, bytes, 0);
    if (n <= 0) {
      if (n < 0 && errno == EINTR) continue;
      throw eckit::Exception("SubprocessPyModelBridge: the model server closed the "
                             "connection. Its stderr is in the log above.", Here());
    }
    p += n;
    bytes -= static_cast<size_t>(n);
  }
}

// -------------------------------------------------------------------------------------------------

std::string SubprocessPyModelBridge::call(const std::string & requestJson,
                                          const Fields & send,
                                          Fields * receive) const {
  // Requests that carry arrays name them, and their sizes, in the header; the payload is the
  // raw doubles in that order.
  std::string header = requestJson;
  if (!send.empty()) {
    std::string arrays = ",\"arrays\":[";
    bool first = true;
    for (const auto & entry : send) {
      if (!first) arrays += ",";
      first = false;
      arrays += "{\"name\":\"" + jsonEscape(entry.first) + "\",\"size\":" +
                std::to_string(entry.second.size()) + "}";
    }
    arrays += "]";
    header.insert(header.size() - 1, arrays);
  }

  const uint64_t headerLength = header.size();
  writeAll(&headerLength, sizeof(headerLength));
  writeAll(header.data(), header.size());
  for (const auto & entry : send) {
    writeAll(entry.second.data(), entry.second.size() * sizeof(double));
  }

  uint64_t responseLength = 0;
  readAll(&responseLength, sizeof(responseLength));
  std::string response(responseLength, '\0');
  readAll(&response[0], responseLength);

  const eckit::Value parsed = eckit::JSONParser::decodeString(response);
  if (valueToString(parsed["status"]) != "ok") {
    std::stringstream msg;
    msg << "SubprocessPyModelBridge: the model adapter failed: "
        << valueToString(parsed["message"]) << "\n"
        << valueToString(parsed["traceback"]);
    throw eckit::Exception(msg.str(), Here());
  }

  if (receive != nullptr && parsed.contains("arrays")) {
    const eckit::Value arrays = parsed["arrays"];
    for (size_t i = 0; i < arrays.size(); ++i) {
      const std::string name = valueToString(arrays[i]["name"]);
      const size_t size = static_cast<size_t>(static_cast<int64_t>(arrays[i]["size"]));
      std::vector<double> buffer(size);
      readAll(buffer.data(), size * sizeof(double));
      (*receive)[name] = std::move(buffer);
    }
  }

  return response;
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::encode(const Fields & inputs, const util::DateTime & validTime) {
  if (!isRoot_) return;
  const std::string request = "{\"op\":\"encode\",\"valid_time\":\"" +
                              validTime.toString() + "\"}";
  call(request, inputs);
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::advance(int steps) {
  if (!isRoot_) return;
  call("{\"op\":\"advance\",\"steps\":" + std::to_string(steps) + "}", {});
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::decode(Fields & outputs) const {
  outputs.clear();
  if (!isRoot_) return;
  call("{\"op\":\"decode\"}", {}, &outputs);
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::reset() {
  if (!isRoot_) return;
  call("{\"op\":\"reset\"}", {});
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::setTrajectory(const Fields & inputs,
                                            const util::DateTime & validTime) {
  if (!isRoot_) return;
  call("{\"op\":\"set_trajectory\",\"valid_time\":\"" + validTime.toString() + "\"}",
       inputs);
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::advanceTL(const Fields & dxIn, Fields & dxOut, int steps) {
  dxOut.clear();
  if (!isRoot_) return;
  call("{\"op\":\"advance_tl\",\"steps\":" + std::to_string(steps) + "}", dxIn, &dxOut);
}

// -------------------------------------------------------------------------------------------------

void SubprocessPyModelBridge::advanceAD(const Fields & dxOut, Fields & dxIn, int steps) {
  dxIn.clear();
  if (!isRoot_) return;
  call("{\"op\":\"advance_ad\",\"steps\":" + std::to_string(steps) + "}", dxOut, &dxIn);
}

// -------------------------------------------------------------------------------------------------

}  // namespace ijedi

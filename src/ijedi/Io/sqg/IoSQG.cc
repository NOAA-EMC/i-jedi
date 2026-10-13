/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ijedi/Io/sqg/IoSQG.h"

#include <netcdf.h>

#include <algorithm>
#include <cstring>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ijedi/Model/sqg/ModelSQG.h"
#include "ijedi/Model/sqg/SQGGrid.h"

namespace ijedi
{
  // -------------------------------------------------------------------------------------------------
  static IoMaker<IoSQG> makerIoSQG_("sqg");

  namespace {
    void ncCheck(int status, const std::string & what) {
      if (status != NC_NOERR) {
        throw eckit::Exception("IoSQG: " + what + ": " + nc_strerror(status), Here());
      }
    }

    // Name of the state variable in the file ("field io names" may override the default "pv")
    std::string fileVarName(const eckit::LocalConfiguration & names) {
      return names.getString(ModelSQG::varname, "pv");
    }

    void checkFields(const atlas::FieldSet & x) {
      for (const auto & field : x) {
        if (field.name() != ModelSQG::varname) {
          throw eckit::BadValue("IoSQG: only " + std::string(ModelSQG::varname)
                                + " is supported, not " + field.name(), Here());
        }
      }
    }
  }  // namespace

  // -------------------------------------------------------------------------------------------------
  IoSQG::IoSQG(const Geometry &geom, const Parameters_ &params)
      : IoBase(geom, params.toConfiguration()),
        geom_(geom),
        filename_(params.filename),
        timeIndex_(params.timeIndex),
        scaleToK_(params.scaleToK)
  {}

  // -------------------------------------------------------------------------------------------------
  void IoSQG::read(atlas::FieldSet &x, const eckit::LocalConfiguration &names,
                   const eckit::LocalConfiguration &) const
  {
    oops::Log::trace() << classname() << " read starting" << std::endl;
    checkFields(x);
    const SQGGrid grid(geom_);
    const int N = grid.N, NN = N * N;

    int ncid, varid, dimids[4];
    ncCheck(nc_open(filename_.c_str(), NC_NOWRITE, &ncid), "open " + filename_);
    const std::string vname = fileVarName(names);
    ncCheck(nc_inq_varid(ncid, vname.c_str(), &varid), "variable " + vname);
    ncCheck(nc_inq_vardimid(ncid, varid, dimids), "dimensions of " + vname);
    size_t len[4];
    for (int d = 0; d < 4; ++d) ncCheck(nc_inq_dimlen(ncid, dimids[d], &len[d]), "dim length");
    if (len[1] != 2 || len[2] != static_cast<size_t>(N) || len[3] != static_cast<size_t>(N)) {
      throw eckit::BadValue("IoSQG: " + filename_ + " is not (t, 2, N, N) with N="
                            + std::to_string(N), Here());
    }
    const int rec = timeIndex_ < 0 ? static_cast<int>(len[0]) + timeIndex_ : timeIndex_;
    if (rec < 0 || rec >= static_cast<int>(len[0])) {
      throw eckit::BadValue("IoSQG: time index out of range for " + filename_, Here());
    }

    float scale;
    if (nc_get_att_float(ncid, varid, "scale_to_K", &scale) != NC_NOERR) {
      float f, theta0, g;
      ncCheck(nc_get_att_float(ncid, NC_GLOBAL, "f", &f), "global attribute f");
      ncCheck(nc_get_att_float(ncid, NC_GLOBAL, "theta0", &theta0), "global attribute theta0");
      ncCheck(nc_get_att_float(ncid, NC_GLOBAL, "g", &g), "global attribute g");
      scale = f * theta0 / g;
    }

    std::vector<float> pv(2 * NN);
    const size_t start[4] = {static_cast<size_t>(rec), 0, 0, 0};
    const size_t count[4] = {1, 2, static_cast<size_t>(N), static_cast<size_t>(N)};
    ncCheck(nc_get_vara_float(ncid, varid, start, count, pv.data()), "read " + vname);
    ncCheck(nc_close(ncid), "close");

    auto view = atlas::array::make_view<double, 2>(x.field(ModelSQG::varname));
    for (size_t p = 0; p < grid.nodes.size(); ++p) {
      for (int k = 0; k < 2; ++k) {
        view(grid.nodes[p], k) = static_cast<double>(pv[k * NN + grid.index[p]]) * scale;
      }
    }
    oops::Log::info() << classname() << " read " << filename_ << " record " << rec
                      << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoSQG::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &names,
                    const eckit::LocalConfiguration &) const
  {
    oops::Log::trace() << classname() << " write starting" << std::endl;
    checkFields(x);
    if (geom_.getComm().size() != 1) {
      throw eckit::NotImplemented("IoSQG::write needs a single MPI task", Here());
    }
    const SQGGrid grid(geom_);
    const int N = grid.N, NN = N * N;

    // File name, with %{valid} -> compact valid time
    std::string fname = filename_;
    const std::string token = "%{valid}";
    const size_t pos = fname.find(token);
    if (pos != std::string::npos) {
      if (!hasValidTime()) {
        throw eckit::BadValue("IoSQG: %{valid} in filename but no valid time", Here());
      }
      std::string vt = validTime().toString();
      vt.erase(std::remove(vt.begin(), vt.end(), '-'), vt.end());
      vt.erase(std::remove(vt.begin(), vt.end(), ':'), vt.end());
      fname.replace(pos, token.size(), vt);
    }

    const float scale = scaleToK_;

    const auto view = atlas::array::make_view<double, 2>(x.field(ModelSQG::varname));
    std::vector<float> pv(2 * NN);
    for (size_t p = 0; p < grid.nodes.size(); ++p) {
      for (int k = 0; k < 2; ++k) {
        pv[k * NN + grid.index[p]] = static_cast<float>(view(grid.nodes[p], k) /
                                                        static_cast<double>(scale));
      }
    }

    int ncid, dimids[4], varid;
    ncCheck(nc_create(fname.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "create " + fname);
    ncCheck(nc_def_dim(ncid, "t", NC_UNLIMITED, &dimids[0]), "def dim t");
    ncCheck(nc_def_dim(ncid, "z", 2, &dimids[1]), "def dim z");
    ncCheck(nc_def_dim(ncid, "y", N, &dimids[2]), "def dim y");
    ncCheck(nc_def_dim(ncid, "x", N, &dimids[3]), "def dim x");
    const std::string vname = fileVarName(names);
    ncCheck(nc_def_var(ncid, vname.c_str(), NC_FLOAT, 4, dimids, &varid), "def var");
    const char units[] = "m/s";
    ncCheck(nc_put_att_text(ncid, varid, "units", std::strlen(units), units), "units");
    ncCheck(nc_put_att_float(ncid, varid, "scale_to_K", NC_FLOAT, 1, &scale), "scale_to_K");
    if (hasValidTime()) {
      const std::string vt = validTime().toString();
      ncCheck(nc_put_att_text(ncid, NC_GLOBAL, "valid_time", vt.size(), vt.c_str()),
              "valid_time");
    }
    ncCheck(nc_enddef(ncid), "enddef");
    const size_t start[4] = {0, 0, 0, 0};
    const size_t count[4] = {1, 2, static_cast<size_t>(N), static_cast<size_t>(N)};
    ncCheck(nc_put_vara_float(ncid, varid, start, count, pv.data()), "write " + vname);
    ncCheck(nc_close(ncid), "close");
    oops::Log::info() << classname() << " wrote " << fname << std::endl;
  }

  // -------------------------------------------------------------------------------------------------
  void IoSQG::print(std::ostream &os) const
  {
    os << classname() << " SQG netCDF file " << filename_;
  }

  // -------------------------------------------------------------------------------------------------
}  // namespace ijedi

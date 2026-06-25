#include <netcdf.h>

#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/fv3-history/IoFV3History.h"

namespace ijedi
{
    // -----------------------------------------------------------------------------
    static IoMaker<IoFV3History> makerIoFV3_("fv3 history");
    // -----------------------------------------------------------------------------
    IoFV3History::IoFV3History(const Geometry &geom, const Parameters_ &params)
        : IoBase(geom, params.toConfiguration()), parameters_(params), geom_(geom)
    {
        util::Timer timer(classname(), "IoFV3History");
        oops::Log::trace() << classname() << " constructor starting" << std::endl;
        oops::Log::trace() << classname() << " constructor done" << std::endl;
    }
    // -----------------------------------------------------------------------------
    IoFV3History::~IoFV3History()
    {
        util::Timer timer(classname(), "~IoFV3History");
        oops::Log::trace() << classname() << " destructor starting" << std::endl;
        oops::Log::trace() << classname() << " destructor done" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void IoFV3History::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &fileioscaling) const
    {
        util::Timer timer(classname(), "read state");
        oops::Log::trace() << classname() << " read state starting" << std::endl;

        const std::string datapath = parameters_.datapath.value();

        // Build the list of file paths from atm_file and sfc_file
        std::vector<std::string> filepaths;
        filepaths.push_back(datapath + "/" + parameters_.atm_file.value());
        filepaths.push_back(datapath + "/" + parameters_.sfc_file.value());

        // Get per-file dimension name overrides (or use defaults)
        const auto &xdimOpt = parameters_.xdim.value();
        const auto &ydimOpt = parameters_.ydim.value();
        const auto &zfdimOpt = parameters_.zfdim.value();
        const auto &tiledimOpt = parameters_.tiledim.value();

        // Get the JEDI field names from fileionames keys.
        // Each key is a JEDI long name; the value is the NetCDF variable name in the file.
        const std::vector<std::string> jediNames = fileionames.keys();

        // Track which fields have been read (to avoid reading from a second file)
        std::vector<bool> fieldRead(jediNames.size(), false);
        oops::Log::info() << classname() << " reading history files for "
                          << jediNames.size() << " requested fields" << std::endl;

        // Loop over each file (e.g. atm file, sfc file)
        for (size_t ifile = 0; ifile < filepaths.size(); ++ifile)
        {
            const std::string &filepath = filepaths[ifile];
            oops::Log::info() << classname() << " reading history file: "
                              << filepath << std::endl;

            // Dimension name defaults for this file
            const std::string xdimName = (xdimOpt && ifile < xdimOpt->size())
                                             ? (*xdimOpt)[ifile]
                                             : "grid_xt";
            const std::string ydimName = (ydimOpt && ifile < ydimOpt->size())
                                             ? (*ydimOpt)[ifile]
                                             : "grid_yt";
            const std::string zfdimName = (zfdimOpt && ifile < zfdimOpt->size())
                                              ? (*zfdimOpt)[ifile]
                                              : "pfull";
            const bool hasTileDim = (tiledimOpt && ifile < tiledimOpt->size())
                                        ? (*tiledimOpt)[ifile]
                                        : true;

            // Open the NetCDF file
            int ncid;
            checkNetCDF(nc_open(filepath.c_str(), NC_NOWRITE, &ncid),
                        "opening " + filepath);

            // Read dimensions
            int dimid;
            size_t nx, ny, nz;

            checkNetCDF(nc_inq_dimid(ncid, xdimName.c_str(), &dimid),
                        "getting " + xdimName + " dimension");
            checkNetCDF(nc_inq_dimlen(ncid, dimid, &nx),
                        "getting " + xdimName + " length");

            checkNetCDF(nc_inq_dimid(ncid, ydimName.c_str(), &dimid),
                        "getting " + ydimName + " dimension");
            checkNetCDF(nc_inq_dimlen(ncid, dimid, &ny),
                        "getting " + ydimName + " length");

            checkNetCDF(nc_inq_dimid(ncid, zfdimName.c_str(), &dimid),
                        "getting " + zfdimName + " dimension");
            checkNetCDF(nc_inq_dimlen(ncid, dimid, &nz),
                        "getting " + zfdimName + " length");

            size_t ntiles = 1;
            if (hasTileDim)
            {
                checkNetCDF(nc_inq_dimid(ncid, "tile", &dimid),
                            "getting tile dimension");
                checkNetCDF(nc_inq_dimlen(ncid, dimid, &ntiles),
                            "getting tile length");
            }

            oops::Log::debug() << classname() << " file dimensions:"
                               << " nx=" << nx << " ny=" << ny << " nz=" << nz
                               << " ntiles=" << ntiles << std::endl;

            // Get local-to-global index mapping from the function space
            const auto &funcSpace = geom_.functionSpace();
            const size_t numNodes = funcSpace.size();
            auto globalIdx = atlas::array::make_view<atlas::gidx_t, 1>(
                funcSpace.global_index());
            const size_t nxy = ny * nx;  // spatial points per tile

            // Loop over requested fields and attempt to read from this file
            for (size_t ifield = 0; ifield < jediNames.size(); ++ifield)
            {
                if (fieldRead[ifield])
                    continue;  // already read from a previous file

                const std::string &jediName = jediNames[ifield];
                const std::string ncVarName = fileionames.getString(jediName);

                // Check if this variable exists in the current file
                int varid;
                if (nc_inq_varid(ncid, ncVarName.c_str(), &varid) != NC_NOERR)
                {
                    continue;  // variable not in this file, try next file
                }

                // Determine number of dimensions
                int ndims;
                checkNetCDF(nc_inq_varndims(ncid, varid, &ndims),
                            "getting ndims for " + ncVarName);

                // Variables are either:
                //   5D: (time, tile, z, y, x)
                //   4D: (time, tile, y, x)
                //   When tile is not a dimension, reduce by 1D each

                const int expected3d = hasTileDim ? 5 : 4;  // with z
                const int expected2d = hasTileDim ? 4 : 3;  // without z

                if (ndims == expected3d)
                {
                    // 3D field — read full data into buffer, then scatter to atlas field
                    const size_t totalSize = ntiles * nz * nxy;
                    std::vector<double> buffer(totalSize);

                    std::vector<size_t> start, count;
                    if (hasTileDim)
                    {
                        start = {0, 0, 0, 0, 0};
                        count = {1, ntiles, nz, ny, nx};
                    } else {
                        start = {0, 0, 0, 0};
                        count = {1, nz, ny, nx};
                    }
                    checkNetCDF(nc_get_vara_double(ncid, varid,
                                                   start.data(), count.data(),
                                                   buffer.data()),
                                "reading " + ncVarName);

                    // Populate existing rank-2 atlas field (nodes, levels)
                    atlas::Field &field = x.field(jediName);
                    auto view = atlas::array::make_view<double, 2>(field);

                    for (size_t jnode = 0; jnode < numNodes; ++jnode)
                    {
                        const size_t gid0 = static_cast<size_t>(globalIdx(jnode)) - 1;
                        const size_t tile0 = gid0 / nxy;
                        const size_t spatialIdx = gid0 % nxy;
                        for (size_t z = 0; z < nz; ++z)
                        {
                            // NetCDF C-order: buffer[tile][z][y][x]
                            const size_t bufIdx = tile0 * (nz * nxy) + z * nxy + spatialIdx;
                            view(jnode, z) = buffer[bufIdx];
                        }
                    }

                    // Apply scaling if provided
                    if (fileioscaling.has(jediName))
                    {
                        const double scale = fileioscaling.getDouble(jediName);
                        for (size_t jnode = 0; jnode < numNodes; ++jnode)
                        {
                            for (size_t z = 0; z < nz; ++z)
                            {
                                view(jnode, z) *= scale;
                            }
                        }
                    }

                    fieldRead[ifield] = true;
                    oops::Log::info() << classname() << " read 3D field: " << jediName
                                      << " (" << ncVarName << ") from " << filepath
                                      << std::endl;
                } else if (ndims == expected2d) {
                    // 2D field — read full data into buffer, then scatter to atlas field
                    const size_t totalSize = ntiles * nxy;
                    std::vector<double> buffer(totalSize);

                    std::vector<size_t> start, count;
                    if (hasTileDim)
                    {
                        start = {0, 0, 0, 0};
                        count = {1, ntiles, ny, nx};
                    } else {
                        start = {0, 0, 0};
                        count = {1, ny, nx};
                    }
                    checkNetCDF(nc_get_vara_double(ncid, varid,
                                                   start.data(), count.data(),
                                                   buffer.data()),
                                "reading " + ncVarName);

                    // Populate existing rank-2 atlas field (nodes, levels=1)
                    atlas::Field &field = x.field(jediName);
                    auto view = atlas::array::make_view<double, 2>(field);

                    for (size_t jnode = 0; jnode < numNodes; ++jnode)
                    {
                        const size_t gid0 = static_cast<size_t>(globalIdx(jnode)) - 1;
                        view(jnode, 0) = buffer[gid0];
                    }

                    // Apply scaling if provided
                    if (fileioscaling.has(jediName))
                    {
                        const double scale = fileioscaling.getDouble(jediName);
                        for (size_t jnode = 0; jnode < numNodes; ++jnode)
                        {
                            view(jnode, 0) *= scale;
                        }
                    }

                    fieldRead[ifield] = true;
                    oops::Log::info() << classname() << " read 2D field: " << jediName
                                      << " (" << ncVarName << ") from " << filepath
                                      << std::endl;
                } else {
                    oops::Log::warning() << classname() << " unexpected ndims=" << ndims
                                         << " for variable " << ncVarName
                                         << ", skipping" << std::endl;
                }
            }

            checkNetCDF(nc_close(ncid), "closing " + filepath);
        }

        // Check that all requested fields were found
        for (size_t ifield = 0; ifield < jediNames.size(); ++ifield)
        {
            if (!fieldRead[ifield])
            {
                throw eckit::Exception(classname() + "::readHistoryFiles: "
                                                     "Variable \"" +
                                       fileionames.getString(jediNames[ifield]) +
                                       "\" (JEDI name: " + jediNames[ifield] +
                                       ") not found in any of the provided files");
            }
        }

        oops::Log::info() << classname() << " finished reading "
                          << x.size() << " fields from "
                          << filepaths.size() << " history file(s)" << std::endl;

        oops::Log::trace() << classname() << " read state done" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void IoFV3History::write(
        const atlas::FieldSet &x,
        const eckit::LocalConfiguration &fileionames,
        const eckit::LocalConfiguration &fileioscaling) const
    {
        util::Timer timer(classname(), "write state");
        oops::Log::trace() << classname() << " write state starting" << std::endl;
        const std::string datapath = parameters_.datapath.value();
        const eckit::mpi::Comm &comm = geom_.getComm();
        const bool isRoot = (comm.rank() == 0);

        // Build the list of file paths from atm_file and sfc_file
        std::vector<std::string> filepaths;
        filepaths.push_back(datapath + "/" + parameters_.atm_file.value());
        filepaths.push_back(datapath + "/" + parameters_.sfc_file.value());

        // Get per-file dimension name overrides (or use defaults)
        const auto &xdimOpt = parameters_.xdim.value();
        const auto &ydimOpt = parameters_.ydim.value();
        const auto &zfdimOpt = parameters_.zfdim.value();
        const auto &tiledimOpt = parameters_.tiledim.value();

        const std::vector<std::string> jediNames = fileionames.keys();

        // Prepare for gather
        const auto &funcSpace = geom_.functionSpace();

        // Loop over each file (e.g. atm file, sfc file)
        for (size_t ifile = 0; ifile < filepaths.size(); ++ifile)
        {
            const std::string &filepath = filepaths[ifile];

            // Dimension name defaults for this file
            const std::string xdimName = (xdimOpt && ifile < xdimOpt->size())
                                             ? (*xdimOpt)[ifile]
                                             : "grid_xt";
            const std::string ydimName = (ydimOpt && ifile < ydimOpt->size())
                                             ? (*ydimOpt)[ifile]
                                             : "grid_yt";
            const std::string zfdimName = (zfdimOpt && ifile < zfdimOpt->size())
                                              ? (*zfdimOpt)[ifile]
                                              : "pfull";
            const bool hasTileDim = (tiledimOpt && ifile < tiledimOpt->size())
                                        ? (*tiledimOpt)[ifile]
                                        : true;

            size_t nx = 0, ny = 0, nz = 0, ntiles = 0;

            auto globalIdx = atlas::array::make_view<atlas::gidx_t, 1>(funcSpace.global_index());
            atlas::gidx_t maxGIdx = 0;
            for (atlas::idx_t jnode = 0; jnode < funcSpace.size(); ++jnode) {
                if (globalIdx(jnode) > maxGIdx) maxGIdx = globalIdx(jnode);
            }
            comm.allReduceInPlace(maxGIdx, eckit::mpi::max());

            if (maxGIdx == 0) throw eckit::Exception("maxGIdx is zero, cannot write history file");

            ntiles = (hasTileDim && (maxGIdx % 6 == 0)) ? 6 : 1;

            // Prefer compute-domain size from owned nodes, which excludes halo-expanded
            // nodes that may be present in the function space indexing.
            size_t nxy = 0;
            if (geom_.fields().has("owned")) {
                const atlas::Field &owned = geom_.fields().field("owned");
                auto ownedView = atlas::array::make_view<int, 2>(owned);
                atlas::idx_t ownedLocal = 0;
                for (atlas::idx_t jnode = 0; jnode < funcSpace.size(); ++jnode) {
                    if (ownedView(jnode, 0) > 0) ++ownedLocal;
                }
                atlas::idx_t ownedGlobal = ownedLocal;
                comm.allReduceInPlace(ownedGlobal, eckit::mpi::sum());
                if (ownedGlobal > 0 && ownedGlobal % static_cast<atlas::idx_t>(ntiles) == 0) {
                    const size_t nOwnedPerTile = static_cast<size_t>(ownedGlobal / ntiles);
                    const size_t side =
                        static_cast<size_t>(std::sqrt(static_cast<double>(nOwnedPerTile)));
                    if (side * side == nOwnedPerTile) {
                        nx = ny = side;
                        nxy = nOwnedPerTile;
                    }
                }
            }

            // Fallback for geometries without an owned mask or non-square owned count.
            if (nxy == 0) {
                size_t n2 = maxGIdx / ntiles;
                nx = ny = std::sqrt(n2);
                nxy = nx * ny;
            }

            if (nxy == 0) throw eckit::Exception("nx*ny is zero, cannot write history file");

            for (const auto & jediName : jediNames) {
                if (x.has(jediName)) {
                    nz = std::max(nz, static_cast<size_t>(x.field(jediName).shape(1)));
                }
            }

            if (isRoot) {
                oops::Log::info() << classname() << " writing history file: "
                                  << filepath << std::endl;
                int ncid;
                checkNetCDF(
                    nc_create(filepath.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
                    "creating " + filepath);

                int xdimid, ydimid, zdimid, tdimid, timedimid;
                checkNetCDF(nc_def_dim(ncid, xdimName.c_str(), nx, &xdimid), "def dim x");
                checkNetCDF(nc_def_dim(ncid, ydimName.c_str(), ny, &ydimid), "def dim y");
                checkNetCDF(nc_def_dim(ncid, zfdimName.c_str(), nz, &zdimid), "def dim z");
                if (hasTileDim) {
                    checkNetCDF(nc_def_dim(ncid, "tile", ntiles, &tdimid), "def dim tile");
                }
                checkNetCDF(nc_def_dim(ncid, "time", NC_UNLIMITED, &timedimid), "def dim time");

                for (const auto & jediName : jediNames) {
                    if (ifile > 0) break;
                    if (!x.has(jediName)) continue;

                    std::string ncVarName = fileionames.getString(jediName);
                    int varid;
                    int curr_nz = x.field(jediName).shape(1);

                    std::vector<int> dimids;
                    dimids.push_back(timedimid);
                    if (hasTileDim) dimids.push_back(tdimid);
                    if (curr_nz > 1) dimids.push_back(zdimid);
                    dimids.push_back(ydimid);
                    dimids.push_back(xdimid);

                    checkNetCDF(
                        nc_def_var(ncid, ncVarName.c_str(), NC_DOUBLE, dimids.size(),
                                   dimids.data(), &varid),
                        "def var " + ncVarName);
                }
                checkNetCDF(nc_enddef(ncid), "enddef");

                for (const auto & jediName : jediNames) {
                    if (ifile > 0) break;
                    if (!x.has(jediName)) continue;

                    std::string ncVarName = fileionames.getString(jediName);
                    const atlas::Field &field = x.field(jediName);
                    int curr_nz = field.shape(1);

                    atlas::Field globalField = funcSpace.createField<double>(
                        atlas::option::name(jediName)
                        | atlas::option::levels(curr_nz)
                        | atlas::option::global());

                    funcSpace.gather(field, globalField);

                    auto globalView = atlas::array::make_view<double, 2>(globalField);

                    size_t totalNodes = globalField.shape(0);
                    std::vector<double> buffer;
                    if (curr_nz > 1) {
                        buffer.assign(ntiles * curr_nz * nxy, 0.0);
                    } else {
                        buffer.assign(ntiles * nxy, 0.0);
                    }

                    atlas::Field gidxDouble =
                        funcSpace.createField<double>(atlas::option::levels(1));
                    auto gidxDoubleView = atlas::array::make_view<double, 2>(gidxDouble);
                    auto localGidxView =
                        atlas::array::make_view<atlas::gidx_t, 1>(
                            funcSpace.global_index());
                    for (atlas::idx_t j = 0; j < funcSpace.size(); ++j) {
                        gidxDoubleView(j, 0) =
                            static_cast<double>(localGidxView(j));
                    }

                    atlas::Field globalGidxDouble =
                        funcSpace.createField<double>(atlas::option::levels(1)
                                                      | atlas::option::global());
                    funcSpace.gather(gidxDouble, globalGidxDouble);
                    auto globalGidxDoubleView =
                        atlas::array::make_view<double, 2>(globalGidxDouble);

                    double scale = 1.0;
                    if (fileioscaling.has(jediName)) {
                        scale = 1.0 / fileioscaling.getDouble(jediName);
                    }

                    for (size_t gn = 0; gn < totalNodes; ++gn) {
                        size_t gid0 = static_cast<size_t>(globalGidxDoubleView(gn, 0)) - 1;
                        size_t tile0 = gid0 / nxy;
                        size_t spatialIdx = gid0 % nxy;
                        if (tile0 >= ntiles) continue;

                        if (curr_nz > 1) {
                            for (size_t z = 0; z < (size_t)curr_nz; ++z) {
                                size_t bufIdx = tile0 * (curr_nz * nxy) + z * nxy + spatialIdx;
                                buffer[bufIdx] = globalView(gn, z) * scale;
                            }
                        } else {
                            size_t bufIdx = tile0 * nxy + spatialIdx;
                            buffer[bufIdx] = globalView(gn, 0) * scale;
                        }
                    }

                    int varid;
                    checkNetCDF(nc_inq_varid(ncid, ncVarName.c_str(), &varid), "inq varid");
                    std::vector<size_t> start, count;
                    start.push_back(0);
                    if (hasTileDim) start.push_back(0);
                    if (curr_nz > 1) start.push_back(0);
                    start.push_back(0);
                    start.push_back(0);

                    count.push_back(1);
                    if (hasTileDim) count.push_back(ntiles);
                    if (curr_nz > 1) count.push_back(curr_nz);
                    count.push_back(ny);
                    count.push_back(nx);

                    checkNetCDF(
                        nc_put_vara_double(ncid, varid, start.data(), count.data(),
                                           buffer.data()),
                        "put vara " + ncVarName);
                }

                checkNetCDF(nc_close(ncid), "closing " + filepath);
            } else {
                for (const auto & jediName : jediNames) {
                    if (ifile > 0) break;
                    if (!x.has(jediName)) continue;

                    const atlas::Field &field = x.field(jediName);
                    int curr_nz = field.shape(1);
                    atlas::Field globalField = funcSpace.createField<double>(
                        atlas::option::name(jediName)
                        | atlas::option::levels(curr_nz)
                        | atlas::option::global());
                    funcSpace.gather(field, globalField);

                    atlas::Field gidxDouble =
                        funcSpace.createField<double>(atlas::option::levels(1));
                    auto gidxDoubleView = atlas::array::make_view<double, 2>(gidxDouble);
                    auto localGidxView =
                        atlas::array::make_view<atlas::gidx_t, 1>(
                            funcSpace.global_index());
                    for (atlas::idx_t j = 0; j < funcSpace.size(); ++j) {
                        gidxDoubleView(j, 0) =
                            static_cast<double>(localGidxView(j));
                    }
                    atlas::Field globalGidxDouble =
                        funcSpace.createField<double>(atlas::option::levels(1)
                                                      | atlas::option::global());
                    funcSpace.gather(gidxDouble, globalGidxDouble);
                }
            }
        }
        oops::Log::trace() << classname() << " write state done" << std::endl;
    }
    // -----------------------------------------------------------------------------
    void IoFV3History::print(std::ostream &os) const
    {
        os << classname() << " IO for Cube Sphere History files";
    }

    // -----------------------------------------------------------------------------
    void IoFV3History::checkNetCDF(int status, const std::string &operation) const
    {
        if (status != NC_NOERR)
        {
            throw eckit::Exception("NetCDF error in " + operation + ": " + nc_strerror(status));
        }
    }
    // -----------------------------------------------------------------------------
}  // namespace ijedi

#include <ostream>
#include <string>
#include <vector>

#include <netcdf.h>
#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Io/fv3/IoFV3.h"

namespace ijedi
{
    // -------------------------------------------------------------------------------------------------
    static IoMaker<IoFV3> makerIoFV3_("fv3");
    // -------------------------------------------------------------------------------------------------
    IoFV3::IoFV3(const Geometry &geom, const Parameters_ &params)
        : IoBase(geom, params.toConfiguration()), parameters_(params), geom_(geom)
    {
        util::Timer timer(classname(), "IoFV3");
        oops::Log::trace() << classname() << " constructor starting" << std::endl;
        oops::Log::trace() << classname() << " constructor done" << std::endl;
    }
    // -------------------------------------------------------------------------------------------------
    IoFV3::~IoFV3()
    {
        util::Timer timer(classname(), "~IoFV3");
        oops::Log::trace() << classname() << " destructor starting" << std::endl;
        oops::Log::trace() << classname() << " destructor done" << std::endl;
    }
    // -------------------------------------------------------------------------------------------------
    void IoFV3::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                     const eckit::LocalConfiguration &fileioscaling) const
    {
        util::Timer timer(classname(), "read state");
        oops::Log::trace() << classname() << " read state starting" << std::endl;

        const std::string source = parameters_.source.value();
        if (source == "history")
        {
            readHistoryFiles(x, fileionames, fileioscaling);
        }
        else if (source == "restart")
        {
            throw eckit::Exception("Reading restart files not yet implemented");
        }
        else
        {
            throw eckit::Exception("Invalid source parameter: " + source);
        }

        oops::Log::trace() << classname() << " read state done" << std::endl;
    }
    // -------------------------------------------------------------------------------------------------
    void IoFV3::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                      const eckit::LocalConfiguration &fileioscaling) const
    {
        util::Timer timer(classname(), "write state");
        oops::Log::trace() << classname() << " write state starting" << std::endl;

        const std::string source = parameters_.source.value();
        if (source == "history")
        {
            writeHistoryFiles(x, fileionames, fileioscaling);
        }
        else if (source == "restart")
        {
            throw eckit::Exception("Writing restart files not yet implemented");
        }
        else
        {
            throw eckit::Exception("Invalid source parameter: " + source);
        }

        oops::Log::trace() << classname() << " write state done" << std::endl;
    }
    // -------------------------------------------------------------------------------------------------
    void IoFV3::readHistoryFiles(atlas::FieldSet &fieldSet,
                                 const eckit::LocalConfiguration &fileionames,
                                 const eckit::LocalConfiguration &fileioscaling) const
    {
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
                    continue; // already read from a previous file

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
                    }
                    else
                    {
                        start = {0, 0, 0, 0};
                        count = {1, nz, ny, nx};
                    }
                    checkNetCDF(nc_get_vara_double(ncid, varid,
                                                   start.data(), count.data(),
                                                   buffer.data()),
                                "reading " + ncVarName);

                    // Populate existing rank-2 atlas field (nodes, levels)
                    atlas::Field &field = fieldSet.field(jediName);
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
                }
                else if (ndims == expected2d)
                {
                    // 2D field — read full data into buffer, then scatter to atlas field
                    const size_t totalSize = ntiles * nxy;
                    std::vector<double> buffer(totalSize);

                    std::vector<size_t> start, count;
                    if (hasTileDim)
                    {
                        start = {0, 0, 0, 0};
                        count = {1, ntiles, ny, nx};
                    }
                    else
                    {
                        start = {0, 0, 0};
                        count = {1, ny, nx};
                    }
                    checkNetCDF(nc_get_vara_double(ncid, varid,
                                                   start.data(), count.data(),
                                                   buffer.data()),
                                "reading " + ncVarName);

                    // Populate existing rank-2 atlas field (nodes, levels=1)
                    atlas::Field &field = fieldSet.field(jediName);
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
                }
                else
                {
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
                          << fieldSet.size() << " fields from "
                          << filepaths.size() << " history file(s)" << std::endl;
    }
    // -------------------------------------------------------------------------------------------------
    void IoFV3::writeHistoryFiles(const atlas::FieldSet &fieldSet,
                                  const eckit::LocalConfiguration &fileionames,
                                  const eckit::LocalConfiguration &fileioscaling) const
    {
        const std::string datapath = parameters_.datapath.value();
        const eckit::mpi::Comm &comm = geom_.comm();
        const size_t myRank = comm.rank();

        // Build the list of file paths from atm_file and sfc_file
        std::vector<std::string> filepaths;
        filepaths.push_back(datapath + "/" + parameters_.atm_file.value());
        filepaths.push_back(datapath + "/" + parameters_.sfc_file.value());

        // Get per-file dimension name overrides (or use defaults)
        const auto &xdimOpt = parameters_.xdim.value();
        const auto &ydimOpt = parameters_.ydim.value();
        const auto &zfdimOpt = parameters_.zfdim.value();
        const auto &tiledimOpt = parameters_.tiledim.value();

        // Grid dimensions from geometry
        const eckit::Configuration &geomVars = geom_.geomVariables();
        const size_t nx = static_cast<size_t>(geomVars.getInt("npx")) - 1;
        const size_t ny = static_cast<size_t>(geomVars.getInt("npy")) - 1;
        const size_t nz = static_cast<size_t>(geomVars.getInt("npz"));
        const size_t ntiles = static_cast<size_t>(geomVars.getInt("ntiles"));
        const size_t nxy = ny * nx;

        // Function space info for global index mapping
        const auto &funcSpace = geom_.functionSpace();
        const size_t numNodes = funcSpace.size();
        auto globalIdx = atlas::array::make_view<atlas::gidx_t, 1>(
            funcSpace.global_index());
        auto ghostView = atlas::array::make_view<int, 1>(funcSpace.ghost());

        // Get the JEDI field names from fileionames keys.
        const std::vector<std::string> jediNames = fileionames.keys();

        // Track which fields have been written
        std::vector<bool> fieldWritten(jediNames.size(), false);
        oops::Log::info() << classname() << " writing history files for "
                          << jediNames.size() << " fields" << std::endl;

        // Loop over each output file
        for (size_t ifile = 0; ifile < filepaths.size(); ++ifile)
        {
            const std::string &filepath = filepaths[ifile];

            // Dimension names for this file
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

            // Determine which fields belong to this file by checking which exist
            // in the source file (same logic as read — try to find var in file).
            // For write we need to know if a field is 2D or 3D. We use the atlas
            // field's number of levels: levels > 1 means 3D, levels == 1 means 2D.

            // Collect fields to write to this file: check if the variable exists
            // in the source file (if it exists) or write all remaining fields to
            // the first file that gets them.
            // Strategy: try opening the source file to probe variable existence.
            // If the source file doesn't exist, write all remaining fields.
            int srcNcid = -1;
            bool srcExists = (nc_open(filepath.c_str(), NC_NOWRITE, &srcNcid) == NC_NOERR);

            // Collect indices of fields to write to this file
            std::vector<size_t> fieldsForFile;
            for (size_t ifield = 0; ifield < jediNames.size(); ++ifield)
            {
                if (fieldWritten[ifield])
                    continue;

                const std::string ncVarName = fileionames.getString(jediNames[ifield]);
                if (srcExists)
                {
                    int varid;
                    if (nc_inq_varid(srcNcid, ncVarName.c_str(), &varid) == NC_NOERR)
                    {
                        fieldsForFile.push_back(ifield);
                    }
                }
                else
                {
                    // No source file — assign all remaining fields here
                    fieldsForFile.push_back(ifield);
                }
            }
            if (srcExists)
            {
                nc_close(srcNcid);
            }

            if (fieldsForFile.empty())
                continue;

            oops::Log::info() << classname() << " writing " << fieldsForFile.size()
                              << " fields to history file: " << filepath << std::endl;

            // For each field, gather data from all ranks into a global buffer
            // and have rank 0 write the NetCDF file.

            // Create the NetCDF file on rank 0
            int ncid = -1;
            if (myRank == 0)
            {
                checkNetCDF(nc_create(filepath.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid),
                            "creating " + filepath);

                // Define dimensions
                int timeDimId, tileDimId, zDimId, yDimId, xDimId;
                checkNetCDF(nc_def_dim(ncid, "time", NC_UNLIMITED, &timeDimId),
                            "defining time dimension");
                if (hasTileDim)
                {
                    checkNetCDF(nc_def_dim(ncid, "tile", ntiles, &tileDimId),
                                "defining tile dimension");
                }
                checkNetCDF(nc_def_dim(ncid, zfdimName.c_str(), nz, &zDimId),
                            "defining " + zfdimName + " dimension");
                checkNetCDF(nc_def_dim(ncid, ydimName.c_str(), ny, &yDimId),
                            "defining " + ydimName + " dimension");
                checkNetCDF(nc_def_dim(ncid, xdimName.c_str(), nx, &xDimId),
                            "defining " + xdimName + " dimension");

                // Define variables
                for (const size_t ifield : fieldsForFile)
                {
                    const std::string &jediName = jediNames[ifield];
                    const std::string ncVarName = fileionames.getString(jediName);
                    const atlas::Field &field = fieldSet.field(jediName);
                    const size_t nlevs = field.levels();

                    int varid;
                    if (nlevs > 1)
                    {
                        // 3D variable
                        if (hasTileDim)
                        {
                            int dimids[] = {timeDimId, tileDimId, zDimId, yDimId, xDimId};
                            checkNetCDF(nc_def_var(ncid, ncVarName.c_str(), NC_DOUBLE,
                                                   5, dimids, &varid),
                                        "defining variable " + ncVarName);
                        }
                        else
                        {
                            int dimids[] = {timeDimId, zDimId, yDimId, xDimId};
                            checkNetCDF(nc_def_var(ncid, ncVarName.c_str(), NC_DOUBLE,
                                                   4, dimids, &varid),
                                        "defining variable " + ncVarName);
                        }
                    }
                    else
                    {
                        // 2D variable
                        if (hasTileDim)
                        {
                            int dimids[] = {timeDimId, tileDimId, yDimId, xDimId};
                            checkNetCDF(nc_def_var(ncid, ncVarName.c_str(), NC_DOUBLE,
                                                   4, dimids, &varid),
                                        "defining variable " + ncVarName);
                        }
                        else
                        {
                            int dimids[] = {timeDimId, yDimId, xDimId};
                            checkNetCDF(nc_def_var(ncid, ncVarName.c_str(), NC_DOUBLE,
                                                   3, dimids, &varid),
                                        "defining variable " + ncVarName);
                        }
                    }
                }

                checkNetCDF(nc_enddef(ncid), "ending define mode for " + filepath);
            }

            // Write each field
            for (const size_t ifield : fieldsForFile)
            {
                const std::string &jediName = jediNames[ifield];
                const std::string ncVarName = fileionames.getString(jediName);
                const atlas::Field &field = fieldSet.field(jediName);
                auto view = atlas::array::make_view<double, 2>(field);
                const size_t nlevs = field.levels();
                const bool is3d = (nlevs > 1);
                const size_t nzWrite = is3d ? nz : 1;

                // Build local contribution to the global buffer
                const size_t globalSize = ntiles * nzWrite * nxy;
                std::vector<double> localBuf(globalSize, 0.0);

                for (size_t jnode = 0; jnode < numNodes; ++jnode)
                {
                    if (ghostView(jnode))
                        continue;  // skip ghost/halo nodes

                    const size_t gid0 = static_cast<size_t>(globalIdx(jnode)) - 1;
                    const size_t tile0 = gid0 / nxy;
                    const size_t spatialIdx = gid0 % nxy;

                    if (is3d)
                    {
                        for (size_t z = 0; z < nz; ++z)
                        {
                            const size_t bufIdx = tile0 * (nz * nxy) + z * nxy + spatialIdx;
                            localBuf[bufIdx] = view(jnode, z);
                        }
                    }
                    else
                    {
                        localBuf[gid0] = view(jnode, 0);
                    }
                }

                // Reduce to rank 0 (sum — each owned point contributed by exactly one rank)
                std::vector<double> globalBuf(globalSize, 0.0);
                comm.reduce(localBuf.data(), globalBuf.data(), globalSize,
                            eckit::mpi::sum(), 0);

                // Apply inverse scaling if provided
                if (myRank == 0 && fileioscaling.has(jediName))
                {
                    const double scale = fileioscaling.getDouble(jediName);
                    if (scale != 0.0)
                    {
                        const double invScale = 1.0 / scale;
                        for (size_t i = 0; i < globalSize; ++i)
                        {
                            globalBuf[i] *= invScale;
                        }
                    }
                }

                // Rank 0 writes the data
                if (myRank == 0)
                {
                    int varid;
                    checkNetCDF(nc_inq_varid(ncid, ncVarName.c_str(), &varid),
                                "finding variable " + ncVarName + " for writing");

                    std::vector<size_t> start, count;
                    if (is3d)
                    {
                        if (hasTileDim)
                        {
                            start = {0, 0, 0, 0, 0};
                            count = {1, ntiles, nz, ny, nx};
                        }
                        else
                        {
                            start = {0, 0, 0, 0};
                            count = {1, nz, ny, nx};
                        }
                    }
                    else
                    {
                        if (hasTileDim)
                        {
                            start = {0, 0, 0, 0};
                            count = {1, ntiles, ny, nx};
                        }
                        else
                        {
                            start = {0, 0, 0};
                            count = {1, ny, nx};
                        }
                    }

                    checkNetCDF(nc_put_vara_double(ncid, varid,
                                                   start.data(), count.data(),
                                                   globalBuf.data()),
                                "writing " + ncVarName);

                    oops::Log::info() << classname() << " wrote "
                                      << (is3d ? "3D" : "2D") << " field: " << jediName
                                      << " (" << ncVarName << ") to " << filepath
                                      << std::endl;
                }

                fieldWritten[ifield] = true;
            }

            // Close the file on rank 0
            if (myRank == 0)
            {
                checkNetCDF(nc_close(ncid), "closing " + filepath);
            }
        }

        // Check that all fields were written
        for (size_t ifield = 0; ifield < jediNames.size(); ++ifield)
        {
            if (!fieldWritten[ifield])
            {
                oops::Log::warning() << classname() << "::writeHistoryFiles: "
                                     << "Variable \"" << fileionames.getString(jediNames[ifield])
                                     << "\" (JEDI name: " << jediNames[ifield]
                                     << ") was not written to any file" << std::endl;
            }
        }

        oops::Log::info() << classname() << " finished writing "
                          << jediNames.size() << " fields to "
                          << filepaths.size() << " history file(s)" << std::endl;
    }
    // -------------------------------------------------------------------------------------------------
    void IoFV3::print(std::ostream &os) const
    {
        os << classname() << " IO for Cube Sphere History and Restart files";
    }

    // -------------------------------------------------------------------------------------------------
    void IoFV3::checkNetCDF(int status, const std::string &operation) const
    {
        if (status != NC_NOERR)
        {
            throw eckit::Exception("NetCDF error in " + operation + ": " + nc_strerror(status));
        }
    }
    // -------------------------------------------------------------------------------------------------
}  // namespace ijedi
#include <ostream>
#include <string>
#include <vector>

#include <netcdf.h>
#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
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
        // CALL WRITE
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
                    continue; // variable not in this file, try next file
                }

                // Determine number of dimensions
                int ndims;
                checkNetCDF(nc_inq_varndims(ncid, varid, &ndims),
                            "getting ndims for " + ncVarName);

                // Variables are either:
                //   5D: (time, tile, z, y, x) → atlas field shape (ntiles, nz, ny, nx)
                //   4D: (time, tile, y, x)    → atlas field shape (ntiles, ny, nx)
                //   When tile is not a dimension, reduce by 1D each

                const int expected3d = hasTileDim ? 5 : 4; // with z
                const int expected2d = hasTileDim ? 4 : 3; // without z

                if (ndims == expected3d)
                {
                    // 3D field (has vertical levels)
                    std::vector<size_t> shape = {ntiles, nz, ny, nx};
                    atlas::Field field(jediName,
                                       atlas::array::make_datatype<double>(),
                                       atlas::array::ArrayShape(shape));
                    auto view = atlas::array::make_view<double, 4>(field);

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
                                                   view.data()),
                                "reading " + ncVarName);

                    // Apply scaling if provided
                    if (fileioscaling.has(jediName))
                    {
                        const double scale = fileioscaling.getDouble(jediName);
                        const size_t total = ntiles * nz * ny * nx;
                        double *raw = view.data();
                        for (size_t i = 0; i < total; ++i)
                        {
                            raw[i] *= scale;
                        }
                    }

                    fieldSet.add(field);
                    fieldRead[ifield] = true;
                    oops::Log::info() << classname() << " read 3D field: " << jediName
                                      << " (" << ncVarName << ") from " << filepath
                                      << std::endl;
                }
                else if (ndims == expected2d)
                {
                    // 2D field (no vertical levels)
                    std::vector<size_t> shape = {ntiles, ny, nx};
                    atlas::Field field(jediName,
                                       atlas::array::make_datatype<double>(),
                                       atlas::array::ArrayShape(shape));
                    auto view = atlas::array::make_view<double, 3>(field);

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
                                                   view.data()),
                                "reading " + ncVarName);

                    // Apply scaling if provided
                    if (fileioscaling.has(jediName))
                    {
                        const double scale = fileioscaling.getDouble(jediName);
                        const size_t total = ntiles * ny * nx;
                        double *raw = view.data();
                        for (size_t i = 0; i < total; ++i)
                        {
                            raw[i] *= scale;
                        }
                    }

                    fieldSet.add(field);
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
} // namespace ijedi
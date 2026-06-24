#pragma once

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

// Forward declarations
namespace eckit
{
    class Configuration;
}

namespace util
{
    class DateTime;
    class Duration;
}

namespace ijedi
{
    extern "C"
    {
        void f_fv3_geom_initialize(const eckit::LocalConfiguration &, const eckit::mpi::Comm *);
        void f_fv3_geom_create(const eckit::Configuration &, const eckit::Configuration &,
                               const eckit::mpi::Comm *);
        void f_fv3_geom_set_and_fill_geometry_fields(const void *, const void *, const char *, int, int,
                                                     const double *, const double *, const double *,
                                                     const double *);
    }  // extern "C"
}  // namespace ijedi

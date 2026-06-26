#pragma once

#include "atlas/field.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{
extern "C"
{
  void ijedi_io_fv3_restart_read_f90(const eckit::Configuration &, const eckit::Configuration &,
                                     atlas::field::FieldSetImpl *,
                                     const eckit::Configuration &, const eckit::Configuration &);
  void ijedi_io_fv3_restart_write_f90(const eckit::Configuration &, const eckit::Configuration &,
                                      const atlas::field::FieldSetImpl *,
                                      const eckit::Configuration &, const eckit::Configuration &);
}
}  // namespace ijedi

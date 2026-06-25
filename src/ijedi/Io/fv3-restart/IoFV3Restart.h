#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "ijedi/Geometry/Geometry.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ijedi/Io/IoBase.h"

namespace ijedi
{

  // -------------------------------------------------------------------------------------------------

  class IoFV3RestartParameters : public IoParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(IoFV3RestartParameters, IoParametersBase)

   public:
      // Are files to be read restarts or not?
      oops::Parameter<bool> is_restart{"is restart", "is restart", true, this};

      // Data path for files being read
      oops::Parameter<std::string> datapath{"datapath",
                                            "path to location of files to be read",
                                            "./", this};

      // Filename to be read or written
      oops::Parameter<std::string> filename_nonrestart{"filename_nonrestart",
                                                       "filename_nonrestart",
                                                       "fms_nonrestart.nc", this};

      // Restart filenames to be read
      oops::Parameter<std::string> filename_core{"filename_core",
                                                 "filename_core",
                                                 "fv_core.res.nc", this};
      oops::Parameter<std::string> filename_trcr{"filename_trcr",
                                                 "filename_trcr",
                                                 "fv_tracer.res.nc", this};
      oops::Parameter<std::string> filename_sfcd{"filename_sfcd",
                                                 "filename_sfcd",
                                                 "sfc_data.nc", this};
      oops::Parameter<std::string> filename_sfcw{"filename_sfcw",
                                                 "filename_sfcw",
                                                 "fv_srf_wnd.res.nc", this};
      oops::Parameter<std::string> filename_cplr{"filename_cplr",
                                                 "filename_cplr",
                                                 "coupler.res", this};
      oops::Parameter<std::string> filename_spec{"filename_spec",
                                                 "filename_spec",
                                                 "null", this};
      oops::Parameter<std::string> filename_phys{"filename_phys",
                                                 "filename_phys",
                                                 "phy_data.nc", this};
      oops::Parameter<std::string> filename_orog{"filename_orog",
                                                 "filename_orog",
                                                 "oro_data.nc", this};
      oops::Parameter<std::string> filename_cold{"filename_cold",
                                                 "filename_cold",
                                                 "gfs_data.nc", this};

      // Input filename may be templated with datetimes
      oops::Parameter<bool> filename_is_datetime_templated{
        "filename is datetime templated",
        "filename is datetime templated",
        false, this};

      // Skip reading/writing the coupler.res file
      oops::Parameter<bool> skip_coupler_file{"skip coupler file",
                                              "skip coupler file",
                                              false, this};

      // Prepend the files with the date
      oops::Parameter<bool> prepend_files_with_date{"prepend files with date",
                                                    "prepend files with date",
                                                    true, this};

      // Force tell the system that surface pressure is in the file
      oops::Parameter<bool> psinfile{"psinfile",
                                     "tell the system surface pressure is in the file",
                                     false, this};

      // Optionally the config may contain member
      oops::OptionalParameter<int> member{"member", "ensemble member number", this};

      // Let user set the calendar type
      oops::Parameter<int> calendar_type{"calendar type", "calendar type", 2, this};

      // Ignore checksum for FMS restarts?
      oops::Parameter<bool> ignore_checksum{"ignore checksum",
                                            "whether to ignore restart checksums",
                                            true, this};

      // Write only a subset of fields?
      oops::OptionalParameter<std::vector<std::string>> fields_to_write{
        "fields to write",
        "names of fields to write",
        this};
  };

  // -------------------------------------------------------------------------------------------------
  class IoFV3Restart : public IoBase, private util::ObjectCounter<IoFV3Restart>
  {
   public:
    static const std::string classname() { return "ijedi::IoFV3Restart"; }

    typedef IoFV3RestartParameters Parameters_;

    IoFV3Restart(const Geometry &, const Parameters_ &);
    ~IoFV3Restart();
    void read(atlas::FieldSet &, const eckit::LocalConfiguration &,
              const eckit::LocalConfiguration &) const override;
    void write(const atlas::FieldSet &, const eckit::LocalConfiguration &,
               const eckit::LocalConfiguration &) const override;

   private:
    void print(std::ostream &) const override;

    const Geometry & geom_;
    Parameters_ parameters_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace ijedi

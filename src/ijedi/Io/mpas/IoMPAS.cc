#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/functionspace/NodeColumns.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "ijedi/Geometry/Geometry.h"
#include "ijedi/Geometry/mpas/GeometryMPAS.h"
#include "ijedi/Io/mpas/IoMPAS.h"
#include "ijedi/Io/mpas/IoMPAS.interface.h"

namespace {

void readMPASStream(const std::string & filepath,
                    atlas::FieldSet & x,
                    const std::vector<std::string> & fileVarNames,
                    const std::vector<double> & scalings,
                    const void * fortranGeom,
                    const std::string & streamName) {
  const auto fs = atlas::functionspace::NodeColumns(x.field(0).functionspace());
  auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(fs.nodes().global_index());
  const int nLocal = static_cast<int>(fs.nodes().size());

  std::vector<int> globalIndices(nLocal, 0);
  for (int n = 0; n < nLocal; ++n) {
    globalIndices[n] = static_cast<int>(localGidx(n));
  }

  int fi = 0;
  for (auto & field : x) {
    const std::string & fileVarName = fileVarNames.at(fi);
    const double scale = scalings.at(fi);
    ++fi;

    if (!fileVarName.empty()) {
      const int nLevelsField = field.shape(1);
      std::vector<double> values(static_cast<size_t>(nLocal) *
                                 static_cast<size_t>(nLevelsField), 0.0);

      ijedi::ijedi_mpas_state_read_field_f90(fortranGeom,
                                             filepath.c_str(),
                                             streamName.c_str(),
                                             fileVarName.c_str(),
                                             nLocal,
                                             globalIndices.data(),
                                             nLevelsField,
                                             values.data(),
                                             scale);

      auto view = atlas::array::make_view<double, 2>(field);
      for (int n = 0; n < nLocal; ++n) {
        for (int k = 0; k < nLevelsField; ++k) {
          view(n, k) = values[static_cast<size_t>(n) * nLevelsField + k];
        }
      }
    }

    field.metadata().set("interp_type", "default");
  }

  for (auto & field : x) {
    atlas::functionspace::NodeColumns(field.functionspace()).haloExchange(field);
  }
}

void writeMPASStream(const std::string & filepath,
                     const atlas::FieldSet & x,
                     const std::vector<std::string> & fileVarNames,
                     const void * fortranGeom,
                     const std::string & streamName) {
  const auto fs = atlas::functionspace::NodeColumns(x.field(0).functionspace());
  auto localGidx = atlas::array::make_view<atlas::gidx_t, 1>(fs.nodes().global_index());
  const int nLocal = static_cast<int>(fs.nodes().size());

  std::vector<int> globalIndices(nLocal, 0);
  for (int n = 0; n < nLocal; ++n) {
    globalIndices[n] = static_cast<int>(localGidx(n));
  }

  int fi = 0;
  for (const auto & field : x) {
    const std::string & varName = fileVarNames.at(fi++);
    if (varName.empty()) {
      continue;
    }

    const int nLevelsField = field.shape(1);
    std::vector<double> values(static_cast<size_t>(nLocal) *
                               static_cast<size_t>(nLevelsField), 0.0);

    auto view = atlas::array::make_view<double, 2>(field);
    for (int n = 0; n < nLocal; ++n) {
      for (int k = 0; k < nLevelsField; ++k) {
        values[static_cast<size_t>(n) * nLevelsField + k] = view(n, k);
      }
    }

    ijedi::ijedi_mpas_state_write_field_f90(fortranGeom,
                                            filepath.c_str(),
                                            streamName.c_str(),
                                            varName.c_str(),
                                            nLocal,
                                            globalIndices.data(),
                                            nLevelsField,
                                            values.data());
  }

  ijedi::ijedi_mpas_state_write_flush_f90(fortranGeom, filepath.c_str(), streamName.c_str());
}

}  // namespace

namespace ijedi
{

static IoMaker<IoMPAS> makerIoMPAS_("mpas");

IoMPAS::IoMPAS(const Geometry &geom, const Parameters_ &params)
    : IoBase(geom, params.toConfiguration()),
      geom_(geom),
      datapath_(params.datapath),
      filename_(params.filename.value().value_or("")),
      streamName_(params.streamName.value().value_or(""))
{
  util::Timer timer(classname(), "IoMPAS");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}

IoMPAS::~IoMPAS()
{
  util::Timer timer(classname(), "~IoMPAS");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

void IoMPAS::read(atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                  const eckit::LocalConfiguration &fileioscaling) const
{
  util::Timer timer(classname(), "read state");
  oops::Log::trace() << classname() << " read state starting" << std::endl;

  const std::string filepath = (!filename_.empty() && filename_[0] == '/')
                                   ? filename_
                                   : datapath_ + "/" + filename_;

  const auto & geomMPAS = dynamic_cast<const GeometryMPAS &>(geom_.geometryImpl());

  std::vector<std::string> fileVarNames;
  std::vector<double> scalings;
  fileVarNames.reserve(x.size());
  scalings.reserve(x.size());

  for (const auto & field : x)
  {
    const std::string jediName = field.name();
    fileVarNames.push_back(fileionames.has(jediName)
                               ? fileionames.getString(jediName)
                               : jediName);
    scalings.push_back(fileioscaling.has(jediName)
                           ? fileioscaling.getDouble(jediName)
                           : 0.0);
  }

  std::string readStream = streamName_;
  if (readStream.empty()) {
    readStream = "restart";
  }

  readMPASStream(filepath, x, fileVarNames, scalings, geomMPAS.fortranGeom(), readStream);

  oops::Log::trace() << classname() << " read state done" << std::endl;
}

void IoMPAS::write(const atlas::FieldSet &x, const eckit::LocalConfiguration &fileionames,
                   const eckit::LocalConfiguration &fileioscaling) const
{
  util::Timer timer(classname(), "write state");
  oops::Log::trace() << classname() << " write state starting" << std::endl;

  const std::string filepath = (!filename_.empty() && filename_[0] == '/')
                                   ? filename_
                                   : datapath_ + "/" + filename_;

  const auto & geomMPAS = dynamic_cast<const GeometryMPAS &>(geom_.geometryImpl());
  (void)fileioscaling;

  std::vector<std::string> fileVarNames;
  fileVarNames.reserve(x.size());
  for (const auto & field : x)
  {
    const std::string jediName = field.name();
    fileVarNames.push_back(fileionames.has(jediName)
                               ? fileionames.getString(jediName)
                               : jediName);
  }

  std::string writeStream = streamName_;
  if (writeStream.empty()) {
    writeStream = "da_state";
  }

  writeMPASStream(filepath, x, fileVarNames, geomMPAS.fortranGeom(), writeStream);

  oops::Log::trace() << classname() << " write state done" << std::endl;
}

void IoMPAS::print(std::ostream &os) const
{
  os << classname() << " Io for MPAS restart files";
}

}  // namespace ijedi

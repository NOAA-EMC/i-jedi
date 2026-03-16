// (C) Copyright 2026- NOAA.
// This software is licensed under the terms of the Creative Commons
// Attribution-NonCommercial-ShareAlike Licence.
// See LICENSE file in the top-level directory for details.

#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "ijedi/Geometry/base/GeometryBase.h"

namespace eckit
{
  class Configuration;
}

namespace ijedi
{

  class GeometryMOM6 : public GeometryBase
  {
  public:
    GeometryMOM6(const eckit::Configuration &, const eckit::mpi::Comm &);
    void print(std::ostream &) const override;
    eckit::LocalConfiguration gridSpecific() const override;

    // Accessors for MOM6-specific (structured) decomposition.
    // Used by ModelMOM6 for state injection — not exposed via GeometryBase.
    const atlas::FunctionSpace & mom6FunctionSpace() const { return mom6FunctionSpace_; }
    const atlas::FieldSet      & mom6Fields()        const { return mom6Fields_; }

    // Scatter/gather map: for each JEDI-owned ocean point, the owning MOM6 rank
    // and the row-major local index on that rank.  Built once at init; used by
    // ModelMOM6 to exchange state before/after MOM6::run().
    struct ScatterMap {
      std::vector<int> mom6Rank;      // [0..ownedCount_) → MOM6 rank
      std::vector<int> mom6LocalIdx;  // [0..ownedCount_) → row-major local idx on that rank
    };
    const ScatterMap & scatterMap() const { return scatterMap_; }

  private:
    // Replicates FMS compute_extent() — integer-division domain partition
    static int computeExtent(int N, int ndivs, int pe);

    // --- Grid read helpers (always return full global arrays) ---
    // lonT/latT: T-cell centres, size [njGlobal * niGlobal], row-major (j, i)
    void readHgridLonLat(const std::string & path,
                         std::vector<double> & lon,
                         std::vector<double> & lat) const;
    // dxT, dyT, areaT: T-cell metrics, same layout
    void readHgridMetrics(const std::string & path,
                          std::vector<double> & dxT,
                          std::vector<double> & dyT,
                          std::vector<double> & areaT) const;
    // lonU/latU: east-face (U) centres; lonV/latV: north-face (V) centres; same layout
    void readHgridUVLonLat(const std::string & path,
                           std::vector<double> & lonU,
                           std::vector<double> & latU,
                           std::vector<double> & lonV,
                           std::vector<double> & latV) const;
    // depth(j, i) and wet(j, i) — MOM6's authoritative land/sea mask, size [njGlobal * niGlobal]
    void readTopog(const std::string & path,
                   std::vector<double> & depth,
                   std::vector<double> & wet) const;

    // --- Function space builders ---
    // MOM6 structured space: replicates the FMS rectangular tile decomposition → mom6FunctionSpace_
    void buildMom6FunctionSpace(const eckit::mpi::Comm & comm,
                                 const std::vector<double> & lonGlobal,
                                 const std::vector<double> & latGlobal);

    // JEDI unstructured space: ocean + land-fringe points, load-balanced → functionSpace_
    void buildJediFunctionSpace(const eckit::mpi::Comm & comm,
                                 const std::vector<double> & wetGlobal,
                                 const std::vector<double> & lonGlobal,
                                 const std::vector<double> & latGlobal);

    // --- Field builders ---
    // Populate mom6Fields_ on mom6FunctionSpace_ from global arrays (extracts local slice)
    void buildMom6Fields(const std::vector<double> & lonGlobal,
                          const std::vector<double> & latGlobal,
                          const std::vector<double> & depthGlobal,
                          const std::vector<double> & wetGlobal,
                          const std::vector<double> & dxTGlobal,
                          const std::vector<double> & dyTGlobal,
                          const std::vector<double> & areaTGlobal,
                          const std::vector<double> & lonUGlobal,
                          const std::vector<double> & latUGlobal,
                          const std::vector<double> & lonVGlobal,
                          const std::vector<double> & latVGlobal);

    // Populate fields_ on functionSpace_ (JEDI unstructured) from global arrays via jediPoints_
    void buildFields(const std::vector<double> & lonGlobal,
                     const std::vector<double> & latGlobal,
                     const std::vector<double> & depthGlobal,
                     const std::vector<double> & wetGlobal,
                     const std::vector<double> & dxTGlobal,
                     const std::vector<double> & dyTGlobal,
                     const std::vector<double> & areaTGlobal,
                     const std::vector<double> & lonUGlobal,
                     const std::vector<double> & latUGlobal,
                     const std::vector<double> & lonVGlobal,
                     const std::vector<double> & latVGlobal);

    // --- Scatter map ---
    void buildScatterMap();

    // --- Output helpers ---
    void saveGrid(const std::string & filename,
                  const eckit::mpi::Comm & comm) const;
    void saveStructuredGrid(const std::string & filename,
                            const eckit::mpi::Comm & comm) const;
    void saveDebugMesh(const std::string & prefix,
                       const eckit::mpi::Comm & comm) const;

    // Grid parameters (from YAML config)
    int niGlobal_ = 0;
    int njGlobal_ = 0;
    int layoutX_  = 1;
    int layoutY_  = 1;
    double minimumDepth_ = 0.0;

    // Local MOM6 compute domain (1-based, inclusive)
    int iStart_ = 1, iCount_ = 0;
    int jStart_ = 1, jCount_ = 0;

    // MOM6 structured function space and fields (local compute domain)
    atlas::FunctionSpace mom6FunctionSpace_;
    atlas::FieldSet      mom6Fields_;

    // JEDI unstructured partition
    int nActiveGlobal_ = 0;   // global active point count (ocean + fringe)
    int ownedCount_    = 0;   // number of active points owned by this rank
    int ownedMin_      = 0;   // min owned count across ranks (load balance)
    int ownedMax_      = 0;   // max owned count across ranks (load balance)

    // (iGlobal, jGlobal) for each JEDI node in functionSpace_ order (owned first, then ghost)
    std::vector<std::pair<int,int>> jediPoints_;

    // Scatter/gather map: JEDI unstructured ↔ MOM6 structured
    ScatterMap scatterMap_;
  };

} // namespace ijedi

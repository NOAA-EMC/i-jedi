/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

// Vendored from saber/src/saber/interpolation/Geometry.cc. Keep in sync with upstream.

#include "ijedi/Geometry/gsibec/GsiGrid.h"

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "atlas/grid/detail/spacing/gaussian/Latitudes.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

namespace ijedi
{
  const std::string GsiGridKey = "custom grid matching gsi";
  const std::string GsiPartitionerKey = "custom partitioner matching gsi";

  std::vector<int> computeS2NCheckerboardPartition(const atlas::RegularGrid & rg,
                                                   const int ntasks, const int nbands) {
    // Number of MPI tasks (=partitions) per band
    if (ntasks % nbands != 0) {
      atlas::throw_Exception("number of bands doesn't divide number of tasks", Here());
    }

    const auto map_1d_point_to_1d_partition = [&](const int npoints,
                                                  const int npartitions) -> std::vector<int> {
      const int npart = npoints / npartitions;
      const int nremain = npoints % npartitions;
      std::vector<int> mapping(npoints);
      int i = 0;
      for (int p = 0; p < npartitions; ++p) {
        const int npoints_this_part = npart + (p < nremain ? 1 : 0);
        for (int count = 0; count < npoints_this_part; ++count) {
          mapping[i] = p;
          ++i;
        }
      }
      ASSERT(i == npoints);
      return mapping;
    };

    const int nx = rg.nx();
    const int ny = rg.ny();
    const size_t parts_per_band = ntasks / nbands;

    const auto band = map_1d_point_to_1d_partition(ny, nbands);
    const auto part_in_band = map_1d_point_to_1d_partition(nx, parts_per_band);

    std::vector<int> partition(rg.size());
    for (atlas::idx_t j = 0; j < ny; ++j) {
      for (atlas::idx_t i = 0; i < nx; ++i) {
        partition[rg.index(i, j)] = band[j] * parts_per_band + part_in_band[i];
      }
    }
    return partition;
  }

  void setupGsiMatchingGrid(const eckit::Configuration & config,
                            const eckit::mpi::Comm & comm,
                            atlas::Grid & grid,
                            atlas::FunctionSpace & functionSpace,
                            atlas::FieldSet & fieldSet) {
    const std::string grid_type = config.getString(GsiGridKey + ".type");
    ASSERT(grid_type == "gaussian" || grid_type == "latlon" || grid_type == "rotated_lonlat");


    const auto require_parameter = [&](const std::string & name) {
      const std::string key = GsiGridKey + "." + name;
      if (!config.has(key)) {
        throw eckit::BadParameter(
          "Missing required parameter \"" + key + "\" for rotated_lonlat GSI-matching grid");
      }
    };

    if (grid_type == "rotated_lonlat") {
      require_parameter("lat_start");
      require_parameter("lat_end");
      require_parameter("lon_start");
      require_parameter("lon_end");
      require_parameter("north_pole_lat");
      require_parameter("north_pole_lon");
    }



    const int nlats = config.getInt(GsiGridKey + ".lats");  // pole to pole
    const int nlons = config.getInt(GsiGridKey + ".lons");
    const bool lendp = config.getInt(GsiGridKey + ".endpoint");
    const double lat_start = config.has(GsiGridKey + ".lat_start") ?
                             config.getDouble(GsiGridKey + ".lat_start") : 0.0;
    const double lat_end = config.has(GsiGridKey + ".lat_end") ?
                           config.getDouble(GsiGridKey + ".lat_end") : 0.0;
    const double lon_start = config.has(GsiGridKey + ".lon_start") ?
                             config.getDouble(GsiGridKey + ".lon_start") : 0.0;
    const double lon_end = config.has(GsiGridKey + ".lon_end") ?
                           config.getDouble(GsiGridKey + ".lon_end") : 0.0;
    const double north_pole_lat = config.has(GsiGridKey + ".north_pole_lat") ?
                                  config.getDouble(GsiGridKey + ".north_pole_lat") : 0.0;
    const double north_pole_lon = config.has(GsiGridKey + ".north_pole_lon") ?
                                  config.getDouble(GsiGridKey + ".north_pole_lon") : 0.0;

    
    const auto gsi_gaussian_points = [](const int N) -> std::vector<double> {
      ASSERT(N % 2 == 0);  // code below would need verification, probably fixing, in odd case
      std::vector<double> result(N);
      // north-to-south order, following atlas's default
      result[0] = 90.0;
      atlas::grid::spacing::gaussian::gaussian_latitudes_npole_spole((N-2)/2, result.data()+1);
      result[N-1] = -90.0;
      // flip sign to obtain south-to-north order, following GSI's default
      for (auto & r : result) {
        r *= 1.0; //-1.0;  keep north-to-south order for surface 
      }
      return result;
    };
    
    std::vector<double> grid_xt_ar = {
    0, 0.9375, 1.875, 2.8125, 3.75, 4.6875, 5.625, 6.5625, 7.5,
    8.4375, 9.375, 10.3125, 11.25, 12.1875, 13.125, 14.0625, 15, 15.9375,
    16.875, 17.8125, 18.75, 19.6875, 20.625, 21.5625, 22.5, 23.4375, 24.375,
    25.3125, 26.25, 27.1875, 28.125, 29.0625, 30, 30.9375, 31.875, 32.8125,
    33.75, 34.6875, 35.625, 36.5625, 37.5, 38.4375, 39.375, 40.3125, 41.25,
    42.1875, 43.125, 44.0625, 45, 45.9375, 46.875, 47.8125, 48.75, 49.6875,
    50.625, 51.5625, 52.5, 53.4375, 54.375, 55.3125, 56.25, 57.1875, 58.125,
    59.0625, 60, 60.9375, 61.875, 62.8125, 63.75, 64.6875, 65.625, 66.5625,
    67.5, 68.4375, 69.375, 70.3125, 71.25, 72.1875, 73.125, 74.0625, 75,
    75.9375, 76.875, 77.8125, 78.75, 79.6875, 80.625, 81.5625, 82.5, 83.4375,
    84.375, 85.3125, 86.25, 87.1875, 88.125, 89.0625, 90, 90.9375, 91.875,
    92.8125, 93.75, 94.6875, 95.625, 96.5625, 97.5, 98.4375, 99.375,
    100.3125, 101.25, 102.1875, 103.125, 104.0625, 105, 105.9375, 106.875,
    107.8125, 108.75, 109.6875, 110.625, 111.5625, 112.5, 113.4375, 114.375,
    115.3125, 116.25, 117.1875, 118.125, 119.0625, 120, 120.9375, 121.875,
    122.8125, 123.75, 124.6875, 125.625, 126.5625, 127.5, 128.4375, 129.375,
    130.3125, 131.25, 132.1875, 133.125, 134.0625, 135, 135.9375, 136.875,
    137.8125, 138.75, 139.6875, 140.625, 141.5625, 142.5, 143.4375, 144.375,
    145.3125, 146.25, 147.1875, 148.125, 149.0625, 150, 150.9375, 151.875,
    152.8125, 153.75, 154.6875, 155.625, 156.5625, 157.5, 158.4375, 159.375,
    160.3125, 161.25, 162.1875, 163.125, 164.0625, 165, 165.9375, 166.875,
    167.8125, 168.75, 169.6875, 170.625, 171.5625, 172.5, 173.4375, 174.375,
    175.3125, 176.25, 177.1875, 178.125, 179.0625, 180, 180.9375, 181.875,
    182.8125, 183.75, 184.6875, 185.625, 186.5625, 187.5, 188.4375, 189.375,
    190.3125, 191.25, 192.1875, 193.125, 194.0625, 195, 195.9375, 196.875,
    197.8125, 198.75, 199.6875, 200.625, 201.5625, 202.5, 203.4375, 204.375,
    205.3125, 206.25, 207.1875, 208.125, 209.0625, 210, 210.9375, 211.875,
    212.8125, 213.75, 214.6875, 215.625, 216.5625, 217.5, 218.4375, 219.375,
    220.3125, 221.25, 222.1875, 223.125, 224.0625, 225, 225.9375, 226.875,
    227.8125, 228.75, 229.6875, 230.625, 231.5625, 232.5, 233.4375, 234.375,
    235.3125, 236.25, 237.1875, 238.125, 239.0625, 240, 240.9375, 241.875,
    242.8125, 243.75, 244.6875, 245.625, 246.5625, 247.5, 248.4375, 249.375,
    250.3125, 251.25, 252.1875, 253.125, 254.0625, 255, 255.9375, 256.875,
    257.8125, 258.75, 259.6875, 260.625, 261.5625, 262.5, 263.4375, 264.375,
    265.3125, 266.25, 267.1875, 268.125, 269.0625, 270, 270.9375, 271.875,
    272.8125, 273.75, 274.6875, 275.625, 276.5625, 277.5, 278.4375, 279.375,
    280.3125, 281.25, 282.1875, 283.125, 284.0625, 285, 285.9375, 286.875,
    287.8125, 288.75, 289.6875, 290.625, 291.5625, 292.5, 293.4375, 294.375,
    295.3125, 296.25, 297.1875, 298.125, 299.0625, 300, 300.9375, 301.875,
    302.8125, 303.75, 304.6875, 305.625, 306.5625, 307.5, 308.4375, 309.375,
    310.3125, 311.25, 312.1875, 313.125, 314.0625, 315, 315.9375, 316.875,
    317.8125, 318.75, 319.6875, 320.625, 321.5625, 322.5, 323.4375, 324.375,
    325.3125, 326.25, 327.1875, 328.125, 329.0625, 330, 330.9375, 331.875,
    332.8125, 333.75, 334.6875, 335.625, 336.5625, 337.5, 338.4375, 339.375,
    340.3125, 341.25, 342.1875, 343.125, 344.0625, 345, 345.9375, 346.875,
    347.8125, 348.75, 349.6875, 350.625, 351.5625, 352.5, 353.4375, 354.375,
    355.3125, 356.25, 357.1875, 358.125, 359.0625 } ;

    const auto build_xspace_config = [&](const std::string & grid_type)
                                     -> eckit::LocalConfiguration {
      eckit::LocalConfiguration lc{};
      if (grid_type == "rotated_lonlat") {
        lc.set("type", "linear");
        lc.set("N", nlons);
        lc.set("start", lon_start);
        lc.set("end", lon_end);
      /*} else if (grid_type == "gaussian") {
        lc.set("type", "custom");
        lc.set("N", nlons);
        lc.set("values", grid_xt_ar); */
      } else {
        lc.set("type", "linear");
        lc.set("N", nlons);
        lc.set("interval", std::vector<double>{{0.0, 360.0}});   //359.0625}}); 
        lc.set("endpoint", lendp);
      }
      return lc;
    };
    
    std::vector<double> grid_yt_ar = {
    89.2842275325136, 88.3570035186649, 87.4243037460699,
    86.4903667662812, 85.5559604848927, 84.6213271076488, 83.6865668165639,
    82.7517284734307, 81.8168387286032, 80.8819133467975, 79.9469622473857,
    79.0119919826722, 78.0770070543037, 77.1420106570536, 76.2070051208607,
    75.2719921848602, 74.336973173452, 73.4019491129572, 72.4669208110036,
    71.5318889118272, 70.5968539356029, 69.6618163069383, 68.7267763758595,
    67.7917344335006, 66.8566907239931, 65.9216454535868, 64.9865987977267,
    64.0515509066005, 63.1165019095293, 62.1814519184735, 61.2464010308544,
    60.3113493318435, 59.3762968962322, 58.4412437899697, 57.5061900714349,
    56.5711357924953, 55.6360809993933, 54.7010257334911, 53.765970031901,
    52.8309139280201, 51.8958574519866, 50.9608006310702, 50.0257434900084,
    49.0906860512962, 48.1556283354368, 47.2205703611601, 46.2855121456128,
    45.3504537045245, 44.4153950523538, 43.4803362024163, 42.5452771669981,
    41.6102179574557, 40.6751585843046, 39.7400990572982, 38.8050393854982,
    37.8699795773369, 36.934919640674, 35.999859582846, 35.064799410712,
    34.1297391306945, 33.1946787488156, 32.2596182707304, 31.3245577017573,
    30.3894970469048, 29.4544363108967, 28.5193754981943, 27.5843146130172,
    26.649253659362, 25.7141926410195, 24.7791315615903, 23.8440704244996,
    22.9090092330105, 21.9739479902356, 21.0388866991495, 20.1038253625982,
    19.1687639833093, 18.2337025639012, 17.2986411068909, 16.3635796147026,
    15.4285180896742, 14.4934565340649, 13.5583949500608, 12.6233333397819,
    11.6882717052868, 10.7532100485791, 9.81814837161151, 8.88308667629159,
    7.94802496448583, 7.01296323802434, 6.07790149870505, 5.14283974829787,
    4.20777798854868, 3.27271622118319, 2.33765444791075, 1.40259267042804,
    0.467530890422768, -0.467530890422768, -1.40259267042804,
    -2.33765444791075, -3.27271622118319, -4.20777798854868,
    -5.14283974829787, -6.07790149870505, -7.01296323802434,
    -7.94802496448583, -8.88308667629159, -9.81814837161151,
    -10.7532100485791, -11.6882717052868, -12.6233333397819,
    -13.5583949500608, -14.4934565340649, -15.4285180896742,
    -16.3635796147026, -17.2986411068909, -18.2337025639012,
    -19.1687639833093, -20.1038253625982, -21.0388866991495,
    -21.9739479902356, -22.9090092330105, -23.8440704244996,
    -24.7791315615903, -25.7141926410195, -26.649253659362,
    -27.5843146130172, -28.5193754981943, -29.4544363108967,
    -30.3894970469048, -31.3245577017573, -32.2596182707304,
    -33.1946787488156, -34.1297391306945, -35.064799410712, -35.999859582846,
    -36.934919640674, -37.8699795773369, -38.8050393854982,
    -39.7400990572982, -40.6751585843046, -41.6102179574557,
    -42.5452771669981, -43.4803362024163, -44.4153950523538,
    -45.3504537045245, -46.2855121456128, -47.2205703611601,
    -48.1556283354368, -49.0906860512962, -50.0257434900084,
    -50.9608006310702, -51.8958574519866, -52.8309139280201,
    -53.765970031901, -54.7010257334911, -55.6360809993933,
    -56.5711357924953, -57.5061900714349, -58.4412437899697,
    -59.3762968962322, -60.3113493318435, -61.2464010308544,
    -62.1814519184735, -63.1165019095293, -64.0515509066005,
    -64.9865987977267, -65.9216454535868, -66.8566907239931,
    -67.7917344335006, -68.7267763758595, -69.6618163069383,
    -70.5968539356029, -71.5318889118272, -72.4669208110036,
    -73.4019491129572, -74.336973173452, -75.2719921848602,
    -76.2070051208607, -77.1420106570536, -78.0770070543037,
    -79.0119919826722, -79.9469622473857, -80.8819133467975,
    -81.8168387286032, -82.7517284734307, -83.6865668165639,
    -84.6213271076488, -85.5559604848927, -86.4903667662812,
    -87.4243037460699, -88.3570035186649, -89.2842275325136 };

    const auto build_yspace_config = [&](const std::string & grid_type) ->
                                     eckit::LocalConfiguration {
      eckit::LocalConfiguration lc{};
      if (grid_type == "rotated_lonlat") {
        lc.set("type", "linear");
        lc.set("N", nlats);
        lc.set("start", lat_start);
        lc.set("end", lat_end);
      } else if (grid_type == "gaussian") {
        lc.set("type", "custom");
        lc.set("N", nlats);
        lc.set("values", gsi_gaussian_points(nlats)); //grid_yt_ar); //
      } else {
        lc.set("type", "linear");
        lc.set("N", nlats);
        lc.set("interval", std::vector<double>{{90.0, -90.0}});
      }
      return lc;
    };

    const auto build_projection_config = [&](const std::string & grid_type) ->
                                         eckit::LocalConfiguration {
      eckit::LocalConfiguration lc{};
      lc.set("type", "rotated_lonlat");
      lc.set("north_pole", std::vector<double>{{north_pole_lon, north_pole_lat}});
      return lc;
    };

    eckit::LocalConfiguration testconfig{};
    testconfig.set("type", "structured");
    testconfig.set("xspace", build_xspace_config(grid_type));
    testconfig.set("yspace", build_yspace_config(grid_type));
    if (grid_type == "rotated_lonlat") {
      testconfig.set("projection", build_projection_config(grid_type));
    }
    grid = atlas::Grid{testconfig};

    const atlas::RegularGrid rg{grid};
    ASSERT(rg);

    const int ntasks = comm.size();
    const int nbands = config.getInt(GsiPartitionerKey + ".bands");
    ASSERT(nbands >= 1 && nbands <= ntasks);
    std::vector<int> partition = computeS2NCheckerboardPartition(rg, ntasks, nbands);

    const atlas::grid::Distribution distribution(ntasks, partition.size(), partition.data());
    const unsigned halo = config.getUnsigned("halo", 1);
    functionSpace = atlas::functionspace::StructuredColumns(grid, distribution,
                                                            atlas::option::halo(halo));

    // Using atlas::mpi::Scope in the call to atlas::functionspace::StructuredColumns
    // may have reverted the default communicator to the world communicator.
    // We set back the default communicator to `comm` to fix this.
    eckit::mpi::setCommDefault(comm.name().c_str());

    fieldSet.clear();  // return empty
  }

}  // namespace ijedi

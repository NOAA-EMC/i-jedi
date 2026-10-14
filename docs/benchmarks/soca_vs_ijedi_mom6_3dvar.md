# SOCA vs i-jedi/MOM6: 3DVar walltime and memory benchmark

This document compares the walltime and memory of the same ocean 3DVar run two ways:

- with SOCA (`soca_var.x`);
- with the i-jedi MOM6 interface (`ijedi.x var`).

It gives the results at 1 degree and everything needed to set the test up again on
HPC at higher resolution.

## Summary

At 1 degree with 6 MPI tasks:

- **i-jedi/MOM6 is about 1.3 - 1.5x faster than SOCA.**
- **It uses about 25 % less memory** (11 GB vs 14 GB in total).

The two runs solve the same problem: their background departures are identical.
Almost all of the cost is the SABER diffusion B. That cost scales with the number of
points each interface hands to SABER:

- SOCA carries the whole grid, land included: 115 200 points per level.
- i-jedi carries only wet columns: 77 852 points per level.

## Test definition

| Item | Setting |
|---|---|
| Grid | MOM6 tripolar, 360 x 320 x 75 (1 deg) |
| Background | MOM6 restart `MOM.res.nc` (`Temp`, `Salt`, `h`, `sfc`), 2020-12-16T00Z |
| Window | 2020-12-13T00Z + P6D |
| Analysis variables | `sea_water_potential_temperature`, `sea_water_salinity`, `sea_surface_height_above_geoid` |
| Observations | SST, Identity operator, 202 obs (`test/data/obs/sst.cdl`); ADT, 101 obs (`test/data/obs/adt.cdl`) |
| B | SABER `diffusion` central block only. T and S: univariate, horizontal + vertical. SSH: horizontal only |
| Diffusion parameters | Horizontal 200 km Gaussian, vertical 20 levels Gaussian, `method: implicit`, 4 iterations, 1000 normalization iterations |
| Minimizer | DRPCG, 1 outer loop, 10 inner iterations |
| MPI | 6 tasks, MOM6 `LAYOUT` 2 x 3, for both applications |
| Software | One bundle build: oops `f0c617bd`, saber `85f1fec4`, soca `develop` `5f3d068a`, ijedi `feature/mom6-uv` |
| Platform | JEDI GNU/OpenMPI singularity container, 16-core / 32 GB workstation |

The diffusion operator is the same SABER code on both sides. It is calibrated
separately on each geometry, because the two interfaces hand SABER different point
sets.

## Results

| | SOCA | i-jedi/MOM6 | SOCA / i-jedi |
|---|---|---|---|
| Wall time, 3DVar | **449 s** | **342 s** (292 s on a rerun) | 1.3 - 1.5 |
| `SABER::multiply` (diffusion TL + AD) | 433 s (97 %) | 330 s (97 %) | 1.31 |
| `Geometry::Geometry` | 0.6 s (+ 4.4 s offline `soca_gridgen.x`) | 2.3 s | |
| `State::State` (read background) | 0.73 s | 0.80 s | ~1 |
| Increment linear algebra (axpy, dot, `=`, `+=`) | 6.6 s | 4.4 s | 1.5 |
| Horizontal points handed to SABER | 115 200 | 77 852 | 1.48 |
| Memory per task after linearization | 1.17 GB | 0.99 GB | 1.18 |
| Memory per task at end (min - max) | 2.09 - 2.47 GB | 1.63 - 2.12 GB | ~1.25 |
| Total memory at end (`OOPS_STATS Run end`) | **14.1 GB** | **11.1 GB** | 1.27 |
| Norm reduction after 10 inner iterations | 0.041 | 0.049 | |
| Final J | 21.14 | 21.21 | |
| Diffusion calibration wall time (one-off) | 260 s | 145 s | 1.8 |

Notes on the numbers:

- **Timing noise.** On the same workstation, i-jedi ran in 342 s and then 292 s on a
  second run, with identical convergence. Treat differences of about 15 % as noise
  on this machine. SOCA was timed once.
- **Same problem on both sides.** The background departures are identical for both
  applications: SST RMS 1.2634 on 201 obs, ADT RMS 0.1240 on 100 obs. The final J
  differs slightly because each B was normalized on its own grid.

### Interpretation

- **B dominates.** At 1 degree, about 97 % of the time is in SABER diffusion. The
  benchmark therefore mostly measures the size of the field set each interface
  exposes to SABER. SOCA masks land but still carries it; i-jedi drops it from the
  function space. That is the source of the 1.3 - 1.5x gap in both B and the
  increment operations, and it is consistent with the calibration files: i-jedi's
  are a third smaller (95 MB vs 140 MB vertical).
- **I/O and geometry are negligible at 1 degree.** They are expected to matter more
  at 1/4 degree and above, and with more tasks.
- **Memory grows during the minimization.** Both applications grow by about
  50 - 100 MB per task per inner iteration. That is DRPCG storing its vectors, and it
  scales with the increment size.

## Setting up the test

### 1. Build

Both applications must come from the same bundle build, so they use identical oops
and saber code. Steps:

1. Clone soca into `bundle/soca` (`https://github.com/JCSDA-internal/soca`, branch
   `develop`) and run `git submodule update --init --recursive` to get MOM6.
2. Enable soca in `bundle/CMakeLists.txt`, after `ufo`:

   ```cmake
   option(BUILD_SOCA "Build soca" ON)
   if (BUILD_SOCA)
     ecbuild_bundle( PROJECT soca    SOURCE "./soca" )
   endif()
   ```

3. Reconfigure and build from the top-level build directory, not `build/ijedi`,
   because the bundle has to be reconfigured to pick up soca:
   `cd build && cmake . && make -j<N>`.

Executables used: `soca_gridgen.x`, `soca_error_covariance_toolbox.x`,
`soca_var.x` and `ijedi.x`.

### 2. Data

| File | Used by | Notes |
|---|---|---|
| `ocean_hgrid.nc` | both | MOM6 supergrid |
| `ocean_topog.nc` | both | variable `depth` |
| `MOM.res.nc` | both | `Temp`, `Salt`, `h`, `sfc`; `Layer` dimension = NZ |
| `sst.nc`, `adt.nc` | both | `ncgen -k nc4 -o sst.nc test/data/obs/sst.cdl`, and the same for `adt.cdl` |
| `obsop_name_map.yml` | both | `test/data/obs/obsop_name_map.yml` |
| `rossrad.nc` | SOCA gridgen | `<build>/soca/test/data_static/rossrad.nc`. It is a global climatology, so it is resolution independent |
| `fields_metadata.yml` | SOCA | copy of `<build>/soca/test/data_static/fields_metadata.yml` with the SSH `io name: ave_ssh` changed to `io name: sfc` |
| `diag_table`, `field_table` | SOCA | any minimal FMS tables, e.g. from a soca test work directory |

### 3. Directory layout

```
bench/
  ijedi/  var.yaml diffusion-init.yaml ocean_hgrid.nc ocean_topog.nc MOM.res.nc
          sst.nc adt.nc obsop_name_map.yml
  soca/   var.yaml diffusion-init.yaml gridgen.yml MOM_input mom6.nml
          diag_table field_table fields_metadata.yml rossrad.nc
          sst.nc adt.nc obsop_name_map.yml
          INPUT/ocean_hgrid.nc INPUT/ocean_topog.nc
          restarts/MOM.res.nc  RESTART/  out/
```

### 4. Run order

All steps use `LX*LY` MPI tasks, and both applications use the same task count.

```bash
# SOCA (from bench/soca)
srun -n $N soca_gridgen.x                  gridgen.yml          # writes soca_gridspec.nc
srun -n $N soca_error_covariance_toolbox.x diffusion-init.yaml  # writes diffusion_{hz,vt}.nc
srun -n $N soca_var.x                      var.yaml             # benchmark

# i-jedi (from bench/ijedi)
srun -n $N ijedi.x errortoolbox diffusion-init.yaml             # writes diffusion_{hz,vt}.nc
srun -n $N ijedi.x var          var.yaml                        # benchmark
```

Set `OOPS_TRAPFPE=0` and `OMP_NUM_THREADS=1` for every step.

The calibration files from SOCA and i-jedi are not interchangeable, because the two
point sets differ. Only the two `var` runs are the benchmark; gridgen and calibration
are one-off costs.

## Configuration files (1 degree)

The resolution-dependent entries are marked `# <- res`:

- `NIGLOBAL`, `NJGLOBAL`, `NZ`/`NK` and `LAYOUT`;
- the vertical diffusion `levels`.

Everything else, including the diffusion length scales, stays the same at every
resolution.

### i-jedi `var.yaml`

```yaml
cost function:
  cost type: 3D-Var
  jb evaluation: false
  time window:
    begin: 2020-12-13T00:00:00Z
    length: P6D
  geometry: &geom
    geometry_type: mom6
    ocean_hgrid: ocean_hgrid.nc
    ocean_topog: ocean_topog.nc
    vertical geometry from: MOM.res.nc
    coarsen factor: 1
    fringe width: 1
    has northern fold: true
    MOM_input:
      NIGLOBAL: 360      # <- res
      NJGLOBAL: 320      # <- res
      NZ: 75             # <- res
      LAYOUT: [2, 3]     # <- res
      MINIMUM_DEPTH: 0.5
  analysis variables:
  - sea_water_potential_temperature
  - sea_water_salinity
  - sea_surface_height_above_geoid
  background:
    date: 2020-12-16T00:00:00Z
    state variables:
    - sea_water_salinity
    - sea_water_potential_temperature
    - sea_surface_height_above_geoid
    - sea_water_cell_thickness
    io:
      filetype: mom6
      ocn_file: MOM.res.nc
      field io names: &ionames
        sea_water_potential_temperature: Temp
        sea_water_salinity: Salt
        sea_surface_height_above_geoid: sfc
        sea_water_cell_thickness: h
  background error:
    covariance model: SABER
    saber central block:
      saber block name: diffusion
      read:
        groups:
        - variables:
          - sea_water_potential_temperature
          - sea_water_salinity
          multivariate strategy: univariate
          horizontal:
            filepath: diffusion_hz
          vertical:
            levels: 75   # <- res
            filepath: diffusion_vt
            method: implicit
            iterations: 4
        - variables: [sea_surface_height_above_geoid]
          horizontal:
            filepath: diffusion_hz
  observations:
    observers:
    - obs space:
        name: SST
        obsdatain:
          engine: {type: H5File, obsfile: sst.nc}
        simulated variables: [seaSurfaceTemperature]
      obs operator:
        name: Identity
        observation alias file: obsop_name_map.yml
    - obs space:
        name: ADT
        obsdatain:
          engine: {type: H5File, obsfile: adt.nc}
        simulated variables: [absoluteDynamicTopography]
      obs operator:
        name: ADT
variational:
  minimizer:
    algorithm: DRPCG
  iterations:
  - geometry: *geom
    ninner: 10
    gradient norm reduction: 1e-10
output:
  io:
    filetype: mom6
    ocn_file: an.nc
    field io names: *ionames
final:
  increment:
    output:
      state component:
        io:
          filetype: mom6
          ocn_file: incr.nc
          field io names: *ionames
    geometry: *geom
final j evaluation: false
```

### i-jedi `diffusion-init.yaml`

```yaml
geometry:
  # same geometry block as in var.yaml
background:
  date: 2020-12-16T00:00:00Z
  state variables: [sea_water_potential_temperature]
  io:
    filetype: mom6
    ocn_file: MOM.res.nc
    field io names:
      sea_water_potential_temperature: Temp
background error:
  covariance model: SABER
  saber central block:
    saber block name: diffusion
    calibration:
      normalization:
        iterations: 1000
      groups:
      - horizontal:
          fixed value: 200000.0
          as gaussian: true
        write:
          filepath: diffusion_hz
      - vertical:
          fixed value: 20.0
          levels: 75   # <- res
          as gaussian: true
          method: implicit
          iterations: 4
        write:
          filepath: diffusion_vt
```

### SOCA `var.yaml`

This is the same as the i-jedi `var.yaml`, with the following blocks changed. The
`analysis variables`, `background error`, `observations` and `variational` blocks are
identical.

```yaml
cost function:
  geometry: &geom
    geom_grid_file: soca_gridspec.nc
    mom6_input_nml: mom6.nml
    fields metadata: fields_metadata.yml
  background:
    read_from_file: 1
    basename: restarts/
    ocn_filename: MOM.res.nc
    date: 2020-12-16T00:00:00Z
    state variables:
    - sea_water_salinity
    - sea_water_potential_temperature
    - sea_surface_height_above_geoid
    - sea_water_cell_thickness
output:
  datadir: out
  exp: bench
  type: an
final:
  increment:
    output:
      state component:
        datadir: out
        exp: bench
        type: incr
```

### SOCA `diffusion-init.yaml`

This uses the SOCA geometry and background blocks, with
`state variables: [sea_water_potential_temperature, sea_water_cell_thickness]`.
The calibration block is the same as the i-jedi one, with one addition:
`mask: "interp_mask"` in both the horizontal and the vertical group. SOCA's function
space contains land, so its calibration needs the mask; i-jedi's does not.

### SOCA `gridgen.yml`

```yaml
geometry:
  geom_grid_file: soca_gridspec.nc
  mom6_input_nml: mom6.nml
  fields metadata: fields_metadata.yml
  rossby file: rossrad.nc
```

### SOCA `mom6.nml`

The namelist must not be named `input.nml`. SOCA copies it to `input.nml` itself and
refuses a source file that already has that name.

```
 &MOM_input_nml
        output_directory = './',
        input_filename = 'n'
        restart_input_dir = 'restarts/',
        restart_output_dir = 'RESTART/',
        parameter_filename = 'MOM_input'
/
 &diag_manager_nml
 /
 &fms_io_nml
      max_files_w=100
      checksum_required=.false.
/
 &fms_nml
       clock_grain='MODULE'
       domains_stack_size = 20000000
       clock_flags='SYNC' /
```

### SOCA `MOM_input`

`MOM_input` is only used by `soca_gridgen.x`, which runs a full `initialize_MOM`. It
needs three things:

- the grid and topography;
- an analytic initial state, because gridgen keeps only the horizontal grid;
- a set of physics parameters that MOM6 requires at init (copied from SOCA's
  72x35x25 test).

None of this affects the 3DVar: SOCA takes the vertical grid from `h` in the restart.

```
USE_REGRIDDING = True
REGRIDDING_COORDINATE_MODE = "Z*"
ALE_COORDINATE_CONFIG = "UNIFORM"
DT = 3600
TRIPOLAR_N = True
NIGLOBAL = 360          ! <- res
NJGLOBAL = 320          ! <- res
NK = 75                 ! <- res
LAYOUT = 2, 3           ! <- res
INPUTDIR = "INPUT"
GRID_CONFIG = "mosaic"
GRID_FILE = "ocean_hgrid.nc"
TOPO_CONFIG = "file"
TOPO_FILE = "ocean_topog.nc"
TOPO_VARNAME = "depth"
MINIMUM_DEPTH = 0.5
MASKING_DEPTH = 0.0
MAXIMUM_DEPTH = 6500.0
USE_TRIPOLAR_GEOLONB_BUG = False
GRID_ROTATION_ANGLE_BUGS = False
RESTART_CHECKSUMS_REQUIRED = False
REPORT_UNUSED_PARAMS = False
! analytic initial state (gridgen only needs the horizontal grid)
THICKNESS_CONFIG = "uniform"
TS_CONFIG = "linear"
T_TOP = 20.0
T_RANGE = 18.0
S_TOP = 34.0
S_RANGE = 1.0
! physics parameters required at init (from soca 72x35x25 test)
USE_MEKE = True
USE_VARIABLE_MIXING = True
CHANNEL_DRAG = True
HBBL = 10.0
KV = 1.0E-04
CORIOLIS_SCHEME = "SADOURNY75_ENSTRO"
BOUND_CORIOLIS = True
MASS_WEIGHT_IN_PRESSURE_GRADIENT = True
LAPLACIAN = True
KH_VEL_SCALE = 0.01
KH_SIN_LAT = 2000.0
SMAGORINSKY_KH = True
SMAG_LAP_CONST = 0.15
AH_VEL_SCALE = 0.01
SMAGORINSKY_AH = True
SMAG_BI_CONST = 0.06
USE_LAND_MASK_FOR_HVISC = True
HMIX_FIXED = 0.5
MAXVEL = 5.0
BOUND_BT_CORRECTION = True
BT_PROJECT_VELOCITY = True
DYNAMIC_SURFACE_PRESSURE = True
BEBT = 0.2
DTBT = -0.9
MIXEDLAYER_RESTRAT = True
FOX_KEMPER_ML_RESTRAT_COEF = 1.0
MLE_FRONT_LENGTH = 200.0
MLE_USE_PBL_MLD = True
MLE_MLD_DECAY_TIME = 2.592E+06
ENERGETICS_SFC_PBL = True
KD = 1.5E-05
KD_MIN = 2.0E-06
HENYEY_IGW_BACKGROUND = True
USE_JACKSON_PARAM = True
MAX_RINO_IT = 205
EPBL_USTAR_MIN = 1.45842E-18
TRACER_ADVECTION_SCHEME = "PPM:H3"
KHTR = 50.0
CHECK_DIFFUSIVE_CFL = True
MAX_TR_DIFFUSION_CFL = 2.0
USE_NEUTRAL_DIFFUSION = True
MAXTRUNC = 1000
THICKNESSDIFFUSE = True
FRAZIL = True
BOUND_SALINITY = True
DT_THERM = 3600.0
```

## Extracting the numbers

All the numbers in the results table come from the `OOPS_STATS` lines in each
`var` log:

```bash
# runtime and total / per-task memory at start and end
grep "OOPS_STATS Run \(start\|end\)" var.log

# memory growth through the minimization
grep "Local Memory" var.log | grep -E "linearize 0|DRPCG (start|iteration|end)"

# per-timer min/max/avg across tasks (first parallel-timing block)
L=$(grep -n "Parallel Timing" var.log | head -1 | cut -d: -f1)
sed -n "$((L+3)),$((L+80))p" var.log | \
  grep -E "SABER::multiply|Diffusion::multiplySqrt|Geometry::Geometry |State::State |Increment::|Run::execute"

# check both runs solve the same problem
grep "Jo Dep" var.log | head -2
grep "Norm reduction" var.log | tail -1
grep "Nonlinear J =" var.log | tail -1
```

On Slurm, also check `sacct -j <jobid> --format=JobID,MaxRSS,AveRSS,Elapsed` to get
the per-step peak RSS. Run each application in its own job, on exclusive nodes.

## Going to higher resolution

- **Resolution settings.** For 1/4 degree, set `NIGLOBAL = 1440`, `NJGLOBAL = 1080`
  and `NK = 75` (with matching `NZ` and diffusion `levels`), and use a 1/4 degree
  grid, topography and restart.
- **Task count.** 1/4 degree has 13.5x the horizontal points of 1 degree. At the
  1 degree cost per point, the B alone would take about 75 min on 6 tasks, so start
  at `LAYOUT` 8 x 12 (96 tasks).
- **Scaling sweep.** Compare at equal task counts, and run a short sweep (for example
  4x6, 8x12 and 16x24) for both applications. How the I/O and geometry share of the
  cost changes with resolution and task count is the interesting part. At 1 degree
  that share was negligible.
- **Calibration cost.** Calibration costs about 1000 applications of the diffusion
  operator. If it becomes too expensive, lower `normalization: iterations`, but use
  the same value for both applications.
- **Isolating interface cost.** To measure interface overheads (I/O, halo exchanges,
  interpolation) without the diffusion dominating, rerun both `var.yaml` with
  `saber block name: ID` as the central block.

help([[
Load environment for running the i-jedi with Intel compilers and MPI on Ursa.
]])

local pkgName    = myModuleName()
local pkgVersion = myModuleVersion()
local pkgNameVer = myModuleFullName()

prepend_path("MODULEPATH", '/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core')

load("stack-oneapi/2024.2.1")
load("stack-intel-oneapi-mpi/2021.13")
load("intel-oneapi-mkl/2024.2.1")

load("git/2.43.5")
load("git-lfs/3.5.1")
load("hdf5/1.14.3")
load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("netcdf-cxx4/4.3.1")
load("udunits/2.2.28")
load("cmake/3.30.2")
load("ecbuild/3.7.2")
load("eigen/3.4.0")
load("openblas/0.3.27")
load("boost/1.84.0")
load("eckit/1.28.3")
load("fckit/0.13.2")
load("atlas/0.40.0")
load("fms/2024.02")
load("parallel-netcdf/1.12.3")
load("parallelio/2.6.2")
load("gsl-lite/0.37.0")
load("nccmp/1.9.1.0")
load("py-pycodestyle/2.11.0")
load("py-netcdf4/1.7.1.post2")

setenv("CC","mpiicx")
setenv("CXX","mpiicpx")
setenv("FC","mpiifort")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_F90", "ifort")

local mpiexec = '/apps/slurm/default/bin/srun'
local mpinproc = '-n'
setenv('MPIEXEC_EXEC', mpiexec)
setenv('MPIEXEC_NPROC', mpinproc)

whatis("Name: ".. pkgName)
whatis("Version: ".. tostring(pkgVersion))
whatis("Category: i-jedi")
whatis("Description: Load all libraries needed for i-jedi")

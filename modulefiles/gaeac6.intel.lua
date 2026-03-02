help([[
Load environment for running the i-jedi with Intel compilers and MPI on GaeaC6.
]])

local pkgName    = myModuleName()
local pkgVersion = myModuleVersion() or "0.0.1"
local pkgNameVer = myModuleFullName()

prepend_path("MODULEPATH", '/ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.3/envs/ue-oneapi-2024.2.1/install/modulefiles/Core')
prepend_path("MODULEPATH", '/ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.2/envs/ue-intel-2023.2.0/install/modulefiles/gcc/12.3.0')

-- below two lines get us access to the spack-stack modules
load("stack-oneapi/2024.2.1")
load("stack-cray-mpich/8.1.32")
load("git/2.42.0")
load("git-lfs/2.11.0")
load("hdf5/1.14.3")
load("netcdf-c/4.9.2")
load("netcdf-fortran/4.6.1")
load("netcdf-cxx4/4.3.1")
load("udunits/2.2.28")
load("cmake/3.27.9")
load("ecbuild/3.7.2")
load("eigen/3.4.0")
load("openblas/0.3.26")
load("boost/1.84.0")
load("eckit/1.28.3")
load("fckit/0.13.2")
load("atlas/0.40.0")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

local mpiexec = '/usr/bin/srun'
local mpinproc = '-n'
setenv('MPIEXEC_EXEC', mpiexec)
setenv('MPIEXEC_NPROC', mpinproc)

whatis("Name: ".. pkgName)
whatis("Version: ".. pkgVersion)
whatis("Category: i-jedi")
whatis("Description: Load all libraries needed for i-jedi")

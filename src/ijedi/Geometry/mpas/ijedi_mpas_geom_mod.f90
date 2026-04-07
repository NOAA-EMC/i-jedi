! (C) Copyright 2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module ijedi_mpas_geom_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use ESMF

use kinds, only: kind_real

use mpas_derived_types
use mpas_kind_types, only: RKIND
use mpas_constants, only: pii
use mpas_dmpar, only: mpas_dmpar_sum_int, mpas_dmpar_exch_halo_field
use mpas_subdriver, only: mpas_init, mpas_finalize
use atm_core
use mpas_pool_routines, only: mpas_pool_get_subpool, mpas_pool_get_dimension, &
                              mpas_pool_get_array, mpas_pool_get_field

implicit none
private

public :: ijedi_mpas_geom, &
          geom_setup, geom_clone, geom_delete

real(kind=kind_real), parameter :: RAD2DEG = 180.0_kind_real / real(pii, kind_real)
real(kind=kind_real), parameter :: HALF_PI = real(pii, kind_real) / 2.0_kind_real

logical, save :: esmf_initialized = .false.

character(len=1024) :: message

type :: ijedi_mpas_geom
  integer :: nCellsGlobal
  integer :: nCells
  integer :: nCellsSolve
  integer :: nVerticesGlobal
  integer :: nVertices
  integer :: nVerticesSolve
  integer :: nVertLevels
  integer :: nVertLevelsP1
  integer :: vertexDegree

  real(kind=RKIND), allocatable :: latCell(:), lonCell(:)
  real(kind=RKIND), allocatable :: areaCell(:)
  real(kind=RKIND), allocatable :: zgrid(:,:)
  integer, allocatable :: cellsOnVertex(:,:)
  integer, allocatable :: bdyMaskVertex(:)
  logical :: is_regional = .false.
  logical :: owns_mpas = .false.

  type(domain_type), pointer :: domain => null()
  type(core_type), pointer :: corelist => null()
  type(fckit_mpi_comm) :: comm

  contains
  procedure, public :: get_num_nodes_and_elements
  procedure, public :: get_coords_and_connectivities
end type ijedi_mpas_geom

contains

! ------------------------------------------------------------------------------

subroutine geom_setup(self, f_conf, comm)

  type(ijedi_mpas_geom), intent(inout) :: self
  type(fckit_configuration), intent(in) :: f_conf
  type(fckit_mpi_comm), intent(in) :: comm

  character(len=512) :: nml_file, streams_file
  character(len=:), allocatable :: str
  type(mpas_pool_type), pointer :: meshPool
  type(block_type), pointer :: block_ptr
  real(kind=RKIND), pointer :: r1d_ptr(:), r2d_ptr(:,:)
  integer, pointer :: i0d_ptr, i1d_ptr(:), i2d_ptr(:,:)
  integer :: ii

  self%comm = comm

  call f_conf%get_or_die("nml_file", str)
  nml_file = str
  call f_conf%get_or_die("streams_file", str)
  streams_file = str

  if (.not. esmf_initialized) then
    call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN)
    esmf_initialized = .true.
  end if

  call mpas_init(self%corelist, self%domain, &
                 external_comm=self%comm%communicator(), &
                 namelistFileParam=trim(nml_file), &
                 streamsFileParam=trim(streams_file))
  self%owns_mpas = .true.

  block_ptr => self%domain%blocklist
  call mpas_pool_get_subpool(block_ptr%structs, 'mesh', meshPool)

  call mpas_pool_get_dimension(block_ptr%dimensions, 'nCells', i0d_ptr)
  self%nCells = i0d_ptr
  call mpas_pool_get_dimension(block_ptr%dimensions, 'nCellsSolve', i0d_ptr)
  self%nCellsSolve = i0d_ptr
  call mpas_dmpar_sum_int(self%domain%dminfo, self%nCellsSolve, self%nCellsGlobal)

  call mpas_pool_get_dimension(block_ptr%dimensions, 'nVertices', i0d_ptr)
  self%nVertices = i0d_ptr
  call mpas_pool_get_dimension(block_ptr%dimensions, 'nVerticesSolve', i0d_ptr)
  self%nVerticesSolve = i0d_ptr
  call mpas_dmpar_sum_int(self%domain%dminfo, self%nVerticesSolve, self%nVerticesGlobal)

  call mpas_pool_get_dimension(block_ptr%dimensions, 'nVertLevels', i0d_ptr)
  self%nVertLevels = i0d_ptr
  call mpas_pool_get_dimension(block_ptr%dimensions, 'nVertLevelsP1', i0d_ptr)
  self%nVertLevelsP1 = i0d_ptr
  call mpas_pool_get_dimension(block_ptr%dimensions, 'vertexDegree', i0d_ptr)
  self%vertexDegree = i0d_ptr

  allocate(self%latCell(self%nCells))
  allocate(self%lonCell(self%nCells))
  allocate(self%areaCell(self%nCells))
  allocate(self%cellsOnVertex(self%vertexDegree, self%nVertices))
  allocate(self%bdyMaskVertex(self%nVertices))

  call mpas_pool_get_array(meshPool, 'latCell', r1d_ptr)
  self%latCell = r1d_ptr(1:self%nCells)
  where (self%latCell > HALF_PI) self%latCell = HALF_PI
  where (self%latCell < -HALF_PI) self%latCell = -HALF_PI

  call mpas_pool_get_array(meshPool, 'lonCell', r1d_ptr)
  self%lonCell = r1d_ptr(1:self%nCells)
  call mpas_pool_get_array(meshPool, 'areaCell', r1d_ptr)
  self%areaCell = r1d_ptr(1:self%nCells)

  call mpas_pool_get_array(meshPool, 'cellsOnVertex', i2d_ptr)
  self%cellsOnVertex = i2d_ptr(1:self%vertexDegree, 1:self%nVertices)

  call mpas_pool_get_array(meshPool, 'bdyMaskVertex', i1d_ptr)
  self%bdyMaskVertex = i1d_ptr(1:self%nVertices)
  ! bdyMaskVertex == 7 is the MPAS convention for boundary-zone vertices in
  ! regional meshes (value 7 = outermost relaxation-zone boundary).
  self%is_regional = any(self%bdyMaskVertex(1:self%nVerticesSolve) == 7)

  allocate(self%zgrid(self%nVertLevelsP1, self%nCells))
  call mpas_pool_get_array(meshPool, 'zgrid', r2d_ptr)
  self%zgrid = r2d_ptr(1:self%nVertLevelsP1, 1:self%nCells)

end subroutine geom_setup

! ------------------------------------------------------------------------------

subroutine geom_clone(self, other)

  type(ijedi_mpas_geom), intent(inout) :: self
  type(ijedi_mpas_geom), intent(in) :: other

  integer :: ii

  self%comm = other%comm

  self%nCellsGlobal = other%nCellsGlobal
  self%nCells = other%nCells
  self%nCellsSolve = other%nCellsSolve
  self%nVerticesGlobal = other%nVerticesGlobal
  self%nVertices = other%nVertices
  self%nVerticesSolve = other%nVerticesSolve
  self%nVertLevels = other%nVertLevels
  self%nVertLevelsP1 = other%nVertLevelsP1
  self%vertexDegree = other%vertexDegree
  self%is_regional = other%is_regional
  self%owns_mpas = .false.

  if (.not. allocated(self%latCell)) allocate(self%latCell(other%nCells))
  if (.not. allocated(self%lonCell)) allocate(self%lonCell(other%nCells))
  if (.not. allocated(self%areaCell)) allocate(self%areaCell(other%nCells))
  if (.not. allocated(self%zgrid)) allocate(self%zgrid(other%nVertLevelsP1, other%nCells))
  if (.not. allocated(self%cellsOnVertex)) &
    allocate(self%cellsOnVertex(other%vertexDegree, other%nVertices))
  if (.not. allocated(self%bdyMaskVertex)) allocate(self%bdyMaskVertex(other%nVertices))

  self%latCell = other%latCell
  self%lonCell = other%lonCell
  self%areaCell = other%areaCell
  self%zgrid = other%zgrid
  self%cellsOnVertex = other%cellsOnVertex
  self%bdyMaskVertex = other%bdyMaskVertex

  self%corelist => other%corelist
  self%domain => other%domain

end subroutine geom_clone

! ------------------------------------------------------------------------------

subroutine geom_delete(self)

  type(ijedi_mpas_geom), intent(inout) :: self

  if (allocated(self%latCell)) deallocate(self%latCell)
  if (allocated(self%lonCell)) deallocate(self%lonCell)
  if (allocated(self%areaCell)) deallocate(self%areaCell)
  if (allocated(self%zgrid)) deallocate(self%zgrid)
  if (allocated(self%cellsOnVertex)) deallocate(self%cellsOnVertex)
  if (allocated(self%bdyMaskVertex)) deallocate(self%bdyMaskVertex)

  if (self%owns_mpas .and. associated(self%corelist) .and. associated(self%domain)) then
    call mpas_finalize(self%corelist, self%domain)
  end if

  self%owns_mpas = .false.
  nullify(self%corelist)
  nullify(self%domain)

end subroutine geom_delete

! ------------------------------------------------------------------------------

subroutine get_num_nodes_and_elements(self, num_nodes, num_tris)

  class(ijedi_mpas_geom), intent(in) :: self
  integer, intent(out) :: num_nodes
  integer, intent(out) :: num_tris

  integer :: nVerticesBdy7

  num_nodes = self%nCells
  num_tris = self%nVerticesSolve

  if (self%is_regional) then
    nVerticesBdy7 = count(self%bdyMaskVertex(1:self%nVerticesSolve) == 7)
    num_tris = self%nVerticesSolve - nVerticesBdy7
  end if

end subroutine get_num_nodes_and_elements

! ------------------------------------------------------------------------------

subroutine get_coords_and_connectivities(self, num_nodes, num_tri_boundary_nodes, &
                                         lons, lats, ghosts, global_indices, &
                                         remote_indices, partition, &
                                         raw_tri_boundary_nodes)

  class(ijedi_mpas_geom), intent(in) :: self
  integer, intent(in) :: num_nodes
  integer, intent(in) :: num_tri_boundary_nodes
  real(kind_real), intent(out) :: lons(num_nodes)
  real(kind_real), intent(out) :: lats(num_nodes)
  integer, intent(out) :: ghosts(num_nodes)
  integer, intent(out) :: global_indices(num_nodes)
  integer, intent(out) :: remote_indices(num_nodes)
  integer, intent(out) :: partition(num_nodes)
  integer, intent(out) :: raw_tri_boundary_nodes(num_tri_boundary_nodes)

  integer :: i, iVertValid
  type(field1DInteger), pointer :: indexToCellID, iTmp

  lons = self%lonCell * RAD2DEG
  lats = self%latCell * RAD2DEG

  ghosts = 1
  ghosts(1:self%nCellsSolve) = 0

  call mpas_pool_get_field(self%domain%blocklist%allFields, 'indexToCellID', indexToCellID)
  call mpas_duplicate_field(indexToCellID, iTmp)

  global_indices(1:num_nodes) = indexToCellID%array(1:num_nodes)

  iTmp%array(:) = -1
  do i = 1, self%nCellsSolve
    iTmp%array(i) = i
  end do
  call mpas_dmpar_exch_halo_field(iTmp)
  remote_indices(1:num_nodes) = iTmp%array(1:num_nodes)

  iTmp%array(:) = -1
  do i = 1, self%nCellsSolve
    iTmp%array(i) = self%comm%rank()
  end do
  call mpas_dmpar_exch_halo_field(iTmp)
  partition(1:num_nodes) = iTmp%array(1:num_nodes)

  call mpas_deallocate_field(iTmp)

  iVertValid = 1
  do i = 1, self%nVerticesSolve
    if (self%bdyMaskVertex(i) == 7) cycle
    raw_tri_boundary_nodes(3*(iVertValid-1)+1) = indexToCellID%array(self%cellsOnVertex(1,i))
    raw_tri_boundary_nodes(3*(iVertValid-1)+2) = indexToCellID%array(self%cellsOnVertex(2,i))
    raw_tri_boundary_nodes(3*(iVertValid-1)+3) = indexToCellID%array(self%cellsOnVertex(3,i))
    iVertValid = iVertValid + 1
  end do

end subroutine get_coords_and_connectivities

! ------------------------------------------------------------------------------

end module ijedi_mpas_geom_mod


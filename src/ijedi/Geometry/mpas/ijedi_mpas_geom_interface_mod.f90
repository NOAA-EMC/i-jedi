! (C) Copyright 2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_setup(c_geom, c_conf, c_comm) &
    bind(c, name='ijedi_mpas_geom_setup_f90')

  use iso_c_binding
  use fckit_configuration_module, only: fckit_configuration
  use fckit_mpi_module, only: fckit_mpi_comm
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), intent(out) :: c_geom
  type(c_ptr), value, intent(in) :: c_conf
  type(c_ptr), value, intent(in) :: c_comm

  type(ijedi_mpas_geom), pointer :: geom
  type(fckit_configuration) :: conf
  type(fckit_mpi_comm) :: comm

  conf = fckit_configuration(c_conf)
  comm = fckit_mpi_comm(c_comm)

  allocate(geom)
  call geom_setup(geom, conf, comm)
  c_geom = c_loc(geom)

  call conf%final()
  call comm%final()

end subroutine c_ijedi_mpas_geom_setup

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_clone(c_geom, c_other) &
    bind(c, name='ijedi_mpas_geom_clone_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), intent(out) :: c_geom
  type(c_ptr), value, intent(in) :: c_other

  type(ijedi_mpas_geom), pointer :: geom, other

  call c_f_pointer(c_other, other)
  allocate(geom)
  call geom_clone(geom, other)
  c_geom = c_loc(geom)

end subroutine c_ijedi_mpas_geom_clone

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_delete(c_geom) &
    bind(c, name='ijedi_mpas_geom_delete_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), intent(inout) :: c_geom

  type(ijedi_mpas_geom), pointer :: geom

  if (.not. c_associated(c_geom)) return

  call c_f_pointer(c_geom, geom)
  call geom_delete(geom)
  deallocate(geom)
  c_geom = c_null_ptr

end subroutine c_ijedi_mpas_geom_delete

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_get_num_nodes_and_elements(c_geom, num_nodes, num_tris) &
    bind(c, name='ijedi_mpas_geom_get_num_nodes_and_elements_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  integer(c_int), intent(out) :: num_nodes
  integer(c_int), intent(out) :: num_tris

  type(ijedi_mpas_geom), pointer :: geom

  call c_f_pointer(c_geom, geom)
  call geom%get_num_nodes_and_elements(num_nodes, num_tris)

end subroutine c_ijedi_mpas_geom_get_num_nodes_and_elements

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_get_coords_and_connectivities(c_geom, &
    num_nodes, lons, lats, ghosts, global_indices, remote_indices, partition, &
    num_tri_nodes, raw_tri_boundary_nodes) &
    bind(c, name='ijedi_mpas_geom_get_coords_and_connectivities_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  integer(c_int), intent(in) :: num_nodes
  real(c_double), intent(out) :: lons(num_nodes)
  real(c_double), intent(out) :: lats(num_nodes)
  integer(c_int), intent(out) :: ghosts(num_nodes)
  integer(c_int), intent(out) :: global_indices(num_nodes)
  integer(c_int), intent(out) :: remote_indices(num_nodes)
  integer(c_int), intent(out) :: partition(num_nodes)
  integer(c_int), intent(in) :: num_tri_nodes
  integer(c_int), intent(out) :: raw_tri_boundary_nodes(num_tri_nodes)

  type(ijedi_mpas_geom), pointer :: geom

  call c_f_pointer(c_geom, geom)
  call geom%get_coords_and_connectivities(num_nodes, num_tri_nodes, &
      lons, lats, ghosts, global_indices, remote_indices, partition, &
      raw_tri_boundary_nodes)

end subroutine c_ijedi_mpas_geom_get_coords_and_connectivities

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_get_vertical_resolution(c_geom, nVertLevels) &
    bind(c, name='ijedi_mpas_geom_get_vertical_resolution_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  integer(c_int), intent(out) :: nVertLevels

  type(ijedi_mpas_geom), pointer :: geom

  call c_f_pointer(c_geom, geom)
  nVertLevels = geom%nVertLevels

end subroutine c_ijedi_mpas_geom_get_vertical_resolution

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_get_local_cell_counts(c_geom, nCells, nCellsSolve) &
    bind(c, name='ijedi_mpas_geom_get_local_cell_counts_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  integer(c_int), intent(out) :: nCells
  integer(c_int), intent(out) :: nCellsSolve

  type(ijedi_mpas_geom), pointer :: geom

  call c_f_pointer(c_geom, geom)
  nCells      = geom%nCells
  nCellsSolve = geom%nCellsSolve

end subroutine c_ijedi_mpas_geom_get_local_cell_counts

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_get_global_cell_count(c_geom, nCellsGlobal) &
    bind(c, name='ijedi_mpas_geom_get_global_cell_count_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  integer(c_int), intent(out) :: nCellsGlobal

  type(ijedi_mpas_geom), pointer :: geom

  call c_f_pointer(c_geom, geom)
  nCellsGlobal = geom%nCellsGlobal

end subroutine c_ijedi_mpas_geom_get_global_cell_count

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_geom_get_area(c_geom, n, area) &
    bind(c, name='ijedi_mpas_geom_get_area_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  integer(c_int), intent(in) :: n
  real(c_double), intent(out) :: area(n)

  type(ijedi_mpas_geom), pointer :: geom

  call c_f_pointer(c_geom, geom)
  area(1:n) = real(geom%areaCell(1:n), c_double)

end subroutine c_ijedi_mpas_geom_get_area


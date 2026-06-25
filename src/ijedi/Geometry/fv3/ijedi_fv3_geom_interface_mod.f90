module fv3jedi_geom_interface_mod

use atlas_module, only: atlas_fieldset, atlas_functionspace
use iso_c_binding

use fckit_mpi_module,           only: fckit_mpi_comm
use fckit_configuration_module, only: fckit_configuration

use ijedi_kinds_mod
use ijedi_fv3_geom_mod

implicit none
private

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine c_fv3_geom_initialize(c_conf, c_comm) bind(c, name='f_fv3_geom_initialize')

! Args
type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_comm

! Locals
type(fckit_mpi_comm)        :: f_comm
type(fckit_configuration)   :: f_conf

! Fortran APIs
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm(c_comm)

! Call impl
call fv3_geom_initialize(f_conf, f_comm)

end subroutine c_fv3_geom_initialize

! --------------------------------------------------------------------------------------------------

subroutine c_fv3_geom_create(c_geom_conf, c_geom_vars, c_comm) bind(c, name='f_fv3_geom_create')

!Arguments
type(c_ptr), value, intent(in) :: c_geom_conf
type(c_ptr), value, intent(in) :: c_geom_vars
type(c_ptr), value, intent(in) :: c_comm

! Locals
type(fckit_configuration) :: f_geom_conf
type(fckit_configuration) :: f_geom_vars
type(fckit_mpi_comm)      :: f_comm

! Fortran APIs
! ------------
f_geom_conf = fckit_configuration(c_geom_conf)
f_geom_vars = fckit_configuration(c_geom_vars)
f_comm = fckit_mpi_comm(c_comm)

! Call implementation
! -------------------
call fv3_geom_create(f_geom_conf, f_geom_vars, f_comm)

end subroutine c_fv3_geom_create

! --------------------------------------------------------------------------------------------------

end module fv3jedi_geom_interface_mod

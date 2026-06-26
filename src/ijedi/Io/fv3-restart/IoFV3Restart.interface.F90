module ijedi_io_fv3_restart_interface_mod

use atlas_module, only: atlas_fieldset
use iso_c_binding

use fckit_configuration_module,  only: fckit_configuration
use ijedi_fv3_restart_mod,       only: ijedi_fv3_restart_read, ijedi_fv3_restart_write

implicit none
private

contains

subroutine c_ijedi_io_fv3_restart_read(c_conf, c_geom_conf, c_fieldset, c_fileionames, &
                                       c_fileioscaling) bind(c, name='ijedi_io_fv3_restart_read_f90')

type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_geom_conf
type(c_ptr), value, intent(in) :: c_fieldset
type(c_ptr), value, intent(in) :: c_fileionames
type(c_ptr), value, intent(in) :: c_fileioscaling

type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_geom_conf
type(fckit_configuration) :: f_fileionames
type(fckit_configuration) :: f_fileioscaling
type(atlas_fieldset)      :: f_fieldset

f_conf = fckit_configuration(c_conf)
f_geom_conf = fckit_configuration(c_geom_conf)
f_fileionames = fckit_configuration(c_fileionames)
f_fileioscaling = fckit_configuration(c_fileioscaling)
f_fieldset = atlas_fieldset(c_fieldset)

call ijedi_fv3_restart_read(f_conf, f_geom_conf, f_fieldset, f_fileionames, f_fileioscaling)

call f_fieldset%final()

end subroutine c_ijedi_io_fv3_restart_read

subroutine c_ijedi_io_fv3_restart_write(c_conf, c_geom_conf, c_fieldset, c_fileionames, &
                                        c_fileioscaling) bind(c, name='ijedi_io_fv3_restart_write_f90')

type(c_ptr), value, intent(in) :: c_conf
type(c_ptr), value, intent(in) :: c_geom_conf
type(c_ptr), value, intent(in) :: c_fieldset
type(c_ptr), value, intent(in) :: c_fileionames
type(c_ptr), value, intent(in) :: c_fileioscaling

type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_geom_conf
type(fckit_configuration) :: f_fileionames
type(fckit_configuration) :: f_fileioscaling
type(atlas_fieldset)      :: f_fieldset

f_conf = fckit_configuration(c_conf)
f_geom_conf = fckit_configuration(c_geom_conf)
f_fileionames = fckit_configuration(c_fileionames)
f_fileioscaling = fckit_configuration(c_fileioscaling)
f_fieldset = atlas_fieldset(c_fieldset)

call ijedi_fv3_restart_write(f_conf, f_geom_conf, f_fieldset, f_fileionames, f_fileioscaling)

call f_fieldset%final()

end subroutine c_ijedi_io_fv3_restart_write

end module ijedi_io_fv3_restart_interface_mod

! (C) Copyright 2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

subroutine c_ijedi_mpas_state_read_field(c_geom, c_filepath, c_stream_name, c_var_name, &
                                         nvals, global_indices, nlevels, values, scaling) &
    bind(c, name='ijedi_mpas_state_read_field_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod, only: ijedi_mpas_geom
  use mpas_derived_types, only: MPAS_streamManager_type, MPAS_Time_type, field1DReal, field2DReal, &
                                field1DInteger, field3DReal, mpas_pool_type
  use mpas_kind_types, only: StrKIND
  use mpas_pool_routines, only: mpas_pool_get_field, mpas_pool_get_subpool, mpas_pool_get_dimension
  use mpas_stream_manager
  use mpas_timekeeping, only: MPAS_NOW, MPAS_START_TIME, mpas_get_clock_time, mpas_get_time, &
                              mpas_set_clock_time

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  character(kind=c_char), intent(in) :: c_filepath(*)
  character(kind=c_char), intent(in) :: c_stream_name(*)
  character(kind=c_char), intent(in) :: c_var_name(*)
  integer(c_int), intent(in) :: nvals
  integer(c_int), intent(in) :: global_indices(nvals)
  integer(c_int), intent(in) :: nlevels
  real(c_double), intent(out) :: values(nlevels, nvals)
  real(c_double), intent(in) :: scaling

  type(ijedi_mpas_geom), pointer :: geom
  type(MPAS_streamManager_type), pointer :: manager
  type(MPAS_Time_type) :: local_time
  type(field1DInteger), pointer :: indexToCellID
  type(field1DReal), pointer :: f1
  type(field2DReal), pointer :: f2
  type(field3DReal), pointer :: scalars
  type(mpas_pool_type), pointer :: statePool, diagPool
  integer, pointer :: index_qv

  character(len=StrKIND) :: filepath, stream_name, var_name, dateTimeString
  integer :: ierr, i, k, nCells, gid, lid
  integer, allocatable :: gid_to_local(:)
  character(len=1024) :: message

  call c_f_pointer(c_geom, geom)
  call c_to_f_string(c_filepath, filepath)
  call c_to_f_string(c_stream_name, stream_name)
  call c_to_f_string(c_var_name, var_name)

  values(:, :) = 0.0_c_double

  manager => geom%domain%streamManager
  local_time = mpas_get_clock_time(geom%domain%clock, MPAS_NOW, ierr)
  if (ierr /= 0) call abor1_ftn('ijedi_mpas_state_read_field_f90: mpas_get_clock_time failed')

  call mpas_get_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
  if (ierr /= 0) call abor1_ftn('ijedi_mpas_state_read_field_f90: mpas_get_time failed')

  call mpas_set_clock_time(geom%domain%clock, local_time, MPAS_NOW)
  call mpas_set_clock_time(geom%domain%clock, local_time, MPAS_START_TIME)

  ierr = 0
  call MPAS_stream_mgr_read(manager, streamID=trim(stream_name), when=trim(dateTimeString), &
                            rightNow=.true., ierr=ierr)
  if (ierr /= 0) then
    write(message, '(a,i0)') 'ijedi_mpas_state_read_field_f90: MPAS_stream_mgr_read failed ierr=', ierr
    call abor1_ftn(message)
  end if

  call mpas_pool_get_field(geom%domain%blocklist%allFields, 'indexToCellID', indexToCellID)
  nCells = size(indexToCellID%array)
  allocate(gid_to_local(geom%nCellsGlobal))
  gid_to_local = 0
  do i = 1, nCells
    gid = indexToCellID%array(i)
    if (gid >= 1 .and. gid <= geom%nCellsGlobal) gid_to_local(gid) = i
  end do

  nullify(f2)
  call mpas_pool_get_field(geom%domain%blocklist%allFields, trim(var_name), f2)
  if (.not. associated(f2)) then
    nullify(statePool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'state', statePool)
    if (associated(statePool)) call mpas_pool_get_field(statePool, trim(var_name), f2)
  end if
  if (.not. associated(f2)) then
    nullify(diagPool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'diag', diagPool)
    if (associated(diagPool)) call mpas_pool_get_field(diagPool, trim(var_name), f2)
  end if
  if (associated(f2)) then
    if (size(f2%array, 1) /= nlevels) then
      call abor1_ftn('ijedi_mpas_state_read_field_f90: vertical levels mismatch for 2D field')
    end if
    do i = 1, nvals
      gid = global_indices(i)
      if (gid < 1 .or. gid > geom%nCellsGlobal) cycle
      lid = gid_to_local(gid)
      if (lid <= 0) cycle
      do k = 1, nlevels
        values(k, i) = real(f2%array(k, lid), kind=c_double)
      end do
    end do
    if (scaling /= 0.0_c_double) values(:, :) = values(:, :) * scaling
    deallocate(gid_to_local)
    return
  end if

  nullify(f1)
  call mpas_pool_get_field(geom%domain%blocklist%allFields, trim(var_name), f1)
  if (.not. associated(f1)) then
    nullify(statePool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'state', statePool)
    if (associated(statePool)) call mpas_pool_get_field(statePool, trim(var_name), f1)
  end if
  if (.not. associated(f1)) then
    nullify(diagPool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'diag', diagPool)
    if (associated(diagPool)) call mpas_pool_get_field(diagPool, trim(var_name), f1)
  end if
  if (.not. associated(f1)) then
    if (trim(var_name) == 'qv') then
      nullify(statePool)
      call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'state', statePool)
      if (associated(statePool)) then
        nullify(scalars)
        nullify(index_qv)
        call mpas_pool_get_field(statePool, 'scalars', scalars)
        call mpas_pool_get_dimension(statePool, 'index_qv', index_qv)
        if (associated(scalars) .and. associated(index_qv)) then
          if (nlevels <= 1) call abor1_ftn('ijedi_mpas_state_read_field_f90: qv expects 3D field')
          do i = 1, nvals
            gid = global_indices(i)
            if (gid < 1 .or. gid > geom%nCellsGlobal) cycle
            lid = gid_to_local(gid)
            if (lid <= 0) cycle
            do k = 1, nlevels
              values(k, i) = real(scalars%array(index_qv, k, lid), kind=c_double)
            end do
          end do
          if (scaling /= 0.0_c_double) values(:, :) = values(:, :) * scaling
          deallocate(gid_to_local)
          return
        end if
      end if
    end if
    call abor1_ftn('ijedi_mpas_state_read_field_f90: field not found: ' // trim(var_name))
  end if

  if (nlevels /= 1) then
    call abor1_ftn('ijedi_mpas_state_read_field_f90: expected 1 level for 1D field')
  end if

  do i = 1, nvals
    gid = global_indices(i)
    if (gid < 1 .or. gid > geom%nCellsGlobal) cycle
    lid = gid_to_local(gid)
    if (lid <= 0) cycle
    values(1, i) = real(f1%array(lid), kind=c_double)
  end do
  if (scaling /= 0.0_c_double) values(:, :) = values(:, :) * scaling

  deallocate(gid_to_local)

end subroutine c_ijedi_mpas_state_read_field

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_state_write_field(c_geom, c_filepath, c_stream_name, c_var_name, &
                                          nvals, global_indices, nlevels, values) &
    bind(c, name='ijedi_mpas_state_write_field_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod, only: ijedi_mpas_geom
  use mpas_derived_types, only: field1DReal, field2DReal, field1DInteger, field3DReal, mpas_pool_type
  use mpas_kind_types, only: RKIND
  use mpas_pool_routines, only: mpas_pool_get_field, mpas_pool_get_subpool, mpas_pool_get_dimension

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  character(kind=c_char), intent(in) :: c_filepath(*)
  character(kind=c_char), intent(in) :: c_stream_name(*)
  character(kind=c_char), intent(in) :: c_var_name(*)
  integer(c_int), intent(in) :: nvals
  integer(c_int), intent(in) :: global_indices(nvals)
  integer(c_int), intent(in) :: nlevels
  real(c_double), intent(in) :: values(nlevels, nvals)

  type(ijedi_mpas_geom), pointer :: geom
  type(field1DInteger), pointer :: indexToCellID
  type(field1DReal), pointer :: f1
  type(field2DReal), pointer :: f2
  type(field3DReal), pointer :: scalars
  type(mpas_pool_type), pointer :: statePool, diagPool
  integer, pointer :: index_qv

  character(len=1024) :: filepath, stream_name, var_name
  integer :: i, k, nCells, gid, lid
  integer, allocatable :: gid_to_local(:)

  call c_f_pointer(c_geom, geom)
  call c_to_f_string(c_filepath, filepath)
  call c_to_f_string(c_stream_name, stream_name)
  call c_to_f_string(c_var_name, var_name)

  call mpas_pool_get_field(geom%domain%blocklist%allFields, 'indexToCellID', indexToCellID)
  nCells = size(indexToCellID%array)
  allocate(gid_to_local(geom%nCellsGlobal))
  gid_to_local = 0
  do i = 1, nCells
    gid = indexToCellID%array(i)
    if (gid >= 1 .and. gid <= geom%nCellsGlobal) gid_to_local(gid) = i
  end do

  nullify(f2)
  call mpas_pool_get_field(geom%domain%blocklist%allFields, trim(var_name), f2)
  if (.not. associated(f2)) then
    nullify(statePool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'state', statePool)
    if (associated(statePool)) call mpas_pool_get_field(statePool, trim(var_name), f2)
  end if
  if (.not. associated(f2)) then
    nullify(diagPool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'diag', diagPool)
    if (associated(diagPool)) call mpas_pool_get_field(diagPool, trim(var_name), f2)
  end if
  if (associated(f2)) then
    if (size(f2%array, 1) /= nlevels) then
      call abor1_ftn('ijedi_mpas_state_write_field_f90: vertical levels mismatch for 2D field')
    end if
    do i = 1, nvals
      gid = global_indices(i)
      if (gid < 1 .or. gid > geom%nCellsGlobal) cycle
      lid = gid_to_local(gid)
      if (lid <= 0) cycle
      do k = 1, nlevels
        f2%array(k, lid) = real(values(k, i), kind=RKIND)
      end do
    end do
    deallocate(gid_to_local)
    return
  end if

  nullify(f1)
  call mpas_pool_get_field(geom%domain%blocklist%allFields, trim(var_name), f1)
  if (.not. associated(f1)) then
    nullify(statePool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'state', statePool)
    if (associated(statePool)) call mpas_pool_get_field(statePool, trim(var_name), f1)
  end if
  if (.not. associated(f1)) then
    nullify(diagPool)
    call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'diag', diagPool)
    if (associated(diagPool)) call mpas_pool_get_field(diagPool, trim(var_name), f1)
  end if
  if (.not. associated(f1)) then
    if (trim(var_name) == 'qv') then
      nullify(statePool)
      call mpas_pool_get_subpool(geom%domain%blocklist%structs, 'state', statePool)
      if (associated(statePool)) then
        nullify(scalars)
        nullify(index_qv)
        call mpas_pool_get_field(statePool, 'scalars', scalars)
        call mpas_pool_get_dimension(statePool, 'index_qv', index_qv)
        if (associated(scalars) .and. associated(index_qv)) then
          if (nlevels <= 1) call abor1_ftn('ijedi_mpas_state_write_field_f90: qv expects 3D field')
          do i = 1, nvals
            gid = global_indices(i)
            if (gid < 1 .or. gid > geom%nCellsGlobal) cycle
            lid = gid_to_local(gid)
            if (lid <= 0) cycle
            do k = 1, nlevels
              scalars%array(index_qv, k, lid) = real(values(k, i), kind=RKIND)
            end do
          end do
          deallocate(gid_to_local)
          return
        end if
      end if
    end if
    call abor1_ftn('ijedi_mpas_state_write_field_f90: field not found: ' // trim(var_name))
  end if

  if (nlevels /= 1) then
    call abor1_ftn('ijedi_mpas_state_write_field_f90: expected 1 level for 1D field')
  end if

  do i = 1, nvals
    gid = global_indices(i)
    if (gid < 1 .or. gid > geom%nCellsGlobal) cycle
    lid = gid_to_local(gid)
    if (lid <= 0) cycle
    f1%array(lid) = real(values(1, i), kind=RKIND)
  end do

  deallocate(gid_to_local)

end subroutine c_ijedi_mpas_state_write_field

! ------------------------------------------------------------------------------

subroutine c_ijedi_mpas_state_write_flush(c_geom, c_filepath, c_stream_name) &
    bind(c, name='ijedi_mpas_state_write_flush_f90')

  use iso_c_binding
  use ijedi_mpas_geom_mod, only: ijedi_mpas_geom
  use mpas_derived_types, only: MPAS_streamManager_type, MPAS_Time_type
  use mpas_kind_types, only: StrKIND
  use mpas_stream_manager
  use mpas_timekeeping, only: MPAS_NOW, mpas_get_clock_time, mpas_get_time

  implicit none

  type(c_ptr), value, intent(in) :: c_geom
  character(kind=c_char), intent(in) :: c_filepath(*)
  character(kind=c_char), intent(in) :: c_stream_name(*)

  type(ijedi_mpas_geom), pointer :: geom
  type(MPAS_streamManager_type), pointer :: manager
  type(MPAS_Time_type) :: local_time

  character(len=StrKIND) :: filepath, stream_name, dateTimeString
  integer :: ierr
  character(len=1024) :: message
  integer, parameter :: STREAM_PROPERTY_FILENAME = 7

  call c_f_pointer(c_geom, geom)
  call c_to_f_string(c_filepath, filepath)
  call c_to_f_string(c_stream_name, stream_name)

  manager => geom%domain%streamManager
  local_time = mpas_get_clock_time(geom%domain%clock, MPAS_NOW, ierr)
  if (ierr /= 0) call abor1_ftn('ijedi_mpas_state_write_flush_f90: mpas_get_clock_time failed')

  call mpas_get_time(local_time, dateTimeString=dateTimeString, ierr=ierr)
  if (ierr /= 0) call abor1_ftn('ijedi_mpas_state_write_flush_f90: mpas_get_time failed')

  call MPAS_stream_mgr_set_property(manager, trim(stream_name), STREAM_PROPERTY_FILENAME, &
                                    trim(filepath), ierr=ierr)
  if (ierr /= 0) call abor1_ftn('ijedi_mpas_state_write_flush_f90: set_property(filename) failed')

  ierr = 0
  call mpas_stream_mgr_write(manager, streamID=trim(stream_name), forceWriteNow=.true., &
                             writeTime=trim(dateTimeString), ierr=ierr)
  if (ierr /= 0) then
    write(message, '(a,i0)') 'ijedi_mpas_state_write_flush_f90: MPAS_stream_mgr_write failed ierr=', ierr
    call abor1_ftn(message)
  end if

end subroutine c_ijedi_mpas_state_write_flush

! ------------------------------------------------------------------------------

subroutine c_to_f_string(c_string, f_string)

  use iso_c_binding

  implicit none

  character(kind=c_char), intent(in) :: c_string(*)
  character(len=*), intent(out) :: f_string

  integer :: i

  f_string = ''
  do i = 1, len(f_string)
    if (c_string(i) == c_null_char) exit
    f_string(i:i) = transfer(c_string(i), 'a')
  end do

end subroutine c_to_f_string

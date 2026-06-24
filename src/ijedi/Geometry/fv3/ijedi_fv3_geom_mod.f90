module ijedi_fv3_geom_mod

use netcdf
use mpi
use string_f_c_mod

! atlas uses
use atlas_module,               only: atlas_field, atlas_fieldset, &
                                      atlas_integer, atlas_real, atlas_functionspace

! fckit uses
use fckit_mpi_module,           only: fckit_mpi_comm
use fckit_configuration_module, only: fckit_configuration

! fms uses
use fms_mod,                    only: fms_init
use mpp_mod,                    only: mpp_exit, mpp_pe, mpp_npes, mpp_error, FATAL, NOTE, &
                                      mpp_set_current_pelist
use mpp_domains_mod,            only: domain2D, mpp_deallocate_domain, mpp_define_layout, &
                                      mpp_define_mosaic, mpp_define_io_domain, mpp_domains_exit, &
                                      mpp_domains_set_stack_size
use ensemble_manager_mod,       only: get_ensemble_id, get_ensemble_size
use field_manager_mod,          only: fm_string_len, field_manager_init
use ensemble_manager_mod,       only: ensemble_manager_init, ensemble_pelist_setup
use ensemble_manager_mod,       only: get_ensemble_pelist

! fv3 uses
use ijedi_fv3_akbk_mod,       only: akbk_gfs_127
use ijedi_fv3_arrays_mod,     only: fv_atmos_type, deallocate_fv_atmos_type
use ijedi_fv3_control_mod,    only: fv_control_init

! ijedi uses
use ijedi_constants_mod,      only: constant
use ijedi_kinds_mod,          only: kind_int, kind_real
use ijedi_netcdf_utils_mod,   only: nccheck
use ijedi_fv3_namelist_mod,   only: ijedi_fmsnamelist

implicit none
private
public :: fv3_geom_initialize
public :: fv3_geom_create
public :: fv3_geom_set_and_fill_geometry_fields
public :: fv3_geom_setup_domain
public :: fv3_geom_write_geom
public :: fv3_geom_getVerticalCoord
public :: fv3_geom_getVerticalCoordLogP
public :: fv3_geom_pedges2pmidlayer
public :: fv3_geom_nodes_to_atlas_nodes

! Create interface for generic nodes to atlas nodes procedure
interface fv3_geom_nodes_to_atlas_nodes
  module procedure fv3_geom_nodes_to_atlas_nodes_r
  module procedure fv3_geom_nodes_to_atlas_nodes_i
end interface fv3_geom_nodes_to_atlas_nodes

! --------------------------------------------------------------------------------------------------


! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_initialize(conf, comm)

type(fckit_configuration), intent(in) :: conf
type(fckit_mpi_comm),      intent(in) :: comm

integer :: stackmax
character(len=1024) :: nml_filename, field_table_filename
character(len=:), allocatable :: str

! Path for input.nml file
call conf%get_or_die("namelist filename",str)
if (len(str) > 1024) call abor1_ftn("Length of fms namelist filename too long")
nml_filename = str
deallocate(str)

! Field table file
call conf%get_or_die("field table filename",str)
if (len(str) > fm_string_len) call abor1_ftn("Length of fms field table filename too long")
field_table_filename = str
deallocate(str)

! Initialize fms, mpp, etc.
call fms_init(localcomm=comm%communicator(), alt_input_nml_path = nml_filename)

! Set max stacksize
call conf%get_or_die("stackmax", stackmax)
call mpp_domains_set_stack_size(stackmax)

! Initialize the tracers
call field_manager_init(table_name = field_table_filename)

end subroutine fv3_geom_initialize

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_create(geom_conf, geom_vars, comm)

!Arguments
type(fckit_configuration), intent(in)    :: geom_conf
type(fckit_configuration), intent(inout) :: geom_vars
type(fckit_mpi_comm),      intent(in)    :: comm

!Locals
character(len=256)                    :: file_akbk
type(fv_atmos_type), allocatable      :: Atm(:)
logical, allocatable                  :: grids_on_this_pe(:)
integer                               :: i, j, jj, this_grid
integer                               :: p_split = 1
integer                               :: ncstat, ncid, akvarid, bkvarid, readdim, dcount
integer, dimension(nf90_max_var_dims) :: dimids, dimlens

integer :: npx, npy, npz, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, ntiles

character(len=:), allocatable :: str
real(kind=kind_real) :: sf, t_lon, t_lat
logical :: do_write_geom = .false.
integer :: iterator_dimension = 2

real(kind=kind_real), allocatable :: ak(:), bk(:)

type(ijedi_fmsnamelist) :: fmsnamelist

integer              :: ensNum
integer, allocatable :: Atm_pelist(:)
integer, allocatable :: Ocean_pelist(:)
integer, allocatable :: Land_pelist(:)
integer, allocatable :: Ice_fast_pelist(:)
integer              :: atmos_npes
integer              :: ocean_npes
integer              :: land_npes
integer              :: ice_npes
integer              :: ensemble_id
integer              :: ens_siz(6), ensemble_size, npes
integer, allocatable :: ensemble_pelist(:, :)

type(domain2D) :: domain

! Things for Atlas connections
integer :: ngrid, num_nodes, num_tri_elements, num_quad_elements
integer :: num_tri_boundary_nodes, num_quad_boundary_nodes

real(kind_real), allocatable :: lons(:)
real(kind_real), allocatable :: lats(:)
integer, allocatable :: ghosts(:)
integer, allocatable :: global_indices(:)
integer, allocatable :: remote_indices(:)
integer, allocatable :: partition(:)
integer, allocatable :: raw_tri_boundary_nodes(:)
integer, allocatable :: raw_quad_boundary_nodes(:)


! Update the fms name list with this Geometry
! -------------------------------------------
call fmsnamelist%replace_namelist(geom_conf)

! Initialize using the model setup routine
! ----------------------------------------
call fv_control_init(Atm, 300.0_kind_real, this_grid, grids_on_this_pe, p_split, &
                     skip_nml_read_in=.true.)

! For convenience copy grid integer things
npx = Atm(1)%npx
npy = Atm(1)%npy
npz = Atm(1)%npz
isc = Atm(1)%bd%isc
iec = Atm(1)%bd%iec
jsc = Atm(1)%bd%jsc
jec = Atm(1)%bd%jec
isd = Atm(1)%bd%isd
ied = Atm(1)%bd%ied
jsd = Atm(1)%bd%jsd
jed = Atm(1)%bd%jed
ntile = Atm(1)%global_tile
ntiles = Atm(1)%flagstruct%ntiles
ngrid = (Atm(1)%bd%iec-Atm(1)%bd%isc+1)*(Atm(1)%bd%jec-Atm(1)%bd%jsc+1)

! Sanity check
if (this_grid .ne. 1) call abor1_ftn("Geometry not ready for this_grid > 1")

! Copy relevant contents of Atm
! -----------------------------
call geom_vars%set("npx", npx)
call geom_vars%set("npy", npy)
call geom_vars%set("nLevels", npz)
call geom_vars%set("isc", isc)
call geom_vars%set("iec", iec)
call geom_vars%set("jsc", jsc)
call geom_vars%set("jec", jec)
call geom_vars%set("isd", isd)
call geom_vars%set("ied", ied)
call geom_vars%set("jsd", jsd)
call geom_vars%set("jed", jed)
call geom_vars%set("ntile", ntile)
call geom_vars%set("ntiles", ntiles)
call geom_vars%set("ngrid", ngrid)
call geom_vars%set("layout_x", Atm(1)%layout(1))
call geom_vars%set("layout_y", Atm(1)%layout(2))


!Allocatable arrays
allocate(ak(npz+1))
allocate(bk(npz+1))

! ak and bk hybrid coordinate coefficients
! ----------------------------------------
if (npz > 1) then

  ! Set path/filename for ak and bk file
  call geom_conf%get_or_die("akbk", str)
  file_akbk = str

  if (trim(file_akbk) == "gfs_127") then

    call akbk_gfs_127(npz, ak, bk)

  else

    !Open file
    call nccheck ( nf90_open(file_akbk, nf90_nowrite, ncid), "nf90_open "//file_akbk )

    !Search for ak in the file
    ncstat = nf90_inq_varid(ncid, "ak", akvarid)
    if(ncstat /= nf90_noerr) call abor1_ftn("Failed to find ak in file "//file_akbk)

    !Search for bk in the file
    ncstat = nf90_inq_varid(ncid, "bk", bkvarid)
    if(ncstat /= nf90_noerr) call abor1_ftn("Failed to find bk in file "//file_akbk)

    ! Check that dimension of ak/bk in the file match vertical levels of model
    dimids = 0
    call nccheck ( nf90_inquire_variable(ncid, akvarid, dimids = dimids), "nf90_inq_var ak" )
    readdim = -1
    dcount = 0
    do i = 1,nf90_max_var_dims
      if (dimids(i) > 0) then
         call nccheck( nf90_inquire_dimension(ncid, dimids(i), len = dimlens(i)), &
                       "nf90_inquire_dimension" )
         if (dimlens(i) == npz+1) then
            readdim = i
         endif
         dcount = dcount + 1
      endif
    enddo
    if (readdim == -1) call abor1_ftn("ak/bk in file does not match dimension of npz from input.nml")

    !Read ak and bk from the file
    call nccheck( nf90_get_var(ncid, akvarid, ak), "ijedi_fv3_geom, nf90_get_var ak" )
    call nccheck( nf90_get_var(ncid, bkvarid, bk), "ijedi_fv3_geom, nf90_get_var bk" )

  endif

else
  ak = 0.0_kind_real
  bk = 0.0_kind_real
endif

! Put ak/bk into the configuration for use in other places

call geom_vars%set("sigma_pressure_hybrid_coordinate_a_coefficient", ak)
call geom_vars%set("sigma_pressure_hybrid_coordinate_b_coefficient", bk)
call geom_vars%set("air_pressure_at_top_of_atmosphere_model", ak(1))

! Save some things later needed in Atlas-based Geometry Fields
! ------------------------------------------------------------
call geom_vars%set("area", reshape(Atm(1)%gridstruct%area_64(isc:iec, jsc:jec), (/ngrid/)))
call geom_vars%set("surface_pressure", reshape(Atm(1)%ps(isc:iec, jsc:jec), (/ngrid/)))
call geom_vars%set("surface_geopotential", reshape(Atm(1)%phis(isc:iec, jsc:jec), (/ngrid/)))

! Ensemble manager
! ----------------
if (.not. geom_conf%get("member_number", ensNum)) then
  ensNum = 0
endif

if( ensNum > 0 ) then
  call ensemble_manager_init()
  ens_siz = get_ensemble_size()
  ensemble_size = ens_siz(1)
  npes = ens_siz(2)

  atmos_npes = npes
  ocean_npes = 0
  land_npes = 0
  ice_npes = 0

  allocate( Atm_pelist  (atmos_npes) )
  allocate( Ocean_pelist(ocean_npes) )
  allocate( Land_pelist (land_npes) )
  allocate( Ice_fast_pelist(ice_npes) )

  call ensemble_pelist_setup(.true., atmos_npes, ocean_npes, land_npes, ice_npes, &
                               Atm_pelist, Ocean_pelist, Land_pelist, Ice_fast_pelist)
  ensemble_id = get_ensemble_id()
  allocate(ensemble_pelist(1:ensemble_size,1:npes))
  call get_ensemble_pelist(ensemble_pelist)
  call mpp_set_current_pelist(ensemble_pelist(ensemble_id,:))
  deallocate( Atm_pelist )
  deallocate( Ocean_pelist )
  deallocate( Land_pelist )
  deallocate( Ice_fast_pelist )
endif

! Place ensemble num in geom_vars
call geom_vars%set("ensNum", ensNum)

! Assert that ntiles is 1 or 6
if ((ntiles /= 6) .and. (ntiles /= 1)) then
  call mpp_error(FATAL, "get_num_nodes_and_elements: ntiles != 1 or 6")
endif

! Create fms domain for communication prior to creating atlas structures
call fv3_geom_setup_domain( domain, npx-1, npy-1, &
                            ntiles, Atm(1)%layout, Atm(1)%io_layout, 3 )

! Atlas connection requirements
! -----------------------------




! Get number of nodes and elements
if (ntiles == 6) then
  call fv3_geom_get_num_nodes_and_elements_global(Atm(1)%global_tile, &
                                                  Atm(1)%bd%isc, Atm(1)%bd%iec, &
                                                  Atm(1)%bd%jsc, Atm(1)%bd%jec, &
                                                  Atm(1)%npx, Atm(1)%npy, &
                                                  num_nodes, num_tri_elements, num_quad_elements)
else if (ntiles == 1) then
  call fv3_geom_get_num_nodes_and_elements_regional(Atm(1)%bd%isc, Atm(1)%bd%iec, &
                                                    Atm(1)%bd%jsc, Atm(1)%bd%jec, &
                                                    Atm(1)%npx, Atm(1)%npy, &
                                                    num_nodes, num_tri_elements, num_quad_elements)
end if

num_tri_boundary_nodes = 3 * num_tri_elements;
num_quad_boundary_nodes = 4 * num_quad_elements;

! Allocate arrays in Atlas form
allocate(lons(num_nodes))
allocate(lats(num_nodes))
allocate(ghosts(num_nodes))
allocate(global_indices(num_nodes))
allocate(remote_indices(num_nodes))
allocate(partition(num_nodes))
allocate(raw_tri_boundary_nodes(num_tri_boundary_nodes))
allocate(raw_quad_boundary_nodes(num_quad_boundary_nodes))

! Get arrays needed for atlas
if (ntiles == 6) then
  call fv3_geom_get_coords_and_connectivities_global(Atm(1)%bd%isc, Atm(1)%bd%iec, &
                                                     Atm(1)%bd%jsc, Atm(1)%bd%jec, &
                                                     Atm(1)%bd%isd, Atm(1)%bd%ied, &
                                                     Atm(1)%bd%jsd, Atm(1)%bd%jed, &
                                                     Atm(1)%npx, Atm(1)%npy, ngrid, &
                                                     Atm(1)%global_tile, &
                                                     Atm(1)%flagstruct%ntiles, &
                                                     Atm(1)%gridstruct%agrid_64(:,:,1), &
                                                     Atm(1)%gridstruct%agrid_64(:,:,2), &
                                                     domain, comm, &
                                                     num_nodes, &
                                                     num_tri_boundary_nodes, &
                                                     num_quad_boundary_nodes, &
                                                     lons, lats, ghosts, global_indices, &
                                                     remote_indices, partition, &
                                                     raw_tri_boundary_nodes, &
                                                     raw_quad_boundary_nodes)

else if (ntiles == 1) then
  call fv3_geom_get_coords_and_connectivities_regional(Atm(1)%bd%isc, Atm(1)%bd%iec, &
                                                       Atm(1)%bd%jsc, Atm(1)%bd%jec, &
                                                       Atm(1)%bd%isd, Atm(1)%bd%ied, &
                                                       Atm(1)%bd%jsd, Atm(1)%bd%jed, &
                                                       Atm(1)%npx, Atm(1)%npy, ngrid, &
                                                       Atm(1)%global_tile, &
                                                       Atm(1)%flagstruct%ntiles, &
                                                       Atm(1)%gridstruct%agrid_64(:,:,1), &
                                                       Atm(1)%gridstruct%agrid_64(:,:,2), &
                                                       domain, comm, &
                                                       num_nodes, &
                                                       num_tri_boundary_nodes, &
                                                       num_quad_boundary_nodes, &
                                                       lons, lats, ghosts, global_indices, &
                                                       remote_indices, partition, &
                                                       raw_tri_boundary_nodes, &
                                                       raw_quad_boundary_nodes)
end if

! Add everything to the geometry variables
call geom_vars%set("num_nodes", num_nodes)
call geom_vars%set("num_tri_elements",  num_tri_elements)
call geom_vars%set("num_quad_elements", num_quad_elements)
call geom_vars%set("lons", lons)
call geom_vars%set("lats", lats)
call geom_vars%set("ghosts", ghosts)
call geom_vars%set("global_indices", global_indices)
call geom_vars%set("remote_indices", remote_indices)
call geom_vars%set("partition", partition)
call geom_vars%set("raw_tri_boundary_nodes", raw_tri_boundary_nodes)
call geom_vars%set("raw_quad_boundary_nodes", raw_quad_boundary_nodes)

! Safe deallocate of the grid structure
! -------------------------------------
call deallocate_fv_atmos_type(Atm(1))
deallocate(Atm)
deallocate(grids_on_this_pe)

! Revert the fms namelist
! -----------------------
call fmsnamelist%revert_namelist

end subroutine fv3_geom_create

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_set_and_fill_geometry_fields(npz, ngrid, ak, bk, area, surface_pressure, &
                                                 surface_geopotential, vertcoord_selector, &
                                                 afunctionspace, afieldset)

!Arguments
integer,                   intent(in) :: npz, ngrid, vertcoord_selector
real(kind=kind_real),      intent(in) :: ak(npz+1), bk(npz+1)
real(kind=kind_real),      intent(in) :: area(ngrid)
real(kind=kind_real),      intent(in) :: surface_pressure(ngrid)
real(kind=kind_real),      intent(in) :: surface_geopotential(ngrid)
type(atlas_functionspace), intent(inout) :: afunctionspace
type(atlas_fieldset),      intent(inout) :: afieldset

!Locals
type(atlas_field) :: afield
integer :: jl, jn
integer, pointer :: int_ptr(:,:)
real(kind=kind_real), pointer :: real_ptr(:,:)
real(kind=kind_real) :: sigmaup, sigmadn, p_mid
real(kind=kind_real), parameter :: grav = 9.80665_kind_real

! Add owned vs halo/BC field
afield = afunctionspace%create_field(name='owned', kind=atlas_integer(kind_int), levels=1)
call afield%data(int_ptr)
int_ptr(1, :) = 0
int_ptr(1, 1:ngrid) = 1
call afieldset%add(afield)

! Add area
afield = afunctionspace%create_field(name='area', kind=atlas_real(kind_real), levels=1)
call afield%data(real_ptr)
real_ptr(1, :) = -1.0_kind_real
real_ptr(1, 1:ngrid) = area(:)
call afieldset%add(afield)

! Add vertical coordinate
if (vertcoord_selector == 1) then
   afield = afunctionspace%create_field(name='vert_coord', kind=atlas_real(kind_real), levels=npz)
   call afield%data(real_ptr)
   real_ptr(:, :) = 0.0_kind_real
   do jl=1,npz
      do jn=1,ngrid
         sigmaup = ak(jl+1)/surface_pressure(jn)+bk(jl+1)
         sigmadn = ak(jl  )/surface_pressure(jn)+bk(jl  )
         real_ptr(jl, jn) = 0.5_kind_real*(sigmaup+sigmadn)
      enddo
   enddo
else if (vertcoord_selector == 2) then
   afield = afunctionspace%create_field(name='vert_coord', kind=atlas_real(kind_real), levels=npz)
   call afield%data(real_ptr)
   real_ptr(:, :) = 0.0_kind_real
   do jl=1,npz
      do jn=1,ngrid
         p_mid = 0.5_kind_real*(ak(jl) + ak(jl+1)) + &
                 0.5_kind_real*(bk(jl) + bk(jl+1))*surface_pressure(jn)
         real_ptr(jl, jn) = log(p_mid)
      enddo
   enddo
else if (vertcoord_selector == 3) then
   afield = afunctionspace%create_field(name='vert_coord', kind=atlas_real(kind_real), levels=1)
   call afield%data(real_ptr)
   real_ptr(:, :) = 0.0_kind_real
   real_ptr(1, 1:ngrid) = surface_geopotential(:) / grav
else
   call abor1_ftn('ijedi_fv3_geom_mod%set_and_fill_geometry_fields: unknown vertical coordinate type')
endif
call afieldset%add(afield)
call afield%final()

end subroutine fv3_geom_set_and_fill_geometry_fields

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_setup_domain(domain, nx, ny, ntiles, layout_in, io_layout, halo)

 type(domain2D),   intent(inout) :: domain
 integer,          intent(in)    :: nx, ny, ntiles
 integer,          intent(in)    :: layout_in(:), io_layout(:)
 integer,          intent(in)    :: halo

 integer                              :: pe, npes, npes_per_tile, tile
 integer                              :: num_contact, num_alloc
 integer                              :: n, layout(2)
 integer, allocatable, dimension(:,:) :: global_indices, layout2D
 integer, allocatable, dimension(:)   :: pe_start, pe_end
 integer, allocatable, dimension(:)   :: tile1, tile2
 integer, allocatable, dimension(:)   :: istart1, iend1, jstart1, jend1
 integer, allocatable, dimension(:)   :: istart2, iend2, jstart2, jend2
 integer, allocatable :: tile_id(:), ensNum
 logical :: is_symmetry

  pe = mpp_pe()
  npes = mpp_npes()
  ensNum = get_ensemble_id()

  if (mod(npes,ntiles) /= 0) then
     call mpp_error(NOTE, "setup_domain: npes can not be divided by ntiles")
     return
  endif
  npes_per_tile = npes/ntiles
  tile = pe/npes_per_tile + 1

  if (layout_in(1)*layout_in(2) == npes_per_tile) then
     layout = layout_in
  else
     call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
  endif

  if (io_layout(1) <1 .or. io_layout(2) <1) call mpp_error(FATAL, &
          "setup_domain: both elements of variable io_layout must be positive integer")
  if (mod(layout(1), io_layout(1)) /= 0 ) call mpp_error(FATAL, &
       "setup_domain: layout(1) must be divided by io_layout(1)")
  if (mod(layout(2), io_layout(2)) /= 0 ) call mpp_error(FATAL, &
       "setup_domain: layout(2) must be divided by io_layout(2)")

  allocate(global_indices(4,ntiles), layout2D(2,ntiles), pe_start(ntiles), pe_end(ntiles) )

  ! select case based off of 1 or 6 tiles
  select case(ntiles)
  case ( 1 ) ! FV3-SAR
    num_contact = 0
  case ( 6 ) ! FV3 global
    num_contact = 12
  case default
    call mpp_error(FATAL, "setup_domain: ntiles != 1 or 6")
  end select

  do n = 1, ntiles
     global_indices(:,n) = (/1,nx,1,ny/)
     layout2D(:,n)       = layout
     pe_start(n)         = (n-1)*npes_per_tile + (ensNum -1) * 6 * npes_per_tile
     pe_end(n)           = n*npes_per_tile-1 + (ensNum -1) * 6 * npes_per_tile
  enddo

  num_alloc = max(1, num_contact)
  ! this code copied from domain_decomp in fv_mp_mod.f90
  allocate(tile1(num_alloc), tile2(num_alloc) )
  allocate(tile_id(ntiles))
  allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
  allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )
  ! select case based off of 1 or 6 tiles
  select case(ntiles)
  case ( 1 ) ! FV3-SAR
    ! No contacts, do nothing
  case ( 6 ) ! FV3 global
    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
    istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
    istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1; tile2(3) = 5
    istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
    istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
    istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2; tile2(5) = 3
    istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
    istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
    istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
    istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
    istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
    istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
    istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
    istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
  end select
  is_symmetry = .true.
  do n = 1, ntiles
     tile_id(n) = n
  enddo

  call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                         pe_start, pe_end, whalo=halo, ehalo=halo, shalo=halo, nhalo=halo,    &
                         symmetry=is_symmetry, tile_id=tile_id, &
                         name='cubic_grid')

  call mpp_define_io_domain(domain, io_layout)

  deallocate(pe_start, pe_end)
  deallocate(layout2D, global_indices)
  deallocate(tile1, tile2, tile_id)
  deallocate(istart1, iend1, jstart1, jend1)
  deallocate(istart2, iend2, jstart2, jend2)

end subroutine fv3_geom_setup_domain

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_write_geom(f_comm, isc, iec, jsc, jec, npx, npy, ntile, grid_lon, grid_lat, egrid_lon, &
                      egrid_lat)

  ! Arguments
  type(fckit_mpi_comm), intent(in) :: f_comm
  integer,              intent(in) :: isc, iec, jsc, jec, npx, npy, ntile
  real(kind=kind_real), intent(in) :: grid_lon(:,:), grid_lat(:,:)
  real(kind=kind_real), intent(in) :: egrid_lon(:,:), egrid_lat(:,:)

  ! Locals
  character(len=255) :: filename
  integer :: ncid, xf_dimid, yf_dimid, xv_dimid, yv_dimid, ti_dimid, pe_dimid
  integer :: mydims(3,3), ijdims(1), ijdimf(1), tmpij(1)
  integer :: varid(8)

  write(filename,"(A9,I0.4,A4)") 'fv3grid_c', npx-1, '.nc4'

  ! Create and open the file for parallel write
  call nccheck( nf90_create( trim(filename), ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                             comm = f_comm%communicator(), info = MPI_INFO_NULL), "nf90_create" )

  !Dimensions
  call nccheck ( nf90_def_dim(ncid, 'fxdim', npx-1   , xf_dimid), "nf90_def_dim fxdim" )
  call nccheck ( nf90_def_dim(ncid, 'fydim', npy-1   , yf_dimid), "nf90_def_dim fydim" )
  call nccheck ( nf90_def_dim(ncid, 'vxdim', npx     , xv_dimid), "nf90_def_dim vxdim" )
  call nccheck ( nf90_def_dim(ncid, 'vydim', npy     , yv_dimid), "nf90_def_dim vydim" )
  call nccheck ( nf90_def_dim(ncid, 'ntile', 6            , ti_dimid), "nf90_def_dim ntile" )
  call nccheck ( nf90_def_dim(ncid, 'nproc', f_comm%size(), pe_dimid), "nf90_def_dim ntile" )

  !Define variables
  call nccheck( nf90_def_var(ncid, "flons", NF90_DOUBLE, (/ xf_dimid, yf_dimid, ti_dimid /), varid(1)), "nf90_def_var flons" )
  call nccheck( nf90_put_att(ncid, varid(1), "long_name", "longitude of faces") )
  call nccheck( nf90_put_att(ncid, varid(1), "units", "degrees_east") )

  call nccheck( nf90_def_var(ncid, "flats", NF90_DOUBLE, (/ xf_dimid, yf_dimid, ti_dimid /), varid(2)), "nf90_def_var flats" )
  call nccheck( nf90_put_att(ncid, varid(2), "long_name", "latitude of faces") )
  call nccheck( nf90_put_att(ncid, varid(2), "units", "degrees_north") )

  call nccheck( nf90_def_var(ncid, "vlons", NF90_DOUBLE, (/ xv_dimid, yv_dimid, ti_dimid /), varid(3)), "nf90_def_var vlons" )
  call nccheck( nf90_put_att(ncid, varid(3), "long_name", "longitude of vertices") )
  call nccheck( nf90_put_att(ncid, varid(3), "units", "degrees_east") )

  call nccheck( nf90_def_var(ncid, "vlats", NF90_DOUBLE, (/ xv_dimid, yv_dimid, ti_dimid /), varid(4)), "nf90_def_var vlats" )
  call nccheck( nf90_put_att(ncid, varid(4), "long_name", "latitude of vertices") )
  call nccheck( nf90_put_att(ncid, varid(4), "units", "degrees_north") )

  call nccheck( nf90_def_var(ncid, "isc", NF90_INT, (/ pe_dimid /), varid(5)), "nf90_def_var isc" )
  call nccheck( nf90_put_att(ncid, varid(5), "long_name", "starting index i direction") )
  call nccheck( nf90_put_att(ncid, varid(5), "units", "1") )

  call nccheck( nf90_def_var(ncid, "iec", NF90_INT, (/ pe_dimid /), varid(6)), "nf90_def_var iec" )
  call nccheck( nf90_put_att(ncid, varid(6), "long_name", "ending index i direction") )
  call nccheck( nf90_put_att(ncid, varid(6), "units", "1") )

  call nccheck( nf90_def_var(ncid, "jsc", NF90_INT, (/ pe_dimid /), varid(7)), "nf90_def_var jsc" )
  call nccheck( nf90_put_att(ncid, varid(7), "long_name", "starting index j direction") )
  call nccheck( nf90_put_att(ncid, varid(7), "units", "1") )

  call nccheck( nf90_def_var(ncid, "jec", NF90_INT, (/ pe_dimid /), varid(8)), "nf90_def_var jec" )
  call nccheck( nf90_put_att(ncid, varid(8), "long_name", "ending index j direction") )
  call nccheck( nf90_put_att(ncid, varid(8), "units", "1") )

  ! End define mode
  call nccheck( nf90_enddef(ncid), "nf90_enddef" )

  ! Write variables
  mydims(1,1) = 1;          mydims(2,1) = npx-1
  mydims(1,2) = 1;          mydims(2,2) = npy-1
  mydims(1,3) = ntile; mydims(2,3) = 1

  call nccheck( nf90_put_var( ncid, varid(1), grid_lon(isc:iec,jsc:jec), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var flons" )

  call nccheck( nf90_put_var( ncid, varid(2), grid_lat(isc:iec,jsc:jec), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var flats" )

  mydims(1,1) = 1;          mydims(2,1) = npx
  mydims(1,2) = 1;          mydims(2,2) = npy
  mydims(1,3) = ntile; mydims(2,3) = 1

  call nccheck( nf90_put_var( ncid, varid(3), egrid_lon(isc:iec+1,jsc:jec+1), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var vlons" )

  call nccheck( nf90_put_var( ncid, varid(4), egrid_lat(isc:iec+1,jsc:jec+1), &
                              start = mydims(1,:), count = mydims(2,:) ), "nf90_put_var vlats" )

  ijdims(1) = f_comm%rank()+1
  ijdimf(1) = 1

  tmpij = isc
  call nccheck( nf90_put_var( ncid, varid(5), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var isc" )

  tmpij = iec
  call nccheck( nf90_put_var( ncid, varid(6), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var iec" )

  tmpij = jsc
  call nccheck( nf90_put_var( ncid, varid(7), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var jsc" )

  tmpij = jec
  call nccheck( nf90_put_var( ncid, varid(8), tmpij, start = ijdims, count = ijdimf ), "nf90_put_var jec" )

  ! Close the file
  call nccheck ( nf90_close(ncid), "nf90_close" )

end subroutine fv3_geom_write_geom

!----------------------------------------------------------------------------
! 1d pressure_edge to pressure_mid
!----------------------------------------------------------------------------

subroutine fv3_geom_pedges2pmidlayer(npz,ptype,pe1d,kappa,p1d)
 integer,              intent(in)  :: npz       !number of model layers
 character(len=*),     intent(in)  :: ptype     !midlayer pressure definition: 'average' or 'Philips'
 real(kind=kind_real), intent(in)  :: pe1d(npz+1) !pressure edge
 real(kind=kind_real), intent(in)  :: kappa
 real(kind=kind_real), intent(out) :: p1d(npz)    !pressure mid

 real(kind=kind_real) :: kap1, kapr

 kap1 = kappa + 1.0_kind_real
 kapr = 1.0_kind_real/kappa

 select case (ptype)
   case('Philips')
     p1d = ((pe1d(2:npz+1)**kap1 - pe1d(1:npz)**kap1)/&
            (kap1*(pe1d(2:npz+1) - pe1d(1:npz))))**kapr
   case default
     p1d = 0.5*(pe1d(2:npz+1) + pe1d(1:npz))
 end select

end subroutine fv3_geom_pedges2pmidlayer

!--------------------------------------------------------------------------------------------------
subroutine fv3_geom_getVerticalCoord(ak, bk, vc, npz, psurf)
  ! returns log(pressure) at mid level of the vertical column with surface
  ! prsssure of psurf
  ! coded using an example from Jeff Whitaker used in GSI ENKF pacakge

  real(kind=kind_real), intent(in) :: ak(npz+1)
  real(kind=kind_real), intent(in) :: bk(npz+1)
  integer,              intent(in) :: npz
  real(kind=kind_real), intent(in) :: psurf
  real(kind=kind_real), intent(out) :: vc(npz)

  real(kind=kind_real) :: plevli(npz+1), p(npz), kappa
  integer :: k

  ! compute interface pressure
  do k=1,npz+1
    plevli(k) = ak(k) + bk(k)*psurf
  enddo

  ! get kappa
  kappa = constant('kappa')

  ! compute presure at mid level and convert it to logp
  call fv3_geom_pedges2pmidlayer(npz,'Philips',plevli,kappa,vc)

end subroutine fv3_geom_getVerticalCoord

!--------------------------------------------------------------------------------------------------
subroutine fv3_geom_getVerticalCoordLogP(ak, bk, vc, npz, psurf)
  ! returns log(pressure) at mid level of the vertical column with surface prsssure of psurf
  ! coded using an example from Jeff Whitaker used in GSI ENKF pacakge

  real(kind=kind_real), intent(in) :: ak(npz+1)
  real(kind=kind_real), intent(in) :: bk(npz+1)
  integer,              intent(in) :: npz
  real(kind=kind_real), intent(in) :: psurf
  real(kind=kind_real), intent(out) :: vc(npz)

  real(kind=kind_real) :: p(npz)

  call fv3_geom_getVerticalCoord(ak, bk, p, npz, psurf)
  vc = - log(p)

end subroutine fv3_geom_getVerticalCoordLogP

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_get_num_nodes_and_elements(ntiles, ntile, isc, iec, jsc, jec, npx, npy, &
                                        num_nodes, num_tris, num_quads)

  integer, intent(in)  :: ntiles, ntile
  integer, intent(in)  :: isc, iec, jsc, jec
  integer, intent(in)  :: npx, npy
  integer, intent(out) :: num_nodes
  integer, intent(out) :: num_tris
  integer, intent(out) :: num_quads

  if (ntiles == 6) then
    call fv3_geom_get_num_nodes_and_elements_global(ntile, isc, iec, jsc, jec, npx, npy, &
                                           num_nodes, num_tris, num_quads)
  else if (ntiles == 1) then
    call fv3_geom_get_num_nodes_and_elements_regional(isc, iec, jsc, jec, npx, npy, &
                                             num_nodes, num_tris, num_quads)
  else
    call mpp_error(FATAL, "get_num_nodes_and_elements: ntiles != 1 or 6")
  end if

end subroutine fv3_geom_get_num_nodes_and_elements

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_get_num_nodes_and_elements_global(ntile, isc, iec, jsc, jec, npx, npy, &
                                               num_nodes, num_tris, num_quads)

  integer, intent(in)  :: ntile
  integer, intent(in)  :: isc, iec, jsc, jec
  integer, intent(in)  :: npx, npy
  integer, intent(out) :: num_nodes
  integer, intent(out) :: num_tris
  integer, intent(out) :: num_quads

  integer :: nx, ny
  logical :: lower_left_corner, upper_left_corner, lower_right_corner

  ! extra +1 from adding the ghost nodes on the lower side of each dimension
  nx = iec - isc + 2
  ny = jec - jsc + 2

  ! default case
  num_nodes = nx * ny
  num_tris = 0
  num_quads = (nx - 1) * (ny - 1)

  lower_left_corner = (isc == 1 .and. jsc == 1)
  upper_left_corner = (isc == 1 .and. jec == npy-1)
  lower_right_corner = (iec == npx-1 .and. jsc == 1)

  ! if at lower-left corner of any tile, then lower-left quad is a tri
  if (lower_left_corner) then
    num_nodes = num_nodes - 1
    num_tris = num_tris + 1
    num_quads = num_quads - 1
  end if

  ! if at upper-left corner of tile #3, then add extra tri in upper-left corner
  if (upper_left_corner .and. ntile == 3) then
    num_nodes = num_nodes + 1
    num_tris = num_tris + 1
  end if

  ! if at lower-right corner of tile #6, then add extra tri in lower-right corner
  if (lower_right_corner .and. ntile == 6) then
    num_nodes = num_nodes + 1
    num_tris = num_tris + 1
  end if

end subroutine fv3_geom_get_num_nodes_and_elements_global

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_get_num_nodes_and_elements_regional(isc, iec, jsc, jec, npx, npy, &
                                                 num_nodes, num_tris, num_quads)

  integer, intent(in)  :: isc, iec, jsc, jec
  integer, intent(in)  :: npx, npy
  integer, intent(out) :: num_nodes
  integer, intent(out) :: num_tris
  integer, intent(out) :: num_quads

  integer :: nx, ny
  logical :: right_bdry, upper_bdry

  ! extra +1 from adding the ghost nodes on the lower side of each dimension
  nx = iec - isc + 2
  ny = jec - jsc + 2

  right_bdry = (iec == npx-1)
  upper_bdry = (jec == npy-1)

  ! if at upper or right edges, need to adjust the nx,ny for a differently-sized rectangle
  if (right_bdry) then
    nx = nx + 1
  end if
  if (upper_bdry) then
    ny = ny + 1
  end if

  num_nodes = nx * ny
  num_tris = 0
  num_quads = (nx - 1) * (ny - 1)

end subroutine fv3_geom_get_num_nodes_and_elements_regional

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_get_coords_and_connectivities(ntiles, ntile, isc, iec, jsc, jec, isd, ied, jsd, jed, &
    npx, npy, ngrid, grid_lon, grid_lat, domain, f_comm, &
    num_nodes, num_tri_boundary_nodes, num_quad_boundary_nodes, &
    lons, lats, ghosts, global_indices, remote_indices, partition, &
    raw_tri_boundary_nodes, raw_quad_boundary_nodes)

  use mpp_domains_mod, only: domain2D

  integer, intent(in) :: ntiles, ntile
  integer, intent(in) :: isc, iec, jsc, jec
  integer, intent(in) :: isd, ied, jsd, jed
  integer, intent(in) :: npx, npy, ngrid
  real(kind_real), intent(in) :: grid_lon(isd:ied, jsd:jed)
  real(kind_real), intent(in) :: grid_lat(isd:ied, jsd:jed)
  type(domain2D), intent(inout) :: domain
  type(fckit_mpi_comm), intent(in) :: f_comm
  integer, intent(in) :: num_nodes
  integer, intent(in) :: num_tri_boundary_nodes
  integer, intent(in) :: num_quad_boundary_nodes
  real(kind_real), intent(out) :: lons(num_nodes)
  real(kind_real), intent(out) :: lats(num_nodes)
  integer, intent(out) :: ghosts(num_nodes)
  integer, intent(out) :: global_indices(num_nodes)
  integer, intent(out) :: remote_indices(num_nodes)
  integer, intent(out) :: partition(num_nodes)
  integer, intent(out) :: raw_tri_boundary_nodes(num_tri_boundary_nodes)
  integer, intent(out) :: raw_quad_boundary_nodes(num_quad_boundary_nodes)

  if (ntiles == 6) then
    call fv3_geom_get_coords_and_connectivities_global(isc, iec, jsc, jec, isd, ied, jsd, jed, &
        npx, npy, ngrid, ntile, ntiles, grid_lon, grid_lat, domain, f_comm, &
        num_nodes, num_tri_boundary_nodes, num_quad_boundary_nodes, &
        lons, lats, ghosts, global_indices, remote_indices, partition, &
        raw_tri_boundary_nodes, raw_quad_boundary_nodes)
  else if (ntiles == 1) then
    call fv3_geom_get_coords_and_connectivities_regional(isc, iec, jsc, jec, isd, ied, jsd, jed, &
        npx, npy, ngrid, ntile, ntiles, grid_lon, grid_lat, domain, f_comm, &
        num_nodes, num_tri_boundary_nodes, num_quad_boundary_nodes, &
        lons, lats, ghosts, global_indices, remote_indices, partition, &
        raw_tri_boundary_nodes, raw_quad_boundary_nodes)
  else
    call mpp_error(FATAL, "fv3_geom_get_coords_and_connectivities: ntiles != 1 or 6")
  end if

end subroutine fv3_geom_get_coords_and_connectivities

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_get_coords_and_connectivities_global(isc, iec, jsc, jec, isd, ied, jsd, jed, &
    npx, npy, ngrid, ntile, ntiles, grid_lon, grid_lat, domain, f_comm, &
    num_nodes, num_tri_boundary_nodes, num_quad_boundary_nodes, &
    lons, lats, ghosts, global_indices, remote_indices, partition, &
    raw_tri_boundary_nodes, raw_quad_boundary_nodes)

  use mpp_domains_mod, only: mpp_update_domains, domain2D

  integer, intent(in) :: isc, iec, jsc, jec
  integer, intent(in) :: isd, ied, jsd, jed
  integer, intent(in) :: npx, npy, ngrid, ntile, ntiles
  real(kind_real), intent(in) :: grid_lon(isd:ied, jsd:jed)
  real(kind_real), intent(in) :: grid_lat(isd:ied, jsd:jed)
  type(domain2D), intent(inout) :: domain
  type(fckit_mpi_comm), intent(in) :: f_comm
  integer, intent(in) :: num_nodes
  integer, intent(in) :: num_tri_boundary_nodes
  integer, intent(in) :: num_quad_boundary_nodes
  real(kind_real), intent(out) :: lons(num_nodes)
  real(kind_real), intent(out) :: lats(num_nodes)
  integer, intent(out) :: ghosts(num_nodes)
  integer, intent(out) :: global_indices(num_nodes)
  integer, intent(out) :: remote_indices(num_nodes)
  integer, intent(out) :: partition(num_nodes)
  integer, intent(out) :: raw_tri_boundary_nodes(num_tri_boundary_nodes)
  integer, intent(out) :: raw_quad_boundary_nodes(num_quad_boundary_nodes)

  integer :: i, j, node_counter, tri_counter, quad_counter
  logical :: lower_left_corner, upper_left_corner, lower_right_corner

  integer :: loc_ghost(isd:ied, jsd:jed)
  integer :: loc_global_index(isd:ied, jsd:jed)
  integer :: loc_remote_index(isd:ied, jsd:jed)
  integer :: loc_partition(isd:ied, jsd:jed)

  lower_left_corner = (isc == 1 .and. jsc == 1)
  upper_left_corner = (isc == 1 .and. jec == npy-1)
  lower_right_corner = (iec == npx-1 .and. jsc == 1)

  ! local 2d array for ghost, no need to exchange
  loc_ghost = 1
  loc_ghost(isc:iec, jsc:jec) = 0

  ! local 2d arrays for global_index, remote_index, and partition for exchanging across tasks
  loc_global_index = -1
  loc_global_index(isc:iec, jsc:jec) = (npx-1) * (npy-1) * (ntile-1)
  do j = jsc, jec
    do i = isc, iec
      ! 1-based index for global index
      loc_global_index(i,j) = loc_global_index(i,j) + (j - 1) * (npx-1) + i
    end do
  end do
  call mpp_update_domains(loc_global_index, domain)

  loc_remote_index = -1
  do j = jsc, jec
    do i = isc, iec
      ! 1-based index
      loc_remote_index(i,j) = (j - jsc) * (iec - isc + 1) + (i - isc) + 1
    end do
  end do
  call mpp_update_domains(loc_remote_index, domain)

  loc_partition = -1
  loc_partition(isc:iec, jsc:jec) = f_comm%rank()
  call mpp_update_domains(loc_partition, domain)

  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, grid_lon, lons)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, grid_lat, lats)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_ghost, ghosts)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_global_index, global_indices)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_remote_index, remote_indices)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_partition, partition)

  lons = constant('rad2deg') * lons
  lats = constant('rad2deg') * lats

  tri_counter = 1
  quad_counter = 1
  do j = jsc-1, jec
    do i = isc-1, iec

      ! if at lower-left corner of any tile, then lower-left quad is a tri => skip a point
      if (lower_left_corner .and. (j == jsc-1) .and. (i == isc-1)) then
        raw_tri_boundary_nodes(tri_counter)   = loc_global_index(i+1, j)
        raw_tri_boundary_nodes(tri_counter+1) = loc_global_index(i+1, j+1)
        raw_tri_boundary_nodes(tri_counter+2) = loc_global_index(i, j+1)
        tri_counter = tri_counter + 3
        cycle
      end if

      if ((j /= jec) .and. (i /= iec)) then
        raw_quad_boundary_nodes(quad_counter)   = loc_global_index(i, j)
        raw_quad_boundary_nodes(quad_counter+1) = loc_global_index(i+1, j)
        raw_quad_boundary_nodes(quad_counter+2) = loc_global_index(i+1, j+1)
        raw_quad_boundary_nodes(quad_counter+3) = loc_global_index(i, j+1)
        quad_counter = quad_counter + 4
      end if
    end do
  end do

  ! at upper-left corner of tile #3, then add extra tri => add extra point
  if (upper_left_corner .and. (ntile == 3)) then
    raw_tri_boundary_nodes(tri_counter)   = loc_global_index(isc-1, jec)
    raw_tri_boundary_nodes(tri_counter+1) = loc_global_index(isc, jec)
    raw_tri_boundary_nodes(tri_counter+2) = loc_global_index(isc, jec+1)
    tri_counter = tri_counter + 3
  end if

  ! if at lower-right corner of tile #6, then add extra tri => add extra point
  if (lower_right_corner .and. (ntile == 6)) then
    raw_tri_boundary_nodes(tri_counter)   = loc_global_index(iec, jsc-1)
    raw_tri_boundary_nodes(tri_counter+1) = loc_global_index(iec+1, jsc)
    raw_tri_boundary_nodes(tri_counter+2) = loc_global_index(iec, jsc)
    tri_counter = tri_counter + 3
  end if

  ! sanity checks: tri_counter-1 == num_tri_boundary_nodes
  if (tri_counter-1 /= num_tri_boundary_nodes) then
    call abor1_ftn('ijedi_fv3_geom_mod: inconsistent tri counter when getting connectivities')
  end if
  ! quad_counter-1 == num_quad_boundary_nodes
  if (quad_counter-1 /= num_quad_boundary_nodes) then
    call abor1_ftn('ijedi_fv3_geom_mod: inconsistent quad counter when getting connectivities')
  end if

end subroutine fv3_geom_get_coords_and_connectivities_global

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_get_coords_and_connectivities_regional(isc, iec, jsc, jec, isd, ied, jsd, jed, &
    npx, npy, ngrid, ntile, ntiles, grid_lon, grid_lat, domain, f_comm, &
    num_nodes, num_tri_boundary_nodes, num_quad_boundary_nodes, &
    lons, lats, ghosts, global_indices, remote_indices, partition, &
    raw_tri_boundary_nodes, raw_quad_boundary_nodes)

  use mpp_domains_mod, only: mpp_update_domains, domain2D

  integer, intent(in) :: isc, iec, jsc, jec
  integer, intent(in) :: isd, ied, jsd, jed
  integer, intent(in) :: npx, npy, ngrid, ntile, ntiles
  real(kind_real), intent(in) :: grid_lon(isd:ied, jsd:jed)
  real(kind_real), intent(in) :: grid_lat(isd:ied, jsd:jed)
  type(domain2D), intent(inout) :: domain
  type(fckit_mpi_comm), intent(in) :: f_comm
  integer, intent(in) :: num_nodes
  integer, intent(in) :: num_tri_boundary_nodes
  integer, intent(in) :: num_quad_boundary_nodes
  real(kind_real), intent(out) :: lons(num_nodes)
  real(kind_real), intent(out) :: lats(num_nodes)
  integer, intent(out) :: ghosts(num_nodes)
  integer, intent(out) :: global_indices(num_nodes)
  integer, intent(out) :: remote_indices(num_nodes)
  integer, intent(out) :: partition(num_nodes)
  integer, intent(out) :: raw_tri_boundary_nodes(num_tri_boundary_nodes)
  integer, intent(out) :: raw_quad_boundary_nodes(num_quad_boundary_nodes)

  integer :: i, j, node_counter, quad_counter
  integer :: imax, jmax
  integer :: counter_local_idx
  logical :: left_bdry, right_bdry, lower_bdry, upper_bdry

  logical :: loc_bc(isd:ied, jsd:jed)
  integer :: loc_ghost(isd:ied, jsd:jed)
  integer :: loc_global_index(isd:ied, jsd:jed)
  integer :: loc_remote_index(isd:ied, jsd:jed)
  integer :: exchange_remote_index(isd:ied, jsd:jed)
  integer :: loc_partition(isd:ied, jsd:jed)

  left_bdry = (isc == 1)
  right_bdry = (iec == npx-1)
  lower_bdry = (jsc == 1)
  upper_bdry = (jec == npy-1)

  imax = iec
  if (right_bdry) then
    imax = imax + 1
  end if

  jmax = jec
  if (upper_bdry) then
    jmax = jmax + 1
  end if

  ! identify which points in first halo layer are this tasks's BC points
  loc_bc = .false.
  if (left_bdry) loc_bc(isc-1, jsc:jec) = .true.
  if (right_bdry) loc_bc(iec+1, jsc:jec) = .true.
  if (lower_bdry) loc_bc(isc:iec, jsc-1) = .true.
  if (upper_bdry) loc_bc(isc:iec, jec+1) = .true.
  if (left_bdry .and. lower_bdry) loc_bc(isc-1, jsc-1) = .true.
  if (left_bdry .and. upper_bdry) loc_bc(isc-1, jec+1) = .true.
  if (right_bdry .and. lower_bdry) loc_bc(iec+1, jsc-1) = .true.
  if (right_bdry .and. upper_bdry) loc_bc(iec+1, jec+1) = .true.

  ! local 2d array for ghost, no need to exchange
  loc_ghost = 1
  loc_ghost(isc:iec, jsc:jec) = 0
  where (loc_bc) loc_ghost = 0

  ! local 2d arrays for global_index, remote_index, and partition for exchanging across tasks

  ! global_index runs over the entire regional "compute" domain +/- 1 point
  loc_global_index = -1
  do j = jsc-1, jec+1
    do i = isc-1, iec+1
      ! 1-based index
      loc_global_index(i,j) = j * (npx + 1) + i + 1
    end do
  end do

  loc_remote_index = -1
  do j = jsc, jec
    do i = isc, iec
      ! 1-based index
      loc_remote_index(i,j) = (j - jsc) * (iec - isc + 1) + (i - isc) + 1
    end do
  end do
  counter_local_idx = maxval(loc_remote_index)
  ! use exchange to fill halo points with neighboring task's index
  call mpp_update_domains(loc_remote_index, domain)
  ! for halo points that are actually a BC, generate new local indices
  do j = jsc-1, jec+1
    do i = isc-1, iec+1
      if (loc_bc(i, j)) then
        counter_local_idx = counter_local_idx + 1
        loc_remote_index(i,j) = counter_local_idx
      end if
    end do
  end do

  loc_partition = -1
  loc_partition(isc:iec, jsc:jec) = f_comm%rank()
  call mpp_update_domains(loc_partition, domain)
  where (loc_bc) loc_partition = f_comm%rank()

  ! special case handling of halo points within the BC region:
  !
  ! to allow JEDI's generic code (using atlas) to perform halo exchanges within the BC region, we
  ! need to pass connectivity information (atlas's partition number and index of each point)
  ! between adjacent processors on the boundary. this is slightly tedious to do, because fv3's own
  ! halo-exchanges do NOT pass around BC information.
  !
  ! our strategy is take the indices that were generated to fill loc_remote_index, copy them into
  ! the *owned* portion of a dummy array, then use fv3's halo exchanges to send them to neighboring
  ! MPI tasks. this adds an extra MPI communication, but avoids the need to implement complicated
  ! logic to recreated the neighboring task's locally-generated indices into its BC regions.
  !
  ! in this illustration,
  !
  ! boundary condition        bc0   bc1   bc2   bc3 | bc4   bc5   bc6   bc7
  !                                                 |
  ! upper boundary of domain  ----------------------+----------------------
  !                                                 |
  ! interior points           x0    x1    x2    x3  | y0    y1    y2    y3
  !                                 task 0          |       task 1
  !
  ! task 0 needs to fill the halo location `bc4` adjacent to its BC location `bc3`. the partition
  ! number is easily obtained by reading the partition of points `y0`, as these must match.
  ! however, the index of `bc4` on task 1 is hard to recompute on task 0, so we follow these steps,
  ! - task 1 copies the generated index of `bc4` into the `y0` slot of a dummy array
  ! - use fv3's halo exchange to send this index into the halo region on task 0
  ! - task 0 reads `bc4`'s local index (from task 1) in the `y0` slot of the dummy array
  exchange_remote_index = -1
  if (left_bdry) then
    if (.not.upper_bdry) then
      ! fill upper-left OWNED point with generated index of corresponding BC point
      exchange_remote_index(isc, jec) = loc_remote_index(isc-1, jec)
    end if
    if (.not.lower_bdry) then
      exchange_remote_index(isc, jsc) = loc_remote_index(isc-1, jsc)
    end if
  end if
  if (right_bdry) then
    if (.not.upper_bdry) then
      exchange_remote_index(iec, jec) = loc_remote_index(iec+1, jec)
    end if
    if (.not.lower_bdry) then
      exchange_remote_index(iec, jsc) = loc_remote_index(iec+1, jsc)
    end if
  end if
  if (lower_bdry) then
    if (.not.left_bdry) then
      exchange_remote_index(isc, jsc) = loc_remote_index(isc, jsc-1)
    end if
    if (.not.right_bdry) then
      exchange_remote_index(iec, jsc) = loc_remote_index(iec, jsc-1)
    end if
  end if
  if (upper_bdry) then
    if (.not.left_bdry) then
      exchange_remote_index(isc, jec) = loc_remote_index(isc, jec+1)
    end if
    if (.not.right_bdry) then
      exchange_remote_index(iec, jec) = loc_remote_index(iec, jec+1)
    end if
  end if
  ! halo-exchange the dummy array
  call mpp_update_domains(exchange_remote_index, domain)
  if (left_bdry) then
    if (.not.upper_bdry) then
      ! fill upper-left BC from neighbor's lower-right OWNED point, using the index already in halo
      loc_remote_index(isc-1, jec+1) = exchange_remote_index(isc, jec+1)
      loc_partition(isc-1, jec+1) = loc_partition(isc, jec+1)
    end if
    if (.not.lower_bdry) then
      loc_remote_index(isc-1, jsc-1) = exchange_remote_index(isc, jsc-1)
      loc_partition(isc-1, jsc-1) = loc_partition(isc, jsc-1)
    end if
  end if
  if (right_bdry) then
    if (.not.upper_bdry) then
      loc_remote_index(iec+1, jec+1) = exchange_remote_index(iec, jec+1)
      loc_partition(iec+1, jec+1) = loc_partition(iec, jec+1)
    end if
    if (.not.lower_bdry) then
      loc_remote_index(iec+1, jsc-1) = exchange_remote_index(iec, jsc-1)
      loc_partition(iec+1, jsc-1) = loc_partition(iec, jsc-1)
    end if
  end if
  if (lower_bdry) then
    if (.not.left_bdry) then
      loc_remote_index(isc-1, jsc-1) = exchange_remote_index(isc-1, jsc)
      loc_partition(isc-1, jsc-1) = loc_partition(isc-1, jsc)
    end if
    if (.not.right_bdry) then
      loc_remote_index(iec+1, jsc-1) = exchange_remote_index(iec+1, jsc)
      loc_partition(iec+1, jsc-1) = loc_partition(iec+1, jsc)
    end if
  end if
  if (upper_bdry) then
    if (.not.left_bdry) then
      loc_remote_index(isc-1, jec+1) = exchange_remote_index(isc-1, jec)
      loc_partition(isc-1, jec+1) = loc_partition(isc-1, jec)
    end if
    if (.not.right_bdry) then
      loc_remote_index(iec+1, jec+1) = exchange_remote_index(iec+1, jec)
      loc_partition(iec+1, jec+1) = loc_partition(iec+1, jec)
    end if
  end if

  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, grid_lon, lons)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, grid_lat, lats)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_ghost, ghosts)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_global_index, global_indices)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_remote_index, remote_indices)
  call fv3_geom_nodes_to_atlas_nodes(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                ntiles, ngrid, loc_partition, partition)

  lons = constant('rad2deg') * lons
  lats = constant('rad2deg') * lats

  quad_counter = 1
  do j = jsc-1, jmax
    do i = isc-1, imax
      if ((j /= jmax) .and. (i /= imax)) then
        raw_quad_boundary_nodes(quad_counter)   = loc_global_index(i, j)
        raw_quad_boundary_nodes(quad_counter+1) = loc_global_index(i+1, j)
        raw_quad_boundary_nodes(quad_counter+2) = loc_global_index(i+1, j+1)
        raw_quad_boundary_nodes(quad_counter+3) = loc_global_index(i, j+1)
        quad_counter = quad_counter + 4
      end if
    end do
  end do

  ! quad_counter-1 == num_quad_boundary_nodes
  if (quad_counter-1 /= num_quad_boundary_nodes) then
    call abor1_ftn('ijedi_fv3_geom_mod: inconsistent quad counter when getting connectivities')
  end if

  ! Avoid compilation warning
  raw_tri_boundary_nodes = 0

end subroutine fv3_geom_get_coords_and_connectivities_regional

! --------------------------------------------------------------------------------------------------

subroutine fv3_geom_nodes_to_atlas_nodes_r(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                      ntiles, ngrid, fv3_data, atlas_data)

  integer, intent(in) :: npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, ntiles, ngrid
  real(kind_real), intent(in) :: fv3_data(isd:ied, jsd:jed)
  real(kind_real), intent(inout) :: atlas_data(:)

  integer :: a, b, ncopy
  logical :: at_lower_left_corner, at_upper_left_corner, at_lower_right_corner
  logical :: at_right_edge, at_upper_edge
  logical :: halo_w, halo_e, halo_s, halo_n, halo_sw, halo_nw, halo_ne, halo_se, halo_nw3, halo_se6

  ! Identify which halos need including
  ! Default case for PEs interior to a tile
  halo_w = .true.
  halo_s = .true.
  halo_sw = .true.
  halo_e = .false.
  halo_n = .false.
  halo_nw = .false.
  halo_ne = .false.
  halo_se = .false.
  halo_nw3 = .false.
  halo_se6 = .false.

  ! Edges and corners depend on specifics...
  if (ntiles == 6) then
    ! Global grid -- handle corners between cubed-sphere tiles
    at_lower_left_corner = (isc == 1 .and. jsc == 1)
    at_upper_left_corner = (isc == 1 .and. jec == npy-1)
    at_lower_right_corner = (iec == npx-1 .and. jsc == 1)

    ! at lower-left corner of any tile, use a triangle => no diagonal point
    if (at_lower_left_corner) then
      halo_sw = .false.
    end if
    ! at upper-left corner of tile #3, place extra tri => add extra point
    if (at_upper_left_corner .and. (ntile == 3)) then
      halo_nw3 = .true.
    end if
    ! at lower-right corner of tile #6, place extra tri => add extra point
    if (at_lower_right_corner .and. (ntile == 6)) then
      halo_se6 = .true.
    end if

  else if (ntiles == 1) then
    ! Regional grid -- handle "boundary condition" points around patch
    at_right_edge = (iec == npx-1)
    at_upper_edge = (jec == npy-1)

    if (at_upper_edge) then
      halo_n = .true.
      halo_nw = .true.
    end if
    if (at_right_edge) then
      halo_e = .true.
      halo_se = .true.
      if (at_upper_edge) then
        halo_ne = .true.
      end if
    end if

  else
    call mpp_error(FATAL, "fv3_nodes_to_atlas_nodes: ntiles != 1 or 6")
  end if

  ! First, copy owned points
  ncopy = ngrid
  a = 1
  b = ncopy
  atlas_data(a:b) = reshape(fv3_data(isc:iec, jsc:jec), (/ncopy/))

  ! Copy west + east edge halos
  ncopy = (jec - jsc + 1)
  if (halo_w) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(isc-1, jsc:jec), (/ncopy/))
  end if
  if (halo_e) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(iec+1, jsc:jec), (/ncopy/))
  end if

  ! Copy south + north edge halos
  ncopy = (iec - isc + 1)
  if (halo_s) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(isc:iec, jsc-1), (/ncopy/))
  end if
  if (halo_n) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(isc:iec, jec+1), (/ncopy/))
  end if

  ! Copy corners
  if (halo_sw) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(isc-1, jsc-1)
  end if
  if (halo_nw) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(isc-1, jec+1)
  end if
  if (halo_ne) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(iec+1, jec+1)
  end if
  if (halo_se) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(iec+1, jsc-1)
  end if

  if (halo_nw3) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(isc, jec+1)
  end if
  if (halo_se6) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(iec+1, jsc)
  end if

  ! sanity check on size: b = size(atlas_data)
  if (b /= size(atlas_data)) then
    call abor1_ftn('fv3jedi_geom_mod%fv3_nodes_to_atlas_nodes: inconsistent atlas_data size')
  end if

end subroutine fv3_geom_nodes_to_atlas_nodes_r

! --------------------------------------------------------------------------------------------------

! displeasing!
! this is a copy of the real interface above with just one replacement real -> integer
subroutine fv3_geom_nodes_to_atlas_nodes_i(npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, &
                                      ntiles, ngrid, fv3_data, atlas_data)

  integer, intent(in) :: npx, npy, isc, iec, jsc, jec, isd, ied, jsd, jed, ntile, ntiles, ngrid
  integer, intent(in) :: fv3_data(isd:ied, jsd:jed)
  integer, intent(inout) :: atlas_data(:)

  integer :: a, b, ncopy
  logical :: at_lower_left_corner, at_upper_left_corner, at_lower_right_corner
  logical :: at_right_edge, at_upper_edge
  logical :: halo_w, halo_e, halo_s, halo_n, halo_sw, halo_nw, halo_ne, halo_se, halo_nw3, halo_se6

  ! Identify which halos need including
  ! Default case for PEs interior to a tile
  halo_w = .true.
  halo_s = .true.
  halo_sw = .true.
  halo_e = .false.
  halo_n = .false.
  halo_nw = .false.
  halo_ne = .false.
  halo_se = .false.
  halo_nw3 = .false.
  halo_se6 = .false.

  ! Edges and corners depend on specifics...
  if (ntiles == 6) then
    ! Global grid -- handle corners between cubed-sphere tiles
    at_lower_left_corner = (isc == 1 .and. jsc == 1)
    at_upper_left_corner = (isc == 1 .and. jec == npy-1)
    at_lower_right_corner = (iec == npx-1 .and. jsc == 1)

    ! at lower-left corner of any tile, use a triangle => no diagonal point
    if (at_lower_left_corner) then
      halo_sw = .false.
    end if
    ! at upper-left corner of tile #3, place extra tri => add extra point
    if (at_upper_left_corner .and. (ntile == 3)) then
      halo_nw3 = .true.
    end if
    ! at lower-right corner of tile #6, place extra tri => add extra point
    if (at_lower_right_corner .and. (ntile == 6)) then
      halo_se6 = .true.
    end if

  else if (ntiles == 1) then
    ! Regional grid -- handle "boundary condition" points around patch
    at_right_edge = (iec == npx-1)
    at_upper_edge = (jec == npy-1)

    if (at_upper_edge) then
      halo_n = .true.
      halo_nw = .true.
    end if
    if (at_right_edge) then
      halo_e = .true.
      halo_se = .true.
      if (at_upper_edge) then
        halo_ne = .true.
      end if
    end if

  else
    call mpp_error(FATAL, "fv3_nodes_to_atlas_nodes: ntiles != 1 or 6")
  end if

  ! First, copy owned points
  ncopy = ngrid
  a = 1
  b = ncopy
  atlas_data(a:b) = reshape(fv3_data(isc:iec, jsc:jec), (/ncopy/))

  ! Copy west + east edge halos
  ncopy = (jec - jsc + 1)
  if (halo_w) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(isc-1, jsc:jec), (/ncopy/))
  end if
  if (halo_e) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(iec+1, jsc:jec), (/ncopy/))
  end if

  ! Copy south + north edge halos
  ncopy = (iec - isc + 1)
  if (halo_s) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(isc:iec, jsc-1), (/ncopy/))
  end if
  if (halo_n) then
    a = b + 1
    b = b + ncopy
    atlas_data(a:b) = reshape(fv3_data(isc:iec, jec+1), (/ncopy/))
  end if

  ! Copy corners
  if (halo_sw) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(isc-1, jsc-1)
  end if
  if (halo_nw) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(isc-1, jec+1)
  end if
  if (halo_ne) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(iec+1, jec+1)
  end if
  if (halo_se) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(iec+1, jsc-1)
  end if

  if (halo_nw3) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(isc, jec+1)
  end if
  if (halo_se6) then
    a = b + 1
    b = b + 1
    atlas_data(a) = fv3_data(iec+1, jsc)
  end if

  ! sanity check on size: b = size(atlas_data)
  if (b /= size(atlas_data)) then
    call abor1_ftn('ijedi_fv3_geom_mod%fv3_nodes_to_atlas_nodes: inconsistent atlas_data size')
  end if

end subroutine fv3_geom_nodes_to_atlas_nodes_i

! --------------------------------------------------------------------------------------------------

end module ijedi_fv3_geom_mod

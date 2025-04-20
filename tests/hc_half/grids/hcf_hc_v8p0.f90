!*******************************************************************************
!
!          --- Grid generation code for a hemisphere cylinder ---
!
!
! This code generates a 3D grid over a hemisphere cylinder.
!
! Version 8.0 (June 29, 2017).
!
!------------------------------------------------------------------------------
!
!-------------
!  Input:  Parameters are set in the input file named as 'input.nml'.
!         See 'Input parameter module' below.
!
!   Note: The geometry of HC used in the experiment Ref[1]: Radius of the hemisphere is 0.5 inches, and the
!         body length is 10 inches, which correspond to x2 = 10.0 (the radius is 0.5 by default). 
!
!         [1] Tsieh, T., An Investigation of Separated Flow About a Hemisphere-Cylinder at 0- to 19-Deg
!             Incidence in the Mach Number Range from 0.6 to 1.5, AEDC-TR-76-112, 1976.
!
!-------------
! Output: 2 files by default:
!
!   (1) Grid file = "hc_'element_type'.1.(b8/l8)ugrid"
!
!       Note: Unstructured grid file (for all grid types).
!             Sometimes, you don't want this file.
!
!   (2) Grid file = "hc_'element_type'.1.p3d/ufmt" and "hc_'element_type'.1.nmf"
!
!       Note: Structured PLOT3D grid file ("igrid_type = 5" only).
!             Sometimes, you don't want these files.
!
!   (3) Structured indices = "hc_'element_type'.1.k" (for all types of elements).
!
!        Set generate_k_file = T to write this file.
!
!       Note: This .k file is required in 'hcf_coarsening_v3p3.f90' for
!             generating a series of coarse grids.
!
!   (4) Node line-info( BL) = "hc_'element_type'.1.lines_fmt"     - Points in lines in the BL region.
!   (5) Node line-info(All) = "hc_'element_type'.1.lines_fmt_all" - Points in lines all the way to outer.
!                        (NOTE: Rename .lines_fmt_all as .lines_fmt, for use in FUN3D.)
!
!   (6) Cell line-info( BL) = "hc_'element_type'.1.lines_fmt_cc"     - Cells in lines in the BL region.
!   (7) Cell line-info(All) = "hc_'element_type'.1.lines_fmt_cc_all" - Cells in lines all the way to outer.
!                        (NOTE: Rename .lines_fmt_all as .lines_fmt, for use in FUN3D.)
!
!       Note: Line information can be used for line relaxation or line agglomeration.
!       Note: These line files are required in 'hcf_coarsening_v3p3.f90' for
!             generating a series of coarse grids.
!       Note: Element(cell) numbers are defined in the order: tetra -> prism -> hex.
!
!   (8)"hc_'element_type'.1.mapbc"        !Boundary condition file for FUN3D
!   (9) Tecplot file = "hc_'element_type'.1.tec_bndary.dat" !Boundary grid for viewing
!  (10) Tecplot file = "hc_'element_type'.1.tec_volume.dat" !Volume grid for viewing
!
!
! [Send comments/bug-report to Hiro at hiro(at)nianet.org.]
!
!*******************************************************************************


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!  Input parameter module
!
!-------------------------------------------------------------------------------
!
!  Sample input file: 'input.nml' to generate a tet grid for HC.
!  ------------------------------------------------------
!  &input_parameters
!          target_reynolds_number = 350000
!                   target_y_plus = 0.5
!                      igrid_type = 5
!                              x2 = 10.0
!                         R_outer = 100
!           nnodes_cylinder_input = 256
!                           nr_gs = 64
!                             nre = 1280
!                      domain_cut = 2
!             generate_ugrid_file = T
!          ugrid_file_unformatted = T
!          generate_p3d_ufmt_file = T
!                 generate_k_file = T
!              k_file_unformatted = T
!           generate_line_file_nc = T
!           generate_line_file_cc = T
!             generate_tec_file_b = T
!             generate_tec_file_v = F
!  /
!  ------------------------------------------------------
!
!  Note: No need to specify all namelist variables.
!        In the above, those not shown are given their default values
!        as defined below.
!
!   Note: 'nodes_cylinder_input', 'nr_gs', and 'nre' determine the final grid size
!          (total # of nodes).
!
!   Note: If you want to regularly coarsen the grid afterwards (by 'hcf_coarsening_v3p3.f90'),
!         make sure that 'nr_gs' and 'nodes_cylinder_input' are 
!         multiples of 2 (e.g., 4, 8, 16, 32, 64, ...), and also set nr = a multiple of 2.
!
!  NOTE: If you wish to generate coarser grids later by regular coarsening, you must set
!
!        generate_ugrid_file   = T
!        generate_line_file_nc = T
!        generate_line_file_cc = T
!        generate_k_file       = T
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
 module input_parameter_module

  implicit none

  integer , parameter ::     dp = selected_real_kind(P=15)

  public

!----------------------------
! Default input values
!----------------------------

!----------------------------
! debug_mode = T: for debug use only.
!              F: Default
!
  logical :: debug_mode = .false.

!----------------------------
! target_reynolds_number = Reynolds number to determine the mesh spacing at wall,
!                                based on the root chord: Re=rho*U*(root chord)/mu.

  real(dp) :: target_reynolds_number = 1.0e7

!----------------------------
! target_y_plus = y-plus to determine the mesh spacing at wall (e.g., 1.0).

  real(dp) :: target_y_plus = 10.0

!----------------------------
!  igrid_type = Element type as listed below.
!                            1 = Prismatic
!                            2 = Tetrahedral
!                            3 = Mixed - Prism in BL and Tetra outside
!                            4 = Mixed - Prism on the hemisphere and Hexa over the cylinder.
!                            5 = Strctured - Hexa except prisms around the center of the rounded tip.

  integer :: igrid_type = 1

!----------------------------
! x2 = Length of the hemisphere-cylinder.
!   Apex of the hemisphere is at x=0, this is the position of
!   the end of the cylinder, which is equal to the length of HC.

  real(dp) :: x2 = 10.0_dp

!----------------------------
! R_outer = Radius of the outer-boundary hemisphere.

  real(dp) :: R_outer = 100.0_dp

!----------------------------
! nodes_cylinder_input = Number of elements along the cylinder part.
!
! Note: If this value is negative, the number of elements is set by
!        --> int( (x2-x1)/s )/abs(nnodes_cylinder_input), s is a spacing over the hemisphere.
! So, if nnodes_cylinder_input=-1, then we get a uniform spacing over the cylinder;
! and if nnodes_cylinder_input=-2, less nodes and stretched to match the spacing with hemisphere grid.

  integer :: nnodes_cylinder_input = 24

!----------------------------
! nr_gs = # of Elements across the rounded tip /2 (-> HC: x=0 to 0.5)

  integer :: nr_gs = 8

!----------------------------
! nre = # of Elements in the radial direction (from HC to farfield)
!   Note: If a non-positive value is given for 'nre' (number of spacings in the radial
!         direction), then it will be automatically determined in the code.

  integer :: nre = 88

!----------------------------
! BL region height factor: BL region is defined as
!  bl_height_factor x (theoretical estimate).

  real(dp) :: bl_height_factor = 0.5_dp

!----------------------------
! Lifting factor for the first-off-the-wall node at TE to take a larger
! distance to avoid vanishing element volumes.

  real(dp) :: lifting_factor_node_te = 3.0_dp

!----------------------------
! To use the same stretching factor at all radial lines (true) or not (false).

  logical :: uniform_stretching_to_outer = .true.

!----------------------------
! spacing_ratio = Transition from hemisphere to cylinder in mesh spacing
!                 E.g., spacing_ratio=2 gives the spacing of the first two
!                 nodes in the cylinder part twice as large as the mesh spacing
!                 in the hemisphere surface grid.
!
! spacing_ratio_outer = the ratio of the radial spacing to the surface spacing of
!                       the outer boundary. Cells will be isotropic if this is 1.0,
!                        but it will create a lot of cells in the farfield.

  real(dp) :: spacing_ratio       = 1.1_dp
  real(dp) :: spacing_ratio_outer = 1.7_dp

!----------------------------
! 1: full geometry, 2: a half domain (y > 0 only).

  integer :: domain_cut = 1

!----------------------------
! generate_ugrid_file = T to write .ugrid file (required by the coarsening program).
!                       F not to write.
! >>> [=T Required by the coarsening program]

  logical :: generate_ugrid_file = .true.

!----------------------------
! ugrid_file_unformatted = T: unformatted, F: formatted

  logical :: ugrid_file_unformatted = .true.

!----------------------------
! generate_p3d_ufmt_file = T = Write .p3d/.umft and .nmf files (igrid_type = 5 only); F: don't write.

  logical :: generate_p3d_ufmt_file = .true.

!----------------------------
! generate_k_file = T = Write a 'k'-file; F: don't write.
! >>> [=T Required by the coarsening program]

  logical :: generate_k_file = .true.

!----------------------------
! k_file_unformatted = T: unformatted, F: formatted

  logical :: k_file_unformatted = .true.

!----------------------------
! generate_line_file_nc = T or F.
! >>> [=T Required by the coarsening program]

  logical :: generate_line_file_nc = .true.

!----------------------------
! generate_line_file_cc = T or F.
! >>> [=T Required by the coarsening program]

  logical :: generate_line_file_cc = .true.

!----------------------------
! generate_tec_file_b = T: Write a Tecplot file, xxx_'element_type'.1.tec_bndary.dat

  logical :: generate_tec_file_b = .false.

!----------------------------
! generate_tec_file_v = T: Write a Tecplot file, xxx_'element_type'.1.tec_volume.dat

  logical :: generate_tec_file_v = .false.

!----------------------------
! Compute grid statistics (may require large memory...)

  logical :: grid_statistics = .false.

!----------------------------
! Output file prefix (optional)
!----------------------------
  character(80) :: output_prefix = "hc"

!----------------------------
! End of Default input values
!----------------------------

  namelist / input_parameters /   &
                      debug_mode, &
          target_reynolds_number, &
                   target_y_plus, &
                      igrid_type, &
                              x2, &
                         R_outer, &
           nnodes_cylinder_input, &
                           nr_gs, &
                             nre, &
                      domain_cut, &
     uniform_stretching_to_outer, &
                   spacing_ratio, &
             spacing_ratio_outer, &
             generate_ugrid_file, &
          ugrid_file_unformatted, &
          generate_p3d_ufmt_file, &
                 generate_k_file, &
              k_file_unformatted, &
           generate_line_file_nc, &
           generate_line_file_cc, &
             generate_tec_file_b, &
             generate_tec_file_v, &
             grid_statistics,     &
             output_prefix

 contains

!*****************************************************************************
!* Read input_parameters in the input file: file name = namelist_file
!*****************************************************************************
  subroutine read_nml_input_parameters(namelist_file)

  implicit none
  character(*), intent(in) :: namelist_file
  integer :: os

  write(*,*) "**************************************************************"
  write(*,*) " List of namelist variables and their values"
  write(*,*)

  open(unit=10,file=namelist_file,form='formatted',status='old',iostat=os)
  if (os /= 0) then
    write(*,*) "Warning: Could not open input file ", trim(namelist_file)
    write(*,*) "Will use default values and command line arguments."
  else
    read(unit=10,nml=input_parameters,iostat=os)
    if (os /= 0) then
      write(*,*) "Warning: Error reading namelist from file ", trim(namelist_file)
      write(*,*) "Will use default values and command line arguments."
    else
      write(*,nml=input_parameters) ! Print the namelist variables.
    end if
    close(10)
  end if

  end subroutine read_nml_input_parameters

 end module input_parameter_module
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!  End of input parameter module
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!*******************************************************************************
!*******************************************************************************
! Main program begins here.
!*******************************************************************************
!*******************************************************************************
program hemisphere_cylinder_grid

 use input_parameter_module

 implicit none

  ! Add the interface to the C functions
  interface
    subroutine edu3d_euler_debug_hooks_init() bind(C, name="edu3d_euler_debug_hooks_init")
    end subroutine edu3d_euler_debug_hooks_init
    
    subroutine edu3d_euler_debug_hooks_cleanup() bind(C, name="edu3d_euler_debug_hooks_cleanup")
    end subroutine edu3d_euler_debug_hooks_cleanup
  end interface

! --- CLI argument parsing ---
integer :: n_args, i_arg
character(256) :: arg, value_str
logical :: show_help
logical :: has_namelist_file
character(256) :: namelist_file

! Parameters
  integer , parameter ::  kd = selected_int_kind(8)
  real(dp), parameter ::  zero = 0.0_dp
  real(dp), parameter ::   one = 1.0_dp
  real(dp), parameter ::   two = 2.0_dp
  real(dp), parameter :: three = 3.0_dp
  real(dp), parameter ::  half = 0.5_dp
  real(dp), parameter ::    pi = 3.14159265358979323846_dp

  integer(kd), parameter :: zero_int = 0_kd
  integer(kd), parameter ::  two_int = 2_kd

! Custom data types
  type tria_data
   integer(kd), dimension(3) :: v    !Vertices (nodes) of the triangle
   integer(kd)               :: type !Type of triangle (upward, downward, cylinder)
  end type tria_data

  type node_data_yz
    integer(kd) :: gnode     !Global node number
   real(dp) :: y, z      !Coordinates in yz-plane.
   integer(kd)  :: i         !1D coordinate index from the apex to the base.
  end type node_data_yz

  type node_data
   real(dp) :: x,y,z       !Coordinates in xyz space
   integer(kd)  :: i           !1D coordinate index from the apex to the base.
   real(dp) :: nx,ny,nz    !Unit vector normal to the surface
   real(dp) :: nx2,ny2,nz2 !Unit vector along the second direction
  end type node_data

  type cell_data
   real(dp) :: x,y,z       !Coordinates in xyz space
  end type cell_data

! Output File

  character(80) :: elmtype

! (1) Default output: Files to be used for CFD computations
  character(80) :: filename_mapbc
  character(80) :: filename_ugrid
  character(80) :: filename_p3d
  character(80) :: filename_nmf

! (2)[Optional] Tecplot files for debugging.
  character(80) :: filename01 = "debug01_generating_sector.dat"
  character(80) :: filename02 = "debug02_disk.dat"
  character(80) :: filename03 = "debug03_hemisphere_surface.dat"
  character(80) :: filename07 = "debug07_hemisphere_cylinder_surface.dat"
  character(80) :: filename09 = "debug09_outer_boundary.dat"
  character(80) :: filename11 = "debug11_cc_line_nm.dat"
  character(80) :: filename12 = "debug12_cc_line_all.dat"

! (3)[Optional] Tecplot files for viewing grids.
  character(80) :: filename_tecplot_b
  character(80) :: filename_tecplot_v

! (4)[Optional] Auxiliary Files
  character(80) :: filename_lines      !Line information within BL
  character(80) :: filename_lines_all !Line information all the way to outer boundary
  
  character(80) :: filename_lines_c      !Line information within BL
  character(80) :: filename_lines_c_all !Line information all the way to outer boundary
  
  character(80) :: filename_stats

  character(80) :: filename_k          !Structured indices.

! Local variables
  integer(kd) :: os           !IO constant
  real(dp)    :: x1           !Left end of the cylinder.
  integer(kd) :: i,j,k,inode

  integer(kd) :: ntrias_gs    !Number of triangles in the generating sector
  integer(kd) :: nnodes_gs    !Number of nodes in the generating sector
  real(dp)    :: Rd           !Radius of the hemisphere
  real(dp)    :: r_gs         !Radius of the generating sector
  real(dp)    :: dr_gs        !Spacing along the side of the generating sector
  real(dp)    :: r2, dtheta, theta

  type(node_data_yz), dimension(:),     pointer :: nodes1, nodes2, node_gs
  integer(kd)       , dimension(:),     pointer :: k1_gs, k2_gs, k3_gs
  type(tria_data)   , dimension(:), allocatable :: tria_gs

  integer(kd)                                       :: nnodes_disk, ntrias_disk
  type(node_data_yz)    , dimension(:),     pointer :: node_disk
  integer(kd)           , dimension(:),     pointer :: k1_disk
  integer(kd)           , dimension(:),     pointer :: k2_disk
  integer(kd)           , dimension(:),     pointer :: k3_disk
  type(tria_data)       , dimension(:), allocatable :: tria_disk

  real(dp)     :: xp, yp, zp
  integer(kd)  :: nnodes
  real(dp)     :: s
  
  type(node_data), dimension(:),     pointer :: node_body
  type(node_data), dimension(:),     pointer :: node_outer
  type(tria_data), dimension(:), allocatable :: tria
  
  integer(kd), dimension(:),     pointer :: k1_body, k2_body, k3_body, k4_body
  integer(kd), dimension(:), allocatable :: node_map
  integer(kd)                            :: node1, node2, node3

  integer(kd)                              :: nnodes_circum, nnodes_cylinder , ntrias
  integer(kd), dimension(:)  , allocatable :: nodes_circum
  integer(kd), dimension(:,:), allocatable :: nodes_cylinder
  real(dp) :: dx
  type(cell_data), dimension(:)    ,  pointer    :: cell
  type(node_data), dimension(:)    ,  pointer    :: node
  integer(kd)    , dimension(:)    , allocatable :: k1, k2, k3, k4, k5
  integer(kd)    , dimension(:,:,:), allocatable :: k2n
  integer(kd)    , dimension(:,:,:), allocatable :: k2n_temp

  integer(kd) :: nr, nm
  real(dp)    :: rmax, drN, dr1, gr
  real(dp)    :: sf
  real(dp)    :: dxnmp1

  real(dp), allocatable, dimension(:) :: vspacing

  integer(kd) :: nnodes_body
  integer(kd) :: node4, node5, node6
  integer(kd) :: ntrias_b, nquads_b, ntet, nprs,   nhex,nquads_cyl

  integer(kd),  dimension(:,:), allocatable :: tri, quad, tet, prs, quad_cyl, hex
  integer(kd),  dimension(:)  , allocatable :: node_above
  integer(kd),  dimension(:)  , allocatable :: cell_above
  integer(kd),  dimension(:)  , allocatable :: cell_body

  real(dp)     :: dirx, diry, dirz, xi
  integer(kd)  :: mglevels
  real(dp)     :: cf
  integer(kd)  :: sk12p, sk23p, sk31p
  integer(kd)  :: sk12m, sk23m, sk31m
  integer(kd)  :: n_sections, n_temp

! New interior grid generation
  real(dp)    :: mag, magp, h_outer, h_hemisphere
  real(dp)    :: xx,yy,zz, xis
  integer(kd) :: kk, km

! Strct grid
  real(dp)    :: thetax, dthetax, thetayz, dthetayz
  integer(kd) :: nquads_hem, n1

  integer(kd),  dimension(:,:), allocatable :: quad_hem

  logical       :: iterate_for_sf
  character(80) :: nm_char
  integer(kd)   :: nh0, nc0, nr0, cglevels, nk1, nk4, nk5

  logical, dimension(3) :: negative_volume_detected
  
  real(dp) :: temp_real
  integer :: temp_int
  logical :: temp_logical

  ! Initialize debug hooks
  call edu3d_euler_debug_hooks_init()

  ! Nullify pointers, otherwise my_alloc_ndyz_ptr may attempt
  ! to deallocate a dangling pointer.
  nodes1 => NULL()
  nodes2 => NULL()

show_help = .false.
has_namelist_file = .false.
namelist_file = 'input.nml'

! Process command line arguments
n_args = command_argument_count()
i_arg = 1

do while (i_arg <= n_args)
  call get_command_argument(i_arg, arg)
  select case (trim(arg))
  case ('-h', '--help')
    show_help = .true.
    i_arg = i_arg + 1
  case ('-nml', '--namelist')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, arg)
      has_namelist_file = .true.
      namelist_file = trim(arg)
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('-o', '--output')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, arg)
      output_prefix = trim(arg)
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--debug')
    debug_mode = .true.
    i_arg = i_arg + 1
  case ('--reynolds', '-re')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_real
      target_reynolds_number = temp_real
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--yplus', '-yp')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_real
      target_y_plus = temp_real
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--grid-type', '-gt')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_int
      igrid_type = temp_int
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--length', '-l')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_real
      x2 = temp_real
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--radius', '-r')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_real
      R_outer = temp_real
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--cylinder-nodes', '-cn')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_int
      nnodes_cylinder_input = temp_int
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--tip-elements', '-te')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_int
      nr_gs = temp_int
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--radial-elements', '-radel')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_int
      nre = temp_int
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--domain', '-d')
    if (i_arg + 1 <= n_args) then
      i_arg = i_arg + 1
      call get_command_argument(i_arg, value_str)
      read(value_str, *) temp_int
      domain_cut = temp_int
      i_arg = i_arg + 1
    else
      write(*,*) 'Error: Missing value after ', trim(arg)
      stop
    endif
  case ('--no-ugrid')
    generate_ugrid_file = .false.
    i_arg = i_arg + 1
  case ('--no-p3d')
    generate_p3d_ufmt_file = .false.
    i_arg = i_arg + 1
  case ('--no-k-file')
    generate_k_file = .false.
    i_arg = i_arg + 1
  case ('--no-lines-nc')
    generate_line_file_nc = .false.
    i_arg = i_arg + 1
  case ('--no-lines-cc')
    generate_line_file_cc = .false.
    i_arg = i_arg + 1
  case ('--tecplot-boundary', '-tb')
    generate_tec_file_b = .true.
    i_arg = i_arg + 1
  case ('--tecplot-volume', '-tv')
    generate_tec_file_v = .true.
    i_arg = i_arg + 1
  case ('--statistics', '-s')
    grid_statistics = .true.
    i_arg = i_arg + 1
  case default
    write(*,*) 'Unrecognized command-line option: ', trim(arg)
    write(*,*) 'Use --help for usage information'
    stop
  end select
end do

if (show_help) then
  write(*,*) 'Usage: hcf_hc [options]'
  write(*,*) 'Options:'
  write(*,*) '  -h, --help                       Show this help message'
  write(*,*) '  -nml, --namelist FILE            Specify input namelist file (default: input.nml)'
  write(*,*) '  -o, --output PREFIX              Set output file prefix (default: hc)'
  write(*,*) '  --debug                          Enable debug mode'
  write(*,*) '  --reynolds, -re VALUE            Set Reynolds number (default: 1.0e7)'
  write(*,*) '  --yplus, -yp VALUE               Set target y-plus (default: 10.0)'
  write(*,*) '  --grid-type, -gt TYPE            Set grid type (1-5, default: 1)'
  write(*,*) '                                    1 = Prismatic'
  write(*,*) '                                    2 = Tetrahedral'
  write(*,*) '                                    3 = Mixed - Prism in BL, Tetra outside'
  write(*,*) '                                    4 = Mixed - Prism on hemisphere, Hexa on cylinder'
  write(*,*) '                                    5 = Structured - Hexa with prisms at tip center'
  write(*,*) '  --length, -l VALUE               Set hemisphere-cylinder length (default: 10.0)'
  write(*,*) '  --radius, -r VALUE               Set outer boundary radius (default: 100.0)'
  write(*,*) '  --cylinder-nodes, -cn VALUE      Set number of elements along cylinder (default: 24)'
  write(*,*) '  --tip-elements, -te VALUE        Set number of elements across tip/2 (default: 8)'
  write(*,*) '  --radial-elements, -re VALUE     Set number of elements in radial direction (default: 88)'
  write(*,*) '  --domain, -d VALUE               Set domain type (1=full, 2=half, default: 1)'
  write(*,*) '  --no-ugrid                       Disable generation of .ugrid file'
  write(*,*) '  --no-p3d                         Disable generation of PLOT3D files'
  write(*,*) '  --no-k-file                      Disable generation of K file'
  write(*,*) '  --no-lines-nc                    Disable generation of node line files'
  write(*,*) '  --no-lines-cc                    Disable generation of cell line files'
  write(*,*) '  --tecplot-boundary, -tb          Enable generation of tecplot boundary file'
  write(*,*) '  --tecplot-volume, -tv            Enable generation of tecplot volume file'
  write(*,*) '  --statistics, -s                 Enable grid statistics computation'
  stop
endif

    ntrias     = 0
    nquads_hem = 0
    nquads_cyl = 0
  
          nhex = 0
          nprs = 0
          ntet = 0
      ntrias_b = 0
      nquads_b = 0

   write(*,*)
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) " Hemisphere-cylinder grid generator"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*) "----------------------------------------------------------------"
   write(*,*)

!*******************************************************************************
! Read the input parameters, defined in the file named as 'input.nml'.
!*******************************************************************************

   write(*,*) "Reading the input file: ", trim(namelist_file), "..... "
   write(*,*)
   call read_nml_input_parameters(namelist_file)
   write(*,*)
   
   ! Display CLI overrides
   write(*,*) "Command-line parameter settings:"
   write(*,*) "  Output prefix           = ", trim(output_prefix)
   write(*,*) "  Debug mode              = ", debug_mode
   write(*,*) "  Reynolds number         = ", target_reynolds_number
   write(*,*) "  Target y-plus           = ", target_y_plus
   write(*,*) "  Grid type               = ", igrid_type
   write(*,*) "  HC Length               = ", x2
   write(*,*) "  Outer radius            = ", R_outer
   write(*,*) "  Cylinder nodes          = ", nnodes_cylinder_input
   write(*,*) "  Tip elements            = ", nr_gs
   write(*,*) "  Radial elements         = ", nre
   write(*,*) "  Domain (1=full, 2=half) = ", domain_cut
   write(*,*)

!*******************************************************************************
! Hemisphere Cylinder Geometry (This is an axisymmetric 3D geometry.)
!
!
!      Outer boundary is a hemisphere of radius "R_outer"
!      centered at the base center of the hemisphere-cylinder.
!
!                                  *
!                          *       |
!                  *               |
!             *                    |
!          *                       |
!        *                         |
!      *                           |     z
!    *                             |     ^  y
!   *                              |     | /
!  *         ______________________|     |/
!  *        (______________________|     -------> x
!  *          Hemisphere cylinder  |
!   *                              |
!    *                             |
!      *                           |
!        *                         |
!          *                       |
!             *                    |
!                  *               |
!                          *       |
!                                  *
!
! Note: The inner hemisphere has the unit diameter (radius = 0.5).
!
!
!            <-Hemisphere-> <--- Cylinder --------->
!                      .   . .......................z=0.5
!                .         .                       |         z
!             .            .                       |         ^
!           .              .                       |         |
!           .......................................|         ------> x
!         x=0.0         x=x1=0.5                  x=x2
!
! Note: the length of HC = x2.
!
!*******************************************************************************
!  x0 = zero ! The apex of the hemisphere (leading edge location).
!  y0 = zero
!  z0 = zero

   x1 =  0.5_dp ! Hemisphere-cylinder junction coordinate (=radius of hemisphere)
   Rd =  x1

!*******************************************************************************
! 0. Input parameters
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 0. Set up input parameters"
  write(*,*) "***********************************************************"

!--------------------------------------------------------------------
! Input parameters
!
   if     (igrid_type==1) then

     elmtype = "prism"

   elseif (igrid_type==2) then

     elmtype = "tetra"

   elseif (igrid_type==3) then

     elmtype = "mixed"

   elseif (igrid_type==4) then

     elmtype = "mixed_ph"

   elseif (igrid_type==5) then

     elmtype = "strct"

   else

    write(*,*) " >>> Invalid input value: igrid_type = ", igrid_type
    stop

   endif

! - ugrid_file_unformatted = Grid file format - T: unformatted, F: formatted

   if (ugrid_file_unformatted) then
    if ( big_endian_io(9999) ) then
      write(*,*) 'The system is big Endian'
      write(*,*) ' Ensure big Endian -> setenv F_UFMTENDIAN big'
    else
      write(*,*) 'The system is little Endian'
      write(*,*) ' Ensure little Endian -> setenv F_UFMTENDIAN little'
    endif
   else
      write(*,*) ' Using ascii type for ugrid file'
   endif

! - nre = Number of elements in the radial direction from a boundary node
!         to the corresponding outer boundary node.

   if (nre <= 0) then
    write(*,*) " Number of nodes in the radial direction will be determined in the code."
   endif

!--------------------------------------------------------------------
! Define file names

   filename_mapbc       = trim(output_prefix) // '_' // trim(elmtype) // '.1.mapbc'
   filename_lines       = trim(output_prefix) // '_' // trim(elmtype) // '.1.lines_fmt'
   filename_lines_all   = trim(output_prefix) // '_' // trim(elmtype) // '.1.lines_fmt_all'
   filename_lines_c     = trim(output_prefix) // '_' // trim(elmtype) // '.1.lines_fmt_cc'
   filename_lines_c_all = trim(output_prefix) // '_' // trim(elmtype) // '.1.lines_fmt_cc_all'
   if ( ugrid_file_unformatted ) then
     if ( big_endian_io(9999) ) then
      filename_ugrid  = trim(output_prefix) // '_' // trim(elmtype) // '.1.b8.ugrid'
      filename_p3d    = trim(output_prefix) // '_' // trim(elmtype) // '.1.ufmt'
     else
      filename_ugrid  = trim(output_prefix) // '_' // trim(elmtype) // '.1.l8.ugrid'
      filename_p3d    = trim(output_prefix) // '_' // trim(elmtype) // '.1.ufmt'
     end if
   else
      filename_ugrid  = trim(output_prefix) // '_' // trim(elmtype) // '.1.ugrid'
      filename_p3d    = trim(output_prefix) // '_' // trim(elmtype) // '.1.p3d'
   endif
   filename_nmf       = trim(output_prefix) // '_' // trim(elmtype) // '.1.nmf'
   filename_k         = trim(output_prefix) // '_' // trim(elmtype) // '.1.k'
   filename_tecplot_b = trim(output_prefix) // '_' // trim(elmtype) // '.1.tec_bndary.dat'
   filename_tecplot_v = trim(output_prefix) // '_' // trim(elmtype) // '.1.tec_volume.dat'

   filename_stats     = trim(output_prefix) // '_' // trim(elmtype) // '.1.stats.txt'

!--------------------------------------------------------------------
! Domain cut: Currently 2 choices. Full or half.

    if     (domain_cut == 1) then
     n_sections = 6 ! Full geometry
    elseif (domain_cut == 2) then
     n_sections = 3 ! 180-degree section grid
    else
     write(*,*) " Not valid: domain_cut must be 1 or 2... Stop..."
     stop
    endif
!   n_sections = 1 !  60-degree section grid [Not correctly implemented yet]

!--------------------------------------------------------------------
! Maximum coarse grid level.
!--------------------------------------------------------------------
   write(*,*)
   write(*,*)
   write(*,*) "   -------------------------------------------------  "
   write(*,*) "   ----   Maximum level of regular coarsening   ----  "
   write(*,*) "   -------------------------------------------------  "

   nh0 = nr_gs
   nc0 = nnodes_cylinder_input
   nr0 = nre

   if (igrid_type == 5) then
    nk4 = nh0 + nc0
   else
    nk4 = nc0
   endif

    nk5 = nr0

   if     (domain_cut == 1) then
    nk1 = 6*nh0
   elseif (domain_cut == 2) then
    nk1 = 3*nh0
   endif

     cglevels = 0
   do

   !Check the current level.
    if ( mod(nk1,two_int)==zero_int .and. mod(nk4,two_int)==zero_int .and. mod(nk5,two_int)==zero_int ) then
     cglevels = cglevels + 1 !can be regualrly coarsened.
    else
     exit  !cannot be regularly coarsend. Stop. The current is the coarsest one.
    endif

    !Coarsen and proceed.
     nk1 = nk1 / 2
     nk4 = nk4 / 2
     nk5 = nk5 / 2

   end do

     write(*,*) "         # of coarse grids that can be generated = ", cglevels
     write(*,*) "   Total grid level (including the target level) = ", cglevels+1

   write(*,*) "   -------------------------------------------------  "
   write(*,*) "  "
   write(*,*) "   -------------------------------------------------  "
   write(*,*)
!--------------------------------------------------------------------
! End of Maximum coarse grid level.
!--------------------------------------------------------------------


!*******************************************************************************
!
! Start of hemisphere surface grid generation.
!
!*******************************************************************************

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! Option 1: A regular triangulation of a hemisphere surface (no polar singularity)
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

 hemisphere_surface_grid : if ( igrid_type < 5  ) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 1. Generate a hemisphere surface grid"
  write(*,*) "***********************************************************"

!*******************************************************************************
! Step 1. Systematic triangulation of the generating sector (isotropic).
!    Resolution of this sector will determine other dimensions.
!
! This is a sector with the central angle 60 degrees.
! It is called here the generating sector.
! It is located in the yz-plane.
! We first triangulate this, and use it to build a triangulation of the
! whole circle (the disk).
!
!       ____________
!       \/\/\/\/\/\/    This is an example corresponding to nr_gs = 6.
!        \/\/\/\/\/
!         \/\/\/\/        z ^
!          \/\/\/           |
!           \/\/            |
!            \/             ----> y
!
! NOTE: The number of triangles = 1 + 3 + 5 + ... + 2*nr_gs-1 = nr_gs*nr_gs
!       The number of nodes  = 1 + 2 + 3 + ... + nr_gs+1
!
! NOTE: Important to distinguish two types of triangles.
!       They will define the way the prism is subdivided into 3 tets.
!       The subdivision must be done carefully to match the triangular faces
!       between two adjacent tetrahedra.
!
!*******************************************************************************

   write(*,*) "  --------------------------------------------------------------"
   write(*,*) "   Step 1. Triangulate the generating sector, the building block"
   write(*,*) "  --------------------------------------------------------------"

   k = 0
   do
    k = k + 1
    if (mod(nr_gs,2**(k-1))==0) then
      mglevels = k
    else
      write(*,*)
      write(*,*) "  Maximum multigrid level would be ", mglevels
      write(*,*)
      exit
    endif
   end do

   r_gs = half*pi*Rd          ! Radius of the isotropic triangle
   write(*,*) "  Radius of the isotropic triangle, r_gs = ", r_gs

  dr_gs = r_gs/real(nr_gs,dp) ! Uniform refinement
   write(*,*) "  dr_gs = ", dr_gs

  nnodes_gs = (nr_gs+1)*((nr_gs+1)+1)/2
   write(*,*) "  nnodes_gs = ", nnodes_gs
  ntrias_gs = nr_gs**2
   write(*,*) "  ntrias_gs = ", ntrias_gs
  allocate(node_gs(nnodes_gs))
  allocate(tria_gs(ntrias_gs))
  allocate(k1_gs(nnodes_gs))
  allocate(k2_gs(nnodes_gs))
  allocate(k3_gs(nnodes_gs))

  nnodes_gs = 0
  ntrias_gs = 0

  triangulate : do i = 1, nr_gs

!   r1 = dr_gs * real(i-1,dp) ! Current arc
    r2 = dr_gs * real(i,dp)   !    Next arc

! Nodes on the current arc (r = r1)
   call my_alloc_ndyz_ptr(nodes1,i)

   if (i ==  1) then

      nodes1(1)%y = zero
      nodes1(1)%z = zero
      nodes1(1)%gnode = 1

        nnodes_gs = 1
     node_gs(1)%y = zero
     node_gs(1)%z = zero
     node_gs(1)%i = 1

     k1_gs(1) = 0
     k2_gs(1) = 0
     k3_gs(1) = - ( k1_gs(1) + k2_gs(1) )

   else

    do k = 1, i
     nodes1(k)%y     = nodes2(k)%y
     nodes1(k)%z     = nodes2(k)%z
     nodes1(k)%gnode = nodes2(k)%gnode
    end do

   endif

! Nodes on the next arc (r = r2): New nodes

    call my_alloc_ndyz_ptr(nodes2,i+1)
    dtheta = (pi/three) / real(i,dp)
    do k = 1, i+1
     theta = dtheta * real(k-1,dp)
     nodes2(k)%y = r2 * cos(theta)
     nodes2(k)%z = r2 * sin(theta)

     nnodes_gs = nnodes_gs + 1
     nodes2(k)%gnode = nnodes_gs
     node_gs(nnodes_gs)%y = nodes2(k)%y
     node_gs(nnodes_gs)%z = nodes2(k)%z
     node_gs(nnodes_gs)%i = i+1


     k1_gs(nnodes_gs) = i + (1 - k)
     k2_gs(nnodes_gs) = k - 1
     k3_gs(nnodes_gs) = - ( k1_gs(nnodes_gs) + k2_gs(nnodes_gs) )

    end do

! Triangulate the region between nodes1 and nodes2.
! Right to left

! NOTE: Nodes are ordered clockwise here. It will be counter-clockwise
!       when the hemisphere surface is seen from the interior domain.

! Type 1 triangle
!
! nodes2(k+1)   nodes2(k)
!      1             2
!       o-----------o
!        \         /  
!         \       /
!          \     /
!           \   /
!            \ /
!             o
!             3
!         nodes1(k)

   downward_tria : do k = 1, i
    ntrias_gs = ntrias_gs + 1
    tria_gs(ntrias_gs)%v(1) = nodes2(k+1)%gnode
    tria_gs(ntrias_gs)%v(2) = nodes2(k  )%gnode
    tria_gs(ntrias_gs)%v(3) = nodes1(k  )%gnode
    tria_gs(ntrias_gs)%type = 1
   end do downward_tria


! Type 2 triangle
!
!         nodes2(k+1)
!             3
!             o
!            / \
!           /   \
!          /     \ 
!         /       \
!        /         \
!       o-----------o
!      2             1
! nodes1(k+1)   nodes1(k)

   if (i > 1) then
    upward_tria : do k = 1, i-1
     ntrias_gs = ntrias_gs + 1
     tria_gs(ntrias_gs)%v(1) = nodes1(k  )%gnode
     tria_gs(ntrias_gs)%v(2) = nodes1(k+1)%gnode
     tria_gs(ntrias_gs)%v(3) = nodes2(k+1)%gnode
     tria_gs(ntrias_gs)%type = 2

    end do upward_tria
   endif

  end do triangulate

! We now have a triangulation defined by tria_gs and node_gs.
! Number of triangles = ntrias_gs
! Number of nodes     = nnodes_gs

!*******************************************************************************
! Write a Tecplot file for the generating sector
!******************************************************************************
 debug_mode_01 : if (debug_mode) then
 open(unit=1, file=filename01, status="unknown", iostat=os)

  write(1,*) 'TITLE = "GRID"'
  write(1,*) 'VARIABLES = "x","y","z","k1+k2","i"'
  write(1,*) 'ZONE  N=', nnodes_gs,',E=', ntrias_gs,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes_gs
    write(1,'(3ES20.10,2i10)') 0.0, node_gs(i)%y, node_gs(i)%z, k1_gs(i)+k2_gs(i),node_gs(i)%i
  end do

! Triangles
  do i = 1, ntrias_gs
   write(1,'(3I10)') tria_gs(i)%v(1), tria_gs(i)%v(2), tria_gs(i)%v(3)
  end do

 close(1)

 write(*,*) "  Tecplot file has been written: ", filename01
 endif debug_mode_01
!*******************************************************************************

!*******************************************************************************
! Step 2. Generate a triangulated disk.
!
!
! Rotate and copy the triangulation generated above (0) onto 5 places
! [1,2,3,4,5 below] to form a triangulation of a whole disk (6 patches in total).
!
!           .  . 
!        \        /                  ^ z
!     .   \   1  /   .               |
!    .     \    /     .              |
!   .   2   \  /   0   .             |
!   _________\/_________             ------> y
!   .        /\        .
!   .   3   /  \   5   .
!    .     /    \     .     Copy the grid 0 to 1, 2, 3, 4, 5.
!     .   /      \   .
!        /   4    \ .
!           .  .
!
! Generate new data: tria_disk and node_disk
! Number of triangles = ntrias_disk
! Number of nodes     = nnodes_disk
!*******************************************************************************

  write(*,*) "  --------------------------------------------------------------"
  write(*,*) "   Step 2. Generate a disk by copying the generating sector"
  write(*,*) "  --------------------------------------------------------------"

  nnodes_disk = nnodes_gs * n_sections !<- More than enough (by the # of overlapping nodes)
  ntrias_disk = ntrias_gs * n_sections
  allocate(node_disk(nnodes_disk))
  allocate(tria_disk(ntrias_disk))
  allocate(k1_disk(nnodes_disk))
  allocate(k2_disk(nnodes_disk))
  allocate(k3_disk(nnodes_disk))

  nnodes_disk = 0
  ntrias_disk = 0

! Copy the data from the generating-sector (0) data to the disk data.

!  Copy the node data
!  This is the first one: theta = 0 to 60 degrees

   do i = 1, nnodes_gs

    node_disk(i)%y = node_gs(i)%y
    node_disk(i)%z = node_gs(i)%z
    node_disk(i)%i = node_gs(i)%i

        k1_disk(i) = k1_gs(i)
        k2_disk(i) = k2_gs(i)
        k3_disk(i) = k3_gs(i)

   end do

   nnodes_disk = nnodes_gs

!  Copy the triangle data
!  This is the first one: theta = 0 to 60 degrees

   do i = 1, ntrias_gs
    tria_disk(i)%v    = tria_gs(i)%v
    tria_disk(i)%type = tria_gs(i)%type
   end do

   ntrias_disk = ntrias_gs

! Now generate other parts of the disk: i=1,2,3,4,5
! 1. theta =  60 to 120 degrees
! 2. theta = 120 to 180 degrees
! 3. theta = 180 to 240 degrees
! 4. theta = 240 to 300 degrees
! 5. theta = 300 to 360 degrees

! Each part has (n+1)*(n+2)/2 nodes: 1+2+3+...+(nr_gs+1) = (n_gs+1)*(n_gs+2)/2
  allocate(node_map((nr_gs + 1)*(nr_gs + 2)/2))

 full_geometry : if (n_sections /= 1) then

 new_sectors : do i = 1, n_sections-1

! (1)Generate new nodes and construct a node map (one to one)
!    inode = local node numbering for node_map(:)

  node_map(1) = 1 !Node at the origin

  do k = 2, nr_gs + 1 !Origin to outer boundary
   do j = 1, k        !Right to left

    inode = (k-1)*k/2 + j !Local node number = Right-most node + increment(1,2,..,k)

!   Right side nodes are existing nodes: Left side nodes of the previous sector
    if (j==1) then

     if     (i == 1) then
      node_map(inode) = (k-1)*k/2 + k         !Left-most node of the original sector
     elseif (i == 2) then
      node_map(inode) = nnodes_gs + (k-1)*k/2 !Left-most node of the second sector
     else
      node_map(inode) = nnodes_gs + (nnodes_gs-(nr_gs+1))*(i-2) + (k-1)*k/2
     endif

!   Left side of the last one is the right side of the original sector
    elseif (i==5 .and. j==k) then

      node_map(inode) = (k-1)*k/2 + 1         !Right-most node of the original sector

!   New nodes: Rotate the original nodes by theta = i*pi/3 (i times 60 degrees).
    else

     theta = (pi/three) * real(i,dp)
     nnodes_disk = nnodes_disk + 1
     node_map(inode) = nnodes_disk
     node_disk(nnodes_disk)%y = cos(theta)*node_gs(inode)%y - sin(theta)*node_gs(inode)%z
     node_disk(nnodes_disk)%z = sin(theta)*node_gs(inode)%y + cos(theta)*node_gs(inode)%z
     node_disk(nnodes_disk)%i = k

     if (i==1) then

      k1_disk(nnodes_disk) = -k2_gs(inode)
      k2_disk(nnodes_disk) = -k3_gs(inode)
      k3_disk(nnodes_disk) = -k1_gs(inode)

     elseif(i==2) then

      k1_disk(nnodes_disk) =  k3_gs(inode)
      k2_disk(nnodes_disk) =  k1_gs(inode)
      k3_disk(nnodes_disk) =  k2_gs(inode)
 
     elseif(i==3) then
 
      k1_disk(nnodes_disk) = -k1_gs(inode)
      k2_disk(nnodes_disk) = -k2_gs(inode)
      k3_disk(nnodes_disk) = -k3_gs(inode)

     elseif(i==4) then
 
      k1_disk(nnodes_disk) =  k2_gs(inode)
      k2_disk(nnodes_disk) =  k3_gs(inode)
      k3_disk(nnodes_disk) =  k1_gs(inode)

     elseif(i==5) then
 
      k1_disk(nnodes_disk) = -k3_gs(inode)
      k2_disk(nnodes_disk) = -k1_gs(inode)
      k3_disk(nnodes_disk) = -k2_gs(inode)

     endif

    endif
   end do
  end do

! (2)Generate triangles on the new sector.

  do k = 1, ntrias_gs
   ntrias_disk = ntrias_disk + 1
   tria_disk(ntrias_disk)%v    = node_map(tria_gs(k)%v)
   tria_disk(ntrias_disk)%type = tria_gs(k)%type
  end do

 end do new_sectors

 endif full_geometry

  write(*,*) "  nnodes_disk = ", nnodes_disk
  write(*,*) "  ntrias_disk = ", ntrias_disk

 deallocate(k1_gs,k2_gs,k3_gs)

!--------------------------------------------------------------------------------

! At this point, we have a triangulation of a disk defined by
! tria_disk and node_disk.
! Number of triangles = ntrias_disk
! Number of nodes     = nnodes_disk

!*******************************************************************************
! Write a Tecplot file for the triangulated disk.
!******************************************************************************
 debug_mode_02 : if (debug_mode) then
 open(unit=2, file=filename02, status="unknown", iostat=os)

  write(2,*) 'TITLE = "GRID"'
  write(2,*) 'VARIABLES = "x","y","z","k1","k2","k3","kd","i"'
  write(2,*) 'ZONE  N=', nnodes_disk,',E=', ntrias_disk,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes_disk

   sk12p = ( 1 + int_sign( k1_disk(i)*k2_disk(i) ) )/2
   sk23p = ( 1 + int_sign( k2_disk(i)*k3_disk(i) ) )/2
   sk31p = ( 1 + int_sign( k3_disk(i)*k1_disk(i) ) )/2

   sk12m = ( 1 - int_sign( k1_disk(i)*k2_disk(i) ) )/2
   sk23m = ( 1 - int_sign( k2_disk(i)*k3_disk(i) ) )/2
   sk31m = ( 1 - int_sign( k3_disk(i)*k1_disk(i) ) )/2

   if     (sk12p > 0) then

     sk23p = 0
     sk23m = 1

     sk31p = 0
     sk31m = 1

   elseif (sk23p > 0) then

     sk12p = 0
     sk12m = 1

     sk31p = 0
     sk31m = 1

   elseif (sk31p > 0) then

     sk12p = 0
     sk12m = 1

     sk23p = 0
     sk23m = 1

   endif

    write(2,'(3ES20.10,5i10)') 0.0,  node_disk(i)%y, node_disk(i)%z, k1_disk(i),k2_disk(i),k3_disk(i), &
    sk31m*sk23m*sk12p*abs(k3_disk(i)) + &
    sk31m*sk12m*sk23p*abs(k1_disk(i)) + &
    sk23m*sk12m*sk31p*abs(k2_disk(i)), node_disk(i)%i

  end do

! Triangles
  do i = 1, ntrias_disk
   write(2,'(3I10)') tria_disk(i)%v(1), tria_disk(i)%v(2), tria_disk(i)%v(3)
  end do

 close(2)

 write(*,*) "  Tecplot file has been written: ", filename02
 endif debug_mode_02
!*******************************************************************************

!*******************************************************************************
! 3. Map the disk triangulation onto the hemisphere surface.
!
! Gently attach the disk onto the hemisphere.
! This is done locally at each node, just by using the node_disk data.
! It doesn't matter how they are ordered.
! Connectivity data are unchanged.
!
! Note: In either full or half geometry case, k2=0 indicates that a node
!       is located at z=0 plane:
! 
!           o  o    
!        o        o            z
! k2=0  o__________o           ^
!       o          o           |
!        o        o            |
!           o  o       y<-------
!

!*******************************************************************************

  write(*,*) "  --------------------------------------------------------------"
  write(*,*) "   Step 3. Place the disk onto the hemisphere "
  write(*,*) "  --------------------------------------------------------------"

  s = (half*pi*Rd) / real(nr_gs,dp)
  write(*,*) "  s = ", s
  nnodes_cylinder = int( (x2-x1)/s )   !Isotropic grid
  nnodes_cylinder = int( (x2-x1)/s )/12

! Use input value if requested
  if (nnodes_cylinder_input > 0) then
    nnodes_cylinder = nnodes_cylinder_input
  endif

  write(*,*) "  nnodes_cylinder (actual) = ", nnodes_cylinder


!  nnodes_circum = 6*(nr_gs)
  nnodes_circum = n_sections*nr_gs

  !Add the last node for a non-closed circumferential grid:
  if (n_sections /= 6) nnodes_circum = nnodes_circum + 1

  nnodes = nnodes_disk + nnodes_circum*(nnodes_cylinder+1)
  allocate(node_body(nnodes))
  allocate(k1_body(nnodes))
  allocate(k2_body(nnodes))
  allocate(k3_body(nnodes))
  allocate(k4_body(nnodes))

   nnodes = 0

  do i = 1, nnodes_disk

   k1_body(i) = -k1_disk(i) ! Flip the sign JFC
   k2_body(i) = -k2_disk(i) ! Flip the sign JFC
   k3_body(i) = -k3_disk(i) ! Flip the sign JFC
   k4_body(i) =  0          ! <- On the hemisphere

! Push the node onto the circle located at x=Rd.
   s = sqrt(node_disk(i)%y**2 + node_disk(i)%z**2)
   if (i==1) then
    yp = zero
    zp = zero
   else
    yp = node_disk(i)%y/s * Rd !Extend it to the circle of radius Rd
    zp = node_disk(i)%z/s * Rd !Extend it to the circle of radius Rd
   endif

!    xp = Rd  The circle is located at x = Rd.

! Now, the node (xp,yp,zp) is located on the perimeter of the circle at x=Rd.
! Rotate the node onto the sphere.
            theta = s/Rd
           nnodes = nnodes + 1
   node_body(nnodes)%y = yp*sin(theta)
   node_body(nnodes)%z = zp*sin(theta)
   node_body(nnodes)%x = Rd - Rd*cos(theta)

!  Surface normal direction along which we go up to generate prismatic elements.
   node_body(nnodes)%nx = node_body(nnodes)%x - Rd
   node_body(nnodes)%ny = node_body(nnodes)%y
   node_body(nnodes)%nz = node_body(nnodes)%z

!  Make it the unit vector (well, probably already a unit vector, though...)
   sf = sqrt(node_body(nnodes)%nx**2 + node_body(nnodes)%ny**2 + node_body(nnodes)%nz**2)
   node_body(nnodes)%nx = node_body(nnodes)%nx / sf
   node_body(nnodes)%ny = node_body(nnodes)%ny / sf
   node_body(nnodes)%nz = node_body(nnodes)%nz / sf

   node_body(nnodes)%i = node_disk(i)%i

  end do

 deallocate(k1_disk,k2_disk,k3_disk)
 write(*,*) "  nnodes, nnodes_disk = ", nnodes, nnodes_disk 

!*******************************************************************************
! Write a Tecplot file for the triangulated hemisphere surface.
!******************************************************************************
 debug_mode_03 : if (debug_mode) then
 open(unit=3, file=filename03, status="unknown", iostat=os)

  write(3,*) 'TITLE = "GRID"'
  write(3,*) 'VARIABLES = "x","y","z","k1+k2","i"'
  write(3,*) 'ZONE  N=', nnodes,',E=', ntrias_disk,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes

    write(3,'(3ES20.10,2i10)') node_body(i)%x,  node_body(i)%y, node_body(i)%z, &
                          k1_body(i)+ k2_body(i), node_body(i)%i
  end do

! Triangles
  do i = 1, ntrias_disk
   write(3,'(3I10)') tria_disk(i)%v(1), tria_disk(i)%v(2), tria_disk(i)%v(3)
  end do

 close(3)

 write(*,*) "  Tecplot file has been written: ", filename03
 endif debug_mode_03
!*******************************************************************************


! Construct the list of nodes around the hemisphere at the shoulder.
! nodes_cicum(1:nnodes_circum)
! The list starts at the node located at (y,z)=(0.5,0)
! and goes in the negative y-direction (rotating y-axis to z-axis, pointing
! the positive x-direction).

  allocate(nodes_circum(nnodes_circum+1))
  write(*,*) "  nnodes around the cylinder = ", nnodes_circum

  nnodes_circum = 0

! Circum of generating sector
  do k = 1, nr_gs+1
   nnodes_circum = nnodes_circum + 1
   nodes_circum(nnodes_circum) = nr_gs*(nr_gs+1)/2 + k
  end do

!-------------------------------------------------------------------------------
 full_geom : if (n_sections == 6) then

! Circum of other sectors
  do i = 1, 4
   do k = 1, nr_gs
    nnodes_circum = nnodes_circum + 1
    nodes_circum(nnodes_circum) = nnodes_gs + (i-1)*(nr_gs+1)*nr_gs/2 &
                                            + (nr_gs-1)*nr_gs/2 + k
   end do
  end do

! i == 5
   do k = 1, nr_gs-1
    nnodes_circum = nnodes_circum + 1
    nodes_circum(nnodes_circum) = nnodes_gs + (5-1)*(nr_gs+1)*nr_gs/2 &
                                            + (nr_gs-2)*(nr_gs-1)/2 + k
   end do

  nodes_circum(nnodes_circum + 1) = nodes_circum(1)

 else

! Circum of other sectors
  do i = 1, n_sections-1
   do k = 1, nr_gs
    nnodes_circum = nnodes_circum + 1
    nodes_circum(nnodes_circum) = nnodes_gs + (i-1)*(nr_gs+1)*nr_gs/2 &
                                            + (nr_gs-1)*nr_gs/2 + k
   end do
  end do

 endif full_geom
!-------------------------------------------------------------------------------

! At this point, we have a triangulation of a hemisphere defined by
! tria_disk and node_body.
! Number of triangles = ntrias_disk
! Number of nodes     = nnodes

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! Option 2: Structured polar grid for hemisphere surface.
!           Note: The apex node has many neighbors - polar singularity.
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

 elseif  ( igrid_type == 5  ) then

   write(*,*)
   write(*,*) "***********************************************************"
   write(*,*) " 1. Generate a structured grid over a hemisphere"
   write(*,*) "***********************************************************"

! (1) Full geometry case:
!         
!                          o  o    
!                       o        o            z
!  k1=1 -------------> o          o           ^
!  k1=nnodes_circum -> o          o           |
!                       o        o            |
!                          o  o       y<-------
!         

!(2) Half geometry case:
!                                             z
!                                             ^
!              o  o                           |
!           o        o                        |
! k1=1 --> o__________o               y<-------
!                     ^
!                     |
!            k1=nnodes_circum
!

  !For now, just like the regular triangulation, the number of nodes
  !around the cylinder is set to be

   if (n_sections == 6) then
    nnodes_circum = 6*nr_gs
   else
    nnodes_circum = n_sections*nr_gs + 1
   endif

  !so that the total of 'nnodes_circum' grid-lines run along the body
  !from the apex to the base.

  !Total number of nodes over the hemisphere:

   !         'nr_gs' per grid-line + apex
    nnodes = nr_gs * nnodes_circum +    1

  !Estimate the total number of nodes over the surface.

   !Grid spacing over the hemisphere.

    s = (half*pi*Rd) / real(nr_gs,dp)
    if (debug_mode) write(*,*) " s = ", s

   !Default number of nodes if 'nnodes_cylinder_input' < 0.

    nnodes_cylinder = int( (x2-x1)/s )   !Isotropic grid
    nnodes_cylinder = int( (x2-x1)/s )/12

   !Use input value if requested

    if (nnodes_cylinder_input > 0) then
     nnodes_cylinder = nnodes_cylinder_input
    endif

   !Total number of nodes over the hemisphere-cylinder surface.

    write(*,*) " nnodes_cylinder (actual) = ", nnodes_cylinder
    nnodes = nnodes + nnodes_circum*(nnodes_cylinder+1)


   ! Allocate arrays
    allocate(node_body(nnodes))
    allocate(k1_body(nnodes))
    allocate(k2_body(nnodes))
    allocate(k3_body(nnodes))
    allocate(k4_body(nnodes))

    ntrias_disk = nnodes_circum
    allocate(tria_disk(ntrias_disk))

    nquads_hem = nr_gs * nnodes_circum
    allocate(quad_hem(nquads_hem,5))

    allocate(nodes_circum(nnodes_circum+1))

  !------------------------------------------------------------------
  ! Generate nodes and elements over the hemisphere surface.

    ntrias_disk = 0
    nquads_hem  = 0

   ! Apex node
     nnodes = 1
     node_body(nnodes)%x = zero
     node_body(nnodes)%y = zero
     node_body(nnodes)%z = zero
     node_body(nnodes)%i = 1  !Index along the HC (apex to base): i=1,2,3,...
         k1_body(nnodes) =  0 !Around the HC
         k2_body(nnodes) =  0 !Not used
         k3_body(nnodes) =  0 !Not used (used to indicate the half geometry).
         k4_body(nnodes) =  0 !Along the HC (apex to base)      : i=0,1,2,3,...
    node_body(nnodes)%nx = -one
    node_body(nnodes)%ny = zero
    node_body(nnodes)%nz = zero

    dthetax  = half*pi / real(        nr_gs,dp)
    
   if     (n_sections == 6) then
    dthetayz =  two*pi / real(nnodes_circum,dp)
   elseif (n_sections == 3) then
    dthetayz =      pi / real(nnodes_circum-1,dp)
   else
     write(*,*) " Not valid: domain_cut must be 1 or 2... Stop..."
     stop
   endif

   ! All other nodes
 
    !Looping from the apex to the shoulder.
    do i = 1, nr_gs

      thetax  = dthetax*real(i,dp)  !Angle from the -ve x-axis towards z-axis

     !Looping from a node at (0,0.5,0.0) counterclockwise in (y,z)-plane.
      n1 = nnodes + 1
     do k = 1, nnodes_circum

      thetayz = dthetayz*real(k-1,dp) !Angle from y-axis towards z-axis
       nnodes = nnodes + 1

       !Store the node list at the last location in 'nodes_circum'
        if (i==nr_gs) nodes_circum(k) = nnodes

       !Save the first node
        if (k==1) n1 = nnodes

     !Generate a node.
       node_body(nnodes)%x =  Rd-Rd*cos(thetax) !x=0 at thetax=0; x=Rd at thetax=pi/2
       node_body(nnodes)%y =     Rd*sin(thetax) * cos(thetayz)
       node_body(nnodes)%z =     Rd*sin(thetax) * sin(thetayz)
       node_body(nnodes)%i =  i+1 !Index along the HC (apex to base): i=1,2,3,...
           k1_body(nnodes) =  k   !Index around the HC: k1=1,2,3,...
           k2_body(nnodes) =  0   !Not used
           k3_body(nnodes) =  0   !Not used (used to indicate the half geometry).
           k4_body(nnodes) =  i   !Index along the HC (apex to base): i=0,1,2,3,...
      node_body(nnodes)%nx = node_body(nnodes)%x - Rd
      node_body(nnodes)%ny = node_body(nnodes)%y
      node_body(nnodes)%nz = node_body(nnodes)%z
      sf = sqrt(node_body(nnodes)%nx**2 + node_body(nnodes)%ny**2 + node_body(nnodes)%nz**2)
      node_body(nnodes)%nx = node_body(nnodes)%nx / sf
      node_body(nnodes)%ny = node_body(nnodes)%ny / sf
      node_body(nnodes)%nz = node_body(nnodes)%nz / sf

     !Construct tria elements around the apex node.
      tria_formed : if (i==1) then

       !Triangle vertices ordered pointing inward (to the interior domain).
       if (k > 1) then
 
                ntrias_disk = ntrias_disk + 1
          tria_disk(ntrias_disk)%v(1) = 1        !apex node
          tria_disk(ntrias_disk)%v(2) = nnodes   !right
          tria_disk(ntrias_disk)%v(3) = nnodes-1 !left

         if (k == nnodes_circum .and. n_sections == 6) then
                ntrias_disk = ntrias_disk + 1
          tria_disk(ntrias_disk)%v(1) = 1        !apex node
          tria_disk(ntrias_disk)%v(2) = n1       !right
          tria_disk(ntrias_disk)%v(3) = nnodes   !left
         endif

       endif

      endif tria_formed

     !Construct quad elements
      quad_formed : if (i > 1) then

       !Quad vertices ordered pointing inward (to the interior domain).
       if (k > 1) then
 
                       nquads_hem = nquads_hem + 1
           quad_hem(nquads_hem,1) = nnodes
           quad_hem(nquads_hem,2) = nnodes-1
           quad_hem(nquads_hem,3) = (nnodes-nnodes_circum)-1
           quad_hem(nquads_hem,4) = (nnodes-nnodes_circum)
           quad_hem(nquads_hem,5) = 1

          if (k == nnodes_circum .and. n_sections == 6) then
                       nquads_hem = nquads_hem + 1
           quad_hem(nquads_hem,1) = n1
           quad_hem(nquads_hem,2) = nnodes
           quad_hem(nquads_hem,3) = (nnodes-nnodes_circum)
           quad_hem(nquads_hem,4) = (    n1-nnodes_circum)
           quad_hem(nquads_hem,5) = 1
          endif

       endif

      endif quad_formed

     end do
    end do

    !Close the list by adding the first node in the full geometry case:
     if (n_sections == 6) nodes_circum(nnodes_circum + 1) = nodes_circum(1)

!*******************************************************************************
! Write a Tecplot file for the hemisphere surface.
!******************************************************************************
 debug_mode_03s : if (debug_mode) then
 open(unit=3, file=filename03, status="unknown", iostat=os)

  write(3,*) 'TITLE = "GRID"'
  write(3,*) 'VARIABLES = "x","y","z","k1","k2","k3","k4","i"'
  write(3,*) 'ZONE  N=', nnodes,',E=', ntrias_disk+nquads_hem, &
             ' ,ET=quadrilateral, F=FEPOINT'

! Nodes
  do i = 1, nnodes

    write(3,'(3ES20.10,5i10)') node_body(i)%x,  node_body(i)%y, node_body(i)%z, &
                  k1_body(i),k2_body(i), k3_body(i),k4_body(i), node_body(i)%i
  end do

! Triangles
  do i = 1, ntrias_disk
     write(3,'(4I10)') tria_disk(i)%v(1), tria_disk(i)%v(2), tria_disk(i)%v(3), &
                                                             tria_disk(i)%v(3)
  end do

! Quads
  do i = 1, nquads_hem
   write(3,'(4I10)') quad_hem(i,1),quad_hem(i,2),quad_hem(i,3),quad_hem(i,4)
  end do

 close(3)

 write(*,*) "Tecplot file has been written: ", filename03
 endif debug_mode_03s
!*******************************************************************************


 endif hemisphere_surface_grid


 ! Finally, for the sake of smoothness in the final grid.
 ! Adjust the surface normal at the shoulder nodes (deflect a bit toward negative x).

    dthetax = half * (half*pi / real(nr_gs,dp))
   do k = 1, nnodes_circum
    node_body(nodes_circum(k))%nx = -Rd*sin(dthetax)
    sf = sqrt(node_body(nnodes)%nx**2 + node_body(nnodes)%ny**2 + node_body(nnodes)%nz**2)
    node_body(nnodes)%nx = node_body(nnodes)%nx / sf
    node_body(nnodes)%ny = node_body(nnodes)%ny / sf
    node_body(nnodes)%nz = node_body(nnodes)%nz / sf
   end do
 
!*******************************************************************************
!
! End of hemisphere surface grid generation.
!
!*******************************************************************************

!*******************************************************************************
!
! Start of cylinder surface grid generation.
!
!*******************************************************************************

!*******************************************************************************
! 2. Generate a cylinder surface grid.
!
!    Nodes are generated based on the list of nodes around the base circle of
!    the hemisphere part: nodes_circum(1:nnodes_circum+1).
!    This list needs to be available at this point for whatever the type of
!    the hemisphere surface grid.
!
!
! NOTE: Stretching is applied to the distribution of nodes along the cylinder
!       for a smooth the spacing variation from the shoulder of the hemisphere.
!
! NOTE: Important to distinguish two types of triangles.
!       They will define the way the prism is subdivided into 3 tets.
!       The subdivision must be done carefully to match the triangular faces
!       between two adjacent tetrahedra.
!
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 2. Generate a cylinder surface grid."
  write(*,*) "***********************************************************"

! OK, generate a surface grid over the cylinder part.

  write(*,*) " nnodes along the cylinder = ", nnodes_cylinder

! Compute the stretching factor for a smooth stretching along the cylinder.
! The stretching is activated only if we don't have enough nodes to generate
! an isotropic grid. The default is the isotropic grid.

!  Length of the edge on the hemisphere surface
    s = (half*pi*Rd) / real(nr_gs,dp)

!  User specified spacing for a smooth transition (default=1.1, 10% increase)
    s = s*spacing_ratio

! Stretching is applied if nnodes_cylinder is less than that
! for a uniform spacing.
  if (nnodes_cylinder < int( (x2-x1)/s )) then

    dx = (x2-x1)/real(nnodes_cylinder,dp)
    sf = 2.0_dp

   do k = 1, 500

    dxnmp1 = (one-exp(sf*dx/(x2-x1)))/(one-exp(sf)) * (x2-x1)

    if (abs(dxnmp1 - s )/s < 1.0e-03_dp) exit

    if (dxnmp1 > s ) then
      sf = sf + 0.01_dp
    else
      sf = sf - 0.01_dp
    endif

   end do

   write(*,*) "Determined stretching factor = ",sf, " at iterations = ", k
   write(*,*) "          First spacing on the cylinder = ", dxnmp1
   write(*,*) "  Edge length on the hemisphere surface = ", s/1.18_dp

  endif

! Generate nodes on the cylinder surface

!-------------------------------------------------------------------------
! - Full geometry

 full_geom2 : if (n_sections == 6) then

  allocate(nodes_cylinder(nnodes_circum+1,nnodes_cylinder+1))

! Copy the data
  do k = 1, nnodes_circum+1
   nodes_cylinder(k,1) = nodes_circum(k)
  end do

! Uniform spacing
  dx = (x2-x1)/real(nnodes_cylinder,dp)

! Move along the cylinder in the positive x-direction
  do i = 2, nnodes_cylinder+1

!  Go around the circumferential direction, counterclockwise in (y,z)-plane.
   do k = 1, nnodes_circum

    nnodes = nnodes + 1
    nodes_cylinder(k,i) = nnodes
    node_body(nnodes)%y = node_body(nodes_circum(k))%y
    node_body(nnodes)%z = node_body(nodes_circum(k))%z

    node_body(nnodes)%i = (nr_gs+1) +(i-1)!Index along the HC (apex to base): i=1,2,3,...

    if (nnodes_cylinder /= int( (x2-x1)/s )) then
!    Stretched
     node_body(nnodes)%x = x1 + (one-exp(sf*real(i-1,dp)*dx/(x2-x1)))/(one-exp(sf))*(x2-x1)
    else
!    Uniform spacing
    node_body(nnodes)%x  = node_body(nodes_circum(k))%x + real(i-1,dp)*dx
    endif

!   Unit vector normal to the cylinder surface.
    node_body(nnodes)%nx = zero
    node_body(nnodes)%ny = node_body(nodes_circum(k))%ny
    node_body(nnodes)%nz = node_body(nodes_circum(k))%nz

    if  ( igrid_type < 5  ) then
     k1_body(nnodes) = k1_body(nodes_circum(k))
     k2_body(nnodes) = k2_body(nodes_circum(k))
     k3_body(nnodes) = k3_body(nodes_circum(k))
     !Index along the HC: i=0,1,2,3,... (0 on the hemisphere)
     k4_body(nnodes) = i-1
    elseif ( igrid_type == 5  ) then
     !Index around the HC: counterclockwise in (y,z)-plane.
     k1_body(nnodes) = k1_body(nodes_circum(k)) 
     k2_body(nnodes) = 0            !Not used
     k3_body(nnodes) = 0            !Not used (used to indicate the half geometry).
     k4_body(nnodes) = nr_gs +(i-1) !Index along the HC (apex to base): i=0,1,2,3,...
    endif

    if (k==1)  nodes_cylinder(nnodes_circum+1,i) = nnodes

   end do

  end do

!-------------------------------------------------------------------------
! - Partial geometry

 else !if (n_sections /= 6) then

  allocate(nodes_cylinder(nnodes_circum,nnodes_cylinder+1))

! Copy the data
  do k = 1, nnodes_circum
   nodes_cylinder(k,1) = nodes_circum(k)
  end do

! Uniform spacing
  dx = (x2-x1)/real(nnodes_cylinder,dp)

! Move along the cylinder in the positive x-direction
  do i = 2, nnodes_cylinder+1

!  Go around the circumferential direction
   do k = 1, nnodes_circum

    nnodes = nnodes + 1
    nodes_cylinder(k,i) = nnodes
    node_body(nnodes)%y = node_body(nodes_circum(k))%y
    node_body(nnodes)%z = node_body(nodes_circum(k))%z

    node_body(nnodes)%i = (nr_gs+1) +(i-1)!Index along the HC (apex to base): i=1,2,3,...

    if (nnodes_cylinder /= int( (x2-x1)/s )) then
!    Stretched
     node_body(nnodes)%x = x1 + (one-exp(sf*real(i-1,dp)*dx/(x2-x1)))/(one-exp(sf))*(x2-x1)
    else
!    Uniform spacing
    node_body(nnodes)%x  = node_body(nodes_circum(k))%x + real(i-1,dp)*dx
    endif

!   Unit vector normal to the cylinder surface.
    node_body(nnodes)%nx = zero
    node_body(nnodes)%ny = node_body(nodes_circum(k))%ny
    node_body(nnodes)%nz = node_body(nodes_circum(k))%nz

    if  ( igrid_type < 5  ) then
     k1_body(nnodes) = k1_body(nodes_circum(k))
     k2_body(nnodes) = k2_body(nodes_circum(k))
     k3_body(nnodes) = k3_body(nodes_circum(k))
    !Index along the HC: i=0,1,2,3,... (0 on the hemisphere)
     k4_body(nnodes) = i-1
    elseif ( igrid_type == 5  ) then
    !Index around the HC: counterclockwise in (y,z)-plane.
     k1_body(nnodes) = k1_body(nodes_circum(k))
     k2_body(nnodes) = 0            !Not used
     k3_body(nnodes) = 0            !Not used (used to indicate the half geometry).
     k4_body(nnodes) = nr_gs +(i-1) !Index along the HC (apex to base): i=0,1,2,3,...
    endif

   end do

  end do

 endif full_geom2

!-------------------------------------------------------------------------

 write(*,*) " nnodes_circum = ", nnodes_circum

!*************************************************************************************
!*************************************************************************************
! Quadrilateral grid over the cylinder.
!*************************************************************************************
  mixed_ph_cylinder : if ( igrid_type == 4 .or. igrid_type == 5 ) then

! Allocate and copy the hemisphere triangulation data to the global array

  ntrias = ntrias_disk ! + nnodes_circum*nnodes_cylinder*2
  allocate(tria(ntrias))

 !-----------------------------------------------------------------
 ! Mixed - Prism/Hex

  if     ( igrid_type == 4 ) then

   nquads_cyl = nnodes_circum * nnodes_cylinder
   allocate(quad_cyl(nquads_cyl,5))

  !Initialization.
   nquads_cyl = 0

 !-----------------------------------------------------------------
 ! Structured (Prism + some prisms)

  elseif ( igrid_type == 5 ) then

   nquads_cyl = nnodes_circum * nnodes_cylinder  + nquads_hem
   allocate(quad_cyl(nquads_cyl,5))

  ! Copy the data of quads on the hemisphere surface.
   do i = 1, nquads_hem
    quad_cyl(i,1) = quad_hem(i,1)
    quad_cyl(i,2) = quad_hem(i,2)
    quad_cyl(i,3) = quad_hem(i,3)
    quad_cyl(i,4) = quad_hem(i,4)
    quad_cyl(i,5) = quad_hem(i,5)
   end do

  !Initialization.
   nquads_cyl = nquads_hem

  endif

! Copy the data for triangles on the hemisphere surface.
  do i = 1, ntrias
   tria(i)%v    = tria_disk(i)%v
   tria(i)%type = tria_disk(i)%type
  end do

! Quadrangulate(?) the cylinder surface

  if     (n_sections==6) then
   n_temp = nnodes_circum
  else !if (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

  do i = 1, nnodes_cylinder
   do k = 1, n_temp

!
!   (k,i)        (k,i+1)
!      1           4
!       o----------o          ------> x
!       |          |          |
!       |          |          |
!       |          |          |
!       |          |          v
!       |          |   Circumferential direction
!       o----------o
!       2          3
!   (k+1,i)    (k+1,i+1)

      node1 = nodes_cylinder(k  ,i  )
      node2 = nodes_cylinder(k+1,i  )
      node3 = nodes_cylinder(k+1,i+1)
      node4 = nodes_cylinder(k  ,i+1)
     nquads_cyl = nquads_cyl + 1
     quad_cyl(nquads_cyl,1) = node1
     quad_cyl(nquads_cyl,2) = node2
     quad_cyl(nquads_cyl,3) = node3
     quad_cyl(nquads_cyl,4) = node4
     quad_cyl(nquads_cyl,5) = 1

   end do
  end do

!*************************************************************************************
! Triangular grid over the cylinder.
!*************************************************************************************
  else

! Allocate and copy the hemisphere triangulation data to the global array

  nquads_cyl = 0
  allocate(quad_cyl(1,5))

  ntrias = ntrias_disk + nnodes_circum*nnodes_cylinder*2
  allocate(tria(ntrias))

! Copy the data
  do i = 1, ntrias_disk
   tria(i)%v    = tria_disk(i)%v
   tria(i)%type = tria_disk(i)%type
  end do

! Triangulate the cylinder surface

  ntrias = ntrias_disk

  if     (n_sections==6) then
   n_temp = nnodes_circum
  else !if (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

  do i = 1, nnodes_cylinder
   do k = 1, n_temp

! Type 3 triangle
!
!   (k,i)        (k,i+1)
!      1            3
!       o----------o          ------> x
!       |        .            |
!       |      .              |
!       |    .                |
!       |  .                  v
!       |.          Circumferential direction
!       o
!       2
!   (k+1,i)

      node1 = nodes_cylinder(k  ,i  )
      node2 = nodes_cylinder(k+1,i  )
      node3 = nodes_cylinder(k  ,i+1)
     ntrias = ntrias + 1
     tria(ntrias)%v(1) = node1
     tria(ntrias)%v(2) = node2
     tria(ntrias)%v(3) = node3
     tria(ntrias)%type = 3

! Type 4 triangle
!
!              (k,i+1)
!                  2
!                  o              ------> x
!                . |              |
!             .    |              |
!           .      |              |
!         .        |              v
!       .          |    Circumferential direction
!      o---------- o
!      3           1
!    (k+1,i)    (k+1,i+1)

      node1 = nodes_cylinder(k+1,i+1)
      node2 = nodes_cylinder(k  ,i+1)
      node3 = nodes_cylinder(k+1,i  )
     ntrias = ntrias + 1
     tria(ntrias)%v(1) = node1
     tria(ntrias)%v(2) = node2
     tria(ntrias)%v(3) = node3
     tria(ntrias)%type = 4

   end do
  end do

  endif mixed_ph_cylinder
!*************************************************************************************
!*************************************************************************************
!*************************************************************************************

!*******************************************************************************
! Write a Tecplot file for the hemisphere cylinder surface triangulation.
!******************************************************************************
 debug_mode_04 : if (debug_mode) then
 open(unit=4, file=filename07, status="unknown", iostat=os)

  write(4,*) 'TITLE = "GRID"'
  write(4,*) 'VARIABLES = "x","y","z","i","k1","k2","k3","k4"'
  write(4,*) 'ZONE  N=', nnodes,',E=', ntrias+nquads_cyl,' , ET=quadrilateral, F=FEPOINT'

! Nodes
  do i = 1, nnodes
    write(4,'(3ES20.10,5i10)') node_body(i)%x, node_body(i)%y, node_body(i)%z, node_body(i)%i, &
                              k1_body(i),k2_body(i), k3_body(i),k4_body(i)
  end do

! Triangles
  do i = 1, ntrias
   write(4,'(4I10)') tria(i)%v(1), tria(i)%v(2), tria(i)%v(3), tria(i)%v(3)
  end do

! Quads
  do i = 1, nquads_cyl
   write(4,'(4I10)') quad_cyl(i,1),quad_cyl(i,2),quad_cyl(i,3),quad_cyl(i,4)
  end do


 close(4)

 write(*,*) "Tecplot file has been written: ", filename07

 endif debug_mode_04
!*******************************************************************************

! At this point, we have a surface grid of a hemisphere-cylinder defined by
! List of triangles   = tria
! List of quads       = quad_cyl
! List of nodes       = node_body.
! Number of triangles = ntrias
! Number of quads     = nquads_cyl
! Number of nodes     = nnodes

 write(*,*)
 write(*,*) "------------------------------------------------------------------"
 write(*,*) " At this point, we have a surface grid of a hemisphere-cylinder:"
 write(*,*) "           Nodes = ", nnodes
 write(*,*) "       Triangles = ", ntrias
 write(*,*) "  Quadrilaterals = ", nquads_cyl + nquads_hem
 write(*,*) "------------------------------------------------------------------"
 write(*,*)

!*******************************************************************************
!
! End of cylinder surface grid generation.
!
!*******************************************************************************

!*******************************************************************************
!
! Start of outer boundary grid generation
!
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 3. Generate a hemisphere-shape outer boundary grid."
  write(*,*) "***********************************************************"

  write(*,*) " Map the surface grid to the outer hemisphere."
  write(*,*) " So, the number of nodes are the same as that of the surface grid."

!--------------------------------------------------------------------------------
! Generate nodes on the outer boundary, which is a large hemisphere of
! radius = 'R_outer'.

  nnodes_body = nnodes
  allocate(node_outer(nnodes_body))

  !Uniform node distribution along the hemisphere-cylinder surface,
  !by using the data, node_body(i)%i, which runs from 1 to nr_gs+nnodes_cylinder,
  !one-dimensionally along the hemisphere-cylinder surface from the apex
  !to the base.

   dtheta  = half*pi/real(nr_gs+nnodes_cylinder,dp)
   h_outer = R_outer*sin(dtheta) ! Save the uniform arc-length spacing.

   write(*,*) " Radius of the outer hemisphere   = ", R_outer
   write(*,*) "                  Number of nodes = ", nnodes_body
   write(*,*) "                          Spacing = ", h_outer

  do i = 1, nnodes_body

   theta = real(node_body(i)%i-1,dp)*dtheta

   node_outer(i)%x = x2 - R_outer*cos(theta)
               mag = sqrt( node_body(i)%ny**2 + node_body(i)%nz**2 )
   node_outer(i)%y = R_outer*sin(theta)*node_body(i)%ny/(mag+1.0e-16)
   node_outer(i)%z = R_outer*sin(theta)*node_body(i)%nz/(mag+1.0e-16)
   node_outer(i)%i = node_body(i)%i

  end do

!*******************************************************************************
! Write a Tecplot file for the outer boundary.
!******************************************************************************
 debug_mode_05 : if (debug_mode) then
 open(unit=5, file=filename09, status="unknown", iostat=os)

  write(5,*) 'TITLE = "GRID"'
  write(5,*) 'VARIABLES = "x","y","z","i","k1","k2","k3","k4"'
  write(5,*) 'ZONE  N=', nnodes_body,',E=', ntrias,' , ET=triangle, F=FEPOINT'

! Nodes
  do i = 1, nnodes_body
     write(5,'(3ES20.10,5i10)') node_outer(i)%x, node_outer(i)%y, &
                                node_outer(i)%z, node_outer(i)%i, &
                     k1_body(i),k2_body(i), k3_body(i),k4_body(i)
  end do

! Triangles
  do i = 1, ntrias
   write(5,'(3I10)') tria(i)%v(1), tria(i)%v(2), tria(i)%v(3)
  end do

 close(5)

 write(*,*) "Tecplot file has been written: ", filename09
 endif debug_mode_05
!*******************************************************************************

!*******************************************************************************
!
! End of outer boundary grid generation
!
!*******************************************************************************

!*******************************************************************************
!
! Start of interior grid generation.
!
!*******************************************************************************

!*******************************************************************************
! 4. We now go up to generate the interior nodes.
!
! Nodes are determined separately in the BL region and outer region.
!
! In BL region, determine the vertical spacing, vspacing: r(i+1) = vspacing(i) + r(i),
! and generate interior nodes by going up in the direction of surface normal.
! Outside BL region, exponential stretching is used.
!
!
!   Geometric sequence       Exponential stretching
! o-o--o---o----o-----o-----x-------x----------x--------------x        -----> r
! 1 2  3              nm                                 nr
!
!   BL region                     Outer region         Outer boundary
!
!
! (1)BL region: First spacing off the wall (dr1) is determined for given
!               target Reynolds number and y-plus values.
!               Last spacing in BL is determined based on dr1 or the surface
!               grid spacing.
!
! (2)Outside  : Nodes are generated by following a blended direction.
!               First, generate nodes along the surface normal, and
!               then gradually switch to a vector directing towards
!               the corresponding outer boundary node.
!               Stretching factor is determined to generate isotropic cells
!               at the farfield.
!
! Note: The number of nodes in each region will be automatically determined.
!
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 4. Go up to the outer boundary and generate interior nodes"
  write(*,*) "***********************************************************"

!------------------------------------------------------------------------
! Some global parameters
!------------------------------------------------------------------------

   write(*,*) " Parameters for interior node generation"

! First vertical spacing off the wall:

     cf = 0.026_dp/target_reynolds_number**(1.0_dp/7.0_dp)
    dr1 = ( sqrt(2.0_dp/cf)/target_reynolds_number) * target_y_plus
    write(*,*)
    write(*,*) ">>> Compute the first off-the-wall spacing dr1:"
    write(*,*) "    dr1 determined by the target Re and y_plus:"
    write(*,*) "    target_reynolds_number = ", target_reynolds_number
    write(*,*) "             target_y_plus = ", target_y_plus
    write(*,*) "    --------->         dr1 = ", dr1

! rmax = thickness of the specified boundary-layer region:

    write(*,*)
    write(*,*) ">>> Boundary-layer thickness estimate, rmax, "
    write(*,*) "    from a turbulent flow over a flat plate:"
    write(*,*) "    Eq.(6-72), page 430, 'Viscous Fluid Flow' 2nd edition (1991) by White."
    write(*,*) "    (Use 50% of the estimate.)"
    rmax = bl_height_factor * ( 0.37_dp*x2/target_reynolds_number**(0.2_dp) )
    write(*,*) "    ------------> rmax = ", rmax
    write(*,*)

    if (rmax < dr1) then
     rmax = dr1 * 10.0_dp
     write(*,*)
     write(*,*) ">>> The estimated rmax is too small (rmax < dr1)."
     write(*,*) "    Set rmax = dr1 * 10.0_dp."
     write(*,*) "    ------------> rmax = ", rmax
    endif

    write(*,*)
    write(*,*) "  NOTE: rmax defines the BL region in the grid."
    write(*,*)

! Spacing over the hemisphere

   h_hemisphere = half*pi*Rd / real(nr_gs,dp)
    write(*,*)
    write(*,*) ">>> Typical spacing over the hemisphere: "
    write(*,*) "     h_hemisphere = ", h_hemisphere

!------------------------------------------------------------------------
! Compute the parameters for geometric sequence.
! Determine the number of nodes in r-direction within BL for a given dr1:
!------------------------------------------------------------------------

  !Spacing of the outer most cell in the layer.
   drN = max( 2.0_dp*rmax/(rmax/dr1), 0.1_dp*h_hemisphere)
    write(*,*)
    write(*,*) ">>> Set a spacing of the outer most cell in the BL region: drN"
    write(*,*) "    ------------> drN = ", drN

    if (drN > rmax) then
     write(*,*)
     write(*,*) ">>> drN > rmax !! Impossible... drN = ", drN
     drN = 0.1_dp*h_hemisphere
     write(*,*) "    ------------> Modified drN = ", drN
    endif

     write(*,*)
     write(*,*) ">>> AR = ", h_hemisphere/dr1

   gr = (rmax-dr1)/(rmax-drN)
   nm = ceiling( log(gr*drN/dr1)/log(gr) )
   gr = (drN/dr1)**( one/(real(nm,dp)-one) )

   write(*,*)
   write(*,*) " NOTE: dr1 is a non-dimensionalized length. "
   write(*,*) "       Reference length is the diameter of the hemisphere (=1.0 in the grid)."
   write(*,*)

   write(*,*) "                       Hemisphere spacing = ", h_hemisphere
   write(*,*) "          rmax (the end of the BR region) = ", rmax
   write(*,*) "                                      dr1 = ", dr1
   write(*,*) "                                      drN = ", drN
   write(*,*) "                                       gr = ", gr
   write(*,*) "                           dr1*gr**(nm-1) = ", dr1*gr**(nm-1)

!  At this point, nm is the number of elements within the layer.
!  Add 1 to make it the number of nodes.
   nm = nm + 1 ! where change the type of cell from prism to tetra for a mixed grid.
   write(*,*)
   write(*,*) ">>> Number of elements in the BL region:"
   write(*,*) "    ---> nm (lements) = ", nm

!  If the requested value of nre is less than nm, then we have a problem.
   if (nre > 0) then
    if (nm > nre) then

     write( nm_char  , '(i0)' ) nm
     write(*,*)
     write(*,*) "*************************************************************************"
     write(*,*) "*************************************************************************"
     write(*,*) " Problem! Code has stopped....."
     write(*,*)

     write(*,*) " Input value: Elements in the radial direction = ", nre, " is not enough..."
     write(*,*)
     write(*,*) " [The number of elements in BL region has been found to be larger than the"
     write(*,*) "  total number of elements in the radial direction requested.]"
     write(*,*)

     write(*,*) " Try one of the followings:"
     write(*,*)
     write(*,*) "  (1)Increase the input value: Elements in the radial direction (>", trim(nm_char),")."
     write(*,*) "     Note: (the input value) - ",trim(nm_char), " is the # of elements outside the BL region."
     write(*,*)
     write(*,*) "  (2)Let the code determines the number of elements in the radial direction by"
     write(*,*) "     giving a negative input value: e.g., "
     write(*,*)
     write(*,*) "               Elements in the radial direction = ", -1
     write(*,*)
     write(*,*) "     Note: # of elements is determined by surface spacing of the outer boundary."
     write(*,*)
     write(*,*)
     write(*,*) "  (3)Increase Target y-plus value to reduce the number of elements in BL region."
     write(*,*)
     write(*,*)
     write(*,*) "  (4)Decrease Target Reynolds number to reduce the number of elements in BL region."
     write(*,*)
     write(*,*)
     write(*,*) "*************************************************************************"
     write(*,*) "*************************************************************************"
     write(*,*)
     write(*,*) " No grids generated... Try again."
     write(*,*)

     stop

    endif
   endif

!  Allocate the vertical spacing array
!  Geometric sequence for the cell spacing inside the layer.

   allocate(vspacing(nm-1))

   do i = 1, nm-1
    if (i == 1) then
     vspacing(i) = dr1
    else
     vspacing(i) = dr1*gr**(i-1)
    endif
   end do

! NOTE: The last node in the prismatic layer is nm-th node.
!       The spacing array vspacing() goes from 1 to nm-1.

!------------------------------------------------------------------------
!  Determine the number of nodes outside the BL-region.
!------------------------------------------------------------------------

! Compute 'magp', the spacing of the last segment of the BL grid
! at a boundary node i = 1 as a reference.

   !Boundary node 1

    i = 1

   !  Second node

    xp = vspacing(1)*node_body(i)%nx + node_body(i)%x
    yp = vspacing(1)*node_body(i)%ny + node_body(i)%y
    zp = vspacing(1)*node_body(i)%nz + node_body(i)%z

   !  Third node to nm-th node
   !  Nodes in the prismatic layer along the surface normal.
   !  Just compute the last spacing (=magp).

    do k = 3, nm
     xx = vspacing(k-1)*node_body(i)%nx + xp
     yy = vspacing(k-1)*node_body(i)%ny + yp
     zz = vspacing(k-1)*node_body(i)%nz + zp
!     if (k==nm) then
      magp = sqrt( (xx-xp)**2 + (yy-yp)**2 +(zz-zp)**2 )
!     endif
     xp = xx
     yp = yy
     zp = zz
    end do

! Figure out 'nr' (the total number of nodes in the radial line)
! to match the spacing at far field, as well as match the spacing with 'magp'
! for a smooth transition from BL grid to outside.
! Note: There is no guarantee that the iteration converges to a unique value.
!       It probably depends on the initial guess.

     !My initial guess
      nr = nm + int(ceiling( (R_outer - rmax) / h_hemisphere )/ &
                (50.0_dp * ( 0.1_dp*real(nr_gs)/5.0_dp + 0.9_dp ) ))
      if (debug_mode) write(364,*)
      if (debug_mode) write(364,*) " Initial nr = ", nr

     !If requested (by a positive input for nre), use 'nre' specified by a user.
      if (nre > 0) then
       nr = nre + 1  !nre = elements, nr = nodes = elements+1.
      endif

   !--------------------------------------------
    itr_nr : do kk = 1, 100

      if (debug_mode) write(364,*)
      if (debug_mode) write(364,*) "--- Iteration for nr = ", kk

     ! Find the best stretching factor 'sf' to match the spacing with 'magp'
     ! for a smooth transition from BL grid to outside.

       k = nm+1
      xi = real( k - (nm) , dp ) / real( nr - (nm) , dp)

     !Determine sf iteratively: (one-exp(sf*xi))/(one-exp(sf))
        sf = 5.0_dp
       if (debug_mode) write(364,*) "        Initial stretching factor = ", sf
       do km = 1, 100
         gr = (one-exp(sf*xi))/(one-exp(sf))
         xx = gr*(node_outer(i)%x - xp)
         yy = gr*(node_outer(i)%y - yp)
         zz = gr*(node_outer(i)%z - zp)
        mag = sqrt( (xx)**2 + (yy)**2 + (zz)**2 )
        if (abs(mag - magp)/magp < 1.0e-01_dp) exit
        if (mag > magp) then
          sf = sf + 0.1_dp
        else
          sf = sf - 0.1_dp
        endif
       end do
       if (debug_mode) write(364,*) "                    sf iterations = ", km
       if (debug_mode) write(364,*) "                stretching factor = ", sf
       if (debug_mode) write(364,*) "                  Last BL spacing = ", magp
       if (debug_mode) write(364,*) "            First outside spacing = ", mag

     ! NOTE: This 'sf' is just an estimate.
     !       Actual stretching factor will be determined later.

     !Check the last spacing to see if it matches h_outer.

        k = nr
       xi = real( k - (nm) , dp ) / real( nr - (nm) , dp)
       gr = (one-exp(sf*xi))/(one-exp(sf))
       xx = gr*(node_outer(i)%x - xp)
       yy = gr*(node_outer(i)%y - yp)
       zz = gr*(node_outer(i)%z - zp)

        k = nr-1
       xi = real( k - (nm) , dp ) / real( nr - (nm) , dp)
       gr = (one-exp(sf*xi))/(one-exp(sf))
       xp = gr*(node_outer(i)%x - xp)
       yp = gr*(node_outer(i)%y - yp)
       zp = gr*(node_outer(i)%z - zp)

      mag = sqrt( (xx-xp)**2 + (yy-yp)**2 +(zz-zp)**2 )

      if (debug_mode) write(364,*) " Farfield spacing = ", h_outer
      if (debug_mode) write(364,*) "     Grid spacing = ", mag

      !If nr is specified, not iteration on 'nr'
      if (nre > 0) then 
        write(*,*)
        write(*,*) ">>> Stop iteration and use the requested value: nre = ", nre
        write(*,*) "    --->  nr(nodes) = nre + 1 = ", nr
        write(*,*) 
        exit itr_nr
      endif

      if (abs(mag - h_outer*spacing_ratio_outer)/(h_outer*spacing_ratio_outer)&
          < 1.0e-01) then
       exit itr_nr
      else
       if (mag < h_outer*spacing_ratio_outer) then
        nr = nr - 5 !If smaller spacing, reduce nr and try again.
       else
        nr = nr + 5 !If larger spacing, increase nr and try again.
       endif
      endif

     end do itr_nr
   !--------------------------------------------

     if (debug_mode) write(364,*) "       Farfield spacing = ", mag
     if (debug_mode) write(364,*) " Outer boundary spacing = ", h_outer
     if (debug_mode) write(364,*) "          Determined nr = ", nr
     if (debug_mode) write(364,*) "          at iterations = ", kk
     if (debug_mode .and. kk==1) write(364,*) " >>> Good initial guess!!!"  
     if (debug_mode) write(364,*)
   
   write(*,*)
   write(*,*) ">>> Total number of nodes in the radial direction:"
   write(*,*) "    --->  nr = ", nr
   write(*,*)
   write(*,*) "     Number of nodes outside the BL region = nr-nm = ", nr-nm
   write(*,*)

!------------------------------------------------------------------------
! Allocate arrays
!------------------------------------------------------------------------

  nnodes_body = nnodes
          nnodes = nnodes*nr
  allocate(node(nnodes))
  allocate(node_above(nr*nnodes_body))
  allocate(cell_above(3*nr*(ntrias+nquads_cyl)))
  allocate(cell_body(       ntrias+nquads_cyl ))
  allocate(k1(nnodes))
  allocate(k2(nnodes))
  allocate(k3(nnodes))
  allocate(k4(nnodes))
  allocate(k5(nnodes))

  write(*,*)
  write(*,*) ">>> Total number of nodes over the hemisphere-cylinder surface:"
  write(*,*) "    ---> nnodes_body =", nnodes_body
  write(*,*)

! Copy the data on the body
  do i = 1, nnodes_body
   node(i)%x  = node_body(i)%x
   node(i)%y  = node_body(i)%y
   node(i)%z  = node_body(i)%z
   node(i)%nx = node_body(i)%nx
   node(i)%ny = node_body(i)%ny
   node(i)%nz = node_body(i)%nz
   k1(i)  = k1_body(i)
   k2(i)  = k2_body(i)
   k3(i)  = k3_body(i)
   k4(i)  = k4_body(i)
   k5(i)  = 0          ! <- On the hemisphere cylinder
  end do

! Vector from a surface node to the corresponding outer-boundary node.

  do i = 1, nnodes_body
   node(i)%nx2 = node_outer(i)%x - node_body(i)%x
   node(i)%ny2 = node_outer(i)%y - node_body(i)%y
   node(i)%nz2 = node_outer(i)%z - node_body(i)%z
            sf = sqrt(node(i)%nx2**2 + node(i)%ny2**2 + node(i)%nz2**2)
   node(i)%nx2 = node(i)%nx2 / sf
   node(i)%ny2 = node(i)%ny2 / sf
   node(i)%nz2 = node(i)%nz2 / sf
  end do

!--------------------------------------------------------------------------------
! Go up and generate interior nodes.
!--------------------------------------------------------------------------------

  write(*,*)
  write(*,*) " Generating interior nodes.........."
  write(*,*)

 !--------------------
  nnodes = nnodes_body
  node_on_body : do i = 1, nnodes_body

 !---------------------------------------------------------------------
 ! Inside the specified boundary layer region, vspacing(1:nm) provides
 ! grid spacing in the radial direction.

 !  Second node
   nnodes = nnodes + 1
   node(nnodes)%x = vspacing(1)*node(i)%nx + node(i)%x
   node(nnodes)%y = vspacing(1)*node(i)%ny + node(i)%y
   node(nnodes)%z = vspacing(1)*node(i)%nz + node(i)%z
   xp = node(nnodes)%x
   yp = node(nnodes)%y
   zp = node(nnodes)%z
   node_above(i) = nnodes

    k1(nnodes)  = k1_body(i)
    k2(nnodes)  = k2_body(i)
    k3(nnodes)  = k3_body(i)
    k4(nnodes)  = k4_body(i)
    k5(nnodes)  = 1

 !  Third node to nm-th node
 !  Nodes in the prismatic layer along the surface normal. 
   do k = 3, nm

    nnodes = nnodes + 1
    node(nnodes)%x = vspacing(k-1)*node(i)%nx + xp
    node(nnodes)%y = vspacing(k-1)*node(i)%ny + yp
    node(nnodes)%z = vspacing(k-1)*node(i)%nz + zp
    node_above(nnodes-1) = nnodes

!    if (k==nm) then
     magp = sqrt( (node(nnodes)%x-xp)**2 + (node(nnodes)%y-yp)**2 +(node(nnodes)%z-zp)**2 )
!    endif

    xp = node(nnodes)%x
    yp = node(nnodes)%y
    zp = node(nnodes)%z

    k1(nnodes)  = k1_body(i)
    k2(nnodes)  = k2_body(i)
    k3(nnodes)  = k3_body(i)
    k4(nnodes)  = k4_body(i)
    k5(nnodes)  = k-1

   end do

 !---------------------------------------------------------------------
 !  Outer region
 !  Gradually deflected towards the second direction
 !  to avoid high-aspect ratio cells near the outer boundary

 ! Re-define the direction: from the edge of the boundary-layer region to outer.
   node(i)%nx2 = node_outer(i)%x - xp
   node(i)%ny2 = node_outer(i)%y - yp
   node(i)%nz2 = node_outer(i)%z - zp
          mag = sqrt( node(i)%nx2**2 + node(i)%ny2**2 + node(i)%nz2**2 )
   node(i)%nx2 = node(i)%nx2/mag
   node(i)%ny2 = node(i)%ny2/mag
   node(i)%nz2 = node(i)%nz2/mag

 ! Compute the actual stretching factor in xi for a smooth transition:

   iterate_for_sf = .true.

  !Set uniform_stretching_to_outer = .true. to avoid a gap which may be
  !created by a slight difference in the stretching factor. It can happen
  !for a large iteration tolerance, e.g., 1.0e-01. Or use a smaller tolerance.
   if (uniform_stretching_to_outer) then
    if (i > 1) iterate_for_sf = .false.
   endif

   if (iterate_for_sf) then

       k = nm + 1 
      xi = real( k - (nm) , dp ) / real( nr - (nm) , dp)
     !Determine sf iteratively (again): (one-exp(sf*xi))/(one-exp(sf))
      sf = 8.0_dp
     !if (debug_mode) write(*,*) " Initial stretching factor = ", sf
     do kk = 1, 100
       gr = (one-exp(sf*xi))/(one-exp(sf))
       xx = gr*(node_outer(i)%x - xp)
       yy = gr*(node_outer(i)%y - yp)
       zz = gr*(node_outer(i)%z - zp)
      mag = sqrt( (xx)**2 + (yy)**2 + (zz)**2 )
      if (abs(mag - magp)/magp < 1.0e-02_dp) exit !<--Smaller the better.
      if (mag > magp) then
        sf = sf + 0.1_dp
      else
        sf = sf - 0.1_dp
      endif
     end do
     if (debug_mode .and. i==1) write(364,*) "Actual stretching factor = ", sf
     if (debug_mode .and. i==1) write(364,*) "           at iterations = ", kk
     if (debug_mode .and. i==1) write(364,*) "           at   a node i = ", i

     if (debug_mode) write(365,*)
     if (debug_mode) write(365,*) "Actual stretching factor = ", sf
     if (debug_mode) write(365,*) "           at iterations = ", kk
     if (debug_mode) write(365,*) "           at   a node i = ", i

   endif


   do k = nm+1, nr

!   Uniform parameter function: e.g., xi = [0, 0.1, 0.2, 0.3, ..., 0.8, 0.9, 1]

     xi = real( k - (nm) , dp ) / real( nr - (nm) , dp)

    !Stretched parameter function
    ! xis = [0,0.01,0.015,0.02,...,0.5,0.8,1.0]
     xis = (one-exp(sf*xi))/(one-exp(sf))

!   Straight distance function
      xx = xis*(node_outer(i)%x - xp)
      yy = xis*(node_outer(i)%y - yp)
      zz = xis*(node_outer(i)%z - zp)
     mag = sqrt( xx**2 + yy**2 + zz**2 )

!   Define the blended direction based on inversely stretched parameter function.
    ! xi = [0,0.2,0.4,0.7,...,0.990,0.995,0.999,1.0]
      xi = (one-exp(-sf*xi))/(one-exp(-sf))
    dirx = (one-xi)*node(i)%nx + xi*node(i)%nx2
    diry = (one-xi)*node(i)%ny + xi*node(i)%ny2
    dirz = (one-xi)*node(i)%nz + xi*node(i)%nz2

!   Go up in the blended direction and generate a new interior node.
    nnodes = nnodes + 1
    node(nnodes)%x = xp + mag*dirx
    node(nnodes)%y = yp + mag*diry
    node(nnodes)%z = zp + mag*dirz
    node_above(nnodes-1) = nnodes

    k1(nnodes)  = k1_body(i)
    k2(nnodes)  = k2_body(i)
    k3(nnodes)  = k3_body(i)
    k4(nnodes)  = k4_body(i)
    k5(nnodes)  = k-1

    if (debug_mode .and. k==nr) then
     write(364,*) " xi = ", xi
     write(364,*) " node(nnodes)%x, node(nnodes)%y, node(nnodes)%z:"
     write(364,*)   node(nnodes)%x, node(nnodes)%y, node(nnodes)%z
     write(364,*) " node_outer(i)%x, node_outer(i)%y, node_outer(i)%z:"
     write(364,*)   node_outer(i)%x, node_outer(i)%y, node_outer(i)%z
     write(364,*)
    endif

   end do

  end do node_on_body
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

 deallocate(k1_body,k2_body,k3_body,k4_body)

! At this point, all nodes have been generated.
! It is now a matter of how to connect them (i.e., type of elements).

!*******************************************************************************
! 4.5. Let us rotate the grid by 90 degrees, so that a full and half domain
!      grids will be compatible: the generated half grid matches the half of
!      a full grid.
!
!      Also, we set y=0 precisely at the symmetry plane.
!
! Note: All y and z in other places in the code refer to the original (y,z),
!       even in the element construction performed after this rotation to
!       avoid any confusion.
!
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 4.5. Rotate the grid at 90 degrees, and enforce y=0.0 "
  write(*,*) "***********************************************************"

  write(*,*) " Performing the switch: y <- z and z <- (-y)....."
  write(*,*) " Also, enforcing y = 0.0 at the symmetry plane....."

 !--------------------------------------------------------------------
 !--------------------------------------------------------------------
  zero_y_nodes: do i = 1, nnodes

   !---------------------------------------------------------------
   !---------------------------------------------------------------
   ! Rotation

    yy =  node(i)%z
    zz = -node(i)%y

    node(i)%y = yy
    node(i)%z = zz

   !---------------------------------------------------------------
   !---------------------------------------------------------------
   ! Precisely zero the coordinates on the symmetry plane at y = 0.

    !**************************************
    !- Prism, tetra, mixed, mixed_ph

     if     (igrid_type  < 5) then

       !Symmetry plane is given by k2 = 0:
        if ( k2(i) == 0 ) then
          node(i)%y = zero
        endif

    !**************************************
    !- Structured

     elseif (igrid_type == 5) then

        !-------------------
        !- Full geometry

        if     (domain_cut == 1) then

         !Symmetry plane nodes are given by k1 = 1, 0, nnodes_circum/2 + 1 .
          if ( k1(i) == 1 .or. k1(i) == 0 .or. k1(i) == nnodes_circum/2 + 1 ) then
            node(i)%y = zero
          endif

        !-------------------
        !- Half geometry

        elseif (domain_cut == 2) then

         !Symmetry plane nodes are given by k1 = 1, 0, nnodes_circum
          if ( k1(i) == 1 .or. k1(i) == 0 .or. k1(i) == nnodes_circum ) then
            node(i)%y = zero
          endif

        endif

        !-------------------

     endif

    !**************************************

   !---------------------------------------------------------------
   !---------------------------------------------------------------

  end do zero_y_nodes
 !--------------------------------------------------------------------
 !--------------------------------------------------------------------

  write(*,*)
  write(*,*) " The grid has been rotated by 90 degrees from +ve z-axis to +ve y-axis."
  write(*,*) " Enforced y=0.0 at symmetry plane."
  write(*,*)

!*******************************************************************************
! 5. Write a map file: boundary marks
!
! There are three boundary parts (full-geometry case):
! 1. Hemisphere cylinder body
! 2. Outflow plane
! 3. Outer boundary
!
! In the case of a half geometry:
! 1. Hemisphere cylinder body
! 2. Outflow plane
! 3. Outer boundary
! 4. Symmetry plane at y = 0
!
! NOTE: The boundary condition numbers (e.g., 4000) are specific to a solver.
!       Appropriate number needs to be assigned for your solver.
!
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 5. Write a boundary info file"
  write(*,*) "***********************************************************"

  open(unit=16, file=filename_mapbc, status="unknown", iostat=os)

 !Full domain
  full_geometry_mapbc : if (n_sections == 6) then

   write(*,*) " Full geometry -> 3 boundary groups."
   write(*,*) "   --- (1)Hemisphere-cylinder."
   write(*,*) "   --- (2)Outflow."
   write(*,*) "   --- (3)Farfield"

   write(16,'(a57)') "3       !Number of boundary parts (boundary conditions)"
   write(16,'(a32)') "1, 4000 !Viscous wall in FUN3D"
   write(16,'(a46)') "2, 5051 !Outflow with back pressure in FUN3D"
   write(16,'(a55)') "3, 5050 !Characteristic-based inflow/outflow in FUN3D"

 !Half-domain
  elseif (n_sections == 3) then

   write(*,*) " Partial geometry -> 5 boundary groups."
   write(*,*) "   --- (1)Hemisphere-cylinder."
   write(*,*) "   --- (2)Outflow."
   write(*,*) "   --- (3)Farfield"
   write(*,*) "   --- (4)y-symmetry"

   ! Here we use boundaries types supported by EDU3D
   write(16,*) '4                        ! Number of boundary parts (boundary conditions)'
   write(16,*) '1, "slip_wall"           ! Viscous wall (4000 in FUN3D)'
   write(16,*) '2, "outflow_subsonic_p0" ! Outflow with back pressure (5051 in FUN3D)'
   write(16,*) '3, "freestream"          ! Characteristic-based inflow/outflow (5050 in FUN3D)'
   write(16,*) '4, "slip_wall"           ! Y-Symmetry (6662 in FUN3D)'

 !Partial domain
  else

   write(*,*) " Partial geometry -> 5 boundary groups."
   write(*,*) "   --- (1)Hemisphere-cylinder."
   write(*,*) "   --- (2)Outflow."
   write(*,*) "   --- (3)Farfield"
   write(*,*) "   --- (4)Symmetry 1"
   write(*,*) "   --- (5)Symmetry 2"

   write(16,'(a57)') "5       !Number of boundary parts (boundary conditions)"
   write(16,'(a32)') "1, 4000 !Viscous wall in FUN3D"
   write(16,'(a46)') "2, 5051 !Outflow with back pressure in FUN3D"
   write(16,'(a55)') "3, 5050 !Characteristic-based inflow/outflow in FUN3D"
   write(16,'(a22)') "4, 6662 !Y-Symmetry?"
   write(16,'(a22)') "5, 6662 !Y-Symmetry?"

  endif full_geometry_mapbc

  close(16)

  write(*,*)
  write(*,*) " Boundary info file = ", filename_mapbc

!*******************************************************************************
! 6. Generate elements
!
!    Different grids will be generated by connecting existing nodes.
!
!    igrid_type = Choice of element type:
!
!      1 = Prismatic grid
!      2 = Tetrahedral grid
!      3 = Mixed grid (Prism/Tet)
!      4 = Mixed grid (Prism/Hex)
!      5 = Strct grid (Hex + some prisms)
!
!*******************************************************************************
!*******************************************************************************
!* (1). Generate a prismatic grid
!*******************************************************************************
 if (igrid_type == 1) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 6. Generate a prismatic grid"
  write(*,*) "***********************************************************"

 call prismatic_grid

!*******************************************************************************
!* (2). Generate a tet grid
!*******************************************************************************
 elseif ( igrid_type == 2 ) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 6. Generate a tetrahedral grid"
  write(*,*) "***********************************************************"

  call tet_grid

!*******************************************************************************
!* (3). Generate a mixed grid  (Prism/Tet)
!*******************************************************************************
 elseif ( igrid_type == 3 ) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 6. Generate a mixed grid: Prism/Tetra"
  write(*,*) "***********************************************************"

 call mixed_grid

 write(*,*) " Mixed has been generated."

!*******************************************************************************
!* (4). Generate a mixed grid (Prism/Hex)
!*******************************************************************************
 elseif ( igrid_type == 4 ) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 6. Generate a mixed grid: Prism/Hex"
  write(*,*) "***********************************************************"

 call mixed_ph_grid

 write(*,*) " Mixed has been generated."

!*******************************************************************************
!* (5). Generate a structured grid (Hex + some prisms)
!*******************************************************************************
 elseif ( igrid_type == 5 ) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 6. Generate a structured grid: Hex + some prisms"
  write(*,*) "***********************************************************"

 call mixed_strct_grid

 write(*,*) " Structured has been generated."

!*******************************************************************************
!* (6). Error
!*******************************************************************************
 else

  write(*,*) "Invalid input: igrid_type = ", igrid_type
  write(*,*) "               igrid_type must be 1, 2, 3, 4, or 5. Try again."
  stop

 endif

!*******************************************************************************
!*******************************************************************************
!
! Negative volume check.
!
! Note: The code will proceed even if negative volumes are detected.
!       Some methods work just fine with zero/negative-volume elements.
!
!       Ref.: H. Nishikawa, "Uses of Zero and Negative Volume Elements for
!             Node-Centered Edge-Based Discretization", AIAA Paper 2017-4295,
!             23rd AIAA CFD Conference, 5 - 9 June 2017, Denver, Colorado.
!             http://hiroakinishikawa.com/My_papers/nishikawa_aiaa2017-4295.pdf
!
!*******************************************************************************
!*******************************************************************************

  write(*,*)
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*)
  write(*,*) ">>> Check element volumes..."
  write(*,*)
  write(*,*) "  # of Tetrahedra = ", ntet
  write(*,*) "  # of Hexahedra  = ", nhex
  write(*,*) "  # of Prisms     = ", nprs
  write(*,*)

  call volume_check(negative_volume_detected)

  !negative_volume_detected(1) = T if negative Tet   detected.
  !negative_volume_detected(2) = T if negative Hexa  detected.
  !negative_volume_detected(3) = T if negative Prism detected.

  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*)

!*******************************************************************************
!*******************************************************************************
!
! Check non-planerity of quadrilateral element-faces.
!
!*******************************************************************************
!*******************************************************************************

  write(*,*)
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*)
  write(*,*) ">>> Check non-planerity of quadrilateral faces..."
  write(*,*)

  if (nhex == 0 .and. nprs == 0) then

   write(*,*) " Well, this is a pure tet grid, and so, "
   write(*,*) " no quadrilateral faces exist. Skip the check! "

  else

   call non_planarity_check

  endif

  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  write(*,*)


!*******************************************************************************
! 7. Generate two line files:
!
!    (1) hemisphere.lines_fmt     - Lines are within the viscous layer.
!    (2) hemisphere.lines_fmt_all - Lines go up to the outer boundary.
!
!  Note: Line files need to be generated if one wants to apply regular coarsening
!        afterwards.
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 7.1 NC: Generate line-info files"
  write(*,*) "***********************************************************"

  if( generate_line_file_nc ) then

! Lines within the viscous layer

  write(*,*)
  write(*,*) "Writing a regular BL line-info file...."
  write(*,*) " - Points within the BL region only."
  write(*,*)
  call line_info_file(nm, filename_lines)

! Lines go all the way up to the outer boundary (all nodes are in lines)

  write(*,*)
  write(*,*) "Writing an all-node line-info file...."
  write(*,*) " - All points in a line from the HC body to the outer boundary."
  write(*,*)
  call line_info_file_all(nr, filename_lines_all)

  else
    write(*,*) "Skipping write of line-info files"
  endif

!------------------------------------------------------------------------------
! Line information for cells:
!------------------------------------------------------------------------------

  if( generate_line_file_cc ) then

  !Compute the cell centroid coordinates.

  allocate(cell(ntet+nprs+nhex))

    if (ntet > 0) then
     do i = 1, ntet
      cell(i)%x = (node(tet(i,1))%x+node(tet(i,2))%x+node(tet(i,3))%x+node(tet(i,4))%x)/4.0_dp
      cell(i)%y = (node(tet(i,1))%y+node(tet(i,2))%y+node(tet(i,3))%y+node(tet(i,4))%y)/4.0_dp
      cell(i)%z = (node(tet(i,1))%z+node(tet(i,2))%z+node(tet(i,3))%z+node(tet(i,4))%z)/4.0_dp
     end do
    endif

    if (nprs > 0) then
     do i = 1, nprs
      cell(ntet+i)%x = (node(prs(i,1))%x+node(prs(i,2))%x+node(prs(i,3))%x + &
                        node(prs(i,4))%x+node(prs(i,5))%x+node(prs(i,6))%x   )/6.0_dp
      cell(ntet+i)%y = (node(prs(i,1))%y+node(prs(i,2))%y+node(prs(i,3))%y + &
                        node(prs(i,4))%y+node(prs(i,5))%y+node(prs(i,6))%y   )/6.0_dp
      cell(ntet+i)%z = (node(prs(i,1))%z+node(prs(i,2))%z+node(prs(i,3))%z + &
                        node(prs(i,4))%z+node(prs(i,5))%z+node(prs(i,6))%z   )/6.0_dp
     end do
    endif

    if (nhex > 0) then
     do i = 1, nhex
      cell(ntet+nprs+i)%x = ( &
                   node(hex(i,1))%x+node(hex(i,2))%x+node(hex(i,3))%x + &
                   node(hex(i,4))%x+node(hex(i,5))%x+node(hex(i,6))%x + &
                   node(hex(i,7))%x+node(hex(i,8))%x                    )/8.0_dp
      cell(ntet+nprs+i)%y = ( &
                   node(hex(i,1))%y+node(hex(i,2))%y+node(hex(i,3))%y + &
                   node(hex(i,4))%y+node(hex(i,5))%y+node(hex(i,6))%y + &
                   node(hex(i,7))%y+node(hex(i,8))%y                    )/8.0_dp
      cell(ntet+nprs+i)%z = ( &
                   node(hex(i,1))%z+node(hex(i,2))%z+node(hex(i,3))%z + &
                   node(hex(i,4))%z+node(hex(i,5))%z+node(hex(i,6))%z + &
                   node(hex(i,7))%z+node(hex(i,8))%z                    )/8.0_dp
     end do
    endif

  write(*,*) "--------------------------------------------"
  write(*,*) " 7-2. CC: Generate lin-info files"
  write(*,*) "--------------------------------------------"

 ! Lines within the viscous layer

  write(*,*)
  write(*,*) "Writing a regular BL line-info file...."
  write(*,*) " - (nm-1) within the BL region only: (nm-1) = ",(nm-1)
  write(*,*)

  if (igrid_type==2) then !Tetra
   call line_info_file_c(3*(nm-1), filename_lines_c)
  else
   call line_info_file_c(nm-1, filename_lines_c)
  endif

  !These lines are tiny, and hard to see...
  if (debug_mode) then
   open(unit=3, file=filename11, status="unknown", iostat=os)
    write(3,*) 'GEOMETRY T=LINE3D'
    write(3,*) min(10,ntrias + nquads_cyl)
    do i = 1, min(10,ntrias + nquads_cyl)
     write(3,*) nm-1
     write(3,*) cell(cell_body(i))%x, cell(cell_body(i))%y, cell(cell_body(i))%z
     n_temp = cell_body(i)
     do k = 2, nm-1
     write(3,*) cell(cell_above(n_temp))%x, cell(cell_above(n_temp))%y, cell(cell_above(n_temp))%z
     n_temp = cell_above(n_temp)
     end do
    end do
   close(3)
  endif

 ! Lines go all the way up to the outer boundary (all nodes are in lines)

  write(*,*)
  write(*,*) "Writing an all-node line-info file...."
  write(*,*) " - All points in a line from the HC to the outer boundary."
  write(*,*)

  if     (igrid_type==2) then !Tetra

   call line_info_file_c( 3*(nr-1), filename_lines_c_all)

  elseif (igrid_type==3) then !Mixed

  !                       prisms        tetra
   call line_info_file_c( (nm-1) + 3*( (nr-1)-(nm-1) ), filename_lines_c_all)

  else

   call line_info_file_c( (nr-1), filename_lines_c_all)

  endif

  !These lines are OK. But Tecplot allows only 50 lines.
  if (debug_mode) then
   open(unit=3, file=filename12, status="unknown", iostat=os)
    write(3,*) 'GEOMETRY T=LINE3D'
    write(3,*) min(50,ntrias + nquads_cyl)
    do i = 100, 99+min(50,ntrias + nquads_cyl)
     write(3,*) nr-1
     write(3,*) cell(cell_body(i))%x, cell(cell_body(i))%y, cell(cell_body(i))%z
     n_temp = cell_body(i)
     do k = 2, nr-1
     write(3,*) cell(cell_above(n_temp))%x, cell(cell_above(n_temp))%y, cell(cell_above(n_temp))%z
     n_temp = cell_above(n_temp)
     end do
    end do
   close(3)
  endif

    deallocate(cell)

  else
    write(*,*) "Skipping write of line-info files"
  endif

!*******************************************************************************
! 8. Write a ugrid file.
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 8. Generate a k file"
  write(*,*) "***********************************************************"

  if (generate_k_file) then
    write(*,*) "Writing a k file...."
    call write_k_file
    write(*,*) "K file written : ", filename_k
  else
    write(*,*) "Skipping write of K file : ", filename_k
  endif

!*******************************************************************************
! 9. Write a ugrid file.
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 9. Generate a .ugrid file for running a flow solver"
  write(*,*) "***********************************************************"

  if (generate_ugrid_file) then
   write(*,*) "Writing a ugrid file...."
   call write_ugrid_file
   write(*,*) "UGRID file written : ", filename_ugrid
  else
    write(*,*) "Skipping write of .ugrid file : ", filename_ugrid
  endif

!*******************************************************************************
! 10. Write a Tecplot file for boundaries.
!*******************************************************************************

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 10. Generate a tecplot file for viewing boundaries"
  write(*,*) "***********************************************************"

  if (generate_tec_file_b) then

   write(*,*) "Writing a Tecplot file for boundaries...."
   call write_tecplot_boundary_file
   write(*,*) "Tecplot file for boundaries written : ", filename_tecplot_b

  else

   write(*,*) "Skipping write of Teplot boundary file : ", filename_tecplot_b

  endif

!*******************************************************************************
! 11. Write a Tecplot file for the volume grid.
!*******************************************************************************

   write(*,*)
   write(*,*) "***********************************************************"
   write(*,*) " 11. Generate a tecplot file for viewing a volume grid"
   write(*,*) "***********************************************************"

  if (generate_tec_file_v) then

   write(*,*) "Writing a Tecplot file for volumes...."
   call write_tecplot_volume_file
   write(*,*) "Tecplot file for volume grid written : ", filename_tecplot_v

  else
   write(*,*)
   write(*,*) "Skipping write of Teplot volume file : ", filename_tecplot_v
   write(*,*)
  endif

!*******************************************************************************
! 12. Write a PLOT3D file (for structured-grid only: igrid_type=5).
!*******************************************************************************

  if (igrid_type == 5 .and. generate_p3d_ufmt_file) then

  write(*,*)
  write(*,*) "***********************************************************"
  write(*,*) " 12. Write a PLOT3D file"
  write(*,*) "***********************************************************"

     write(*,*) "Writing a PLOT3D file...."
     call write_plot3d_file
     write(*,*) "PLOT3D file written : ", filename_p3d
     write(*,*) "  .nmf file written : ", filename_nmf

  else

   if (igrid_type == 5) then
    write(*,*) "Skipping write of p3d/ufmt file : ", filename_p3d
   endif

  endif


 write(*,*)
 write(*,*) "Congratulations!"
 write(*,*) "Grid generation successfully completed."

!*******************************************************************************
! Report negative volumes if any...
!*******************************************************************************

  if (negative_volume_detected(1)) then
   write(*,*)
   write(*,*) "  Well, sorry to tell you that negative Tet volume elements exist..."
   write(*,*) "   - Go up and check the output message..."
   write(*,*) "   - See also fort.8933 for details......."
   write(*,*) "   - First 20 negative volumes can be viewed: negative_vol_tet_xx_tec.dat"
   write(*,*)
  endif

  if (negative_volume_detected(2)) then
   write(*,*)
   write(*,*) "  Well, sorry to tell you that negative Hex volume elements exist..."
   write(*,*) "   - Go up and check the output message..."
   write(*,*) "   - See also fort.8933 for details......."
   write(*,*) "   - First 20 negative volumes can be viewed: negative_vol_hex_xx_tec.dat"
   write(*,*)
  endif

  if (negative_volume_detected(3)) then
   write(*,*)
   write(*,*) "  Well, sorry to tell you that negative Prism volume elements exist..."
   write(*,*) "   - Go up and check the output message..."
   write(*,*) "   - See also fort.8933 for details......."
   write(*,*) "   - First 20 negative volumes can be viewed: negative_vol_prs_xx_tec.dat"
   write(*,*)
  endif

!*******************************************************************************
! (Debug) Some additional information about a tetrahedral grid.
!*******************************************************************************

  if (grid_statistics .and. igrid_type==2) then
   write(*,*) " Computing some statistics..."
   call tet_statistics(filename_stats)
   write(*,*) " Finished. See ", trim(filename_stats)
  endif

  ! Cleanup debug hooks
  call edu3d_euler_debug_hooks_cleanup()

contains 



!*******************************************************************************
!* Prismatic Grid Generation
!*******************************************************************************
 subroutine prismatic_grid
 implicit none
 integer(kd) :: cell_below, ntet0, nhex0, nprs0, ntrias_b0, nquads_b0
!*******************************************************************************
! Generate prismatic grid
!*******************************************************************************
  ntrias_b0 = ntrias + ntrias        !Inner and outer boundaries
  allocate(tri(ntrias_b0,5))

  if (n_sections==6) then
   nquads_b0 = nnodes_circum * (nr-1) !Outflow boundary
   allocate(quad(nquads_b0,5))
  else !if (n_sections==1) then
   nquads_b0 = nnodes_circum * (nr-1) + 2*(nr_gs+nnodes_cylinder)*(nr-1)
   allocate(quad(nquads_b0,5))
  endif

  ntet0 = 0
  nprs0 = ntrias*(nr-1)
  allocate(prs(nprs0,6) )
  nhex0 = 0

  write(*,*) "-- Prismatic grid : Expected values..."
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs0
  write(*,*) "    tetrahedra = ", ntet0
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b0
  write(*,*) "         quads = ", nquads_b0
  write(*,*)

  nprs = 0
  ntrias_b = 0
  nquads_b = 0

! Copy the triangulation on the body

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate prisms

  do i = 1, ntrias

   node1 = tria(i)%v(1)
   node2 = tria(i)%v(2)
   node3 = tria(i)%v(3)
   node4 = node_above(tria(i)%v(1))
   node5 = node_above(tria(i)%v(2))
   node6 = node_above(tria(i)%v(3))

     nprs = nprs + 1

     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_body(i) = nprs
     cell_below   = nprs

   do k = 2, nr-1

    node1 = node4
    node2 = node5
    node3 = node6
    node4 = node_above(node1)
    node5 = node_above(node2)
    node6 = node_above(node3)

    nprs = nprs + 1
    prs(nprs,1) = node1
    prs(nprs,2) = node2
    prs(nprs,3) = node3
    prs(nprs,4) = node4
    prs(nprs,5) = node5
    prs(nprs,6) = node6

    cell_above(cell_below) = nprs
    cell_below             = nprs

   end do

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,5) = 3
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

  end do

! Outflow boundary
   nquads_b = 0

  if     (n_sections==6) then
   n_temp = nnodes_circum
  else !if (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

   do i = 1, n_temp

    node1 = nodes_cylinder(i  , nnodes_cylinder+1)
    node2 = nodes_cylinder(i+1, nnodes_cylinder+1)
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 2
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 2
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 section1 : if (n_sections/=6) then

 ! Symmetry boundary 1: z=0 and positive y.

   do i = 1, nr_gs

    node1 = (i-1)*i/2 + 1 !Left-most node of the original sector
    node2 = (i+1)*i/2 + 1 !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do


   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(1, i  ) !Left-most node of the original sector
    node2 = nodes_cylinder(1, i+1) !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 ! Symmetry boundary 2: z=0 and negative y.

   do i = 1, nr_gs

    node1 = (i+1)*i/2 + i +1 !Right-most node of the original sector
    node2 = (i-1)*i/2 + i    !Right-most node of the original sector

   !Map onto the last sector:
    node1 = node_map((i+1)*i/2 + i +1)
    node2 = node_map((i-1)*i/2 + i)

    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(nnodes_circum, i+1)
    node2 = nodes_cylinder(nnodes_circum, i  )
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 endif section1

  write(*,*) " -----------------------------------------------------"
  write(*,*) " Generated elements....."
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    hex        = ", nhex
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*) " -----------------------------------------------------"
  write(*,*)
  write(*,*) " check(tetra ) = ", ntet0-ntet
  write(*,*) " check(prisms) = ", nprs0-nprs
  write(*,*) " check(hexa  ) = ", nhex0-nhex
  write(*,*) " check(tria  ) = ", ntrias_b0-ntrias_b
  write(*,*) " check(quads ) = ", nquads_b0-nquads_b
  write(*,*)

 end subroutine prismatic_grid
!********************************************************************************


!*******************************************************************************
!* Tet grid Generation
!
!*******************************************************************************
 subroutine tet_grid
 implicit none
 integer(kd) :: tria_type
 integer(kd) :: cell_below, ntet0, nhex0, nprs0, ntrias_b0, nquads_b0
!*******************************************************************************
! Generate tet grid
!*******************************************************************************

  if (n_sections==6) then
   ntrias_b0 = ntrias + ntrias + 2*nnodes_circum * (nr-1) !Inner and outer boundaries
   allocate(tri(ntrias_b0,5))
  else !if (n_sections==1) then
  !Note: nnodes_circum is the number of nodes in case of half geometry.
   ntrias_b  = ntrias + ntrias + 2*(nnodes_circum-1) * (nr-1) !Inner and outer boundaries
   ntrias_b0 = ntrias_b + 2*2*(nr_gs+nnodes_cylinder)*(nr-1)
   allocate(tri(ntrias_b0,5))
  endif

!  ntrias_b = ntrias + ntrias + 2*nnodes_circum * (nr-1) !Inner and outer boundaries
!  allocate(tri(ntrias_b,5))

  nquads_b0 = 0

  ntet0 = ntrias*(nr-1)*3
  allocate(tet(ntet0,6) )
  nprs0 = 0
  nhex0 = 0

  write(*,*) "-- Tetrahedral grid : Expected values... "
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs0
  write(*,*) "    tetrahedra = ", ntet0
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b0
  write(*,*) "         quads = ", nquads_b0

  nprs = 0
  ntet = 0
  ntrias_b = 0
  nquads_b = 0

! Copy the triangulation on the body.

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate tetrahedra

  do i = 1, ntrias

     tria_type = tria(i)%type

   do k = 1, nr-1

    if (k == 1) then

     node1 = tria(i)%v(1)
     node2 = tria(i)%v(2)
     node3 = tria(i)%v(3)
     node4 = node_above(tria(i)%v(1))
     node5 = node_above(tria(i)%v(2))
     node6 = node_above(tria(i)%v(3))

    else

     node1 = node4
     node2 = node5
     node3 = node6
     node4 = node_above(node1)
     node5 = node_above(node2)
     node6 = node_above(node3)

    endif

 !--------------------------------------------
     t_type :if (tria_type == 1) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      if (k==1) then
       cell_body(i) = ntet
      else
       cell_above(cell_below) = ntet
      endif
       cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      cell_above(cell_below) = ntet
      cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node4

      cell_above(cell_below) = ntet
      cell_below   = ntet

 !--------------------------------------------
     elseif (tria_type == 2) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      if (k==1) then
       cell_body(i) = ntet
      else
       cell_above(cell_below) = ntet
      endif
       cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below   = ntet

 !--------------------------------------------
     elseif (tria_type == 20) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      if (k==1) then
       cell_body(i) = ntet
      else
       cell_above(cell_below) = ntet
      endif
       cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below   = ntet

 !--------------------------------------------
     elseif (tria_type == 3) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node5

      if (k==1) then
       cell_body(i) = ntet
      else
       cell_above(cell_below) = ntet
      endif
       cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node4
      tet(ntet,3) = node3
      tet(ntet,4) = node5

      cell_above(cell_below) = ntet
      cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node1
      tet(ntet,3) = node3
      tet(ntet,4) = node5

      cell_above(cell_below) = ntet
      cell_below   = ntet

 !--------------------------------------------
     elseif (tria_type == 4) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node6

      if (k==1) then
       cell_body(i) = ntet
      else
       cell_above(cell_below) = ntet
      endif
       cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node5
      tet(ntet,3) = node2
      tet(ntet,4) = node6

      cell_above(cell_below) = ntet
      cell_below   = ntet

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node2
      tet(ntet,3) = node1
      tet(ntet,4) = node6

      cell_above(cell_below) = ntet
      cell_below   = ntet

     endif t_type

   end do

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,5) = 3
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

  end do

! Outflow boundary
   nquads_b = 0

  if     (n_sections==6) then
   n_temp = nnodes_circum
  else !if (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

   do i = 1, n_temp

    do k = 1, nr-1

     if (k == 1) then

      node1 = nodes_cylinder(i  , nnodes_cylinder+1)
      node2 = nodes_cylinder(i+1, nnodes_cylinder+1)
      node3 = node_above(node1)
      node4 = node_above(node2)

     else

      node1 = node3
      node2 = node4
      node3 = node_above(node1)
      node4 = node_above(node2)

     endif

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 2
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 2
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    end do

   end do

  section1 : if (n_sections/=6) then

 ! Symmetry boundary 1: z=0 and positive y.

   do i = 1, nr_gs

    node1 = (i-1)*i/2 + 1
    node2 = (i+1)*i/2 + 1
    node3 = node_above(node1)
    node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    end do

   end do

   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(1, i  ) !Left-most node of the original sector
    node2 = nodes_cylinder(1, i+1) !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node3

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node2
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node3

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node2
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    end do

   end do
   
 ! Symmetry boundary 2: z=0 and negative y.

   do i = 1, nr_gs

    node1 = (i+1)*i/2 + i+1 !Right-most node of the original sector
    node2 = (i-1)*i/2 + i   !Right-most node of the original sector

   !Map onto the last sector:
    node1 = node_map((i+1)*i/2 + i +1)
    node2 = node_map((i-1)*i/2 + i)

    node3 = node_above(node1)
    node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node3

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node2
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node3

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node2
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    end do

   end do

   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(nnodes_circum, i+1) !Left-most node of the original sector
    node2 = nodes_cylinder(nnodes_circum, i  ) !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

    end do

   end do

  endif section1

  write(*,*) " -----------------------------------------------------"
  write(*,*) " Generated elements....."
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    hex        = ", nhex
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*) " -----------------------------------------------------"
  write(*,*)
  write(*,*) " check(tetra ) = ", ntet0-ntet
  write(*,*) " check(prisms) = ", nprs0-nprs
  write(*,*) " check(hexa  ) = ", nhex0-nhex
  write(*,*) " check(tria  ) = ", ntrias_b0-ntrias_b
  write(*,*) " check(quads ) = ", nquads_b0-nquads_b
  write(*,*)

 end subroutine tet_grid
!********************************************************************************


!*******************************************************************************
!* Mixed grid Generation
!
!*******************************************************************************
 subroutine mixed_grid
 implicit none
 integer(kd) :: tria_type, ntet0, nhex0, nprs0, ntrias_b0, nquads_b0
 integer(kd) :: cell_below
!*******************************************************************************
! Generate mixed grid
!*******************************************************************************

  ntrias_b0 = 2*ntrias + 2*nnodes_circum * ( (nr-1)-(nm-1) )  !All boundaries
!  allocate(tri(ntrias_b,5))
  nquads_b0 = nnodes_circum * (nm-1) !A part of outflow boundary
!  allocate(quad(nquads_b,5))

  ntet0 = ntrias*((nr-1)-(nm-1))*3   !Outer region
  allocate(tet(ntet0,6) )
  nprs0 = ntrias*(nm-1)              !Boundary layer
  allocate(prs(nprs0,6) )
  nhex0 = 0

  if (n_sections/=6) then
   ntrias_b0 = ntrias_b0 + 2*2*(nr_gs+nnodes_cylinder)*(nr-1)
   nquads_b0 = nquads_b0 +   2*(nr_gs+nnodes_cylinder)*(nr-1)
  endif

  allocate(tri(ntrias_b0,5))
  allocate(quad(nquads_b0,5))

  write(*,*) "-- Mixed grid: Expected values... "
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet0
  write(*,*) "    prisms     = ", nprs0
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b0
  write(*,*) "         quads = ", nquads_b0

  nprs = 0
  ntet = 0
  ntrias_b = 0
  nquads_b = 0

! Copy the triangulation on the body

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate prisms

  do i = 1, ntrias

   node1 = tria(i)%v(1)
   node2 = tria(i)%v(2)
   node3 = tria(i)%v(3)
   node4 = node_above(tria(i)%v(1))
   node5 = node_above(tria(i)%v(2))
   node6 = node_above(tria(i)%v(3))
   tria_type = tria(i)%type

     nprs = nprs + 1

     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_body(i) = ntet0 + nprs
     cell_below   = ntet0 + nprs

   do k = 2, nr-1

    node1 = node4
    node2 = node5
    node3 = node6
    node4 = node_above(node1)
    node5 = node_above(node2)
    node6 = node_above(node3)

    prs_or_tet : if (k < nm) then

     nprs = nprs + 1
     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_above(cell_below) = ntet0 + nprs
     cell_below             = ntet0 + nprs

    else

     t_type :if (tria_type == 1) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node4

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node4

      cell_above(cell_below) = ntet
      cell_below             = ntet

     elseif (tria_type == 2) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below             = ntet

     elseif (tria_type == 20) then

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node6
      tet(ntet,3) = node5
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node6
      tet(ntet,3) = node3
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node5
      tet(ntet,2) = node3
      tet(ntet,3) = node2
      tet(ntet,4) = node1

      cell_above(cell_below) = ntet
      cell_below             = ntet

     elseif (tria_type == 3) then

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node2
      tet(ntet,3) = node3
      tet(ntet,4) = node5

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node4
      tet(ntet,3) = node3
      tet(ntet,4) = node5

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node4
      tet(ntet,2) = node1
      tet(ntet,3) = node3
      tet(ntet,4) = node5

      cell_above(cell_below) = ntet
      cell_below             = ntet

     elseif (tria_type == 4) then

      ntet = ntet + 1
      tet(ntet,1) = node6
      tet(ntet,2) = node5
      tet(ntet,3) = node4
      tet(ntet,4) = node2

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node3
      tet(ntet,3) = node6
      tet(ntet,4) = node2

      cell_above(cell_below) = ntet
      cell_below             = ntet

      ntet = ntet + 1
      tet(ntet,1) = node1
      tet(ntet,2) = node6
      tet(ntet,3) = node4
      tet(ntet,4) = node2

      cell_above(cell_below) = ntet
      cell_below             = ntet

     endif t_type

    endif prs_or_tet

   end do

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,5) = 3
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

  end do

! Outflow boundary
   nquads_b = 0

  if     (n_sections==6) then
   n_temp = nnodes_circum
  else !if (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

   do i = 1, n_temp

    node1 = nodes_cylinder(i  , nnodes_cylinder+1)
    node2 = nodes_cylinder(i+1, nnodes_cylinder+1)
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 2
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     if (k < nm) then

      nquads_b = nquads_b + 1
      quad(nquads_b,5) = 2
      quad(nquads_b,1) = node1
      quad(nquads_b,2) = node2
      quad(nquads_b,3) = node4
      quad(nquads_b,4) = node3

     else

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 2
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 2
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

     endif

    end do

   end do


  section1 : if (n_sections/=6) then

 ! Symmetry boundary 1: z=0 and positive y.

   do i = 1, nr_gs

    node1 = (i-1)*i/2 + 1
    node2 = (i+1)*i/2 + 1
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     if (k < nm) then

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

     else

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

     endif

    end do

   end do


   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(1, i  ) !Left-most node of the original sector
    node2 = nodes_cylinder(1, i+1) !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     if (k < nm) then

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

     else

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node3

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node2
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

     endif

    end do

   end do

 ! Symmetry boundary 2: z=0 and negative y.

   do i = 1, nr_gs

    node1 = (i+1)*i/2 + i+1 !Right-most node of the original sector
    node2 = (i-1)*i/2 + i   !Right-most node of the original sector

   !Map onto the last sector:
    node1 = node_map((i+1)*i/2 + i +1)
    node2 = node_map((i-1)*i/2 + i)

    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     if (k < nm) then

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

     else

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node3

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node2
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

     endif

    end do

   end do


   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(nnodes_circum, i+1) !Left-most node of the original sector
    node2 = nodes_cylinder(nnodes_circum, i  ) !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     if (k < nm) then

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

     else

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node2
      tri(ntrias_b,3) = node4

      ntrias_b = ntrias_b + 1
      tri(ntrias_b,5) = 5
      if (n_sections==3) tri(ntrias_b,5) = 4
      tri(ntrias_b,1) = node1
      tri(ntrias_b,2) = node4
      tri(ntrias_b,3) = node3

     endif

    end do

   end do

  endif section1

  write(*,*) " -----------------------------------------------------"
  write(*,*) " Generated elements....."
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    hex        = ", nhex
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*) " -----------------------------------------------------"
  write(*,*)
  write(*,*) " check(tetra ) = ", ntet0-ntet
  write(*,*) " check(prisms) = ", nprs0-nprs
  write(*,*) " check(hexa  ) = ", nhex0-nhex
  write(*,*) " check(tria  ) = ", ntrias_b0-ntrias_b
  write(*,*) " check(quads ) = ", nquads_b0-nquads_b
  write(*,*)

 end subroutine mixed_grid
!********************************************************************************


!*******************************************************************************
!* Mixed grid Generation: Prism and Hex
!
!*******************************************************************************
 subroutine mixed_ph_grid
 implicit none
 integer(kd) :: node7, node8, ntet0, nhex0, nprs0, ntrias_b0, nquads_b0
 integer(kd) :: cell_below
!*******************************************************************************
! Generate mixed grid
!*******************************************************************************
  ntrias_b0 = 2*ntrias !All boundaries
  allocate(tri(ntrias_b0,5))

  if (n_sections==6) then
   nquads_b  = nnodes_circum * (nr-1)   ! All over outflow boundary
   nquads_b0 = nquads_b + 2*nquads_cyl  ! + All over cylinder surface and outer
   allocate(quad(nquads_b0,5))
  else !if (n_sections==1) then
   nquads_b = nnodes_circum * (nr-1)   ! All over outflow boundary
   nquads_b0 = nquads_b + 2*nquads_cyl  ! + All over cylinder surface and outer
   nquads_b0 = nquads_b0 + 2*(nr_gs+nnodes_cylinder)*(nr-1)
   allocate(quad(nquads_b0,5))
  endif

  ntet0 = 0
  allocate(tet(ntet0,6) )
  nprs0 = ntrias*(nr-1)              ! All over the hemisphere.
  allocate(prs(nprs0,6) )
  nhex0 = nquads_cyl*(nr-1)          ! All over the cylinder part
  allocate(hex(nhex0,8))

  write(*,*) "-- Mixed grid: Expected values... "
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    prisms     = ", nprs0
  write(*,*) "    tetrahedra = ", ntet0
  write(*,*) "    hex        = ", nhex0
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b0
  write(*,*) "         quads = ", nquads_b0

  nprs = 0
  ntet = 0
  nhex = 0
  ntrias_b = 0
  nquads_b = 0

!------------------------------------------------------------
! Copy the triangulation on the body
!------------------------------------------------------------

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate prisms

  do i = 1, ntrias

   node1 = tria(i)%v(1)
   node2 = tria(i)%v(2)
   node3 = tria(i)%v(3)
   node4 = node_above(tria(i)%v(1))
   node5 = node_above(tria(i)%v(2))
   node6 = node_above(tria(i)%v(3))

     nprs = nprs + 1

     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_body(i) = nprs
     cell_below   = nprs

   do k = 2, nr-1

    node1 = node4
    node2 = node5
    node3 = node6
    node4 = node_above(node1)
    node5 = node_above(node2)
    node6 = node_above(node3)

     nprs = nprs + 1
     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_above(cell_below) = nprs
     cell_below             = nprs

   end do

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,5) = 3
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

  end do

!------------------------------------------------------------
! Quads on the cylinder part
!------------------------------------------------------------

  nquads_b = 0

  do i = 1, nquads_cyl

     node1 = quad_cyl(i,1)
     node2 = quad_cyl(i,2)
     node3 = quad_cyl(i,3)
     node4 = quad_cyl(i,4)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 1      ! Quads on the cylinder surface
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node3
     quad(nquads_b,4) = node4

  end do

! Generate Hex

  do i = 1, nquads_cyl

     node1 = quad_cyl(i,1)
     node2 = quad_cyl(i,2)
     node3 = quad_cyl(i,3)
     node4 = quad_cyl(i,4)

   node5 = node_above(node1)
   node6 = node_above(node2)
   node7 = node_above(node3)
   node8 = node_above(node4)

     nhex = nhex + 1
     hex(nhex,1) = node1
     hex(nhex,2) = node2
     hex(nhex,3) = node3
     hex(nhex,4) = node4
     hex(nhex,5) = node5
     hex(nhex,6) = node6
     hex(nhex,7) = node7
     hex(nhex,8) = node8

     cell_body(ntrias+i) = nprs + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)
     cell_below          = nprs + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)

   do k = 2, nr-1

    node1 = node5
    node2 = node6
    node3 = node7
    node4 = node8
    node5 = node_above(node1)
    node6 = node_above(node2)
    node7 = node_above(node3)
    node8 = node_above(node4)

     nhex = nhex + 1
     hex(nhex,1) = node1
     hex(nhex,2) = node2
     hex(nhex,3) = node3
     hex(nhex,4) = node4
     hex(nhex,5) = node5
     hex(nhex,6) = node6
     hex(nhex,7) = node7
     hex(nhex,8) = node8

     cell_above(cell_below) = nprs + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)
     cell_below             = nprs + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)

   end do

    ! Quads on the outer boundary
    ! Reverse the ordering, so that it points inward to the interior.

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 3
     quad(nquads_b,1) = node8
     quad(nquads_b,2) = node7
     quad(nquads_b,3) = node6
     quad(nquads_b,4) = node5

  end do

! Outflow boundary

  if     (n_sections==6) then
   n_temp = nnodes_circum
  else !if (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

   do i = 1, n_temp

    node1 = nodes_cylinder(i  , nnodes_cylinder+1)
    node2 = nodes_cylinder(i+1, nnodes_cylinder+1)
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 2   ! Quads on the outflow boundary
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

      nquads_b = nquads_b + 1
      quad(nquads_b,5) = 2   ! Quads on the outflow boundary
      quad(nquads_b,1) = node1
      quad(nquads_b,2) = node2
      quad(nquads_b,3) = node4
      quad(nquads_b,4) = node3

    end do

   end do

 section1 : if (n_sections/=6) then

 ! Symmetry boundary 1: z=0 and positive y.

   do i = 1, nr_gs

    node1 = (i-1)*i/2 + 1
    node2 = (i+1)*i/2 + 1
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do


   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(1, i  ) !Left-most node of the original sector
    node2 = nodes_cylinder(1, i+1) !Left-most node of the original sector
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 ! Symmetry boundary 2: z=0 and negative y.

   do i = 1, nr_gs

    node1 = (i+1)*i/2 + i+1 !Right-most node of the original sector
    node2 = (i-1)*i/2 + i   !Right-most node of the original sector

   !Map onto the last sector:
    node1 = node_map((i+1)*i/2 + i +1)
    node2 = node_map((i-1)*i/2 + i)

    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
     if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(nnodes_circum, i+1)
    node2 = nodes_cylinder(nnodes_circum, i  )
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 endif section1

  write(*,*) " -----------------------------------------------------"
  write(*,*) " Generated elements....."
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    hex        = ", nhex
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*) " -----------------------------------------------------"
  write(*,*)
  write(*,*) " check(tetra ) = ", ntet0-ntet
  write(*,*) " check(prisms) = ", nprs0-nprs
  write(*,*) " check(hexa  ) = ", nhex0-nhex
  write(*,*) " check(tria  ) = ", ntrias_b0-ntrias_b
  write(*,*) " check(quads ) = ", nquads_b0-nquads_b
  write(*,*)

 end subroutine mixed_ph_grid
!********************************************************************************



!*******************************************************************************
!* Mixed grid Generation: Prism and Hex
!
!*******************************************************************************
 subroutine mixed_strct_grid
 implicit none
 integer(kd) :: node7, node8
 integer(kd) :: cell_below, ntet0, nhex0, nprs0, ntrias_b0, nquads_b0
!*******************************************************************************
! Generate mixed grid
!*******************************************************************************
  ntrias_b0 = 2*ntrias !All boundaries
  allocate(tri(ntrias_b0,5))

  if (n_sections==6) then
   nquads_b  = nnodes_circum * (nr-1)   ! All over outflow boundary
   nquads_b0 = nquads_b + 2*nquads_cyl  ! + All over the wing/hc and outer
   allocate(quad(nquads_b0,5))
  else !if (n_sections==1) then
   nquads_b  = nnodes_circum * (nr-1)   ! All over outflow boundary
   nquads_b0 = nquads_b + 2*nquads_cyl  ! + All over the HC surface and outer
   nquads_b0 = nquads_b0 + 2*(nr_gs+nnodes_cylinder)*(nr-1)
   allocate(quad(nquads_b0,5))
  endif

  nprs0 = ntrias*(nr-1)              ! All over the apex region
  allocate(prs(nprs0,6) )
  nhex0 = nquads_cyl*(nr-1)          ! All over the wing/hc part
  allocate(hex(nhex0,8))
  ntet0 = 0

  write(*,*) "-- Structured grid  Expected values... "
  write(*,*) "    nodes      = ", nnodes
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet0
  write(*,*) "    prisms     = ", nprs0
  write(*,*) "    hex        = ", nhex0
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b

  nprs = 0
  ntet = 0
  nhex = 0
  ntrias_b = 0
  nquads_b = 0

!------------------------------------------------------------
! Copy the triangulation on the body
!------------------------------------------------------------

  do i = 1, ntrias
          ntrias_b = ntrias_b + 1
   tri(ntrias_b,5) = 1
   tri(ntrias_b,1) = tria(i)%v(1)
   tri(ntrias_b,2) = tria(i)%v(2)
   tri(ntrias_b,3) = tria(i)%v(3)
  end do

! Generate prisms

  do i = 1, ntrias

   node1 = tria(i)%v(1)
   node2 = tria(i)%v(2)
   node3 = tria(i)%v(3)
   node4 = node_above(tria(i)%v(1))
   node5 = node_above(tria(i)%v(2))
   node6 = node_above(tria(i)%v(3))

     nprs = nprs + 1

     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_body(i) = nprs
     cell_below   = nprs

   do k = 2, nr-1

    node1 = node4
    node2 = node5
    node3 = node6
    node4 = node_above(node1)
    node5 = node_above(node2)
    node6 = node_above(node3)

     nprs = nprs + 1
     prs(nprs,1) = node1
     prs(nprs,2) = node2
     prs(nprs,3) = node3
     prs(nprs,4) = node4
     prs(nprs,5) = node5
     prs(nprs,6) = node6

     cell_above(cell_below) = nprs
     cell_below             = nprs

   end do

     ntrias_b = ntrias_b + 1
     tri(ntrias_b,5) = 3
     tri(ntrias_b,1) = node6
     tri(ntrias_b,2) = node5
     tri(ntrias_b,3) = node4

  end do

!------------------------------------------------------------
! Quads on the hemisphere-cylinder
!------------------------------------------------------------

  nquads_b = 0

  do i = 1, nquads_cyl

     node1 = quad_cyl(i,1)
     node2 = quad_cyl(i,2)
     node3 = quad_cyl(i,3)
     node4 = quad_cyl(i,4)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 1      ! Quads on the hemisphere-cylinder surface
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node3
     quad(nquads_b,4) = node4

  end do

! Generate Hex

  do i = 1, nquads_cyl

     node1 = quad_cyl(i,1)
     node2 = quad_cyl(i,2)
     node3 = quad_cyl(i,3)
     node4 = quad_cyl(i,4)

   node5 = node_above(node1)
   node6 = node_above(node2)
   node7 = node_above(node3)
   node8 = node_above(node4)

     nhex = nhex + 1
     hex(nhex,1) = node1
     hex(nhex,2) = node2
     hex(nhex,3) = node3
     hex(nhex,4) = node4
     hex(nhex,5) = node5
     hex(nhex,6) = node6
     hex(nhex,7) = node7
     hex(nhex,8) = node8

     cell_body(ntrias+i) = nprs0 + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)
     cell_below          = nprs0 + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)

   do k = 2, nr-1

    node1 = node5
    node2 = node6
    node3 = node7
    node4 = node8
    node5 = node_above(node1)
    node6 = node_above(node2)
    node7 = node_above(node3)
    node8 = node_above(node4)

     nhex = nhex + 1
     hex(nhex,1) = node1
     hex(nhex,2) = node2
     hex(nhex,3) = node3
     hex(nhex,4) = node4
     hex(nhex,5) = node5
     hex(nhex,6) = node6
     hex(nhex,7) = node7
     hex(nhex,8) = node8

    cell_above(cell_below) = nprs0 + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)
    cell_below             = nprs0 + nhex ! Cell # = 1,2,3,...,(ntet+nprs+nhex)

   end do

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 3      ! Quads on the outer boundary
     quad(nquads_b,1) = node8
     quad(nquads_b,2) = node7
     quad(nquads_b,3) = node6
     quad(nquads_b,4) = node5

  end do

! Outflow boundary

  if     (n_sections==6) then
   n_temp = nnodes_circum
  elseif (n_sections==1) then
   n_temp = nnodes_circum - 1
  endif

   do i = 1, n_temp

    node1 = nodes_cylinder(i  , nnodes_cylinder+1)
    node2 = nodes_cylinder(i+1, nnodes_cylinder+1)
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 2   ! Quads on the outflow boundary
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

      nquads_b = nquads_b + 1
      quad(nquads_b,5) = 2   ! Quads on the outflow boundary
      quad(nquads_b,1) = node1
      quad(nquads_b,2) = node2
      quad(nquads_b,3) = node4
      quad(nquads_b,4) = node3

    end do

   end do

 section1 : if (n_sections/=6) then

 ! Symmetry boundary 1: z=0 and positive y.

  !Left-most node in the prismatic region.

     node1 = 1        !apex node
     node2 = 2        !First node off the apex (left node)

   do k = 1, nr-1

     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

     node1 = node3
     node2 = node4

   end do

  !Left-most nodes in the hex region over hemisphere.
   do i = 2, nr_gs

    node1 = 2 + (i-2)*nnodes_circum
    node2 = node1 + nnodes_circum
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

  !Hex region over the cylinder.
   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(1, i  )
    node2 = nodes_cylinder(1, i+1)
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 ! Symmetry boundary 2: z=0 and negative y.

  !Right-most node in the prismatic region.

     node2 = 1                 !apex node
     node1 = 1 + nnodes_circum !Last node

   do k = 1, nr-1

     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

     node1 = node3
     node2 = node4

   end do

  !Right-most nodes in the hex region over hemisphere.
   do i = 2, nr_gs

    node2 = (1 + nnodes_circum) + (i-2)*nnodes_circum
    node1 = node2 + nnodes_circum
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 4
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

  !Hex region over the cylinder.
   do i = 1, nnodes_cylinder

    node1 = nodes_cylinder(nnodes_circum, i+1)
    node2 = nodes_cylinder(nnodes_circum, i  )
    node3 = node_above(node1)
    node4 = node_above(node2)

    nquads_b = nquads_b + 1
    quad(nquads_b,5) = 5
    if (n_sections==3) quad(nquads_b,5) = 4
    quad(nquads_b,1) = node1
    quad(nquads_b,2) = node2
    quad(nquads_b,3) = node4
    quad(nquads_b,4) = node3

    do k = 2, nr-1

     node1 = node3
     node2 = node4
     node3 = node_above(node1)
     node4 = node_above(node2)

     nquads_b = nquads_b + 1
     quad(nquads_b,5) = 5
     if (n_sections==3) quad(nquads_b,5) = 4
     quad(nquads_b,1) = node1
     quad(nquads_b,2) = node2
     quad(nquads_b,3) = node4
     quad(nquads_b,4) = node3

    end do

   end do

 endif section1

  write(*,*) " -----------------------------------------------------"
  write(*,*) " Generated elements....."
  write(*,*) "  Volume elements:"
  write(*,*) "    tetrahedra = ", ntet
  write(*,*) "    prisms     = ", nprs
  write(*,*) "    hex        = ", nhex
  write(*,*) "  Boundary elements:"
  write(*,*) "     triangles = ", ntrias_b
  write(*,*) "         quads = ", nquads_b
  write(*,*) " -----------------------------------------------------"
  write(*,*)
  write(*,*) " check(tetra ) = ", ntet0-ntet
  write(*,*) " check(prisms) = ", nprs0-nprs
  write(*,*) " check(hexa  ) = ", nhex0-nhex
  write(*,*) " check(tria  ) = ", ntrias_b0-ntrias_b
  write(*,*) " check(quads ) = ", nquads_b0-nquads_b
  write(*,*)

 end subroutine mixed_strct_grid
!********************************************************************************

!*******************************************************************************
!* Write out a file containing line info (only within the viscous layer)):
!* This information is for the implicit line relaxations or line agglomeration.
!* This format is used by FUN3D.
!*******************************************************************************
 subroutine line_info_file(n_points, file_name)
 implicit none

 integer(kd)      , intent(in) :: n_points
 character(80), intent(in) :: file_name

 integer(kd) :: node_below
 integer(kd) :: i, k, n_lines, n_total_points, max_points, min_points, i_count, inode,os

  open(unit=1, file=file_name, status="unknown", iostat=os)
  max_points = n_points
  min_points = n_points
  n_lines = nnodes_body
  n_total_points = n_lines * n_points

  write(*,*) " Total number of lines = ", n_lines
  write(*,*) "            max points = ", max_points
  write(*,*) "            min points = ", min_points

  write(1,*) n_lines, n_total_points, "Total lines and points"
  write(1,*) min_points, max_points, "Min and max points in line"

  i_count = 0

! Go up and write out the node info for every node on the body.
! NOTE: Go up to the given limit: n_points

  do i = 1, nnodes_body

     i_count = i_count + 1

!  1. First node (on the body)
     write(1,*) n_points, " Points in line for line = ", i_count
     write(1,'(i10,a10,3es30.20)') i, " x/y/z= ", node(i)%x, node(i)%y, node(i)%z
     node_below = i

!  2. Second node to the node below the last one.
   do k = 2, n_points-1
     write(1,'(i10)') node_above(node_below)
     node_below = node_above(node_below)
   end do

!  3. Last node
     inode = node_above(node_below)
     write(1,'(i10,a10,3es30.20)') inode, " x/y/z= ", node(inode)%x, node(inode)%y, node(inode)%z

  end do

 if (i_count /= n_lines) write(*,*) "Error: i_count /= n_lines"
 write(*,*)
 write(*,*) "lines_fmt file has been written: ", file_name

 close(1)
 end subroutine line_info_file
!********************************************************************************


!*******************************************************************************
!* Write out a file containing line info (Lines go all the way up to the outer):
!* This information is for the implicit line relaxations or line agglomeration.
!* This format is used by FUN3D.
!*******************************************************************************
 subroutine line_info_file_all(n_points, file_name)
 implicit none

 integer(kd)      , intent(in) :: n_points
 character(80), intent(in) :: file_name

 integer(kd) :: node_below
 integer(kd) :: i, k, n_lines, n_total_points, max_points, min_points, i_count, inode,os

  open(unit=1, file=file_name, status="unknown", iostat=os)
  max_points = n_points
  min_points = n_points
  n_lines = nnodes_body
  n_total_points = n_lines * n_points

  write(*,*) " Total number of lines = ", n_lines
  write(*,*) "            max points = ", max_points
  write(*,*) "            min points = ", min_points

  write(1,*) n_lines, n_total_points, "Total lines and points"
  write(1,*) min_points, max_points, "Min and max points in line"

  i_count = 0

! Go up and write out the node info for every node on the body.
! NOTE: Go up to the given limit: n_points

  do i = 1, nnodes_body

     i_count = i_count + 1

!  1. First node (on the body)
     write(1,*) n_points, " Points in line for line = ", i_count
     write(1,'(i10,a10,3es30.20)') i, " x/y/z= ", node(i)%x, node(i)%y, node(i)%z
     node_below = i

!  2. Second node to the node below the last one.
   do k = 2, n_points-1
     write(1,'(i10)') node_above(node_below)
     node_below = node_above(node_below)
   end do

!  3. Last node
     inode = node_above(node_below)
     write(1,'(i10,a10,3es30.20)') inode, " x/y/z= ", node(inode)%x, &
                                       node(inode)%y, node(inode)%z

  end do

 if (i_count /= n_lines) write(*,*) "Error: i_count /= n_lines"
 write(*,*)
 write(*,*) "lines_fmt file has been written: ", file_name

 close(1)
 end subroutine line_info_file_all
!********************************************************************************


!*******************************************************************************
! This subroutine writes a Tecplot file for the volume grid.
!*******************************************************************************
 subroutine write_tecplot_volume_file

 open(unit=8, file=filename_tecplot_v, status="unknown", iostat=os)
 write(8,*) 'TITLE = "GRID"'
 write(8,*) 'VARIABLES = "x","y","z","k1","k2","k3","k4","k5"'

! Tetra Zone
  if (ntet > 0) then

   write(8,*) 'zone  n=', nnodes,',e=', ntet,' , et=tetrahedron, f=fepoint'
   do i = 1, nnodes
     write(8,'(3es20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, ntet
    write(8,'(4i10)') tet(i,1), tet(i,2), tet(i,3), tet(i,4)
   end do

  endif

! Prism zone
  if (nprs > 0) then

   write(8,*) 'zone  n=', nnodes,',e=', nprs,' , et=brick, f=fepoint'
   do i = 1, nnodes
     write(8,'(3es20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nprs
    write(8,'(8i10)') prs(i,1), prs(i,2), prs(i,3), prs(i,3), &
                      prs(i,4), prs(i,5), prs(i,6), prs(i,6)
   end do

  endif

! Hex zone
  if (nhex > 0) then

   write(8,*) 'zone  n=', nnodes,',e=', nhex,' , et=brick, f=fepoint'
   do i = 1, nnodes
     write(8,'(3es20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nhex
    write(8,'(8i10)') hex(i,1), hex(i,2), hex(i,3), hex(i,4), &
                      hex(i,5), hex(i,6), hex(i,7), hex(i,8)
   end do

  endif

 close(8)

 end subroutine write_tecplot_volume_file
!********************************************************************************


!*******************************************************************************
!* CC line data:
!*
!* Write out a file containing line info (only within the viscous layer)):
!* This information is for the implicit line relaxations or line agglomeration.
!*
!*******************************************************************************
 subroutine line_info_file_c(n_cells, file_name)
 implicit none

 integer(kd)      , intent(in) :: n_cells
 character(80), intent(in) :: file_name

 integer(kd) :: cell_below, icell,os
 integer(kd) :: i, k, n_lines, n_total_cells, max_cells, min_cells

  open(unit=1, file=file_name, status="unknown", iostat=os)
  max_cells     = n_cells
  min_cells     = n_cells
  n_lines       = ntrias + nquads_cyl
  n_total_cells = n_lines * n_cells

  write(*,*)
  write(*,*) " Cell-line information: "
  write(*,*) " Total number of lines = ", n_lines
  write(*,*) "            max cells = ", max_cells
  write(*,*) "            min cells = ", min_cells

!  write(1,*) n_lines, n_total_cells, "Total lines and cells"
!  write(1,*) min_cells, max_cells, "Min and max cells in line"

! Go up and write out the node info for every node on the body.
! NOTE: Go up to the given limit: n_cells

  write(1,*) n_lines, n_total_cells, "Total lines and cells"
  write(1,*)   min_cells, max_cells, "Min and max cells in line"

  do i = 1, ntrias + nquads_cyl

!  1. First cell (on the body)
     write(1,*) n_cells, " Cells in line for line = ", i
     icell = cell_body(i)
     write(1,'(i12,a10,3es30.20)') icell, " x/y/z= ", cell(icell)%x, &
                                       cell(icell)%y, cell(icell)%z

     cell_below = icell

!  2. Second cell to the cell below the last one.
   do k = 2, n_cells-1
     write(1,'(i12)') cell_above(cell_below)
         cell_below = cell_above(cell_below)
   end do

!  3. Last cell
     if (n_cells > 1) then
      icell = cell_above(cell_below)
      write(1,'(i12,a10,3es30.20)') icell, " x/y/z= ", cell(icell)%x, &
                                        cell(icell)%y, cell(icell)%z
     endif

  end do

 write(*,*)
 write(*,*) "lines_cc_fmt file has been written: ", file_name

 close(1)
 end subroutine line_info_file_c
!********************************************************************************




!*******************************************************************************
! This subroutine writes  a Tecplot file for boundaries.
!******************************************************************************
 subroutine write_tecplot_boundary_file

 integer(kd) :: nelms
 integer(kd)                               :: nnodes_loc, i, j, k, i_boundary
 integer(kd),  dimension(:,:), allocatable :: g2l
 integer(kd),  dimension(:  ), allocatable :: l2g

 open(unit=7, file=filename_tecplot_b, status="unknown", iostat=os)
 write(7,*) 'TITLE = "GRID"'
 write(7,*) 'VARIABLES = "x","y","z","k1","k2","k3","k4","k5"'

 allocate(g2l(nnodes,2))
 allocate(l2g(nnodes  ))

!--------------------------------------------------------------
! Triangles on the hemisphere-cylinder surface: i_boundary = 1.

  i_boundary = 1

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, ntrias_b
   if ( tri(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 3
     if (g2l(tri(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(tri(i,k),1) = nnodes_loc
      g2l(tri(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = tri(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="HC body: triangles"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, ntrias_b
    if ( tri(i,5) == i_boundary ) write(7,'(4I10)') g2l(tri(i,1),1), g2l(tri(i,2),1), &
                                                    g2l(tri(i,3),1), g2l(tri(i,3),1)
   end do
  endif

!--------------------------------------------------------------
! Triangles on the outflow plane: i_boundary = 2.

  i_boundary = 2

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, ntrias_b
   if ( tri(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 3
     if (g2l(tri(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(tri(i,k),1) = nnodes_loc
      g2l(tri(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = tri(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Outflow: triangles"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, ntrias_b
    if ( tri(i,5) == i_boundary ) write(7,'(4I10)') g2l(tri(i,1),1), g2l(tri(i,2),1), &
                                                    g2l(tri(i,3),1), g2l(tri(i,3),1)
   end do
  endif

!--------------------------------------------------------------
! Triangles on the outer boundary: i_boundary = 3.

  i_boundary = 3

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, ntrias_b
   if ( tri(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 3
     if (g2l(tri(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(tri(i,k),1) = nnodes_loc
      g2l(tri(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = tri(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Farfield: triangles"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, ntrias_b
    if ( tri(i,5) == i_boundary ) write(7,'(4I10)') g2l(tri(i,1),1), g2l(tri(i,2),1), &
                                                    g2l(tri(i,3),1), g2l(tri(i,3),1)
   end do
  endif

!--------------------------------------------------------------
! Triangles on the symmetry plane 1: i_boundary = 4.

  i_boundary = 4

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, ntrias_b
   if ( tri(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 3
     if (g2l(tri(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(tri(i,k),1) = nnodes_loc
      g2l(tri(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = tri(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Symmetry 1: triangles"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, ntrias_b
    if ( tri(i,5) == i_boundary ) write(7,'(4I10)') g2l(tri(i,1),1), g2l(tri(i,2),1), &
                                                    g2l(tri(i,3),1), g2l(tri(i,3),1)
   end do
  endif

!--------------------------------------------------------------
! Triangles on the symmetry plane 2: i_boundary = 5.

  i_boundary = 5

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, ntrias_b
   if ( tri(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 3
     if (g2l(tri(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(tri(i,k),1) = nnodes_loc
      g2l(tri(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = tri(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Symmetry 2: triangles"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, ntrias_b
    if ( tri(i,5) == i_boundary ) write(7,'(4I10)') g2l(tri(i,1),1), g2l(tri(i,2),1), &
                                                    g2l(tri(i,3),1), g2l(tri(i,3),1)
   end do
  endif

!--------------------------------------------------------------
! Quad elements
!--------------------------------------------------------------

 !--------------------------------------
 !Body: i_boundary = 1.

  i_boundary = 1

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, nquads_b
   if ( quad(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 4
     if (g2l(quad(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(quad(i,k),1) = nnodes_loc
      g2l(quad(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = quad(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="HC body: quads"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nquads_b
    if ( quad(i,5) == i_boundary ) write(7,'(4I10)') g2l(quad(i,1),1), g2l(quad(i,2),1), &
                                                     g2l(quad(i,3),1), g2l(quad(i,4),1)
   end do
  endif

 !--------------------------------------------------------------
 !Outflow boundary: i_boundary = 2.

  i_boundary = 2

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, nquads_b
   if ( quad(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 4
     if (g2l(quad(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(quad(i,k),1) = nnodes_loc
      g2l(quad(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = quad(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Outflow: quads"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nquads_b
    if ( quad(i,5) == i_boundary ) write(7,'(4I10)') g2l(quad(i,1),1), g2l(quad(i,2),1), &
                                                     g2l(quad(i,3),1), g2l(quad(i,4),1)
   end do
  endif

 !--------------------------------------------------------------
 !Outer boundary: i_boundary = 3.

  i_boundary = 3

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, nquads_b
   if ( quad(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 4
     if (g2l(quad(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(quad(i,k),1) = nnodes_loc
      g2l(quad(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = quad(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Farfield: quads"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nquads_b
    if ( quad(i,5) == i_boundary ) write(7,'(4I10)') g2l(quad(i,1),1), g2l(quad(i,2),1), &
                                                     g2l(quad(i,3),1), g2l(quad(i,4),1)
   end do
  endif

 !--------------------------------------------------------------
 !Symmetry plane 1: i_boundary = 4

  i_boundary = 4

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, nquads_b
   if ( quad(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 4
     if (g2l(quad(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(quad(i,k),1) = nnodes_loc
      g2l(quad(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = quad(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Symmetry 1: quads"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nquads_b
    if ( quad(i,5) == i_boundary ) write(7,'(4I10)') g2l(quad(i,1),1), g2l(quad(i,2),1), &
                                                     g2l(quad(i,3),1), g2l(quad(i,4),1)
   end do
  endif

 !--------------------------------------------------------------
 !Symmetry plane 1: i_boundary = 5

  i_boundary = 5

       nelms = 0
  nnodes_loc = 0
         g2l = -1
  do i = 1, nquads_b
   if ( quad(i,5) == i_boundary ) then
    nelms = nelms + 1
    do k = 1, 4
     if (g2l(quad(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(quad(i,k),1) = nnodes_loc
      g2l(quad(i,k),2) = 100 !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = quad(i,k)
     endif
    end do
   endif
  end do

  if (nelms > 0) then
   write(7,*) 'ZONE T="Symmetry 2: quads"  N=', nnodes_loc,',E=', nelms,&
              ' , ET=quadrilateral, F=FEPOINT'

   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(3ES20.10,5i13)') node(i)%x, node(i)%y, node(i)%z, &
                                k1(i),k2(i),k3(i),k4(i),k5(i)
   end do

   do i = 1, nquads_b
    if ( quad(i,5) == i_boundary ) write(7,'(4I10)') g2l(quad(i,1),1), g2l(quad(i,2),1), &
                                                     g2l(quad(i,3),1), g2l(quad(i,4),1)
   end do
  endif

  !--------------------------------------

 close(7)

 end subroutine write_tecplot_boundary_file
!********************************************************************************

!*******************************************************************************
! This subroutine writes a ugrid file.
!*******************************************************************************
 subroutine write_ugrid_file

  if ( ugrid_file_unformatted ) then
    open(unit=9, file=filename_ugrid, form='unformatted',access="stream",&
                                      status='unknown', iostat=os )
    write(9) nnodes,   ntrias_b,    nquads_b,   ntet,    0, nprs, nhex
  else
    open(unit=9, file=filename_ugrid, status="unknown", iostat=os)
    !                    #nodes, #tri_faces, #quad_faces, #tetra, #pyr, #prz,
    !                    #hex
    write(9,'(7I20)') nnodes,   ntrias_b,    nquads_b,   ntet,    0, nprs, nhex
  endif

!---------------------------------------------------------------
!(1) Unformatted grid file
!---------------------------------------------------------------

  if ( ugrid_file_unformatted ) then

  ! Nodes
    do i = 1, nnodes
     write(9) node(i)%x, node(i)%y, node(i)%z
    end do

  ! Triangular faces = ntri
    if (ntrias_b > 0) then
     do i = 1, ntrias_b
      write(9) tri(i,1), tri(i,2), tri(i,3)
     end do
    endif

  ! Quad faces = nquad
    if (nquads_b > 0) then
     do i = 1, nquads_b
      write(9) quad(i,1), quad(i,2), quad(i,3), quad(i,4)
     end do
    endif

  ! Face tag
    if (ntrias_b > 0) then
     do i = 1, ntrias_b
      write(9)  tri(i,5)
     end do
    endif

    if (nquads_b > 0) then
     do i = 1, nquads_b
      write(9) quad(i,5)
     end do
    endif

  ! tet
    if (ntet > 0) then
     do i = 1, ntet
      write(9) tet(i,1), tet(i,2), tet(i,3), tet(i,4)
     end do
    endif

  ! Prism
    if (nprs > 0) then
     do i = 1, nprs
      write(9) prs(i,1), prs(i,2), prs(i,3), &
               prs(i,4), prs(i,5), prs(i,6)
     end do
    endif

  ! Hex
    if (nhex > 0) then
     do i = 1, nhex
      write(9) hex(i,1), hex(i,2), hex(i,3), &
               hex(i,4), hex(i,5), hex(i,6), hex(i,7), hex(i,8)
     end do
    endif

!---------------------------------------------------------------
!(2) UFormatted grid file
!---------------------------------------------------------------

  else

  ! Nodes
    do i = 1, nnodes
     write(9,'(3ES26.15)') node(i)%x, node(i)%y, node(i)%z
    end do

  ! Triangular faces = ntri
    if (ntrias_b > 0) then
     do i = 1, ntrias_b
      write(9,'(3I20)') tri(i,1), tri(i,2), tri(i,3)
     end do
    endif

  ! Quad faces = nquad
    if (nquads_b > 0) then
     do i = 1, nquads_b
      write(9,'(4I20)') quad(i,1), quad(i,2), quad(i,3), quad(i,4)
     end do
    endif

  ! Face tag
    if (ntrias_b > 0) then
     do i = 1, ntrias_b
      write(9,'(I110)')  tri(i,5)
     end do
    endif

    if (nquads_b > 0) then
     do i = 1, nquads_b
      write(9,'(I110)') quad(i,5)
     end do
    endif

  ! tet
    if (ntet > 0) then
     do i = 1, ntet
      write(9,'(4I20)') tet(i,1), tet(i,2), tet(i,3), tet(i,4)
     end do
    endif

  ! Prism
    if (nprs > 0) then
     do i = 1, nprs
      write(9,'(6I20)') prs(i,1), prs(i,2), prs(i,3), &
                        prs(i,4), prs(i,5), prs(i,6)
     end do
    endif

  ! Hex
    if (nhex > 0) then
     do i = 1, nhex
      write(9,'(8I20)') hex(i,1), hex(i,2), hex(i,3), &
                        hex(i,4), hex(i,5), hex(i,6), hex(i,7), hex(i,8)
     end do
    endif

  endif

!---------------------------------------------------------------
! End of Unformatted or Formatted
!---------------------------------------------------------------

  close(9)

 end subroutine write_ugrid_file
!********************************************************************************

!*******************************************************************************
! This subroutine writes a PLOT3D file (for structured grids only: igrid_type=5).
!*******************************************************************************
 subroutine write_plot3d_file

  integer(kd) :: i, j, k
  integer(kd) :: k1min, k4min, k5min
  integer(kd) :: k1max, k4max, k5max
  logical :: full_gmtry

  integer(kd) :: imax, jmax, kmax

   if (n_sections == 6) then
    full_gmtry = .true.
   else
    full_gmtry = .false.
   endif

  !Original data

   k1min = minval(k1)
   k4min = minval(k4)
   k5min = minval(k5)

   k1max = maxval(k1)
   k4max = maxval(k4)
   k5max = maxval(k5)

  !Original data: start from 0, k1=0 is used for the apex node only.

   write(*,*)
   write(*,*) " --- Original:"
   write(*,*) "  i = ", k4min, ",", k4max
   write(*,*) "  j = ", k5min, ",", k5max
   write(*,*) "  k = ", k1min, ",", k1max
   write(*,*)

    allocate( k2n_temp(k4min:k4max, k5min:k5max, k1min:k1max) )

   do i = 1, nnodes
    k2n_temp( k4(i), k5(i), k1(i) ) = i
   end do

  !Modified data: start from 1, k1=0 is distributed over k1=1,k1max.

   k4max = k4max + 1
   k5max = k5max + 1
   k1max = k1max

   if (full_gmtry) then
    k1max = k1max + 1
   endif

   allocate( k2n(k4max, k5max, k1max) )

    k2n = -1

   do i = 1, nnodes
    if (k1(i)==0) cycle !<- Exclude the apex node.
    k2n( k4(i)+1, k5(i)+1, k1(i) ) = k2n_temp( k4(i), k5(i), k1(i) )
   end do

   do k = 1, k1max
     k2n( 1, :, k ) = k2n_temp( 0, 0:k5max-1, 0 ) !Nodes above the apex.
   end do

   if (full_gmtry) then
    !Copy the first plane to close.
     k2n( :, :, k1max ) = k2n( :, :, 1 )
   endif

   write(*,*)
   write(*,*) " --- Modified for PLOT3D:"
   write(*,*) " i = ", 1, ",", k4max
   write(*,*) " j = ", 1, ",", k5max
   write(*,*) " k = ", 1, ",", k1max
   write(*,*)

   deallocate( k2n_temp )

  if ( ugrid_file_unformatted ) then
    !Note: access="stream" removed (04-04-2017). It was causing some problems...
    open(unit=9, file=filename_p3d, form='unformatted', status='unknown', iostat=os )
    write(9) 1 ! Single block
    write(9) k4max, k5max, k1max ! ni, nj, nk
  else
    open(unit=9, file=filename_p3d, status="unknown", iostat=os)
    write(9,'(I20)') 1 ! Single block
    write(9,'(3I20)') k4max, k5max, k1max ! ni, nj, nk
  endif

   if ( ugrid_file_unformatted ) then
    write(9)                                                            &
         ((( node( k2n(i,j,k) )%x , i=1,k4max), j=1,k5max), k=1,k1max), &
         ((( node( k2n(i,j,k) )%y , i=1,k4max), j=1,k5max), k=1,k1max), &
         ((( node( k2n(i,j,k) )%z , i=1,k4max), j=1,k5max), k=1,k1max)
   else
    write(9,'(3ES26.15)')                                               &
         ((( node( k2n(i,j,k) )%x , i=1,k4max), j=1,k5max), k=1,k1max), &
         ((( node( k2n(i,j,k) )%y , i=1,k4max), j=1,k5max), k=1,k1max), &
         ((( node( k2n(i,j,k) )%z , i=1,k4max), j=1,k5max), k=1,k1max) 
   endif

  close(9)


  deallocate( k2n )

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
! Write .nmf file.
!
! 1 -> kmin,   i, j
! 2 -> kmax,   i, j
!
! 3 -> imin,   j, k
! 4 -> imax,   j, k
!
! 5 -> jmin,   k, i
! 6 -> jmax,   k, i
!----------------------------------------------------------------------------------

! imin = 1
  imax = k4max
! jmin = 1
  jmax = k5max
! kmin = 1
  kmax = k1max

 open(unit=9, file=filename_nmf, status='unknown', iostat=os )

 write(9,'(a)') "# ============= NASA Langley Geometry Laboratory TOG Neutral Map File ==============="
 write(9,'(a)') "# ==================================================================================="
 write(9,'(a)') "# Block#   IDIM   JDIM   KDIM"
 write(9,'(a)') "# ==================================================================================="
 write(9,'(i3)')     1
 write(9,*)
 write(9,'(4i15)')   1,  imax, jmax, kmax
 write(9,*)
 write(9,'(a)') "# ==================================================================================="
 write(9,'(a)') "# Type         B1  F1     S1   E1     S2   E2    B2  F2     S1   E1     S2   E2  Swap"
 write(9,'(a)') "#------------------------------------------------------------------------------------"
 write(9,'(a,2i3,i3,i15,i3,i15,i3)') "         pole ", 1, 3, 1, jmax, 1, kmax, 2
 write(9,'(a,2i3,i3,i15,i3,i15)')    "back_pressure ", 1, 4, 1, jmax, 1, kmax
 write(9,'(a,2i3,i3,i15,i3,i15)')    "viscous_solid ", 1, 5, 1, kmax, 1, imax
 write(9,'(a,2i3,i3,i15,i3,i15)')    "     farfield ", 1, 6, 1, kmax, 1, imax

 !(1)Full geometry: Just to tell two boundaries match and form an interior plane.
 if (full_gmtry) then
  write(9,'(a,2i3,i3,i15,i3,i15, 2i3, i3,i15,i3,i15, a)') "   one_to_one "      , &
                                                          1, 1, 1, imax, 1, jmax, &
                                                          1, 2, 1, imax, 1, jmax, "  FALSE"
 !(2)Half geometry -> Two symmetry planes (one for positive z; the other for negative z).
 else
  write(9,'(a,2i3,i3,i15,i3,i15)')    "   symmetry_y ", 1, 1, 1, imax, 1, jmax 
  write(9,'(a,2i3,i3,i15,i3,i15)')    "   symmetry_y ", 1, 2, 1, imax, 1, jmax 
 endif

 close(9)
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------


 end subroutine write_plot3d_file
!********************************************************************************



!*******************************************************************************
! This subroutine writes an k file.
! Note: All grids are structured in the k-coordinate system.
!******************************************************************************
 subroutine write_k_file

  integer(kd) :: k1min, k2min, k3min, k4min, k5min
  integer(kd) :: k1max, k2max, k3max, k4max, k5max

 !----------------------------------------------------------------------------------
  if ( k_file_unformatted ) then
    open(unit=10, file=filename_k, form='unformatted', status='unknown', iostat=os )
 !----------------------------------------------------------------------------------
  else
    open(unit=10, file=filename_k, status="unknown", iostat=os)
  endif
 !----------------------------------------------------------------------------------

 if ( igrid_type == 5 .and. n_sections == 3) then
  k3(1) = -1 !To indicate that this is a half-geometry structured grid.
  write(*,*)
  write(*,*) " Set k3(1) = -1 to indicate that this is a half-geometry grid."
 endif

! Write node_number, k1, k2, k3, k4, k5

 !----------------------------------------------------------------------------------
  if ( k_file_unformatted ) then
   write(10) nnodes
   do i = 1, nnodes
    write(10) i, k1(i),k2(i),k3(i),k4(i),k5(i)
   end do
 !----------------------------------------------------------------------------------
  else
   write(10,*) nnodes
   do i = 1, nnodes
    write(10,'(6i13)') i, k1(i),k2(i),k3(i),k4(i),k5(i)
   end do
  endif
 !----------------------------------------------------------------------------------

 close(10)

   k1min = minval(k1)
   k2min = minval(k2)
   k3min = minval(k3)
   k4min = minval(k4)
   k5min = minval(k5)

   k1max = maxval(k1)
   k2max = maxval(k2)
   k3max = maxval(k3)
   k4max = maxval(k4)
   k5max = maxval(k5)

   write(*,*)


  if     ( igrid_type < 5 ) then

   write(*,*) "k1 min max = ", k1min, k1max, " k1-coord of the hemisphere grid points."
   write(*,*) "k2 min max = ", k2min, k2max, " k2-coord of the hemisphere grid points."
   write(*,*) "k3 min max = ", k3min, k3max, " k3 = -(k1+k2)"
   write(*,*) "k4 min max = ", k4min, k4max, " Index along the cylinder (=0 on hemisphere)."
   write(*,*) "k5 min max = ", k5min, k5max, " Index in the radial direction"
   write(*,*)

  elseif ( igrid_type == 5 ) then

   write(*,*) "k1 min max = ", k1min, k1max, " Around HC."
   write(*,*) "k2 min max = ", k2min, k2max, " Not used in structured grid."
   write(*,*) "k3 min max = ", k3min, k3max, " Used to indicate half-geometry (-1)."
   write(*,*) "k4 min max = ", k4min, k4max, " Index along the HC (apex to base)."
   write(*,*) "k5 min max = ", k5min, k5max, " Index in the radial direction"
   write(*,*)

  endif

   if ( igrid_type == 5 .and. n_sections == 3) then
    write(*,*) " NOTE: k1max must be odd for regular coarsening in the"
    write(*,*) "       case of a structured grid for a half geometry."
   endif

   write(*,*)

 end subroutine write_k_file
!********************************************************************************


!********************************************************************************
!* This subroutine is useful to expand or shrink type(node_data_yz) arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*                     expanded to the requested dimension, n, if n > dim(x).
!*                     shrunk to the requested dimension, n, if n < dim(x).
!*
!********************************************************************************
  subroutine my_alloc_ndyz_ptr(x,n)
  implicit none

  integer(kd),                    intent(in   ) :: n
  type(node_data_yz), dimension(:), pointer :: x

  integer(kd) :: i
  type(node_data_yz), dimension(:), pointer :: temp

  if (n <= 0) then
   write(*,*) "my_alloc_ndyz_ptr received non-positive dimension. Stop."
   stop
  endif

! If initial, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n))
    temp(n)%gnode = 0
    temp(n)%y     = zero
    temp(n)%z     = zero

! (1) Expand the array dimension
  if ( n > size(x) ) then

   do i = 1, size(x) 
    temp(i)%gnode = x(i)%gnode
    temp(i)%y     = x(i)%y
    temp(i)%z     = x(i)%z
  end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
  else

   do i = 1, n
    temp(i)%gnode = x(i)%gnode
    temp(i)%y     = x(i)%y
    temp(i)%z     = x(i)%z
   end do

  endif

! Destroy the target of x
  deallocate(x)

! Re-assign the pointer
   x => temp

  return

  end subroutine my_alloc_ndyz_ptr
!********************************************************************************

!********************************************************************************
! Just to get the sign of an integer...
!********************************************************************************
 function int_sign(n)
 implicit none
 integer(kd), intent(in) :: n
 integer(kd) :: int_sign

  if (n >= 0) then
   int_sign = 1
  else
   int_sign = -1
  endif

 end function int_sign
!********************************************************************************

!********************************************************************************
! Find out big_endian_io.
!********************************************************************************
 function big_endian_io( opt_unit )

 integer, intent(in) :: opt_unit
 logical             :: big_endian_io
! one-byte integer
 integer, parameter :: i1 = selected_int_kind(2)
! two-byte integer
 integer, parameter :: i2 = selected_int_kind(4)
 integer(i1) :: byte_one, byte_two
! 00000000 00000001 big-endian binary
 integer(i2) :: two_byte_int = 1_i2

    open(opt_unit,status='scratch',form='unformatted')
      write( opt_unit) two_byte_int
      rewind(opt_unit)
      read(  opt_unit) byte_one, byte_two
    close(opt_unit)
    big_endian_io = ( byte_one == 0 .and. byte_two == 1 )

 end function big_endian_io
!********************************************************************************

!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************
!********************************************************************************


!********************************************************************************
!* Check the volume of all types of elements.
!*
!********************************************************************************
  subroutine volume_check(negative_volume)
  implicit none

  logical, dimension(3), intent(out) :: negative_volume

  real(dp) :: x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4
  real(dp) :: x5,x6,x7,x8, y5,y6,y7,y8, z5,z6,z7,z8
  real(dp) :: volume, vol_max, vol_min, vol_ave

  integer(kd)               :: i, n_negative_vol, itet, ihex, iprs
  integer(kd), dimension(8) :: v
  character(80)         :: i_char, filename

  itet = 1
  ihex = 2
  iprs = 3
  negative_volume = .false.

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! (1)Check Tetrahedral elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 tet_vol_check : if (ntet > 0) then

   n_negative_vol =  0
          vol_max = -1.0_dp
          vol_min =  1.0e+15_dp
          vol_ave =  0.0_dp

  !-------------------------------------------------------------
  ! Loop over tetrahedra and check the volumne one by one.
  !-------------------------------------------------------------
   do i = 1, ntet

    x1 = node( tet(i,1) )%x
    y1 = node( tet(i,1) )%y
    z1 = node( tet(i,1) )%z

    x2 = node( tet(i,2) )%x
    y2 = node( tet(i,2) )%y
    z2 = node( tet(i,2) )%z

    x3 = node( tet(i,3) )%x
    y3 = node( tet(i,3) )%y
    z3 = node( tet(i,3) )%z

    x4 = node( tet(i,4) )%x
    y4 = node( tet(i,4) )%y
    z4 = node( tet(i,4) )%z

     volume = tet_volume(x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4)

    vol_max = max(vol_max,volume)
    vol_min = min(vol_min,volume)
    vol_ave = vol_ave + volume

    if (volume < 0.0_dp) then

     write(8933,*) " Negative volume tet at element = ", i

     !Write a Tecplot file for the negative-volume tetrahedra,
     !only for the first 20 tetra (to avoid too many files).
     if (n_negative_vol < 21) then

      write( i_char  , '(i0)' ) i
      filename = "negative_vol_tet_" // trim(i_char) // "_tec.dat"
      open(unit=3, file=filename, status="unknown", iostat=os)
      write(3,*) 'TITLE = ', "negative_vol_tet_" // trim(i_char)
      write(3,*) 'VARIABLES = "x","y","z"'
      write(3,*) 'ZONE  N=', 4,',E=', 1,' , et=tetrahedron, F=FEPOINT'
      write(3,'(3ES20.10)') x1, y1, z1
      write(3,'(3ES20.10)') x2, y2, z2
      write(3,'(3ES20.10)') x3, y3, z3
      write(3,'(3ES20.10)') x4, y4, z4
      write(3,'(4i10)'    )  1, 2, 3, 4
      close(3)

     endif

     n_negative_vol = n_negative_vol + 1

     ! Display a WARNING MESSAGE
     if (n_negative_vol == 1) then
      write(*,*)
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) " Negative tet volume deteceted. See fort.8933..."
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*)
     endif

    endif

   end do
  !-------------------------------------------------------------
  ! End of loop over tetrahedra
  !-------------------------------------------------------------

   vol_ave = vol_ave / real(ntet,dp)

   if (n_negative_vol == 0) then
     write(*,*)
     write(*,*) " Negative tetra volume not detected."
     write(*,*)
   endif

   write(*,*)
   write(*,'(a,es26.16)') " Minimum tetrahedral volume = ", vol_min
   write(*,'(a,es26.16)') " Maximum tetrahedral volume = ", vol_max
   write(*,'(a,es26.16)') " Average tetrahedral volume = ", vol_ave
   write(*,*)

   if (n_negative_vol > 0) then
     write(*,*) " Tetrahedra with negative volume created....."
     write(*,*) "  # of tets having negative volume = ", n_negative_vol
     negative_volume(itet) = .true.
   endif

 endif tet_vol_check
!--------------------------------------------------------------------
!--------------------------------------------------------------------
! End fo (1)Check Tetrahedral elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! (2)Check Hexahedral elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 hex_vol_check : if (nhex > 0) then

   n_negative_vol = 0
   vol_max = -1.0_dp
   vol_min =  1.0e+15_dp
   vol_ave =  0.0_dp

  !-------------------------------------------------------------
  ! Loop over hex and check the volumne one by one.
  !-------------------------------------------------------------
   do i = 1, nhex

    v(1:8) =  hex(i,1:8)

    x1 = node( v(1) )%x
    y1 = node( v(1) )%y
    z1 = node( v(1) )%z

    x2 = node( v(2) )%x
    y2 = node( v(2) )%y
    z2 = node( v(2) )%z

    x3 = node( v(3) )%x
    y3 = node( v(3) )%y
    z3 = node( v(3) )%z

    x4 = node( v(4) )%x
    y4 = node( v(4) )%y
    z4 = node( v(4) )%z

    x5 = node( v(5) )%x
    y5 = node( v(5) )%y
    z5 = node( v(5) )%z

    x6 = node( v(6) )%x
    y6 = node( v(6) )%y
    z6 = node( v(6) )%z

    x7 = node( v(7) )%x
    y7 = node( v(7) )%y
    z7 = node( v(7) )%z

    x8 = node( v(8) )%x
    y8 = node( v(8) )%y
    z8 = node( v(8) )%z

  !Cut into 6 tetrahedra and compute the volume of hex:

     volume = 0.0_dp
     volume = volume + tet_volume( x1, x4, x5, x6, &
                                   y1, y4, y5, y6, &
                                   z1, z4, z5, z6  )

     volume = volume + tet_volume( x4, x8, x5, x6, &
                                   y4, y8, y5, y6, &
                                   z4, z8, z5, z6  )

     volume = volume + tet_volume( x3, x2, x6, x4, &
                                   y3, y2, y6, y4, &
                                   z3, z2, z6, z4  )

     volume = volume + tet_volume( x3, x6, x7, x4, &
                                   y3, y6, y7, y4, &
                                   z3, z6, z7, z4  )

     volume = volume + tet_volume( x1, x2, x4, x6, &
                                   y1, y2, y4, y6, &
                                   z1, z2, z4, z6  )

     volume = volume + tet_volume( x4, x7, x8, x6, &
                                   y4, y7, y8, y6, &
                                   z4, z7, z8, z6  )

    vol_max = max(vol_max,volume)
    vol_min = min(vol_min,volume)
    vol_ave = vol_ave + volume

    if (volume < 0.0_dp) then

     write(8933,*) " Negative volume hex at element = ", i

     !Write a Tecplot file for the negative-volume tetrahedra,
     !only for the first 20 tetra (to avoid too many files).
     if (n_negative_vol < 21) then

      write( i_char  , '(i0)' ) i
      filename = "negative_vol_hex_" // trim(i_char) // "_tec.dat"
      open(unit=3, file=filename, status="unknown", iostat=os)
      write(3,*) 'TITLE = ', "negative_vol_tet_" // trim(i_char)
      write(3,*) 'VARIABLES = "x","y","z"'
      write(3,*) 'ZONE  N=', 8,',E=', 1,' , et=tetrahedron, F=FEPOINT'
      write(3,'(3ES20.10)') x1, y1, z1
      write(3,'(3ES20.10)') x2, y2, z2
      write(3,'(3ES20.10)') x3, y3, z3
      write(3,'(3ES20.10)') x4, y4, z4
      write(3,'(3ES20.10)') x5, y5, z5
      write(3,'(3ES20.10)') x6, y6, z6
      write(3,'(3ES20.10)') x7, y7, z7
      write(3,'(3ES20.10)') x8, y8, z8
      write(3,'(4i10)'    )  1, 2, 3, 4, 5, 6, 7, 8
      close(3)

     endif

     n_negative_vol = n_negative_vol + 1

     ! Display a WARNING MESSAGE
     if (n_negative_vol == 1) then
      write(*,*)
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) " Negative Hex volume deteceted. See fort.8933..."
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*)
     endif

    endif

   end do
  !-------------------------------------------------------------
  ! End of loop over hex
  !-------------------------------------------------------------

   vol_ave = vol_ave / real(nhex,dp)

   if (n_negative_vol == 0) then
     write(*,*)
     write(*,*) " Negative hex volume not detected. "
     write(*,*)
   endif

   write(*,*)
   write(*,'(a,es26.16)') " Minimum hexahedral volume = ", vol_min
   write(*,'(a,es26.16)') " Maximum hexahedral volume = ", vol_max
   write(*,'(a,es26.16)') " Average hexahedral volume = ", vol_ave
   write(*,*)

   if (n_negative_vol > 0) then
     write(*,*) " Hexahedra with negative volume created....."
     write(*,*) "  # of hex having negative volume = ", n_negative_vol
     negative_volume(ihex) = .true.
   endif

 endif hex_vol_check

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! End fo (2)Check Hexahedral elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!--------------------------------------------------------------------
! (3)Check Prismatic elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 prs_vol_check : if (nprs > 0) then

   n_negative_vol = 0
   vol_max = -1.0_dp
   vol_min =  1.0e+15_dp
   vol_ave =  0.0_dp

  !-------------------------------------------------------------
  ! Loop over prism and check the volumne one by one.
  !-------------------------------------------------------------
   do i = 1, nprs

    v(1:6) =  prs(i,1:6)

    x1 = node( v(1) )%x
    y1 = node( v(1) )%y
    z1 = node( v(1) )%z

    x2 = node( v(2) )%x
    y2 = node( v(2) )%y
    z2 = node( v(2) )%z

    x3 = node( v(3) )%x
    y3 = node( v(3) )%y
    z3 = node( v(3) )%z

    x4 = node( v(4) )%x
    y4 = node( v(4) )%y
    z4 = node( v(4) )%z

    x5 = node( v(5) )%x
    y5 = node( v(5) )%y
    z5 = node( v(5) )%z

    x6 = node( v(6) )%x
    y6 = node( v(6) )%y
    z6 = node( v(6) )%z

  !Cut prism into 3 tetrahedra and compute the volume of prism:

     volume = 0.0_dp
     volume = volume + tet_volume( x1, x2, x3, x4, &
                                   y1, y2, y3, y4, &
                                   z1, z2, z3, z4  )

     volume = volume + tet_volume( x5, x6, x3, x4, &
                                   y5, y6, y3, y4, &
                                   z5, z6, z3, z4  )

     volume = volume + tet_volume( x5, x3, x2, x4, &
                                   y5, y3, y2, y4, &
                                   z5, z3, z2, z4  )

    vol_max = max(vol_max,volume)
    vol_min = min(vol_min,volume)
    vol_ave = vol_ave + volume

    if (volume < 0.0_dp) then

     write(8933,*) " Negative volume prism at element = ", i

     !Write a Tecplot file for the negative-volume tetrahedra,
     !only for the first 20 tetra (to avoid too many files).
     if (n_negative_vol < 21) then

      write( i_char  , '(i0)' ) i
      filename = "negative_vol_prism_" // trim(i_char) // "_tec.dat"
      open(unit=3, file=filename, status="unknown", iostat=os)
      write(3,*) 'TITLE = ', "negative_vol_tet_" // trim(i_char)
      write(3,*) 'VARIABLES = "x","y","z"'
      write(3,*) 'ZONE  N=', 6,',E=', 1,' , et=tetrahedron, F=FEPOINT'
      write(3,'(3ES20.10)') x1, y1, z1
      write(3,'(3ES20.10)') x2, y2, z2
      write(3,'(3ES20.10)') x3, y3, z3
      write(3,'(3ES20.10)') x4, y4, z4
      write(3,'(3ES20.10)') x5, y5, z5
      write(3,'(3ES20.10)') x6, y6, z6
      write(3,'(4i10)'    )  1, 2, 3, 4, 5, 6
      close(3)

     endif

     n_negative_vol = n_negative_vol + 1

     ! Display a WARNING MESSAGE
     if (n_negative_vol == 1) then
      write(*,*)
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) " Negative Prism volume deteceted. See fort.8933..."
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(*,*)
     endif

    endif

   end do
  !-------------------------------------------------------------
  ! End of loop over prism
  !-------------------------------------------------------------

   vol_ave = vol_ave / real(nprs,dp)

   if (n_negative_vol == 0) then
     write(*,*)
     write(*,*) " Negative prism volume not detected."
     write(*,*)
   endif

   write(*,*)
   write(*,'(a,es26.16)') " Minimum prismatic volume = ", vol_min
   write(*,'(a,es26.16)') " Maximum prismatic volume = ", vol_max
   write(*,'(a,es26.16)') " Average prismatic volume = ", vol_ave
   write(*,*)

   if (n_negative_vol > 0) then
     write(*,*) " prisms with negative volume created....."
     write(*,*) "  # of prisms having negative volume = ", n_negative_vol
     negative_volume(iprs) = .true.
   endif

 endif prs_vol_check

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! End fo (3)Check Prismatic elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------


  end subroutine volume_check

!*******************************************************************************
! Compute the volume of a tetrahedron defined by 4 vertices:
!
!       (x1,y1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4),
!
! which are ordered as follows:
!
!            1
!            o
!           /| .
!          / |   .
!         /  |     .
!        /   |       .
!     2 o----|-------o 3
!        \   |     .
!         \  |    .
!          \ |  .
!           \|.
!            o
!            4
!
! Note: Volume = volume integral of 1 = 1/3 * volume integral of div(x,y,z) dV
!              = surface integral of (x,y,z)*dS
!              = sum of [ (xc,yc,zc)*area_vector ] over triangular faces.
!
! where the last step is exact because (x,y,z) vary linearly over the triangle.
! There are other ways to compute the volume, of course.
!
!*******************************************************************************
 function tet_volume(x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4)

 implicit none

 integer, parameter :: dp = selected_real_kind(15) !Double precision

 !Input
 real(dp), intent(in)   :: x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4
 !Output
 real(dp)               :: tet_volume

 real(dp)               :: xc, yc, zc
 real(dp), dimension(3) :: area
 integer(kd)                :: ix=1, iy=2, iz=3


 tet_volume = 0.0_dp

! Triangle 1-3-2

   !Centroid of the triangular face
      xc = (x1+x3+x2)/3.0_dp
      yc = (y1+y3+y2)/3.0_dp
      zc = (z1+z3+z2)/3.0_dp
   !Outward normal surface vector
   area = triangle_area_vector(x1,x3,x2, y1,y3,y2, z1,z3,z2)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

! Triangle 1-4-3

   !Centroid of the triangular face
      xc = (x1+x4+x3)/3.0_dp
      yc = (y1+y4+y3)/3.0_dp
      zc = (z1+z4+z3)/3.0_dp
   !Outward normal surface vector
   area = triangle_area_vector(x1,x4,x3, y1,y4,y3, z1,z4,z3)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

! Triangle 1-2-4

   !Centroid of the triangular face
      xc = (x1+x2+x4)/3.0_dp
      yc = (y1+y2+y4)/3.0_dp
      zc = (z1+z2+z4)/3.0_dp
   !Outward normal surface vector
   area = triangle_area_vector(x1,x2,x4, y1,y2,y4, z1,z2,z4)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

! Triangle 2-3-4

   !Centroid of the triangular face
      xc = (x2+x3+x4)/3.0_dp
      yc = (y2+y3+y4)/3.0_dp
      zc = (z2+z3+z4)/3.0_dp
   !Outward normal surface vector
   area = triangle_area_vector(x2,x3,x4, y2,y3,y4, z2,z3,z4)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

   tet_volume = tet_volume / 3.0_dp

 end function tet_volume

!*******************************************************************************
! Compute the area of a triangle in 3D defined by 3 vertices:
!
!       (x1,y1,z1), (x2,y2,z2), (x3,y3,z3),
!
! which is assumed to be ordered clockwise.
!
!     1             2
!      o------------o
!       \         .
!        \       . --------->
!         \    .
!          \ .
!           o
!           3
!
! Note: Area is a vector based on the right-hand rule: 
!       when wrapping the right hand around the triangle with the fingers in the
!       direction of the vertices [1,2,3], the thumb points in the positive
!       direction of the area.
!
! Note: Area vector is computed as the cross product of edge vectors [31] and [32].
!
!*******************************************************************************
 function triangle_area_vector(x1,x2,x3, y1,y2,y3, z1,z2,z3) result(area_vector)
 
 implicit none
 integer , parameter :: dp = selected_real_kind(15) !Double precision

 !Input
  real(dp), intent(in)   :: x1,x2,x3, y1,y2,y3, z1,z2,z3
 !Output
  real(dp), dimension(3) :: area_vector

  integer(kd) :: ix=1, iy=2, iz=3

   area_vector(ix) = 0.5_dp*( (y1-y3)*(z2-z3)-(z1-z3)*(y2-y3) )
   area_vector(iy) = 0.5_dp*( (z1-z3)*(x2-x3)-(x1-x3)*(z2-z3) )
   area_vector(iz) = 0.5_dp*( (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3) )

 end function triangle_area_vector


!*******************************************************************************
! Compute the unit area vector of a triangle in 3D.
!
!*******************************************************************************
 function unit_triangle_area_vector(x1,x2,x3, y1,y2,y3, z1,z2,z3) result(area_vector)
 
 implicit none
 integer , parameter :: dp = selected_real_kind(15) !Double precision

 !Input
  real(dp), intent(in)   :: x1,x2,x3, y1,y2,y3, z1,z2,z3
 !Output
  real(dp), dimension(3) :: area_vector

  real(dp) :: mag
  integer(kd)  :: ix=1, iy=2, iz=3

   area_vector(ix) = 0.5_dp*( (y1-y3)*(z2-z3)-(z1-z3)*(y2-y3) )
   area_vector(iy) = 0.5_dp*( (z1-z3)*(x2-x3)-(x1-x3)*(z2-z3) )
   area_vector(iz) = 0.5_dp*( (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3) )

   mag = sqrt(   area_vector(ix)**2 &
               + area_vector(iy)**2 &
               + area_vector(iz)**2 )

   if (mag < 1.0e-18) then
    write(*,*) "!!!!!!!!!!!!!!!!!!! Very small area vector...: Area = ", mag
    stop
   endif

   area_vector = area_vector / mag

 end function unit_triangle_area_vector

!********************************************************************************
!* Check non-planarity of quad faces.
!*
!* Note: All faces are oriented clockwise to give an outward area vector.
!*
!********************************************************************************
  subroutine non_planarity_check
  implicit none

  real(dp) :: x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4
  real(dp) :: cos_max, cos_min, cos_ave

  integer(kd)               :: i, k, nqface
  integer(kd), dimension(8) :: v

  !3 quad-faces for prism, 6 quad-faces for hex.
  integer(kd), dimension(6,4) :: face

  real(dp), dimension(3) :: unit_norm01, unit_norm02
  real(dp)               :: cos_p1, cos_p2, cos_p

   cos_max = -1.0_dp
   cos_min =  1.0e+15_dp
   cos_ave =  0.0_dp
    nqface = 0

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! (2)Loop over Hexahedral elements
!
!      8 ______________ 7 
!       /             /|
!      /             / |
!     /____________ /  |
!    5|  |         |6  |
!     |  |         |   |
!     |  |         |   |
!     |  |4________|___|3
!     |  /         |  /
!     | /          | /
!     |/___________|/
!     1            2
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 hex_check : if (nhex > 0) then

  !-------------------------------------------------------------
  ! Loop over hex and check quad faces  one by one.
  !-------------------------------------------------------------
   do i = 1, nhex

    v(1:8) =  hex(i,1:8)

    face(1, 1) = v(1)
    face(1, 2) = v(4)
    face(1, 3) = v(3)
    face(1, 4) = v(2)

    face(2, 1) = v(5)
    face(2, 2) = v(6)
    face(2, 3) = v(7)
    face(2, 4) = v(8)

    face(3, 1) = v(1)
    face(3, 2) = v(2)
    face(3, 3) = v(6)
    face(3, 4) = v(5)

    face(4, 1) = v(3)
    face(4, 2) = v(4)
    face(4, 3) = v(8)
    face(4, 4) = v(7)

    face(5, 1) = v(1)
    face(5, 2) = v(5)
    face(5, 3) = v(8)
    face(5, 4) = v(4)

    face(6, 1) = v(2)
    face(6, 2) = v(3)
    face(6, 3) = v(7)
    face(6, 4) = v(6)

    do k = 1, 6

     nqface = nqface + 1

     x1 = node( face(k,1) )%x
     y1 = node( face(k,1) )%y
     z1 = node( face(k,1) )%z

     x2 = node( face(k,2) )%x
     y2 = node( face(k,2) )%y
     z2 = node( face(k,2) )%z

     x3 = node( face(k,3) )%x
     y3 = node( face(k,3) )%y
     z3 = node( face(k,3) )%z

     x4 = node( face(k,4) )%x
     y4 = node( face(k,4) )%y
     z4 = node( face(k,4) )%z

   !Compute the dot product of unit normals of triangular faces within the quad:
   !There are two possible pairs of triangular faces.
   !The dot product is equal to cos(angle) because the normals are unit vectors.

     unit_norm01 = unit_triangle_area_vector(x1,x2,x3, y1,y2,y3, z1,z2,z3)
     unit_norm02 = unit_triangle_area_vector(x1,x3,x4, y1,y3,y4, z1,z3,z4)

     cos_p1 = unit_norm01(1)*unit_norm02(1) &
            + unit_norm01(2)*unit_norm02(2) &
            + unit_norm01(3)*unit_norm02(3)

     unit_norm01 = unit_triangle_area_vector(x1,x2,x4, y1,y2,y4, z1,z2,z4)
     unit_norm02 = unit_triangle_area_vector(x2,x3,x4, y2,y3,y4, z2,z3,z4)

     cos_p2 = unit_norm01(1)*unit_norm02(1) &
            + unit_norm01(2)*unit_norm02(2) &
            + unit_norm01(3)*unit_norm02(3)

     cos_p = max(cos_p1, cos_p2)

    cos_max = max(cos_max,cos_p)
    cos_min = min(cos_min,cos_p)
    cos_ave = cos_ave + cos_p

    end do

   end do
  !-------------------------------------------------------------
  ! End of loop over hex
  !-------------------------------------------------------------

 endif hex_check

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! End fo (2)Loop over Hexahedral elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! (2)Loop over Prismatic elements
!
!        6
!       / .
!      /      .
!     /________ 
!    4|  |      |5
!     |  |      |
!     |  |      |
!     |  |3     |
!     |  /  .   |
!     | /     . |
!     |/________|
!     1            2
!
!--------------------------------------------------------------------
!--------------------------------------------------------------------
 prs_check : if (nprs > 0) then

  !-------------------------------------------------------------
  ! Loop over prisms and check quad faces one by one.
  !-------------------------------------------------------------
   do i = 1, nprs

    v(1:6) =  prs(i,1:6)

    face(1, 1) = v(1)
    face(1, 2) = v(2)
    face(1, 3) = v(5)
    face(1, 4) = v(4)

    face(2, 1) = v(3)
    face(2, 2) = v(1)
    face(2, 3) = v(4)
    face(2, 4) = v(6)

    face(3, 1) = v(2)
    face(3, 2) = v(3)
    face(3, 3) = v(6)
    face(3, 4) = v(5)

    do k = 1, 3

     nqface = nqface + 1

     x1 = node( face(k,1) )%x
     y1 = node( face(k,1) )%y
     z1 = node( face(k,1) )%z

     x2 = node( face(k,2) )%x
     y2 = node( face(k,2) )%y
     z2 = node( face(k,2) )%z

     x3 = node( face(k,3) )%x
     y3 = node( face(k,3) )%y
     z3 = node( face(k,3) )%z

     x4 = node( face(k,4) )%x
     y4 = node( face(k,4) )%y
     z4 = node( face(k,4) )%z

   !Compute the dot product of unit normals of triangular faces within the quad:
   !There are two possible pairs of triangular faces.
   !The dot product is equal to cos(angle) because the normals are unit vectors.

     unit_norm01 = unit_triangle_area_vector(x1,x2,x3, y1,y2,y3, z1,z2,z3)
     unit_norm02 = unit_triangle_area_vector(x1,x3,x4, y1,y3,y4, z1,z3,z4)

     cos_p1 = unit_norm01(1)*unit_norm02(1) &
            + unit_norm01(2)*unit_norm02(2) &
            + unit_norm01(3)*unit_norm02(3)

     unit_norm01 = unit_triangle_area_vector(x1,x2,x4, y1,y2,y4, z1,z2,z4)
     unit_norm02 = unit_triangle_area_vector(x2,x3,x4, y2,y3,y4, z2,z3,z4)

     cos_p2 = unit_norm01(1)*unit_norm02(1) &
            + unit_norm01(2)*unit_norm02(2) &
            + unit_norm01(3)*unit_norm02(3)

     cos_p = max(cos_p1, cos_p2)

    cos_max = max(cos_max,cos_p)
    cos_min = min(cos_min,cos_p)
    cos_ave = cos_ave + cos_p

    end do

   end do
  !-------------------------------------------------------------
  ! End of loop over hex
  !-------------------------------------------------------------

 endif prs_check

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! End of (2)Loop over Prismatic elements
!--------------------------------------------------------------------
!--------------------------------------------------------------------


   cos_ave = cos_ave / real(nqface,dp)

   write(*,*)
   write(*,*) " Non-planarity:"
   write(*,*) "   Each quad face is split into two triangles: T1 and T2."
   write(*,*) "   Dot product of unit tria-face normals computed: n(T1)*n(T2)."
   write(*,*) "   Two possible splits considered, and the larger is taken."
   write(*,*)
   write(*,'(a,es26.16)') " Minimum  = ", cos_min
   write(*,'(a,es26.16)') " Maximum  = ", cos_max
   write(*,'(a,es26.16)') " Average  = ", cos_ave
   write(*,*)

 end subroutine non_planarity_check

!********************************************************************************
!* Compute mesh spacing... Many ways to define the spacing!
!********************************************************************************
! subroutine effective_mesh_spacing

! implicit none

! real(dp), dimension(3) :: h_vol
! real(dp), dimension(3) :: h_cell, h_node

! end subroutine effective_mesh_spacing

!********************************************************************************
!* Statistics for tetra grid.
!*
!********************************************************************************

!*******************************************************************************
 subroutine tet_statistics(filename)

 character(80), intent(in) :: filename

 integer(kd), allocatable, dimension(:)   :: n_n2tet
 integer(kd), allocatable, dimension(:,:) :: n2tet
 integer(kd), allocatable, dimension(:,:) :: n2tet_k
 
 integer(kd)                              :: nedges
 integer(kd), allocatable, dimension(:,:) :: edge

 integer(kd)                              :: nfaces_i, nfaces_b, nfaces
 integer(kd), allocatable, dimension(:,:) :: face

 integer(kd), allocatable, dimension(:)   :: nc_nnghbrs
 integer(kd), allocatable, dimension(:,:) :: nc_nghbr
 integer(kd), allocatable, dimension(:)   :: nc_nnghbrs2
 integer(kd), allocatable, dimension(:,:) :: nc_nghbr2

 integer(kd), allocatable, dimension(:)   :: tet_nnghbrs
 integer(kd), allocatable, dimension(:,:) :: tet_nghbr
 integer(kd), allocatable, dimension(:)   :: tet_nnghbrs2
 integer(kd), allocatable, dimension(:,:) :: tet_nghbr2

 integer(kd) :: i, j, k, m, kv
 integer(kd) :: n1, n2
 integer(kd), dimension(4,3) :: nghbr

 integer(kd), allocatable, dimension(:)   :: bmark
 integer(kd), allocatable, dimension(:)   :: bmark_tet

 logical :: found
 integer(kd) :: itet
 integer(kd) :: v1, v2, v3, v4
 integer(kd) :: euler_residual

 character(80) :: temp_char

 open(unit=77, file=filename, status="unknown", iostat=os)

  write(77,*)
  write(77,*)
  write(77,*) "*********************************************************"
  write(77,*) "*********************************************************"
  write(77,*) "*********************************************************"
  write(77,*) " Some statistics for a tetra grid"
  write(77,*) "*********************************************************"
  write(77,*) "*********************************************************"
  write(77,*) "*********************************************************"
  write(77,*)

 !Neighbor nodes within a tetrahedron.

  nghbr(1,1) = 2
  nghbr(1,2) = 3
  nghbr(1,3) = 4

  nghbr(2,1) = 1
  nghbr(2,2) = 4
  nghbr(2,3) = 3

  nghbr(3,1) = 1
  nghbr(3,2) = 2
  nghbr(3,3) = 4

  nghbr(4,1) = 1
  nghbr(4,2) = 3
  nghbr(4,3) = 2

!--------------------------------------------------------------
!--------------------------------------------------------------
! Boundary mark: 0 for interior, 1 for boundary
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing bmark..."
  write(77,*)

  allocate(bmark(nnodes))
  bmark = 0

  do i = 1, ntrias_b
   do k = 1, 3
    bmark( tri(i,k) )  = -1
   end do
  end do

  do i = 1, nquads_b
   do k = 1, 4
    bmark( quad(i,k) ) = -1
   end do
  end do

!--------------------------------------------------------------
!--------------------------------------------------------------
! List of tets around each node.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of tets around each node..."
  write(77,*) " - Tetrahedral elements sharing the same node."
  write(77,*)

  allocate( n_n2tet(nnodes)  )

! Count the number of tetrahedra sharing the same node:
    n_n2tet = 0
  do i = 1, ntet
   do k = 1, 4
    n_n2tet( tet(i,k) ) = n_n2tet( tet(i,k) ) + 1 !# of tets around the node tet(i,k).
   end do
  end do

! (TBI) Add the number of prismatic elements.

! (TBI) Add the number of hexahedral  elements.


! Allocate the array to store the neighboring elments at each node.
! Use the maximum n_node2tet.

  allocate( n2tet  (nnodes, maxval(n_n2tet) ) )
  allocate( n2tet_k(nnodes, maxval(n_n2tet) ) )

! Construct the list of elements around each node.

  n_n2tet = 0

  do i = 1, ntet
   do k = 1, 4
    n_n2tet( tet(i,k) ) = n_n2tet( tet(i,k) ) + 1   !# of tets around the node tet(i,k).
      n2tet( tet(i,k)   , n_n2tet( tet(i,k) ) ) = i !Add tet i to the node tet(i,k).
    n2tet_k( tet(i,k)   , n_n2tet( tet(i,k) ) ) = k !Which vertex.
   end do
  end do

 !-------------------------------------------------
 ! Compute and display the data:
 !-------------------------------------------------

   temp_char= "Tets around node"
   call histogram(temp_char,nnodes,n_n2tet,bmark,-100000)

   temp_char= "Tets around node (interior)"
   call histogram(temp_char,nnodes,n_n2tet,bmark,-1)

   temp_char= "Tets around node (boundary)"
   call histogram(temp_char,nnodes,n_n2tet,bmark,0)

!--------------------------------------------------------------
!--------------------------------------------------------------
! List of nodes around each node.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of nodes around each node..."
  write(77,*) " - Edge-connected neighbor nodes for each node."
  write(77,*)

! real(dp) :: ave_n2tet, min_nc_nghbr, max_nc_nghbr
! allocate(hist_nc_nghbr)

  allocate( nc_nnghbrs(nnodes) )
  nc_nnghbrs = 0
  allocate( nc_nghbr(nnodes, maxval(n_n2tet) ) )

 !-------------------------------------------------
 ! Find edge-connected neghbr nodes from tets.
  do i = 1, ntet

  !-------------------------------------------------
  ! Loop over the vertices.
   do k = 1, 4

     n1 = tet(i,k)
     
   !-------------------------------------------------
   ! Loop over the rest of vertices.
    do j = 1, 3
     n2 = tet(i,nghbr(k,j))

     if ( nc_nnghbrs(n1) == 0 ) then

      nc_nnghbrs(n1) = 1
      nc_nghbr(n1,1) = n2

     else

     !Check if the neighbor is already added.
      found = .false.
      do m = 1, nc_nnghbrs(n1)
       if (n2 == nc_nghbr(n1,m)) then
        found = .true.
        exit
       endif
      end do

     !If not added yet, add it to the list.
      if (.not.found) then
       nc_nnghbrs(n1) = nc_nnghbrs(n1) + 1
       nc_nghbr(n1,nc_nnghbrs(n1)) = n2
      endif

     endif

    end do
   !-------------------------------------------------
   ! End of Loop over the rest of vertices.

   end do
  !-------------------------------------------------
  ! End of Loop over the vertices.

  end do
 ! End of Find edge-connected neghbr nodes from tets.
 !-------------------------------------------------
 
 !-------------------------------------------------
 ! (TBI) Find edge-connected neghbr nodes from prism.

 !-------------------------------------------------
 ! (TBI) Find edge-connected neghbr nodes from hexa.

 !-------------------------------------------------
 ! Compute and display the data:
 !-------------------------------------------------

   temp_char= "Node neighbors"
   call histogram(temp_char,nnodes,nc_nnghbrs,bmark,-100000)

   temp_char= "Node neighbors (interior)"
   call histogram(temp_char,nnodes,nc_nnghbrs,bmark,-1)

   temp_char= "Node neighbors (boundary)"
   call histogram(temp_char,nnodes,nc_nnghbrs,bmark,0)

!--------------------------------------------------------------
!--------------------------------------------------------------
! List of level2-nodes around each node.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of lvl2-nodes around each node..."
  write(77,*) " - Edge-connected neighbors and their neighbors."
  write(77,*)

  allocate( nc_nnghbrs2(nnodes) )
  nc_nnghbrs2 = 0
  allocate( nc_nghbr2(nnodes, 10*maxval(n_n2tet) ) )

 !-------------------------------------------------
  do i = 1, nnodes

  !Begin with the edge-connected neighbors.
    nc_nnghbrs2(i) = nc_nnghbrs(i)
   do k = 1, nc_nnghbrs(i)
    nc_nghbr2(i,k) = nc_nghbr(i,k)
   end do

  !-------------------------------------------------
  ! Loop over the neighbors
   do k = 1, nc_nnghbrs(i)

     n1 = nc_nghbr(i,k)
     
   !-------------------------------------------------
   ! Loop over the neighbors of the neighbor.
    do j = 1, nc_nnghbrs(n1)
     n2 = nc_nghbr(n1,k)

     !Check if the neighbor is already added.
      found = .false.
      do m = 1, nc_nnghbrs2(i)
       if (n2 == nc_nghbr2(i,m) .or. n2 == i) then
        found = .true.
        exit
       endif
      end do

     !If not added yet, add it to the list.
      if (.not.found) then
       nc_nnghbrs2(i) = nc_nnghbrs2(i) + 1
       nc_nghbr2(i,nc_nnghbrs2(i)) = n2
      endif

    end do
   !-------------------------------------------------

   end do
  !-------------------------------------------------

  end do
 !-------------------------------------------------

 !-------------------------------------------------
 ! Compute and display the data:
 !-------------------------------------------------

   temp_char= "Level2 node-nghbrs"
   call histogram(temp_char,nnodes,nc_nnghbrs2,bmark,-100000)

   temp_char= "Level2 node-nghbrs (interior)"
   call histogram(temp_char,nnodes,nc_nnghbrs2,bmark,-1)

   temp_char= "Level2 node-nghbrs (boundary)"
   call histogram(temp_char,nnodes,nc_nnghbrs2,bmark,0)

!--------------------------------------------------------------
!--------------------------------------------------------------
! Edge list
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of edges..."
  write(77,*)

 !---------------------------------------
 ! First, count the number of edges.

  nedges = 0
  do i = 1, nnodes
   do k = 1, nc_nnghbrs(i)
    if ( i < nc_nghbr(i,k) ) then
      nedges = nedges + 1
    endif
   end do
  end do

 !---------------------------------------
 ! Allocate and fill edge array.

  allocate(edge(nedges,2))

  nedges = 0
  do i = 1, nnodes
   do k = 1, nc_nnghbrs(i)
    if ( i < nc_nghbr(i,k) ) then
      nedges = nedges + 1
      edge(nedges,1) = i
      edge(nedges,2) = nc_nghbr(i,k)
    endif
   end do
  end do

 !-------------------------------------------------
 ! Compute and display the data:
 !-------------------------------------------------

  write(77,*) " Number of edges = ", nedges
  write(77,*)
  write(77,*) " P2 element nodes = ", nnodes, "+", nedges, "=", nnodes + nedges
  write(77,*)

!--------------------------------------------------------------
!--------------------------------------------------------------
! Construct tet neighbor list.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of tet neighbors..."
  write(77,*) " - Face-neighbors of each tetrahedral element."
  write(77,*)

  allocate(tet_nnghbrs(ntet))
  tet_nnghbrs = 0
  allocate(tet_nghbr(ntet,4))

  do i = 1, ntet

   !--------------------------------------------------------
   ! Face 123
   !--------------------------------------------------------
    v1 = tet(i,1)
    v2 = tet(i,2)
    v3 = tet(i,3)
    v4 = 4

    found = .false.

    do k = 1, n_n2tet( v1    )

     itet = n2tet(     v1, k )
       kv = n2tet_k(   v1, k )

     if ( tet(itet,kv) /= v1 ) then 
      write(77,*) " Error in the surrounding tets..."
      stop
     endif

     if     ( tet(itet, nghbr(kv,2)) == v2 .and. tet(itet, nghbr(kv,1) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,3)) == v2 .and. tet(itet, nghbr(kv,2) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,1)) == v2 .and. tet(itet, nghbr(kv,3) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     endif

     if (found) exit

    end do

    !--------------------------------------------------------
    !If not found, this must be a boundary face. Assign -1.
     if (.not.found) then
      tet_nghbr(i,v4) = -1
     endif

   !--------------------------------------------------------
   ! Face 141
   !--------------------------------------------------------
    v1 = tet(i,1)
    v2 = tet(i,4)
    v3 = tet(i,2)
    v4 = 3

    found = .false.

    do k = 1, n_n2tet( v1    )

     itet = n2tet(     v1, k )
       kv = n2tet_k(   v1, k )

     if ( tet(itet,kv) /= v1 ) then 
      write(77,*) " Error in the surrounding tets..."
      stop
     endif

     if     ( tet(itet, nghbr(kv,2)) == v2 .and. tet(itet, nghbr(kv,1) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,3)) == v2 .and. tet(itet, nghbr(kv,2) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,1)) == v2 .and. tet(itet, nghbr(kv,3) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     endif

     if (found) exit

    end do

    !--------------------------------------------------------
    !If not found, this must be a boudnary face. Assign -1.
     if (.not.found) then
      tet_nghbr(i,v4) = -1
     endif

   !--------------------------------------------------------
   ! Face 134
   !--------------------------------------------------------
    v1 = tet(i,1)
    v2 = tet(i,3)
    v3 = tet(i,4)
    v4 = 2

    found = .false.

    do k = 1, n_n2tet( v1    )

     itet = n2tet(     v1, k )
       kv = n2tet_k(   v1, k )

     if ( tet(itet,kv) /= v1 ) then 
      write(77,*) " Error in the surrounding tets..."
      stop
     endif

     if     ( tet(itet, nghbr(kv,2)) == v2 .and. tet(itet, nghbr(kv,1) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,3)) == v2 .and. tet(itet, nghbr(kv,2) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,1)) == v2 .and. tet(itet, nghbr(kv,3) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     endif

     if (found) exit

    end do

    !--------------------------------------------------------
    !If not found, this must be a boudnary face. Assign -1.
     if (.not.found) then
      tet_nghbr(i,v4) = -1
     endif

   !--------------------------------------------------------
   ! Face 243
   !--------------------------------------------------------
    v1 = tet(i,2)
    v2 = tet(i,4)
    v3 = tet(i,3)
    v4 = 1

    found = .false.

    do k = 1, n_n2tet( v1    )

     itet = n2tet(     v1, k )
       kv = n2tet_k(   v1, k )

     if ( tet(itet,kv) /= v1 ) then 
      write(77,*) " Error in the surrounding tets..."
      stop
     endif

     if     ( tet(itet, nghbr(kv,2)) == v2 .and. tet(itet, nghbr(kv,1) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,3)) == v2 .and. tet(itet, nghbr(kv,2) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     elseif ( tet(itet, nghbr(kv,1)) == v2 .and. tet(itet, nghbr(kv,3) ) == v3 ) then

                found = .true.
      tet_nnghbrs(i)  = tet_nnghbrs(i) + 1
      tet_nghbr(i,v4) = itet

     endif

     if (found) exit

    end do

    !--------------------------------------------------------
    !If not found, this must be a boudnary face. Assign -1.
     if (.not.found) then
      tet_nghbr(i,v4) = -1
     endif

  end do


  allocate(bmark_tet(ntet))
  bmark_tet = 0

  do i = 1, ntet
   if ( minval(tet_nghbr(i,:)) == -1 ) then
    bmark_tet( i ) = -1
   endif
  end do

 !-------------------------------------------------
 ! Compute and display the data:
 !-------------------------------------------------

   temp_char= "cc-Tet-nghbrs"
   call histogram(temp_char,ntet,tet_nnghbrs,bmark_tet,-100000)

   temp_char= "cc-Tet-nghbrs (interior)"
   call histogram(temp_char,ntet,tet_nnghbrs,bmark_tet,-1)

   temp_char= "cc-Tet-nghbrs (boundary)"
   call histogram(temp_char,ntet,tet_nnghbrs,bmark_tet,0)


!--------------------------------------------------------------
!--------------------------------------------------------------
! Construct lvl2-tet neighbor list.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of lvl2-tet neighbors..."
  write(77,*) " - Face-neighbors and their neighbors."
  write(77,*)

  allocate(tet_nnghbrs2(ntet))
  tet_nnghbrs2 = 0
  allocate(tet_nghbr2(ntet,50))

 !Loop over tets
  do i = 1, ntet

    tet_nnghbrs2(i) = tet_nnghbrs(i)
   do k = 1, 4
    tet_nghbr2(i,k) = tet_nghbr(i,k)
   end do

  !Loop over tet neighbors
   do k = 1, tet_nnghbrs(i)
    n1 = tet_nghbr(i,k)
    if (n1 < 0) cycle

   !Loop over tet neighbors of the tet neighbors
    do j = 1, tet_nnghbrs(n1)
    n2 = tet_nghbr(n1,j)
    if (n2 < 0) cycle

     !Check if the neighbor is already added.
      found = .false.
      do m = 1, tet_nnghbrs2(i)
       if (n2 == tet_nghbr2(i,m) .or. n2 == i) then
        found = .true.
        exit
       endif
      end do

     !If not added yet, add it to the list.
      if (.not.found) then
       tet_nnghbrs2(i) = tet_nnghbrs2(i) + 1
       tet_nghbr2(i,tet_nnghbrs2(i)) = n2
      endif

    end do

   end do

  end do


 !-------------------------------------------------
 ! Compute and display the data:
 !-------------------------------------------------

   temp_char= "Level2 cc-tet-nghbrs"
   call histogram(temp_char,ntet,tet_nnghbrs2,bmark_tet,-100000)

   temp_char= "Level2 cc-tet-nghbrs (interior)"
   call histogram(temp_char,ntet,tet_nnghbrs2,bmark_tet,-1)

   temp_char= "Level2 cc-tet-nghbrs (boundary)"
   call histogram(temp_char,ntet,tet_nnghbrs2,bmark_tet,0)


!--------------------------------------------------------------
!--------------------------------------------------------------
! Construct face list.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Constructing a list of faces..."
  write(77,*)

  nfaces_i = 0
  nfaces_b = 0

  do i = 1, ntet
   do k = 1, 4
    if ( tet_nghbr(i,k) > 0 ) then
     if ( i < tet_nghbr(i,k) ) nfaces_i = nfaces_i + 1
    else
     nfaces_b = nfaces_b + 1
    endif
   end do
  end do

   nfaces = nfaces_i + nfaces_b

  allocate( face(nfaces,2) )

  nfaces = 0

 !Loop over tets
  do i = 1, ntet

  !Loop over tet neighbors
   do k = 1, 4
    n1 = tet_nghbr(i,k)

    if ( i < n1 ) then
     nfaces = nfaces + 1
     face(nfaces,1) = i
     face(nfaces,2) = n1
    endif

    if (n1 < 0) then
     nfaces = nfaces + 1
     face(nfaces,1) = i
     face(nfaces,2) = n1
    endif

   end do

  end do

  write(77,*) " Number of interior faces = ", nfaces_i
  write(77,*) " Number of boundary faces = ", nfaces_b
  write(77,*) " ---------------------------------------"
  write(77,*) " Number of    total faces = ", nfaces, nfaces_i+nfaces_b
  write(77,*)


!--------------------------------------------------------------
!--------------------------------------------------------------
! Check with Euler's formula: nnodes-nedges+nfaces+ncells=0.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Factors:"
  write(77,*)
  write(77,*) "   nnodes = ", nnodes
  write(77,*) "     ntet = ", ntet  , ": /nnodes = ", real  (ntet,dp)/real(nnodes,dp)
  write(77,*) "   nedges = ", nedges, ": /nnodes = ", real(nedges,dp)/real(nnodes,dp)
  write(77,*) "   nfaces = ", nfaces, ": /nnodes = ", real(nfaces,dp)/real(nnodes,dp)
  write(77,*) " P2 nodes = ", nnodes+nedges, ": /nnodes = ", &
              real(nnodes+nedges,dp)/real(nnodes,dp)
  write(77,*)
  write(77,*)

!--------------------------------------------------------------
!--------------------------------------------------------------
! Check with Euler's formula: nnodes-nedges+nfaces+ncells=0.
!--------------------------------------------------------------
!--------------------------------------------------------------

  write(77,*) "*********************************************************"
  write(77,*) " Check with Euler's formula: nnodes-nedges+nfaces-ncells=0,"
  write(77,*) " NOTE: ncells includes outside the domain: ncells <- ncells+1"
  write(77,*)
    euler_residual = nnodes - nedges + nfaces - (ntet + 1)
  write(77,*) "   Euler's formula (n-e+f-(c+1))  = ", euler_residual
  write(77,*)


  close(77)


 end subroutine tet_statistics
!*******************************************************************************
!*******************************************************************************
 subroutine histogram(dataname,ndata,data,flag,skip_value)
 implicit none

 character(80)                , intent(in) :: dataname
 integer(kd),                   intent(in) :: ndata
 integer      ,                 intent(in) :: skip_value
 integer(kd), dimension(ndata), intent(in) :: data
 integer(kd), dimension(ndata), intent(in) :: flag

 real(dp) :: ave_q
 integer(kd)  :: min_q, max_q, ndata2
 integer(kd), dimension(:), allocatable :: hist_q
 
 ndata2 = 0
  ave_q =  0.0_dp
  min_q =  100000000
  max_q = -100000000

  write(77,*)
  write(77,*) " --------------------------------------------------------"
  write(77,*) " --------------------------------------------------------"
  write(77,'(10x,2a)') " Histogram for ", trim(dataname)
  write(77,*) " --------------------------------------------------------"
  write(77,*) " --------------------------------------------------------"

  do i = 1, ndata
   if (flag(i)==skip_value) cycle
   ndata2 = ndata2 + 1
    ave_q =     ave_q + real(data(i),dp)
    min_q = min(min_q,  data(i) )
    max_q = max(max_q,  data(i) )
  end do
    ave_q = ave_q / real(ndata2,dp)

  write(77,*) "  average = ", ave_q
  write(77,*) "  minimum = ", min_q
  write(77,*) "  maximum = ", max_q

  allocate( hist_q(max_q) )

   hist_q = 0
  do i = 1, ndata
   if (flag(i)==skip_value) cycle
   hist_q( data(i) ) = hist_q( data(i) ) + 1
  end do

   write(77,*) " --------------------------------------------------------"
   write(77,'(3x,2a)') trim(dataname),"   :  Occurence   "
  do k = 1, max_q
   if ( hist_q(k) == 0 ) cycle
   write(77,'(i20,a7,i10)') k, " : ", hist_q(k)
  end do
   write(77,*) " --------------------------------------------------------"

   write(77,*) " --------------------------------------------------------"
   write(77,*) " --------------------------------------------------------"
   write(77,*)

  deallocate( hist_q )

 end subroutine histogram
!*******************************************************************************



end program hemisphere_cylinder_grid

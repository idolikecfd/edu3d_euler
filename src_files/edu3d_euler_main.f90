!********************************************************************************
! Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
! This program reads an unstructured grid (.ugrid) and solves the Euler equations
! in 3D for steady inviscid flow problems.
!
!
!  - This is a serial code. Think about how it can be parallelized.
!
!  - Node-centered edge-based discretization
!    -- Solutions stored at nodes, flux computed at edge-midpoints.
!    -- Numerical solutions are point values, not cell averages.
!    -- 1st-order accurate on arbitray grids.
!    -- 2nd-order accurate on arbitrary tetrahedra, and regular prism/hex/mixed.
!    -- 3rd-order accurate on arbitrary tetrahedra (not implemented).
!    -- It can be considered as finite-volume based on edge-based flux quadrature.
!    -- It can be considered as finite-element based on edge-based flux quadrature.
!    -- It can be considered as conservative finite-difference scheme.
!
!  - Steady solver (nonlinear solver) can be chosen from
!    -- Explicit forward-Euler pseudo-time iteration
!    -- Explicit 2-stage RK pseudo-time iteration
!    -- Implicit defect-correction solver with exact first-order Jacobians
!    -- Jacobian-Free Newton-Krylov solver with GCR
!
!  - Boundary conditions: all implemented in the weak form through a numerical flux.
!    -- Free stream
!    -- Subsonic outflow (back pressure)
!    -- Supersonic outflow
!    -- Slip wall
!
!  - Input parameters controlled by namelist variables set in 'input.nml'.
!
!
!  - Perform a numerical truncation error analysis to verify the residual implementation.
!
!
!  - Latest reference on the node-centered edge-based discretization:
!
!    H. Nishikawa and Y. Liu, " Accuracy-Preserving Source Term Quadrature for
!    Third-Order Edge-Based Discretization", Journal of Computational Physics,
!    Volume 344, September 2017, Pages 595-622.
!    https://doi.org/10.1016/j.jcp.2017.04.075
!    http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_liu_jcp2017_preprint.pdf
!
!  - Example cases are included in the package with grid generation codes:
!    (1)Free stream in a cube.
!    (2)Subsonic flow over a hemisphere cylinder
!    (3)Subsonic flow over an ONERA M6 wing
!
!  - Things that are missing and you might want to implement by yourself.
!
!    -- Computation of the lift/drag coefficient.
!    -- Symmetry/periodic boundary conditions.
!    -- CFL ramping (gradually increase CFL from the start of iteration)
!    -- Order-sequencing: begin with 1st-order, and then switch to 2nd-order.
!    -- Other inviscid flux functions (e.g., HLLC, AUSM, LDFSS, SLAU2, etc.)
!    -- Kappa reconstruction scheme (kappa=1/2 or kappa=1/3?)
!    -- Extend to unsteady problem (e.g., add BDF2 term).
!    -- Extend to the Navier-Stokes equations (add a viscous flux, and no-slip bc).
!       See http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_AIAA-2011-3044.pdf
!       for a simple and practical viscous numerical flux that can be added easily.
!    -- Extend to other types of elements, prism/pyramid/hexahedra.
!    -- Limiter functions for shock waves.
!    -- Mechanisms to cure low-Mach accuracy/instability problems.
!    -- HANIM for a better control of the JFNK solver.
!    -- Node reordering for efficiency
!    -- Grid coloring for multi-color Gauss-Seidel
!    -- Local-preconditioning to accelerate convergence (and low-Mach cure)
!    -- Other grid file formats (other than just .ugrid).
!    -- Other output data file formats (other than just Tecplot)
!    -- Parallelization for faster computation
!    -- Other nonlinear solvers (e.g., multigrid)
!    -- Subdivision of prism/pyramid/hexahedra into tetrahedra.
!    -- 3rd-order edge-based method (quadratic-LSQ/extrapolation/source).
!    -- Cell-centered FV/DG solver (I think this is very easy to do).
!    -- Turbulence models for RANS (or implicit LES or DNS?)
!    -- Grid adaptation capability
!
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! This is Version 4.4 (2020).
!
! 02/28/21: Implemented write/read solution files.
! 02/27/21: Modified slip BC and boundary flux quadrature. Thanks to Li Juan Chan.
! 02/16/21: Corrected the normal vectors at nodes thanks to Li Juan Chan.
! 01/23/21: Added an option to run a numerical TE analysis (use test_mms_te=T).
! 01/23/21: Added the kappa recnstruction scheme (UMUSCL). Specify kappa.
! 01/08/21: Fixed another bug (about bmark) reported by Li Juan Chan.
! 12/26/20: Removed bugs thanks to Li Juan Chan.
!
!
! This F90 program is written and made available for an educational purpose.
!
! This file may be updated in future.
!
! Katate Masatsuka, February 2021. http://www.cfdbooks.com
!
! Notes:
!
!  The purpose of this code is to give a beginner an opportunity to learn how to
!  write an unstructured CFD code. Hence, the focus here is on the simplicity.
!  The code is not optimized for efficiency.
!
!  If the code is not simple enough to understand, please send questions to Hiro
!  at sunmasen(at)hotmail.com. He'll greatly appreciate it and revise the code.
!
!  If the code helps you understand how to write your own code that is more
!  efficient and has more features, it'll have served its purpose.
!
!********************************************************************************

 program edu3d_euler

 ! This main program contains the following subroutines:
 !
 !  - set_filenames                ! define file names
 !  - allocate_nceb_arrays         ! allocate arrays for the edge-based scheme
 !  - set_initial_solution         ! set up an initial solution
 !  - read_grid                    ! read a grid file
 !  - write_tecplot_boundary_file  ! write a tecplot boundary grid file
 !  - write_tecplot_volume_file    ! write a tecplot volume grid file
 !  - big_endian_io                ! detect endian

  use input_parameter, only : read_nml_input_parameters , &
                                    generate_tec_file_b , &
                                    generate_tec_file_v , &
                                    generate_soln_file  , &
                                    read_soln_file      , &
                                    test_mms_te

  use data_module    , only : filename_soln

  use nceb_data      , only : construct_nceb_data

  use steady_solver  , only : solve_steady
  use data_module    , only : nnodes, x, y, z, w

  implicit none

  integer , parameter :: p2 = selected_real_kind(p=15)

  integer  :: i, os
  real(p2) :: xp, yp, zp 

   write(*,*)
   write(*,*) "----------------------------------------------------------------"
   write(*,*)
   write(*,*) "  This is EDU3D Euler, version 4.1........"
   write(*,*)
   write(*,*) "----------------------------------------------------------------"
   write(*,*)

!*******************************************************************************
! Read the input parameters, defined in the file named as 'input_coarsen.nml'.
!*******************************************************************************

   write(*,*) "Reading the input file: input.nml..... "
   write(*,*)
   call read_nml_input_parameters("input.nml")
   write(*,*)

!-------------------------------------------------------------------------------
! Perform a TE test by MMS and stop if you add "test_mms_te=T" in the input file.
!
!  (1) Point evaluation of MMS source: Add mms_s_quadrature = "point" - Default
!  (2) Special source discretization : Add mms_s_quadrature = "kappa_s" and kappa_s=0.3.
!      Kapa_s is a parameter. Third-order accurate if also add oukappa=0.5 for the
!      kappa reconstruction scheme and kappa_s=0.3. This is UMUSCL-SSQ, See
!
!       "Resolving Confusion Over Third-Order Accuracy of U-MUSCL", AIAA2021-0056.
!

  if (test_mms_te) then

   call perform_te_test_with_mms
   stop

  endif

!-------------------------------------------------------------------------------
! 1. Set file names

  call set_filenames

!-------------------------------------------------------------------------------
! 2. Read a grid file (.ugrid) and a bc file (.bcmap)

  call read_grid

!-------------------------------------------------------------------------------
! 3. Construct grid data required for the node-centered edge-based discretization.

  call construct_nceb_data

!-------------------------------------------------------------------------------
! 4. Allocate arrays: solution, residual, jacobian, etc.

  call allocate_nceb_arrays

!-------------------------------------------------------------------------------
! 5. Set up the initial solution with/without random perturbation or read 
!    a solution file.

  if (read_soln_file) then

   open(unit=18, file=filename_soln, status="unknown", iostat=os)
   do i = 1, nnodes
    read(8,*) xp,yp,zp, w(i,1:5)
    if ( abs(xp-x(i))+abs(yp-y(i))+abs(zp-z(i)) > 1.0e-10_p2 ) then
     write(*,*) " Reading a soln file: Error = ", abs(xp-x(i))+abs(yp-y(i))+abs(zp-z(i))
     write(*,*) " Read (x,y,z) = ", xp  ,yp  ,zp
     write(*,*) " Grid (x,y,z) = ", x(i),y(i),z(i)
     stop
    endif
   end do
   close(18)

  else

   call set_initial_solution

  endif

!-------------------------------------------------------------------------------
! 6. Compute the solution.

  !Compute the steady soution by solving Res(U)=0.

      call  solve_steady

  !Compute the unsteady soution by solving dU/dt+Res(U)=0:

  !   call  solve_unsteady (not implemented...)

!-------------------------------------------------------------------------------
! 7. Write output files for visualization

  !Tecplot boundary grid data file if requested:
   if (generate_tec_file_b) call write_tecplot_boundary_file

  !Tecplot volume grid data file if requested:
   if (generate_tec_file_v) call write_tecplot_volume_file

  !Solution data file if requested:
   if (generate_soln_file) call write_solution_file

   write(*,*) 
   write(*,*) " See fort.1000 for GS relaxation convergence history (implicit/jfnk)."
   write(*,*) " See fort.2000 for GCR convergence history (jfnk)."
   write(*,*)

   write(*,*)
   write(*,*) "----------------------------------------------------------------"
   write(*,*)
   write(*,*) "  End of EDU3D Euler. "
   write(*,*)
   write(*,*) "----------------------------------------------------------------"
   write(*,*)

 stop

 contains
 
!********************************************************************************
!
! Define file names.
!
!********************************************************************************
 subroutine set_filenames

  use input_parameter, only : project, ugrid_unformatted
  use data_module    , only : filename_ugrid    , filename_bcmap
  use data_module    , only : filename_tecplot_b, filename_tecplot_v
  use data_module    , only : filename_soln

 implicit none

 !-----------------------------------------------------------------------
 ! Input grid file (.ugrid)

   !(1) Binary file
    if (ugrid_unformatted) then

      if ( big_endian_io(9999) ) then
       filename_ugrid     = trim(project) // '.b8.ugrid'
      else
       filename_ugrid     = trim(project) // '.l8.ugrid'
      end if

   !(2) ASCII file
    else

       filename_ugrid     = trim(project) // '.ugrid'

    endif

 !-----------------------------------------------------------------------
 ! Input boundary condition file (ASCII file)

      filename_bcmap     = trim(project) // '.bcmap'

 !-----------------------------------------------------------------------
 ! Output: Tecplot boundary grid file (ASCII file)

      filename_tecplot_b = trim(project) // '_boundary_tec.dat'

 !-----------------------------------------------------------------------
 ! Output: Tecplot volume grid file (ASCII file)

      filename_tecplot_v = trim(project) // '_volume_tec.dat'

 !-----------------------------------------------------------------------
 ! Output: Tecplot volume grid file (ASCII file)

      filename_soln      = trim(project) // '_soln.data'

 !-----------------------------------------------------------------------

 end subroutine set_filenames


!********************************************************************************
!
! Allocate arrays required for the node-centered edge-based discretization,
! and nonlinear solvers.
!
!********************************************************************************
 subroutine allocate_nceb_arrays

 use data_module, only : nnodes, nghbrs_total, u, w, gradw, dtau, wsn, res
 use data_module, only : du, u0, jac_diag, jac_dinv, jac_off

 use input_parameter, only : solver_type

 implicit none

! Solution, gradients, and residual arrays.

  allocate(     u(nnodes,5  ) )
  allocate(    du(nnodes,5  ) )
  allocate(     w(nnodes,5  ) )
  allocate( gradw(nnodes,5,3) )
  allocate(   res(nnodes,5  ) )

  allocate(  dtau(nnodes    ) )
  allocate(   wsn(nnodes    ) )

! Extra arrays required by nonlinear solvers

  if (     trim(solver_type) == "fwd_euler" ) then

  ! No extra arrays are needed.

  elseif ( trim(solver_type) == "tvdrk" ) then

   allocate( u0(nnodes,5) )               !to save the solution for 2-stage RK.

  elseif ( trim(solver_type) == "implicit" ) then

   allocate( jac_diag(      nnodes,5,5) ) !jacobian: diagonal blocks
   allocate( jac_dinv(      nnodes,5,5) ) !jacobian: inverse of diagonal blocks
   allocate( jac_off( nghbrs_total,5,5) ) !jacobian: off-diagonal blocks 

  elseif ( trim(solver_type) == "jfnk" ) then

   allocate( jac_diag(      nnodes,5,5) ) !jacobian: diagonal blocks
   allocate( jac_dinv(      nnodes,5,5) ) !jacobian: inverse of diagonal blocks
   allocate( jac_off( nghbrs_total,5,5) ) !jacobian: off-diagonal blocks 

  else

   write(*,*) " Invalid input: solver_type = ", trim(solver_type)
   write(*,*) "   Available options: fwd_euler, tvdrk, implciit, jfnk. Stop..."
   stop

  endif

 end subroutine allocate_nceb_arrays

!********************************************************************************
!
! Set the free stream values, and initial solution.
!
!********************************************************************************
 subroutine set_initial_solution

 use data_module           , only : p2, one, zero, pi
 use data_module           , only : nnodes, u, w, bmark

 use euler_eqns     , only : rho_inf, u_inf, v_inf, w_inf, p_inf, gamma
 use euler_eqns     , only : ir, iu, iv, iw, ip, w2u

 use input_parameter, only : M_inf, aoa, perturb_initial

 implicit none

 integer  :: i
 real(p2) :: rn

  write(*,*)
  write(*,*) "---------------------------------------"
  write(*,*) " Set initial solution"
  write(*,*)

  if (perturb_initial==1) then
    write(*,*) "     perturb_initial = ", perturb_initial
    write(*,*) "  -> Add random perturbation in the y-velocity. "
  endif

! Set free stream values.
! Angle of attack (aoa) is taken from the positive x-axis to the positive z-axis.

  rho_inf = one
    u_inf = M_inf*cos(aoa*pi/180_p2)
    v_inf = zero
    w_inf = M_inf*sin(aoa*pi/180_p2)
    p_inf = one/gamma

! Let the initial solution be the free stream values.

  loop_node : do i = 1, nnodes

  !Primitive variables
   w(i,ir) = rho_inf
   w(i,iu) =   u_inf
   w(i,iv) =   v_inf
   w(i,iw) =   w_inf
   w(i,ip) =   p_inf

  !Add random perturbation in the y-velocity if perturb_initial=1.
   if (perturb_initial==1) then
    call random_number(rn)
    if ( bmark(i) == 0 ) w(i,iv) =  0.1_p2 * rn ! interior nodes only
   endif

  !Conservative variables
   u(i,:) = w2u( w(i,:) )

  end do loop_node

 end subroutine set_initial_solution

!********************************************************************************
!* Read the grid and boundary condition file.
!*
!* Note: Currently, the EDU3D-Euler works only for pure tetrahedral grids.
!*       I think you can modify the code to make it work for other elements.
!*
!********************************************************************************
!*
!* 1. "datafile_grid_in": .ugrid file name. 
!*
!*
!* 2. "datafile_bcmap_in" is assumed have been written in the following format:
!*
!*   -----------------------------------------------------------------------
!*    write(*,*) nb (# of boundaries)
!*    do i = 1, nb
!*     write(*,*) i, bc_type
!*    end do
!*   -----------------------------------------------------------------------
!*
!*   NOTE: bc_type is the name of the boundary condition, e.g.,
!*
!*         1. "freestream"
!*             Roe flux with freestream condition on the right state.
!*
!*         2. "outflow_supersonic"
!*             Roe flux with the interior state as the right state.
!*             (equivalent to the interior-extrapolation condition.)
!*
!*         3. "slip_wall"
!*             Right state with zero mass flux through the boundary.
!*
!*         4. "outflow_subsonic_p0"
!*             Fix the back pressure. This should work for subsonic flows in a
!*             large enough domain.
!*
!*         !You can implement more BCs.
!*
!*
!********************************************************************************
!* Data to be read and stored:
!*
!* 1. Some numbers
!*    nnodes        = number of nodes
!*    ntria         = number of triangular boundary elements
!*    nquad         = number of quadrilateral boundary elements
!*    ntet          = number of tetrahedra
!*    npyr          = number of pyramids
!*    nprs          = number of prims
!*    nhex          = number of hexahedra
!*
!* 2. Boundary element data:
!*    tria = list of vertices of each triangular boundary element
!*    quad = list of vertices of each quadrilateral boundary element
!*
!* 3. Volume element data:
!*    tet  = list of vertices of each tetrahedron
!*    pyr  = list of vertices of each pyramid
!*    prs  = list of vertices of each prism
!*    hex  = list of vertices of each hexehedron
!*
!* 4. x, y, z coordinates of ndoes
!*    x    = x-coordinate of the nodes
!*    y    = y-coordinate of the nodes
!*    z    = z-coordinate of the nodes
!*
!* 5. Boundary Data:
!*    nb      = number of boundary groups
!*    bc_type = boundary condition name
!*
!********************************************************************************
 subroutine read_grid

 !Input
  use input_parameter, only : ugrid_unformatted

 !Fill in these data:
  use data_module    , only : nnodes            , & !# of nodes
                              x, y, z           , & !nodal coords
                              ntet, tet         , & !elements
                              ntria, tria       , & !boundary elements
                              nb, bc_type           !# of boundaries, and BCs

 !Input file names
  use data_module    , only : filename_ugrid
  use data_module    , only : filename_bcmap

 implicit none

 !Data not used in the code at the moment (tetrahedra only).

 integer                           :: nquad
 integer , dimension(:,:), pointer :: quad

 integer                           :: nhex
 integer , dimension(:,:), pointer :: hex

 integer                           :: nprs
 integer , dimension(:,:), pointer :: prs

 integer                           :: npyr
 integer , dimension(:,:), pointer :: pyr

 integer :: i, dummy_int, os

  write(*,*)
  write(*,*) "---------------------------------------"
  write(*,*) " Reading : ", trim(filename_ugrid)
  write(*,*)

 !------------------------------------------------------------------------------
 !(1)Read an unformatted file.

  unformatted : if ( ugrid_unformatted ) then

   open(unit=2, file=filename_ugrid, form='unformatted',access="stream", &
                                               status="unknown", iostat=os )
   read(2) nnodes, ntria, nquad, ntet, npyr, nprs, nhex

  write(*,*) "          nnodes = ", nnodes
  write(*,*) "          prisms = ", nprs
  write(*,*) "           tetra = ", ntet
  write(*,*) "            hexa = ", nhex
  write(*,*) "         pyramid = ", npyr
  write(*,*) "       triangles = ", ntria
  write(*,*) "  quadrilaterals = ", nquad

                allocate(tria(ntria,4))
  if (nquad >0) allocate(quad(nquad,5))
  if (ntet  >0) allocate(  tet(ntet,4))
  if (npyr  >0) allocate(  pyr(nprs,5))
  if (nprs  >0) allocate(  prs(nprs,6))
  if (nhex  >0) allocate(  hex(nhex,8))
                allocate(x(nnodes), y(nnodes), z(nnodes))

   do i = 1, nnodes
    read(2) x(i), y(i), z(i)
   end do

  ! Triangular boundary faces
     do i = 1, ntria
      read(2) tria(i,1), tria(i,2), tria(i,3)
     end do

  ! Quad boundary faces
    if (nquad>0) then
     do i = 1, nquad
      read(2) quad(i,1), quad(i,2), quad(i,3), quad(i,4)
     end do
    endif

  ! Face tag: Boundary group number
     do i = 1, ntria
      read(2) tria(i,4)
     end do

     do i = 1, nquad
      read(2) quad(i,5)
     end do

  ! Tets
    if (ntet>0) then
     do i = 1, ntet
     read(2) tet(i,1), tet(i,2), tet(i,3), tet(i,4)
     end do
    endif

  ! Pyramids
    if (npyr>0) then
     do i = 1, npyr
     read(2) pyr(i,1), pyr(i,2), pyr(i,3), pyr(i,4), pyr(i,5)
     end do
    endif

  ! Prisms
    if (nprs>0) then
     do i = 1, nprs
     read(2) prs(i,1), prs(i,2), prs(i,3), prs(i,4), prs(i,5), prs(i,6)
     end do
    endif

  ! Hexahedra
    if (nhex>0) then
     do i = 1, nhex
     read(2) hex(i,1), hex(i,2), hex(i,3), hex(i,4), hex(i,5), hex(i,6), hex(i,7), hex(i,8)
     end do
    endif

 !------------------------------------------------------------------------------
 !(2)Read a formatted file.

  else

   open(unit=2, file=filename_ugrid, status="unknown", iostat=os)
   read(2,*) nnodes, ntria, nquad, ntet, npyr, nprs, nhex

  write(*,*) "          nnodes = ", nnodes
  write(*,*) "          prisms = ", nprs
  write(*,*) "           tetra = ", ntet
  write(*,*) "            hexa = ", nhex
  write(*,*) "         pyramid = ", npyr
  write(*,*) "       triangles = ", ntria
  write(*,*) "  quadrilaterals = ", nquad

                allocate(tria(ntria,4))
  if (nquad >0) allocate(quad(nquad,5))
  if (ntet  >0) allocate(  tet(ntet,4))
  if (npyr  >0) allocate(  pyr(nprs,5))
  if (nprs  >0) allocate(  prs(nprs,6))
  if (nhex  >0) allocate(  hex(nhex,8))
                allocate(x(nnodes), y(nnodes), z(nnodes))

   do i = 1, nnodes
    read(2,*) x(i), y(i), z(i)
   end do

  ! Triangular boundary faces
     do i = 1, ntria
      read(2,*) tria(i,1), tria(i,2), tria(i,3)
     end do

  ! Quad boundary faces
    if (nquad>0) then
     do i = 1, nquad
      read(2,*) quad(i,1), quad(i,2), quad(i,3), quad(i,4)
     end do
    endif

  ! Face tag: Boundary group number
     do i = 1, ntria
      read(2,*) tria(i,4)
     end do

     do i = 1, nquad
      read(2,*) quad(i,5)
     end do

  ! Tets
    if (ntet>0) then
     do i = 1, ntet
     read(2,*) tet(i,1), tet(i,2), tet(i,3), tet(i,4)
     end do
    endif

  ! Pyramids
    if (npyr>0) then
     do i = 1, npyr
     read(2,*) pyr(i,1), pyr(i,2), pyr(i,3), pyr(i,4), pyr(i,5)
     end do
    endif

  ! Prisms
    if (nprs>0) then
     do i = 1, nprs
     read(2,*) prs(i,1), prs(i,2), prs(i,3), prs(i,4), prs(i,5), prs(i,6)
     end do
    endif

  ! Hexahedra
    if (nhex>0) then
     do i = 1, nhex
     read(2,*) hex(i,1), hex(i,2), hex(i,3), hex(i,4), hex(i,5), hex(i,6), hex(i,7), hex(i,8)
     end do
    endif

  endif unformatted
 !------------------------------------------------------------------------------

  close(2)

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! 2. Read the boundary condition data file

   write(*,*)
   write(*,*) "---------------------------------------"
   write(*,*) " Reading the boundary condition file....", trim(filename_bcmap)
   write(*,*)

! Open the input file.
  open(unit=2, file=filename_bcmap, status="unknown", iostat=os)

    read(2,*) nb

  if (nquad == 0) then
   allocate(quad(1,5))
   quad = 0
  endif

  if ( nb /= max( maxval(tria(:,4)), maxval(quad(:,5)) ) ) then
   write(*,*) " Error in bcmap file..." 
   stop
  endif

  allocate(bc_type(nb))

! READ: Read the boundary condition type. bc_type is character(80).
   do i = 1, nb
    read(2,*) dummy_int, bc_type(i)
   end do

!  Print the data
    write(*,*)
   do i = 1, nb
    write(*,'(a10,i3,a12,a35)') " boundary", i, "  bc_type = ", trim(bc_type(i))
   end do

    write(*,*)

  close(2)

! End of Read the boundary condition data file
!--------------------------------------------------------------------------------

 end subroutine read_grid

!*******************************************************************************
!
! This subroutine writes a Tecplot file for boundaries.
!
!*******************************************************************************
 subroutine write_tecplot_boundary_file

  use data_module      , only : x, y, z     , & !nodal coords
                                w           , & !primitive variables
                                ntria, tria , & !boundary elements
                                nnodes      , & !number of nodes
                                nb          , & !number of boundaries
                                bc_type         !BC type name

 use euler_eqns , only : ir, iu, iv, iw, ip

! Output tecplot boundary filename
  use data_module               , only : filename_tecplot_b, bmark, bmarks


 Implicit none

!Input

 character(80)                         :: boundary_name
 integer                               :: nelms
 integer                               :: nnodes_loc, i, j, k, ib, os
 integer,  dimension(:,:), allocatable :: g2l
 integer,  dimension(:  ), allocatable :: l2g

 write(*,*)
 write(*,*) ' Tecplot boundary file = ', trim(filename_tecplot_b)
 write(*,*)

!Allocate global-to-local and local-to-global arrays.

 allocate(g2l(nnodes,2))
 allocate(l2g(nnodes  ))

 open(unit=7, file=filename_tecplot_b, status="unknown", iostat=os)
 write(7,*) 'TITLE = "GRID"'
 write(7,*) 'VARIABLES = "x","y","z","rho","u","v","w","p","bmark","bmarks"'

!Write boundaries separately.

 tria_boundary_loop : do ib = 1, nb

!---------------------------------------------------------------------------
! Triangular boundaries.

       nelms = 0
  nnodes_loc = 0
         g2l = -1
         l2g = 0

 !To save memory, we write the boundary data based on the data
 !local to each boundary. Local numbers are assigned to each
 !node on the boundary, which is accessed by g2l (global to local).

  do i = 1, ntria

   if ( tria(i,4) == ib ) then

    nelms = nelms + 1

    do k = 1, 3
     if (g2l(tria(i,k),2) == -1) then
      nnodes_loc = nnodes_loc + 1
      g2l(tria(i,k),1) = nnodes_loc
      g2l(tria(i,k),2) = 100         !Flag to indicate the node is recorded.
      l2g(nnodes_loc)  = tria(i,k)
     endif
    end do

   endif

  end do

 !So, there are nnodes_loc nodes and nelms elements on this boundary (ib).

  if (nelms > 0) then
   boundary_name = trim(bc_type(ib)) // " (triangles)"
   write(7,*) 'ZONE T="', trim(boundary_name), '"  N=', nnodes_loc, &
              ', E=', nelms, ' , ET=quadrilateral, F=FEPOINT'

  !Write only the local nodes.
   do j = 1, nnodes_loc
     i = l2g(j)
     write(7,'(8es25.15,2i5)') x(i), y(i), z(i), &
                           w(i,ir),w(i,iu),w(i,iv),w(i,iw),w(i,ip), bmark(i),bmarks(i)

   end do

  !Write only the local elements.
   do i = 1, ntria
    if ( tria(i,4) == ib ) write(7,'(4i10)') g2l(tria(i,1),1), g2l(tria(i,2),1), &
                                             g2l(tria(i,3),1), g2l(tria(i,3),1)
   end do
  endif

 end do tria_boundary_loop

!---------------------------------------------------------------------------

 close(7)

 end subroutine write_tecplot_boundary_file
!********************************************************************************

!*******************************************************************************
! This subroutine writes a Tecplot file for the volume grid.
! In this program, only tetrahedral grids are considered.
!*******************************************************************************
 subroutine write_tecplot_volume_file

  use data_module       , only : x, y, z     , & !nodal coords
                                 w           , & !primitive variables
                                 ntet, tet   , & !number of tetra, and tetra info
                                 nnodes          !number of nodes

  use euler_eqns , only : ir, iu, iv, iw, ip

! Output tecplot volume filename
  use data_module               , only : filename_tecplot_v, bmark, bmarks

  implicit none

  integer :: i, os

 write(*,*)
 write(*,*) ' Tecplot volume file = ', trim(filename_tecplot_v)
 write(*,*)

  open(unit=8, file=filename_tecplot_v, status="unknown", iostat=os)
  write(8,*) 'TITLE = "GRID"'
  write(8,*) 'VARIABLES = "x","y","z","rho","u","v","w","p","bmark","bmarks"'

!---------------------------------------------------------------------------
! Tetrahedral zone

  if (ntet > 0) then

   write(8,*) 'zone  n=', nnodes,',e=', ntet,' , et=tetrahedron, f=fepoint'
   do i = 1, nnodes
     write(8,'(8es25.15,2i5)') x(i), y(i), z(i), &
                           w(i,ir),w(i,iu),w(i,iv),w(i,iw),w(i,ip), bmark(i),bmarks(i)
   end do

   do i = 1, ntet
    write(8,'(4i10)') tet(i,1), tet(i,2), tet(i,3), tet(i,4)
   end do

  endif

!---------------------------------------------------------------------------

 close(8)

 end subroutine write_tecplot_volume_file
!********************************************************************************


!*******************************************************************************
! This subroutine writes a solution data file, containing, node numbers, and
! solution values.
!*******************************************************************************
 subroutine write_solution_file

  use data_module       , only : x, y, z     , & !nodal coords
                                 w           , & !primitive variables
                                 nnodes          !number of nodes

! Output tecplot volume filename
  use data_module               , only : filename_soln

  implicit none

  integer :: i, os

 write(*,*)
 write(*,*) ' Tecplot volume file = ', trim(filename_soln)
 write(*,*)

  open(unit=8, file=filename_soln, status="unknown", iostat=os)

!---------------------------------------------------------------------------
! Write soution data.

     write(8,*) nnodes

   do i = 1, nnodes
     write(8,'(8es25.15)') x(i), y(i), z(i), w(i,1:5)
   end do

!---------------------------------------------------------------------------

 close(8)

 end subroutine write_solution_file
!********************************************************************************


!********************************************************************************
! Find out big_endian_io.
!********************************************************************************
 function big_endian_io( opt_unit )

 integer, intent(in) :: opt_unit
 logical             :: big_endian_io

! one-byte integer(kind)
 integer, parameter  :: i1 = selected_int_kind(2)

! two-byte integer(kind)
 integer, parameter  :: i2 = selected_int_kind(4)
 integer(i1)         :: byte_one, byte_two

! 00000000 00000001 big-endian binary
 integer(i2)         :: two_byte_int = 1_i2

    open(opt_unit,status='scratch',form='unformatted')
      write( opt_unit) two_byte_int
      rewind(opt_unit)
      read(  opt_unit) byte_one, byte_two
    close(opt_unit)
    big_endian_io = ( byte_one == 0 .and. byte_two == 1 )

 end function big_endian_io
!********************************************************************************




!********************************************************************************
!*
!* This subroutine performs a numerical TE analysis to verify the order of
!* accuracy of the truncation error.
!*
!********************************************************************************
 subroutine perform_te_test_with_mms

  use input_parameter, only : generate_tec_file_b,  generate_tec_file_v
  use input_parameter, only : project, ugrid_unformatted, mms_s_quadrature, kappa_s

  use data_module    , only : filename_ugrid    , filename_bcmap
  use data_module    , only : filename_tecplot_b, filename_tecplot_v

  use nceb_data      , only : construct_nceb_data
  use residual       , only : compute_residual

  use data_module    , only : nnodes, res, vol, x, y, z, res, w, u, bmarks
  use data_module    , only : edge, nedges, n12, n12_mag, e12, e12_mag

  use edu3d_mms_euler, only : compute_manufactured_sol_and_source
 implicit none
 integer , parameter :: p2 = selected_real_kind(p=15)

 real(p2), parameter ::       zero = 0.0_p2
 real(p2), parameter ::       half = 0.5_p2
 real(p2), parameter ::        one = 1.0_p2
 real(p2), parameter ::      three = 3.0_p2

 real(p2)     , dimension(  10) :: heff
 real(p2)     , dimension(5,10) :: norm1, ooa

 integer                          :: i, k, n, igrid, ngrids
 character(80)                :: i_char
 character(80), dimension(10)   :: mms_grid_filename

 integer                         :: node1, node2
 real(p2), dimension(3) :: da    , ev
 real(p2)                        :: da_mag, ev_mag, vjk
 real(p2), dimension(5) :: s1, sx1, sy1, sz1, ds1
 real(p2), dimension(5) :: s2, sx2, sy2, sz2, ds2
 
 real(p2) , dimension(:,:), allocatable  :: forcing, sx, sy, sz

  write(*,*)
  write(*,*) " --- MMS Truncation Error Test ---"
  write(*,*)

 !--------------------------------------------------------------------
 ! Boundary condition file: all must be "mms_dirichlet":

   filename_bcmap = trim(project) // '.bcmap'

 !--------------------------------------------------------------------
 ! Define grid files to be read.
 
   ngrids = 3 !<- 3 grids will be used.

   mms_grid_filename(1) = "./grids/tetgrid_08x08x08.ugrid" !1st grid
   mms_grid_filename(2) = "./grids/tetgrid_16x16x16.ugrid" !2nd grid
   mms_grid_filename(3) = "./grids/tetgrid_32x32x32.ugrid" !3rd grid

 !--------------------------------------------------------------------
 ! Loop over grids
 !--------------------------------------------------------------------
  do igrid = 1, ngrids

   !Filename for a boundary grid Tecplot file.
    write(i_char,'(i0)') igrid
    filename_tecplot_b    = trim(project) // '_boundary_tec' // trim(i_char) // '.dat'
    filename_tecplot_v    = trim(project) // '_volume_tec'   // trim(i_char) // '.dat'


   !Filename for igrid-th grid.
    filename_ugrid = mms_grid_filename(igrid)

    write(*,*) " Grid ", igrid, trim(filename_ugrid)

    call read_grid             !Read a grid.
    call construct_nceb_data   !Construct grid quantities for NCEB.
    call allocate_nceb_arrays  !Allocate solution arrays and others.

  !------------------------------------------------------------
  ! Compute the residual with the exact soltuion.

    allocate( forcing(nnodes,5) )
    allocate(      sx(nnodes,5) )
    allocate(      sy(nnodes,5) )
    allocate(      sz(nnodes,5) )

   !Compute and store the source and its derivativs.

    do i = 1, nnodes
     call compute_manufactured_sol_and_source(  &
             x(i),y(i),z(i), w(i,:),u(i,:),forcing(i,:),sx(i,:),sy(i,:),sz(i,:) )
    end do

   !Compute the spatial residual of the Euler.

    call compute_residual
    
   !Add the source term: two formulas.

    !--------------------------------------------------------------------------------
    !(1) Simple point evaluation

     if     (trim(mms_s_quadrature)=="point") then

      do i = 1, nnodes
       res(i,:) = res(i,:) - forcing(i,:)*vol(i)
      end do

    !--------------------------------------------------------------------------------
    !(2) Special source quadrature

     elseif (trim(mms_s_quadrature)=="kappa_s") then
 
      write(*,*) "nedges=", nedges 
      
      do i = 1, nedges

         node1 = edge(i,1)   ! Left node of the edge
         node2 = edge(i,2)   ! Right node of the edge
               da = n12(i,:)    ! This is the directed area vector (unit vector)
       da_mag = n12_mag(i)  ! Magnitude of the directed area vector
                ev = e12(i,:)    ! Edge vector (unit vector)
       ev_mag = e12_mag(i)  ! Magnitude of the edge vector

           s1  = forcing(node1,:)
           sx1 =      sx(node1,:)
           sy1 =      sy(node1,:)
           sz1 =      sz(node1,:)

           s2  = forcing(node2,:)
           sx2 =      sx(node2,:)
           sy2 =      sy(node2,:)
           sz2 =      sz(node2,:)
           
           vjk = ( ev(1)*da(1) + ev(2)*da(2) + ev(3)*da(3) )*da_mag*ev_mag/6.0_p2

                ds1 =  half*( sx1*ev(1) + sy1*ev(2) + sz1*ev(3) )*ev_mag
       res(node1,:) = res(node1,:) - ( kappa_s*half*(s1+s2)+(one-kappa_s)*(s1+ds1) )*vjk
       
                ds2 = -half*( sx2*ev(1) + sy2*ev(2) + sz2*ev(3) )*ev_mag
       res(node2,:) = res(node2,:) - ( kappa_s*half*(s2+s1)+(one-kappa_s)*(s2+ds2) )*vjk

      end do

    !--------------------------------------------------------------------------------
     else

      write(*,*) " Invalid input: mms_s_quadrature = ", trim(mms_s_quadrature)
      stop

     endif
    !--------------------------------------------------------------------------------
    
    deallocate( forcing, sx, sy, sz )
 
  !------------------------------------------------------------
  ! Compute the residual norm: res/vol = truncation error.

                 n = 0

      heff( igrid) = zero
    norm1(:,igrid) = zero

    do i = 1, nnodes

     !Skip boundary nodes,     their neighbors
     if (  bmarks(i) == 1 .or. bmarks(i) == 2 ) cycle

                  n = n + 1

        heff(igrid) = heff(igrid) + ( vol(i) )**(one/three)
     norm1(:,igrid) = norm1(:,igrid) + abs( res(i,:)/vol(i) )

    end do 

        heff(igrid) =    heff(igrid)/real(n,p2)
     norm1(:,igrid) = norm1(:,igrid)/real(n,p2)

  !------------------------------------------------------------
  ! Write a boundary grid Tecplot file if requested.

   if (generate_tec_file_b) call write_tecplot_boundary_file
   if (generate_tec_file_v) call write_tecplot_volume_file

   call deallocate_all

  end do

 !--------------------------------------------------------------------
 ! End of Loop over grids
 !--------------------------------------------------------------------


 !------------------------------------------------------------
 ! Print the results
 !------------------------------------------------------------

   do igrid = 2, ngrids
    ooa(:,igrid) = ( log10(norm1(:,igrid)) - log10(norm1(:,igrid-1)) ) / ( log10(heff(igrid)) - log10(heff(igrid-1)) ) 
   end do


    write(*,*)
    write(*,*) " ------------------------------------------------------------------------ "
    write(*,*) " ----- Results of TE analysis: "
    write(*,*) " ------------------------------------------------------------------------ "
    write(*,*) "                                Res(1)    Res(2)    Res(3)    Res(4)    Res(5)     heff"

   do igrid = 1, ngrids

    write(*,'(a,i5,a,6es10.3)') " Grid level = ", igrid, ": L1(TE) = ", ( norm1(k,igrid), k=1,5), heff(igrid)

   end do

    write(*,*)
    
   do igrid = 2, ngrids

    write(*,'(a,i5,a,5es10.3)') " Grid level = ", igrid, ": Order = ", ( ooa(k,igrid), k=1,5 ) 

   end do

    write(*,*)

 !------------------------------------------------------------
 ! End of Print the results
 !------------------------------------------------------------


 end subroutine perform_te_test_with_mms

!********************************************************************************
!
! Allocate arrays required for the node-centered edge-based discretization,
! and nonlinear solvers.
!
!********************************************************************************
 subroutine deallocate_all

 use data_module    , only : u, w, gradw, dtau, wsn, res, bc_type, nghbr, k0, nnghbrs
 use data_module    , only : du, u0, jac_diag, jac_dinv, jac_off, edge, kth_nghbr_of_1, kth_nghbr_of_2, vol
 use data_module    , only : n12, n12_mag, e12, e12_mag, bmark, bmarks, bn, bn_mag, bnj, bnj_mag
 use input_parameter, only : solver_type
 use data_module    , only : x, y, z, tet, tria
 use data_module    , only : linear_lsq3x3_cx, linear_lsq3x3_cy, linear_lsq3x3_cz

 implicit none

! Arrays allocated in "subroutine read_grid".

  deallocate( x, y, z, tet, tria, bc_type )

! Arrays for the LSQ method

  deallocate( linear_lsq3x3_cx, linear_lsq3x3_cy, linear_lsq3x3_cz )

! Arrays allocated in "subroutine construct_nceb_data".

  deallocate( u,du, w,gradw, res, dtau,wsn, nghbr, k0, nnghbrs, edge, kth_nghbr_of_1, kth_nghbr_of_2, vol )
  deallocate( n12, n12_mag, e12, e12_mag, bmark, bmarks, bn, bn_mag, bnj, bnj_mag )

  ! Extra arrays required by nonlinear solvers

  if (     trim(solver_type) == "fwd_euler" ) then

  ! No extra arrays are needed.

  elseif ( trim(solver_type) == "tvdrk" ) then

   deallocate( u0 )

  elseif ( trim(solver_type) == "implicit" ) then

   deallocate( jac_diag )
   deallocate( jac_dinv )
   deallocate( jac_off  )

  elseif ( trim(solver_type) == "jfnk" ) then

   deallocate( jac_diag )
   deallocate( jac_dinv )
   deallocate( jac_off  )

  else

   write(*,*) " Invalid input: solver_type = ", trim(solver_type)
   write(*,*) "   Available options: fwd_euler, tvdrk, implciit, jfnk. Stop..."
   stop

  endif

 end subroutine deallocate_all

 end program edu3d_euler


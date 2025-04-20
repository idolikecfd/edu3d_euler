!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes grid, solution, and other data used throughout
!* the whole program. These data are accessed by the use statement.
!*
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* 
!* This F90 program is written and made available for an educational purpose.
!*
!* This file may be updated in future.
!*
!* Katate Masatsuka, October 2017. http://www.cfdbooks.com
!********************************************************************************

 module data_module

  implicit none

  public
 
 !--------------------------------------------------------------------
 ! Constants

  integer , parameter :: p2 = selected_real_kind(P=15) !Double precision
  real(p2), parameter ::   zero = 0.0_p2
  real(p2), parameter ::   half = 0.5_p2
  real(p2), parameter ::    one = 1.0_p2
  real(p2), parameter ::    two = 2.0_p2
  real(p2), parameter ::  three = 3.0_p2
  real(p2), parameter ::   four = 4.0_p2
  real(p2), parameter ::   five = 5.0_p2
  real(p2), parameter ::  third = 1.0_p2/3.0_p2
  real(p2), parameter :: fourth = 1.0_p2/4.0_p2
  real(p2), parameter ::  sixth = 1.0_p2/6.0_p2
  real(p2), parameter :: my_eps = epsilon(one)

  real(p2), parameter :: pi = 3.141592653589793238_p2

  integer, parameter :: ix = 1, iy = 2, iz = 3

 !--------------------------------------------------------------------
 !--------------------------------------------------------------------
 ! file names

  character(80) :: filename_ugrid     ! input grid filename (.ugrid)
  character(80) :: filename_mapbc     ! input bc   filename (.mapbc)
  character(80) :: filename_tecplot_b ! output tecplot boundary filename (.dat)
  character(80) :: filename_tecplot_v ! output tecplot volume   filename (.dat)
  character(80) :: filename_soln      ! output solution filename (.data)

 !--------------------------------------------------------------------
 ! Data used in the edge-based discretization

  !------------------------------------------
  !>> Node data
  integer                                 :: nnodes
  real(p2), dimension(:  )  , allocatable :: x, y, z
  real(p2), dimension(:  )  , allocatable :: vol

  real(p2), dimension(:  )  , allocatable :: dtau
  real(p2), dimension(:  )  , allocatable :: wsn

  real(p2), dimension(:,: ) , allocatable :: u
  real(p2), dimension(:,:)  , allocatable :: w
  real(p2), dimension(:,:)  , allocatable :: du
  real(p2), dimension(:,:)  , allocatable :: res
  real(p2), dimension(:,:,:), allocatable :: gradw

  integer                                 :: nghbrs_total
  integer , dimension(:    ), allocatable :: k0
  integer , dimension(:    ), allocatable :: nnghbrs
  integer , dimension(:    ), allocatable :: nghbr    !i=1,nnodes*nnghbrs

  integer , dimension(:    ), allocatable :: bmark  !Boundary node mark: 0 interior, >0 boundary.
  integer , dimension(:    ), allocatable :: bmarks !Boundary node mark

  !LSQ gradient coefficients
  real(p2), dimension(:)    , allocatable :: linear_lsq3x3_cx
  real(p2), dimension(:)    , allocatable :: linear_lsq3x3_cy
  real(p2), dimension(:)    , allocatable :: linear_lsq3x3_cz

  !------------------------------------------
  !>> Edge data
  integer                                 :: nedges
  integer , dimension(:,:)  , allocatable :: edge
  real(p2), dimension(:,:)  , allocatable :: n12
  real(p2), dimension(:  )  , allocatable :: n12_mag
  real(p2), dimension(:,:)  , allocatable :: e12
  real(p2), dimension(:  )  , allocatable :: e12_mag
  integer , dimension(:    ), allocatable :: kth_nghbr_of_1
  integer , dimension(:    ), allocatable :: kth_nghbr_of_2

  !------------------------------------------
  !>> Boundary element data
  integer                                 :: ntria
  integer , dimension(:,:)  , allocatable ::  tria
  real(p2), dimension(:,:)  , allocatable ::  bn      !Unit normal of the tria face.
  real(p2), dimension(:  )  , allocatable ::  bn_mag  !Magnitude of bn
  real(p2), dimension(:,:,:), allocatable ::  bnj     !Unit normal at nodes of the tria.
  real(p2), dimension(:,:)  , allocatable ::  bnj_mag !Magnitude of bnj

  integer                                  :: nb
  character(80), dimension(:), allocatable :: bc_type !type of boundary condition

 !--------------------------------------------------------------------
 !--------------------------------------------------------------------
 ! Temporary data

  !>> Element data
  integer                               :: ntet
  integer , dimension(:,:), allocatable ::  tet

 !--------------------------------------------------------------------
 !--------------------------------------------------------------------
 ! Steady solver

  integer  :: i_iteration

  real(p2), dimension(:,: ) , allocatable :: u0
  real(p2), dimension(:,:,:), allocatable :: jac_diag !i=1,nnodes
  real(p2), dimension(:,:,:), allocatable :: jac_off  !i=1,nnodes*nnghbrs
  real(p2), dimension(:,:,:), allocatable :: jac_dinv !i=1,nnodes

  integer  :: n_residual_evaluation !#of residual evaluations


 !--------------------------------------------------------------------
 !--------------------------------------------------------------------
 ! Linear solver

  integer  :: lrelax_sweeps_actual
  real(p2) :: lrelax_roc

  integer  :: gcr_proj_actual
  integer  :: status_pc
  integer  :: status_gcr
  real(p2) :: gcr_rkr0

 !--------------------------------------------------------------------
 !--------------------------------------------------------------------


 end module data_module


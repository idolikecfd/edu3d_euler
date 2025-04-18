!********************************************************************************
! Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
! This module defines input parameters, and set the default values.
! These parameters are specified in the file named 'input.nml', and read
! in the main program. 
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! 
! This F90 program is written and made available for an educational purpose.
!
! This file may be updated in future.
!
! Katate Masatsuka, October 2017. http://www.cfdbooks.com
!********************************************************************************

 module input_parameter

  implicit none

 !This module contains the following subroutine:

  public :: read_nml_input_parameters ! read the input file
  public

  integer , parameter ::     p2 = selected_real_kind(P=15)

!----------------------------
! Default input values
!----------------------------

  character(80) :: project             = "default"    ! project name

  logical       :: test_mms_te         = .false.      ! to preform TE test to verify the residual.

  real(p2)      :: M_inf               = 0.3_p2       ! free stream Mach number to positive x.

  real(p2)      :: aoa                 = 0.0_p2       ! free stream angle of attack: x -> z

  character(80) :: inviscid_flux       = "roe"        ! numerical flux: roe, hll, rusanov, rhll

  character(80) :: inviscid_jac        = "roe"        ! flux Jacobian : roe, hll, rusanov, rhll

  integer       :: accuracy_order      = 1            ! 1=1st-order, 2=2nd-order

  real(p2)      :: lsq_weight_power    = 0.0_p2       ! lsq weight: 1/distance^lsq_weight_power

  real(p2)      ::            kappa    = 0.0_p2       ! U-MUSCL parameter.

  integer       :: perturb_initial     = 1            ! 0=unperturbed, 1=perturbed

  character(80) :: mms_s_quadrature    = "point"      ! point or kappa: source quadrature
  real(p2)      ::          kappa_s    = 0.0_p2       ! parameter in kappa-source quadrature
                                                      ! kappa_s=3/10 for 3rd-order UMUSCL-SSQ

! Nonlinear solver (steady solver)
  character(80) :: solver_type         = "implicit"   ! solver = "fwd_euler"
                                                      !          "tvdrk"
                                                      !          "implicit"
                                                      !          "jfnk"

  real(p2)      :: solver_tolerance    = 1.0e-10_p2   ! solver: tolerance for nonlinear solve

  integer       :: solver_max_itr      = 10000        ! solver: max nonlinear iterations

  real(p2)      :: CFL                 = 0.9_p2       ! CFL for pseudo time

! JFNK-GCR
  integer       :: gcr_projection      = 10           ! max GCR projections
  real(p2)      :: gcr_mu              = 0.1_p2       ! tolerance for GCR (reduction)

! Linear relaxation (preconditioner)
  character(80) :: lrelax_scheme       = "gs"    ! preconditioner scheme type

  integer       :: lrelax_sweeps       = 500     ! preconditioner max relaxation

  real(p2)      :: lrelax_tolerance    = 0.1_p2  ! preconditioner tolerance (reduction)

! Data files
  logical       :: ugrid_unformatted   = .true.  ! ugrid_file_unformatted = T: unformatted, F: formatted
  logical       :: generate_tec_file_b = .true.  ! generate_tec_file_b = T
  logical       :: generate_tec_file_v = .false. ! generate_tec_file_v = T
  logical       :: generate_soln_file  = .false. !
  logical       :: read_soln_file      = .false. !

!----------------------------
! End of Default input values
!----------------------------

 namelist / input_parameters /   &
!
   project            , &
   test_mms_te        , &
   M_inf              , &
   aoa                , &
   inviscid_flux      , &
   inviscid_jac       , &
   accuracy_order     , &
   lsq_weight_power   , &
   kappa              , &
   perturb_initial    , &
   mms_s_quadrature   , &
   kappa_s            , &
   solver_type        , &
   solver_tolerance   , &
   solver_max_itr     , &
   CFL                , &
   gcr_projection     , &
   gcr_mu             , &
   lrelax_scheme      , &
   lrelax_sweeps      , &
   lrelax_tolerance   , & 
   ugrid_unformatted  , &
   generate_tec_file_b, &
   generate_tec_file_v, &
   generate_soln_file , &
   read_soln_file 

 contains

!*****************************************************************************
!
! This subroutine reads input_parameters specified in the input file:
!
!    file name = namelist_file
!
! prints the content on screen.
!
! In the main program, we set: namelist_file = "input.nml"
!
!*****************************************************************************

  subroutine read_nml_input_parameters(namelist_file)

  implicit none
  character(9), intent(in) :: namelist_file
  integer :: os

  write(*,*) "**************************************************************"
  write(*,*) " List of namelist variables and their values"
  write(*,*)

  open(unit=10,file=trim(namelist_file),form='formatted',status='old',iostat=os)
  read(unit=10,nml=input_parameters)

  write(*,nml=input_parameters) ! Print the namelist variables.
  close(10)

  end subroutine read_nml_input_parameters
!*****************************************************************************

 end module input_parameter
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!  End of input parameter module

!********************************************************************************
! Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
! This module containes subroutines that compute the right state based on
! various boundary conditions. This provides a boundary state for computing
! a boundary flux for implementing weak boundary conditions.
!
! All subroutines are written in terms of the derivative data type
! (automatic differentiation). So, it returns the right state as well as
! the derivative dwR/dwj, where wR is the right state and wj is the interior
! state.
!
! Note: In the node-centered edge-based scheme, the boundary flux is computed
!       at a boundary node. So, we compute a flux with the nodal value, and
!       a right state at the same location. The right state is determined by
!       a boundary condition. That is exactly what is determined here.
!
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

 module bc_states

  implicit none

 !This module contains the following subroutines:

  public :: get_right_state     !get the right state by bc
  public :: dw1dw2_to_du1du2    !transformation
  public :: freestream          !use free stream values
  public :: outflow_subsonic_p0 !fix the pressure value
  public :: slip_wall           !eliminate the normal velocity component
  public :: outflow_supersonic  !supersonic outflow (use the interior state)

 contains

!********************************************************************************
! Compute the right state that incorporates a boundary condition.
!
! - Compute wR for given wj.
!
! - nj is the unit normal at the node j, which is used for slip condition.
!
! - This subroutines return both wR and duR/duj. The latter is requried to
!   compute the residual Jacobian exactly.
!
!********************************************************************************

 subroutine get_right_state(bc_state_type,wj,nj,x,y,z,  duRduj,wR)

 use edu3d_mms_euler, only : compute_manufactured_sol_and_source

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 real(p2), dimension(   5), intent( in) :: wj
 real(p2), dimension(   3), intent( in) :: nj
 real(p2),                  intent( in) :: x, y, z
 character(80),             intent( in) :: bc_state_type

!Output
 real(p2), dimension(5,5), intent(out) :: duRduj
 real(p2), dimension(  5), intent(out) :: wR

!Local variables
 real(p2), dimension(  5) :: wL
 real(p2), dimension(5,5) :: dwRdwL

 real(p2), dimension(  5) :: dummy5,dummy5x,dummy5y,dummy5z

  wL = wj

!-----------------------------------------------------------------------
! Below, input is wL = the primitive variabes [rho,u,v,p].
! Output is the right state in wR = [rho,u,v,p] and the Jacobian duR/duL.

  select case(trim(bc_state_type))

   case('freestream')

     call freestream(wL, dwRdwL,wR)

   case('outflow_subsonic_p0')

     call outflow_subsonic_p0(wL, dwRdwL,wR)

   case('slip_wall')

     call slip_wall(wL,nj, dwRdwL,wR)

   case('outflow_supersonic')

     call outflow_supersonic(wL, dwRdwL,wR)

   case('mms_dirichlet')

     call compute_manufactured_sol_and_source(x,y,z, wR,wR,dummy5,dummy5x,dummy5y,dummy5z )

   case default

     write(*,*) "Boundary condition=",trim(bc_state_type),"  not implemented."
     write(*,*) " --- Stop at get_right_state() in bc_states.f90..."
     stop

  end select
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Output 1: Return the derviative duR/duL = (duRdwR)*dwRdwL*(dwLduL)

  duRduj = dw1dw2_to_du1du2(wR,wL,dwRdwL)

!-----------------------------------------------------------------------
! Output 2: Return the right state in the primitive variables: [rho,u,v,w,p]

      wR = wR

 end subroutine get_right_state

!********************************************************************************
!* Transformation: du1/du2 = (du1/dw1)*dw1/dw2*(dw2/du2)
!*
!*  Input: w1,w2,dw1dw2
!* Output:       du1du2
!*
!* Note: w  = [rho,     u,     v,     p, x,x,x,...]
!* Note: u  = [rho, rho*u, rho*v, rho*E, x,x,x,...]
!********************************************************************************
 function dw1dw2_to_du1du2(w1,w2,dw1dw2) result(du1du2)

 implicit none

 integer , parameter ::    p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter ::  zero = 0.0_p2
 real(p2), parameter ::  half = 0.5_p2
 real(p2), parameter ::   one = 1.0_p2
 real(p2), parameter :: gamma = 1.4_p2

 real(p2), dimension(5)  , intent( in) :: w1, w2
 real(p2), dimension(5,5), intent( in) :: dw1dw2

 real(p2), dimension(5,5)              :: du1du2 !output

 real(p2), dimension(5,5) :: du1dw1, dw2du2
 real(p2) :: q2

!-------------------------------------------------------------
! Primitive variables = [rho, u, v, p]

! 1. Compute du1/dw1.

            q2 = w1(2)*w1(2) + w1(3)*w1(3) + w1(4)*w1(4)
   du1dw1      = zero

   du1dw1(1,1) =  one

   du1dw1(2,1) =   w1(2); du1dw1(2,2) = w1(1);
   du1dw1(3,1) =   w1(3); du1dw1(3,3) = w1(1);
   du1dw1(4,1) =   w1(4); du1dw1(4,4) = w1(1);

   du1dw1(5,1) = half*q2;
   du1dw1(5,2) = w1(1)*w1(2); 
   du1dw1(5,3) = w1(1)*w1(3); 
   du1dw1(5,4) = w1(1)*w1(4); 
   du1dw1(5,5) = one/(gamma-one);

! 2. Compute dw2/du2.

            q2 = w2(2)*w2(2) + w2(3)*w2(3) + w2(4)*w2(4)
   dw2du2      = zero

   dw2du2(1,1) = one

   dw2du2(2,1) = -w2(2)/w2(1); dw2du2(2,2) = one/w2(1); 
   dw2du2(3,1) = -w2(3)/w2(1); dw2du2(3,3) = one/w2(1); 
   dw2du2(4,1) = -w2(4)/w2(1); dw2du2(4,4) = one/w2(1); 

   dw2du2(5,1) =  (gamma-one)*half*q2;
   dw2du2(5,2) = -(gamma-one)*w2(2)
   dw2du2(5,3) = -(gamma-one)*w2(3)
   dw2du2(5,4) = -(gamma-one)*w2(4)
   dw2du2(5,5) =  (gamma-one)

! 3. Compute du1/du2 = (du1/dw1)*dw1/dw2*(dw2/du2).

   du1du2 = matmul(du1dw1,matmul(dw1dw2,dw2du2))

 end function dw1dw2_to_du1du2

!

!********************************************************************************
! Freesteam (Compressible)
!
!********************************************************************************
 subroutine freestream(wL_real,dwRdwL,wR_real)

 use derivative_data_df5
 use euler_eqns  , only : rho_inf, u_inf, v_inf, w_inf, p_inf

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

 real(p2)                      , dimension(5)  , intent( in) :: wL_real
 real(p2)                      , dimension(5,5), intent(out) :: dwRdwL
 real(p2)                      , dimension(5)  , intent(out) :: wR_real

 type(derivative_data_type_df5), dimension(5) :: wL, wR

 integer :: i, j

  wL = wL_real
  call ddt_seed(wL)

   wR(1)    = rho_inf
   wR(2)    =   u_inf
   wR(3)    =   v_inf
   wR(4)    =   w_inf
   wR(5)    =   p_inf

!-------------------------------------------------------------------------
!-------------------------
!--------------
! Output

 !(1)Right state

     wR_real = wR%f

 !(2)Derivative

    do i = 1, 5
     do j = 1, 5
      dwRdwL(i,j) = wR(i)%df(j)
     end do
    end do

 end subroutine freestream


!********************************************************************************
! Back Pressure (Compressible): Just fix pressure.
!
!********************************************************************************
 subroutine outflow_subsonic_p0(wL_real,dwRdwL,wR_real)

 use derivative_data_df5
 use euler_eqns  , only : p_inf

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

 real(p2)                      , dimension(5)  , intent( in) :: wL_real
 real(p2)                      , dimension(5,5), intent(out) :: dwRdwL
 real(p2)                      , dimension(5)  , intent(out) :: wR_real

 type(derivative_data_type_df5), dimension(5) :: wL, wR
 integer                                      :: i, j, ip

   ip = 5
   wL = wL_real

   call ddt_seed(wL)

   wR     = wL
   wR(ip) = p_inf

!-------------------------------------------------------------------------
!-------------------------
!--------------
! Output

 !(1)Right state

     wR_real = wR%f

 !(2)Derivative

    do i = 1, 5
     do j = 1, 5
      dwRdwL(i,j) = wR(i)%df(j)
     end do
    end do

 end subroutine outflow_subsonic_p0

 
!********************************************************************************
!* Slip wall (Compressible)
!********************************************************************************
 subroutine slip_wall(wL_real,nj, dwRdwL,wR_real)

 use derivative_data_df5

 implicit none

 integer , parameter ::  p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter :: two = 2.0_p2

 real(p2)                       , dimension(5)  , intent( in) :: wL_real
 real(p2)                       , dimension(3)  , intent( in) :: nj
 real(p2)                       , dimension(5,5), intent(out) :: dwRdwL
 real(p2)                       , dimension(5)  , intent(out) :: wR_real

 type(derivative_data_type_df5), dimension(5) :: wL, wR
 type(derivative_data_type_df5)               :: dqdn

 real(p2) :: nx, ny, nz
 integer  :: i, j

  nx = nj(1)
  ny = nj(2)
  nz = nj(3)

  wL = wL_real
  call ddt_seed(wL)

! We begin by copying the left state.

  wR = wL

!----------------------------------------------------------------
! Veloicty components
!----------------------------------------------------------------

! Zero normal velocity.

    dqdn = wL(2)*nx + wL(3)*ny + wL(4)*nz

   wR(2) = wL(2) - dqdn * nx
   wR(3) = wL(3) - dqdn * ny
   wR(4) = wL(4) - dqdn * nz


!-------------------------------------------------------------------------
!-------------------------
!--------------
! Output

 !(1)Right state

     wR_real = wR%f

 !(2)Derivative

    do i = 1, 5
     do j = 1, 5
      dwRdwL(i,j) = wR(i)%df(j)
     end do
    end do

 end subroutine slip_wall


!********************************************************************************
!* Outflow supersonic(Compressible)
!*
!* ------------------------------------------------------------------------------
!*  Input: 
!*
!* Output: node(:)%res = the residual computed by the current solution.
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine outflow_supersonic(wL_real,dwRdwL,wR_real)

 use derivative_data_df5

 implicit none

 integer , parameter ::   p2 = selected_real_kind(15) ! Double precision

 real(p2)                      , dimension(5)  , intent( in) :: wL_real
 real(p2)                      , dimension(5,5), intent(out) :: dwRdwL
 real(p2)                      , dimension(5)  , intent(out) :: wR_real

 type(derivative_data_type_df5), dimension(5) :: wL, wR

 integer :: i, j

  wL = wL_real
  call ddt_seed(wL)

!-------------------------
! Right state = left state:

  wR = wL

!-------------------------------------------------------------------------
!-------------------------
!--------------
! Output

 !(1)Right state

     wR_real = wR%f

 !(2)Derivative

    do i = 1, 5
     do j = 1, 5
      dwRdwL(i,j) = wR(i)%df(j)
     end do
    end do

 end subroutine outflow_supersonic

 

 end module bc_states
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

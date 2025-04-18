!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes data and functions relevant to the Euler equations.
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

 module euler_eqns

  implicit none

  integer , parameter :: p2 = selected_real_kind(P=15)

  public :: w2u !primitive to conservative variables
  public :: u2w !conservative to primitive variables
  public

 !--------------------------------------------------------------------
 ! Indices of primitive variables.

  integer             :: ir = 1
  integer             :: iu = 2
  integer             :: iv = 3
  integer             :: iw = 4
  integer             :: ip = 5

 !--------------------------------------------------------------------
 ! Ratio of specific heats for the air.

  real(p2), parameter :: gamma = 1.4_p2

 !--------------------------------------------------------------------
 ! To store free stream values

  real(p2) :: rho_inf, u_inf, v_inf, w_inf, p_inf


 contains


!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  w =    primitive variables (rho,     u,     v,     p)
!* Output:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function w2u(w) result(u)

 use data_module    , only : p2, one, half

 implicit none

 real(p2), dimension(5), intent(in) :: w

!Local variables
 real(p2), dimension(5)             :: u !output

  u(1) = w(1)
  u(2) = w(1)*w(2)
  u(3) = w(1)*w(3)
  u(4) = w(1)*w(4)
  u(5) = w(5)/(gamma-one)+half*w(1)*(w(2)*w(2)+w(3)*w(3)+w(4)*w(4))

 end function w2u
!--------------------------------------------------------------------------------

!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* Output:  w =    primitive variables (rho,     u,     v,     p)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function u2w(u) result(w)

 use data_module    , only : p2, one, half

 implicit none

 real(p2), dimension(5), intent(in) :: u

!Local variables
 real(p2), dimension(5)             :: w !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = u(3)/u(1)
  w(4) = u(4)/u(1)
  w(5) = (gamma-one)*( u(5) - half*w(1)*(w(2)*w(2)+w(3)*w(3)+w(4)*w(4)) )

 end function u2w
!--------------------------------------------------------------------------------


 end module euler_eqns
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

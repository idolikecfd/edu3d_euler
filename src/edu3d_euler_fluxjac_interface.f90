!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes subroutines that comptute numerical fluxes and flux
!* Jacobians by using flux functions written based on automatic differentiation.
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
!********************************************************************************

 module fluxjac_interface

  implicit none

 !This module contains the following subroutines:

  public :: interface_flux   ! compute a numerical flux at interior face
  public :: interface_flux_b ! compute a numerical flux at boundary face

  public :: interface_jac    ! compute a flux Jacobian at interior face
  public :: interface_jac_b  ! compute a flux Jacobian at boundary face

 contains

!********************************************************************************
! Compute a numerical flux at an internal face.
!
! The numerical flux is computed at the edge midpoint between node j and k.
!
!               
!   ------------------> ejk(1:3)*ejk_mag:  This is the edge vector.
!
!   o--------x--------o
!   j       mid       k
!
!
! Left and right states linearly extrpolated from nodes to edge midpoint
! for 2nd-order accuracy. No reconstruction in 1st-order scheme.
!
!              |           . |           |
!              |        .  L |.R         |
!              |      .      |   .       |
!              |    .        |      .    |
!              | .           |        .  |
!              |             |           |
!         ------------o------x------o-------
!                     j     mid     k
!
! Note: the location of the face is always at the midpoint of the edge.
!
!********************************************************************************
 subroutine interface_flux(wj,wk,gradwj,gradwk,njk,ejk,ejk_mag, flux, wsn)

 use derivative_data_df5

 use data_module        , only : p2, ix, iy ,iz
 use input_parameter    , only : inviscid_flux
 use flux_functions_ddt , only :        roe_ddt, &
                                    rusanov_ddt, &
                                        hll_ddt, &
                                       rhll_ddt

 use euler_eqns         , only : w2u, ir, ip
 use input_parameter    , only : accuracy_order, kappa

 implicit none

 real(p2) ::    zero = 0.0_p2
 real(p2) ::     one = 1.0_p2
 real(p2) ::    half = 0.5_p2

!Input
 real(p2), dimension(5)  , intent( in) ::     wj  !primitive variables at j
 real(p2), dimension(5)  , intent( in) ::     wk  !primitive variables at k
 real(p2), dimension(5,3), intent( in) :: gradwj  !gradient at j
 real(p2), dimension(5,3), intent( in) :: gradwk  !gradient at k
 real(p2), dimension(3)  , intent( in) :: njk     !unit drected area vector
 real(p2), dimension(3)  , intent( in) :: ejk     !unit edge vector
 real(p2)                , intent( in) :: ejk_mag !magnitude of the edge vector

!Output
 real(p2), dimension(5)  , intent(out) :: flux !output: numerical flux
 real(p2)                , intent(out) :: wsn  !output: eigenvalue to be used for dtau

!Left and right states
 real(p2)                      , dimension(5  ) :: wL ,wR, dwj, dwk
 type(derivative_data_type_df5), dimension(5  ) :: uL_ddt, uR_ddt

!Dummy array to store Jacobian (computed but not used here).
 real(p2),                       dimension(5,5) :: dummy5x5

!-------------------------------------------------------------------------------
! Left and right solutions at the midpoint:

  if     (accuracy_order == 1) then

   !No reconstruction for 1st-order accuracy
    wL = wj
    wR = wk

  elseif (accuracy_order == 2) then

   !U-MSUCL reconstruction scheme at the edge midpoint for 2nd-order accuracy:
   !See "Resolving Confusion Over Third-Order Accuracy of U-MUSCL", AIAA2021-0056.

   ! U-MUSCL scheme: kappa=1 is central (no dissipation), kappa=0 is linear extrapolation.

     dwj = ( gradwj(:,ix)*ejk(ix) + gradwj(:,iy)*ejk(iy) + gradwj(:,iz)*ejk(iz) )*ejk_mag*half
     dwk = ( gradwk(:,ix)*ejk(ix) + gradwk(:,iy)*ejk(iy) + gradwk(:,iz)*ejk(iz) )*ejk_mag*half

      wL = kappa*half*( wj + wk ) + (one-kappa)*( wj + dwj )
      wR = kappa*half*( wj + wk ) + (one-kappa)*( wk - dwk )

   !If negative density/pressure arises, locally and temporarily revert to first-order.
   !You may want to see how often this is enforced for smooth flows, or how it is affected
   !by the type of solver, the magnitude of the CFL number, initial solution, etc.
   !Is there a better way to control this?

    if ( wL(ir) < zero .or. wL(ip) < zero ) wL = wj
    if ( wR(ir) < zero .or. wR(ip) < zero ) wR = wk

  else

   write(*,*) " Invalid input: accuracy_order = ", accuracy_order
   write(*,*) "  accuracy_order must be 1 or 2. Stop..."
   stop

  endif

!-------------------------------------------------------------------------------
! Covert them to conservative variables in ddt.
! Here, ddt is not really used, but needed since numerical fluxes
! are written based on ddt variables for easy computation of Jacobian.

   uL_ddt = w2u(wL)
   uR_ddt = w2u(wR)

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!  Compute inviscid flux by 3D flux subroutines.
!
!  Note: These flux subroutines are written based on automatic differentiation,
!        and thus return the flux derivative also. But the derivative is not
!        used here, and so it is received in 'dummy5x5', and not used.
!
!  Note: Input variables to each flux function must be derivative-data-type (ddt),
!        which contains derivative information as defined in the module
!        'derivative_data_df5'.
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

 !------------------------------------------------------------
 !  (1) Roe flux
 !------------------------------------------------------------
   if         (trim(inviscid_flux)=="roe") then

     call roe_ddt(uL_ddt,uR_ddt,njk, flux,dummy5x5,wsn)

 !------------------------------------------------------------
 !  (2) Rusanov flux
 !------------------------------------------------------------
   elseif     (trim(inviscid_flux)=="rusanov") then

     call rusanov_ddt(uL_ddt,uR_ddt,njk, flux,dummy5x5,wsn)

 !------------------------------------------------------------
 !  (3) HLL flux
 !------------------------------------------------------------
   elseif     (trim(inviscid_flux)=="hll") then

     call hll_ddt(uL_ddt,uR_ddt,njk, flux,dummy5x5,wsn)

 !------------------------------------------------------------
 !  (4) RHLL flux
 !------------------------------------------------------------
   elseif     (trim(inviscid_flux)=="rhll") then

     call rhll_ddt(uL_ddt,uR_ddt,njk, flux,dummy5x5,wsn,.false.)

 !------------------------------------------------------------
 !  Others...
 !------------------------------------------------------------
   else

    write(*,*) " Invalid input for inviscid_flux = ", trim(inviscid_flux)
    write(*,*) " Choose roe, hll, rusanov, or rhll, and try again."
    write(*,*) " ... Stop."
    stop

   endif

 end subroutine interface_flux
!********************************************************************************


!********************************************************************************
! Compute a numerical flux at a boundary node.
!
! In node-centered edge-based method, a boundary flux is computed at a node.
! This subroutine computes the flux based on the nodal value wj, and a boundary
! state given in wk (conceptually located at the same location as j).
! There is no edge here.
!
! The Roe flux is used to take into account characteristics at boundary.
!
!********************************************************************************
 subroutine interface_flux_b(wj,wk,njk, flux, wsn)

 use derivative_data_df5
 use data_module           , only : p2
 use flux_functions_ddt    , only :     roe_ddt
 use euler_eqns     , only : w2u

 implicit none

!Input
 real(p2), dimension(5)  , intent( in) :: wj  !primitive variables at j
 real(p2), dimension(5)  , intent( in) :: wk  !primitive variables at k
 real(p2), dimension(3)  , intent( in) :: njk !unit boundary element normal vector

!Output
 real(p2), dimension(5)  , intent(out) :: flux !output
 real(p2)                , intent(out) :: wsn  !output

!Left and right states
 real(p2)                      , dimension(5  ) :: wL ,wR
 type(derivative_data_type_df5), dimension(5  ) :: uL_ddt, uR_ddt

!Dummy array to store Jacobian (computed but not used here).
 real(p2),                       dimension(5,5) :: dummy5x5

 !No reconstruction at a boundary node. (edge doesn't exist at a boundary node).

   wL = wj
   wR = wk

 !Covert them to conservative variables in ddt.
 !Here, ddt is not really used, but needed since numerical fluxes
 !are written based on ddt variables for easy computation of Jacobian.

   uL_ddt = w2u(wL)
   uR_ddt = w2u(wR)

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
! Use the Roe flux at boundary to implement a characteristic boundary flux.
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

     call roe_ddt(uL_ddt,uR_ddt,njk, flux,dummy5x5,wsn)

 end subroutine interface_flux_b
!********************************************************************************



!********************************************************************************
! Compute two Jacobians for an internal face:
!
!   dFnduL = d(numerical flux)/d(uj)
!   dFnduR = d(numerical flux)/d(uk)
!
! Flux Jacobians are computed at the edge midpoint between node j and k.
!
!               
!   ------------------> ejk(1:3)*ejk_mag:  This is the edge vector.
!
!   o--------x--------o
!   j       mid       k
!
! No reconstruction for Jacobian, whcih we compute exactly for 1st-order scheme.
!
!********************************************************************************
 subroutine interface_jac(wj,wk,njk, dFnduL, dFnduR )

 use derivative_data_df5
 use data_module           , only : p2
 use input_parameter       , only : inviscid_jac
 use flux_functions_ddt    , only :     roe_ddt, &
                                    rusanov_ddt, &
                                        hll_ddt, &
                                       rhll_ddt

 use euler_eqns            , only : w2u

 implicit none

!Input
 real(p2), dimension(5)  , intent( in) ::     wj,     wk
 real(p2), dimension(3)  , intent( in) :: njk

!Output
 real(p2), dimension(5,5), intent(out) :: dFnduL, dFnduR

!Left and right states

 real(p2), dimension(5)   :: wL ,wR
 real(p2), dimension(5,5) :: dfndu
 real(p2), dimension(5  ) :: dummy5
 real(p2)                 :: wsn

 integer :: i
 type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt

  jac_L_R : do i = 1, 2

  !No reconstruction here.

    wL = wj
    wR = wk

  !Covert them to conservative variables in ddt.
  !Here, ddt is not really used, but needed since numerical fluxes
  !are written based on ddt variables for easy computation of Jacobian.

    uL_ddt = w2u(wL)
    uR_ddt = w2u(wR)

  !This step determines with respect to which state the flux derivative is computed.

    if (i==1) then

     !This sets the derivative of uL_ddt = 1.
      call ddt_seed(uL_ddt)

    else

     !This sets the derivative of uR_ddt = 1.
      call ddt_seed(uR_ddt)

    endif

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
!  Compute inviscid Jacobian by 3D flux subroutines.
!
!  Note: These flux subroutines are written based on automatic differentiation,
!        and thus return the flux and its derivative. Here, we only want the
!        derivative. So, the flux is received in 'dummy5', and not used.
!
!  Note: Input variables to each flux function must be derivative-data-type (ddt),
!        which contains derivative information as defined in the module
!        'derivative_data_df5'.
!
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

 !------------------------------------------------------------
 !  (1) Roe flux
 !------------------------------------------------------------
   if         (trim(inviscid_jac)=="roe") then

     call roe_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)

 !------------------------------------------------------------
 !  (2) Rusanov flux
 !------------------------------------------------------------
   elseif     (trim(inviscid_jac)=="rusanov") then

     call rusanov_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)

 !------------------------------------------------------------
 !  (3) HLL flux
 !------------------------------------------------------------
   elseif     (trim(inviscid_jac)=="hll") then

     call hll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)

 !------------------------------------------------------------
 !  (4) RHLL flux: the last argumant -> exact_jac = .false.
 !                  so that the jac = a1*HLL_jac+a2*Roe_jac
!                   with a1 and a2 not differentiated.
 !------------------------------------------------------------
   elseif     (trim(inviscid_jac)=="rhll") then

     call rhll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn,.false.)

 !------------------------------------------------------------
 !  Others...
 !------------------------------------------------------------
   else

    write(*,*) " Invalid input for inviscid_jac = ", trim(inviscid_jac)
    write(*,*) " Choose roe or rhll, and try again."
    write(*,*) " ... Stop."
    stop

   endif

   if (i==1) then
    dFnduL = dfndu
   else
    dFnduR = dfndu
   endif

  end do jac_L_R

 end subroutine interface_jac
!********************************************************************************


!********************************************************************************
! Compute two Jacobians for the flux defined at at a boundary node:
!
!   dFnduL = d(numerical flux)/d(uj)
!   dFnduR = d(numerical flux)/d(uk)
!
! In node-centered edge-based method, a boundary flux is computed at a node.
! This subroutine computes the flux based on the nodal value wj, and a boundary
! state given in wk (conceptually located at the same location as j).
! There is no edge here.
!
! The Roe flux is used to take into account characteristics at boundary.
!
!********************************************************************************
 subroutine interface_jac_b(wj,wk,njk, dFnduL, dFnduR )

 use derivative_data_df5
 use data_module           , only : p2
 use flux_functions_ddt    , only : roe_ddt
 use euler_eqns            , only : w2u

 implicit none

!Input
 real(p2), dimension(5)  , intent( in) :: wj  !primitive variables at j
 real(p2), dimension(5)  , intent( in) :: wk  !primitive variables at k
 real(p2), dimension(3)  , intent( in) :: njk !unit boundary element normal vector

!Output
 real(p2), dimension(5,5), intent(out) :: dFnduL, dFnduR !Left and right Jacobians

!Left and right states

 real(p2), dimension(5)   :: wL ,wR
 real(p2), dimension(5,5) :: dfndu
 real(p2), dimension(5  ) :: dummy5
 real(p2)                 :: wsn

 integer :: i
 type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt

  jac_L_R : do i = 1, 2

  !No reconstruction here (edge doesn't exist at a boundary node).

    wL = wj
    wR = wk

  !Covert them to conservative variables in ddt.
  !Here, ddt is not really used, but needed since numerical fluxes
  !are written based on ddt variables for easy computation of Jacobian.

    uL_ddt = w2u(wL)
    uR_ddt = w2u(wR)

    if (i==1) then
     !This sets the derivative of uL_ddt = 1.
     call ddt_seed(uL_ddt)
    else
     !This sets the derivative of uR_ddt = 1.
     call ddt_seed(uR_ddt)
    endif

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------
! Use the Roe flux at boundary to implement a characteristic boundary flux.
!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

     call roe_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndu,wsn)


   if (i==1) then
    dFnduL = dfndu
   else
    dFnduR = dfndu
   endif

  end do jac_L_R

 end subroutine interface_jac_b
!********************************************************************************


 end module fluxjac_interface


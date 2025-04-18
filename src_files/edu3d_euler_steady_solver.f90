!********************************************************************************
! Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
!
! This module containes subroutines that solve the global system of nonlinear
! residual equations:
!
!     Res(U) = 0.
!
! - This is a nonlinear system solver, which some people call a steady solver.
!
! - So, this can solve unsteady residuals as well, which incorporate physical
!   time derivatives as source terms, or any nonlinear system as long as
!   a Jacobian is provided.
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

 module steady_solver

 !This module contains the following subroutines:

  public :: solve_steady                 ! solve residual equations (nonlinear system)

  public :: compute_residual_norm        ! compute residual norms
  public :: compute_local_time_step_dtau ! compute pseudo time step

 !Nonlinear solvers:

  public :: explicit_pseudo_time_forward_euler  !Explicit forward Euler
  public :: explicit_pseudo_time_tvdrk          !Explicit 2-stage TVD RK
  public :: implicit                            !Implicit defect-correction
  public :: jfnk                                !Jacobian-Free Newton-Krylov

 contains

!*******************************************************************************
!
! This subroutine solves the global system of nonlinear residual eqs: Res(U) = 0
!
! Four nonlinear solvers are implemented:
!
! (1)Forward Euler  explicit time stepping. ! many iterations
! (2)2-stage TVD RK explicit time stepping. ! many iterations
! (3)Implicit defect-correction solver.     ! much less iterations
! (4)Jacobian-Free Newton-Krylov with GCR.  ! can be much much less iterations
!
! Note: In (4), the solver (3) is used as a preconditioner.
!
!*******************************************************************************
 subroutine solve_steady

 use input_parameter, only : solver_type, solver_max_itr, solver_tolerance
 use input_parameter, only : inviscid_flux, inviscid_jac, CFL, accuracy_order
 use data_module    , only : i_iteration
 use data_module    , only : lrelax_sweeps_actual, lrelax_roc
 use residual       , only : compute_residual
 use data_module    , only : gcr_proj_actual, gcr_rkr0, n_residual_evaluation

 integer ,      parameter :: p2 = selected_real_kind(15) ! Double precision
 real(p2), dimension(5,3) :: res_norm
 real(p2), dimension(5,3) :: res_norm_initial
 integer                  :: L1 = 1

  write(*,*)
  write(*,*) "---------------------------------------"
  write(*,*) " Begin nonlinear iteration"
  write(*,*)

  write(*,*)
  write(*,'(2a)')   "      solver_type =  ", trim(solver_type)
  write(*,'(a,i1)') "   accuracy_order =  ", accuracy_order
  write(*,'(2a)')   "    inviscid_flux =  ", trim(inviscid_flux)
  write(*,'(2a)')   "    inviscid_jac  =  ", trim(inviscid_jac)
  write(*,'(a,es9.2)') "             CFL  = ", CFL
  write(*,*)

! Warning...

  if (trim(inviscid_flux) == "rhll") then
    write(*,*)
    write(*,*) "     Note: RHLL may not converge to machine zero unless n1 and n2 are frozen..."
    write(*,*)
  endif

! Initialization

           i_iteration = 0
  lrelax_sweeps_actual = 0
            lrelax_roc = 0.0_p2
 n_residual_evaluation = 0

! Print the title.

   if ( trim(solver_type) == "jfnk" ) then

    write(*,'(a92)') "|     GCR solver"
    write(*,'(2a)') " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy", &
                   "    |   proj     reduction"

   elseif ( trim(solver_type) == "implicit" ) then

    write(*,'(a90)') "|     GS relax"
    write(*,'(2a)') " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy", &
                   "    |  relax    reduction"

   else

    write(*,*)
    write(*,'(2a)') " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy"

   endif


! Begin nonlinear iteration.

  iteration : do

    call compute_residual

    call compute_residual_norm(res_norm)

  ! Print the current residual norms on screen.

    if ( trim(solver_type) == "implicit" ) then

     write(*,'(i10,5es13.3,a,i6,es12.1)') i_iteration, res_norm(:,L1),  " | ",   &
                                              lrelax_sweeps_actual, lrelax_roc

    elseif ( trim(solver_type) == "jfnk" ) then

     write(*,'(i10,5es13.3,a,i6,es12.1)') i_iteration, res_norm(:,L1),  " | ",   &
                                              gcr_proj_actual, gcr_rkr0
    else

     write(*,'(i10,5es13.3)') i_iteration, res_norm(:,L1)

    endif


  ! Check convergence and exit if converged.

    if (i_iteration == 0) then !<- No iteration has been performed yet.

    ! Save off the initial residual norm.
      res_norm_initial = res_norm

    else

    ! Exit if the residual is reduced by a user-specified order.
     if (maxval(res_norm(:,1)/(res_norm_initial(:,1)+1.0e-07_p2) ) < solver_tolerance) then
      write(*,*)
      write(*,*) " Converged. Res reduction achieved = ", maxval(res_norm(:,1)/(res_norm_initial(:,1)+1.0e-07_p2))
      write(*,*)
      exit iteration
     endif

    ! Update the reference residual norm if increased.
      res_norm_initial(:,1) = max( res_norm_initial(:,1), res_norm(:,1) )

    endif

  ! Exit if the iteration reaches a user-specified maximum number.
    if (i_iteration == solver_max_itr) exit iteration


  ! Proceed to the next iteration.

    call compute_local_time_step_dtau

  !----------------------------------------------------------
  ! Perform one iteration: u <- u + du

     method : if ( trim(solver_type) == "fwd_euler") then

      ! Forward Euler explicit iteration
        call explicit_pseudo_time_forward_euler

     elseif ( trim(solver_type) == "tvdrk" ) then

      ! TVD RK explicit iteration
        call explicit_pseudo_time_tvdrk
 
     elseif ( trim(solver_type) == "implicit" ) then

      ! Implicit Defect-Correction solver
        call implicit

     elseif ( trim(solver_type) == "jfnk" ) then

      ! Jacobian-Free Newton-Krylov with GCR
        call jfnk

     endif method
  !----------------------------------------------------------

   i_iteration = i_iteration + 1

  end do iteration

  write(*,*)
  write(*,*) " Total number of residual evaluations = ", n_residual_evaluation
  write(*,*)

 end subroutine solve_steady


!********************************************************************************
!
! This subroutine computes the residual norms: L1, L2, L_infty
!
!********************************************************************************
 subroutine compute_residual_norm(res_norm)

 use data_module, only : nnodes, res

 implicit none

 integer , parameter ::   p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter :: zero = 0.0_p2
 real(p2), parameter ::  one = 1.0_p2

 real(p2), dimension(5,3), intent(out) :: res_norm

!Local variables
 real(p2), dimension(5) :: residual
 integer                :: i

   res_norm(:,1) =  zero
   res_norm(:,2) =  zero
   res_norm(:,3) = - one


  loop_nodes : do i = 1, nnodes

   residual = abs( res(i,:) )

   res_norm(:,1) = res_norm(:,1)    + residual    !L1   norm
   res_norm(:,2) = res_norm(:,2)    + residual**2 !L2   norm
   res_norm(:,3) = max(res_norm(:,3), residual)   !Linf norm

  end do loop_nodes


   res_norm(:,1) =      res_norm(:,1)/real(nnodes,p2)  !L1 norm
   res_norm(:,2) = sqrt(res_norm(:,2)/real(nnodes,p2)) !L2 norm

 end subroutine compute_residual_norm
!********************************************************************************

!********************************************************************************
!
! This subroutine computes the local pseudo-time step at each node.
!
!********************************************************************************
 subroutine compute_local_time_step_dtau

 use data_module    , only : nnodes, dtau, wsn, vol
 use input_parameter, only : CFL

 implicit none

 integer , parameter ::   p2 = selected_real_kind(P=15)
 real(p2), parameter :: half = 0.5_p2
 integer             :: i


  loop_nodes : do i = 1, nnodes

  ! Local time step: dtau = CFL*volume/sum(0.5*max_wave_speed*face_area).

   dtau(i) = CFL*vol(i) / ( half*wsn(i) )

  end do loop_nodes


 end subroutine compute_local_time_step_dtau
!********************************************************************************


!********************************************************************************
!
! Perform 1 nonlinear iteration by the forward-Euler pseudo-time stepping.
!
!   u^{n+1} = u^n - (dt/vol)*Res(u^n)
!
! Note: This is very slow to converge to steady state.
!
!*******************************************************************************
 subroutine explicit_pseudo_time_forward_euler

 use data_module, only : nnodes, u, w, du, dtau, res, vol
 use euler_eqns , only : u2w

 implicit none

 integer :: i


  loop_nodes : do i = 1, nnodes

  !This is the correction suggested by the forward Euler scheme:

   du(i,:) = - dtau(i)/vol(i) * res(i,:)

  !Update the solution by applying the correction.

    u(i,:) = u(i,:) + du(i,:)

  !Update the primitive variables.

    w(i,:) = u2w( u(i,:) )

  end do loop_nodes


 end subroutine explicit_pseudo_time_forward_euler
!********************************************************************************

!********************************************************************************
!
! Perform 1 nonlinear iteration by an explicit 2-stage RK pseudo-time stepping.
!
!  Two-stage Runge-Kutta scheme: u^n is saved as u0(:,:)
!   1. u^*     = u^n - (dt/vol)*Res(u^n)
!   2. u^{n+1} = 1/2*u^n + 1/2*[u^* - (dt/vol)*Res(u^*)]
!
! Note: This is very slow to converge to steady state.
!
!*******************************************************************************
 subroutine explicit_pseudo_time_tvdrk

 use data_module, only : nnodes, u0, u, w, du, dtau, res, vol
 use euler_eqns , only : u2w
 use residual   , only : compute_residual

 implicit none

 integer , parameter :: p2 = selected_real_kind(P=15)
 real(p2), parameter :: half = 0.5_p2
 integer             :: i

 !-----------------------------
 !- 1st stage of Runge-Kutta:
 !-----------------------------

   u0 = u

  loop_nodes : do i = 1, nnodes

  !This is the correction suggested by the forward Euler scheme:

   du(i,:) = - dtau(i)/vol(i) * res(i,:)

  !Update the solution by applying the correction.

    u(i,:) = u(i,:) + du(i,:)

  !Update the primitive variables.

    w(i,:) = u2w( u(i,:) )

  end do loop_nodes

 !-----------------------------
 !- 2nd stage of Runge-Kutta:
 !-----------------------------

  call compute_residual

  u = half*( u + u0 )

  loop_nodes2 : do i = 1, nnodes

  !This is the correction suggested by the forward Euler scheme:

   du(i,:) = - half*dtau(i)/vol(i) * res(i,:)

  !Update the solution by applying the correction.

    u(i,:) = u(i,:) + du(i,:)

  !Update the primitive variables.

    w(i,:) = u2w( u(i,:) )

  end do loop_nodes2

 end subroutine explicit_pseudo_time_tvdrk
!********************************************************************************

!********************************************************************************
!
! Perform 1 nonlinear iteration by the implicit defect-correction solver.
!
! (1)Compute the residual Jacobian:
!
!    - Exact derivative of 1st-order residual (no LSQ gradients).
!    - This is Newton's method for 1st-order scheme (with a large CFL).
!
! (2)Solve the linear system: (V/dtau+dRes1/du)*du = -Res, where dRes1/du is
!    the exact derivative of the 1st-order residual (no LSQ gradients).
!
!    - Gauss-Seidel relaxation is implemented.
!    - Max relaxations and tolerance are specified by input parameters.
!
! (3)Update the solution: u <- u + du.
!
!********************************************************************************
 subroutine implicit

 use data_module  , only : nnodes, u, w, du
 use euler_eqns   , only : u2w
 use jacobian     , only : compute_jacobian
 use linear_solver, only : linear_relaxation
 use data_module  , only : p2

 implicit none

 integer  :: i
 real(p2) :: omegan !under-relaxation factor for nonlinear iteration

  !Compute the residual Jacobian

   call compute_jacobian

  !Compute du by relaxing (approximately solving) the linear system.

   call linear_relaxation

  !Update the solution by applying the correction with a safety factor.

   loop_nodes : do i = 1, nnodes

    !Check negative desnity/pressure and compute a safety factor if needed.
     omegan = safety_factor( u(i,:), du(i,:) )

    !Update the conservative variables
     u(i,:) = u(i,:) + omegan*du(i,:)

    !Update the primitive variables.
     w(i,:) = u2w( u(i,:) )

   end do loop_nodes


 end subroutine implicit
!********************************************************************************


!********************************************************************************
!
! Perform 1 nonlinear iteration by the implicit defect-correction solver.
!
! (1)Compute the residual Jacobian:
!
!    - Exact derivative of 1st-order residual.
!    - This is used as a preconditioner matrix, Jac.
!
! (2)Solve the exact linearized system: (V/dtau+dRes/du)*du = -Res:
!
!    - Frechet derivative is used to evaluate dRes/dU*du.
!    - Gauss-Seidel relaxation is implemented to approximately invert Jac.
!    - Max relaxations and tolerance are specified by input parameters.
!    - Max GCR projection and tolerance are specified by input parameters.
!
! (3)Update the solution: u <- u + du.
!
! Note: This is Newton's method for both 1st- and 2nd-order schemes if
!       CFL is sufficiently large.
!
! Note: Jacobian-Free means that the exact Jacobian dRes/du is not stored.
!       But it requires a preconditioner matrix, which is given by
!       the exact derivative of the 1st-order scheme.
!
!********************************************************************************
 subroutine jfnk

 use input_parameter, only : gcr_projection !Max GCR projection
 use input_parameter, only : gcr_mu         !GCR tolerance (reduction)

 use data_module    , only : gcr_proj_actual      !Actual # of projections
 use data_module    , only : lrelax_sweeps_actual !Last GS relaxations
 use data_module    , only : gcr_rkr0             !Actual GCR-res reduction
 use data_module    , only : status_pc            !0=success in preconditioner
 use data_module    , only : status_gcr           !0=success in GCR

 use data_module    , only : nnodes, u, w, du, i_iteration
 use euler_eqns     , only : u2w
 use jacobian       , only : compute_jacobian
 use gcr            , only : gcr_solver

 use data_module    , only : p2

 implicit none

 integer  :: i
 real(p2) :: omegan !under-relaxation factor for nonlinear iteration

  !Compute the residual Jacobian used as a preconditioner matrix.
  !This matrix is unchanged during the GCR solve.

   call compute_jacobian

  !Compute du by JFNK-GCR solver by solving the exactly linearized system.

   call gcr_solver(gcr_projection, gcr_mu, gcr_proj_actual, lrelax_sweeps_actual, &
                                      gcr_rkr0, status_pc, status_gcr,i_iteration+1 )

  !Update the solution by applying the correction with a safety factor.

   loop_nodes : do i = 1, nnodes

    !Check negative desnity/pressure and compute a safety factor if needed.
     omegan = safety_factor( u(i,:), du(i,:) )

    !Update the conservative variables
     u(i,:) = u(i,:) + omegan*du(i,:)

    !Update the primitive variables.
     w(i,:) = u2w( u(i,:) )

   end do loop_nodes


 end subroutine jfnk
!********************************************************************************


!********************************************************************************
! Compute a safety factor (under relaxation for nonlinear update) to make sure
! the updated density and pressure are postive.
!
! This is just a simple exmaple.
! Can you come up with a better and more efficient way to control this?
!
!********************************************************************************
 function safety_factor(u,du)

 use data_module    , only : p2
 use euler_eqns     , only : u2w, ir, ip
 use input_parameter, only : M_inf

 implicit none

 real(p2) ::    zero = 0.00_p2

 real(p2), dimension(5), intent(in) :: u, du
 real(p2)                           :: safety_factor

 real(p2), dimension(5)             :: u_updated
 real(p2), dimension(5)             :: w_updated
 real(p2)                           :: p_updated

 ! Default safety_factor

   safety_factor = 1.0_p2

 ! Temporarily update the solution:

       u_updated = u + safety_factor*du
       w_updated = u2w(u_updated)

 !-----------------------------
 ! Return if both updated density and pressure are positive

   if ( w_updated(ir) > zero .and. w_updated(ip) > zero ) then
 
   !Good. Keep safety_factor = 1.0_p2, and return.

    return

   endif

 !-----------------------------
 ! Negative density fix

   if ( w_updated(ir) <= zero) then ! meaning du(ir) < zero, reducing the density.

     safety_factor = -u(ir)/du(ir) * 0.25_p2 ! to reduce the density only by half.

   endif

 !-----------------------------
 ! Negative pressure fix
 !
 ! Note: Take into account a case of density < 0 and pressure > 0.
 !       Then, must check if the safety factor computed above for density
 !       will give positive pressure also, and re-compute it if necessary.

   if ( w_updated(ip) <= zero .or. w_updated(ir) <= zero) then

   !Note: Limiting value of safety_factor is zero, i.e., no update and pressure > 0.
    do

          u_updated = u + safety_factor*du
          w_updated = u2w(u_updated)
          p_updated = w_updated(ip)

     ! For low-Mach flows, theoretically, pressure = O(Mach^2).
     ! We require the pressure be larger than 1.0e-05*(free stream Mach)^2.
     ! Is there a better estimate for the minimum pressure?

      if (p_updated > 1.0e-05_p2*M_inf**2) exit

      safety_factor = 0.75_p2*safety_factor ! reduce the factor by 25%, and continue.

    end do

   endif


 end function safety_factor


 end module steady_solver
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


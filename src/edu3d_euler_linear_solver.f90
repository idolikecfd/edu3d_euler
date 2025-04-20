!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!*
!* This module containes a linear relaxtion solver. Currently, the sequential
!* Gauss-Seidel scheme is implemented. Other schemes may be implemented and
!* explored: e.g., multi-color GS, symmetric GS, Jacobi, etc.
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

 module linear_solver

  implicit none

 !This module contains the following subroutine:

  public :: linear_relaxation ! perform linear relaxations.

 contains

!********************************************************************************
!
! This subroutine relaxes the linear system by Sequential Gauss-Seidel method.
!
! - Other methods are possible: multi-color GS, symmetric GS, etc.
! - This subroutine is used in the preconditioner for the JFNK solver.
! - Convergence results will be printed in the file: fort.1000.
!
!********************************************************************************
 subroutine linear_relaxation

 use data_module    , only : p2, zero, one, my_eps
 use data_module    , only : nnodes, k0, nnghbrs, nghbr
 use data_module    , only : du, res, jac_off, jac_dinv
 use data_module    , only : i_iteration
 use data_module    , only : lrelax_sweeps_actual, lrelax_roc
 use data_module    , only : gcr_proj_actual

 use input_parameter, only : solver_type
 use input_parameter, only : lrelax_sweeps, lrelax_tolerance, lrelax_scheme

 implicit none

!Local variables
 real(p2), dimension(5)   :: b
 integer                  :: i, k

 integer  :: ii

 !Linear system residual
 real(p2), dimension(5  ) :: linear_res

 !Residual norms(L1,L2,Linf) for the linear system
 real(p2), dimension(5,3) :: linear_res_norm
 real(p2), dimension(5,3) :: linear_res_norm_initial
 real(p2), dimension(5,3) :: linear_res_norm_previous

 !Rate of convergence. Residual reduction from the initial.
 real(p2) :: roc
 real(p2) :: roc_previous

 !Under-relaxation parameter.
 real(p2) :: omega_lrelax

      lrelax_sweeps_actual = 0      ! to count actual number of relaxations
              omega_lrelax = 1.0_p2 ! possible under-relaxation parameter


 if    ( trim(solver_type) == "implicit" ) then

 write(1000,*) "-------------------------------------------------------------------------------"
 write(1000,*) " Nonlinear Iteration = ", i_iteration
 write(1000,*) "-------------------------------------------------------------------------------"

 elseif( trim(solver_type) == "jfnk" ) then

 write(1000,*) "-------------------------------------------------------------------------------"
 write(1000,*) " Nonlinear Iteration = ", i_iteration, " GCR projection = ", gcr_proj_actual
 write(1000,*) "-------------------------------------------------------------------------------"

 endif


 write(1000,*)
 write(1000,*) " lrelax_scheme = ", trim(lrelax_scheme)
 write(1000,*)

!--------------------------------------------------------------------------------
! 1. Initialize the correction

     du = zero

!--------------------------------------------------------------------------------
! 2. Linear Relaxation (Sweep)

  relax : do ii = 1, lrelax_sweeps

   linear_res_norm(:,1) = zero

 !---------------------------------------------------------
 ! Sequential Gauss-Seidel Relaxation(sweep)

   if (trim(lrelax_scheme) == "gs") then

    gs_node_loop : do i = 1, nnodes

    ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]!

        b = -res(i,:)
       gs_node_nghbrs : do k = 1, nnghbrs(i)
        b = b - matmul( jac_off( k0(i)+k, :,:), du(nghbr(k0(i)+k),:) )
       end do gs_node_nghbrs

    ! Update du by the GS relaxation:
    !
    ! e.g., for 3 nghbrs, perform the relaxation in the form:
    !
    !                     diagonal block        sum of off-diagonal block contributions
    !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }

        linear_res = matmul( jac_dinv(i,:,:), b ) - du(i,:)

           du(i,:) = du(i,:) + omega_lrelax * linear_res

    ! To compute the linear residual norm
       linear_res_norm(:,1) = linear_res_norm(:,1) + abs( linear_res )

    end do gs_node_loop

   else

    write(*,*) " Sorry, only 'gs' is available at the moment..."
    write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop."
    stop 127

   endif
 !---------------------------------------------------------

    linear_res_norm(:,1) = linear_res_norm(:,1) / real(nnodes, p2)

!--------------------------------------------------------------------------------
! 3. Check the linear residual.

 !  After the first relaxation

    if (ii==1) then

     !-----------------------------------------------------------------
     ! Print the initial linear residual norm in the file 'fort.1000'.

     linear_res_norm_initial = linear_res_norm
     write(1000,'(a,i10,a,es12.5)') "  after relax ", ii, &
                " max(L1 norm) = ", maxval(linear_res_norm(:,1))

     if ( maxval(linear_res_norm(:,1)) < my_eps ) then
      write(1000,*) " Machine zero res reached. Exit GS relaxation. Total sweeps = ", ii
      write(1000,*) " tolerance_linear = ", lrelax_tolerance
      lrelax_sweeps_actual = ii
      exit relax
     endif

 !  After the second relaxation

    else

      roc = maxval(linear_res_norm(1:5,1)/linear_res_norm_initial(1:5,1))

     !-----------------------------------------------------------------
     ! Print the current linear residual norm in the file 'fort.1000'.

       if (roc < one) then

        write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,f8.3)')   "  after relax ", ii,          &
                  " max(L1 norm) = ", maxval(linear_res_norm(1:5,1)), omega_lrelax, roc
       else

        write(1000,'(a,i10,a,es12.5,2x,f6.3,2x,f8.3,a)') "  after relax ", ii,          &
                  " max(L1 norm) = ", maxval(linear_res_norm(1:5,1)), omega_lrelax, roc," <- diverge"
       endif


     !-----------------------------------------------------------------
     ! Tolerance met: Exit

      if (roc < lrelax_tolerance) then
       write(1000,*)
       write(1000,*) " Tolerance met. Exit GS relaxation. Total sweeps = ", ii
       write(1000,*) " tolerance_linear = ", lrelax_tolerance
       lrelax_sweeps_actual = ii
       exit relax

     !-----------------------------------------------------------------
     ! If tolerance is NOT met.

      else

       !---------------------------------------------
       ! Stop if the linear residual is too small.

        if ( maxval(linear_res_norm(1:5,1)) < my_eps ) then
         write(1000,*)
         write(1000,*) " Residuals too small. Exit GS relaxation. Total sweeps = ", ii
         write(1000,*) " maxval(linear_res_norm(1:5,1)) = ", maxval(linear_res_norm(1:5,1))
         lrelax_sweeps_actual = ii
         exit relax
        endif

       !-------------------------------------------------
       ! Stop if the maximum number of sweeps is reached.

        if (ii == lrelax_sweeps) then
         write(1000,*)
         write(1000,*) " Tolerance not met... sweeps = ", lrelax_sweeps
         lrelax_sweeps_actual = lrelax_sweeps
        endif

      endif

     ! End of Tolerance met or NOT met.
     !-----------------------------------------------------------------

     !-----------------------------------------------------------------
     ! Do something if we continue.
     ! Reduce/increase the relaxation factor if it's going worse/better.

       roc_previous = maxval(linear_res_norm_previous(1:5,1)/linear_res_norm_initial(1:5,1))

      !--------------------------------------------------------------
      ! If the linear residual goes beyond the initial one,
       if (roc > one) then

        ! REDUCE if the linear residual increases from the previous relaxation.

        if (roc > roc_previous) omega_lrelax = max(0.95_p2*omega_lrelax, 0.05_p2)

      !--------------------------------------------------------------
      ! If the linear residual is smaller than the initial one,
       else

        ! INCREASE if doing fine.
        omega_lrelax = min(1.05_p2*omega_lrelax, one)
        if (roc < roc_previous) omega_lrelax = min(1.25_p2*omega_lrelax, one)

       endif

     ! End of Do something if we continue.
     !-----------------------------------------------------------------

    endif


    linear_res_norm_previous = linear_res_norm

  end do relax

! End of 3. Linear Relaxation (Sweep)
!--------------------------------------------------------------------------------

    lrelax_roc = roc
    write(1000,*)

    return

 end subroutine linear_relaxation
!********************************************************************************



 end module linear_solver
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

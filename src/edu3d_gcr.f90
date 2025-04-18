!********************************************************************************
! Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
! This module containes subroutines for performing a Jacobian-Free Newton-Krylov
! solver by using a general GCR solver module.
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

 module gcr

 ! This module contains the following subroutine:

 ! General GCR sovler (version 3) and associated subroutines:
  public :: gcr_solver
  public :: gcr_dimension
  public :: gcr_rhs_vector
  public :: gcr_preconditioner
  public :: gcr_compute_Ap
  public :: gcr_store_x_in_your_array
  public :: gcr_dot_product
  public :: gcr_rms
  public :: gcr_length
  public :: rms
  public :: length

 ! These are subroutines from the ED3-Euler program,
 ! plugged into the GCR solver for an easy installation:
 ! [If you want to use GCR, just replace these subroutines by yours.]
  public :: your_dimension
  public :: your_rhs_vector
  public :: your_preconditioner
  public :: your_compute_Ap
  public :: vector_dt_term
  public :: your_store_x_in_your_array

! Constants used in GCR
  integer , parameter ::   p2 = selected_real_kind(p=15)
  real(p2), parameter :: zero = 0.0_p2
  real(p2), parameter ::  one = 1.0_p2

 contains

!********************************************************************************
!
!      Generalized Conjugate Residual (GCR) with a variable preconditioner.
!
! This subroutine 'gcr_solver' solves the linear system: Ax = b, up to a specified
! tolerance by the GCR method with a variable preconditioner that approximately
! computes M^{-1} (e.g., by a relaxation scheme), where M is a preconditioning matrix.
!
! References:
!
!  [1] K. Abe and S.-L. Zhang,
!     A Variable Preconditioning Using The SOR Method for GCR-Like Methods,
!     International Journal of Numerical Analysis and Modeling , Vol.2, No.2,
!     pp.147-161, 2005.
!
!  [2] Mohagna J. Pandya, Boris Diskin, James L. Thomas, Neal T. Frink,
!      Improved Convergence and Robustness of USM3D Solutions on Mixed-Element Grids,
!      AIAA Journal, 2016, Vol.54: 2589-2610, 10.2514/1.J054545
!
!--------------------------------------------------------------------------------
!
!  Input: max_projection_gcr = Maximum GCR projections (search directions)
!                     mu_gcr = Tolerance for the linear residual reduction.
!
! Output: actual_projections = Actual projections performed.
!              actual_sweeps = Total # of relaxations (sweeps) to compute M^{-1}.
!                   rkr0_gcr = Final linear residual reduction
!                  status_pc = success(0), diverged(1 or 2)
!                 status_gcr = success(0), not converged (2)
!                          x = Solution copied to your code before exit.
!
!  Note: Success means that the tolerance is met within max_projection_gcr.
!
!--------------------------------------------------------------------------------
!
! To use 'gcr_solver' in your code:
!
!  (1)Write the following subroutines using your data in your code:
!
!      your_dimension()            : Return total # of discrete unknowns
!      your_rhs_vector()           : Compute b of Ax=b.
!      your_preconditioner()       : Compute p=M^{-1}*r for a given vector r.
!      your_compute_Ap()           : Compute Ap for a given vector p.
!      your_dot_product()          : Compute dot product of two vectors.
!      your_store_x_in_your_array(): Copy x to your own soluition array.
!
!     and add them at the end of this module or make them accesible via use
!     statement where they are called.
!
!     These subroutines are the ones that use specific data structure and
!     subroutines/functions in your code. Only you can write them!
!
!  (2)Then, call gcr_solver() inside your code, e.g.,
!
!       call gcr_solver(10,0.5,.....)
!
!--------------------------------------------------------------------------------
!
!
!  Note: To solve a linear system well, e.g., mu_gcr=1.0e-06, but then you'll
!        need to set a large value for max_projection_gcr such as 50-100, or more.
!
!  Note: To solve a linear system in Newton's method for nonlinear equations,
!        it's perhaps enouigh to set mu_gcr=0.9 with max_projection_gcr = 4, for
!        example.
!
!  Note: Preconditioner is assumed to be an approximate matrix M approximately
!        inverted by a relaxation scheme (e.g., Gauss Seidel). 'actual_sweeps'
!        is the actual nummber of relaxations performed; it is accumulated over
!        the entire GCR projections. The number of relaxations can change every
!        time the preconditioning is performed; thus a variable preconditioner.
!
!  Note: Preconditioner does not have to be variable, and can be fixed.
!        If the preconditioner is of iLU-type, we'll always have status_pc = 0
!        and actual_sweeps=0.
!
!  Note: In CFD codes, a variable preconditioner is typically available as
!        an implicit solver: approximately solve M*dq = -Res(q), where M is an
!        approximate Jacobian. If you already have it, you can use your implicit
!        solver (as it is) as the preconditioner.
!
!  Note: Parallel version may be wrriten by you. Think about which part
!        of the code needs to be parallelized.
!
!--------
! Notes:
!
!  The purpose of this code is to give a beginner an opportunity to learn.
!  Hence, the focus here is on the simplicity. The code is not optimized for
!  efficiency; and no consideration for parallelization is given.
!
!  If the code is not simple enough to understand, please send questions to Hiro
!  at sunmasen(at)hotmail.com. He'll greatly appreciate it and revise the code.
!
!  If the code helps you understand how to write your own code that is more
!  efficient and has more features (possibly in a different programming language),
!  it'll have served its purpose.
!
!
!        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!
! the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!
! This is Version 3 (October 2017).
! This F90 code is written and made available for an educational purpose.
! This file may be updated in future
!
!********************************************************************************
 subroutine gcr_solver(max_projection_gcr, mu_gcr, actual_projections, &
                        actual_sweeps, rkr0_gcr, status_pc, status_gcr,i_iteration )

 implicit none

 integer , intent( in) :: max_projection_gcr !Maximum search directions
 real(p2), intent( in) :: mu_gcr             !Tolerance for GCR
 integer , intent( in) :: i_iteration        !Nonlinear iteration counter

 integer , intent(out) :: actual_projections !Actual GCR projections performed
 integer , intent(out) :: actual_sweeps      !Actual preconditioning relaxations
 real(p2), intent(out) :: rkr0_gcr           !Actual residual reduction
 integer , intent(out) :: status_pc          !Preconditioner status
 integer , intent(out) :: status_gcr         !GCR status

 real(p2), dimension(:)  , allocatable :: x    !Solution vector
 real(p2), dimension(:)  , allocatable :: r    !Residual vector
 real(p2), dimension(:)  , allocatable :: r0   !Initial residual vector
 real(p2), dimension(:)  , allocatable :: a, b !Temporary arrays
 real(p2), dimension(:,:), allocatable :: p    !Search directions
 real(p2), dimension(:,:), allocatable :: Ap   !Ap vectors
 real(p2), dimension(:)  , allocatable :: Apip1!Temporary

 real(p2) :: my_eps, alpha, beta, rms_initial, length_p, length_r
 integer  :: i, k, istat_preconditioner, istat_preconditioner_0
 integer  :: ndim, sweeps
 logical  :: stall_condition

             my_eps  = epsilon(one)
          status_pc  = 0
          status_gcr = 0
              alpha  = one
       actual_sweeps = 0
  actual_projections = 0
     stall_condition = .false.

 write(2000,*) "-------------------------------------------------------------------------------"
 write(2000,*) " Nonlinear Iteration = ", i_iteration
 write(2000,*) "-------------------------------------------------------------------------------"

!------------------------------------------------------------------
! Give me the problem dimension, ndim = total number of discrete unknowns.

  call gcr_dimension(ndim)

!------------------------------------------------------------------
! Allocate gcr arrays:

  allocate(  x(ndim) )
  allocate(  r(ndim) )
  allocate( r0(ndim) )
  allocate(  a(ndim) )
  allocate(  b(ndim) )
  allocate(  p(ndim,max_projection_gcr) )
  allocate( Ap(ndim,max_projection_gcr) )

!------------------------------------------------------------------
!  Initialization

   x = zero
   r = zero
   p = zero

!  Initial residual r0 is b since x=0 (r = b-A*x = b).
!  If your linear system is (dRes/dq)*x=-Res(qk), then r0=-Res(qk).

   call gcr_rhs_vector(r0,ndim)
             r = r0
   rms_initial = gcr_rms(r0,ndim)

!   Preconditioning: p_1 = M^{-1}*r0.
    call gcr_preconditioner(r0,ndim, p(:,1),sweeps,istat_preconditioner_0)
    actual_sweeps = actual_sweeps + sweeps !Accumulate sweeps. 

    if (istat_preconditioner_0 == 100) then
     write(2000,*) ">>>>> Exit GCR: Preconditioner diverged at the beginning(roc > 100.0)."
     status_pc          = 2
     status_gcr         = 0
     actual_projections = 0
     actual_sweeps      = 0
     return
    endif

!   Compute Ap
    call gcr_compute_Ap(p(:,1),r0,ndim, Ap(:,1))

!-----------------------------------------------------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Beginning of GCR Projection
!---------------------------------------------------------------------
!---------------------------------------------------------------------
   i = 0
  gcr_projection : do

   i = i + 1
   actual_projections = i

!-----------------------------------------
!  Update x and r to minimize the residual in L2 norm.
!   Minimize (ri,ri): -> (r-alpha*Ap_i)*d(r-alpha*Ap_i)/d(alpha)
!                      = 0 -> alpha = (r,Ap_i)/(Ap_i,Ap_i).
!
!      alpha = (r,Ap_i)/(Ap_i,Ap_i)  <- inner products.
!            = (r/|r|,Ap_i/|Ap_i|)*|r|*|Ap_i| / |Ap_i|^2
!            = (r/|r|,Ap_i/|Ap_i|)*|r|/|Ap_i|
!          x = x + alpha*p_i
!          r = r - alpha*Ap_i
!
    length_r = gcr_length(r,ndim)
    length_p = gcr_length(Ap(:,i),ndim)
           a =       r/max(my_eps,length_r)
           b = Ap(:,i)/max(my_eps,length_p)
       alpha = gcr_dot_product(a,b,ndim)*( length_r/max(my_eps,length_p) )

        !---------------------------------------------------------------
        ! Check GCR stall for i > 1: Always perform the first iteration.
          if (i > 1) then

           !A stall detection.
            stall_condition = maxval(alpha*Ap(:,i)) < gcr_rms(r,ndim)
 
          endif
        !---------------------------------------------------------------

      if ( .not.stall_condition ) then
        x(:) = x(:) + alpha* p(:,i)
        r(:) = r(:) - alpha*Ap(:,i)
      endif

!-----------------------------------------
!  Compute the convergence rate: rms(current r)/rms(r0).

      rkr0_gcr = gcr_rms(r,ndim) / rms_initial

      write(2000,'(a14,i5, a19,f8.3,a27,i8)') " gcr proj = ", actual_projections, &
                    ": res reduction = ", rkr0_gcr, " - preconditioner relax ", sweeps


!-----------------------------------------
!  Exit GCR if any of the following is true.

  !Exit if tolerance is met.
     if (rkr0_gcr < mu_gcr) then
      status_pc          = 0
      status_gcr         = 0
      write(2000,*) " >>>>> Tolerance met... Exit GCR."
      exit gcr_projection
     endif

  !Exit if stalled.
     if ( stall_condition ) then
      status_pc          = 0
      status_gcr         = 2
      write(2000,*) " >>>>> GCR stall detected... Exit GCR."
      x = p(:,1) !Revert to the implicit-dc solver.
      exit gcr_projection
     endif

  !Exit if max_projection is reached.
     if (i == max_projection_gcr) then
      status_pc          = 0
      status_gcr         = 2
      write(2000,*) " >>>>> Maximum projection reached... Exit GCR."
      exit gcr_projection
     endif

!-----------------------------------------
!  Construct the next search direction.

!   Set the next search direction to be p_{i+1} = M^{-1}*r(i).
    call gcr_preconditioner(r,ndim, p(:,i+1),sweeps,istat_preconditioner)
    actual_sweeps = actual_sweeps + sweeps !Accumulate sweeps. 

    if (istat_preconditioner == 100) then
     status_pc          = 2
     status_gcr         = 0
     exit gcr_projection
    endif

!   Compute Ap
    call gcr_compute_Ap(p(:,i+1),r0,ndim, Ap(:,i+1))

  !---------------------------------------------------------
  ! Orthogonalize p_{i+1} against all previous vectors: p_1, p_2, ..., p_{i}:
  !  p_{i+1} = A_approx^{-1}*r + sum_{k=0,i} beta*p_{k}
  !     beta = -(Ap_{i+1},Ap_{k})/(Ap_{k},Ap_{k})

    Apip1 = Ap(:,i+1)

    do k = 1, i

     length_p = gcr_length(Ap(:,k),ndim)
            a = Apip1  /max(my_eps,length_p)
            b = Ap(:,k)/max(my_eps,length_p)
        beta = -gcr_dot_product(a,b,ndim)
     Ap(:,i+1) = Ap(:,i+1) + Ap(:,k)*beta
      p(:,i+1) =  p(:,i+1) +  p(:,k)*beta

    end do
  !---------------------------------------------------------

  end do gcr_projection
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! End of GCR Projection
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!------------------------------------------------------------------
! Store the solution x(:) in your code.

  call gcr_store_x_in_your_array(x,ndim)

!------------------------------------------------------------------
! Deallocate GCR arrays

  deallocate( x, r, r0, a, b, p, Ap )

 end subroutine gcr_solver
!********************************************************************************






!********************************************************************************
! Assign ndim, which is the total number of discrete unknowns.
! For example, if you have a cell-centered code for 3D Navier-Stokes equations,
! ndim = (number of cells) x 5(equatins) = ncells*5.
!********************************************************************************
 subroutine gcr_dimension(ndim)

 implicit none

 integer, intent(out) :: ndim

 !Compute the total number of discrete unknowns in your code, and return it to GCR.

  call your_dimension(ndim)

 end subroutine gcr_dimension
!********************************************************************************

!********************************************************************************
! GCR works with a global 1D array of residuals. So, please store your RHS
! in a 1D array, r(:), if not in 1D array.
!
! Note: If you solve dRes/dq*dq = -Res, then the rhs is -Res.
!********************************************************************************
 subroutine gcr_rhs_vector(r,ndim)

 implicit none

 integer                  , intent( in) :: ndim  !dimension of r
 real(p2), dimension(ndim), intent(out) :: r

 !Store your rhs in the 1D array r(:), and return it to GCR.

  call your_rhs_vector(r,ndim)

 end subroutine gcr_rhs_vector
!********************************************************************************

!********************************************************************************
! Perform preconditioning: compute p = M^{-1}*r.
!
! The only important output is p = M^{-1}*r.
!
! 'sweeps' may be 0 if the preconditioner is not based on relaxations.
! 'istat' may always be 0 if failure is not defined for a preconditioner.
!
! If your preconditioner consists of
!
! (1)Approximate Jacobian, M.
! (2)Approximate inversion of M by a relaxation scheme (e.g., Gauss Seidel).
!
! You can return p and the residual reduction 'roc', which is the
! ratio of the final linear residual to the initial one for M.p=r,
! and also the number of relaxations 'sweeps'. In this case, istat can be nonzero
! to indicate that the relaxation diverged.
!
! If you use iLU-type preconditioner, just compute p and return sweeps=0, and 
! istat = 0.
!
!********************************************************************************
 subroutine gcr_preconditioner(r,ndim, p, sweeps,istat)

 implicit none

 integer                  , intent( in) :: ndim  !dimension of r and p
 real(p2), dimension(ndim), intent( in) :: r     !residual vector
 real(p2), dimension(ndim), intent(out) :: p     !M^{-1}.r
 integer ,                  intent(out) :: istat !status integer used in HANIM
 integer ,                  intent(out) :: sweeps

!Local
 real(p2) :: roc   !r(final)/r(0)

! Your preconditioner, which computes p as p = M^{-1}.r.

  call your_preconditioner(r,ndim, p,roc,sweeps)

! Status is defined as follows.

  if (    roc > 100.0_p2) then !<- This is bad.
   istat = 100
  elseif (roc > 10.0_p2) then  !<- This is acceptable.
   istat = 10
  elseif (roc > 1.0_p2) then   !<- This is acceptable.
   istat = 1
  else                         !<- This is acceptable.
   istat = 0
  endif

 end subroutine gcr_preconditioner
!********************************************************************************

!********************************************************************************
! Compute Ap.
!
! (1)If you have A, then perform a matrix-vector multiplicaiton.
! (2)If you don't, then Frechet: Ap = dRes/dU*p = [Res(u+epsilon*p)-Res(u)]/eps.
!
!********************************************************************************
 subroutine gcr_compute_Ap(p,r0,ndim, Ap)

 implicit none

 integer                  , intent( in) :: ndim
 real(p2), dimension(ndim), intent( in) :: p
 real(p2), dimension(ndim), intent( in) :: r0
 real(p2), dimension(ndim), intent(out) :: Ap

 !Compute Ap and return it to GCR.

  call your_compute_Ap(p,r0,ndim, Ap)

 end subroutine gcr_compute_Ap
!********************************************************************************

!********************************************************************************
! Copy the GCR solution x to your solution data.
!********************************************************************************
 subroutine gcr_store_x_in_your_array(x,ndim)

 implicit none

 integer                  , intent(in) :: ndim  !dimension of x
 real(p2), dimension(ndim), intent(in) :: x

 !Given a GCR solution x, copy it to the solution array in your code.

  call your_store_x_in_your_array(x,ndim)

 !Note: Do this in your code. No need to return anything to GCR.

 end subroutine gcr_store_x_in_your_array
!********************************************************************************

!********************************************************************************
!* This function returns the dot (or inner) product of two vectors.
!********************************************************************************
 function gcr_dot_product(a,b,n) result(dot_prod)

 implicit none

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: a, b

!Output
 real(p2)                           :: dot_prod

 integer                            :: i

   dot_prod = zero

  do i = 1, n

   dot_prod = dot_prod + a(i)*b(i)

  end do

 end function gcr_dot_product
!********************************************************************************

!********************************************************************************
!* This function returns the root-mean-square of a vector. 
!********************************************************************************
 function gcr_rms(vec,n) result(rms)

 implicit none

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec

!Output
 real(p2)                           :: rms

 integer                            :: i

   rms = zero

  do i = 1, n

   rms = rms + vec(i)**2

  end do

   rms = sqrt( rms / real(n,p2) )

 end function gcr_rms
!********************************************************************************

!********************************************************************************
!* This function returns the length of a vector.
!********************************************************************************
 function gcr_length(vec,n) result(length)

 implicit none

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec
!Output
 real(p2)                           :: length

 integer                            :: i

   length = zero

  do i = 1, n

   length = length + vec(i)**2

  end do

   length = sqrt( length )

 end function gcr_length
!********************************************************************************

!********************************************************************************
!* This function returns a global vector of residual.
!********************************************************************************
 function rms(vec,n)

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter :: zero = 0.0_p2

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec
!Output
 real(p2)                           :: rms
!Local variable
 integer :: i

   rms = zero

  do i = 1, n
   rms = rms + vec(i)**2
  end do

   rms = sqrt( rms / real(n,p2) )

 end function rms
!********************************************************************************

!********************************************************************************
!* This function returns a global vector of residual.
!********************************************************************************
 function length(vec,n)

 implicit none

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter :: zero = 0.0_p2

!Input
 integer               , intent(in) :: n
 real(p2), dimension(n), intent(in) :: vec
!Output
 real(p2)                           :: length
!Local variable
 integer :: i

   length = zero

  do i = 1, n
   length = length + vec(i)**2
  end do

   length = sqrt( length )

 end function length
!********************************************************************************

!********************************************************************************
!********************************************************************************
!********************************************************************************
! End of General GCR subroutines.
!********************************************************************************
!********************************************************************************
!********************************************************************************


!--------------------------------------------------------------------------------
!
! Below, subroutine/functions are written based on the data specific to
! the edu3d_euler program. This file serves as an example of how to use
! the general GCR module. So, by replacing the subroutines below, you can use
! the above GCR subroutines to solve your problem in your code.
!
!--------------------------------------------------------------------------------

!********************************************************************************
!* To return the total dimension of the unknown vector.
!********************************************************************************
 subroutine your_dimension(ndim)

 use data_module, only : nnodes

 implicit none
 integer, intent(out) :: ndim

  ndim = nnodes*5

 end subroutine your_dimension

!********************************************************************************
!* To store the residuals times (-1) into the 1D array, r(1:ndim).
!********************************************************************************
 subroutine your_rhs_vector(r,ndim)

 use data_module, only : nnodes, res

 implicit none
 integer                  , intent( in) :: ndim  !dimension of r
 real(p2), dimension(ndim), intent(out) :: r

 integer                                :: i, k

  do i = 1, nnodes
   do k = 1, 5
    r( 5*(i-1) + k ) = -res(i,k)
   end do
  end do

 end subroutine your_rhs_vector

!********************************************************************************
!* Preconditioner:
!*
!* Perform p = Jac^{-1}*r, approximately by relaxing the linear system
!* Jac*p=r by the GS relaxation scheme to a specified tolerance.
!* This is a variable preconditioner, which can be employed in GCR (or flexible GMRES).
!*
!********************************************************************************
 subroutine your_preconditioner(r,ndim, p,roc,sweeps)

 use data_module  , only : nnodes, res, du
 use linear_solver, only : linear_relaxation
 use data_module  , only : lrelax_sweeps_actual, lrelax_roc, nnodes

!Input
 integer                  , intent( in) :: ndim  !dimension of r and p
 real(p2), dimension(ndim), intent( in) :: r     !residual vector

!Output
 real(p2), dimension(ndim), intent(out) :: p
 real(p2),                  intent(out) :: roc
 integer ,                  intent(out) :: sweeps

 real(p2), dimension(nnodes,5) :: r_temp
 integer                       :: i, k

 !Save off the current residual.

  r_temp = res

 !Store the given r into res. Note: need minus sign since my preconditioner will put
 !the minus sign in the right hand side.

  do i = 1, nnodes
   do k = 1, 5
    res(i,k) = -r( 5*(i-1) + k )
   end do
  end do

 !In edu3d_euler program, this subroutine performs the GS relaxation to aproximately
 !solve Jac*du=-res=r, and stores the solution in du, where Jac is the preconditioner matrix.

  call linear_relaxation

 !Store the solution into p(:):
 ! ->   p(:) = -(Jac)^{-1}*res = (Jac)^{-1}*r, where the inverse is approximate.

  do i = 1, nnodes
   do k = 1, 5
    p( 5*(i-1) + k ) = du(i,K)
   end do
  end do

! Copy back the current residual.

  res = r_temp

! Return the convergence rate (linear residual reduction), and the total
! number of linear relaxations that has been performed in 'linear_relaxation'.

     roc = lrelax_roc
  sweeps = lrelax_sweeps_actual

 end subroutine your_preconditioner

!********************************************************************************
!* Compute the Frechet derivative
!*
!* The matrix vector product, dR/du*p, is performed by using the Frechet
!* derivative as follows:
!*
!*    dR/du*p = [Res(u+epsilon*p) - Res(u)]/epsilon
!*            = [Res(u+epsilon*p/|p|) - Res(u)]/epsilon*|p|
!*
!* and pseudo time terms are added.
!*
!********************************************************************************
 subroutine your_compute_Ap(p,r0,ndim, Ap)

 use data_module, only : res, u, w, nnodes, my_eps
 use euler_eqns , only : u2w
 use residual   , only : compute_residual

!Input
 integer                  , intent( in) :: ndim  !dimension of x
 real(p2), dimension(ndim), intent( in) :: p
 real(p2), dimension(ndim), intent( in) :: r0
!Output
 real(p2), dimension(ndim), intent(out) :: Ap

 real(p2), dimension(nnodes,5)  :: r_temp, w_temp, u_temp
 integer                        :: i, k

 real(p2), dimension(ndim)      :: r_ep, v_over_dt
 real(p2)                       :: eps, length_p

!--------------------------------------------------------------------------------
! Save the current solution and residual

  r_temp = res
  u_temp = u
  w_temp = w

!--------------------------------------------------------------------------------
! Store (u + epsilon*p) into the solution vector u.

   length_p = length(p,ndim)

  !Store the solution temporarily in r_ep(:).
   do i = 1, nnodes
    do k = 1, 5
     r_ep( 5*(i-1) + k ) = u(i,k)
    end do
   end do  

   eps = max( one, rms(r_ep,ndim) )*sqrt(my_eps) ! See Ref.[2].

  !Store perturbed solutions into u and w for residual calculation.
   do i = 1, nnodes
    do k = 1, 5
     u(i,k) = u(i,k) + eps* p( 5*(i-1) + k )/max(my_eps,length_p)
    end do
     w(i,:) = u2w(u(i,:))
   end do

!--------------------------------------------------------------------------------
! Compute the residual vector Res(u + epsilon*p).

   call compute_residual

!--------------------------------------------------------------------------------
! Store the computed residual into the residual vector, r_ep.

   do i = 1, nnodes
    do k = 1, 5
     r_ep( 5*(i-1) + k ) = res(i,k)
    end do
   end do

!--------------------------------------------------------------------------------
! Compute the Frechet derivative (finite-difference approximation).
! Note: r0 = -Res, and so use (-r0) for Res.

   Ap = ( r_ep - (-r0) ) / eps * length_p ! = [Res(u+epsilon*p/|p|) - Res(u)]/epsilon*|p|

!--------------------------------------------------------------------------------
! Compute the pseudo-time term

   call vector_dt_term(v_over_dt,ndim)

   do i = 1, ndim
    Ap(i) = v_over_dt(i)*p(i) + Ap(i)  ! Ap <- (V/dt + A)*p
   end do

!--------------------------------------------------------------------------------
! Copy back the current solution and residual

   res = r_temp
     u = u_temp
     w = w_temp

 end subroutine your_compute_Ap

!********************************************************************************
!* This function returns a global vector of residual.
!********************************************************************************
 subroutine vector_dt_term(vector,ndim)

 use data_module           , only : vol, dtau, nnodes

 implicit none

 integer ,      parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 integer                  , intent( in) :: ndim  !dimension of x

!Output
 real(p2), dimension(ndim), intent(out) :: vector

!local variables
 integer                                :: i, k

  do i = 1, nnodes

   do k = 1, 5
    vector( 5*(i-1) + k ) = vol(i) / dtau(i)
   end do

  end do  

 end subroutine vector_dt_term
!********************************************************************************

!********************************************************************************
!* To store the computed GCR solution in the correction vector, du(1:nnodes,:)
!********************************************************************************
 subroutine your_store_x_in_your_array(x,ndim)

 use data_module, only : du, nnodes

 implicit none

 integer                  , intent(in) :: ndim  !dimension of x
 real(p2), dimension(ndim), intent(in) :: x     !residual vector

 integer :: i, k

  do i = 1, nnodes
   do k = 1, 5
    du(i,k) = x( 5*(i-1) + k ) 
   end do
  end do

 end subroutine your_store_x_in_your_array



 end module gcr


!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes subroutines that compute the linear LSQ coefficients
!* and compute the LSQ gradients at nodes.
!*
!*  - LSQ gradients are computed at a node by looping over its nghbrs, and
!*    accumulating [cx,cy,cz]*(wk-wj), where cx, cy, cz are the LSQ coefficients.
!*
!*  - It can be probably computed also by looping over edges (not implemented).
!*
!*  - To achieve 3rd-order accuracy, a quadratic LSQ fit is required.
!*    This is not implemented yet.
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

 module gradients

 !This module contains the following subroutine:

  public :: compute_gradient         ! compute the gradietns at nodes
  public :: compute_lsq_coefficeints ! compute the LSQ coefficients
  public :: linear_lsq3x3_coeff      ! compute the LSQ coefficient at a node
  public :: lsq_weight               ! LSQ weight function
  public :: qr_factorization         ! QR factorization to solve LSQ problem
  public :: gewp_solve               ! Gauss elimination used in QR

 contains

!********************************************************************************
!
! Compute the gradient at nodes.
!
! - Compute the gradient by [wx,wy,wz] = sum_nghbrs [cx,cy,cz]*(w_nghbr - wj),
!   where [cx,cy,cz] are the LSQ coefficients.
!
! - LSQ coefficients are stored in 2D arrays, e.g., for x-derivative,
!
!    node1 with 5 nghbrs : cx(1, 1), cx(1, 2), ..., cx(1, 5),
!    node2 with 6 nghbrs : cx(1, 6), cx(1, 7), ..., cx(1,11),
!    node3 with 4 nghbrs : cx(1,12), cx(1,13), ..., cx(1,15),
!    etc.
!
!   Loop over nghbrs can be performed by using nnghrs array, which has the
!   number of nghbr nodes for each node: e.g.,
!
!     nnghbrs(node1)=5, nnghbrs(node2)=6, nnghbrs(node3)=4, etc,
!
!   and k0 array, which has the starting index for each node:
!
!     k0(node1)=0, k0(node2)=5, k0(node3)=11, etc.
! 
!   The LSQ coefficeints can be accessed at a node i as
!
!   do k = 1, nnghbrs(i)
!    cx(i,k0(i)+k)
!    cy(i,k0(i)+k)
!    cz(i,k0(i)+k)
!   end do
!
!********************************************************************************
 subroutine compute_gradient

 use data_module, only : w, nnodes, gradw, nghbr, nnghbrs, k0
 use data_module, only : linear_lsq3x3_cx, linear_lsq3x3_cy, linear_lsq3x3_cz

 integer , parameter    :: p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter    :: zero = 0.0_p2

 integer                :: i, k, inghbr
 integer                :: ix=1, iy=2, iz=3
 real(p2), dimension(5) :: dw

  nodes : do i = 1, nnodes

   gradw(i,:,:) = zero

  !loop over the nghbrs
   do k = 1, nnghbrs(i)

           inghbr = nghbr(k0(i)+k)
               dw = w(inghbr,:) - w(i,:)
    gradw(i,:,ix) = gradw(i,:,ix) + linear_lsq3x3_cx(k0(i)+k)*dw
    gradw(i,:,iy) = gradw(i,:,iy) + linear_lsq3x3_cy(k0(i)+k)*dw
    gradw(i,:,iz) = gradw(i,:,iz) + linear_lsq3x3_cz(k0(i)+k)*dw

   end do

  end do nodes

 end subroutine compute_gradient


!********************************************************************************
! Compute the lsq coefficients for the gradient at nodes.
!
! - LSQ coefficients are stored in 2D arrays, e.g., for x-derivative,
!
!    node1 with 5 nghbrs : cx(1, 1), cx(1, 2), ..., cx(1, 5),
!    node2 with 6 nghbrs : cx(1, 6), cx(1, 7), ..., cx(1,11),
!    node3 with 4 nghbrs : cx(1,12), cx(1,13), ..., cx(1,15),
!    etc.
!
!   Loop over nghbrs can be performed by using nnghrs array, which has the
!   number of nghbr nodes for each node: e.g.,
!
!     nnghbrs(node1)=5, nnghbrs(node2)=6, nnghbrs(node3)=4, etc,
!
!   and k0 array, which has the starting index for each node:
!
!     k0(node1)=0, k0(node2)=5, k0(node3)=11, etc.
! 
!   The LSQ coefficeints can be accessed at a node i as
!
!   do k = 1, nnghbrs(i)
!    cx(i,k0(i)+k)
!    cy(i,k0(i)+k)
!    cz(i,k0(i)+k)
!   end do
!
!********************************************************************************
 subroutine compute_lsq_coefficeints

 use data_module, only : nghbrs_total, nnodes
 use data_module, only : linear_lsq3x3_cx, linear_lsq3x3_cy, linear_lsq3x3_cz

 implicit none

 integer :: i

  allocate( linear_lsq3x3_cx(nghbrs_total) )
  allocate( linear_lsq3x3_cy(nghbrs_total) )
  allocate( linear_lsq3x3_cz(nghbrs_total) )

  nodes : do i = 1, nnodes

   call linear_lsq3x3_coeff(i)

  end do nodes

 end subroutine compute_lsq_coefficeints


!********************************************************************************
! --- LSQ coefficients for Least-Squares Gradient Reconstruction ---
!
! Construct a matrix for the linear least-squares(LSQ) gradient reconstruction.
! The normal matrix, A^t*A, is 3x3. But we do not form the normal matrix.
! The nnghbrsx3 rectangular LSQ matrix is inverted by QR factorization.
! This method is known to be more robust than the normal matrix approach.
!
! ------------------------------------------------------------------------------
!  Input:  inode = node number
!
! Output:  [cx,cy,cz] = LSQ coefficietns, which allows us to compute the gradient
!          of u at j as [ux,uy,uz] = sum_nghbrs [cx,cy,cz]*(u_nghbr - uj).
! ------------------------------------------------------------------------------
!
! - LSQ coefficients are stored in 2D arrays, e.g., for x-derivative,
!
!    node1 with 5 nghbrs : cx(1, 1), cx(1, 2), ..., cx(1, 5),
!    node2 with 6 nghbrs : cx(1, 6), cx(1, 7), ..., cx(1,11),
!    node3 with 4 nghbrs : cx(1,12), cx(1,13), ..., cx(1,15),
!    etc.
!
!   Loop over nghbrs can be performed by using nnghrs array, which has the
!   number of nghbr nodes for each node: e.g.,
!
!     nnghbrs(node1)=5, nnghbrs(node2)=6, nnghbrs(node3)=4, etc,
!
!   and k0 array, which has the starting index for each node:
!
!     k0(node1)=0, k0(node2)=5, k0(node3)=11, etc.
! 
!   The LSQ coefficeints can be accessed at a node i as
!
!   do k = 1, nnghbrs(i)
!    cx(i,k0(i)+k)
!    cy(i,k0(i)+k)
!    cz(i,k0(i)+k)
!   end do
!
!********************************************************************************
 subroutine linear_lsq3x3_coeff(inode)

 use data_module, only : x, y, z, nnghbrs, k0, nghbr
 use data_module, only : linear_lsq3x3_cx, linear_lsq3x3_cy, linear_lsq3x3_cz

 implicit none

 integer, intent(in) :: inode

!Local variables

 integer , parameter :: p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter :: zero = 0.0_p2

 real(p2) :: dx, dy, dz, w
 integer  :: k, inghbr, ix=1,iy=2, iz=3

 real(p2), pointer, dimension(:,:) :: a
 real(p2), pointer, dimension(:,:) :: rinvqt
 integer :: m, n

  m = nnghbrs(inode)  ! Number of nghbrs
  n = 3               ! # gradient components: ux, uy, uz.

  allocate(a(m,n))      !mxn = nnghbrsx3 rectangular LSQ matrix.
  allocate(rinvqt(n,m)) !This is the inverse of the LSQ matrix.
  a = zero

!  Loop over the neighbor nodes to construct the mx3 rectangular LSQ matrix.
!
!      weight * ( uk = uj + ux*(xk-xj) + uy*(yk-yj) + uz*(zk-zj) )
!   -> weight*(xk-xj)*ux + weight*(yk-yj)*uy + weight*(zk-zj)*uz = weight*(uk-uj) for all k.

   do k = 1, m
    inghbr = nghbr( k0(inode)+k ) ! kth nghbr of inode.

          dx = x(inghbr) - x(inode)
          dy = y(inghbr) - y(inode)
          dz = z(inghbr) - z(inode)
           w = lsq_weight(dx, dy, dz)
     a(k,ix) = w*dx 
     a(k,iy) = w*dy
     a(k,iz) = w*dz

   end do

!  Solve the mx3 LSQ problem by QR factorization: A*x=b -> x = Ainv*b.
!  Ainv = Pseudo inverse = R^{-1}*Q^t.

   call qr_factorization(a,rinvqt,m,n)


!  Now compute the coefficients for neighbors.
!  A*x=b, where b = [ weight1*(u1-uj), weight2*(u2-uj), ...  ].
!  We compute and store: [cx,cy,cz] = Ainv*[ weight1, weight2, ... ],
!  so that [ux,uy,uz] = sum_k [cx,cy,cz]*(uk - uj).

     loop_nghbr : do k = 1, m
       inghbr = nghbr( k0(inode)+k ) ! kth nghbr of inode.

          dx = x(inghbr) - x(inode)
          dy = y(inghbr) - y(inode)
          dz = z(inghbr) - z(inode)
           w = lsq_weight(dx, dy, dz)

      linear_lsq3x3_cx(k0(inode)+k)  = rinvqt(ix,k) * w
      linear_lsq3x3_cy(k0(inode)+k)  = rinvqt(iy,k) * w
      linear_lsq3x3_cz(k0(inode)+k)  = rinvqt(iz,k) * w

     end do loop_nghbr

  deallocate(a,rinvqt)

 end subroutine linear_lsq3x3_coeff
!********************************************************************************

!********************************************************************************
! Compute the LSQ weight
!
! - Inverse distance to the power of "lsq_weight_power".
! - Unweighted lsq corresponds to lsq_weight_power = 0. 
!
!********************************************************************************
 function lsq_weight(dx, dy, dz)

 use input_parameter, only : lsq_weight_power

 implicit none

 integer , parameter ::  p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter :: one = 1.0_p2

 real(p2), intent(in) :: dx, dy, dz
 real(p2)             :: lsq_weight, distance

     distance = sqrt(dx*dx + dy*dy + dz*dz)

   lsq_weight = one / distance**lsq_weight_power

 end function lsq_weight

!****************************************************************************
! ------------------ QR Factorization ---------------------
!
!  This subroutine solves the LSQ problem: A*x=b, A=mxn matrix.
!
!  IN :       a = (m x n) LSQ matrix.  (m >= n)
!
! OUT :  rinvqt = R^{-1}*Q^t, which gives the solution as x = R^{-1}*Q^t*b.
!
!*****************************************************************************
 subroutine qr_factorization(a,rinvqt,m,n)

 implicit none

  integer , parameter ::   p2 = selected_real_kind(P=15) !Double precision
  real(p2), parameter :: zero = 0.0_p2
  real(p2), parameter ::  one = 1.0_p2
  real(p2), parameter ::  two = 2.0_p2

!Input
 integer ,                 intent( in) :: m, n
 real(p2), dimension(m,n), intent( in) :: a

!Output
 real(p2), dimension(n,m), intent(out) :: rinvqt

! Local variables
!
! Note: Think if you can reduce the number of
!       variables below to save memory.

 integer                  :: i, j, k

 real(p2), dimension(m,n) :: r
 real(p2)                 :: abs_rk, sign_rkk, wTw
 real(p2), dimension(m)   :: w, rk
 real(p2), dimension(m,m) :: qt, wwT

 real(p2), dimension(n,n) ::  r_nxn
 real(p2), dimension(n)   ::   y, b
 real(p2)                 ::    rhs

 if (m < n) then
  write(*,*) " Underdetermined system detected... m < n: "
  write(*,*) "   m =  ", m
  write(*,*) "   n =  ", n
  write(*,*) " qr_factorization() not designed to solve such a problem... Stop. "
  stop
 endif

!-------------------------------------------------------
! Initialization: R = A

   r = a

!-------------------------------------------------------
! Initialization: Qt = I

       qt = zero

  do i = 1, m
   qt(i,i) = one
  end do

!-------------------------------------------------------
! Apply reflection to each column of R, and generate
! the final upper triangular matrix R and the transpose
! Qt of the orthonormal matrix Q.

 column_loop : do k = 1, n

  !Our target are the elements below the (k,k) element
  !in k-th column, i.e., r(k:m).
  !So, rk gets shorter as we move on (as k increases).

        rk      = zero
        rk(k:m) = r(k:m,k)

  !Reflector Hk will zero out all the elements below r(k).

  !Compute the length of rk and the sign of the kth element.

         abs_rk = sqrt( dot_product(rk,rk) )
       sign_rkk = sign( one, rk(k) )

  !Define the reflecting vector w:   w = |rk|*(1,0,0,...,0)-rk
  !                               or w =-|rk|*(1,0,0,...,0)-rk
  !We switch the reflection (there are two possible ones)
  !to avoid w = 0 that can happen if rk=(1,0,0,...,0).

         w      = zero
         w(k)   = -sign_rkk*abs_rk
         w(k:m) = w(k:m) - rk(k:m)

  !Compute the length^2 of w: wt*w = [x,x,...,x]|x| = dot product = scalar.
  !                                             |x|
  !                                             |.|
  !                                             |.|
  !                                             |x|

    wTw = dot_product(w,w)

  !Compute the dyad of w: w*wt = |x|[x,x,...,x] = mxm matrix.
  !                              |x|
  !                              |.|
  !                              |.|
  !                              |x|

    do i = 1, m
     do j = 1, m
      wwT(i,j) = w(i)*w(j)
     end do
    end do

  !We now apply the reflector matrix Hk = I-2*wwt/wTw,
  !and update R and Qt.

  !Update  R:  R = Hk*R  = (I-2*wwt/wTw)*R  = R-2*(wwt*R)/wTw

   r  =  r - two*matmul(wwT,r)/wTw

  !Update Qt: Qt = Hk*Qt = (I-2*wwt/wTw)*Qt = Qt-2*(wwt*Qt)/wTw

   qt = qt - two*matmul(wwT,qt)/wTw

 end do column_loop

!-------------------------------------------------------
! Compute rinvqt(1:n,1:m) = R_{nxn}^{-1} * Q_{nxm}^t by
! solving R_{nxn} * rinvqt(1:n,k) = Q_{nxm}^t(1:n,k)
! for k=1,n. We can solve it easily by back substitution
! since R_{nxn} is upper triangular.

   r_nxn =  r(1:n,1:n)

   do k = 1, m

    !Solve r*y = b, where y is the k-th column of rinvqt.

        b = qt(1:n,k)

    !Solve the lower right equation.

      rhs = b(n)
     y(n) = rhs/r_nxn(n,n)

    !Go up and solve.
     do i = n-1, 1, -1

     !Take all known parts (j=i+1,n) to the rhs.

      !RHS is known, of course.
        rhs = b(i)
      !Below are all known since the solutions y(j=i+1,n) has already been solved.
       do j = i+1, n
        rhs = rhs - r_nxn(i,j)*y(j)
       end do

     !Divide the rhs by the coefficient of the (i,i) part.
      y(i) = rhs/r_nxn(i,i)

     end do

   !The soluton x is the k-th column of rinvqt.
    rinvqt(:,k) = y(:)

   end do

 end subroutine qr_factorization

!****************************************************************************
!* ------------------ GAUSS ELIMINATION WITH PIVOTING ---------------------
!*
!*  This computes the inverse of an (nm)x(nm) matrix "ai" and also
!*  computes the solution to a given lienar system.
!*
!*  IN :       ai = An (nm)x(nm) matrix whoise inverse is sought.
!*             nm = The size of the matrix "ai"
!*
!* OUT :
!*        inverse = the inverse of "ai".
!*       idetstat = 0 -> inverse successfully computed
!*                  1 -> THE INVERSE DOES NOT EXIST (det=0).
!*                  2 -> No unique solutions exist.
!*****************************************************************************
  subroutine gewp_solve(ai,inverse,idetstat,nm)

  implicit none

 integer , parameter ::    p2 = selected_real_kind(15) ! Double precision
 real(p2), parameter ::  zero = 0.0_p2
 real(p2), parameter ::   one = 1.0_p2

 integer ,                   intent( in) :: nm
 real(p2), dimension(nm,nm), intent( in) :: ai

 real(p2), dimension(nm,nm), intent(out) :: inverse
 integer ,                   intent(out) :: idetstat

 real(p2), dimension(nm,nm+1) :: a
 real(p2), dimension(nm)      :: x
 integer , dimension(nm)      :: nrow
 integer                      :: I,J,K,pp,m

  do m = 1, nm
!*****************************************************************************
!* Set up the matrix a
!*****************************************************************************

       do J=1,nm
        do I=1,nm
          a(I,J) = ai(I,J)
        end do
       end do

       do k=1,nm
          a(k,nm+1)=zero; nrow(k)=k
       end do
          a(m,nm+1)=one

!*****************************************************************************
!* HONA IKOKA..... 
!*****************************************************************************
    do j=1,nm-1
!*****************************************************************************
!* FIND SMALLEST pp FOR a(pp,j) IS MAXIMUM IN JTH COLUMN.
!***************************************************************************** 
      call findmax(nm,j,pp,a,nrow)
!*****************************************************************************
!* IF a(nrow(p),j) IS zero, THERE'S NO UNIQUE SOLUTIONS      
!*****************************************************************************
      if (abs(a(nrow(pp),j)) < epsilon(one)) then
       write(6,*) 'THE INVERSE DOES NOT EXIST.'
        idetstat = 1
        return
      endif
!*****************************************************************************
!* IF THE MAX IS NOT A DIAGONAL ELEMENT, SWITCH THOSE ROWS       
!*****************************************************************************
      if (nrow(pp) .ne. nrow(j)) then
      call switch(nm,j,pp,nrow)
      else
      endif

!*****************************************************************************
!* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
!***************************************************************************** 
      call eliminate_below(nm,j,a,nrow)

    end do
!*****************************************************************************
!* CHECK IF a(nrow(N),N)=0.0 .
!*****************************************************************************
      if (abs(a(nrow(nm),nm)) < epsilon(one)) then
        write(6,*) 'NO UNIQUE SOLUTION EXISTS!'
        idetstat = 2
        return
      else
      endif
!*****************************************************************************
!* BACKSUBSTITUTION!
!*****************************************************************************
      call backsub(nm,x,a,nrow)
!*****************************************************************************
!* STORE THE SOLUTIONS, YOU KNOW THEY ARE INVERSE(i,m) i=1...
!*****************************************************************************
      do i=1,nm
         inverse(i,m)=x(i)
      end do
!*****************************************************************************
  end do

      idetstat = 0

    return

!*****************************************************************************
 end subroutine gewp_solve

!*****************************************************************************
!* Four subroutines below are used in gewp_solve() above.
!*****************************************************************************
!* FIND MAXIMUM ELEMENT IN jth COLUMN 
!***************************************************************************** 
      subroutine findmax(nm,j,pp,a,nrow)

      implicit none

      integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
      integer , intent( in) :: nm
      real(p2), intent( in) :: a(nm,nm+1)
      integer , intent( in) :: j,nrow(nm)
      integer , intent(out) :: pp
      real(p2)              :: max
      integer               :: i

            max=abs(a(nrow(j),j)); pp=j

           do i=j+1,nm

             if (max < abs(a(nrow(i),j))) then

                  pp=i; max=abs(a(nrow(i),j))

             endif

           end do

      return

      end subroutine findmax
!*****************************************************************************
!* SWITCH THOSE ROWS       
!*****************************************************************************
      subroutine switch(nm,j,pp,nrow)

      implicit none

      integer, intent(   in) :: nm,j,pp
      integer, intent(inout) :: nrow(nm)
      integer                :: ncopy

      if (nrow(pp).ne.nrow(j)) then

         ncopy=nrow(j)
         nrow(j)=nrow(pp)
         nrow(pp)=ncopy

      endif

      return

      end subroutine switch
!*****************************************************************************
!* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
!*(Give me j, the column you are working on now)
!***************************************************************************** 
      subroutine eliminate_below(nm,j,a,nrow)

      implicit none

      integer , parameter     :: p2 = selected_real_kind(15) ! Double precision
      real(p2), parameter     :: zero = 0.0_p2
      integer , intent(   in) :: nm
      real(p2), intent(inout) :: a(nm,nm+1)
      integer , intent(   in) :: j,nrow(nm)
      real(p2)                :: m
      integer                 :: k,i

      do i=j+1,nm

        m=a(nrow(i),j)/a(nrow(j),j)
        a(nrow(i),j)=zero

          do k=j+1,nm+1
            a(nrow(i),k)=a(nrow(i),k)-m*a(nrow(j),k)
          end do

      end do

      return

      end subroutine eliminate_below
!*****************************************************************************
!* BACKSUBSTITUTION!
!*****************************************************************************
      subroutine backsub(nm,x,a,nrow)

      implicit none

      integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
      real(p2), parameter   :: zero = 0.0_p2

      integer , intent( in) :: nm
      real(p2), intent( in) :: a(nm,nm+1)
      integer , intent( in) :: nrow(nm)
      real(p2), intent(out) :: x(nm)
      real(p2)              :: sum
      integer               :: i,k

      x(nm)=a(nrow(nm),nm+1)/a(nrow(nm),nm)

      do i=nm-1,1,-1

         sum=zero

           do k=i+1,nm

              sum=sum+a(nrow(i),k)*x(k)

           end do

      x(i)=(a(nrow(i),nm+1)-sum)/a(nrow(i),i)

      end do

      return

      end subroutine backsub
!*********************************************************************

 end module gradients
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

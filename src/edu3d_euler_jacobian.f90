!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes subroutines that comptutes the residual Jacobian
!* for the node-centered edge-based discretization. The Jacobian is exact
!* for the first-order node-centered edge-based discretization. So, Newton's
!* method can be constructed for the first-order scheme. This Jacobian can serve
!* as a preconditioning matrix for Jacobian-Free Newton-Krylov solver for
!* higher-order edge-based schemes.
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

 module jacobian

  implicit none

 !This module contains the following subroutine:

  public :: compute_jacobian ! compute the jacobian
  public :: kth_nghbr        ! find a nghbr index
  public :: gewp_solve       ! Gauss elimination to invert the diagonal blocks

 contains

!********************************************************************************
! Compute the residual Jacobian exactly for the first-order residual.
!
! - Jac = [V/dtau+dRes1/dU] for the linear system: Jac*dU = -Residual.
!
!    where Res1     = 1st-order accurate residual (compact stencil),
!          Residual = 1st/2nd-order residual.
!          V/dtau   = global diagonal matrix of dual_vol/dtau.
!
! - This is the linear system that needs to be solved per nonlinear iteration.
! - This is very similar to the residual subroutine.
! - Exact flux Jacobian will be computed by Automatic Differentiation.
! - Diagonal blocks are inverted and stored at the end.
!
! - Off-diagonal blocks are stored as follows:
!
!     node1 with 5 nghbrs : jac_off( 1,:,:), jac_off( 2,:,:), ..., jac_off( 5,:,:),
!     node2 with 6 nghbrs : jac_off( 6,:,:), jac_off( 7,:,:), ..., jac_off(11,:,:),
!     node3 with 4 nghbrs : jac_off(12,:,:), jac_off(13,:,:), ..., jac_off(15,:,:),
!     etc.
!
!    Loop over nghbrs can be performed by using nnghrs array, which has the
!    number of nghbr nodes for each node: e.g.,
!
!      nnghbrs(node1)=5, nnghbrs(node2)=6, nnghbrs(node3)=4, etc,
!
!    and k0 array, which has the starting index for each node:
!
!      k0(node1)=0, k0(node2)=5, k0(node3)=11, etc.
! 
!    The off-diagonal blocks can be accessed at a node j as
!
!     do k = 1, nnghbrs(j)
!      jac_off(i,k0(j)+k,:,:)
!     end do
!
! - The Jacobian matirx is stored in two parts: diagonal and off-diagonal blocks.
!
!       For example, consider the residual at node j, Rj, which is a vector,
!       and suppose j has 5 neighbor nodes k=(1, 2, 3, 4, 5).
!       For implicit defect-correction solver, the residual is taken as 1st-order scheme,
!       so that Res_j is a function of Uj, U1, U2, U3, U4, and U5 (compact stencil).
!
!       Then, the j-th component of [V/dt+dR/dU]*dU = -Residual is
!
!       [        dRj/dU1]*dU1 + [dRj/dU2]*dU2 + [dRj/dU3]*dU3 +
!       [V/dtauj+dRj/dUj]*dUj + [dRj/dU4]*dU4 + [dRj/dU5]*dU5 = - Rj
!
!       where dtauj is the local time step at j, dRj/dUj is the derivative of Rj
!       with respect to the solution at j, i.e., Uj, and similarly for others.
!       The Gauss-Seidel relaxation can be written as
!
!       dUj = [V/dtj+dR/dUj]^{-1}*(- [dRj/dU1]*dU1 - [dRj/dU2]*dU2
!                                  - [dRj/dU3]*dU3 - [dRj/dU4]*dU4
!                                  - [dRj/dU5]*dU5 - Res_j )
!
!       To perform this operation, we store two types of data:
!
!        1.  5x5     Diagonal block:  jac_diag(j, :,:) = [V/dtauj+dRj/dUj]
!        2.  5x5 Off-Diagonal Block:  jac_off(j,k,:,:) = [        dRj/dUk], k = 1,2,3,4,5
!
!       so, we can perform the relaxation at every node by looping over the neighbors:
!
!         b = -Res_j
!        do k = 1, nnghbrs(j)
!         b = b - jac_off(i,k0(i)+k,:,:)*dU( nghbr(k0(i)+k), :)
!        end do
!         dUj = jac_diag(j, :,:)^{-1}*b,
!
!       where nghbr(k0(i)+k) stores the actual node number for kth nghbr.
!
!       Actual relaxation will be performed in edu3d_euler_linear_solver.f90.
!
!       To compute dR/dUj and dR/dUk, just like the residual is computed as a sum of
!       numerical fluxes over edges, we compute the Jacobian blocks as a sum of
!       derivatives of the numerical flux over edges. So, we loop over edges,
!       compute the derivatives of the numerical flux with respect to each end node
!       of the edge, and accumulate the contribution at the two nodes.
!       At the end of the edge loop, at every node j, we have
!
!               dRj/dUj = sum_over_edges ( dflux/dUj )
!               dRj/dUk = dflux/dUk
!
!       which forms the j-th row of the global Jacobian. Here, 'flux' is the
!       numerical flux between nodes j and k.
!
! - The flux Jacobian, dflux/dUj, is computed by automatic differentiation.
!   Numerical flux subroutines in edu3d_flux_functions_ddt.f90 are all written
!   by automatic differentiation. They return the flux vector as well as
!   the derivative. This is very convenient, and it makes it so easy to try out
!   a new numerical flux without coding the Jacobian subroutine.
!   See edu3d_flux_functions_ddt.f90 and edu3d_euler_fluxjac.f90 to learn
!   how it is implemented.
!
!********************************************************************************
 subroutine compute_jacobian

 use data_module     , only : p2, zero, half, third, fourth
 use data_module     , only : nnodes, x, y, z, nedges, edge, w
 use data_module     , only : n12, n12_mag, e12, e12_mag, jac_diag, jac_off, jac_dinv
 use data_module     , only : ntria, tria, bn, bn_mag, bnj, bc_type
 use data_module     , only : k0, kth_nghbr_of_1, kth_nghbr_of_2
 use data_module     , only : vol, dtau

 use fluxjac_interface  , only : interface_jac, interface_jac_b
 use bc_states, only : get_right_state

 implicit none

!Local variables

 real(p2), dimension(3)   :: da    , ev
 real(p2)                 :: da_mag, ev_mag

 real(p2), dimension(5)   :: wb
 real(p2), dimension(5,5) :: duRduL

 integer :: i, j, k, ib
 integer :: node1, node2

 integer , dimension(3,2) :: loc_nghbr
 real(p2), dimension(5,5) :: dFnduL, dFnduR

 integer                  :: n1, n2, n3
 real(p2), dimension(5,5) :: dF1ndu1, dF2ndu2, dF3ndu3

 integer                  :: idestat
 real(p2), dimension(5,5) :: inverse

!     1             3
!      o------------o
!       \         .
!        \       . 
!         \    .
!          \ .
!           o
!           2
!
!
! Neighbor nodes of a given node within a triangle (boundary elements).

  loc_nghbr(1,1) = 2 !nghbr1 of vertex 1.
  loc_nghbr(1,2) = 3 !nghbr2 of vertex 1.

  loc_nghbr(2,1) = 3 !nghbr1 of vertex 2.
  loc_nghbr(2,2) = 1 !nghbr2 of vertex 2.

  loc_nghbr(3,1) = 1 !nghbr1 of vertex 3.
  loc_nghbr(3,2) = 2 !nghbr2 of vertex 3.

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Jacobian computation: interior edges

! Initialization

   jac_diag = zero
   jac_off  = zero

!--------------------------------------------------------------------------------
! Flux computation across internal edges (to be accumulated in res(:))
!
!   node2              1. Extrapolate the solutions to the edge-midpoint
!       o                 from the nodes, n1 and n2.
!        \             2. Compute the numerical flux
!         \            3. Add it to the residual for n1, and subtract it from
!          \              the residual for n2.
!           \
!            o         Directed area is the sum of all dual faces defined
!          node1       by the edge-midpoint, centroids of adjacent tetrahedra
!                      and the centroids of the adjacent triangular faces.
!                      Directed area is positive in n1 -> n2
!--------------------------------------------------------------------------------
  loop_edges : do i = 1, nedges


 ! Left and right nodes of the i-th edge

     node1 = edge(i,1)  ! Left node of the edge
     node2 = edge(i,2)  ! Right node of the edge
       da = n12(i,:)    ! This is the directed area vector (unit vector)
   da_mag = n12_mag(i)  ! Magnitude of the directed area vector
       ev = e12(i,:)    ! Edge vector (unit vector)
   ev_mag = e12_mag(i)  ! Magnitude of the edge vector

 !  Compute the flux Jacobians for given w1 and w2 (nodal values)

   call interface_jac( w(node1,:), w(node2,:), da, dFnduL, dFnduR  )

 !  Add the Jacobian multiplied by the magnitude of the directed area vector to node1.

     jac_diag(node1,:,:)        = jac_diag(node1,:,:)        + (dFnduL) * da_mag
     k = kth_nghbr_of_1(i)
     jac_off( k0(node1)+k ,:,:) = jac_off( k0(node1)+k, :,:) + (dFnduR) * da_mag

 !  Subtract the Jacobian multiplied by the magnitude of the directed area vector from node2.
 !  NOTE: Subtract because the outward face normal is -da for the node2.

     jac_diag(node2,:,:)        = jac_diag(node2,:,:)        - (dFnduR) * da_mag
     k = kth_nghbr_of_2(i)
     jac_off( k0(node2)+k, :,:) = jac_off( k0(node2)+k, :,:) - (dFnduL) * da_mag

  end do loop_edges
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Close the Jacobian by the derivative of the linearly-exact boundary flux quadrature.
!
!     l3             2
!      o------------o
!       \         .
!        \       .
!         \    .
!          \ .
!           o
!           1
!
! - Linearly-exact quadrature formulas for 3D tet/prs/pyr/hex grids can be found in:
!
! H. Nishikawa, "Beyond Interface Gradient: A General Principle for
! Constructing Diffusion Schemes", AIAA Paper 2010-5093, 40th AIAA Fluid
! Dynamics Conference and Exhibit, Chicago, 2010.
! http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_AIAA-2010-5093.pdf
! (NOTE: See Appendix B for boundary flux quadratures for 2D/3D elements.)
!
! Here, the formula for tetrahedral grids is implemented.
!
!--------------------------------------------------------------------------------

  boundary_element_loop : do i = 1, ntria

    ib = tria(i,4)

    !----------------------------------------------------------------
    ! Node 1: Phi_1(bn): BC is implemented here.

     n1 = tria( i, 1 )
     call get_right_state(  bc_type(ib),w(n1,:),bn(i,:),x(n1),y(n1),z(n1),  duRduL,wb )
     call interface_jac_b( w(n1,:), wb, bn(i,:), dFnduL, dFnduR )

     dF1ndu1 = dFnduL + matmul(dFnduR,duRduL)

    !----------------------------------------------------------------
    ! Node 2: Phi_2(bn): BC is implemented here.

     n2 = tria( i, 2 )
     call get_right_state(  bc_type(ib),w(n2,:),bn(i,:),x(n2),y(n2),z(n2),  duRduL,wb )
     call interface_jac_b( w(n2,:), wb, bn(i,:), dFnduL, dFnduR )

     dF2ndu2 = dFnduL + matmul(dFnduR,duRduL)

    !----------------------------------------------------------------
    ! Node 3: Phi_3(bn): BC is implemented here.

     n3 = tria( i, 3 )
     call get_right_state(  bc_type(ib),w(n3,:),bn(i,:),x(n3),y(n3),z(n3),  duRduL,wb )
     call interface_jac_b( w(n3,:), wb, bn(i,:), dFnduL, dFnduR )

     dF3ndu3 = dFnduL + matmul(dFnduR,duRduL)

    !----------------------------------------------------------------
    !----------------------------------------------------------------
    ! Contribution to node 1 

     !Diagonal block    : dRes(n1)/dU(n1)
      jac_diag(n1,:,:) = jac_diag(n1,:,:) + (6.0_p2/8.0_p2)*(dF1ndu1) * (third*bn_mag(i))

     !Off-diagonal block: dRes(n1)/dU(n2)
      k = kth_nghbr(n1,n2)
      jac_off( k0(n1)+k ,:,:) = jac_off( k0(n1)+k, :,:) + (1.0_p2/8.0_p2)*(dF2ndu2) * (third*bn_mag(i))

     !Off-diagonal block: dRes(n1)/dU(n3)
      k = kth_nghbr(n1,n3)
      jac_off( k0(n1)+k ,:,:) = jac_off( k0(n1)+k, :,:) + (1.0_p2/8.0_p2)*(dF3ndu3) * (third*bn_mag(i))

    !----------------------------------------------------------------
    !----------------------------------------------------------------
    ! Contribution to node 2 

     !Diagonal block    : dRes(n2)/dU(n2)
      jac_diag(n2,:,:) = jac_diag(n2,:,:) + (6.0_p2/8.0_p2)*(dF2ndu2) * (third*bn_mag(i))

     !Off-diagonal block: dRes(n2)/dU(n3)
      k = kth_nghbr(n2,n3)
      jac_off( k0(n2)+k ,:,:) = jac_off( k0(n2)+k, :,:) + (1.0_p2/8.0_p2)*(dF3ndu3) * (third*bn_mag(i))

     !Off-diagonal block: dRes(n2)/dU(n1)
      k = kth_nghbr(n2,n1)
      jac_off( k0(n2)+k ,:,:) = jac_off( k0(n2)+k, :,:) + (1.0_p2/8.0_p2)*(dF1ndu1) * (third*bn_mag(i))

    !----------------------------------------------------------------
    !----------------------------------------------------------------
    ! Contribution to node3

     !Diagonal block    : dRes(n3)/dU(n3)
      jac_diag(n3,:,:) = jac_diag(n3,:,:) + (6.0_p2/8.0_p2)*(dF3ndu3) * (third*bn_mag(i))

     !Off-diagonal block: dRes(n3)/dU(n1)
      k = kth_nghbr(n3,n1)
      jac_off( k0(n3)+k ,:,:) = jac_off( k0(n3)+k, :,:) + (1.0_p2/8.0_p2)*(dF1ndu1) * (third*bn_mag(i))

     !Off-diagonal block: dRes(n3)/dU(n2)
      k = kth_nghbr(n3,n2)
      jac_off( k0(n3)+k ,:,:) = jac_off( k0(n3)+k, :,:) + (1.0_p2/8.0_p2)*(dF2ndu2) * (third*bn_mag(i))

    !----------------------------------------------------------------
    !----------------------------------------------------------------

  end do boundary_element_loop

!--------------------------------------------------------------------------------
! Add pseudo time term vol/dtau to the diagonals of the diagonal blocks
! to form: Jacobian = V/dtau + dRes1/dU, where V/dtau is a global diagonal
! matrix having vol(i)/dtau(i) for each node, and Res1 is the first-order
! residual.

  do i = 1, nnodes
   do k = 1, 5
    jac_diag(i,k,k) = jac_diag(i,k,k) + vol(i)/dtau(i)
   end do
  end do

!--------------------------------------------------------------------------------
! Invert the diagonal blocks.
!  - I will invert the diagonal blocks here, and store them.
!  - Is it inefficient? Think about it.

  invert_at_node : do i = 1, nnodes

     inverse = zero
     idestat = 0

 !  Invert the diagonal block at node i by Gauss elimination with pivoting.
 !  Note: gewp_solve() actually solves a linear system, Ax=b, but here
 !        we use it only to obtain the inverse of A. So, b is a dummy.

 !                   A                dim  A^{-1}
    call gewp_solve( jac_diag(i,:,:), 5  , inverse, idestat )

 !  Report errors

    if (idestat/=0) then
     write(*,*) " Error in inverting the diagonal block... Stop"
     write(*,*) "  Node number = ", i
     write(*,*) "  Location    = ", x(i), y(i), z(i)
     do k = 1, 5
      write(*,'(12(es8.1))') ( jac_diag(i, k,j), j=1,5 )
     end do
     stop 127
    endif

 !  Save the inverse in jac_dinv.

    jac_dinv(i,:,:) = inverse

  end do invert_at_node

 end subroutine compute_jacobian
!********************************************************************************


!********************************************************************************
! This function is useful to find 'k', such that 
! node 'n2' is the k-th nghbr of the node 'n1'.
!********************************************************************************
 function kth_nghbr(n1,n2)

 use data_module, only : nghbr, k0, nnghbrs

 implicit none

 integer, intent(in) :: n1, n2

 integer :: k
 integer :: kth_nghbr
 logical :: found

  kth_nghbr = 0
      found = .false.

 !Find k, such that n2 is the k-th neghbor of n1.

  n1_nghbr_loop : do k = 1, nnghbrs(n1)

   if ( nghbr(k0(n1)+k) == n2 ) then

    kth_nghbr = k
    found = .true.
    exit n1_nghbr_loop

   endif

  end do n1_nghbr_loop

  if (.not. found) then
   write(*,*) " Neighbor not found.. Error. Stop"
   stop 127
  endif

 end function kth_nghbr


!****************************************************************************
!* ------------------ GAUSS ELIMINATION WITH PIVOTING ---------------------
!*
!*  This computes the inverse of an (nm)x(nm) matrix "ai" and also
!*  computes the solution to a given lienar system.
!*
!*  IN :       ai = An (nm)x(nm) matrix whoise inverse is sought.
!*             bi = A vector of (nm): Right hand side of the linear sytem
!*             nm = The size of the matrix "ai"
!*
!* OUT :
!*            sol = Solution to the linear system: ai*sol=bi
!*        inverse = the inverse of "ai".
!*       idetstat = 0 -> inverse successfully computed
!*                  1 -> THE INVERSE DOES NOT EXIST (det=0).
!*                  2 -> No unique solutions exist.
!*****************************************************************************
  subroutine gewp_solve(ai,nm, inverse,idetstat)

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


 end module jacobian
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

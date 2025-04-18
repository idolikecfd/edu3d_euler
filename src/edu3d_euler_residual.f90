!********************************************************************************
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D Euler
!*
!* This module containes a subroutine that comptutes the residuals at nodes for
!* the node-centered edge-based discretization. A special boundary flux quadrature
!* is required to preserve second- and third-order accuracy. Numerical flux
!* is computed at the edge-midpoint: edge-based quadrature.
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

 module residual

 !This module contains the following subroutine:
 
  public :: compute_residual ! compute the residual

 contains

!********************************************************************************
!
! This subtourine computes the residual at all nodes by the node-centered 
! edge-based discretization.
!
! The residual is computed in a loop over edges with a numerical flux computed
! at the edge midpoint, followed by the boundary closure.
!
! - The method is second-order accurate on arbitrary tetrahedral grids.
! - Second-order on regular prism/hex and mixed grids.
! - Requires a special boundary quadrature to preserve the design accuracy.
!
! Once you understand this subroutine, you'll see that a cell-centered method
! would be a bit simpler than this, especially the boundary closure part.
!
!********************************************************************************
 subroutine compute_residual

 use data_module      , only : p2, zero, half, third, fourth, ix, iy, iz
 use data_module      , only : x, y, z, nedges, edge, w, gradw
 use data_module      , only : res, wsn, n12, n12_mag, e12, e12_mag
 use data_module      , only : ntria, tria, bn, bn_mag, bnj, bc_type
 use data_module      , only : n_residual_evaluation

 use fluxjac_interface, only : interface_flux, interface_flux_b
 use bc_states        , only : get_right_state
 use gradients        , only : compute_gradient

 use input_parameter, only : accuracy_order

 implicit none

!Local variables
 real(p2), dimension(5)   :: num_flux
 real(p2), dimension(5)   :: num_flux_1, num_flux_2, num_flux_3

 real(p2), dimension(3)   :: da    , ev
 real(p2)                 :: da_mag, ev_mag

 real(p2)                 :: eig

 real(p2), dimension(5)   :: wb
 real(p2), dimension(5,5) :: dummy5x5

 integer :: i, j, ib
 integer :: node1, node2
 integer :: n1, n2, n3

 integer, dimension(3,2)  :: loc_nghbr

! Counting the number of residual evalutions.

   n_residual_evaluation = n_residual_evaluation + 1

!     1             3
!      o------------o
!       \         .
!        \       . 
!         \    .
!          \ .
!           o
!           2
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
! Residual computation: interior edges

! Initialization

   res = zero
   wsn = zero

! Compute gradients at all nodes for 2nd-order accuracy.

   if (accuracy_order == 2) call compute_gradient

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

 !  Compute the numerical flux for given w1 and w2.
 !  Reconstruction to the midpoint is performed in 'interface_flux' for 2nd-order.

   call interface_flux( w(node1,:), w(node2,:), gradw(node1,:,:), gradw(node2,:,:), &
                                                        da, ev, ev_mag, num_flux, eig   )

 !  Add the flux multiplied by the magnitude of the directed area vector to node1,
 !  and accumulate the max wave speed quantity for use in the pseudo-time step calculation.

     res(node1,:) = res(node1,:)  +  num_flux * da_mag
     wsn(node1  ) = wsn(node1  )  +       eig * da_mag

 !  Subtract the flux multiplied by the magnitude of the directed area vector from node2,
 !  and accumulate the max wave speed quantity for use in the pseudo-time step calculation.
 !  NOTE: Subtract because the outward face normal is -da for the node2.

     res(node2,:) = res(node2,:)  -  num_flux * da_mag
     wsn(node2  ) = wsn(node2  )  +       eig * da_mag


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
     call get_right_state(  bc_type(ib),w(n1,:),bn(i,:),x(n1),y(n1),z(n1),  dummy5x5,wb )
     call interface_flux_b( w(n1,:), wb, bn(i,:), num_flux_1, eig )

     wsn(n1) = wsn(n1)  +  eig * bn_mag(i)

    !----------------------------------------------------------------
    ! Node 2: Phi_2(bn): BC is implemented here.

     n2 = tria( i, 2 )
     call get_right_state(  bc_type(ib),w(n2,:),bn(i,:),x(n2),y(n2),z(n2),  dummy5x5,wb )
     call interface_flux_b( w(n2,:), wb, bn(i,:), num_flux_2, eig )

     wsn(n2) = wsn(n2)  +  eig * bn_mag(i)

    !----------------------------------------------------------------
    ! Node 3: Phi_3(bn): BC is implemented here.

     n3 = tria( i, 3 )
     call get_right_state(  bc_type(ib),w(n3,:),bn(i,:),x(n3),y(n3),z(n3),  dummy5x5,wb )
     call interface_flux_b( w(n3,:), wb, bn(i,:), num_flux_3, eig )

     wsn(n3) = wsn(n3)  +  eig * bn_mag(i)

    !----------------------------------------------------------------
    ! Add the boundary fluxes to residuals by the classical linearly-exact quadrature formula.
    
     res(n1,:) = res(n1,:) + ( (6.0_p2/8.0_p2)*num_flux_1 + (1.0_p2/8.0_p2)*(num_flux_2+num_flux_3 ) )*(third*bn_mag(i))
     
     res(n2,:) = res(n2,:) + ( (6.0_p2/8.0_p2)*num_flux_2 + (1.0_p2/8.0_p2)*(num_flux_3+num_flux_1 ) )*(third*bn_mag(i))

     res(n3,:) = res(n3,:) + ( (6.0_p2/8.0_p2)*num_flux_3 + (1.0_p2/8.0_p2)*(num_flux_1+num_flux_2 ) )*(third*bn_mag(i))

  end do boundary_element_loop


!  write(*,*) " i = 276, (x,y,z)=  ", x(276), y(276), z(276)
!  write(*,*) " i = 276, res(1) =  ", res(276,1)
!  write(*,*) " i = 276, res(2) =  ", res(276,2)
!  write(*,*) " i = 276, res(3) =  ", res(276,3)
!  write(*,*) " i = 276, res(4) =  ", res(276,4)
!  write(*,*) " i = 276, res(5) =  ", res(276,5)

 end subroutine compute_residual
!********************************************************************************

 end module residual
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

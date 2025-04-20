!********************************************************************************
!  Educationally-Designed Unstructured 3D (EDU3D) Code
!
!  --- EDU3D Euler
!
!* This module containes a subroutine that comptutes a set of grid data required
!* for the node-centered edge-based discretization for a tetrahedral grid.
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

 module nceb_data

  implicit none

 !This module contains the following subroutine:

  public :: construct_nceb_data  ! construct node-centered edge-based data
  public :: my_alloc_int_ptr     ! function for allocation
  public :: tet_volume           ! compute a volume of tetrahedron
  public :: triangle_area_vector ! compute an area vector for a triangular element

 contains

!*******************************************************************************
! Compute the data required for the edge-based scheme, and check the data.
!
! The following data will be constructed:
!
!
!   --             nghbrs_total = sum of # of nghbrs over all nodes
!   --        nnghbrs(1:nnodes) = # of nghbrs
!   --             k0(1:nnodes) = starting index of nghbrs for each node
!   --    nghbr(1:nghbrs_total) = stored nghbr numbers (to access nghbrs)
!
!   --                 nedges   = # of edges
!   --       edge(1:nedges,1:2) = edge list containing two end nodes
!   -- kth_nghbr_of_1(1:nedges) = node2 of an edge is k-th nghbr of node1
!   -- kth_nghbr_of_2(1:nedges) = node1 of an edge is k-th nghbr of node2
!
!   --            n12(1:nedges) = unit directed area vector at edge
!   --        n12_mag(1:nedges) = magnitude of directed area vector at edge
!
!   --            e12(1:nedges) = unit edge vector from node1 to node2 at edge
!   --        e12_mag(1:nedges) = magnitude of the edge vector
!
!   --            vol(1:nnodes) = dual volume at nodes
!   --          bmark(1:nnodes) = boundary number of boundary nodes
!
!   --          bn(1:ntria,1:3) = unit outward normal vector of boundary element
!   --      bn_mag(1:ntria    ) = magnitude of outward normal vector of boundary element
!
!   --     bnj(1:ntria,1:3,1:3) = unit outward normal vector at boundary node
!   -- bnj_mag(1:ntria,1:3    ) = magnitude of outward normal vector at boundary node
!
!
! Note: 'ntria' is the number of boundary triangular elements read from
!       an input grid file in 'read_grid' in the main program.
!
! Note: The normal at a boundary node, bnj is stored for each node in 
!       a boundary element, considering the possibility for defining multiple 
!       normals at the same boundary node. E.g., bnj(3,1,1), bnj(3,1,2), bnj(3,1,3)
!       are the x-, y-, z-components of the normal vector at node1 in boundary element 3.
!
! Node: Edge vector, ev and ev_mag, can be easily computed on the fly.
!       Then, no need to store them. It's up to you.
!
!
!
! References:
!
! [1] H. Nishikawa and Y. Liu, " Accuracy-Preserving Source Term Quadrature for
!     Third-Order Edge-Based Discretization", Journal of Computational Physics,
!     Volume 344, September 2017, Pages 595-622.
!     https://doi.org/10.1016/j.jcp.2017.04.075
!     http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_liu_jcp2017_preprint.pdf
!
! [2] H. Nishikawa, "Beyond Interface Gradient: A General Principle for
!     Constructing Diffusion Schemes", AIAA Paper 2010-5093, 40th AIAA Fluid
!     Dynamics Conference and Exhibit, Chicago, 2010.
!     http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_AIAA-2010-5093.pdf
!     (NOTE: See Appendix B for dual volumes and directed area vectors, boundary.)
!
!*******************************************************************************
 subroutine construct_nceb_data

!------------------------------------------------------------------------------------

 !Input
  use data_module    , only : nnodes                , & !nodal coords
                              x, y, z               , & !nodal coords
                              ntria, tria           , & !boundary elements
                              ntet, tet             , & !elements
                              ix, iy, iz            , & !ix=1, iy=2, iz=3
                              zero, half, third, fourth, sixth, p2, two

 !Output: data to be filled.
  use data_module    , only : vol                   , & !nodal coords, dual vol
                              nedges, edge, n12     , & !edge data
                              n12, n12_mag          , & !directed area vector
                              nghbrs_total          , & !Total dimension of nghbr array
                              k0, nnghbrs, nghbr    , & !node nghbr data
                              kth_nghbr_of_1        , & !nghbr index
                              kth_nghbr_of_2        , & !nghbr index
                              bn, bn_mag            , & !boundary normals
                              bnj, bnj_mag          , & !boundary normals at nodes
                              e12, e12_mag, bmark, bmarks  !Edge vector, boundary mark

  use gradients      , only : compute_lsq_coefficeints

!------------------------------------------------------------------------------------

  implicit none

 !Special data structure (to be used temporarily)
  type node_nghbr_type
   integer                         :: nnghbrs !number of neighbors
   integer, dimension(:), pointer  :: nghbr   !list of neighbors
  end type node_nghbr_type

! Temprary data
  type(node_nghbr_type), dimension(:), pointer  :: temp

! Dual volume at nodes
  real(p2), allocatable, dimension(:)           :: vol2

! Node neighbors of nodes.
  integer , allocatable, dimension(:,:)         :: edge_number

! Cell centroid coordinates.
  real(p2)                 :: xm , ym , zm
  real(p2)                 :: xc , yc , zc
  real(p2)                 :: xcl, ycl, zcl
  real(p2)                 :: xcr, ycr, zcr
  integer                  :: i_edge
  integer                  :: k1, k2
  integer                  :: n1, n2
  integer                  :: kv1, kv2, kv3
  integer                  :: v1, v2, v3

! Local pointers
  integer , dimension(4,3)   :: loc_nghbr
  integer , dimension(6,2)   :: loc_edge
  integer , dimension(6,2,3) :: loc_face
  integer                    :: left, rght ! Left and right face of an edge in tet.

  integer                    :: max_nnghbrs
  integer                    :: i, j, k, m, iv, itemp
  logical                    :: found


  real(p2), dimension(3)                :: edge_vector
  real(p2), dimension(3)                :: directed_area_contribution
  real(p2)                              :: dxjk_dot_njk

! Two dual face normals of edge in a tetrahedron.
  real(p2), dimension(3)                :: normal1, normal2

! Boundary normal
  real(p2), dimension(3)                :: normal_b

! For data check.
  real(p2), allocatable, dimension(:)   :: diff_dual
  real(p2)                              :: total_vol1
  real(p2)                              :: total_vol2
  real(p2), allocatable, dimension(:,:) :: sum_directed_area_vector

  real(p2)   :: x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4

  logical                               :: error_in_data
  integer                               :: nbmarks(4)

  real(p2), allocatable, dimension(:,:) :: bnj_temp
  real(p2), allocatable, dimension(:)   :: bnj_mag_temp

   write(*,*)
   write(*,*) "---------------------------------------"
   write(*,*) " Edge based data construction and check." 
   write(*,*)

   left = 1
   rght = 2

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Preparation (1):
!
! Define some local pointers within a tetrahedron ordered as below.
!
!            4
!            o
!           /| .
!          / |   .
!         /  |     .
!        /   |       .
!     3 o----|-------o 2
!        \   |     .
!         \  |    .
!          \ |  .
!           \|.
!            o
!            1
!

 !Neighbor nodes of a given node within a tetrahedron.

  loc_nghbr(1,1) = 2 !nghbr1 of vertex 1.
  loc_nghbr(1,2) = 3 !nghbr2 of vertex 1.
  loc_nghbr(1,3) = 4 !nghbr3 of vertex 1.

  loc_nghbr(2,1) = 1 !nghbr1 of vertex 2.
  loc_nghbr(2,2) = 4 !nghbr2 of vertex 2.
  loc_nghbr(2,3) = 3 !nghbr3 of vertex 2.

  loc_nghbr(3,1) = 1 !nghbr1 of vertex 3.
  loc_nghbr(3,2) = 2 !nghbr2 of vertex 3.
  loc_nghbr(3,3) = 4 !nghbr3 of vertex 3.

  loc_nghbr(4,1) = 1 !nghbr1 of vertex 4.
  loc_nghbr(4,2) = 3 !nghbr2 of vertex 4.
  loc_nghbr(4,3) = 2 !nghbr3 of vertex 4.

 !Local edges in a tetrahedron: pointing from node1 to node2.

  loc_edge(1,1) = 1 !node1 of local edge 1
  loc_edge(1,2) = 2 !node2 of local edge 1

  loc_edge(2,1) = 2 !node1 of local edge 2
  loc_edge(2,2) = 3 !node2 of local edge 2

  loc_edge(3,1) = 3 !node1 of local edge 3
  loc_edge(3,2) = 1 !node2 of local edge 3

  loc_edge(4,1) = 1 !node1 of local edge 4
  loc_edge(4,2) = 4 !node2 of local edge 4

  loc_edge(5,1) = 2 !node1 of local edge 5
  loc_edge(5,2) = 4 !node2 of local edge 5

  loc_edge(6,1) = 3 !node1 of local edge 6
  loc_edge(6,2) = 4 !node2 of local edge 6

 !Local adjacent faces of a local edge in a tetrahedron.
 ! - Left and right face seen from outside of the tetrahedron w.r.t. the edge orientation.
 ! - Nodes ordered in each face to point outward.

  loc_face(1,left,1) = 1 !left face of edge 1 (edge nodes: 1-2)
  loc_face(1,left,2) = 2 !left face of edge 1 (edge nodes: 1-2)
  loc_face(1,left,3) = 4 !left face of edge 1 (edge nodes: 1-2)

  loc_face(1,rght,1) = 1 !rght face of edge 1 (edge nodes: 1-2)
  loc_face(1,rght,2) = 3 !rght face of edge 1 (edge nodes: 1-2)
  loc_face(1,rght,3) = 2 !rght face of edge 1 (edge nodes: 1-2)


  loc_face(2,left,1) = 2 !left face of edge 2 (edge nodes: 2-3)
  loc_face(2,left,2) = 3 !left face of edge 2 (edge nodes: 2-3)
  loc_face(2,left,3) = 4 !left face of edge 2 (edge nodes: 2-3)

  loc_face(2,rght,1) = 2 !rght face of edge 2 (edge nodes: 2-3)
  loc_face(2,rght,2) = 1 !rght face of edge 2 (edge nodes: 2-3)
  loc_face(2,rght,3) = 3 !rght face of edge 2 (edge nodes: 2-3)


  loc_face(3,left,1) = 1 !left face of edge 3 (edge nodes: 3-1)
  loc_face(3,left,2) = 4 !left face of edge 3 (edge nodes: 3-1)
  loc_face(3,left,3) = 3 !left face of edge 3 (edge nodes: 3-1)

  loc_face(3,rght,1) = 1 !rght face of edge 3 (edge nodes: 3-1)
  loc_face(3,rght,2) = 3 !rght face of edge 3 (edge nodes: 3-1)
  loc_face(3,rght,3) = 2 !rght face of edge 3 (edge nodes: 3-1)


  loc_face(4,left,1) = 1 !left face of edge 4 (edge nodes: 1-4)
  loc_face(4,left,2) = 4 !left face of edge 4 (edge nodes: 1-4) 
  loc_face(4,left,3) = 3 !left face of edge 4 (edge nodes: 1-4)

  loc_face(4,rght,1) = 1 !rght face of edge 4 (edge nodes: 1-4)
  loc_face(4,rght,2) = 2 !rght face of edge 4 (edge nodes: 1-4)
  loc_face(4,rght,3) = 4 !rght face of edge 4 (edge nodes: 1-4)


  loc_face(5,left,1) = 1 !left face of edge 5 (edge nodes: 2-4)
  loc_face(5,left,2) = 2 !left face of edge 5 (edge nodes: 2-4)
  loc_face(5,left,3) = 4 !left face of edge 5 (edge nodes: 2-4)

  loc_face(5,rght,1) = 2 !rght face of edge 5 (edge nodes: 2-4)
  loc_face(5,rght,2) = 3 !rght face of edge 5 (edge nodes: 2-4)
  loc_face(5,rght,3) = 4 !rght face of edge 5 (edge nodes: 2-4)


  loc_face(6,left,1) = 2 !left face of edge 6 (edge nodes: 3-4)
  loc_face(6,left,2) = 3 !left face of edge 6 (edge nodes: 3-4)
  loc_face(6,left,3) = 4 !left face of edge 6 (edge nodes: 3-4)

  loc_face(6,rght,1) = 1 !rght face of edge 6 (edge nodes: 3-4)
  loc_face(6,rght,2) = 4 !rght face of edge 6 (edge nodes: 3-4)
  loc_face(6,rght,3) = 3 !rght face of edge 6 (edge nodes: 3-4)

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Preparation (2):
!
! Construct a list of neighbor nodes of each node.

  allocate( temp(nnodes) ) !<- will be deallocated later.

  do i = 1, nnodes
   temp(i)%nnghbrs = 0
   call my_alloc_int_ptr(temp(i)%nghbr, 1)
  end do

 !-------------------------------------------------
 ! Find edge-connected neghbr nodes within a tetrahedron.
  do i = 1, ntet

  !-------------------------------------------------
  ! Loop over the nodes of tetrahedron i.
   do k = 1, 4

     n1 = tet(i,k)

   !-------------------------------------------------
   ! Loop over the rest of nodes.
    do j = 1, 3
     n2 = tet(i,loc_nghbr(k,j))

     if ( temp(n1)%nnghbrs == 0 ) then

      temp(n1)%nnghbrs  = 1
      temp(n1)%nghbr(1) = n2

     else

     !Check if the neighbor is already added.
      found = .false.
      do m = 1, temp(n1)%nnghbrs
       if (n2 == temp(n1)%nghbr(m)) then
        found = .true.
        exit
       endif
      end do

     !If not added yet, add it to the list.
      if (.not.found) then
       temp(n1)%nnghbrs = temp(n1)%nnghbrs + 1
       call my_alloc_int_ptr(temp(n1)%nghbr, temp(n1)%nnghbrs)
       temp(n1)%nghbr( temp(n1)%nnghbrs ) = n2
      endif

     endif

    end do
   !-------------------------------------------------
   ! End of Loop over the rest of nodes.

   end do
  !-------------------------------------------------
  ! End of Loop over the nodes of tetrahedron i.

  end do
 ! End of Find edge-connected neghbr nodes from tets.
 !-------------------------------------------------


 !-------------------------------------------------
 ! Store the nghbr info into 1D array.
 !
 ! nghbr(1,2,3,4,....,) stores all nghbr information.
 ! Loop over nghbrs at a node i, therefore, can be
 ! performed as follows:
 !
 !  do k = 1, nnghbrs(i)
 !   nghbr_node = nghbr( k0(i) + k )
 !  end do
 !

  !Count the total number of nghbrs over all nodes.
    nghbrs_total = 0
   do i = 1, nnodes
    nghbrs_total = nghbrs_total + temp(i)%nnghbrs
   end do

  !Allocate the new arrays.
   allocate(nghbr(nghbrs_total)) ! store all nghbrs
   allocate(k0(         nnodes)) ! starting index
   allocate(nnghbrs(    nnodes)) ! number of nghbrs for each node

        itemp = 0
        k0(1) = 0

  !----------------------------------------------
  ! Loop over nodes
    do i = 1, nnodes

     nnghbrs(i) = temp(i)%nnghbrs

   !----------------------------------------------
   ! Loop over nghbr nodes by 'temp' and transfer the infor to nghbr(:).

    do k = 1, nnghbrs(i)
           itemp  = k0(i)+k
     nghbr(itemp) = temp(i)%nghbr(k)
    end do

   ! End of loop over nghbr nodes
   !----------------------------------------------

     !This is the last one for node i, and used as a starting index for node i+1.
       if (i < nnodes) k0(i+1) = itemp

    end do
  ! End of loop over nodes
  !----------------------------------------------

  !We don't use temp any more. Thanks.
   deallocate(temp)

  !From now on, we use 'k0', 'nnghbrs', and 'nghbr' to loop oever nghbr nodes.

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Preparation (3):
!
! Define the edges. Edge points from node1 to node2,
! where node2 > node1 (i.e., smaller node number to larger node number).
! So, the edges are oriented from a node number to a larger node number.

 !This temprary data will be used later to identify the edge number from two node numbers.
  max_nnghbrs = maxval(nnghbrs)
  allocate(  edge_number(nnodes, max_nnghbrs ) ) !<- will be deallocated later.
  edge_number = 0

 !---------------------------------------
 ! First, count the number of total edges.

  nedges = 0

  do i = 1, nnodes

   do k = 1, nnghbrs(i)

    if ( i < nghbr(k0(i)+k) ) then

      nedges = nedges + 1

    endif

   end do

  end do

 !---------------------------------------
 ! Allocate and fill the edge array.

  allocate(edge(nedges,2))

  nedges     = 0

  do i = 1, nnodes

   do k = 1, nnghbrs(i)

    if ( i < nghbr(k0(i)+k) ) then

              nedges = nedges + 1
      edge(nedges,1) = i
      edge(nedges,2) = nghbr(k0(i)+k)

      edge_number(i,k) = nedges

    endif

   end do

  end do

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Set kth_nghbr data (used for Jacobian off-diagonal block construction and LSQ)

  allocate( kth_nghbr_of_1(nedges) )
  allocate( kth_nghbr_of_2(nedges) )

  do j = 1, nedges
  
   n1 = edge(j,1)
   n2 = edge(j,2)

  !Loop over nghbrs of n2, and find n1.
   do k = 1, nnghbrs(n2)

    if ( n1 == nghbr(k0(n2)+k) ) then
     kth_nghbr_of_2(j) = k !Save k: n1 is kth nghbr of n2.
    endif

   end do

  !Loop over nghbrs of n2, and find n1.
   do k = 1, nnghbrs(n1)

    if ( n2 == nghbr(k0(n1)+k) ) then
     kth_nghbr_of_1(j) = k !Save k: n2 is kth nghbr of n1.
    endif

   end do

  end do

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
! Now, we can compute the directed area vector at each edge and the dual volume at
! each node.

  write(*,*) " >>> Compute directed area and dual volume..."

 !Dual volume is at nodes, and computed here in two different ways for checking.
  allocate(vol( nnodes)) !<- This is used in the solver.
  allocate(vol2(nnodes)) !<- will be deallocated later.

 !Directed area vector is at edge.
  allocate(n12(nedges,3))

 !Initialize them.
    vol  = zero
    vol2 = zero
    n12  = zero

 !Compute them via a loop over tetrahedra.

  tet_loop : do i = 1, ntet

   !Compute the tetra centroid coordinates.
     xc = fourth*( x(tet(i,1)) + x(tet(i,2)) + x(tet(i,3)) + x(tet(i,4)) )
     yc = fourth*( y(tet(i,1)) + y(tet(i,2)) + y(tet(i,3)) + y(tet(i,4)) )
     zc = fourth*( z(tet(i,1)) + z(tet(i,2)) + z(tet(i,3)) + z(tet(i,4)) )

   local_edge_loop : do j = 1, 6 ! Each tetra has 6 local edges.

    !Local edge orientation: k1 -> k2.
        k1 = loc_edge(j,1)        !  Local node number of node1 of edge j.
        k2 = loc_edge(j,2)        !  Local node number of node2 of edge j.

    !Global edge orientation not necessarily n1 -> n2. It is min(n1,n2)-> max(n1,n2).
        n1 = tet(i,k1)            ! Global node number of node1 of edge j.
        n2 = tet(i,k2)            ! Global node number of node2 of edge j.

   !Edge midpoint coordinates.
    xm = half*( x(n1) + x(n2) )
    ym = half*( y(n1) + y(n2) )
    zm = half*( z(n1) + z(n2) )
 
   !Centroid coordinates of the left face of the edge j.
    kv1 = loc_face(j,left,1)
    kv2 = loc_face(j,left,2)
    kv3 = loc_face(j,left,3)
     v1 = tet(i,kv1)
     v2 = tet(i,kv2)
     v3 = tet(i,kv3)
    xcl = third*( x(v1) + x(v2) + x(v3) )
    ycl = third*( y(v1) + y(v2) + y(v3) )
    zcl = third*( z(v1) + z(v2) + z(v3) )

   !Centroid coordinates of the right face of the edge j.
    kv1 = loc_face(j,rght,1)
    kv2 = loc_face(j,rght,2)
    kv3 = loc_face(j,rght,3)
     v1 = tet(i,kv1)
     v2 = tet(i,kv2)
     v3 = tet(i,kv3)
    xcr = third*( x(v1) + x(v2) + x(v3) )
    ycr = third*( y(v1) + y(v2) + y(v3) )
    zcr = third*( z(v1) + z(v2) + z(v3) )

   !-------------------------------------------------------------------------------
   ! Dual face normal vectors. Here, we compute them on the assumption that the edge
   ! orientation is n1 -> n2. The sum is the contributon to directed area vector.
   !-------------------------------------------------------------------------------

    ! Dual face 1: triangle with (edge midpoint)-(tetra centroid)-(left face centrpoid).

     normal1 = triangle_area_vector(xm,xc,xcl, ym,yc,ycl, zm,zc,zcl)

    ! Dual face 2: triangle with (edge midpoint)-(right face centrpoid)-(tetra centroid).

     normal2 = triangle_area_vector(xm,xcr,xc, ym,ycr,yc, zm,zcr,zc)

    ! Sum of the two.

     directed_area_contribution = normal1 + normal2

   !-------------------------------------------------------------------------------
   ! To compute the dual volume by
   !        dual volume = sum_over_edges [directed area vector]*[edge vector]/6 in 3D.
   ! Accumulate the contribution to the dual volume.
   !-------------------------------------------------------------------------------

   ! Compute the inner product: directed_area_vector * edge_vector.
   ! This is independent of the edge orientation. So, keep assuming n1 -> n2,
   ! for consistency with the directed area vector contribution.

      edge_vector(ix) = x(n2) - x(n1)
      edge_vector(iy) = y(n2) - y(n1)
      edge_vector(iz) = z(n2) - z(n1)

      dxjk_dot_njk = dot_product( directed_area_contribution, edge_vector )

      vol(n1) = vol(n1) + sixth * dxjk_dot_njk
      vol(n2) = vol(n2) + sixth * dxjk_dot_njk

   !-------------------------------------------------------------------------------
   ! To compute the directed-area vector of the edge.
   ! Add the contribution to the directed-area vector.
   !-------------------------------------------------------------------------------

      i_edge = 0

   ! Case(1): Edge orientations match: local k1 -> k2, global n1 -> n2.
     if (n1 < n2) then

     !Find the global edge number for n1 -> n2.
      do k = 1, nnghbrs(n1)
       if (nghbr(k0(n1)+k) == n2) i_edge =  edge_number(n1,k)
      end do
      if (i_edge == 0) stop 127 ! Error

     !Add the contribution to the directed area vector to 'i_edge'.
      n12(i_edge,:) = n12(i_edge,:) + directed_area_contribution

   ! Case(2): Edge orientations are opposite: local k1 -> k2, global n2 -> n1.
     else

     !Find the global edge number for n2 -> n1.
      do k = 1, nnghbrs(n2)
       if (nghbr(k0(n2)+k) == n1) i_edge =  edge_number(n2,k)
      end do
      if (i_edge == 0) stop 127 ! Error

     !Add the contribution to the directed area vector to 'i_edge' with minus sign.
      n12(i_edge,:) = n12(i_edge,:) - directed_area_contribution

     endif

   end do local_edge_loop

   !----------------------------------------------------------------------------------
   ! Alternative: Compute the dual volume as a sum of 1/4 of the tetrahedral volume.

    x1 = x(tet(i,1))
    x2 = x(tet(i,2))
    x3 = x(tet(i,3))
    x4 = x(tet(i,4))

    y1 = y(tet(i,1))
    y2 = y(tet(i,2))
    y3 = y(tet(i,3))
    y4 = y(tet(i,4))

    z1 = z(tet(i,1))
    z2 = z(tet(i,2))
    z3 = z(tet(i,3))
    z4 = z(tet(i,4))

    vol2(tet(i,:)) = vol2(tet(i,:)) + fourth*tet_volume(x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4)
   !----------------------------------------------------------------------------------

  end do tet_loop

  !-------------------------------------------------------------------------------
  ! Normalized directed area vector and edge vector.
  !-------------------------------------------------------------------------------

   allocate( n12_mag(nedges  ) )
   allocate( e12_mag(nedges  ) )
   allocate( e12(    nedges,3) )

   edge_loop : do j = 1, nedges

    n12_mag(j) = sqrt( n12(j,ix)**2 + n12(j,iy)**2 + n12(j,iz)**2 )
    n12(j,:)   = n12(j,:) / n12_mag(j)

            n1 = edge(j,1)
            n2 = edge(j,2)
     e12(j,ix) = x(n2) - x(n1)
     e12(j,iy) = y(n2) - y(n1)
     e12(j,iz) = z(n2) - z(n1)
    e12_mag(j) = sqrt( e12(j,ix)**2 + e12(j,iy)**2 + e12(j,iz)**2 )
    e12(j,:)   = e12(j,:) / e12_mag(j)

   end do edge_loop

   deallocate(edge_number)

  !-------------------------------------------------------------------------------
  ! Boundary marks:
  !
  !  bmark = 0 for interior nodes.
  !        = 1 for a node on a boundary.
  !-------------------------------------------------------------------------------

     allocate(bmark(nnodes))
     bmark = 0
 
    do i = 1, ntria
     bmark( tria(i,1:3) ) = max( bmark( tria(i,1:3) ),  1 )
    end do

  !-------------------------------------------------------------------------------
  ! Boundary marks:
  !
  !  bmarks = 1 for boundary nodes
  !         = 2 for interior nodes adjacent to a boundary node.
  !         = 3 for interior nodes adjacent to a node with bmarks=2.
  !         = 4 for interior nodes adjacent to a node with bmarks=3.
  !-------------------------------------------------------------------------------

  !(0)Initialize bmarks:

    allocate(bmarks(nnodes))
    bmarks = 0

    nbmarks = 0

  !(1)Assign bmarks = 1 for boundary nodes.

    do i = 1, nnodes
     if (bmark(i) == 1) then 
        bmarks(i) = 1 !=1 for all boundary nodes.
       nbmarks(1) = nbmarks(1) + 1 
     endif
    end do

  !(2)Find neighbors of boundary nodes and assign bmarks = 2.

    do i = 1, nnodes
    if (bmarks(i) == 1) then

    !Loop over neighbors of a boundary node i.
     do k = 1, nnghbrs(i)

      if ( bmarks(nghbr(k0(i)+k)) == 0 ) then
           bmarks(nghbr(k0(i)+k)) = 2
                       nbmarks(2) = nbmarks(2) + 1
      endif

     end do

    endif
    end do

  !(3)Find neighbors of bmarks=2 nodes and assign bmarks=3.

    do i = 1, nnodes
    if (bmarks(i) == 2) then

    !Loop over neighbors of a node i of bmarks=2.
     do k = 1, nnghbrs(i)

      if ( bmarks(nghbr(k0(i)+k)) == 0 ) then
           bmarks(nghbr(k0(i)+k)) = 3
                       nbmarks(3) = nbmarks(3) + 1 
      endif

     end do

    endif
    end do

  !(4)Find neighbors of bmarks=3 nodes and assign bmarks=4.

    do i = 1, nnodes
    if (bmarks(i) == 3) then

    !Loop over neighbors of a node i of bmarks=3.
     do k = 1, nnghbrs(i)

      if ( bmarks(nghbr(k0(i)+k)) == 0 ) then
           bmarks(nghbr(k0(i)+k)) = 4
                       nbmarks(4) = nbmarks(4) + 1 
      endif

     end do

    endif
    end do

    write(*,*)
    write(*,*) " # of    boundary nodes = ", nbmarks(1), &
               " remaining interior = ", nnodes-nbmarks(1)

    write(*,*) " # of nghbrs of bmarks1 = ", nbmarks(2), &
               " remaining interior = ", nnodes-nbmarks(1)-nbmarks(2)

    write(*,*) " # of nghbrs of bmarks2 = ", nbmarks(3), &
               " remaining interior = ", nnodes-nbmarks(1)-nbmarks(2)-nbmarks(3)

    write(*,*) " # of nghbrs of bmarks3 = ", nbmarks(4), &
               " remaining interior = ", nnodes-nbmarks(1)-nbmarks(2)-nbmarks(3)-nbmarks(4)
    write(*,*)

  !-------------------------------------------------------------------------------
  ! Compute and store the boundary element normals.
  !-------------------------------------------------------------------------------

  write(*,*) " >>> Compute boundary normals..."

  allocate( bn(ntria,3), bn_mag(ntria) )

  boundary_element : do i = 1, ntria

  !Three nodes of the triangle i.
   v1 = tria(i,1)
   v2 = tria(i,2)
   v3 = tria(i,3)

  !Boundary triangle is oriented inward. So, reverse the order to compute the outward vector.
   normal_b  = triangle_area_vector( x(v3),x(v2),x(v1), y(v3),y(v2),y(v1), z(v3),z(v2),z(v1) )

   bn_mag(i) = sqrt( normal_b(ix)**2 + normal_b(iy)**2 + normal_b(iz)**2 )
   bn(i,ix)  = normal_b(ix) / bn_mag(i)
   bn(i,iy)  = normal_b(iy) / bn_mag(i)
   bn(i,iz)  = normal_b(iz) / bn_mag(i)

  end do boundary_element

  !-------------------------------------------------------------------------------
  ! Compute and store the normals at boundary nodes.
  !-------------------------------------------------------------------------------


  !------------------------------------------------------------------
  ! Compute the lumped normals at boundary nodes (valid only up to second order).
  ! 
  ! We first compute the normals at boundary nodes in temporal arrays that
  ! runs over all nodes. That's too large (as normals are stored at boundary
  ! nodes only), but it is convenient for normal computations. At the end,
  ! we copy the boundary normals at nodes to a smaller array.

   allocate( bnj_temp(nnodes,3), bnj_mag_temp(nnodes) ) !Large arrays

   bnj_temp = zero !Initialize the array.

   boundary_element2 : do i = 1, ntria

    do k = 1, 3
 
    !Accumlate the boundary-element normals at nodes:

     bnj_temp(tria(i,k),:) = bnj_temp(tria(i,k),:) + bn(i,:)*bn_mag(i)

    end do

   end do boundary_element2

  !------------------------------------------------------------------
  ! Normalize them. Just to make them unit vectors.

   do i = 1, nnodes

     if ( bmark(i) /= 0) then !Only boundary nodes.

     !Magnitude of the normal at a boundary node iv.
      bnj_mag_temp(i) = sqrt( bnj_temp(i,ix)**2 + bnj_temp(i,iy)**2 + bnj_temp(i,iz)**2 )

     !Make it a unit vector
      bnj_temp(  i,:) = bnj_temp(i,:) / bnj_mag_temp(i)

     endif

   end do

  !------------------------------------------------------------------
  ! Now let's copy the data to smaller arrays, which can be accesed
  ! easily in a boundary element loop.

    allocate( bnj(ntria,3,3), bnj_mag(ntria,3) ) !Smaller arrays.

   boundary_element3 : do i = 1, ntria

    do k = 1, 3

                  iv = tria(i,k)

      bnj(    i,k,:) = bnj_temp(    iv,:)
      bnj_mag(i,k  ) = bnj_mag_temp(iv)

    end do

   end do boundary_element3

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  ! For debugging, visualize the normals at nodes as geometry
  ! lines in Tecplot. Load the boundary grid data and then
  ! this file, and you'll see the normal vectors at boundary nodes.

    open(unit=3, file="normals_at_nodes_tec_geometry_lines.dat", status="unknown")

    do i = 1, nnodes
     if (bmark(i)/=0) then !Only the boundary nodes.

      write(3,*) " GEOMETRY T=LINE3D LT=0.35 " !LT=thickness
      write(3,*) 1 ! 1 line
      write(3,*) 2 ! 2 points per line

     !First point is the boundary node location.
      write(3,*) x(i), y(i), z(i)

     !Second point is a small distance from the node location in the normal direction.
      write(3,*) x(i) + bnj_temp(i,ix)*0.1, &
                 y(i) + bnj_temp(i,iy)*0.1, &
                 z(i) + bnj_temp(i,iz)*0.1
     endif
    end do

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  !---------------------------------------------------------------

  !Deallocate temporary arrays.
    deallocate(bnj_temp, bnj_mag_temp)


! Edge-based data construction is complete. Now, we verify them below.

!------------------------------------------------------------------------------------
! Verify the data.
!------------------------------------------------------------------------------------

  write(*,*) " >>> Verify the dual volumes and directed area vectors..."

! 1. Check the dual volume sum, which must be equal to the total volume.
!    Also, check if the dual volumes computed in diffrent ways match.

   error_in_data = .false.

   allocate(diff_dual(nnodes)) !<- Temporary array. It will be deallocated later.

   total_vol1 = zero
   total_vol2 = zero

   do i = 1, nnodes
       total_vol1 = total_vol1 + vol(i)
       total_vol2 = total_vol2 + vol2(i)
     diff_dual(i) = ( vol(i) - vol2(i) ) / vol(i) ! Relative error
   end do

   if ( abs(total_vol1-total_vol2)/total_vol1 > 1.0e-12_p2 ) error_in_data = .true.

   if ( error_in_data ) then
    write(*,*)
    write(*,*) " Total sum of vol1 and that of vol2 must be the same: "
    write(*,*) "           sum of the dual_vol1 = ", total_vol1
    write(*,*) "           sum of the dual_vol2 = ", total_vol2
    write(*,*)
    write(*,*) " vol1 and vol2 must be the same at all nodes: see relative differences"
    write(*,*) " Max of sum_at_node (dual_vol1-dual_vol2)/total_vol1 ", maxval(diff_dual)/total_vol1
    write(*,*) " Min of sum_at_node (dual_vol1-dual_vol2)/total_vol1 ", minval(diff_dual)/total_vol1
    write(*,*)
    write(*,*) " These must be zero. Are these small enough? ERROR in construct_nceb_data? Stop..."
    write(*,*)
    stop 127
   endif

   deallocate(diff_dual)
   deallocate(vol2)

! 2. Check the directed area: Sum must be zero at all nodes.

 !Accumulate the directed area at nodes: all outward from the node.

   error_in_data = .false.

  allocate(sum_directed_area_vector(nnodes,3)) !<- Temporary array. It will be deallocated later.

   sum_directed_area_vector = zero

  edge_loop2 : do j = 1, nedges

   n1 = edge(j,1)
   n2 = edge(j,2)

  !directed_area_vector points from n1 to n2. So, add to n1.
   sum_directed_area_vector(n1,:) =  sum_directed_area_vector(n1,:) + n12(j,:)*n12_mag(j)

  !directed_area_vector points from n1 to n2. So, add to n2 with minus sign.
   sum_directed_area_vector(n2,:) =  sum_directed_area_vector(n2,:) - n12(j,:)*n12_mag(j)

  end do edge_loop2

 !Close them with boundary outward normals.

  boundary_element4 : do i = 1, ntria

  !Three nodes of the triangle i.
   v1 = tria(i,1)
   v2 = tria(i,2)
   v3 = tria(i,3)

  !Dual area for each node is exactly 1/3 of the triangle area.
   sum_directed_area_vector(v1,:) = sum_directed_area_vector(v1,:) + third*bn(i,:)*bn_mag(i)
   sum_directed_area_vector(v2,:) = sum_directed_area_vector(v2,:) + third*bn(i,:)*bn_mag(i)
   sum_directed_area_vector(v3,:) = sum_directed_area_vector(v3,:) + third*bn(i,:)*bn_mag(i)

  end do boundary_element4

 ! Check the max magnitude over all nodes in the grid, which must be zero.

   if ( maxval(abs(sum_directed_area_vector(:,ix)/vol(:)**(two*third))) > 1.0e-12_p2 ) error_in_data = .true.
   if ( maxval(abs(sum_directed_area_vector(:,iy)/vol(:)**(two*third))) > 1.0e-12_p2 ) error_in_data = .true.
   if ( maxval(abs(sum_directed_area_vector(:,iz)/vol(:)**(two*third))) > 1.0e-12_p2 ) error_in_data = .true.

   if ( error_in_data ) then
    write(*,*)
    write(*,*) " Max of sum_at_node(directed area vector(x)) = "                ,  &
                 maxval(abs(sum_directed_area_vector(:,ix)/vol(:)**(two*third)))
    write(*,*) " Max of sum_at_node(directed area vector(y)) = "                ,  &
                 maxval(abs(sum_directed_area_vector(:,iy)/vol(:)**(two*third)))
    write(*,*) " Max of sum_at_node(directed area vector(z)) = "                ,  &
                 maxval(abs(sum_directed_area_vector(:,iz)/vol(:)**(two*third)))
    write(*,*)
    write(*,*) " These must be zero. Are these small enough? ERROR in construct_nceb_data? Stop..."
    write(*,*)
    stop 127
   endif

  deallocate(sum_directed_area_vector)

 !Compute LSQ coefficients
  call compute_lsq_coefficeints

   write(*,*) " >>> Verified."

 end subroutine construct_nceb_data
!*******************************************************************************


!********************************************************************************
!* This subroutine is useful to expand or shrink integer arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
!*  Array, x, will be shrinked to the requested dimension, n, if (n < dim(x)).
!*
!********************************************************************************
  subroutine my_alloc_int_ptr(x,n)

  implicit none
  integer, intent(in) :: n
  integer, dimension(:), pointer :: x
  integer, dimension(:), pointer :: temp
  integer :: i

  if (n <= 0) then
   write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
   stop 127
  endif

! If not allocated, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n))
   temp = 0

! (1) Expand the array dimension
  if ( n > size(x) ) then

   do i = 1, size(x)
    temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
  else

   do i = 1, n
    temp(i) = x(i)
   end do

  endif

! Destroy the target of x
!  deallocate(x)

! Re-assign the pointer
   x => temp

  return

  end subroutine my_alloc_int_ptr
!********************************************************************************

!*******************************************************************************
! Compute the volume of a tetrahedron defined by 4 vertices:
!
!       (x1,y1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4),
!
! which are ordered as follows:
!
!            1
!            o
!           /| .
!          / |   .
!         /  |     .
!        /   |       .
!     2 o----|-------o 3
!        \   |     .
!         \  |    .
!          \ |  .
!           \|.
!            o
!            4
!
! Note: Volume = volume integral of 1 = 1/3 * volume integral of div(x,y,z) dV
!              = surface integral of (x,y,z)*dS
!              = sum of [ (xc,yc,zc)*area_vector ] over triangular faces.
!
! where the last step is exact because (x,y,z) vary linearly over the triangle.
! There are other ways to compute the volume, of course.
!
!*******************************************************************************
 function tet_volume(x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4)

 implicit none

 integer , parameter :: dp = selected_real_kind(15) !Double precision

 !Input
 real(dp), intent(in)   :: x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4
 !Output
 real(dp)               :: tet_volume

 real(dp)               :: xc, yc, zc
 real(dp), dimension(3) :: area
 integer                :: ix=1, iy=2, iz=3


 tet_volume = 0.0_dp

! Triangle 1-3-2

   !Centroid of the triangular face
      xc = (x1+x3+x2)/3.0_dp
      yc = (y1+y3+y2)/3.0_dp
      zc = (z1+z3+z2)/3.0_dp
   !Outward normal surface vector
   area = triangle_area_vector(x1,x3,x2, y1,y3,y2, z1,z3,z2)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

! Triangle 1-4-3

   !Centroid of the triangular face
      xc = (x1+x4+x3)/3.0_dp
      yc = (y1+y4+y3)/3.0_dp
      zc = (z1+z4+z3)/3.0_dp

   !Outward normal surface vector
   area = triangle_area_vector(x1,x4,x3, y1,y4,y3, z1,z4,z3)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

! Triangle 1-2-4

   !Centroid of the triangular face
      xc = (x1+x2+x4)/3.0_dp
      yc = (y1+y2+y4)/3.0_dp
      zc = (z1+z2+z4)/3.0_dp

   !Outward normal surface vector
   area = triangle_area_vector(x1,x2,x4, y1,y2,y4, z1,z2,z4)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

! Triangle 2-3-4

   !Centroid of the triangular face
      xc = (x2+x3+x4)/3.0_dp
      yc = (y2+y3+y4)/3.0_dp
      zc = (z2+z3+z4)/3.0_dp

   !Outward normal surface vector
   area = triangle_area_vector(x2,x3,x4, y2,y3,y4, z2,z3,z4)

   tet_volume = tet_volume + ( xc*area(ix) + yc*area(iy) + zc*area(iz) )

   tet_volume = tet_volume / 3.0_dp

 end function tet_volume

!*******************************************************************************
! Compute the area of a triangle in 3D defined by 3 vertices:
!
!       (x1,y1,z1), (x2,y2,z2), (x3,y3,z3),
!
! which is assumed to be ordered clockwise.
!
!     1             2
!      o------------o
!       \         .
!        \       . --------->
!         \    .
!          \ .
!           o
!           3
!
! Note: Area is a vector based on the right-hand rule: 
!       when wrapping the right hand around the triangle with the fingers in the
!       direction of the vertices [1,2,3], the thumb points in the positive
!       direction of the area.
!
! Note: Area vector is computed as the cross product of edge vectors [31] and [32].
!
!*******************************************************************************
 function triangle_area_vector(x1,x2,x3, y1,y2,y3, z1,z2,z3) result(area_vector)
 
 implicit none
 integer , parameter :: dp = selected_real_kind(15) !Double precision

 !Input
  real(dp), intent(in)   :: x1,x2,x3, y1,y2,y3, z1,z2,z3
 !Output
  real(dp), dimension(3) :: area_vector

  integer :: ix=1, iy=2, iz=3

  !x-component of the area vector
   area_vector(ix) = 0.5_dp*( (y1-y3)*(z2-z3)-(z1-z3)*(y2-y3) )

  !y-component of the area vector
   area_vector(iy) = 0.5_dp*( (z1-z3)*(x2-x3)-(x1-x3)*(z2-z3) )

  !z-component of the area vector
   area_vector(iz) = 0.5_dp*( (x1-x3)*(y2-y3)-(y1-y3)*(x2-x3) )

 end function triangle_area_vector


 end module nceb_data


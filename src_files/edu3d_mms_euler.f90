!********************************************************************************
!* Educationally-Designed Unstructured 2D (EDU2D) Code
!*
!*  --- EDU2D_MMS_TWODNS
!*
!*
!* This program shows how to compute the source terms arising from the method
!* of manufactured solutions for the 2D compressible Navier-Stokes equations.
!*
!*
!* Method of manufactured solutions:
!*
!* (0) Consider a target equation, df(u)/dx + dg(u)/dy = 0.
!* (1) Choose a function, v(x,y).
!* (2) Substitute v into the tartget equation: s = df(v)/dx + dg(v)/dy.
!* (3) Then, v(x,y) is the exact solution for the following eq.:
!*
!*      df(u)/dx + dg(u)/dy = s.
!*
!* (4) Numerically solve this equation for u, then you can
!*     measure the error: error = |u - v|.
!*
!* This program shows how to evaluate the source term 's' numerically
!* for the 2D Navier-Stokes equation, for a chosen function v(x,y),
!* which is given here by
!*
!*     v(x,y) = c0 + cs*sin(cx*x+cy*y),
!*
!* where the coefficients, c0, cs, cx, cy, can be chosen differently
!* for different variables. You can change this function as you like.
!*
!*
!* Note: The main program calls the subroutine:
!*
!*        'compute_manufactured_sol_and_source'
!*
!*       for a specified (x,y), and prints on screen the solutions and
!*       source terms.       
!*
!* Note: You can use the subroutine 'compute_manufactured_sol_and_source'
!*       to compute the manufactured solution and the source term in
!*       your code.
!*
!* Note: The source (or forcing) term, 's', is a function of space only, not
!*       of the solution. Therefore, it is computed only once and stored.
!*
!* Note: It is quite straightforward to extend the subroutine to generate
!*       source terms for the 3D Navier-Stokes equations.
!*
!* Note: The viscosity is computed by Sutherland's law with a suitable
!*       scaling for the Navier-Stokes equations nondimensionalized with
!*       the freestream speed of sound: (u,v)=(u/a_inf,v/a_inf),
!*       p=p/(rho_inf*a_inf^2). Modify the factor if other (or no) 
!*       non-dimensionalization is used. See "I do like CFD, VOL.1" for
!*       details, which is available in PDF at http://www.cfdbooks.com.
!*
!*
!* This is Version 1 (January 23, 2021).
!*
!* This F90 code is written and made available for an educational purpose.
!* This file may be updated in future.
!*
!* Katate Masatsuka, http://www.cfdbooks.com
!********************************************************************************

 module edu3d_mms_euler

 
  implicit none

 !This module contains the following subroutine:

  public :: compute_manufactured_sol_and_source

 contains

!********************************************************************************
!*
!* Educationally-Designed Unstructured 3D (EDU3D) Code
!*
!*  --- EDU3D_MMS_Euler
!*
!*
!* This subroutine computes, for a given position (x,y,z), the manufactured
!* solution values and the source terms. The manufactured solution takes
!* the following form:
!*
!*    a0 + as*exp(sx*x + sy*y + sz*z)
!*
!* with the coefficients, a0, as, sx, sy may be chosen differently for different
!* variable.
!*
!* ------------------------
!* Input:
!*
!*   (x,y,z) = Position at which the solution and source terms are sought.
!*
!* ------------------------
!* Output:
!*         w(1:11) =  [rho,u,v,w,p]
!*         u(1:11) =  [rho*u,rho*u,rho*v,rho*w,p]
!*
!*
!*           s(1:5) = Source terms for continuity, x/y/z-momentum, and energy eqns.
!*          sx(1:5) = x-derivative of s.
!*          sy(1:5) = y-derivative of s.
!*          sz(1:5) = z-derivative of s.
!*
!* ------------------------
!*
!* Note: You can extend this to the Navier-Stokes equations by adding the viscous terms.
!*       See EDU2D_MMS.
!*
!* Note: You can extend this also to the unsteady equations by including
!*       time derivatives when generating the source term: s=dv/dt+df(v)/dx+dg(v)/dy.
!*
!* Note: For high-order schemes, you can use this subroutine to compute the
!*       source terms at quadrature points.
!*
!*
!* This is Version 1 (January 23, 2021).
!*
!* This F90 code is written and made available for an educational purpose.
!* This file may be updated in future.
!*
!* Katate Masatsuka, http://www.cfdbooks.com
!********************************************************************************

 subroutine compute_manufactured_sol_and_source( x,y,z, up,uc,s,sx,sy,sz )

 implicit none
 integer , parameter :: p2 = selected_real_kind(p=15)

 real(p2), parameter ::       half = 0.5_p2
 real(p2), parameter ::        one = 1.0_p2
 real(p2), parameter ::      gamma = 1.4_p2

!Input 
 real(p2),               intent( in) :: x, y, z

!Output
 real(p2), dimension(5), intent(out) :: up  !Primitive variables    (rho,u,v,w,p)
 real(p2), dimension(5), intent(out) :: uc !Conservative variables (rho,rho*u,rho*v,rho*w,rho*E)
 real(p2), dimension(5), intent(out) :: s, sx, sy, sz

!Local variables



 real(p2) ::  rho, rhox, rhoy,rhoz, rhoxx,rhoyy,rhozz,rhoxy,rhoyz,rhozx
 real(p2) ::  u, ux, uy, uz, uxx,uyy,uzz,uxy,uyz,uzx
 real(p2) ::  v, vx, vy, vz, vxx,vyy,vzz,vxy,vyz,vzx
 real(p2) ::  w, wx, wy, wz, wxx,wyy,wzz,wxy,wyz,wzx
 real(p2) ::  p, px, py, pz, pxx,pyy,pzz,pxy,pyz,pzx

 real(p2) ::  a, ax, ay, az, axx,ayy,azz,axy,ayz,azx
 real(p2) ::  b, bx, by, bz, bxx,byy,bzz,bxy,byz,bzx
 real(p2) ::  c, cx, cy, cz, cxx,cyy,czz,cxy,cyz,czx
 real(p2) ::  k, kx, ky, kz, kxx,kyy,kzz,kxy,kyz,kzx

 real(p2) :: u2, u2x, u2y, u2z, u2xx,u2yy,u2zz,u2xy,u2yz,u2zx
 real(p2) :: v2, v2x, v2y, v2z, v2xx,v2yy,v2zz,v2xy,v2yz,v2zx
 real(p2) :: w2, w2x, w2y, w2z, w2xx,w2yy,w2zz,w2xy,w2yz,w2zx

 real(p2) :: au, aux, auy, auz, auxx,auyy,auzz,auxy,auyz,auzx
 real(p2) :: av, avx, avy, avz, avxx,avyy,avzz,avxy,avyz,avzx
 real(p2) :: aw, awx, awy, awz, awxx,awyy,awzz,awxy,awyz,awzx

 real(p2) :: bu, bux, buy, buz, buxx,buyy,buzz,buxy,buyz,buzx
 real(p2) :: bv, bvx, bvy, bvz, bvxx,bvyy,bvzz,bvxy,bvyz,bvzx
 real(p2) :: bw, bwx, bwy, bwz, bwxx,bwyy,bwzz,bwxy,bwyz,bwzx

 real(p2) :: cu, cux, cuy, cuz, cuxx,cuyy,cuzz,cuxy,cuyz,cuzx
 real(p2) :: cv, cvx, cvy, cvz, cvxx,cvyy,cvzz,cvxy,cvyz,cvzx
 real(p2) :: cw, cwx, cwy, cwz, cwxx,cwyy,cwzz,cwxy,cwyz,cwzx

 real(p2) :: rhoH , rhoHx , rhoHy ,  rhoHz, rhoHxx ,rhoHyy ,rhoHzz ,rhoHxy ,rhoHyz ,rhoHzx
 real(p2) :: rhouH, rhouHx, rhouHy, rhouHz, rhouHxx,rhouHyy,rhouHzz,rhouHxy,rhouHyz,rhouHzx
 real(p2) :: rhovH, rhovHx, rhovHy, rhovHz, rhovHxx,rhovHyy,rhovHzz,rhovHxy,rhovHyz,rhovHzx
 real(p2) :: rhowH, rhowHx, rhowHy, rhowHz, rhowHxx,rhowHyy,rhowHzz,rhowHxy,rhowHyz,rhowHzx

 integer  :: ix, iy, iz
 integer  :: irho, iu, iv, iw, ip
 integer  :: iconti, ixmom, iymom, izmom, ienrgy

 real(p2) :: cr_a0, cr_as, cr_sx, cr_sy, cr_sz
 real(p2) :: cu_a0, cu_as, cu_sx, cu_sy, cu_sz
 real(p2) :: cv_a0, cv_as, cv_sx, cv_sy, cv_sz
 real(p2) :: cw_a0, cw_as, cw_sx, cw_sy, cw_sz
 real(p2) :: cp_a0, cp_as, cp_sx, cp_sy, cp_sz

!-----------------------------------------------------------
! Define some indices

        ix =  1 ! x component
        iy =  2 ! y component
        iz =  3 ! z component

      irho =  1 ! density
      iu   =  2 ! x-velocity
      iv   =  3 ! y-velocity
      iw   =  4 ! w-velocity
      ip   =  5 ! pressure

    iconti =  1 ! continuity equation
    ixmom  =  2 ! x-momentum equation
    iymom  =  3 ! y-momentum equation
    izmom  =  4 ! z-momentum equation
    ienrgy =  5 !     energy equation

!-----------------------------------------------------------
! Constants for the exact solution: a0 + as*sin(sx*x+sy*y+sz*z).
!
! Note: Make sure the density and pressure are positive.

 !-----------------------------------------
 ! Density    = cr_a0 + cr_as*exp(cr_sz*x+cr_sy*y+cr_sz*z)

  cr_a0 =  0.5_p2
  cr_as =  1.0_p2
  cr_sx =  0.525_p2
  cr_sy =  0.550_p2
  cr_sz =  0.575_p2

 !-----------------------------------------
 ! X-velocity = cu_a0 + cu_as*exp(cu_sz*x+cu_sy*y+cu_sz*z)

  cu_a0 =  0.1_p2
  cu_as =  0.1_p2
  cu_sx =  0.125_p2
  cu_sy =  0.150_p2
  cu_sz =  0.175_p2

 !-----------------------------------------
 ! Y-velocity = cv_a0 + cv_as*exp(cv_sz*x+cv_sy*y+cv_sz*z)

  cv_a0 =  0.2_p2
  cv_as =  0.2_p2
  cv_sx =  0.225_p2
  cv_sy =  0.250_p2
  cv_sz =  0.275_p2

 !-----------------------------------------
 ! Z-velocity = cw_a0 + cw_as*exp(cw_sz*x+cw_sy*y+cw_sz*z)

  cw_a0 =  0.3_p2
  cw_as =  0.3_p2
  cw_sx =  0.325_p2
  cw_sy =  0.350_p2
  cw_sz =  0.375_p2
  
 !-----------------------------------------
 ! Pressure = cp_a0 + cp_as*exp(cp_sz*x+cp_sy*y+cp_sz*z)

  cp_a0 =  0.4_p2
  cp_as =  1.0_p2/1.4_p2
  cp_sx =  0.425_p2
  cp_sy =  0.450_p2
  cp_sz =  0.475_p2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Part I: Compute w = [rho,u,v,w,p] and derivatives.
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

 !------------------------------------------------------------------------
 ! rho: Density and its derivatives
 
              rho = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 0,0,0,x,y,z)
             rhox = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 1,0,0,x,y,z)
             rhoy = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 0,1,0,x,y,z)
             rhoz = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 0,0,1,x,y,z)

            rhoxx = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 2,0,0,x,y,z)
            rhoyy = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 0,2,0,x,y,z)
            rhozz = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 0,0,2,x,y,z)

            rhoxy = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 1,1,0,x,y,z)
            rhoyz = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 0,1,1,x,y,z)
            rhozx = manufactured_sol(cr_a0,cr_as,cr_sx,cr_sy,cr_sz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! u: x-velocity and its derivatives

              u = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 0,0,0,x,y,z)
             ux = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 1,0,0,x,y,z)
             uy = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 0,1,0,x,y,z)
             uz = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 0,0,1,x,y,z)

            uxx = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 2,0,0,x,y,z)
            uyy = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 0,2,0,x,y,z)
            uzz = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 0,0,2,x,y,z)

            uxy = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 1,1,0,x,y,z)
            uyz = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 0,1,1,x,y,z)
            uzx = manufactured_sol(cu_a0,cu_as,cu_sx,cu_sy,cu_sz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! v: y-velocity and its derivatives

              v = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 0,0,0,x,y,z)
             vx = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 1,0,0,x,y,z)
             vy = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 0,1,0,x,y,z)
             vz = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 0,0,1,x,y,z)

            vxx = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 2,0,0,x,y,z)
            vyy = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 0,2,0,x,y,z)
            vzz = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 0,0,2,x,y,z)

            vxy = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 1,1,0,x,y,z)
            vyz = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 0,1,1,x,y,z)
            vzx = manufactured_sol(cv_a0,cv_as,cv_sx,cv_sy,cv_sz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! w: z-velocity and its derivatives

              w = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 0,0,0,x,y,z)
             wx = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 1,0,0,x,y,z)
             wy = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 0,1,0,x,y,z)
             wz = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 0,0,1,x,y,z)

            wxx = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 2,0,0,x,y,z)
            wyy = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 0,2,0,x,y,z)
            wzz = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 0,0,2,x,y,z)

            wxy = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 1,1,0,x,y,z)
            wyz = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 0,1,1,x,y,z)
            wzx = manufactured_sol(cw_a0,cw_as,cw_sx,cw_sy,cw_sz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! p: pressure and its derivatives

              p = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 0,0,0,x,y,z)
             px = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 1,0,0,x,y,z)
             py = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 0,1,0,x,y,z)
             pz = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 0,0,1,x,y,z)

            pxx = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 2,0,0,x,y,z)
            pyy = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 0,2,0,x,y,z)
            pzz = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 0,0,2,x,y,z)

            pxy = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 1,1,0,x,y,z)
            pyz = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 0,1,1,x,y,z)
            pzx = manufactured_sol(cp_a0,cp_as,cp_sx,cp_sy,cp_sz, 1,0,1,x,y,z)

 !-----------------------------------------------------------------------------
 ! Store the exact solution for return.

  !Primitive variables

      up(irho) = rho 
      up(iu)   = u
      up(iv)   = v
      up(iw)   = w
      up(ip)   = p

  !Conservative variables

         uc(1) = rho
         uc(2) = rho*u
         uc(3) = rho*v
         uc(4) = rho*w
         uc(5) = p/(gamma-one) + half*rho*(u*u + v*v + w*w) ! =rho*E

!-----------------------------------------------------------------------------
!
! Exact (manufactured) solutons and derivatives have been computed.
! We move onto the source terms, which are the governing equations we wish
! to solve (the compressible Euler) evaluated with the exact solution (i.e., 
! the function we wish to make exact with the source terms).
!
! Navier-Stokes: dF(w)/dx + dG(w)/dy + dH(w)/dz = 0.
!
! Our manufactured solutions, wm, are the exact soltuions to the following:
!
!   dF(w)/dx + dG(w)/dy + dH(w)/dz = S,
!
! where S = dF(wm)/dx + dG(wm)/dy + dH(wm)/dz.
!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Part II: Compute the source terms.
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!  Inviscid terms
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

! Derivatives of u^2
  call derivatives_ab2( u, ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uzx,  &
                        u, ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uzx,  & 
                       u2,u2x,u2y,u2z,u2xx,u2yy,u2zz,u2xy,u2yz,u2zx   ) ! a=u^2

! Derivatives of v^2
  call derivatives_ab2( v, vx, vy, vz, vxx, vyy, vzz, vxy, vyz, vzx,  &
                        v, vx, vy, vz, vxx, vyy, vzz, vxy, vyz, vzx,  & 
                       v2,v2x,v2y,v2z,v2xx,v2yy,v2zz,v2xy,v2yz,v2zx   ) ! b=v^2

! Derivatives of w^2
  call derivatives_ab2( w, wx, wy, wz, wxx, wyy, wzz, wxy, wyz, wzx,  &
                        w, wx, wy, wz, wxx, wyy, wzz, wxy, wyz, wzx,  & 
                       w2,w2x,w2y,w2z,w2xx,w2yy,w2zz,w2xy,w2yz,w2zx   ) ! c=w^2

! Derivatives of k=(u^2+v^2+w^2)/2

  k     = half*(u*u   + v*v   + w*w)
  kx    = half*(u2x   + v2x   + w2x)
  ky    = half*(u2y   + v2y   + w2y)
  kz    = half*(u2z   + v2z   + w2z)

  kxx   = half*(u2xx   + v2xx   + w2xx)
  kyy   = half*(u2yy   + v2yy   + w2yy)
  kzz   = half*(u2zz   + v2zz   + w2zz)

  kxy   = half*(u2xy   + v2xy   + w2xy)
  kyz   = half*(u2yz   + v2yz   + w2yz)
  kzx   = half*(u2zx   + v2zx   + w2zx)

! Derivatives of rho*k = rho*(u^2+v^2)/2
  call derivatives_ab2( rho, rhox, rhoy,rhoz, rhoxx,rhoyy,rhozz,rhoxy,rhoyz,rhozx, &
                          k,   kx,   ky,  kz,   kxx,  kyy,  kzz,  kxy,  kyz,  kzx, &
                          a,   ax,   ay,  az,   axx,  ayy,  azz,  axy,  ayz,  azx  ) !a=rho*(u^2+v^2+w^2)/2

! Derivatives of rho*H = gamma/(gamma-1)*p + rho*k

  rhoH     = gamma/(gamma-one)*p     + a
  rhoHx    = gamma/(gamma-one)*px    + ax
  rhoHy    = gamma/(gamma-one)*py    + ay
  rhoHz    = gamma/(gamma-one)*pz    + az

  rhoHxx   = gamma/(gamma-one)*pxx   + axx
  rhoHyy   = gamma/(gamma-one)*pyy   + ayy
  rhoHzz   = gamma/(gamma-one)*pzz   + azz
  rhoHxy   = gamma/(gamma-one)*pxy   + axy
  rhoHyz   = gamma/(gamma-one)*pyz   + ayz
  rhoHzx   = gamma/(gamma-one)*pzx   + azx

!-----------------------------------------------------------------------------

! Compute derivatives of (rho*u)
  call derivatives_ab2( rho, rhox, rhoy,rhoz, rhoxx,rhoyy,rhozz,rhoxy,rhoyz,rhozx, &
                          u,   ux,   uy,  uz,   uxx,  uyy,  uzz,  uxy,  uyz,  uzx, &
                          a,   ax,   ay,  az,   axx,  ayy,  azz,  axy,  ayz,  azx  ) ! a=(rho*u)

! Compute derivatives of (rho*v)
  call derivatives_ab2( rho, rhox, rhoy,rhoz, rhoxx,rhoyy,rhozz,rhoxy,rhoyz,rhozx, &
                          v,   vx,   vy,  vz,   vxx,  vyy,  vzz,  vxy,  vyz,  vzx, &
                          b,   bx,   by,  bz,   bxx,  byy,  bzz,  bxy,  byz,  bzx  ) ! b=(rho*v)

! Compute derivatives of (rho*w)
  call derivatives_ab2( rho, rhox, rhoy,rhoz, rhoxx,rhoyy,rhozz,rhoxy,rhoyz,rhozx, &
                          w,   wx,   wy,  wz,   wxx,  wyy,  wzz,  wxy,  wyz,  wzx, &
                          c,   cx,   cy,  cz,   cxx,  cyy,  czz,  cxy,  cyz,  czx  ) ! c=(rho*w)

!-----------------------------------------------------------------------------

! Compute derivatives of (a*u)=(rho*u*u)

  call derivatives_ab2( a,   ax,   ay,  az,   axx,  ayy,  azz,  axy,  ayz,  azx, &
                        u,   ux,   uy,  uz,   uxx,  uyy,  uzz,  uxy,  uyz,  uzx, &
                       au,  aux,  auy, auz,  auxx, auyy, auzz, auxy, auyz, auzx  )

! Compute derivatives of (a*v)=(rho*u*v)

  call derivatives_ab2( a,   ax,   ay,  az,   axx,  ayy,  azz,  axy,  ayz,  azx, &
                        v,   vx,   vy,  vz,   vxx,  vyy,  vzz,  vxy,  vyz,  vzx, &
                       av,  avx,  avy, avz,  avxx, avyy, avzz, avxy, avyz, avzx  )

! Compute derivatives of (a*v)=(rho*u*w)

  call derivatives_ab2( a,   ax,   ay,  az,   axx,  ayy,  azz,  axy,  ayz,  azx, &
                        w,   wx,   wy,  wz,   wxx,  wyy,  wzz,  wxy,  wyz,  wzx, &
                       aw,  awx,  awy, awz,  awxx, awyy, awzz, awxy, awyz, awzx  )

!---

! Compute derivatives of (b*u)=(rho*v*u)


  call derivatives_ab2( b,   bx,   by,  bz,   bxx,  byy,  bzz,  bxy,  byz,  bzx, &
                        u,   ux,   uy,  uz,   uxx,  uyy,  uzz,  uxy,  uyz,  uzx, &
                       bu,  bux,  buy, buz,  buxx, buyy, buzz, buxy, buyz, buzx  )

! Compute derivatives of (b*v)=(rho*v*v)

  call derivatives_ab2( b,   bx,   by,  bz,   bxx,  byy,  bzz,  bxy,  byz,  bzx, &
                        v,   vx,   vy,  vz,   vxx,  vyy,  vzz,  vxy,  vyz,  vzx, &
                       bv,  bvx,  bvy, bvz,  bvxx, bvyy, bvzz, bvxy, bvyz, bvzx  )

! Compute derivatives of (b*w)=(rho*v*w)

  call derivatives_ab2( b,   bx,   by,  bz,   bxx,  byy,  bzz,  bxy,  byz,  bzx, &
                        w,   wx,   wy,  wz,   wxx,  wyy,  wzz,  wxy,  wyz,  wzx, &
                       bw,  bwx,  bwy, bwz,  bwxx, bwyy, bwzz, bwxy, bwyz, bwzx  )

!---

! Compute derivatives of (c*u)=(rho*w*u)


  call derivatives_ab2( c,   cx,   cy,  cz,   cxx,  cyy,  czz,  cxy,  cyz,  czx, &
                        u,   ux,   uy,  uz,   uxx,  uyy,  uzz,  uxy,  uyz,  uzx, &
                       cu,  cux,  cuy, cuz,  cuxx, cuyy, cuzz, cuxy, cuyz, cuzx  )

! Compute derivatives of (c*v)=(rho*w*v)

  call derivatives_ab2( c,   cx,   cy,  cz,   cxx,  cyy,  czz,  cxy,  cyz,  czx, &
                        v,   vx,   vy,  vz,   vxx,  vyy,  vzz,  vxy,  vyz,  vzx, &
                       cv,  cvx,  cvy, cvz,  cvxx, cvyy, cvzz, cvxy, cvyz, cvzx  )

! Compute derivatives of (c*w)=(rho*w*w)

  call derivatives_ab2( c,   cx,   cy,  cz,   cxx,  cyy,  czz,  cxy,  cyz,  czx, &
                        w,   wx,   wy,  wz,   wxx,  wyy,  wzz,  wxy,  wyz,  wzx, &
                       cw,  cwx,  cwy, cwz,  cwxx, cwyy, cwzz, cwxy, cwyz, cwzx  )

!-----------------------------------------------------------------------------

! Compute derivatives of (rho*u*h) = (u*rH)

  call derivatives_ab2(     u,     ux,     uy,     uz,     uxx,    uyy,    uzz,    uxy,    uyz,    uzx, &
                         rhoH,  rhoHx,  rhoHy,  rhoHz,  rhoHxx, rhoHyy, rhoHzz, rhoHxy, rhoHyz, rhoHzx, &
                        rhouH, rhouHx, rhouHy, rhouHz, rhouHxx,rhouHyy,rhouHzz,rhouHxy,rhouHyz,rhouHzx  )

! Compute derivatives of (rho*v*h) = (v*rH)

  call derivatives_ab2(     v,     vx,     vy,     vz,     vxx,    vyy,    vzz,    vxy,    vyz,    vzx, &
                         rhoH,  rhoHx,  rhoHy,  rhoHz,  rhoHxx, rhoHyy, rhoHzz, rhoHxy, rhoHyz, rhoHzx, &
                        rhovH, rhovHx, rhovHy, rhovHz, rhovHxx,rhovHyy,rhovHzz,rhovHxy,rhovHyz,rhovHzx  )

! Compute derivatives of (rho*w*h) = (w*rH)

  call derivatives_ab2(     w,     wx,     wy,     wz,     wxx,    wyy,    wzz,    wxy,    wyz,    wzx, &
                         rhoH,  rhoHx,  rhoHy,  rhoHz,  rhoHxx, rhoHyy, rhoHzz, rhoHxy, rhoHyz, rhoHzx, &
                        rhowH, rhowHx, rhowHy, rhowHz, rhowHxx,rhowHyy,rhowHzz,rhowHxy,rhowHyz,rhowHzx  )

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!
! Compute the inviscid terms in the source term array, source(:).
!
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 !----------------------------------------------------------------------------
 ! Continuity:              (rho*u)_x   +   (rho*v)_y       + (rho*w)_z
         s(iconti)  = (rhox*u + rho*ux) + (rhoy*v + rho*vy) + (rhoz*w + rho*wz)

 !----------------------------------------------------------------------------
 ! Momentum:          (rho*u*u)_x + (rho*u*v)_y + (rho*u*v)_z + px
         s(ixmom)   =     aux     +    buy      + cuz         + px

 !----------------------------------------------------------------------------
 ! Momentum:          (rho*v*u)_x + (rho*v*v)_y + (rho*v*w)_z + py
         s(iymom)   =     avx     +    bvy      +  cvz        + py

 !----------------------------------------------------------------------------
 ! Momentum:          (rho*w*u)_x + (rho*w*v)_y + (rho*w*w)_z + pz
         s(izmom)   =     awx     +    bwy      +  cwz        + pz

 !----------------------------------------------------------------------------
 ! Energy:            (rho*u*H)_x + (rho*v*H) + (rho*w*H)
         s(ienrgy)  =    rhouHx   +   rhovHy  +   rhowHz

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!
! Compute derivatives of the source terms.
!
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 !----------------------------------------------------------------------------
 ! Derivatives of Continuity:
 !       s(iconti)  = (rhox*u + rho*ux) + (rhoy*v + rho*vy) + (rhoz*w + rho*wz)

        sx(iconti)  = (rhoxx*u + rhox*ux + rhox*ux + rho*uxx) & != d(rhox*u + rho*ux)/dx
                    + (rhoxy*v + rhoy*vx + rhox*vy + rho*vxy) & != d(rhoy*v + rho*vy)/dx
                    + (rhozx*w + rhoz*wx + rhox*wz + rho*wzx)   != d(rhoz*w + rho*wz)/dx

        sy(iconti)  = (rhoxy*u + rhox*uy + rhoy*ux + rho*uxy) & != d(rhox*u + rho*ux)/dy
                    + (rhoyy*v + rhoy*vy + rhoy*vy + rho*vyy) & != d(rhoy*v + rho*vy)/dy
                    + (rhoyz*w + rhoz*wy + rhoy*wz + rho*wyz)   != d(rhoz*w + rho*wz)/dy

        sz(iconti)  = (rhozx*u + rhox*uz + rhoz*ux + rho*uzx) & != d(rhox*u + rho*ux)/dz
                    + (rhoyz*v + rhoy*vz + rhoz*vy + rho*vyz) & != d(rhoy*v + rho*vy)/dz
                    + (rhozz*w + rhoz*wz + rhoz*wz + rho*wzz)   != d(rhoz*w + rho*wz)/dz

 !----------------------------------------------------------------------------
 ! Derivatives of the x-momentum equation:
 !        s(ixmom)   =     aux     +    buy       + cuz         + px

         sx(ixmom)   =     auxx    +    buxy      + cuzx        + pxx
         sy(ixmom)   =     auxy    +    buyy      + cuyz        + pxy
         sz(ixmom)   =     auzx    +    buyz      + cuzz        + pzx

 !----------------------------------------------------------------------------
 ! Derivatives of the y-momentum equation:
 !        s(iymom)   =     avx     +    bvy       +  cvz        + py

         sx(iymom)   =     avxx    +    bvxy      + cvzx        + pxy
         sy(iymom)   =     avxy    +    bvyy      + cvyz        + pyy
         sz(iymom)   =     avzx    +    bvyz      + cvzz        + pyz

 !----------------------------------------------------------------------------
 ! Derivatives of the z-momentum equation:
 !        s(izmom)   =     awx     +    bwy      +  cwz         + pz

         sx(izmom)   =     awxx    +    bwxy      + cwzx        + pzx
         sy(izmom)   =     awxy    +    bwyy      + cwyz        + pyz
         sz(izmom)   =     awzx    +    bwyz      + cwzz        + pzz

 !----------------------------------------------------------------------------
 ! Derivatives of the energy equation:
 !        s(ienrgy)  =    rhouHx   +   rhovHy  +   rhowHz

         sx(ienrgy)  =    rhouHxx  +   rhovHxy +  rhowHzx
         sy(ienrgy)  =    rhouHxy  +   rhovHyy +  rhowHyz
         sz(ienrgy)  =    rhouHzx  +   rhovHyz +  rhowHzz

 end subroutine compute_manufactured_sol_and_source

!********************************************************************************
!********************************************************************************


!********************************************************************************
!
! This subroutine computes 1st and 2nd derivatives of a quadratic term
!
!  Input: a, ax, ay, az !Function value a, and its derivatives.
!         b, bx, by, bz !Function value b, and its derivatives.
!
! Output: ab = a*b, grad(ab), and second derivatives.
!
!********************************************************************************
 subroutine derivatives_ab2( a, ax, ay, az, axx, ayy, azz, axy, ayz, azx, &
                             b, bx, by, bz, bxx, byy, bzz, bxy, byz, bzx, &
                            ab,abx,aby,abz,abxx,abyy,abzz,abxy,abyz,abzx  )

 implicit none
 integer , parameter ::  p2 = selected_real_kind(p=15)
 real(p2), parameter :: two = 2.0_p2

!Input
 real(p2), intent( in) ::  a,  ax,  ay, az, axx, ayy, azz, axy, ayz, azx
 real(p2), intent( in) ::  b,  bx,  by, bz, bxx, byy, bzz, bxy, byz, bzx

!Output
 real(p2), intent(out) :: ab, abx, aby, abz, abxx, abyy, abzz, abxy, abyz, abzx

  ab    = a*b

  abx   = ax*b + a*bx
  aby   = ay*b + a*by
  abz   = az*b + a*bz

  abxx  = axx*b +   two*ax*bx   + a*bxx
  abyy  = ayy*b +   two*ay*by   + a*byy
  abzz  = azz*b +   two*az*bz   + a*bzz

  abxy  = axy*b + ax*by + ay*bx + a*bxy
  abyz  = ayz*b + ay*bz + az*by + a*byz
  abzx  = azx*b + az*bx + ax*bz + a*bzx

 end subroutine derivatives_ab2

!********************************************************************************
!*
!* This function computes the exponential function:
!*
!*       f = a0 + a*exp( sx*x + sy*y + sz*z )
!*
!* and its derivatives:
!*
!*     df/dx^nx/dy^ny/dz^nz = d^{nx+ny_nz}( a0 + a*exp( sx*x + sy*y + sz*z ) )/(dx^nx*dy^ny*dz^nz)
!*
!* depending on the input parameters:
!*
!*
!* Input:
!*
!*  a0,a,sx,sy,sz = coefficients in the function: f = a0 + a*exp( sx*x + sy*y + sz*z )
!*            x = x-coordinate at which the function/derivative is evaluated.
!*            y = y-coordinate at which the function/derivative is evaluated.
!*            z = z-coordinate at which the function/derivative is evaluated.
!*           nx = nx-th derivative with respect to x (nx >= 0).
!*           ny = ny-th derivative with respect to y (ny >= 0).
!*           nz = nz-th derivative with respect to y (ny >= 0).
!*
!* Output: The function value.
!*
!********************************************************************************
 function manufactured_sol(a0,a,sx,sy,sz, nx,ny,nz,x,y,z) result(fval)

 implicit none
 integer , parameter  :: p2 = selected_real_kind(p=15)

!Input
 real(p2), intent(in) :: a0, a, sx, sy, sz, x, y, z
 integer , intent(in) :: nx, ny, nz

!Output
 real(p2)             :: fval

  !- Function value.

  if     ( nx == 0 .and. ny == 0 .and. nz == 0) then

   fval = a0 + a*exp( sx*x + sy*y + sz*z )

  !- Derivatives w.r.t. a single coordinate.

  elseif ( nx  > 0 .and. ny == 0 .and. nz == 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(nx*sx)

  elseif ( nx == 0 .and. ny  > 0 .and. nz == 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(ny*sy)

  elseif ( nx == 0 .and. ny == 0 .and. nz  > 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(nz*sz)

  !- Derivatives w.r.t. two coordinates.

  elseif ( nx == 0 .and. ny  > 0 .and. nz  > 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(ny*sy)*(nz*sz)

  elseif ( nx  > 0 .and. ny == 0 .and. nz  > 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(nx*sx)*(ny*sy)

  elseif ( nx  > 0 .and. ny  > 0 .and. nz == 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(nx*sx)*(ny*sy)

  !- Derivatives w.r.t. three coordinates.

  elseif ( nx  > 0 .and. ny  > 0 .and. nz  > 0) then

   fval = a*exp( sx*x + sy*y + sz*z )*(nx*sx)*(ny*sy)*(nz*sz)

  !-

  else

   write(*,*) " This cannot happen. Stop."
   stop

  endif

 end function manufactured_sol
!********************************************************************************

 end module edu3d_mms_euler


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
!* This is Version 2 (December 16, 2017).
!*
!* This F90 code is written and made available for an educational purpose.
!* This file may be updated in future.
!*
!* Katate Masatsuka, http://www.cfdbooks.com
!********************************************************************************

 module edu3d_mms_euler


 contains

!********************************************************************************
!*
!* Educationally-Designed Unstructured 2D (EDU2D) Code
!*
!*  --- EDU2D_MMS
!*
!*
!* This subroutine computes, for a given position (x,y), the manufactured
!* solution values and the source terms. The manufactured solution takes
!* the following form:
!*
!*    c0 + cs*sin(cx*x+cy*y)
!*
!* with the coefficients, c0, cs, cx, cy chosen differently for different
!* variable.
!*
!* ------------------------
!* Input:
!*
!*   (x,y) = Position at which the solution and source terms are sought.
!*
!* ------------------------
!* Output:
!*                       1 2 3 4     5     6     7  8  9   10 11
!*         w(1:11) =  [rho,u,v,p,tauxx,tauyy,tauxy,qx,qy,vort, T]
!*
!*   gradw(1:11,1) = d[rho,u,v,p,tauxx,tauyy,tauxy,qx,qy,vort, T]/dx
!*
!*   gradw(1:11,2) = d[rho,u,v,p,tauxx,tauyy,tauxy,qx,qy,vort, T]/dy
!*
!*     source(1:5) = Source terms for continuity, x/y-momentum, and energy eqns
!*
!* ------------------------
!*
!* Note: Second derivatives are computed only for [rho,u,v,p]
!*
!*         uxx, uyy, uxy, vxx, vyy, vxy, rhoxx, rhoyy, rhoxy, pxx, pyy, pxy.
!*
!*       just because these are required anyway to compute the first-derivatives
!*       of [tauxx,tauyy,tauxy,qx,qy,vort].
!*
!* Note: You can extend this subroutine to compute higher-order derivatives
!*       of the solution and the source terms as well. Personally, I have
!*       an extended version that computes up to second derivatives of the
!*       source terms.
!*
!* Note: You can extend this also to the unsteady NS equations by including
!*       time derivatives when generating the source term: s=dv/dt+df(v)/dx+dg(v)/dy.
!*
!* Note: For high-order schemes, you can use this subroutine to compute the
!*       source terms at quadrature points.
!*
!*
!* This is Version 2 (December 16, 2017).
!*
!* This F90 code is written and made available for an educational purpose.
!* This file may be updated in future.
!*
!* Katate Masatsuka, http://www.cfdbooks.com
!********************************************************************************

 subroutine compute_manufactured_sol_and_source(x,y,M_inf,gamma, w,u,s,sx,sy,sz )

 implicit none
 integer , parameter :: p2 = selected_real_kind(p=15)

 real(p2), parameter ::         pi = 3.141592653589793238_p2
 real(p2), parameter ::       half = 0.5_p2
 real(p2), parameter ::        one = 1.0_p2
 real(p2), parameter ::        two = 2.0_p2

!Input 
 real(p2),               intent( in) :: x, y, z
 real(p2),               intent( in) :: M_inf, gamma

!Output
 real(p2), dimension(5), intent(out) :: w !Primitive variables    (rho,u,v,w,p)
 real(p2), dimension(5), intent(out) :: u !Conservative variables (rho,rho*u,rho*v,rho*w,rho*E)
 real(p2), dimension(5), intent(out) :: s, sx, sy, sz

!Local variables

 real(p2) :: rho , u , v , p
 real(p2) :: rhox, ux, vx, px
 real(p2) :: rhox, ux, vx, px
 real(p2) :: rhoz, uz, vz, pz

 real(p2) ::  a,  ax,  ay, az
 real(p2) ::  b,  bx,  by, bz
 real(p2) ::  c,  cx,  cy, cz
 real(p2) ::  k,  kx,  ky, kz

 real(p2) :: u2   , u2x   , u2y
 real(p2) :: v2   , v2x   , v2y
 real(p2) :: w2   , w2x   , w2y

 real(p2) :: au   , aux   , auy
 real(p2) :: av   , avx   , avy
 real(p2) :: aw   , awx   , awy

 real(p2) :: bu   , bux   , buy
 real(p2) :: bv   , bvx   , bvy
 real(p2) :: bw   , bwx   , bwy

 real(p2) :: cu   , cux   , cuy
 real(p2) :: cv   , cvx   , cvy
 real(p2) :: cw   , cwx   , cwy

 real(p2) :: rhoH , rhoHx , rhoHy ,  rhoHz
 real(p2) :: rhouH, rhouHx, rhouHy, rhouHz
 real(p2) :: rhovH, rhovHx, rhovHy, rhovHz
 real(p2) :: rhowH, rhowHx, rhowHy, rhowHz

 real(p2) :: rhoxx, rhoyy, rhozz, rhoxy, rhoyz, rhozx
 real(p2) ::   uxx,   uyy,   uzz,   uxy,   uyz,   uzx
 real(p2) ::   vxx,   vyy,   vzz,   vxy,   vyz,   vzx
 real(p2) ::   wxx,   wyy,   wzz,   wxy,   wyz,   wzx
 real(p2) ::   pxx,   pyy,   pzz,   pxy,   pyz,   pzx

 real(p2) ::   axx,   ayy,   axy, ayz, azx
 real(p2) ::   bxx,   byy,   bxy, byz, bzx
 real(p2) ::   cxx,   cyy,   cxy, cyz, czx

 integer  :: ix, iy, iz
 integer  :: irho, iu, iv, iw, ip
 integer  :: iconti, ixmom, iymom, izmom, ienrgy

 real(p2) :: cr0, crs, crx, cry, crz
 real(p2) :: cu0, cus, cux, cuy, cuz
 real(p2) :: cv0, cvs, cvx, cvy, cvz
 real(p2) :: cw0, cws, cwx, cwy, cwz
 real(p2) :: cp0, cps, cpx, cpy, cpz

!-----------------------------------------------------------
! Define some indices

        ix =  1 ! x component
        iy =  2 ! y component
        iz =  3 ! y component

      irho =  1 ! density
      iu   =  2 ! x-velocity
      iv   =  3 ! y-velocity
      iw   =  4 ! w-velocity
      ip   =  5 ! pressure

    iconti =  1 ! continuity equation
    ixmom  =  2 ! x-momentum equation
    iymom  =  3 ! y-momentum equation
    izmom  =  4 ! y-momentum equation
    ienrgy =  5 !     energy equation

!-----------------------------------------------------------
! Constants for the exact solution: c0 + cs*sin(cx*x+cy*y).
!
! Note: Make sure the density and pressure are positive.
! Note: These values are passed to the subroutine, 
!         manufactured_sol(c0,cs,cx,cy, nx,ny,x,y),
!       whcih returns the solution value or derivative.

 !-----------------------------------------
 ! Density    = cr0 + crs*exp(crx*x+cry*y+crz*z)

  cr0 =  1.0_p2
  crs =  0.2_p2
  crx =  0.5_p2
  cry =  0.5_p2
  crz =  0.5_p2

 !-----------------------------------------
 ! X-velocity = cu0 + cus*exp(cux*x+cuy*y+cuz*z)

  cu0 =  0.3_p2
  cus =  0.2_p2
  cux =  0.5_p2
  cuy =  0.5_p2
  cuz =  0.5_p2

 !-----------------------------------------
 ! Y-velocity = cv0 + cvs*exp(cvx*x+cvy*y+cvz*z)

  cv0 =  0.2_p2
  cvs =  0.2_p2
  cvx =  0.5_p2
  cvy =  0.5_p2
  cvz =  0.5_p2

 !-----------------------------------------
 ! Z-velocity = cw0 + cws*exp(cwx*x+cwy*y+cwz*z)

  cw0 =  0.1_p2
  cws =  0.2_p2
  cwx =  0.5_p2
  cwy =  0.5_p2
  cwz =  0.5_p2

 !-----------------------------------------
 ! Pressure  = cp0 + cps*exp(cpx*x+cpy*y+cpz*z)

  cp0 =  1.0_p2
  cps =  0.2_p2
  cpx =  0.5_p2
  cpy =  0.5_p2
  cpz =  0.5_p2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Part I: Compute w = [rho,u,v,p,tauxx,tauyy,tauxy,qx,qy,vort,T] and grad(w).
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

 !------------------------------------------------------------------------
 ! rho: Density and its derivatives
 
              rho = manufactured_sol(cr0,crs,crx,cry,crz, 0,0,0,x,y,z)
             rhox = manufactured_sol(cr0,crs,crx,cry,crz, 1,0,0,x,y,z)
             rhoy = manufactured_sol(cr0,crs,crx,cry,crz, 0,1,0,x,y,z)
             rhoz = manufactured_sol(cr0,crs,crx,cry,crz, 0,0,1,x,y,z)

            rhoxx = manufactured_sol(cr0,crs,crx,cry,crz, 2,0,0,x,y,z)
            rhoyy = manufactured_sol(cr0,crs,crx,cry,crz, 0,2,0,x,y,z)
            rhozz = manufactured_sol(cr0,crs,crx,cry,crz, 0,0,2,x,y,z)

            rhoxy = manufactured_sol(cr0,crs,crx,cry,crz, 1,1,0,x,y,z)
            rhoyz = manufactured_sol(cr0,crs,crx,cry,crz, 0,1,1,x,y,z)
            rhozx = manufactured_sol(cr0,crs,crx,cry,crz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! u: x-velocity and its derivatives

              u = manufactured_sol(cu0,cus,cux,cuy,cuz, 0,0,0,x,y,z)
             ux = manufactured_sol(cu0,cus,cux,cuy,cuz, 1,0,0,x,y,z)
             uy = manufactured_sol(cu0,cus,cux,cuy,cuz, 0,1,0,x,y,z)
             uz = manufactured_sol(cu0,cus,cux,cuy,cuz, 0,0,1,x,y,z)

            uxx = manufactured_sol(cu0,cus,cux,cuy,cuz, 2,0,0,x,y,z)
            uyy = manufactured_sol(cu0,cus,cux,cuy,cuz, 0,2,0,x,y,z)
            uzz = manufactured_sol(cu0,cus,cux,cuy,cuz, 0,0,2,x,y,z)

            uxy = manufactured_sol(cu0,cus,cux,cuy,cuz, 1,1,0,x,y,z)
            uyz = manufactured_sol(cu0,cus,cux,cuy,cuz, 0,1,1,x,y,z)
            uzx = manufactured_sol(cu0,cus,cux,cuy,cuz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! v: y-velocity and its derivatives

              v = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 0,0,0,x,y,z)
             vx = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 1,0,0,x,y,z)
             vy = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 0,1,0,x,y,z)
             vz = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 0,0,1,x,y,z)

            vxx = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 2,0,0,x,y,z)
            vyy = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 0,2,0,x,y,z)
            vzz = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 0,0,2,x,y,z)

            vxy = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 1,1,0,x,y,z)
            vyz = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 0,1,1,x,y,z)
            vzx = manufactured_sol(cv0,cvs,cvx,cvy,cvz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! w: z-velocity and its derivatives

              w = manufactured_sol(cw0,cws,cwx,cwy,cwz, 0,0,0,x,y,z)
             wx = manufactured_sol(cw0,cws,cwx,cwy,cwz, 1,0,0,x,y,z)
             wy = manufactured_sol(cw0,cws,cwx,cwy,cwz, 0,1,0,x,y,z)
             wz = manufactured_sol(cw0,cws,cwx,cwy,cwz, 0,0,1,x,y,z)

            wxx = manufactured_sol(cw0,cws,cwx,cwy,cwz, 2,0,0,x,y,z)
            wyy = manufactured_sol(cw0,cws,cwx,cwy,cwz, 0,2,0,x,y,z)
            wzz = manufactured_sol(cw0,cws,cwx,cwy,cwz, 0,0,2,x,y,z)

            wxy = manufactured_sol(cw0,cws,cwx,cwy,cwz, 1,1,0,x,y,z)
            wyz = manufactured_sol(cw0,cws,cwx,cwy,cwz, 0,1,1,x,y,z)
            wzx = manufactured_sol(cw0,cws,cwx,cwy,cwz, 1,0,1,x,y,z)

 !------------------------------------------------------------------------
 ! p: pressure and its derivatives

              p = manufactured_sol(cp0,cps,cpx,cpy,cpz, 0,0,0,x,y,z)
             px = manufactured_sol(cp0,cps,cpx,cpy,cpz, 1,0,0,x,y,z)
             py = manufactured_sol(cp0,cps,cpx,cpy,cpz, 0,1,0,x,y,z)
             pz = manufactured_sol(cp0,cps,cpx,cpy,cpz, 0,0,1,x,y,z)

            pxx = manufactured_sol(cp0,cps,cpx,cpy,cpz, 2,0,0,x,y,z)
            pyy = manufactured_sol(cp0,cps,cpx,cpy,cpz, 0,2,0,x,y,z)
            pzz = manufactured_sol(cp0,cps,cpx,cpy,cpz, 0,0,2,x,y,z)

            pxy = manufactured_sol(cp0,cps,cpx,cpy,cpz, 1,1,0,x,y,z)
            pyz = manufactured_sol(cp0,cps,cpx,cpy,cpz, 0,1,1,x,y,z)
            pzx = manufactured_sol(cp0,cps,cpx,cpy,cpz, 1,0,1,x,y,z)

 !-----------------------------------------------------------------------------
 ! Store the exact solution for return.

  !Primitive variables

      w(irho) = rho 
      w(iu)   = u
      w(iv)   = v
      w(iw)   = w
      w(ip)   = p

  !Conservative variables

         u(1) = rho
         u(2) = rho*u
         u(3) = rho*v
         u(4) = rho*w
         u(5) = p/(gamma-one) + half*rho*(u*u + v*v + w*w)

!-----------------------------------------------------------------------------
!
! Exact (manufactured) solutons and derivatives have been computed.
! We move on the source terms, which are the governing equations we wish
! to solve (the compressible NS) evaluated by the exact solution (i.e., 
! the function we wish to make exact with the source terms).
!
! Navier-Stokes: dF(w)/dx + dG(w)/dy = 0.
!
! Our manufactured solutions, wm, are the exact soltuions to the following:
!
!   dF(w)/dx + dG(w)/dy = S,
!
! where S = dF(wm)/dx + dG(wm)/dy.
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
  call derivatives_ab(u,ux,uy,uz, u,ux,uy,uz,  u2,u2x,u2y,u2z) ! a=u^2

! Derivatives of v^2
  call derivatives_ab(v,vx,vy,vz, v,vx,vy,vz,  v2,v2x,v2y,v2z) ! b=v^2

! Derivatives of w^2
  call derivatives_ab(w,wx,wy,wz, w,wx,wy,wz,  w2,w2x,w2y,w2z) ! c=w^2

! Derivatives of k=(u^2+v^2+w^2)/2

  k     = half*(u*u   + v*v   + w*w)
  kx    = half*(u2x   + v2x   + w2x)
  ky    = half*(u2y   + v2y   + w2y)
  kz    = half*(u2z   + v2z   + w2z)

! Derivatives of rho*k = rho*(u^2+v^2)/2
  call derivatives_ab(rho,rhox,rhoy,rhoz, k,kx,ky,kz,  a,ax,ay,az) !a=rho*(u^2+v^2)/2

! Derivatives of rho*H = gamma/(gamma-1)*p + rho*k

  rhoH    = gamma/(gamma-one)*p    + a
  rhoHx   = gamma/(gamma-one)*px   + ax
  rhoHy   = gamma/(gamma-one)*py   + ay
  rhoHz   = gamma/(gamma-one)*pz   + az

!-----------------------------------------------------------------------------

! Compute derivatives of (rho*u)
  call derivatives_ab(rho,rhox,rhoy,rhoz, u,ux,uy,uz,   a,ax,ay,az) !a=(rho*u)

! Compute derivatives of (rho*v)
  call derivatives_ab(rho,rhox,rhoy,rhoz, v,vx,vy,vz,   b,bx,by,bz) !ab=(rho*v)

! Compute derivatives of (rho*w)
  call derivatives_ab(rho,rhox,rhoy,rhoz, w,wx,wy,wz,   c,cx,cy,cz) !ab=(rho*w)

!-----------------------------------------------------------------------------

! Compute derivatives of (a*u)=(rho*u*u)

  call derivatives_ab(a,ax,ay,az, u,ux,uy,uz,     au,aux,auy,auz)

! Compute derivatives of (a*v)=(rho*u*v)

  call derivatives_ab(a,ax,ay,az, v,vx,vy,vz,     av,avx,avy,avz)

! Compute derivatives of (a*v)=(rho*u*w)

  call derivatives_ab(a,ax,ay,az, w,wx,wy,wz,     aw,awx,awy,awz)

!---

! Compute derivatives of (b*u)=(rho*v*u)

  call derivatives_ab(b,bx,by,bz, u,ux,uy,uz,     bu,bux,buy,buz)

! Compute derivatives of (b*v)=(rho*v*v)

  call derivatives_ab(b,bx,by,bz, v,vx,vy,vz,     bv,bvx,bvy,bvz)

! Compute derivatives of w=(rho*v*w)

  call derivatives_ab(b,bx,by,bz, w,wx,wy,wz,     bw,bwx,bwy,bwz)

!---

! Compute derivatives of (c*u)=(rho*w*u)

  call derivatives_ab(c,cx,cy,cz, u,ux,uy,uz,     cu,cux,cuy,cuz)

! Compute derivatives of (c*v)=(rho*w*v)

  call derivatives_ab(c,cx,cy,cz, v,vx,vy,vz,     cv,cvx,cvy,cvz)

! Compute derivatives of (c*w)=(rho*w*w)

  call derivatives_ab(c,cx,cy,cz, w,wx,wy,wz,     cw,cwx,cwy,cwz)

!-----------------------------------------------------------------------------

! Compute derivatives of (rho*u*h) = (u*rH)

  call derivatives_ab( u,ux,uy,uz, rhoH,rhoHx,rhoHy,rhoHz,  rhouH,rhouHx,rhouHy,rhouHz)

! Compute derivatives of (rho*v*h) = (v*rH)

  call derivatives_ab( v,vx,vy,vz, rhoH,rhoHx,rhoHy,rhoHz,  rhovH,rhovHx,rhovHy,rhovHz)

! Compute derivatives of (rho*w*h) = (w*rH)

  call derivatives_ab( w,wx,wy,wz, rhoH,rhoHx,rhoHy,rhoHz,  rhowH,rhowHx,rhowHy,rhowHz)

!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Store the inviscid terms in the source term array, source(:).
!---------------------------------------------------------------------

 !--------------------------------------------------------------------------
 ! Continuity:              (rho*u)_x   +   (rho*v)_y       + (rho*w)_z
    source(iconti)  = (rhox*u + rho*ux) + (rhoy*v + rho*vy) + (rhoz*v + rho*vz)

 !--------------------------------------------------------------------------
 ! Momentum:          (rho*u*u)_x + (rho*u*v)_y  + (rho*u*v)_z + px
    source(ixmom)   =     aux     +    buy       + cuz         + px

 !--------------------------------------------------------------------------
 ! Momentum:          (rho*v*u)_x + (rho*v*v)_y + (rho*v*w)_z + py
    source(iymom)   =     avx     +    bvy      +  cvz        + py

 !--------------------------------------------------------------------------
 ! Momentum:          (rho*w*u)_x + (rho*w*v)_y + (rho*w*w)_z + pz
    source(izmom)   =     awx     +    bwy      +  cwz        + pz

 !--------------------------------------------------------------------------
 ! Energy:            (rho*u*H)_x + (rho*v*H) + (rho*w*H)
    source(ienrgy)  =    rhouHx   +   rhovHy  +   rhowHz


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
 real(p2), intent( in) ::  a,  ax,  ay, axx, ayy, azz, axy, ayz, azx
 real(p2), intent( in) ::  b,  bx,  by, bxx, byy, bzz, bxy, byz, bzx

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
  abxy  = ayz*b + ay*bz + az*by + a*byz
  abxy  = azx*b + az*bx + ax*bz + a*bzx

 end subroutine derivatives_ab2

!********************************************************************************
!
! This subroutine computes first derivatives of a quadratic term
!
!  Input: a, ax, ay, az !Function value a, and its derivatives, (ax,ay,az).
!         b, bx, by, bz !Function value b, and its derivatives, (bx,by,bz).
!
! Output: ab = a*b, abx = d(a*b)/dx, aby = d(a*b)/dy, abz = d(a*b)/dz.
!
!********************************************************************************
 subroutine derivatives_ab(a,ax,ay,az,  b,bx,by,bz, ab,abx,aby,abz)

 implicit none
 integer , parameter  :: p2 = selected_real_kind(p=15)

!Input
 real(p2), intent( in) ::  a,  ax,  ay, az
 real(p2), intent( in) ::  b,  bx,  by, bz

!Output
 real(p2), intent(out) :: ab, abx, aby, abz

  ab    = a*b
  abx   = ax*b + a*bx
  aby   = ay*b + a*by
  abz   = az*b + a*bz

 end subroutine derivatives_ab

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
 integer , intent(in) :: nx, ny, zn

!Output
 real(p2)             :: fval

  if ( nx*ny*nz == 0) then

   fval = a0 + a*exp( sx*x + sy*y + sz*z )

  else

   fval = a*exp( sx*x + sy*y + sz*z )*(nx*sx)*(ny*sy)*(nz*sz)

  endif

 end function manufactured_sol
!********************************************************************************

 end module edu3d_mms_euler


c************************************************************************
c
c     common blocks for cas.f
c
c     nm: maximum number of contours allowed
c     npm: maximum total number of nodes on all contours
c
c     x,y,z position of each node [cartesian co-ordinates]
c     u,v,w velocity at each node [u=dx/dt, etc.]
c
c************************************************************************
      parameter (pi=3.1415926535897932385d00,pifac=pi/180.)
      parameter (npm=100000,nm=50)
      logical corn(npm)
      common /xyuv/ x(npm),y(npm),z(npm),u(npm),v(npm),w(npm)
      common /runga/ xi(npm),yi(npm),zi(npm),xf(npm),yf(npm),zf(npm)
      common /extr1/ o(npm),eu(npm),alp(npm),bet(npm),gam(npm),es(npm)
      common /extr2/ a(npm),b(npm),c(npm),d(npm)
      common /extr3/ e(npm),f(npm),g(npm),h(npm)
      common /d01/ dx(npm),dy(npm),dz(npm),ax(npm),ay(npm),az(npm)
      common /d02/ sx(npm),sy(npm),sz(npm),corn
      common /cind/ np(nm),mb(nm),ind(nm),om(nm)
      common /dunno/ amu,ell,emin,dm,n,npt
      common /other/ t,dt,dt2,dt3,dt6,small,small3
      common /vel0/ u0(289,289),v0(289,289),w0(289,289)
      common /vel1/ u1(289,289),v1(289,289),w1(289,289)
      common /grid/ xlat(289,289),xlong(289,289)
      common /dim/dlat,dlong,ms,ns,nlat,nlong
      common /ttime/ t0,t1,deltat
c
      common /c3/ along(npm),alat(npm)
c

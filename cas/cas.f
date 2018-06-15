ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Contour Advection (with Surgery) code
c
c     Winds from NMC; read from daily files with NH data only in winds/
c     Change setting of prefix variable to read from different
c     directory
c
c     Code which advects the material contours (represented by
c     a series of nodes) using velocities stored in a grid.
c     The velocity at each node is determined by linear interpolation
c     (in space and time) from the grid.
c     Contour representation and surgery is the same as contour
c     surgery code of D.G.Dritschel.
c
c     Link with "uvw.f" - contains subroutines "r_grid" & "uvw"
c
c     Original CA code written by D.W.Waugh, 1992.  (based on contour
c     surgery code of D.G.Dritschel)
c     Modified by D.W. Waugh, Sept. 1997.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      program cas
c**********************************************************************
c
c     main program
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      character*44 titlea,titleb,titlec
      character*25 infile,outfile,logfile
      character*6 date,date1,date2
      character*5 daten, daten2
c      character*48 prefix
c      character*59 ufile
      character*6 prefix
c      character*17 ufile
      character*12 ufile
      character*3 surf
c
c     *** READ FROM SCREEN ***
c
c     read in input/output filenames
c     ------------------------------
      write(6,1040)
 1040 format('Input contour filename ',$)
      read(5,1060)infile
 1060 format(a25)
      write(6,1050)
 1050 format('Output contour filename ',$)
      read(5,1060)outfile
      write(6,1055)
 1055 format('Log filename ',$)
      read(5,1060)logfile
c
c     read in start date and surface
c     ------------------------------
      write(6,1000)
c 1000 format('Enter initial date (iy im id) ',$)
c      read(5,*)iyr,imon,iday
 1000 format('Enter initial day',$)
      read(5,*)iday
c      write(6,1002)
c 1002 format('Enter surface (a3) ',$)
c      read(5,'(a3)')surf
c
c     OPEN  FILES
c     ----------
      open(11,file=infile)
      open(26,file=outfile)
      open(99,file=logfile)
c
c     READ INPUT CONTOUR FILE
c     ------------------------
      read(11,5) titlea
      read(11,5) titleb
      read(11,5) titlec
      read(11,*) ndays,ns,deltat,dsave,amu,dm
c
c     ndays  = duration of calculation (in days)
c     ns     = number of timesteps per day
c     deltat = time between gridded data (in days)
c     dsave  = time between output data (in days)
c     amu,dm = contour representation parameters
c     dt = timestep in advection scheme (in days)
c
      dt=1.0d0/dfloat(ns)
      nss=int(ns*deltat)
      nsave=int(dsave*ns)
      nloops=int(ndays/deltat)
c
c     read in type of node redistribution
c     -----------------------------------
c 4002 write(6,3010)
c 3010 format('Renode: Dritschel (1), simple (2), or none (3): ',$)
c      read(5,*)irenode
c      if (irenode.lt.1 .or. irenode.gt.3) goto 4002
      irenode=1
c
c     if surgery required then read in frequency of surgery in days
c     sfreq=1.0 -> surgery done once per day
c     --------------------------------------
      if (dm.ne.0.) then
         write(6,1003)
 1003    format('Frequency of surgery (days) ',$)
         read(5,*)sfreq
      end if
      nsurg=int(sfreq*ns)
c
      if (deltat.ne.0.25 .and. deltat.ne.0.5 .and.
     &    deltat.ne.1.0) then
         write(*,*)'deltat must equal .25,.5, or 1'
         STOP
      else
         idh=int(24.0*deltat)
      end if
      write(*,*)nloops,nss,nsave,amu,dm,idh,nsurg
      write(99,*)nloops,nss,nsave,amu,dm,idh,nsurg
c
c     read in contours
c     ----------------
      read(11,*) n,npt,t
      do 11 j=1,n
         read(11,*) np(j),mb(j),ind(j),om(j)
 11   continue
      do 10 l=1,npt
         read(11,*) along(l),alat(l)
 10   continue
c
c     convert (lat,long) (to x,y,z)
c     -----------------------------
      call lat_to_xyz
c
5     format(a44)
7     format(f10.7,1x,i3)
15    format(i3,1x,i3,2(f6.3,1x),2(f14.10,1x))
17    format(2(i5,1x),f14.10)
18    format(3(i5,1x),f14.10)
19    format(f14.10,1x,f14.10,1x,f14.10)
c
c     save header info. in output file
c     ---------------------------------
      titlea='BOB  winds, ??? renode'
c      write(titlea(13:15),'(a3)')surf
      if (irenode.eq.1) then
         write(titlea(13:15),'(a3)')'DGD'
      else if (irenode.eq.2) then
         write(titlea(13:15),'(a3)')'MRS'
      else
         write(titlea(13:15),'(a3)')' NO'
      end if

c      titleb='start=yymmdd12, anals until yymmdd12'
c      call form_date(iyr,imon,iday,date)
c      write(titleb(7:12),'(a6)')date
       titleb='start=?????, anals until ?????'
       call form_num_date(iday,daten)
       write(titleb(7:11),'(a5)')daten
c       write(titleb(26:30),'(a5)')
c
c     check how many days of analyses
c     -------------------------------
         prefix='winds/'
c      jyr=iyr
c      jmon=imon
c      jday=iday
c      do i=1,nloops
c         call next_day(jyr,jmon,jday)
c         call form_date(jyr,jmon,jday,date1)
c         ufile=prefix//'u'//date1//'.'//surf
c         open(51,file=ufile,status='old',err=7010)
c         close(51)
c         date2=date1
c      end do
c 7010 write(titleb(29:34),'(a6)')date2
c      write(*,*)date2
c
c      iday = 1
c      do i=1,50
c        call next_num_day(iday)
c        call form_num_date(iday,daten)
c        write(*,*)daten
c      end do
c
      write(26,5) titlea
      write(26,5) titleb
      write(26,5) titlec
      write(26,15) ndays,ns,deltat,dsave,amu,dm
c
c     set-up velocity grid arrays
c     ---------------------------
      tt=-99.
c         prefix='winds/'
         call r_grid1(tt)
c
c     read in velocity grid at initial t
c     ----------------------------------
c         call form_date(iyr,imon,iday,date)
c         ufile=prefix//'u'//date//'.'//surf
c         open(31,file=ufile,status='old')
c         ufile=prefix//'v'//date//'.'//surf
c         open(32,file=ufile,status='old')
c         write(99,*)t,' ',ufile
c         write(*,*)t,' ',ufile
          call form_num_date(iday,daten)
          ufile=prefix//'u'//daten
          open(31,file=ufile,status='old')
          ufile=prefix//'v'//daten
          open(32,file=ufile,status='old')
          write(99,*)t,' ',ufile
          write(*,*)t,' ',ufile
c
         tt=0
         call r_grid1(tt)
c
         close(31)
         close(32)
c
c     contour representation parameters
c     ---------------------------------
      ell = 1.0
      emin=(dm/3.)*(1.-1./(1.+(24.*amu)**(-1./3.)))
      dmsq=dm**2
      cosdm=dcos(dm)
      small=1.d-13
      small3=small*small*small
      dt2=dt/2.
      dt3=dt/3.
      dt6=dt/6.
c
c     redistribute the nodes on all
c     contours according to curvature.
c     --------------------------------
      if (irenode.eq.0) then
         call renode(irenode)
      else if (irenode.eq.1) then
         call renode(irenode)
      else if (irenode.eq.2) then
         call add_nodes
      end if
c
c     dump initial contours to output file
c     -------------------------------------
      call dump
c
c     ********************
c     START OF CALCULATION
c     ********************
c
      t1=t
      idump=0
      isurg=0
      do 1001 loop=1,nloops
c
c     u grid at t0 equals grid at t1
c     ------------------------------
         t0=t1
         do 1006 i=1,nlat
            do 1006 j=1,nlong
               u0(i,j)=u1(i,j)
               v0(i,j)=v1(i,j)
               w0(i,j)=w1(i,j)
 1006    continue
c
c     read in velocity grid at t=t1
c     -----------------------------
c
c         call next_day(iyr,imon,iday)
c         call form_date(iyr,imon,iday,date)
c
c         ufile=prefix//'u'//date//'.'//surf
c         open(31,file=ufile,status='old')
c         ufile=prefix//'v'//date//'.'//surf
c         open(32,file=ufile,status='old')
c         write(99,*)t,' ',ufile
c         write(*,*)t,' ',ufile
c
         call next_num_day(iday)
         call form_num_date(iday,daten)
c
         ufile=prefix//'u'//daten
         open(31,file=ufile,status='old')
         ufile=prefix//'v'//daten
         open(32,file=ufile,status='old')
         write(99,*)t,' ',ufile
         write(*,*)t,' ',ufile
c
         call r_grid1(tt)
         close(31)
         close(32)
c
         if (t1.eq.-999) goto 1001
         t1=t0+deltat
c
         do 900 istep=1,nss
            idump=idump+1
            isurg=isurg+1
c
c     integrates the equations of motion from time t to time t + dt
c     using a runga-kutta integration scheme.  'advect' returns the
c     new nodal positions for the contours.
c     -------------------------------------------------------------
            call advect(idata)
c
c     perform contour surgery
c     -----------------------
            if (dm.ne.0) then
               if (isurg.eq.nsurg) then
                  isurg = 0
                  call surgery
                  write(99,*)'Perform surgery at t =',t
                  if (n.gt.nm) then
                     write(99,*)'N EXCEEDS NM'
                     stop
                  end if
               end if
            end if
c
c     redistribute nodes
c     ------------------
            if (irenode.eq.0) then
               call renode(irenode)
            else if (irenode.eq.1) then
               call renode(irenode)
            else if (irenode.eq.2) then
               call add_nodes
            end if
            write(99,*)'t = ',t,' npt = ',npt
            if (npt.gt.npm) then
               write(99,*)'NPT EXCEEDS NPM'
               STOP
            end if
c
c     save output to file
c     -------------------
            if (idump.eq.nsave) then
               idump=0
               call dump
            end if
c
 900     continue
 1001 continue
c
      stop
      end
c
c
      subroutine advect(idata)
c**********************************************************************
c
c     advects all of the nodes by a 4th order runga-kutta scheme.
c     dt2=dt/2, dt3=dt/3, dt6=dt/6
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      do 5 j=1,n
      do 5 i=mb(j)+1,mb(j)+np(j)
5     o(i)=om(j)
c
      do 10 i=1,npt
      xi(i)=x(i)
      yi(i)=y(i)
      zi(i)=z(i)
      xf(i)=x(i)
      yf(i)=y(i)
10    zf(i)=z(i)

      call uvw1
c
      t = t + dt2
      do 20 i=1,npt
      x(i)=xi(i)+dt2*u(i)
      y(i)=yi(i)+dt2*v(i)
      z(i)=zi(i)+dt2*w(i)
      xf(i)=xf(i)+dt6*u(i)
      yf(i)=yf(i)+dt6*v(i)
20    zf(i)=zf(i)+dt6*w(i)
      call adjust
      call uvw1
c
      do 30 i=1,npt
      x(i)=xi(i)+dt2*u(i)
      y(i)=yi(i)+dt2*v(i)
      z(i)=zi(i)+dt2*w(i)
      xf(i)=xf(i)+dt3*u(i)
      yf(i)=yf(i)+dt3*v(i)
30    zf(i)=zf(i)+dt3*w(i)
      call adjust
      call uvw1
c
      t = t + dt2
      do 40 i=1,npt
      x(i)=xi(i)+dt*u(i)
      y(i)=yi(i)+dt*v(i)
      z(i)=zi(i)+dt*w(i)
      xf(i)=xf(i)+dt3*u(i)
      yf(i)=yf(i)+dt3*v(i)
40    zf(i)=zf(i)+dt3*w(i)
      call adjust
      call uvw1
c
      do 50 i=1,npt
      x(i)=xf(i)+dt6*u(i)
      y(i)=yf(i)+dt6*v(i)
50    z(i)=zf(i)+dt6*w(i)
      call adjust
c
      return
      end
c
c
      subroutine dump
c**********************************************************************
c
c     dumps output to file (unit=26)
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
      real*4 ctime,tarray(2)
c
      ctime=etime(tarray)
      write(99,*)'cpu = ',tarray(1)
c
      write(26,13) n,npt,t
      write(99,13) n,npt,t
c
      do 10 j=1,n
         write(99,15) np(j),mb(j),ind(j),om(j)
         write(26,15) np(j),mb(j),ind(j),om(j)
 10   continue
c
c
c     convert to x,y,x
c
      call xyz_to_lat
c
      do 20 i=1,npt
         write(26,9) along(i),alat(i)
 20   continue
c
13    format(2(i5,1x),f14.10)
14    format(f14.10,1x,f14.10,1x,f14.10)
15    format(3(i5,1x),f14.10)
 9    format(f14.10,1x,f14.10)
c
      return
      end
c
c
      subroutine lat_to_xyz
c********************************************************************
c
c     convert latitude,longitude to x,y,z
c
c********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      do k=1,npt
         alat(k)=alat(k)*pifac
         along(k)=along(k)*pifac
         x(k)=dcos(along(k))*dcos(alat(k))
         y(k)=dsin(along(k))*dcos(alat(k))
         z(k)=dsin(alat(k))
      end do
c
      return
      end
c
c
      subroutine xyz_to_lat
c********************************************************************
c
c     convert x,y,z to latitude,longitude
c
c********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      do k=1,npt
         if (z(k).gt.1.0)  z(k)=1.0
         if (z(k).lt.-1.0) z(k)=-1.0
         alat(k)=asin(z(k))/pifac
         along(k)=atan2(y(k),x(k))/pifac
         if (along(k).lt.0.) along(k)=360.+along(k)
      end do
c
      return
      end
c
cCCCCCCC  Routines to form dates  CCCCCCCCCCCCCC
c
      subroutine next_day(iy,im,id)
c***************************************************************
c     determines next day
c***************************************************************
c
      implicit real*8(a-h,o-z)
c
      if (iy.eq.92 .and. im.eq.2 .and. id.eq.29) then
         id=1
         im=3
      else if (iy.ne.92. .and. im.eq.2 .and. id.eq.28) then
         id=1
         im=3
      else if (id.lt.30) then
         id=id+1
      else if (id.eq.30) then
         if (im.eq.1 .or. im.eq.3 .or. im.eq.5 .or. im.eq.7
     %  .or. im.eq.8 .or. im.eq.10 .or. im.eq.12) then
            id=id+1
         else
            id=1
            im=im+1
         end if
      else if (id.eq.31) then
         id=1
         if (im.eq.12) then
            iy=iy+1
            im=1
         else
            im=im+1
         end if
      end if
c
      return
      end
c
      subroutine form_date(iyr,imon,iday,date1)
c***************************************************************
c     forms 6 character date
c***************************************************************
c
      implicit real*8(a-h,o-z)
      character*6 date1
c
      write(date1(1:2),'(i2)')iyr
      if (imon.lt.10) then
         write(date1(3:3),'(i1)')0
         write(date1(4:4),'(i1)')imon
      else
         write(date1(3:4),'(i2)')imon
      end if
      if (iday.lt.10) then
         write(date1(5:5),'(i1)')0
         write(date1(6:6),'(i1)')iday
      else
         write(date1(5:6),'(i2)')iday
      end if
c
      return
      end
c
      subroutine next_num_day(id)
c***************************************************************
c     determines next day
c***************************************************************
      id=id+1
      return
      end
c
      subroutine form_num_date(iday,daten)
c***************************************************************
c     forms date number with preceeding zeros
c***************************************************************
c
      character*5 daten
c
      write(daten(1:5),'(i5.5)')iday
c
      return
      end
c
cCCCCCCC  SIMPLE ROUTINES TO ADD NODES - CAN USE INSTEAD OF DRITSCHEL SCHEME CCCCCCCCCCCCCC
c
      subroutine add_nodes
c************************************************************************
c
c      adds a node if d>amu
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
      dimension np_t(nm),mb_t(nm)
c
      l=0
      do i=1,n
         l=l+1
         np_t(i)=1
         mb_t(i)=l-1
         u(l)=x(mb(i)+1)
         v(l)=y(mb(i)+1)
         w(l)=z(mb(i)+1)
         do k=mb(i)+1,mb(i)+np(i)-1
c
            xd=x(k+1)-x(k)
            yd=y(k+1)-y(k)
            zd=z(k+1)-z(k)
            dd=sqrt(xd**2+yd**2+zd**2)
c
            if (dd.gt.amu) then
               l=l+1
               np_t(i)=np_t(i)+1
               u(l)=x(k)+0.5*xd
               v(l)=y(k)+0.5*yd
               w(l)=z(k)+0.5*zd
               l=l+1
               np_t(i)=np_t(i)+1
               u(l)=x(k+1)
               v(l)=y(k+1)
               w(l)=z(k+1)
            else
               l=l+1
               np_t(i)=np_t(i)+1
               u(l)=x(k+1)
               v(l)=y(k+1)
               w(l)=z(k+1)
            end if
         end do
         km=mb(i)+1
         k=mb(i)+np(i)
         xd=x(km)-x(k)
         yd=y(km)-y(k)
         zd=z(km)-z(k)
         dd=sqrt(xd**2+yd**2+zd**2)
c
         if (dd.gt.amu) then
            l=l+1
            np_t(i)=np_t(i)+1
            u(l)=x(k)+0.5*xd
            v(l)=y(k)+0.5*yd
            w(l)=z(k)+0.5*zd
         end if
      end do
c
      do i=1,n
         mb(i)=mb_t(i)
         np(i)=np_t(i)
      end do
      npt=l
      do i=1,npt
         x(i)=u(i)
         y(i)=v(i)
         z(i)=w(i)
      end do
c
      return
      end


cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     SUBROUTINES BELOW ARE AS IN CONTOUR SURGERY CODE OF D.G. DRITSCHEL
c
c
      subroutine adjust
c************************************************************************
c
c      adjusts nodes so that x**2+y**2+z**2=1.
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      do 20 i=1,npt
      f(i)=1./dsqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
      x(i)=f(i)*x(i)
      y(i)=f(i)*y(i)
20    z(i)=f(i)*z(i)
      return
      end
c
c
      subroutine cubic
c************************************************************************
c
c      calculates the interpolation coefficients for every node
c      on every contour.  the interpolation approximately enforces
c      continuity of curvature (except at corners which have
c      effectively infinite curvature).
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      do 10 i=1,npt-1
      dx(i)=x(i+1)-x(i)
      dy(i)=y(i+1)-y(i)
      dz(i)=z(i+1)-z(i)
      ax(i)=y(i)*z(i+1)-y(i+1)*z(i)
      ay(i)=z(i)*x(i+1)-z(i+1)*x(i)
10    az(i)=x(i)*y(i+1)-x(i+1)*y(i)
      do 11 j=1,n
      i1=mb(j)+1
      i2=mb(j)+np(j)
      dx(i2)=x(i1)-x(i2)
      dy(i2)=y(i1)-y(i2)
      dz(i2)=z(i1)-z(i2)
      ax(i2)=y(i2)*z(i1)-y(i1)*z(i2)
      ay(i2)=z(i2)*x(i1)-z(i1)*x(i2)
11    az(i2)=x(i2)*y(i1)-x(i1)*y(i2)
      do 30 i=1,npt
      es(i)=dx(i)**2+dy(i)**2+dz(i)**2+small
      eu(i)=1./es(i)
      e(i)=dsqrt(es(i))
      f(i)=e(i)/dsqrt(ax(i)**2+ay(i)**2+az(i)**2+small)
      ax(i)=f(i)*ax(i)
      ay(i)=f(i)*ay(i)
      az(i)=f(i)*az(i)
      sx(i)=(dy(i)*az(i)-dz(i)*ay(i))/e(i)
      sy(i)=(dz(i)*ax(i)-dx(i)*az(i))/e(i)
30    sz(i)=(dx(i)*ay(i)-dy(i)*ax(i))/e(i)
      do 40 i=2,npt
      a(i)=es(i-1)
      u(i)=x(i-1)-x(i)
      v(i)=y(i-1)-y(i)
40    w(i)=z(i-1)-z(i)
      do 41 j=1,n
      i1=mb(j)+1
      i2=mb(j)+np(j)
      a(i1)=es(i2)
      u(i1)=x(i2)-x(i1)
      v(i1)=y(i2)-y(i1)
41    w(i1)=z(i2)-z(i1)
      do 50 i=1,npt
      f(i)=u(i)*es(i)-dx(i)*a(i)
      g(i)=v(i)*es(i)-dy(i)*a(i)
      h(i)=w(i)*es(i)-dz(i)*a(i)
50    b(i)=(x(i)*(dy(i)*w(i)-v(i)*dz(i))+
     1      y(i)*(dz(i)*u(i)-w(i)*dx(i))+
     2      z(i)*(dx(i)*v(i)-u(i)*dy(i)))/
     3     dsqrt(f(i)*f(i)+g(i)*g(i)+h(i)*h(i)+small3)
c      set curvature to zero at corners
      do 51 i=1,npt
      if (corn(i)) then
      b(i)=0.
      endif
51    continue
      do 60 i=1,npt-1
      c(i)=e(i)*(b(i+1)+b(i))
60    d(i)=e(i)*(b(i+1)-b(i))
      do 61 j=1,n
      i1=mb(j)+1
      i2=mb(j)+np(j)
      c(i2)=e(i2)*(b(i1)+b(i2))
61    d(i2)=e(i2)*(b(i1)-b(i2))
c      calculate the cubic interpolation coefficients:
      do 80 i=1,npt
      es(i)=.5*e(i)
      alp(i)=d(i)/6.-.5*c(i)
      bet(i)=.5*(c(i)-d(i))
80    gam(i)=d(i)/3.
c      cancel the average curvature along the segments adjacent to
c      a corner to reduce node density around corners (see renode):
      do 52 i=1,npt
      if (corn(i)) c(i)=0.
52    continue
      do 53 j=1,n
      do 53 i=1,np(j)
      if (corn(mb(j)+i)) c(mb(j)+mod(i-2+np(j),np(j))+1)=0.
53    continue
      return
      end
c
c
      subroutine prop(k)
c************************************************************************
c
c      gets centroid and maximum radius for contour k
c      (used exclusively to speed up surgery --- see 'dyn' and 'renode')
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      i1=mb(k)+1
      i2=mb(k)+np(k)
      do 10 i=i1,i2-1
      dx(i)=x(i+1)-x(i)
      dy(i)=y(i+1)-y(i)
      dz(i)=z(i+1)-z(i)
      b(i)=dx(i)*dx(i)+dy(i)*dy(i)+dz(i)*dz(i)
10    e(i)=dsqrt(b(i))
      dx(i2)=x(i1)-x(i2)
      dy(i2)=y(i1)-y(i2)
      dz(i2)=z(i1)-z(i2)
      b(i2)=dx(i2)*dx(i2)+dy(i2)*dy(i2)+dz(i2)*dz(i2)
      e(i2)=dsqrt(b(i2))
      return
      end
c
c
      subroutine purge(j)
c************************************************************************
c
c      purges contour j from the set of contours
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      n=n-1
      npto=npt
      npt=mb(j)
      if (j .gt. n) return
      npo=np(j)
      do 10 i=mb(j+1)+1,npto
      u(i-npo)=x(i)
      v(i-npo)=y(i)
      w(i-npo)=z(i)
      sx(i-npo)=dx(i)
      sy(i-npo)=dy(i)
      sz(i-npo)=dz(i)
10    c(i-npo)=b(i)
      do 20 js=j,n
      np(js)=np(js+1)
      mb(js)=mb(js+1)-npo
      ind(js)=ind(js+1)
20    om(js)=om(js+1)
      npt=mb(n)+np(n)
      do 30 i=mb(j)+1,npt
      x(i)=u(i)
      y(i)=v(i)
      z(i)=w(i)
      dx(i)=sx(i)
      dy(i)=sy(i)
      dz(i)=sz(i)
30    b(i)=c(i)
      return
      end
c
c
c
      subroutine renode(irenode)
c************************************************************************
c
c      re-nodes each contour while preserving corner locations.
c      the node density along a contour is given by 1/(dm + amu*d)
c      where 1/d gives a measure of the local and induced curvature.
c      the curvature could be induced by surrounding nodes on the
c      same or other contours that have high curvatures.
c      effects of contours which are sufficiently distant are ignored
c      (see below; see also 'uvw').
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
      logical tcorn(npm),last
c
      do 901 j=1,n
      do 901 i=1,np(j)
      im=mb(j)+i
      ia=mb(j)+mod(i,np(j))+1
      ib=mb(j)+mod(i-2+np(j),np(j))+1
901   corn(im)=(x(ia)-x(im))*(x(ib)-x(im))+(y(ia)-y(im))*(y(ib)-y(im))+
     @    (z(ia)-z(im))*(z(ib)-z(im)) .gt. 0.
c      first, get the updated cubic interpolation coefficients:
      call cubic
c      determine the fractional number of nodes to be placed
c      between each pair of nodes using the node density fn:
      do 70 i=1,npt
      a(i)=0.
      b(i)=0.
      h(i)=dsqrt(e(i)**2+c(i)**2)
      xi(i)=x(i)+.5*dx(i)
      yi(i)=y(i)+.5*dy(i)
70    zi(i)=z(i)+.5*dz(i)
      do 71 j=1,n
      do 71 i=mb(j)+1,mb(j)+np(j)
71    f(i)=dabs(om(j))
c
      if (irenode.eq.0) then
        do 77 j=1,npt
            do 77 i=1,npt
               xx=x(i)-xi(j)
               yy=y(i)-yi(j)
               zz=z(i)-zi(j)
               dd=f(j)/(xx*xx+yy*yy+zz*zz+small)
               a(i)=a(i)+h(j)*dd
 77      b(i)=b(i)+e(j)*dd
         do 79 i=1,npt
 79      g(i)=a(i)/b(i)
      else
        do i=1,npt
           g(i)=h(i)/e(i)
        end do
      end if
c
      do 81 i=1,npt-1
81    h(i)=.5*(g(i)+g(i+1))
      do 82 j=1,n
      i1=mb(j)+1
      i2=mb(j)+np(j)
82    h(i2)=.5*(g(i2)+g(i1))
      do 83 i=1,npt
83    f(i)=(h(i)*ell)**0.66666666666667/(amu*ell)+h(i)
      do 84 i=1,npt
84    o(i)=f(i)*e(i)/(f(i)*emin+1.)
      npt=0
c      this variable tracks the total number of points on all renoded contours
      do 50 j=1,n
c      renode between corners, if they exist.
      l=1
c      this variable tracks the total number of points on a renoded contour
      mbo=mb(j)
      mb(j)=npt
      i1=mbo+1
40    tcorn(npt+l)=corn(i1)
      u(npt+l)=x(i1)
      v(npt+l)=y(i1)
      w(npt+l)=z(i1)
      sum=0.
      i=i1
44    sum=sum+o(i)
      i=i+1
      last=i .gt. mbo+np(j)
      if (last) goto 48
      if (corn(i)) goto 47
      goto 44
47    if (sum .lt. small) then
      i1=i
      goto 40
      else
      i2=i-1
      goto 49
      endif
48    if (sum .lt. small) then
      l=l-1
      goto 39
      else
      i2=i-1
      endif
49    m=nint(sum)+2
c      m-1 is the number of nodes to be placed between i1 and i2.
      fac=m/sum
      do 20 i=i1,i2
20    o(i)=fac*o(i)
c      now, the sum of the o's is equal to m.
c      the first node along a contour (segment) is fixed
c      find the adjusted nodal positions:
      q=0.
      i=i1-1
      do 30 lm=npt+l+1,npt+l+m-1
      if (q .ge. 1.) goto 34
33    q=q+o(i+1)
      i=i+1
      if (q .lt. 1.) goto 33
34    q=q-1.
      p=1.-q/o(i)
      eta=p*(alp(i)+p*(bet(i)+p*gam(i)))
      del=es(i)*p*(1.-p)
      tcorn(lm)=.false.
      u(lm)=x(i)+p*dx(i)+eta*ax(i)+del*sx(i)
      v(lm)=y(i)+p*dy(i)+eta*ay(i)+del*sy(i)
30    w(lm)=z(i)+p*dz(i)+eta*az(i)+del*sz(i)
      if (last) then
      l=l+m-1
      goto 39
      else
      l=l+m
      i1=i2+1
      goto 40
      endif
39    np(j)=l
      npt=npt+l
50    continue
c      switch arrays around again:
      do 60 i=1,npt
      corn(i)=tcorn(i)
      x(i)=u(i)
      y(i)=v(i)
60    z(i)=w(i)
      return
      end
c
c
c
      subroutine surgery
c************************************************************************
c
c     performs contour surgery (see dritschel(1988,1989) for details).
c
c************************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
      logical logic(npm),keep1,keep2
c
      do 70 j=1,n
70    call prop(j)
c
ccc    ***vortex splitting and merger***
c*     consider first the case when a single contour is involved.
c      in this case, the contour could break into two, either in the
c      form of two separate vortices or an annular vortex.  both cases
c      are treated the same.
      j=0
100   j=j+1
c      contour "j"
      if (j .gt. n) goto 899
      do 120 i=mb(j)+1,mb(j)+np(j)-1
      do 119 is=i+1,mb(j)+np(j)
c      we will see if the node i gets within a distance dm to the
c      line segment between is and is+1.
      ax(is)=x(is)-x(i)
      ay(is)=y(is)-y(i)
      az(is)=z(is)-z(i)
      a(is)=ax(is)*dx(is)+ay(is)*dy(is)+az(is)*dz(is)
119   logic(is)=a(is)*(a(is)+b(is))+small .gt. 0.
      do 120 is=i+1,mb(j)+np(j)
      if (logic(is)) goto 120
c      passing this condition, a perpendicular can be dropped onto the
c      line segment from is to is+1.  check distance to line segment:
c      dis=dabs((ax(is)*dy(is)-ay(is)*dx(is)))/e(is)
      dis=dsqrt(((ay(is)*dz(is)-dy(is)*az(is))**2
     @         +(az(is)*dx(is)-dz(is)*ax(is))**2
     @         +(ax(is)*dy(is)-dx(is)*ay(is))**2)/b(is))
      if (dis*(dis-dm) .lt. 0.) goto 121
120   continue
      goto 199
c      now, break the contour into two at node i.  we temporarily
c      store the two pieces in the (xi,yi,zi) and (xf,yf,zf) arrays, get
c      rid of contour j, and then introduce two new contours.
c      first piece:
121   p=-a(is)/b(is)
      x(i)=x(i)+.5*(ax(is)+p*dx(is))
      y(i)=y(i)+.5*(ay(is)+p*dy(is))
      z(i)=z(i)+.5*(az(is)+p*dz(is))
      fac=1./dsqrt(x(i)**2+y(i)**2+z(i)**2)
      x(i)=fac*x(i)
      y(i)=fac*y(i)
      z(i)=fac*z(i)
      np1=is-i+1
      keep1=np1 .gt. 4
      if (keep1) then
      om1=om(j)
      ind1=ind(j)
      do 135 ih=1,np1
      xi(ih)=x(i+ih-1)
      yi(ih)=y(i+ih-1)
135   zi(ih)=z(i+ih-1)
      endif
c      second piece:
      np2=np(j)-is+i
      keep2=np2 .gt. 4
      if (keep2) then
      om2=om(j)
      ind2=ind(j)
      xf(1)=x(i)
      yf(1)=y(i)
      zf(1)=z(i)
      ii=is-mb(j)-2
      do 140 ih=2,np2
      im=mb(j)+mod(ii+ih,np(j))+1
      xf(ih)=x(im)
      yf(ih)=y(im)
140   zf(ih)=z(im)
      endif
c      get rid of contour j
      call purge(j)
c      reassign temporary arrays:
c--    1st contour:
      if (keep1) then
      n=n+1
      np(n)=np1
      om(n)=om1
      ind(n)=ind1
      mb(n)=npt
      do 180 ih=1,np1
      x(mb(n)+ih)=xi(ih)
      y(mb(n)+ih)=yi(ih)
180   z(mb(n)+ih)=zi(ih)
      jj=n
      call prop(jj)
      npt=npt+np1
      endif
c--    2nd contour:
      if (keep2) then
      n=n+1
      np(n)=np2
      om(n)=om2
      ind(n)=ind2
      mb(n)=npt
      do 195 ih=1,np2
      x(mb(n)+ih)=xf(ih)
      y(mb(n)+ih)=yf(ih)
195   z(mb(n)+ih)=zf(ih)
      jj=n
      call prop(jj)
      npt=npt+np2
      endif
      j=j-1
      goto 100
c
c*     next consider the merger of two distinct vortices (or the
c      breakup of an annular vortex).   this process replaces two
c      contours with one.
199   js=j
201   js=js+1
      if (js .gt. n) goto 100
      if (ind(j) .ne. ind(js)) goto 201
      if (om(j) .ne. om(js)) goto 201
c      these are the mergeability conditions.
      do 220 i=mb(j)+1,mb(j)+np(j)
      do 219 is=mb(js)+1,mb(js)+np(js)
c      we will see if the node i gets within a distance dm to the
c      line segment between is and is+1.
      ax(is)=x(is)-x(i)
      ay(is)=y(is)-y(i)
      az(is)=z(is)-z(i)
      a(is)=ax(is)*dx(is)+ay(is)*dy(is)+az(is)*dz(is)
219   logic(is)=a(is)*(a(is)+b(is))+small .gt. 0.
      do 220 is=mb(js)+1,mb(js)+np(js)
      if (logic(is)) goto 220
c      passing this condition, a perpendicular can be dropped onto the
c      line segment from is to is+1.  check distance to line segment:
c      dis=dabs((ax(is)*dy(is)-ay(is)*dx(is)))/e(is)
      dis=dsqrt(((ay(is)*dz(is)-dy(is)*az(is))**2
     @         +(az(is)*dx(is)-dz(is)*ax(is))**2
     @         +(ax(is)*dy(is)-dx(is)*ay(is))**2)/b(is))
      if (dis*(dis-dm) .lt. 0.) goto 221
220   continue
      goto 201
c      now we paste two contours together at node i:
221   p=-a(is)/b(is)
      x(i)=x(i)+.5*(ax(is)+p*dx(is))
      y(i)=y(i)+.5*(ay(is)+p*dy(is))
      z(i)=z(i)+.5*(az(is)+p*dz(is))
      fac=1./dsqrt(x(i)**2+y(i)**2+z(i)**2)
      x(i)=fac*x(i)
      y(i)=fac*y(i)
      z(i)=fac*z(i)
c      we will use the xi,yi,zi arrays as temporary storage of the
c      combined contour, get rid of contours j and js, and finally
c      introduce a new contour, the combined contour:
      ii=i-mb(j)-2
      do 245 ih=1,np(j)+1
      im=mb(j)+mod(ii+ih,np(j))+1
      xi(ih)=x(im)
      yi(ih)=y(im)
245   zi(ih)=z(im)
      ii=is-mb(js)-1
      do 250 ih=1,np(js)
      im=mb(js)+mod(ii+ih,np(js))+1
      xi(np(j)+1+ih)=x(im)
      yi(np(j)+1+ih)=y(im)
250   zi(np(j)+1+ih)=z(im)
      np1=np(j)+np(js)+1
      om1=om(j)
      ind1=ind(j)
      call purge(js)
      call purge(j)
      n=n+1
      om(n)=om1
      ind(n)=ind1
      np(n)=np1
      mb(n)=npt
      do 265 ih=1,np1
      x(npt+ih)=xi(ih)
      y(npt+ih)=yi(ih)
265   z(npt+ih)=zi(ih)
      jj=n
      call prop(jj)
      npt=npt+np1
      j=j-1
      goto 100
c
899   continue
c
      return
      end
c

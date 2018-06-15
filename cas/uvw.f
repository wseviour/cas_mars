c
      subroutine r_grid1(tt)
c**********************************************************************
c
c     NMC ANALYSES (2x5), N.H. data
c
c     Reads in velocity grid
c     If tt=-99
c              read in header information:
c              dlat, dlong (grid spacing)
c              dt, ns (time between velocity data = ns*dt
c     otherwise
c              read in velocity grid:
c              [ u(i,j) , v(i,j) ]=velocity at
c              lat=(i-1)*dlat, i=1,..,46 : dlat=2
c              long=(j-1)*dlong, j=1,..,144 : dlong=5
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      common /trig1/ clong(289,289),slong(289,289)
      common /trig2/ clat(289,289),slat(289,289)
      dimension uin(144,92),vin(144,92)
      dimension uu(289,289),vv(289,289)
c
c     scaling factor for velocities (t=1 equals 1 day)
c      Earth value:
c      fac=2.0*pi/464.57
c      Mars value
       fac = 2.0*pi/240.36

c
      if (tt.eq.-99.0) then
c
c     form latitude - longitude grid
c     ------------------------------
         dlat=1.
         dlong=2.5
         nlat=92
         nlong=145
         do 74 i=1,nlat
            do 75 j=1,nlong
               xlat(i,j)=(i-1)*dlat
               xlong(i,j)=(j-1)*dlong
               slong(i,j)=dsin(xlong(i,j)*pifac)
               clong(i,j)=dcos(xlong(i,j)*pifac)
               clat(i,j)=dcos(xlat(i,j)*pifac)
               slat(i,j)=dsin(xlat(i,j)*pifac)
 75         continue
 74      continue
      else
c
c     read velocities (spherical components)
c     --------------------------------------
         do j=1,92
            read(31,*)(uin(i,j),i=1,144)
         end do
         do j=1,92
            read(32,*)(vin(i,j),i=1,144)
         end do
c
         do 90 j=1,nlat
            do 91 i=1,nlong
               uu(j,i)=fac*uin(i,nlat+1-j)
               vv(j,i)=fac*vin(i,nlat+1-j)
 91         continue
 90      continue
         do 92 j=1,nlat
            uu(j,nlong)=uu(j,1)
            vv(j,nlong)=vv(j,1)
 92      continue
c
c
c     convert to velocity cartesian components
c     ----------------------------------------
         do i=1,nlat
            do j=1,nlong
               u1(i,j)=-slong(i,j)*uu(i,j)
     &                 -clong(i,j)*slat(i,j)*vv(i,j)
               v1(i,j)= clong(i,j)*uu(i,j)
     &                 -slong(i,j)*slat(i,j)*vv(i,j)
               w1(i,j)= clat(i,j)*vv(i,j)
            end do
         end do
c
c     velocity at pole
c     ----------------
         usum=0.
         vsum=0.
         do j=1,nlong-1
            usum=usum+u1(nlat-1,j)
            vsum=vsum+v1(nlat-1,j)
         end do
         usum=usum/float(nlong-1)
         vsum=vsum/float(nlong-1)
         do j=1,nlong
            u1(nlat,j)=usum
            v1(nlat,j)=vsum
         end do
c
      end if
      tt=0
      return
c
 1000 tt=-999
      return
      end
c
c
      subroutine uvw1
c**********************************************************************
c
c     NMC ANALYSES (2x5), N.H. data
c
c     calculates the velocity at each node by interpolating from
c     the input velocity grid.
c     linear interpolation in space and time.
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      dimension uu(289,289),vv(289,289),ww(289,289)
c
c     linear interpolation in time: t=t0+s*deltat
c     -------------------------------------------
      s = (t - t0) / deltat
      do 70 i=1,nlat
         do 80 j=1,nlong
            uu(i,j)=(1.0-s)*u0(i,j) + s*u1(i,j)
            vv(i,j)=(1.0-s)*v0(i,j) + s*v1(i,j)
            ww(i,j)=(1.0-s)*w0(i,j) + s*w1(i,j)
 80      continue
 70   continue
c
      do 60 k=1,npt
c
c     calculate latitude,longitude of k'th node
c     -----------------------------------------
         if (z(k).gt.1.0)  z(k)=1.0
         if (z(k).lt.-1.0) z(k)=-1.0
         zlat=asin(z(k))/pifac
         zlong=atan2(y(k),x(k))/pifac
         if (zlong.lt.0.) zlong=360.+zlong
c
c     calculate velocity at k'th node
c     -------------------------------
         if (zlat.le.0) then
c     no motion in SH
            u(k)=0.
            v(k)=0.
            w(k)=0.
         else
            i=int( (zlat)/dlat ) + 1
            j=int(zlong/dlong) + 1
            p=(zlat-xlat(i,1))/dlat
            q=(zlong-xlong(1,j))/dlong
            ps=1.-p
            qs=1.-q
c
            u(k) = ps*qs*uu(i,j) + p*qs*uu(i+1,j) + q*ps*uu(i,j+1)
     @          + p*q*uu(i+1,j+1)
            v(k) = ps*qs*vv(i,j) + p*qs*vv(i+1,j) + q*ps*vv(i,j+1)
     @          + p*q*vv(i+1,j+1)
            w(k) = ps*qs*ww(i,j) + p*qs*ww(i+1,j) + q*ps*ww(i,j+1)
     @          + p*q*ww(i+1,j+1)
         end if
c
 60   continue
c
      return
      end
c
c
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c
      subroutine r_grid3(tt)
c**********************************************************************
c
c     NMC ANALYSES (2x5), Global data
c
c     Reads in velocity grid
c     If tt=-99
c              read in header information:
c              dlat, dlong (grid spacing)
c              dt, ns (time between velocity data = ns*dt
c     otherwise
c              read in velocity grid:
c              [ u(i,j) , v(i,j) ]=velocity at
c              lat=-90.+(i-1)*dlat, i=1,..,91 : dlat=2
c              long=(j-1)*dlong, j=1,..,144 : dlong=5
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      common /trig1/ clong(289,289),slong(289,289)
      common /trig2/ clat(289,289),slat(289,289)
      dimension uin(144,91),vin(144,91)
      dimension uu(289,289),vv(289,289)
c
c     scaling factor for velocities (t=1 equals 1 day)
c
c      fac=2.0*pi/464.57
      fac = 2.0*pi/240.36
c
      if (tt.eq.-99.0) then
c
c     form latitude - longitude grid
c     ------------------------------
         dlat=2.
         dlong=5.
         nlat=91
         nlong=145
         do 74 i=1,nlat
            do 75 j=1,nlong
               xlat(i,j)=-90.+(i-1)*dlat
               xlong(i,j)=(j-1)*dlong
               slong(i,j)=dsin(xlong(i,j)*pifac)
               clong(i,j)=dcos(xlong(i,j)*pifac)
               clat(i,j)=dcos(xlat(i,j)*pifac)
               slat(i,j)=dsin(xlat(i,j)*pifac)
 75         continue
 74      continue
      else
c
c     read velocities (spherical components)
c     --------------------------------------
         do j=1,91
            read(31,*)(uin(i,j),i=1,144)
         end do
         do j=1,91
            read(32,*)(vin(i,j),i=1,144)
         end do
c
         do 90 j=1,nlat
            do 91 i=1,nlong
               uu(j,i)=fac*uin(i,j)
               vv(j,i)=fac*vin(i,j)
 91         continue
 90      continue
         do 92 j=1,nlat
            uu(j,nlong)=uu(j,1)
            vv(j,nlong)=vv(j,1)
 92      continue
c
c
c     convert to velocity cartesian components
c     ----------------------------------------
         do i=1,nlat
            do j=1,nlong
               u1(i,j)=-slong(i,j)*uu(i,j)
     &                 -clong(i,j)*slat(i,j)*vv(i,j)
               v1(i,j)= clong(i,j)*uu(i,j)
     &                 -slong(i,j)*slat(i,j)*vv(i,j)
               w1(i,j)= clat(i,j)*vv(i,j)
            end do
         end do
c
c     velocity at poles
c     -----------------
         usum=0.
         vsum=0.
         do j=1,nlong-1
            usum=usum+u1(2,j)
            vsum=vsum+v1(2,j)
         end do
         usum=usum/float(nlong-1)
         vsum=vsum/float(nlong-1)
         do j=1,nlong
            u1(1,j)=usum
            v1(1,j)=vsum
         end do
c
         usum=0.
         vsum=0.
         do j=1,nlong-1
            usum=usum+u1(nlat-1,j)
            vsum=vsum+v1(nlat-1,j)
         end do
         usum=usum/float(nlong-1)
         vsum=vsum/float(nlong-1)
         do j=1,nlong
            u1(nlat,j)=usum
            v1(nlat,j)=vsum
         end do
c
      end if
      tt=0
      return
c
 1000 tt=-999
      return
      end
c
c
      subroutine uvw3
c**********************************************************************
c
c     NMC Global data
c
c     calculates the velocity at each node by interpolating from
c     the input velocity grid.
c     linear interpolation in space and time.
c
c**********************************************************************
      implicit real*8(a-h,o-z)
      include 'cas.h'
c
      dimension uu(289,289),vv(289,289),ww(289,289)
c
c     linear interpolation in time: t=t0+s*deltat
c     -------------------------------------------
      s = (t - t0) / deltat
      do 70 i=1,nlat
         do 80 j=1,nlong
            uu(i,j)=(1.0-s)*u0(i,j) + s*u1(i,j)
            vv(i,j)=(1.0-s)*v0(i,j) + s*v1(i,j)
            ww(i,j)=(1.0-s)*w0(i,j) + s*w1(i,j)
 80      continue
 70   continue
c
      do 60 k=1,npt
c
c     calculate latitude,longitude of k'th node
c     -----------------------------------------
         if (z(k).gt.1.0)  z(k)=1.0
         if (z(k).lt.-1.0) z(k)=-1.0
         zlat=asin(z(k))/pifac
         zlong=atan2(y(k),x(k))/pifac
         if (zlong.lt.0.) zlong=360.+zlong
c
c     calculate velocity at k'th node
c     -------------------------------
         i=int( (zlat+90.)/dlat ) + 1
         j=int(zlong/dlong) + 1
         p=(zlat-xlat(i,1))/dlat
         q=(zlong-xlong(1,j))/dlong
         ps=1.-p
         qs=1.-q
c
         u(k) = ps*qs*uu(i,j) + p*qs*uu(i+1,j) + q*ps*uu(i,j+1)
     @        + p*q*uu(i+1,j+1)
         v(k) = ps*qs*vv(i,j) + p*qs*vv(i+1,j) + q*ps*vv(i,j+1)
     @        + p*q*vv(i+1,j+1)
         w(k) = ps*qs*ww(i,j) + p*qs*ww(i+1,j) + q*ps*ww(i,j+1)
     @        + p*q*ww(i+1,j+1)
c
 60   continue
c
      return
      end
c

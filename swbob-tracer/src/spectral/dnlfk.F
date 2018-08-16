      subroutine dnlfk (m,n,cp)

c-------------------------------------------------------------
c     computes the  "a_m,n,k" coefs in
c     P Swarztrauber, The Vector Harmonic Transform Method for Solving 
c     Partial Differential Equations in Spherical Geometry
c     MWR, vol 121, no 12, pp 3416-3437
c     using eqs (4.5), (4.6), (4.9), (4.10).
c-------------------------------------------------------------

#include <type.h>
c
      REAL_TYPE cp(*),fnum,fden,fnmh,a1,b1,c1,cp2,fnnp1,fnmsq,fk,
     1          t1,t2,pm1,sc10,sc20,sc40
      parameter (sc10=1024.d0)
      parameter (sc20=sc10*sc10)
      parameter (sc40=sc20*sc20)
c
      cp(1) = 0.0d0
      ma = iabs(m)
      if(ma .gt. n) return
      if(n-1) 2,3,5
    2 cp(1) = dsqrt(2.d0)
      return
    3 if(ma .ne. 0) go to 4
      cp(1) = dsqrt(1.5d0)
      return
    4 cp(1) = dsqrt(.75d0)
      if(m .eq. -1) cp(1) = -cp(1)
      return
    5 if(mod(n+ma,2) .ne. 0) go to 10
      nmms2 = (n-ma)/2
      fnum = n+ma+1
      fnmh = n-ma+1
      pm1 = 1.d0
      go to 15
   10 nmms2 = (n-ma-1)/2
      fnum = n+ma+2
      fnmh = n-ma+2
      pm1 = -1.d0
c      t1 = 1.
c      t1 = 2.d0**(n-1)
c      t1 = 1.d0/t1
 15   t1 = 1.d0/sc20
      nex = 20
      fden = 2.d0
      if(nmms2 .lt. 1) go to 20
      do 18 i=1,nmms2
      t1 = fnum*t1/fden
      if(t1 .gt. sc20) then
      t1 = t1/sc40
      nex = nex+40
      end if
      fnum = fnum+2.0d0
      fden = fden+2.0d0
   18 continue
   20 t1 = t1/2.d0**(n-1-nex)
      if(mod(ma/2,2) .ne. 0) t1 = -t1
      t2 = 1.0d0
      if(ma .eq. 0) go to 26
      do 25 i=1,ma
      t2 = fnmh*t2/(fnmh+pm1)
      fnmh = fnmh+2.0d0
   25 continue
   26 cp2 = t1*dsqrt((n+.5d0)*t2)
      fnnp1 = n*(n+1)
      fnmsq = fnnp1-2.d0*ma*ma
      l = (n+1)/2
      if(mod(n,2) .eq. 0 .and. mod(ma,2) .eq. 0) l = l+1
      cp(l) = cp2
      if(m .ge. 0) go to 29
      if(mod(ma,2) .ne. 0) cp(l) = -cp(l)
   29 if(l .le. 1) return
      fk = n
      a1 = (fk-2.0d0)*(fk-1.0d0)-fnnp1
      b1 = 2.0d0*(fk*fk-fnmsq)
      cp(l-1) = b1*cp(l)/a1
   30 l = l-1
      if(l .le. 1) return
      fk = fk-2.0d0
      a1 = (fk-2.0d0)*(fk-1.0d0)-fnnp1
      b1 = -2.0d0*(fk*fk-fnmsq)
      c1 = (fk+1.)*(fk+2.0d0)-fnnp1
      cp(l-1) = -(b1*cp(l)+c1*cp(l+1))/a1
      go to 30
      end


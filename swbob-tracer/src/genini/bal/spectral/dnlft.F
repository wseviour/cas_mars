      subroutine dnlft (m,n,theta,cp,pb)
#include <type.h>

      REAL_TYPE cp(*)
      REAL_TYPE pb,theta,cdt,sdt,cth,sth,chh
      cdt = dcos(theta+theta)
      sdt = dsin(theta+theta)
      nmod=mod(n,2)
      mmod=mod(abs(m),2)
      if(nmod)1,1,2
1     if(mmod)3,3,4
c
c     n even, m even
c
3     kdo=n/2
      pb = 0.50d0*cp(1)
      if(n .eq. 0) return
      cth = cdt
      sth = sdt
      do 170 k=1,kdo
c     pb = pb+cp(k+1)*dcos(2*k*theta)
      pb = pb+cp(k+1)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  170 continue
      return
c
c     n even, m odd
c
    4 kdo = n/2
      pb = 0.0d0
      cth = cdt
      sth = sdt
      do 180 k=1,kdo
c     pb = pb+cp(k)*dsin(2*k*theta)
      pb = pb+cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  180 continue
      return
2     if(mmod)13,13,14
c
c     n odd, m even
c
13    kdo = (n+1)/2
      pb = 0.0d0
      cth = dcos(theta)
      sth = dsin(theta)
      do 190 k=1,kdo
c     pb = pb+cp(k)*dcos((2*k-1)*theta)
      pb = pb+cp(k)*cth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  190 continue
      return
c
c     n odd, m odd
c
   14 kdo = (n+1)/2
      pb = 0.0d0
      cth = dcos(theta)
      sth = dsin(theta)
      do 200 k=1,kdo
c     pb = pb+cp(k)*dsin((2*k-1)*theta)
      pb = pb+cp(k)*sth
      chh = cdt*cth-sdt*sth
      sth = sdt*cth+cdt*sth
      cth = chh
  200 continue
      return
      end

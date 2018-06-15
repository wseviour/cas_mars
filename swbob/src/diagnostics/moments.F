      subroutine moments(phy,
     $                   da,
     $                   lat,
     $                   lon)
      implicit none

#include <coord.h>
#include <modeltime.h>
#include <stptrs.h>
#include <stats.h>

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE da(nlat2)
      REAL_TYPE lat(nlat2)
      REAL_TYPE lon(nlon)

      REAL_TYPE phi
      REAL_TYPE phi0(nlon,nlat2,0:1)
      REAL_TYPE phit(nlon,nlat2,0:1)

      REAL_TYPE cslat(nlat2)
      REAL_TYPE snlat(nlat2)

      INT_TYPE i,j,k,ins

      do j=1,nlat2
         cslat(j)=DCOS(lat(j))
         snlat(j)=DSIN(lat(j))
      end do

      call initphi(phit,
     $             time_current,
     $             Alpha_coord,
     $             cslat,
     $             snlat,
     $             lon)

      call initphi(phi0,
     $             0.0D0,
     $             Alpha_coord,
     $             cslat,
     $             snlat,
     $             lon)

      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               do i=1,nlon
                  phi=phy(i,j,phi_p,ins,k)
                  dphisq(k)=dphisq(k)+(phi-phiave(k))**2
                  dphitsq(k)=dphitsq(k)+(phit(i,j,ins)-phitave(k))**2
                  dphi0sq(k)=dphi0sq(k)+(phi0(i,j,ins)-phi0ave(k))**2
               end do
            end do
         end do
      end do


      end



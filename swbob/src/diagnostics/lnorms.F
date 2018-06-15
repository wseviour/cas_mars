      subroutine lnorms(phy,
     $                  da,
     $                  lat,
     $                  lon,
     $                  jeg)
      implicit none

#include <coord.h>
#include <modeltime.h>
#include <stptrs.h>
#include <stats.h>
#include <io.h>
#include <planet.h>
c      REAL_TYPE Gravity (parameter)

      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE da(nlat2)
      REAL_TYPE lat(nlat2)
      REAL_TYPE lon(nlon)
      INT_TYPE  jeg

      REAL_TYPE phi
      REAL_TYPE phi0(nlon,nlat2,0:1)
      REAL_TYPE phit(nlon,nlat2,0:1)

      INT_TYPE i,j,k,ins
      REAL_TYPE cslat(nlat2)
      REAL_TYPE snlat(nlat2)

c----------------------------------------
c     Executable Statements
c----------------------------------------

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
                  phi=phy(i,j,phi_p,ins,k)/Gravity

c Numerator of l1 norm
                  l1num(k)  = l1num(k)+DABS(phi-phit(i,j,ins)/Gravity)*da(j)
c Denominator of l1 norm                  
                  l1den(k)  = l1den(k)+DABS(phit(i,j,ins)/Gravity)*da(j)

                  l2num(k)  = l2num(k)+((phi-phit(i,j,ins)/Gravity)**2)*da(j)
                  l2den(k)  = l2den(k)+((phit(i,j,ins)/Gravity)**2)*da(j)
                  linfnum(k)= DMAX1(DABS(phi-phit(i,j,ins)/Gravity),linfnum(k))
                  linfden(k)= DMAX1(DABS(phit(i,j,ins)/Gravity),linfden(k))

                  phiave(k) = phiave(k)+phi*da(j)
                  phitave(k)= phitave(k)+phit(i,j,ins)*da(j)
                  phi0ave(k)= phi0ave(k)+phi0(i,j,ins)*da(j)

#if 0
                  if ((jeg.eq.jelemg).and.(ins.eq.0).and.(k.eq.1).and.j.eq.nlat2) then
                     write(MPI_UNIT,10)lon(i),lat(j),phi,phit(i,j,ins),ABS(phi-phit(i,j,ins))*da(j)
 10                  format(f14.8,1x,f14.8,2x,f20.14,2x,f20.14,2x,f20.14)
                  end if
#endif
               end do
            end do
         end do
      end do


      end


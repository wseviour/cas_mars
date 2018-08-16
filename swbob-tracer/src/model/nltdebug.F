      subroutine nltdebug(nltbug)
      implicit none

#include <dims.h>
#include <decomp.h>
#include <mpcluster.h>
#include <stptrs.h>
#include <nltptrs.h>
#include <planet.h>
#include <physical.h>
#include <io.h>

      REAL_TYPE nltbug(nlon,nlat2*nnlt*2*plev,ielemd,jelemd)               

      REAL_TYPE trma(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE trmb(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE trmc(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE trmd(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE trme(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE trmf(nlon,nlat2,0:1,plev,ielemd,jelemd)

      REAL_TYPE zmax  ! compute zonal maximum
      external  zmax

      REAL_TYPE amax,bmax,cmax,dmax,emax
      INT_TYPE  je,jeg,j,ins

      ! extract fields

      call nltunpack(nltbug,trma,a_p)
      call nltunpack(nltbug,trmb,b_p)
      call nltunpack(nltbug,trmc,c_p)
      call nltunpack(nltbug,trmd,d_p)
      call nltunpack(nltbug,trme,e_p)

      write(*,*)
      write(*,*)"In nltdebug: NLT VALUES:"

      ins=0  ! north
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do j=1,nlat2
            amax=zmax(trma,j,je,jeg,ins)
            bmax=zmax(trmb,j,je,jeg,ins)
            cmax=zmax(trmc,j,je,jeg,ins)
            dmax=zmax(trmd,j,je,jeg,ins)
            emax=zmax(trme,j,je,jeg,ins)
            write(*,*)lat(j,jeg),bmax
         end do
      end do

      ins=1 ! south
      do jeg=jeend(rank),jebeg(rank),-1
         je=jeg-jebeg(rank)+1
         do j=nlat2,1,-1
            amax=zmax(trma,j,je,jeg,ins)
            bmax=zmax(trmb,j,je,jeg,ins)
            cmax=zmax(trmc,j,je,jeg,ins)
            dmax=zmax(trmd,j,je,jeg,ins)
            emax=zmax(trme,j,je,jeg,ins)
            write(*,*)-lat(j,jeg),bmax
         end do
      end do
      write(MPI_UNIT,3)

      end















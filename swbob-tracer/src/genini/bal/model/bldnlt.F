      subroutine bldnlt(phy,
     $                  nlt,
     $                  phitopo,
     $                  lat,
     $                  lon)
      implicit none

#include <dims.h>
#include <stptrs.h>
#include <nltptrs.h>
#include <constants.h>
#include <planet.h>


      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE nlt(nlon,nlat2,nnlt,0:1,plev)
      REAL_TYPE phitopo(nlon,nlat2,0:1,plev)
      REAL_TYPE lat(nlat2)
      REAL_TYPE lon(nlon)

      INT_TYPE i,j,k,ins
      REAL_TYPE cslat(nlat2),snlat(nlat2)
      REAL_TYPE ramp, latsign

c----------------------------------------
c     Executable Statements
c---------------------------------------- 

c      ramp=1.d0
      ramp=0.d0

      do j=1,nlat2
         cslat(j)=DCOS(lat(j))
         snlat(j)=DSIN(lat(j))
      end do

      do k=1,plev
         do ins=0,1
            latsign=(1.d0-2.d0*real(ins))
            do j=1,nlat2
               do i=1,nlon
                  nlt(i,j,a_p,ins,k)=phy(i,j,u_p,ins,k)*phy(i,j,vor_p,ins,k)

                  nlt(i,j,b_p,ins,k)=phy(i,j,v_p,ins,k)*phy(i,j,vor_p,ins,k)

                  nlt(i,j,c_p,ins,k)=phy(i,j,u_p,ins,k)*(phy(i,j,phi_p,ins,k)+ramp*phitopo(i,j,ins,k))

                  nlt(i,j,d_p,ins,k)=phy(i,j,v_p,ins,k)*(phy(i,j,phi_p,ins,k)+ramp*phitopo(i,j,ins,k))

                  nlt(i,j,e_p,ins,k)=(phy(i,j,u_p,ins,k)*phy(i,j,u_p,ins,k)+
     $                                phy(i,j,v_p,ins,k)*phy(i,j,v_p,ins,k))
     $                                /(2.0d0*cslat(j)*cslat(j))

               end do
            end do
         end do
      end do

      end 

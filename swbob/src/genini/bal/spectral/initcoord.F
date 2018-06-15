#include <type.h>
      subroutine initcoord(colatitude,gausswt,plonmax,Radius)
      implicit none

!--------------------------------------------------------------------
!      Include Files
#include <constants.h>
#include <gaussian.h>
#include <physical.h>
#include <nltptrs.h>
#include <decomp.h>
#include <mpcluster.h>
!--------------------------------------------------------------------

      REAL_TYPE colatitude(plat/2)    ! colatitudes 
      REAL_TYPE gausswt(plat/2)       ! gaussian weights 
      INT_TYPE  plonmax(jelemg)       ! max longitude for each latitude
      REAL_TYPE Radius                ! Planetary Radius

      INT_TYPE  j,je,jeg,jptr
      INT_TYPE  i,ie,iptr

      REAL_TYPE dx

      ! compute latitude coordinate, gaussian wt factors.

      jptr=1
      do jeg=1,jelemg
         ielem(jeg)    = plonmax(jeg)/nlon
         do j=1,nlat2
            lat(j,jeg)    = (Pi/2.0d0) - colatitude(jptr)
            cslat(j,jeg)  = DCOS(lat(j,jeg))
            snlat(j,jeg)  = DSIN(lat(j,jeg))
            da(j,jeg)     = gausswt(jptr)/(2.0d0*plonmax(jeg))
            gen_gw(j,jeg) = gausswt(jptr)
            gw(j,a_p,jeg) = gausswt(jptr)/(Radius*cslat(j,jeg)*cslat(j,jeg))
            gw(j,b_p,jeg) = gausswt(jptr)/(Radius*cslat(j,jeg)*cslat(j,jeg))
            gw(j,c_p,jeg) = gausswt(jptr)/(Radius*cslat(j,jeg)*cslat(j,jeg))
            gw(j,d_p,jeg) = gausswt(jptr)/(Radius*cslat(j,jeg)*cslat(j,jeg))
            gw(j,e_p,jeg) = gausswt(jptr)
            jptr=jptr+1
         end do
      end do

      ! compute nodal longitude coordinate (assume equally spaced)

      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         dx=(2.0d0*Pi)/plonmax(jeg)       ! tile longitude spacing
         iptr=1
         do ie=1,ielem(jeg)
            do i=1,nlon
               lon(i,ie,je)=(iptr-1)*dx
               iptr=iptr+1
            end do
         end do
      end do


      end








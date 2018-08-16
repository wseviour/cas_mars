
      subroutine initcoord(colatitude,gausswt,plonmax,Radius)

c--------------------------------------------------------------------
c     Description: computes divers lat and weight related values:
c     lat, coslat, sinlat, da, gen_gw, gw, lon 
c
c     Input: colatitude, gausswt,plonmax, Radius
c--------------------------------------------------------------------
      implicit none
c--------------------------------------------------------------------
c      Include Files
c--------------------------------------------------------------------
#include <type.h>
#include <constants.h>
c     REAL_TYPE   Pi
c     parameter   (Pi=3.141592653589793D0) ! mathematical constant
#include <gaussian.h>
c     REAL_TYPE gen_gw(nlat2,jelemg) ! distributed gaussian weights used in gen_analysis
c     REAL_TYPE gw(nlat2,nnlt,jelemg) ! distributed gaussian weights associated
c                                     ! to the nonlinear terms
c     common /gaussian/gw,gen_gw
#include <physical.h>
c     REAL_TYPE lon(nlon,ielemd,jelemd)  ! nodal longitudes in rads
c     REAL_TYPE lat  (nlat2,jelemg) ! latitudes of gaussian grids
c     REAL_TYPE cslat(nlat2,jelemg) ! cosine of latitude
c     REAL_TYPE snlat(nlat2,jelemg) ! sine of latitude
c     REAL_TYPE da   (nlat2,jelemg) ! area weight
c     common /physical/ ielem,lat,cslat,snlat,da,lon
#include <nltptrs.h>
#include <decomp.h>
c      INT_TYPE jebeg(0:P_NODE-1),jeend(0:P_NODE-1) ! bounds of global latitudinal index for a tile element on a given PE 
c      common /jdec/ jbeg,jend,jebeg,jeend
#include <mpcluster.h>

      REAL_TYPE colatitude(plat/2)    ! colatitudes 
      REAL_TYPE gausswt(plat/2)       ! gaussian weights 
      INT_TYPE  plonmax(jelemg)       ! longitude spacing
      REAL_TYPE Radius                ! Planetary Radius

      INT_TYPE  j,je,jeg,i,ie,iptr
      INT_TYPE  jptr ! index running over the global latitudes in gausswt(plat/2)

      REAL_TYPE dx

c-------------------------------------------------------------
c Begin executable statements
c-------------------------------------------------------------

c-------------------------------------------------------------
c     Compute latitude coordinate, gaussian wt factors.
c-------------------------------------------------------------
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
            gw(j,f_p,jeg) = gausswt(jptr)
            gw(j,g_p,jeg) = gausswt(jptr)
c     rks - weights for tracer terms:
            gw(j,m_p,jeg) = gausswt(jptr)/(Radius*cslat(j,jeg)*cslat(j,jeg))
            gw(j,n_p,jeg) = gausswt(jptr)/(Radius*cslat(j,jeg)*cslat(j,jeg))
            gw(j,l_p,jeg) = gausswt(jptr)
            jptr=jptr+1
         end do
      end do

c-------------------------------------------------------------
c     Compute nodal longitude coordinate (assume equally spaced)
c-------------------------------------------------------------
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         dx=(2.0d0*Pi)/plonmax(jeg)  ! compute longitude spacing
         iptr=1
         do ie=1,ielem(jeg)
            do i=1,nlon
               lon(i,ie,je)=(iptr-1)*dx ! compute nodal longitudes in rad
               iptr=iptr+1
            end do
         end do
      end do


      end








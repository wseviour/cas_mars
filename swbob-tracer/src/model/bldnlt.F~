      subroutine bldnlt(phy,
     $                  nlt,
     $                  pkvcoef,
     $                  pktcoef,
     $                  rcoef,
     $                  pfmechu,
     $                  pfmechv,
c     $                  pfum1,
     $                  pfum1n,
     $                  pfum1t,
     $                  lat,
     $                  lon,
     $                  umax)
      implicit none

#include <dims.h>
#include <stptrs.h>
#include <nltptrs.h>
#include <modeltime.h>
#include <constants.h>
#include <planet.h>


      REAL_TYPE phy(nlon,nlat2,nphy,0:1,plev)
      REAL_TYPE nlt(nlon,nlat2,nnlt,0:1,plev)
      REAL_TYPE pkvcoef(nlat2,0:1,plev)
      REAL_TYPE pktcoef(nlat2,0:1,plev)
      REAL_TYPE rcoef(nlon,nlat2,0:1,plev)
      REAL_TYPE pfmechu(nlon,nlat2,0:1,plev)
      REAL_TYPE pfmechv(nlon,nlat2,0:1,plev)
c      REAL_TYPE pfum1(nlon,nlat2,0:1,plev)
      REAL_TYPE pfum1n(nlat2,0:1,plev)   ! zonal sum on each node/tile
      REAL_TYPE pfum1t(nlat2,0:1,plev)   ! overall zonal mean
      REAL_TYPE lat(nlat2)
      REAL_TYPE lon(nlon)

      INT_TYPE i,j,k,ins
      REAL_TYPE cslat(nlat2),snlat(nlat2)
      REAL_TYPE umax, latsign, tmpq

#ifdef FTOPO
      REAL_TYPE rampdays,ramp,rmpc
      common /toporamp/ ramp,rmpc
#endif

c----------------------------------------
c     Executable Statements
c---------------------------------------- 

c     Compute local zonal CFL
#if 1
      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               do i=1,nlon
                  umax=DMAX1(umax,DABS( phy(i,j,u_p,ins,k)/DCOS(lat(j)) ))
               end do
            end do
         end do
      end do
#endif

c     t=0 assume phy(i,j,u_p,...) zonally symmetric, use value at i=1:
#ifdef FURLX
      if (itime.eq.itime0) then 
      do k=1,plev
         do ins=0,1
            do j=1,nlat2
               pfum1t(j,ins,k)=phy(1,j,u_p,ins,k)
            enddo
         enddo
      enddo
      endif
#endif

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

                  nlt(i,j,c_p,ins,k)=phy(i,j,u_p,ins,k)*phy(i,j,phi_p,ins,k)

                  nlt(i,j,d_p,ins,k)=phy(i,j,v_p,ins,k)*phy(i,j,phi_p,ins,k)

                  nlt(i,j,e_p,ins,k)=(phy(i,j,u_p,ins,k)*phy(i,j,u_p,ins,k)+
     $                                phy(i,j,v_p,ins,k)*phy(i,j,v_p,ins,k))
     $                                /(2.0d0*cslat(j)*cslat(j))

#ifdef FDRAG
                  nlt(i,j,f_p,ins,k)= pkvcoef(j,ins,k)*(phy(i,j,vorm1_p,ins,k)
     $                 -latsign*2.d0*Omega*snlat(j))
#endif
#ifdef FRLAX
                  nlt(i,j,g_p,ins,k)= pktcoef(j,ins,k)*phy(i,j,them1_p,ins,k)
#endif
#ifdef RTERM
                  tmpq=rcoef(i,j,ins,k)
#ifdef FRLAX
                  tmpq=tmpq-nlt(i,j,g_p,ins,k)
#endif
                 tmpq=max(tmpq,0.d0)*cslat(j)/(phy(i,j,phi_p,ins,k)+phibar)

c     R=-u*Q/h  ; need A->A-Ru*coslat, B->B+Rv*coslat
c                 (see Hack&Jacob, 1992, appendix A)

c                 nlt(i,j,ru_p,ins,k)=tmpq*phy(i,j,u_p,ins,k)/phy(i,j,phi_p,ins,k)
c                 nlt(i,j,rv_p,ins,k=)tmpq*phy(i,j,v_p,ins,k)/phy(i,j,phi_p,ins,k)
                  nlt(i,j,a_p,ins,k)=nlt(i,j,a_p,ins,k)+tmpq*phy(i,j,v_p,ins,k)
                  nlt(i,j,b_p,ins,k)=nlt(i,j,b_p,ins,k)-tmpq*phy(i,j,u_p,ins,k)
#endif
c     oct23, 2012, NYC: implement FMECH here rather than in advance.
#ifdef FMECH
c     constant force (can multiply above by rmpc, defined in advance):
                  nlt(i,j,a_p,ins,k)=nlt(i,j,a_p,ins,k)-pfmechv(i,j,ins,k)
                  nlt(i,j,b_p,ins,k)=nlt(i,j,b_p,ins,k)+pfmechu(i,j,ins,k)
#endif
#ifdef FURLX
c     relax zonal wind to prescribed profile fmechu defined in initmomfor:
                  nlt(i,j,b_p,ins,k)=nlt(i,j,b_p,ins,k)
     .                 +pkvcoef(j,ins,k)*(pfmechu(i,j,ins,k)-pfum1t(j,ins,k))
c     .                 +pkvcoef(j,ins,k)*(pfmechu(i,j,ins,k))
c                  pfum1(i,j,ins,k)=phy(i,j,u_p,ins,k)
                  pfum1n(j,ins,k)=pfum1n(j,ins,k)+phy(i,j,u_p,ins,k)
c     (fum1 is u*cosphi at previous timestep, stored in fordamp.h)
#endif

               end do
            end do
         end do
      end do

c#ifdef FRELAX
c      print*,'in bldnlt;',nlt(1,1,g_p,0,1)
c#endif
      end 

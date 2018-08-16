
      subroutine initscg(ounit)

      implicit none


#include <dims.h>
#include <stptrs.h>
#include <scptrs.h>
#include <coord.h>
#include <decomp.h>
#include <physical.h>
#include <mpcluster.h>
#include <planet.h>
#include <io.h>

      REAL_TYPE phi(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE vor(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE div(nlon,nlat2,0:1,plev,ielemd,jelemd)
      REAL_TYPE sig(nlon,nlat2,0:1,plev,ielemd,jelemd)

      REAL_TYPE uini(plat,plev),vorini(plat,plev),thini(plat,plev)
      REAL_TYPE sigini(plat,plev)

      REAL_TYPE vorsc(0:1,nsc*plev)
      REAL_TYPE divsc(0:1,nsc*plev)
      REAL_TYPE phisc(0:1,nsc*plev)
      REAL_TYPE sigsc(0:1,nsc*plev)

      REAL_TYPE scg(0:1,nscg)
      REAL_TYPE fld_file,level_file,sqrt2

      INT_TYPE  k,ie,je,jeg
      INT_TYPE  i,j,jptr,ins
      INT_TYPE  sc_p,mch,m,n
      INT_TYPE  itime
      INT_TYPE  ounit               ! file output unit number

c-------------------------------------------------------------
c     Executable Statements
c-------------------------------------------------------------

      sqrt2=DSQRT(2.0D0)
c
c     define uini and calculate vorini by diff wrt phi,
c     calculate dthetadphi from uini_z and integrate wrt phi -> thini,
c     convert vorini to ABSOLUTE vorticity:

      call inituth(uini,vorini,thini,sigini)
c
c     distribute vorini and thini among nodes and tiles
      do k=1,plev
         jptr=jbeg(rank)
         do jeg=jebeg(rank),jeend(rank) ! loop over elements
            je=jeg-jebeg(rank)+1
            do j=1,nlat2
               do ie=1,ielem(jeg)
                  do i=1,nlon
                     vor(i,j,0,k,ie,je)=vorini(jptr,k)
                     vor(i,j,1,k,ie,je)=vorini(plat+1-jptr,k)
                     div(i,j,0,k,ie,je)=0.d0
                     div(i,j,1,k,ie,je)=0.d0
                     phi(i,j,0,k,ie,je)=thini(jptr,k)
                     phi(i,j,1,k,ie,je)=thini(plat+1-jptr,k)
                     sig(i,j,0,k,ie,je)=sigini(jptr,k)
                     sig(i,j,1,k,ie,je)=sigini(plat+1-jptr,k)
                  enddo
               enddo
               jptr=jptr+1
            end do
         end do
      end do


c     Add perturbation to theta field (bump localized in longitude)
      do jeg=jebeg(rank),jeend(rank)
         je=jeg-jebeg(rank)+1
         do ie=1,ielem(jeg)
            do k=1,plev
               call initpert(vor(1,1,0,k,ie,je),
     $                      lat(1,jeg),
     $                      cslat(1,jeg),
     $                      snlat(1,jeg),
     $                      lon(1,ie,je),
     $                      k)
            end do
         end do
      end do


      call physical_to_spectral(vor,vorsc,plev)
      call physical_to_spectral(div,divsc,plev)
      call physical_to_spectral(phi,phisc,plev)
      call physical_to_spectral(sig,sigsc,plev)

c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do itime=0,1                  ! store two identical timelevels

         fld_file=vor_p             ! vorticity field
         do k=1,plev
            call sc_pgather(vorsc,scg,plev,k)
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file
               call writereal(scg,2*nscg,ounit)
            end if
         enddo

         fld_file=div_p             ! div field
         do k=1,plev
            call sc_pgather(divsc,scg,plev,k)
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file
               call writereal(scg,2*nscg,ounit)
            end if
         enddo

         fld_file=theta_p           ! theta/phi field
         do k=1,plev
            call sc_pgather(phisc,scg,plev,k)
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file
               call writereal(scg,2*nscg,ounit)
            end if
         enddo

         fld_file=sig_p             ! tracer field
         do k=1,plev
            call sc_pgather(sigsc,scg,plev,k)
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file
               call writereal(scg,2*nscg,ounit)
            end if
         enddo

      end do

      end

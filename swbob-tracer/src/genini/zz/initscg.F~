
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

      REAL_TYPE vorini(plat,plev),thini(plat,plev)

      REAL_TYPE vorsc(0:1,nsc*plev)  
      REAL_TYPE divsc(0:1,nsc*plev)  
      REAL_TYPE phisc(0:1,nsc*plev)
      REAL_TYPE forsc(0:1,nsc*plev)  

      REAL_TYPE scg(0:1,nscg)
      REAL_TYPE fld_file,level_file

      INT_TYPE  k,ie,je,jeg
      INT_TYPE  i,j,jptr,ins
      INT_TYPE  sc_p,mch,m,n
      INT_TYPE  itime
      INT_TYPE  ounit               ! file output unit number

c-------------------------------------------------------------
c     Executable Statements
c-------------------------------------------------------------

      call initvorth(vorini,thini)
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
                  enddo
               enddo
               jptr=jptr+1
            end do
         end do
      end do
 
      call physical_to_spectral(vor,vorsc,plev)
      call physical_to_spectral(div,divsc,plev)
      call physical_to_spectral(phi,phisc,plev)

#ifdef FSPEC
      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            do k=1,plev
               sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
               do n=m,nn
                  forsc(0,sc_p+n-m)=0.d0
                  forsc(1,sc_p+n-m)=0.d0
               enddo
            enddo
         enddo
      enddo
#endif
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! transform phy to sc

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

      enddo

#ifdef FSPEC
         fld_file=nsp+1             ! forcing field
         do k=1,plev
            call sc_pgather(forsc,scg,plev,k)
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file   
               call writereal(scg,2*nscg,ounit)
            end if
         enddo
#endif

      
      end







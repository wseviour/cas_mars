      subroutine out_resc(lrevorsc,lredivsc,lrephisc,reounit)

c     Calls: spectral/sc_pgather, utils/fldio/writereal
c------------------------------------------------
      implicit none
#include <type.h>
#include <dims.h>
c      parameter(plev=PLEV)
c      parameter(mm=(plon-1)/3,nn=mm)
c      INT_TYPE nscg  ! global spectral coefficient counts
c      parameter(nscg=(nn+1)*(mm+1)-((mm+1)*mm)/2)
c      parameter(nvsc=mm2d*(nn+4),nsc=mm2d*(nn+2))
#include <scptrs.h>
c      INT_TYPE scgptr(0:mm),scptr(0:mm),vscptr(0:mm)
c      common /scptrs/scgptr,scptr,vscptr
#include <stptrs.h>
c     vor_p ... etc
#include <mpcluster.h>
c      INT_TYPE nThreads,rank,root
c-----------------------------------------------
c     Local declaration
c-----------------------------------------------

      INT_TYPE k,m,n,ic
      INT_TYPE  sc_p            
      REAL_TYPE scgrep(0:1,nscg)
      REAL_TYPE lrevorsc(0:1,nsc*plev)
      REAL_TYPE lredivsc(0:1,nsc*plev)
      REAL_TYPE lrephisc(0:1,nsc*plev)
      REAL_TYPE fld_file,level_file
      INT_TYPE reounit

c-----------------------------------------------
c     executable statements
c-----------------------------------------------

         ! vorticity field
         fld_file=vor_p
c         write(*,*)'vorticity field'
         do k=1,plev
            call sc_pgather(lrevorsc,scgrep,plev,k)

            if(k.eq.10) then
               print*,'vorticity field'
               do n=1,nscg
                  print*,n,scgrep(0,n),scgrep(1,n)
               enddo
            endif

            level_file=k
            if (rank.eq.root) then
c               write(*,*)'before the write'
               write(reounit)level_file,fld_file
               ! 2 comes from Re/Im
               ! 1 field and 1 level at a time
               call writereal(scgrep,2*nscg,reounit)
            end if
         end do

         ! div field
         fld_file=div_p

         do k=1,plev
            call sc_pgather(lredivsc,scgrep,plev,k)

            if(k.eq.10) then
               print*,'divergence field'
               do n=1,nscg
                  print*,n,scgrep(0,n),scgrep(1,n)
               enddo
            endif

            level_file=k
            if (rank.eq.root) then
               write(reounit)level_file,fld_file     
               call writereal(scgrep,2*nscg,reounit)
            end if
         end do

         ! theta/phi field
         fld_file=theta_p

         do k=1,plev
            call sc_pgather(lrephisc,scgrep,plev,k)

            if(k.eq.10) then
               print*,'theta field'
               do n=1,nscg
                  print*,n,scgrep(0,n),scgrep(1,n)
               enddo
            endif

            level_file=k
            if (rank.eq.root) then
               write(reounit)level_file,fld_file     
               call writereal(scgrep,2*nscg,reounit)
            end if
         end do

      end

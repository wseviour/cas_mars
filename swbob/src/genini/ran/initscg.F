
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
      do n=0,nn
         an=float(n)
         edensity(n)=amp*an**(gamma/2.d0)/(an+n0)**gamma
         vornsq(n)=an*(an+1)/(radius**2*(2.d0*an+1))*edensity(n)
      enddo
c
      if (rank.eq.root) then
         iseed=1
         call random_seed(iseed)
         do n=0,nn
            an=float(n)
            edensity(n)=amp*an**(gamma/2.d0)/(an+n0)**gamma
            vornsq(n)=an*(an+1)/(radius**2*(2.d0*an+1))*edensity(n)
            call random_number(uni)
            call random_number(pm)
            tralph(i)=2.d0*uni-1.d0
            trbeta(i)=sign(1.d0,2.d0*pm-1.d0)*sqrt(1.d0-tralph(i)**2)
            write(*,1)"In getinputs, i,tralph,trbeta = ", i,tralph(i),trbeta(i)
         enddo
      endif

 1    format(a40,1x,i4,1x,2(1x,f8.4))

      call broadcastreal(tralph,nn)
      call broadcastreal(trbeta,nn)

      do mch=0,1
         do m=mbeg(mch,rank),mend(mch,rank)
            do k=1,plev
               sc_p=1+plev*scptr(m)+(nn+1-m)*(k-1)
               do n=m,nn
                  vorsc(0,sc_p+n-m)=rnda(n)
                  vorsc(1,sc_p+n-m)=rndb(n)
               enddo
            enddo
         enddo
      enddo
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do itime=0,1                  ! store two identical timelevels

         fld_file=vor_p             ! vorticity field
         do k=1,plev
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file   
               call writereal(vorscg,2*nscg,ounit)
            end if
         enddo

         fld_file=div_p             ! div field
         do k=1,plev
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file     
               call writereal(divscg,2*nscg,ounit)
            end if
         enddo

         fld_file=theta_p           ! theta/phi field
         do k=1,plev
            level_file=k
            if (rank.eq.root) then
               write(ounit)level_file,fld_file     
               call writereal(phiscg,2*nscg,ounit)
            end if
         enddo

      end do
      
      end







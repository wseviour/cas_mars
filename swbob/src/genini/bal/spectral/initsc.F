       subroutine initsc(sc,ilev,ifld)
       implicit none

c     Returns sc(local coeffs) from scg(global coeffs)

#include <dims.h>
#include <stptrs.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
#include <io.h>

       REAL_TYPE   sc(2*(plev*nsp)*nsc)  ! time level n state variable spectral coefficients

       REAL_TYPE   scg(0:1,nscg)    !global spectral coefficients
       REAL_TYPE   scbuf(0:1,0:nsc) ! comm buffer for spectral coefs + two slots for field info
       REAL_TYPE   rfld,rlev
       INT_TYPE    ifld,ilev
       INT_TYPE    fld_p,level
       INT_TYPE    sc_p
       INT_TYPE    ierr

       INT_TYPE    boundint
       external    boundint

       INT_TYPE    m,mch,ipe

#if 1
      REAL_TYPE sqrt2
      INT_TYPE  n
#endif
     
       if (rank.eq.root) then

c     define spectral coefficients for vorticity and distribute:

          call initscg(scg,ilev,ifld)

          rlev=float(ilev)
          rfld=float(ifld)

          do ipe=0,M_NODE-1
             
             do mch=0,1
                do m=mbeg(mch,ipe),mend(mch,ipe)
                   call copyreal(scg(0,scgptr(m)+1),scbuf(0,scptr(m)+1),2*(nn+1-m))
                end do
             end do

             scbuf(0,0)=rlev
             scbuf(1,0)=rfld

c     send scbuf to process ipe:

             if (ipe.ne.root) then
                call sendreal(scbuf,2*(nsc+1),ipe,ipe)
             else
                level=NINT(scbuf(0,0))
                fld_p=NINT(scbuf(1,0))
                do mch=0,1
                   do m=mbeg(mch,rank),mend(mch,rank)
                      sc_p=plev*nsp*2*scptr(m)+1
                      call unpacksc(scbuf(0,scptr(m)+1),sc(sc_p),m,fld_p,level)
                   end do
                end do
             endif

          end do                ! loop on ipe

       else

          ! receive spectral coefficient slices decode them and stow them
          call recvreal(scbuf,2*(nsc+1),root,rank)
          level=NINT(scbuf(0,0))
          fld_p=NINT(scbuf(1,0))
          do mch=0,1
             do m=mbeg(mch,rank),mend(mch,rank)
                sc_p=plev*nsp*2*scptr(m)+1
                call unpacksc(scbuf(0,scptr(m)+1),sc(sc_p),m,fld_p,level)
             end do
          end do

       end if


       end 





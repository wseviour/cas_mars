       subroutine initsc(sc,iunit)
       implicit none

c     Returns sc(local coeffs) from scg(global coeffs)

#include <dims.h>
#include <stptrs.h>
#include <scptrs.h>
#include <decomp.h>
#include <mpcluster.h>
#include <io.h>

       REAL_TYPE   sc(2*(plev*nsp)*nsc)  ! time level n state variable spectral coefficients
       INT_TYPE    iunit                 ! unit number for file input

       REAL_TYPE   scg(0:1,nscg)    !global spectral coefficients
       REAL_TYPE   scbuf(0:1,0:nsc) ! comm buffer for spectral coefs + two slots for field info
       REAL_TYPE   level_file
       REAL_TYPE   fld_file
       INT_TYPE    fld_p
       INT_TYPE    level
       INT_TYPE    sc_p
       INT_TYPE    ierr

       INT_TYPE    boundint
       external    boundint

       INT_TYPE    m,mch,ipe

#if 0
      REAL_TYPE sqrt2
      INT_TYPE  n
#endif
     
       if (rank.eq.root) then

          ! read spectral coefficients for vor,div,theta
          ! and distribute

          read(iunit)level_file,fld_file

          level=NINT(level_file)
          fld_p=NINT(fld_file)

          ierr=0
          ierr=ierr+boundint(level,1,plev,'level  ')
          ierr=ierr+boundint(fld_p,1,nsp ,'fld ptr')

          if (ierr.eq.0) then

             call readreal(scg,2*nscg,iunit)
#if 0
             write(MPI_UNIT,*)"------------------------------"
             write(MPI_UNIT,*)"INITIAL GLOBAL SPECTRAL COEFFICIENTS:",fld_p,level
             write(MPI_UNIT,*)"------------------------------"
             write(MPI_UNIT,*)"  m   n        RE           IM"
             write(MPI_UNIT,*)"------------------------------"
#endif

             do ipe=0,M_NODE-1

                do mch=0,1
                   do m=mbeg(mch,ipe),mend(mch,ipe)
                      call copyreal(scg(0,scgptr(m)+1),scbuf(0,scptr(m)+1),2*(nn+1-m))
                   end do
                end do

                scbuf(0,0)=level_file
                scbuf(1,0)=fld_file
#if 0
                sqrt2=DSQRT(2.0)
                do mch=0,1
                   do m=mbeg(mch,ipe),mend(mch,ipe)
                      write(MPI_UNIT,*)
                      do n=m,nn
                         write(MPI_UNIT,10)m,n,scbuf(0,scptr(m)+1+n-m)/sqrt2,scbuf(1,scptr(m)+1+n-m)/sqrt2
 10                      format(i4,2x,i4,2x,e22.15,2x,e22.15)
                      end do
                   end do
                end do
#endif

                ! send scbuf to process ipe

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

             end do
#if 0
             write(MPI_UNIT,*)"------------------------------"
#endif
          end if
              
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





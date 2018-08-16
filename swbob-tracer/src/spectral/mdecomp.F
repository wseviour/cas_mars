      subroutine mdecomp()
      implicit none

c CALLED BY initsht itself by initmodel itself by main
c CALLS: decomp.F (in spectral)
c GOAL: sets mbeg(0:1,ipe),mend(0:1,ipe),vscsum(ipe),scsum(ipe)
       
c in params.h #define M_NODE 8 for ex actual number of Pes

#include <type.h>  
#include <dims.h>
#include <decomp.h>
c      INT_TYPE mbeg(0:1,0:M_NODE-1),mend(0:1,0:M_NODE-1),mchmax(0:M_NODE-1)
c      0:1 for the 2 chunks small m and big m
c      INT_TYPE jbeg(0:P_NODE-1),jend(0:P_NODE-1)     
c      INT_TYPE jebeg(0:P_NODE-1),jeend(0:P_NODE-1)     
#include <mpcluster.h>
#include <io.h>
c     MPI_UNIT

      INT_TYPE ipe,mch,m
      INT_TYPE mbase(0:M_NODE) ! built in decomp.F
      INT_TYPE vscsum(0:M_NODE-1) ! sum of viscous spec coeffs
      INT_TYPE scsum(0:M_NODE-1) ! sum of spec coefs 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call decomp(mbase,mm/2+1,M_NODE)

#if 0
      write(MPI_UNIT,*)" ---------------------------- "
      write(MPI_UNIT,*)" mbase(0:M_NODE) in mdecomp.F "
      do ipe=0,M_NODE
         write(MPI_UNIT,*)ipe, mbase(ipe)
      enddo
      write(MPI_UNIT,*)" ---------------------------- "
#endif

      do ipe=0,M_NODE-1
         mbeg(0,ipe)=mbase(ipe )
         mend(0,ipe)=mbase(ipe+1)-1
         mbeg(1,ipe)=mm+2-mend(0,ipe)
         mend(1,ipe)=mm+2-mbeg(0,ipe)
         do mch=0,1
            mbeg(mch,ipe)=mbeg(mch,ipe)-1
            mend(mch,ipe)=mend(mch,ipe)-1
         end do
         mchmax(ipe)=1
      end do

      ! correct pairing algorithm for folding for even truncation 
      ! (which results in an odd number of wave numbers)

      if (MOD(mm,2).eq.0) then
         if (mbeg(0,M_NODE-1).ne.mend(0,M_NODE-1)) then
            mbeg(1,M_NODE-1)=mbeg(1,M_NODE-1)+1
         else
            mchmax(M_NODE-1)=0
         end if
      end if

#if 1
      write(MPI_UNIT,*)
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" Domain Decomposition by Wave Numbers "
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" pe     range      range     number of vector"
      write(MPI_UNIT,*)"        (low)      (high)     spectral coefs "
      write(MPI_UNIT,*)" mbeg(0:1,ipe),mend(0:1,ipe),vscsum(ipe)     "
      write(MPI_UNIT,*)"============================================="
#endif

      do ipe=0,M_NODE-1
         vscsum(ipe)=0
         scsum(ipe)=0
         do mch=0,1
            do m=mbeg(mch,ipe),mend(mch,ipe)
               vscsum(ipe)=vscsum(ipe)+(nn+2-m)
                scsum(ipe)=scsum(ipe)+(nn+1-m)
            end do
         end do
#if 1
         if (mchmax(ipe).eq.1) then
            write(MPI_UNIT,20)ipe,mbeg(0,ipe),mend(0,ipe),mbeg(1,ipe),mend(1,ipe),vscsum(ipe)
 20         format(i4,"  ",i4,"-",i4,"  ",i4,"-",i4,"    ",i7)
         else
            write(MPI_UNIT,30)ipe,mbeg(0,ipe),mend(0,ipe),vscsum(ipe)
 30         format(i4,"  ",i4,"-",i4,"               ",i7)
         end if
#endif
      end do

      end

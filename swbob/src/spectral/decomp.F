      subroutine decomp(beg,len,npe)
      implicit none

c CALLED BY: spectral/mdecomp.F : call decomp(mbase,mm/2+1,M_NODE)

c GOAL: returns beg(0:npe), decomposes len in npe segments as
c       evenly as possible ie if npe doesn't divide len exactly,
c       the mod(len/npe) first segment receive one more indices
c       so that the segemnt cover all the indices.

#include <type.h> 
c  #define INT_TYPE integer
#include <io.h>
c     MPI_UNIT

      INT_TYPE len,npe,n
      INT_TYPE beg(0:npe)

ccccccccccccccccccccccccccccccccccccccccccc

#if 0
      write(MPI_UNIT,*)" ---------------------------- "
      write(MPI_UNIT,*)" In decomp.F: len,npe "
      write(MPI_UNIT,*)len, npe
      write(MPI_UNIT,*)" ---------------------------- "
#endif

      beg(0)=1
      do n=1,npe-1

         ! if npe doesn't divide len exactly, add
         ! and indice for the mod(len,npe) first segemnts

         if (n.le.mod(len,npe)) then 
            beg(n)=beg(n-1)+(len-1)/npe+1
         else
            beg(n)=beg(n-1)+len/npe
         end if

      end do
      beg(npe)=len+1

      end

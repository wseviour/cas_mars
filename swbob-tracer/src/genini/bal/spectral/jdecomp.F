      subroutine jdecomp()
      implicit none

#include <decomp.h>
#include <io.h>

      INT_TYPE ipe
      INT_TYPE jbase(0:P_NODE)

      call decomp(jbase,plat/2,P_NODE)

      do ipe=0,P_NODE-1
         jbeg(ipe)=jbase(ipe)
         jend(ipe)=jbase(ipe+1)-1
      end do

      write(MPI_UNIT,*)
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" Domain Decomposition by Latitude Pair Chunks"
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" pe     latitude                             " 
      write(MPI_UNIT,*)"         range                               "
      write(MPI_UNIT,*)"============================================="

      do ipe=0,P_NODE-1
         write(MPI_UNIT,20)ipe,jbeg(ipe),jend(ipe)
 20      format(i4,"   ",i4,"-",i4,"           ",i7)
      end do

      write(MPI_UNIT,*)"============================================="

      ! compute decomposition of plat/2 into elements of length nlat2
      ! for spectral transform

      call decomp(jbase,jelemg,P_NODE)

      do ipe=0,P_NODE-1
         jebeg(ipe)=jbase(ipe)
         jeend(ipe)=jbase(ipe+1)-1
      end do

      end





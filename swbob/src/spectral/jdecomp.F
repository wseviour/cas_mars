      subroutine jdecomp()
      implicit none

c GOAL: get jbeg(ipe), jend(ipe)
c           jebeg(ipe), jeend(ipe)

c     P_NODE in params.h
c     include/dims.h:  parameter(nlat2=4,jelemg=(plat/2)/nlat2)  ! Note: plat/2 assumed divisible by nlat2

#include <decomp.h>
c      INT_TYPE jbeg(0:P_NODE-1),jend(0:P_NODE-1)
c      INT_TYPE jebeg(0:P_NODE-1),jeend(0:P_NODE-1)
#include <io.h>

      INT_TYPE ipe
      INT_TYPE jbase(0:P_NODE) ! decomposition of plat/2 over P_NODE

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call decomp(jbase,plat/2,P_NODE)

      do ipe=0,P_NODE-1
         jbeg(ipe)=jbase(ipe)
         jend(ipe)=jbase(ipe+1)-1
      end do

#if 1
      write(MPI_UNIT,*)
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" Domain Decomposition by Latitude Pair Chunks"
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" pe     latitude                             " 
      write(MPI_UNIT,*)"         range                               "
      write(MPI_UNIT,*)" jbeg(ipe),jend(ipe)                         "
      write(MPI_UNIT,*)"============================================="

      do ipe=0,P_NODE-1
         write(MPI_UNIT,20)ipe,jbeg(ipe),jend(ipe)
 20      format(i4,"   ",i4,"-",i4,"           ",i7)
      end do

      write(MPI_UNIT,*)"============================================="
#endif

      ! compute decomposition of plat/2 into elements of length nlat2
      ! for spectral transform

      call decomp(jbase,jelemg,P_NODE)

      do ipe=0,P_NODE-1
         jebeg(ipe)=jbase(ipe)
         jeend(ipe)=jbase(ipe+1)-1
      end do

#if 1
      write(MPI_UNIT,*)
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" Domain Decomposition by elements            "
      write(MPI_UNIT,*)"============================================="
      write(MPI_UNIT,*)" pe    jebeg(ipe), jeend(ipe)                 "
      write(MPI_UNIT,*)"============================================="

      do ipe=0,P_NODE-1
         write(MPI_UNIT,21)ipe,jebeg(ipe),jeend(ipe)
 21      format(i4,"   ",i4,"-",i4)
      end do

      write(MPI_UNIT,*)"============================================="
#endif

      end





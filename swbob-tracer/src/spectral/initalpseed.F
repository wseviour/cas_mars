      subroutine initalpseed(colatitude)
      implicit none

c-------------------------------------------------------------------
c     Compute Pseed: Seed values to initiate the recurrence generating
c     the Legendre Polynomials
c     CF P Swarztrauber, The Vector Harmonic Transform Method for Solving 
c     Partial Differential Equations in Spherical Geometry
c     MWR, vol 121, no 12, pp 3416-3437 especially here pp 3428-3429
c-------------------------------------------------------------------

c--------------------------------------------------------------------
c      Include Files
c--------------------------------------------------------------------
#include <alp.h>     
c      INT_TYPE  mm,nn
c      parameter(mm=(plon-1)/3,nn=mm)
c      REAL_TYPE Pseed(nlat2,0:nn+1,jelemg,0:1,0:1)
c      common /alp/ Pseed, genp, P
#include <decomp.h>

#include <mpcluster.h>


      REAL_TYPE colatitude(plat/2)    ! colatitudes 

      INT_TYPE mch,ioe,m,n
      INT_TYPE jptr,je,j

      REAL_TYPE a((nn+1)/2+1)
      REAL_TYPE Ptmp

c-------------------------------------------------------------
c     Executable Statements
c-------------------------------------------------------------
      do mch=0,1
         do ioe=0,1
            m=mbeg(mch,rank)+ioe
            do n=m,nn+1

               call dnlfk(m,n,a) ! computes "a" coefs of Swarztrauber 1993

               jptr=1
               do je=1,jelemg
                  do j=1,nlat2
                     call dnlft(m,
     $                          n,
     $                          colatitude(jptr),
     $                          a,
     $                          Ptmp)
                     Pseed(j,n,je,ioe,mch)=Ptmp ! properly normalize Pmn's interval -1 to 1
                     jptr=jptr+1
                  end do
               end do

            end do
         end do
      end do
      
      end




























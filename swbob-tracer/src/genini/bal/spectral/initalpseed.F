      subroutine initalpseed(colatitude)
      implicit none

!--------------------------------------------------------------------
!      Include Files
!--------------------------------------------------------------------
#include <alp.h>
!--------------------------------------------------------------------
#include <decomp.h>
!--------------------------------------------------------------------
#include <mpcluster.h>
!--------------------------------------------------------------------

      REAL_TYPE colatitude(plat/2)    ! colatitudes 

      INT_TYPE mch,ioe,m,n
      INT_TYPE jptr,je,j

      REAL_TYPE a((nn+1)/2+1)
      REAL_TYPE Ptmp

      do mch=0,1
         do ioe=0,1
            m=mbeg(mch,rank)+ioe
            do n=m,nn+1

               call dnlfk(m,n,a)

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




























      subroutine initgenalp(genp,
     $                      m)
      implicit none

#include <dims.h>
      
      INT_TYPE m
      REAL_TYPE genp(3,m:nn+1)

      INT_TYPE n
      REAL_TYPE fac1,fac2,fac3

      if (m.lt.2) then
         do n=m,nn+1
            genp(1,n)=0.0d0
            genp(2,n)=0.0d0
            genp(3,n)=0.0d0
         end do
      else
         do n=m,nn+1

            ! ==================================================
            !
            ! See equations 4.20 and 5.14 in Paul Swarztrauber's
            ! MWR December 1993 article, V 121, no. 12, 
            ! pages 3415-3437.
            !
            ! ==================================================

            fac1=DSQRT( (2.0d0*n+1.0d0)*(n+m-2.0d0)*(n+m-3.0d0) / 
     $                 ((2.0d0*n-3.0d0)*(n+m)*(n+m-1.0d0)) )
         
            fac2=DSQRT( (n-m+2.0d0)*(n-m+1.0d0) / 
     $                 ((n+m)*(n+m-1.0d0)) )

            fac3=DSQRT( (2.0d0*n+1.0d0)*(n-m)*(n-m-1.0d0) /
     $                 ((2.0d0*n-3.0d0)*(n+m)*(n+m-1.0d0)) )

            genp(1,n) =  fac1
            genp(2,n) = -fac2
            genp(3,n) =  fac3

         end do

      end if

      end







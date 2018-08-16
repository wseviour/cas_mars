      subroutine myrand(uni)
      implicit none

#include <type.h>

      REAL_TYPE uni
      INT_TYPE iseed

#ifdef AIX
      INT_TYPE seed(100),m,im
#endif

#ifdef AIX
      if (uni.eq.-1.d0) then
         call random_seed(size=m)
         do im=1,m
            seed(im)=myid
         enddo
         call random_seed(generator=2)
         call random_seed(put=seed(1:m))
      else
         call random_number(harvest=uni)
         uni=2.d0*uni-1.d0
      endif
#endif

#if PG
      if (uni.eq.-1.d0) then
         call random_seed(iseed)
      else
         call random_number(uni)
         uni=2.d0*uni-1.d0
      endif
#endif

#if MPG
      if (uni.eq.-1.d0) then
         call random_seed(iseed)
      else
         call random_number(uni)
         uni=2.d0*uni-1.d0
      endif
#endif

#if GNU
      if (uni.eq.-1.d0) then
         call srand(0)
      else
         uni=2.d0*rand(0)-1.d0
      endif
#endif

      return
      end




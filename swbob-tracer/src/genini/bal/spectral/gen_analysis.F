      subroutine gen_analysis(fc,sc,P,gw,lat,m,nmax,nlev)
      implicit none

#include <dims.h>
#include <io.h>

      INT_TYPE  m
      INT_TYPE  nmax
      INT_TYPE  nlev
      REAL_TYPE fc(0:1,nlat2,0:1,nlev)
      REAL_TYPE sc(0:1,m:nmax,nlev)
      REAL_TYPE P(nlat2,0:nn+1)
      REAL_TYPE gw(nlat2)
      REAL_TYPE lat(nlat2)

      INT_TYPE j,k,n,ioe

      do k=1,nlev

         call gen_symf(fc(0,1,0,k),gw)

         do ioe=0,1
            do n=m+ioe,nmax,2
               do j=1,nlat2
                  sc(0,n,k) = sc(0,n,k)+fc(0,j,ioe,k)*P(j,n)
                  sc(1,n,k) = sc(1,n,k)+fc(1,j,ioe,k)*P(j,n)
               end do
            end do
         end do

      end do

      end 

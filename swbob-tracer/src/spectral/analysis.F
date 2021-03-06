#define K_UNROLL 6
      subroutine Analysis(fc,nltsc,P,gw,m)
      implicit none

#include <dims.h>
#include <nltptrs.h>

      INT_TYPE m
      REAL_TYPE fc(0:1,nlat2,nnlt,0:1,plev)
      REAL_TYPE nltsc(0:1,nnlt,m:nn+1,plev)
      REAL_TYPE P(nlat2,0:nn+1)
      REAL_TYPE gw(nlat2,nnlt)

      INT_TYPE j,k,n,ioe

      do k=1,plev

         call symf(fc(0,1,1,0,k),gw)

         do ioe=0,1
            do n=m+ioe,nn+1,2
               do j=1,nlat2
                  nltsc(0,a_p,n,k) = nltsc(0,a_p,n,k)+fc(0,j,a_p,ioe,k)*P(j,n)
                  nltsc(1,a_p,n,k) = nltsc(1,a_p,n,k)+fc(1,j,a_p,ioe,k)*P(j,n)

                  nltsc(0,b_p,n,k) = nltsc(0,b_p,n,k)+fc(0,j,b_p,ioe,k)*P(j,n)
                  nltsc(1,b_p,n,k) = nltsc(1,b_p,n,k)+fc(1,j,b_p,ioe,k)*P(j,n)

                  nltsc(0,c_p,n,k) = nltsc(0,c_p,n,k)+fc(0,j,c_p,ioe,k)*P(j,n)
                  nltsc(1,c_p,n,k) = nltsc(1,c_p,n,k)+fc(1,j,c_p,ioe,k)*P(j,n)

                  nltsc(0,d_p,n,k) = nltsc(0,d_p,n,k)+fc(0,j,d_p,ioe,k)*P(j,n)
                  nltsc(1,d_p,n,k) = nltsc(1,d_p,n,k)+fc(1,j,d_p,ioe,k)*P(j,n)

                  nltsc(0,e_p,n,k) = nltsc(0,e_p,n,k)+fc(0,j,e_p,ioe,k)*P(j,n)
                  nltsc(1,e_p,n,k) = nltsc(1,e_p,n,k)+fc(1,j,e_p,ioe,k)*P(j,n)

                  nltsc(0,f_p,n,k) = nltsc(0,f_p,n,k)+fc(0,j,f_p,ioe,k)*P(j,n)
                  nltsc(1,f_p,n,k) = nltsc(1,f_p,n,k)+fc(1,j,f_p,ioe,k)*P(j,n)

                  nltsc(0,g_p,n,k) = nltsc(0,g_p,n,k)+fc(0,j,g_p,ioe,k)*P(j,n)
                  nltsc(1,g_p,n,k) = nltsc(1,g_p,n,k)+fc(1,j,g_p,ioe,k)*P(j,n)

c     rks - added tracer field, nonlinear terms L, M, N:
                  
                  nltsc(0,l_p,n,k) = nltsc(0,l_p,n,k)+fc(0,j,l_p,ioe,k)*P(j,n)
                  nltsc(1,l_p,n,k) = nltsc(1,l_p,n,k)+fc(1,j,l_p,ioe,k)*P(j,n)

                  nltsc(0,m_p,n,k) = nltsc(0,m_p,n,k)+fc(0,j,m_p,ioe,k)*P(j,n)
                  nltsc(1,m_p,n,k) = nltsc(1,m_p,n,k)+fc(1,j,m_p,ioe,k)*P(j,n)

                  nltsc(0,n_p,n,k) = nltsc(0,n_p,n,k)+fc(0,j,n_p,ioe,k)*P(j,n)
                  nltsc(1,n_p,n,k) = nltsc(1,n_p,n,k)+fc(1,j,n_p,ioe,k)*P(j,n)


               end do
            end do
         end do
      end do

c      if (m.eq.0) print*,'in analysis;',fc(0,1,g_p,1,1)


      end 

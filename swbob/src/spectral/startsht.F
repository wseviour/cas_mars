       subroutine startsht(iunit)
       implicit none

#include <spectral.h>
#include <fields.h>
#include <decomp.h>
#include <mpcluster.h>
#include <scptrs.h>

       INT_TYPE  iunit
       INT_TYPE  layer
       INT_TYPE  m, sc_p
       
       print *,"header file checks out, reading in spectral coefs..."
       print *
       
c     distribute spectral coefficients into sc array from initial file
c     sc and scm1 stored in spectral.h
       do layer=1,plev*nsp
          call initsc(sc,iunit)
c          print *
       end do
       
c     next get spectral coefficients from time level t-1 from file.
       do layer=1,plev*nsp
          call initsc(scm1,iunit)
c          print *
       end do

c     TODO: write a suitable initsc1
c     get forcing coefficients for Markovian forcing:
c       do layer=1,plev
c          call initsc1(fvorsc,iunit)
cc          print *
c       end do

c       print *,"spectral coefs distributed..."

c     transform sc -> sfc to finish up

c       print *,"inverse fourier transforming..."

c     blduv and go to fourier soace
       call spectral_to_fourier(nltsc,sc,scm1,sfcx,sfc)
 
c       print *,"ready to begin timestepping..."

       end 





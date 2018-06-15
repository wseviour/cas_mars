
      program GenGrid

c--------------------------------------------------------
c This is the main program that generates the grid:
c     It computes the gaussian lats (theta)
c     and gaussian weights (gw) and outputs them.
c--------------------------------------------------------

      implicit none

#include <dims.h>

      INT_TYPE lgwork
      parameter(lgwork= plat*(plat/plat+6) ) 

      INT_TYPE  j
      INT_TYPE  ierror,jeg

      INT_TYPE  lonmax(jelemg)
      REAL_TYPE theta(plat)
      REAL_TYPE gw(plat)
      REAL_TYPE gwork(lgwork)

      REAL_TYPE gwsum ! Sum of the gaussian weight, should be one

      INT_TYPE lenchr
      external lenchr

c Variables used to construct the name of the output containing the GRID
      character*1 trunc 
      character*80 stub 
      character*5 tail,tail_lj
      character*80 grdfile

      REAL_TYPE Pi
      parameter   (Pi=3.141592653589793D0)

c--------------------------------------------------------
c Executable statements
c--------------------------------------------------------

c--------------------------------------------------------
c     Build grid file name
c--------------------------------------------------------
      stub="GRID"
      trunc='T'
      write(tail,1)mm
 1    format(i4)
      call ljust(tail,tail_lj)
      grdfile=stub(1:lenchr(stub))//'.'//trunc//
     $     tail_lj(1:lenchr(tail_lj))

c--------------------------------------------------------
c     Open grid file
c--------------------------------------------------------
      open(unit=7,file=grdfile,form="FORMATTED",status="UNKNOWN")

c--------------------------------------------------------
c     Compute gaussian weights 
c--------------------------------------------------------
      call Pgaqd(plat,theta,gw,gwork,lgwork,ierror)

c--------------------------------------------------------
c     Output colatitudes/gaussian weights
c--------------------------------------------------------
      do j=1,plat/2
         write(7,10)j,theta(j),gw(j)
      end do
 10   format(i4,1x,e27.17,1x,e27.17)

c--------------------------------------------------------
c     Check gaussian weights
c--------------------------------------------------------
      gwsum=0.0d0
      do j=1,plat/2
         gwsum=gwsum+gw(j)
      end do
      print *,"sum of gaussian wts:",gwsum

c--------------------------------------------------------
c     Output longitude maxima of grid tiles
c--------------------------------------------------------
      do jeg=1,jelemg
         lonmax(jeg)=plon
         write(7,20)jeg,lonmax(jeg)
 20      format(i4,1x,i4)
      end do

      stop
      end

#include <type.h>
      subroutine gridfile()
      implicit none

c GOAL: On root node: open grdfile,
c       read jj,colatitude(j),gausswt(j)  j=1,plat/2
c       read jj,plonmax(jeg) jeg=1,jelemg     
c     then
c       broadcastreal(colatitude,plat/2)
c       broadcastreal(gausswt,plat/2)
c       broadcastint (plonmax,jelemg)
c     and
c       initialize cosine latitude
c         do j=1,plat/2
c            coslat(j)=DCOS(Pi/2-colatitude(j))
c            coslat(plat+1-j)=coslat(j)
c         end do

!--------------------------------------------------------------------
!      Include Files
!--------------------------------------------------------------------
#include <planet.h>
#include <constants.h>
#include <coord.h>
#include <mpcluster.h>
#include <io.h>
!--------------------------------------------------------------------

      INT_TYPE  ounit
      INT_TYPE  iunit

      INT_TYPE  navu
      external  navu

      INT_TYPE  lenchr
      external  lenchr

      INT_TYPE  ierr        ! file io error flag
      character*80 msg      ! file io error message

      INT_TYPE  j,jj,jeg
      REAL_TYPE gauss_sum

      character*1 trunc
      character*80 stub

      character*5 tail,tail_lj
      character*80 grdfile

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (rank.eq.root) then

         iunit=navu()      

         if (iunit.gt.0) then

            ! build grid file name
            stub="GRID"
            trunc='T'
            write(tail,1)mm
 1          format(i4)
            call ljust(tail,tail_lj)
            grdfile=stub(1:lenchr(stub))//'.'//trunc//tail_lj(1:lenchr(tail_lj))

            print *,"reading in ",grdfile,"...."

            open(UNIT   = iunit,
     $           FILE   = grdfile,
     $           FORM   = 'formatted',
     $           STATUS = 'old',
     $           IOSTAT = ierr)

            if (ierr.eq.0) then
               do j=1,plat/2
                  read(iunit,10)jj,colatitude(j),gausswt(j)
 10               format(i4,1x,e27.17,1x,e27.17,1x,i6)
               end do

               do jeg=1,jelemg
                  read(iunit,20)jj,plonmax(jeg)
 20               format(i4,1x,i4)
               end do

               close(iunit)
            else
               ierr=OPEN_FAIL
            end if

         else

            ierr=UNIT_ASGN_FAIL

         end if

      end if

      call broadcastint(ierr,1)
      
      if (ierr.eq.OPEN_FAIL) then

         msg='ABORT IN GRIDFILE: error opening/reading '//grdfile(1:lenchr(grdfile))
         call haltmp(msg)

      else if (ierr.eq.UNIT_ASGN_FAIL) then

         msg='ABORT IN GRIDFILE: could not assign unit for '//grdfile(1:lenchr(grdfile))
         call haltmp(msg)

      else

         call broadcastreal(colatitude,plat/2)
         call broadcastreal(gausswt,plat/2)
         call broadcastint (plonmax,jelemg)

         ! initialize cosine latitude

         do j=1,plat/2
            coslat(j)=DCOS(Pi/2-colatitude(j))
            coslat(plat+1-j)=coslat(j)
         end do

#if 0
         write(MPI_UNIT,9000)
 9000    format("")
         write(MPI_UNIT,9010)
 9010    format("---------------------------grid data---------------")
         write(MPI_UNIT,9020)
 9020    format("   j   Colatitude          Gaussian Wt             ")
         write(MPI_UNIT,9030)
 9030    format("---------------------------------------------------")

         gauss_sum=0.0d0
         do j=1,plat/2
            write(MPI_UNIT,9040)j,colatitude(j),gausswt(j)
 9040       format(i4,2x,f18.15,2x,f18.15,5x,i4)
            gauss_sum=gauss_sum+gausswt(j)
         end do
         write(MPI_UNIT,9050)
 9050    format("---------------------------------------------------")
         write(MPI_UNIT,9060)gauss_sum
 9060    format(" Sum of Gaussian Wts = ",f18.15)
         write(MPI_UNIT,9070)
 9070    format("---------------------------------------------------")
#endif

      end if

      end








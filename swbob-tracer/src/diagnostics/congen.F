c---------------------------------------------------------------------------
c
c            Routine which takes a given gridded input data field
c            and produces contours in the format used in dqgcs.F
c            (i.e. like in pvcont.dat)
c
c            Adapted from dqgcs.F on 11 Dec 2003 by D.G. Dritschel
c            @ St Andrews
c
c            Adapted for BOB on 31 Dec 2003 by R.K. Scott @ St Andrews
c
c            Called by sf_diag
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine congen(pvgrid,icountp)
c
#include <congen.h>
#include <io.h>
#include <coord.h>
c
      parameter (ccurv=one)
c
c      ccurv controls the dependence on curvature in renode;
c      *** use ccurv = 0 only for very noisy gridded data.
c
c
      real*4 pvgrid(plon,plat,plev)
      real*4 pvgridr(plon,plat)
      real*8 yg2(nh+1)
      real*8 rlat(plat)
c
c
c      Define fixed arrays and constants:
      call init
c
c---------------------------------------------------------------------
c      Process each layer in turn:
      do lz=1,plev
c
c     reverse latitude grid to go from SP to NP:
         do i=1,plon
            do j=1,plat
               pvgridr(i,j)=pvgrid(i,plat+1-j,lz)
            enddo
         enddo
         do iy=1,nh+1
            yg2(iy)=yg(iy)/two
         enddo
         do j=1,plat/2
            rlat(j)=colatitude(j)-Pi/2.d0
            rlat(plat+1-j)=-rlat(j)
         enddo
c
c     For each ix, interpolate onto regular iy grid:
         nh2=nh/2
         qamax=-1e6
         qamin=1e6
         do ix=1,nh
            if (ix.le.nh2) ix180=ix+nh2
            if (ix.gt.nh2) ix180=ix-nh2
            do iy=1,nh+1
               ixm=ix
               ixp=ix
               do j=1,plat
                  if(rlat(j).gt.yg2(iy)) then
                     jp=j
                     jm=j-1
                     if(jp.eq.1) then
                        jm=jp
                        ixm=ix180
                     endif
                     goto 10
                  endif
                  if(j.eq.plat) then
                     jp=plat
                     jm=plat
                     ixp=ix180
                  endif
               enddo
 10            continue
               if(ixp.eq.ixm) then
                  ap=rlat(jp)-yg2(iy)
                  am=yg2(iy)-rlat(jm)
               else if (ixp.eq.ix) then
                  ap=rlat(jp)-yg2(iy)
                  am=yg2(iy)+rlat(jm)+Pi
               else if (ixm.eq.ix) then
                  ap=Pi-(rlat(jp)+yg2(iy))
                  am=yg2(iy)-rlat(jm)
               endif
               qa(iy,ix)=(am*pvgridr(ixp,jp)+ap*pvgridr(ixm,jm))/(ap+am)
               qamax=max(qa(iy,ix),qamax)
               qamin=min(qa(iy,ix),qamin)
            enddo
         enddo
        qapolen=0.d0
        qapoles=0.d0
        do ix=1,nh
           qapolen=qapolen+qa(nh+1,ix)/nh
           qapoles=qapoles+qa(1,ix)/nh
        enddo
        do ix=1,nh
           qa(nh+1,ix)=qapolen
           qa(1,ix)=qapoles
        enddo
        do iy=1,nh+1
           qa(iy,nh+1)=qa(iy,1)
        enddo

        print*,'qamax,qamin',qamax,qamin
c
c        Generate contours (xd,yd):
        call grid2con(ccurv,icountp)
c
c        Store the contours (xd,yd) into (xa,ya):
        call storecon(lz)
c
      enddo
c
c      Write contours in congen.dat:
      call writecon(icountp)
c
      end
c
c=======================================================================
      subroutine init
c
c      Defines various fixed arrays and constants.
c
#include <congen.h>
c
      dqi=one/dq
      qoff=dq*dble(nlevm)
c      qoff: should be a large integer multiple of the 
c            contour interval, dq.  The multiple should exceed 
c            the maximum expected number of contour levels.
c
      do lev=1,2*nlevm
        qlev(lev)=(dble(lev)-half)*dq-qoff
      enddo
c

c--------------------------------------
c      Grid lengths and their inverses:
      gl=twopi/dble(nh)
      gli=one/gl
c      gl: the grid length
c
c--------------------------------------------------
c      Surgical and node redistribution parameters:
      dm=gl/10.d0
      amu=min(two/sqrt(dble(nh)),0.2d0)
      ell=four*dm/amu**2
c      write(*,*)
c      write(*,'(2(a,f9.7))') ' Using mu = ',amu,' and L = ',ell
c
      dmsq=four*dm**2
      dmi=two/dm
      elf=one/ell**2
      densf=one/(amu*sqrt(ell))
c
c-----------------------
c      Used in GRID2CON:
      do ix=1,nh+1
        xg(ix)=gl*dble(ix-1)-pi
      enddo
      do iy=1,nh+1
        yg(iy)=gl*dble(iy-1)-pi
      enddo
c
      do ix=1,nh
        ibx(ix,1)=nh*(ix-1)
      enddo
      do ix=2,nh
        ibx(ix,0)=ibx(ix-1,1)
      enddo
      ibx(1,0)=ibx(nh,1)
c
      do iy=1,nh
        iby(iy,1)=iy
      enddo
      do iy=2,nh
        iby(iy,0)=iby(iy-1,1)
      enddo
      iby(1,0)=iby(nh,1)
c
      na=0
      npta=0
c
      return
      end
c
c=======================================================================
      subroutine grid2con(ccurv,icountp)
c
c      Generates new contours (xd,yd) from the gridded data (qa) 
c      on the horizontal grid (nh,nh).
c
c      ccurv controls the dependence on curvature in renode.
c
#include <congen.h>
c
c      Local Grid -> Contour arrays:
      dimension ycr(ncrm),xcr(ncrm)
      dimension lcr(ncrm),kib(ncrm),kob(ncrm),ipo(ncrm)
      dimension lgr(nh+1),isi(nh)
      dimension icrtab(nh*nh,2)
      dimension ncrpl(2*nlevm),ibpl(2*nlevm)
      integer*1 noctab(nh*nh)
      logical cro(nh),free(ncrm)
c
c      Counter for total number of grid line crossings:
      ncr=0
c      Counters for each level:
      do lev=1,2*nlevm
        ncrpl(lev)=0
      enddo
c
c      Find x grid line crossings first:
      do ix=1,nh
        xgt=xg(ix)
c
        do iy=1,nh+1
          lgr(iy)=nint((qoff+qa(iy,ix))*dqi)
        enddo
c
        do iy=1,nh
          lgrd=lgr(iy+1)-lgr(iy)
          cro(iy)=(lgrd .ne. 0)
          isi(iy)=sign(1,lgrd)
        enddo
c
        do iy=1,nh
          if (cro(iy)) then
            lgrt=lgr(iy)
            jump=isi(iy)
            inc=(1+jump)/2
            kibt=iy+ibx(ix,inc)
            kobt=iy+ibx(ix,1-inc)
            fac=gl/(qa(iy+1,ix)-qa(iy,ix))
            do while (lgrt .ne. lgr(iy+1))
              ncr=ncr+1
              xcr(ncr)=xgt
              lev=lgrt+inc
              lcr(ncr)=lev
              kib(ncr)=kibt
              kob(ncr)=kobt
              ycr(ncr)=yg(iy)+fac*(qlev(lev)-qa(iy,ix))
              ncrpl(lev)=ncrpl(lev)+1
              lgrt=lgrt+jump
            enddo
          endif
        enddo
c
      enddo
c
c      Above, kib = grid box into which the contour (containing icr) is going
c             kob =   "   "  out of "    "     "         "       "    " coming
c        [kob -> icr -> kib:  icr lies at the boundary between kob & kib]
c
c      Find y grid line crossings next:
      do iy=1,nh
        ygt=yg(iy)
c
        do ix=1,nh+1
          lgr(ix)=nint((qoff+qa(iy,ix))*dqi)
        enddo
c
        do ix=1,nh
          lgrd=lgr(ix+1)-lgr(ix)
          cro(ix)=(lgrd .ne. 0)
          isi(ix)=sign(1,lgrd)
        enddo
c
        do ix=1,nh
          if (cro(ix)) then
            lgrt=lgr(ix)
            jump=isi(ix)
            inc=(1+jump)/2
            kibt=ibx(ix,1)+iby(iy,1-inc)
            kobt=ibx(ix,1)+iby(iy,inc)
            fac=gl/(qa(iy,ix+1)-qa(iy,ix))
            do while (lgrt .ne. lgr(ix+1))
              ncr=ncr+1
              ycr(ncr)=ygt
              lev=lgrt+inc
              lcr(ncr)=lev
              kib(ncr)=kibt
              kob(ncr)=kobt
              xcr(ncr)=xg(ix)+fac*(qlev(lev)-qa(iy,ix))
              ncrpl(lev)=ncrpl(lev)+1
              lgrt=lgrt+jump
            enddo
          endif
        enddo
c
      enddo
c
c      Find the first and last level having non-zero number of crossings:
      lev=1
      do while (ncrpl(lev) .eq. 0) 
        lev=lev+1
      enddo
      levbeg=lev
      do while (ncrpl(lev) .ne. 0) 
        lev=lev+1
      enddo
      levend=lev-1
c
      ibpl(levbeg)=1
      do lev=levbeg+1,levend
        ibpl(lev)=ibpl(lev-1)+ncrpl(lev-1)
      enddo
c
      do icr=1,ncr
        lev=lcr(icr)
        ibt=ibpl(lev)
        ipo(ibt)=icr
        ibpl(lev)=ibt+1
      enddo
c
      do lev=levbeg,levend
        ibeg=ibpl(lev)-ncrpl(lev)
        iend=ibpl(lev)-1
c
c        Initialise number of crossings per box:
        do i=ibeg,iend
          noctab(kob(ipo(i)))=0
        enddo
c
        do i=ibeg,iend
          icr=ipo(i)
c          icr is the index of the current crossing at level lev.
          k=kob(icr)
c          accumulate number of crossings in this box:
          noctab(k)=noctab(k)+1
c          assign crossing to box, permitting 2 crossings:
          icrtab(k,noctab(k))=icr
        enddo
c
        do i=ibeg,iend
          icr=ipo(i)
          k=kib(icr)
          noc=noctab(k)
c          Use last crossing in this box as the next node:
          kob(icr)=icrtab(k,noc)
c          kob(icr) now gives the next point after icr
          noctab(k)=noc-1
c          This will normally be zero, except for boxes with 2 crossings;
c          this allows a second use of this box.
        enddo
c
      enddo
c
c-----------------
c      Now re-build contours:
      j=0
      i=0
c
      do icr=1,ncr
        free(icr)=.true.
      enddo
c
      do icr=1,ncr
        if (free(icr)) then
c          A new contour (j) starts here:
          i=i+1
          j=j+1
          i1d(j)=i
          indd(j)=lcr(icr)
          ipo(i)=icr
          icrn=kob(icr)
c
          do while (icrn .ne. icr)
c            Find remaining points on contour j:
            i=i+1
            ipo(i)=icrn
            free(icrn)=.false.
            icrn=kob(icrn)
          enddo
          i2d(j)=i
          npd(j)=i2d(j)-i1d(j)+1
        endif
      enddo
c
      nptd=0
      nrem=0
      ndt=j
c
      do j=1,ndt
        if (npd(j) .lt. 5) then
          nrem=nrem+1
        else
          ibeg=nptd+1
          ioff=ibeg-i1d(j)
          do i=i1d(j),i2d(j)
            icr=ipo(i)
            xd(ioff+i)=xcr(icr)
            yd(ioff+i)=ycr(icr)
          enddo
          nd=j-nrem
          npd(nd)=npd(j)
          indd(nd)=indd(j)
          i1d(nd)=ibeg
          nptd=nptd+npd(j)
          i2d(nd)=nptd
        endif
      enddo
c
c      Done rebuilding contours.
c
c-------------------------------------------------------
c      Remove points that are extremely close together:
      do i=1,nptd-1
        xx=xd(i+1)-xd(i)
        dx(i)=xx-twopi*int(xx*hlxi)
        yy=yd(i+1)-yd(i)
        dy(i)=yy-twopi*int(yy*hlyi)
      enddo
c
      do j=1,nd
        ibeg=i1d(j)
        iend=i2d(j)
        xx=xd(ibeg)-xd(iend)
        dx(iend)=xx-twopi*int(xx*hlxi)
        yy=yd(ibeg)-yd(iend)
        dy(iend)=yy-twopi*int(yy*hlyi)
      enddo
c
      do i=1,nptd
        v(i)=dx(i)**2+dy(i)**2
      enddo
c
      ndt=nd
      nd=0
      nptd=0
c
      do j=1,ndt
        nptbeg=nptd
        do i=i1d(j),i2d(j)
          if (v(i) .gt. small) then
            nptd=nptd+1
            xd(nptd)=xd(i)
            yd(nptd)=yd(i)
          endif
        enddo
        npdiff=nptd-nptbeg
        if (npdiff .lt. 5) then
          nptd=nptbeg
        else
          nd=nd+1
          i1d(nd)=nptbeg+1
          i2d(nd)=nptd
          npd(nd)=npdiff
        endif
      enddo
c
c------------------------------------------------------------------------
c      Redistribute points on the contours:
      call renode(ccurv)
c
      return
      end
c
c=======================================================================
      subroutine storecon(lz)
c
c      Stores the shifted contours xd, yd, i1d, i2d ... 
c                  into the arrays xa, ya, i1a, i2a ....  
c      Also updates na & npta and the node & contours layer indices 
c      il1a, il2a, jl1a & jl2a.
c
#include <congen.h>
c
c-----------------------------------------------------
c      Store starting layer indices:
      jl1a(lz)=na+1
      il1a(lz)=npta+1
c
c      Store nodes:
      do i=1,nptd
        xa(npta+i)=xd(i)
        ya(npta+i)=yd(i)
      enddo
c
c      Store contour indices:
      do j=1,nd
        jn=na+j
        i1a(jn) =npta+i1d(j)
        i2a(jn) =npta+i2d(j)
        npa(jn) =npd(j)
        inda(jn)=indd(j)
        laya(jn)=lz
      enddo
c
c      Augment number of contours and nodes:
      na=na+nd
      npta=npta+nptd
c
c      Store ending layer indices:
      jl2a(lz)=na
      il2a(lz)=npta
c
      return
      end
c
c=======================================================================
      subroutine renode(ccurv)
c      Re-nodes each contour while preserving corner locations.
c      Uses square-root dependence on a weighted sum of nearby 
c      curvature values with an additional term proportional
c      to ccurv*curvature.  It is recommended that ccurv = 0
c      when this is called after grid2con, since the latter
c      generates relatively noisy contours.
c
#include <congen.h>
      logical last
c
c------------------------------------------------------------------------
c      Set up next array:
c
      do i=1,nptd-1
        next(i)=i+1
      enddo
c
      do j=1,nd
        next(i2d(j))=i1d(j)
      enddo
c
c------------------------------------------------------------------------
c      Get the updated cubic interpolation coefficients:
      call cubic
c
c------------------------------------------------------------------------
c      Use the spherical curvature expression (radius of the sphere = ell)
c      to ensure an adequate node density in low curvature regions.
      do i=1,nptd
        ww=one/(v(i)+dmsq)
        u(i)=ww*sqrt(elf*v(i)+u(i)**2)
        v(i)=ww*d(i)
      enddo
c      NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; d(i)=sqrt{v(i)};
c          u(i)/d(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2
c
c      Re-assign curvature at a node from weighted average on either side
c      (v above is the weight):
      do ib=1,nptd
        i=next(ib)
        q(i)=(u(ib)+u(i))/(v(ib)+v(i))
      enddo
c
c      Re-average to get interval value (effectively, four curvature
c      values go into getting the final interval value, u(i)):
      do i=1,nptd
        ia=next(i)
        u(i)=half*(q(i)+q(ia))
      enddo
c
c      Compute fractional number of nodes to be placed between old
c      nodes i and i+1:
      do i=1,nptd
        d(i)=d(i)*min(dmi,densf*sqrt(u(i))+ccurv*u(i))
      enddo
c      NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})
c
c------------------------------------------------------------------------
c      Now begin the redistribution of nodes contour by contour,
c      making sure to preserve corner locations:
      nptd=0
      do j=1,nd
        inew=1
        i1t=i1d(j)
        i1d(j)=nptd+1
300       u(nptd+inew)=xd(i1t)
          v(nptd+inew)=yd(i1t)
          sum=zero
          i=i1t
310         sum=sum+d(i)
            i=i+1
            last=i .gt. i2d(j)
            if (last) goto 330
            if (corn(i)) goto 320
            goto 310
320       if (sum .lt. small) then
            i1t=i
            goto 300
          else
            i2t=i-1
            goto 340
          endif
330       if (sum .lt. small) then
            inew=inew-1
            goto 390
          else
            i2t=i-1
          endif
340       npseg=nint(sum)+1
c          npseg-1 is the number of nodes to be placed on the contour segment.
          fac=dble(npseg)/sum
          do i=i1t,i2t
            d(i)=fac*d(i)
          enddo
c          Now, the sum of d(i) is equal to npseg.
c          The first node along a contour (segment) is fixed;
c          find the new node positions:
          acc=zero
          i=i1t-1
          do im=nptd+inew+1,nptd+inew+npseg-1
            if (acc .ge. one) goto 370
360           acc=acc+d(i+1)
              i=i+1
              if (acc .lt. one) goto 360
370         acc=acc-one
            p=one-acc/d(i)
            eta=p*(a(i)+p*(b(i)+p*c(i)))
            u(im)=xd(i)+p*dx(i)-eta*dy(i)
            v(im)=yd(i)+p*dy(i)+eta*dx(i)
          enddo
          if (last) then
            inew=inew+npseg-1
            goto 390
          else
            inew=inew+npseg
            i1t=i2t+1
            goto 300
          endif
390     npd(j)=inew
        nptd=nptd+inew
      enddo
c
c------------------------------------------------------------------------
c      Switch arrays around again:
      do i=1,nptd
        xd(i)=u(i)-twopi*int(u(i)*hlxi)
        yd(i)=v(i)-twopi*int(v(i)*hlyi)
      enddo
c
c      Reset ending contour indices:
      do j=1,nd
        i2d(j)=i1d(j)+npd(j)-1
      enddo
c
      return
      end
c
c=======================================================================
      subroutine cubic
c      Calculates the interpolation coefficients between the nodes 
c      [xd(i),yd(i)] and [xd(next(i)),yd(next(i))], i = 1, ..., nptd.
c
c      The interpolation approximately enforces continuity of curvature 
c      (except at corners which have effectively infinite curvature).
c
#include <congen.h>
c
c----------------------------------------------------------------------
      do i=1,nptd
        ia=next(i)
        xx=xd(ia)-xd(i)
        dx(i)=xx-twopi*int(xx*hlxi)
        yy=yd(ia)-yd(i)
        dy(i)=yy-twopi*int(yy*hlyi)
        v(i)=dx(i)*dx(i)+dy(i)*dy(i)+small
        d(i)=sqrt(v(i))
      enddo
c 
      do ib=1,nptd
        i=next(ib)
        u(i)=v(ib)
        c(i)=-dx(ib)
        q(i)=-dy(ib)
      enddo
c 
      do i=1,nptd
        corn(i)=dx(i)*c(i)+dy(i)*q(i) .gt. zero
        if (corn(i)) then
c      Set curvature to zero at corners:
          b(i)=zero
          q(i)=zero
        else
          b(i)=(dx(i)*q(i)-c(i)*dy(i))/
     .     sqrt((c(i)*v(i)-dx(i)*u(i))**2+
     .          (q(i)*v(i)-dy(i)*u(i))**2+small3)
          q(i)=one
        endif
      enddo
c
      do i=1,nptd
        ia=next(i)
        u(i)=d(i)*(b(ia)+b(i))
        c(i)=d(i)*(b(ia)-b(i))
      enddo
c
c      Calculate the cubic interpolation coefficients:
      do i=1,nptd
        a(i)=sixth*c(i)-half*u(i)
        b(i)=half*(u(i)-c(i))
        c(i)=third*c(i)
      enddo
c
      return
      end
c
c=======================================================================
      subroutine writecon(icountp)
c      Writes contours to congen.dat and contour levels to levels.dat:
c
#include <congen.h>
      character*7 ofile
      character*5  tail
      INT_TYPE  ounit

      INT_TYPE  navu
      external  navu

c
c----------------------------------------------------------------------
c      Adjust the levels so that - values have negative indices
c                            and + values have positive indices:
      do j=1,na
        lev=inda(j)
        inda(j)=lev-nlevm+(lev-1)/nlevm-1
      enddo
c
c----------------------------------------------------------------------
c
      write(tail,'(i5.5)') icountp
      ofile="c."//tail
      ounit=navu()
      open(ounit,file=ofile,form='formatted',status='unknown')

      write(ounit,'(i6,1x,i7,1x,f7.2)') na,npta,real(icountp)
      do j=1,na
        write(ounit,'(i5,1x,i7,2(1x,i4))') npa(j),i1a(j),inda(j),laya(j)
      enddo

c      do j=1,na
c        i=i1a(j)
c        rlat=ya(i)/two
c        write(*,*) sin(rlat)
c      enddo
c      stop

      do i=1,npta
c        write(33,'(f8.5,1x,f8.5)') xa(i),ya(i)
        rlat=ya(i)/two
        rlon=xa(i)

c     spherical coordinates (for use with gpackage contouring?):
c        r=cos(rlat)
c        write(33,'(f9.6,1x,f9.6,1x,f9.6)') r*cos(rlon),r*sin(rlon),
c     .                                     sin(rlat)

c     Lambert equal area (for plotting and mlm diagnostics with idl):
        r=two*sin((pi/two-rlat)/two)

        write(ounit,'(f9.6,1x,f9.6)') r*cos(rlon),r*sin(rlon)

      enddo

      close(ounit)
c
      return
      end

      subroutine accumulate(iflag)
      implicit none

      INT_TYPE iflag,icount,nacc


      nfor=(nf1+nf0)/2
      do i=1,4
         Eglb(i)=0.d0
         Zglb(i)=0.d0
         Ecnt(i)=0.d0
         Enrm(i)=0.d0
         do n=0,nn
            Eglb(i)=Eglb(i)+Etot(n,i)
            Zglb(i)=Zglb(i)+Ztot(n,i)
#ifdef GAUSS
            Ecnt(i)=Ecnt(i)+Etot(n,i)*float(n)
            Enrm(i)=Enrm(i)+Etot(n,i)
#else
            if ((abs(n-nfor)).gt.(nf1-nf0)) Ecnt(i)=Ecnt(i)+Etot(n,i)*float(n)
            if ((abs(n-nfor)).gt.(nf1-nf0)) Enrm(i)=Enrm(i)+Etot(n,i)
#endif
            Erotacc(n,i)=Erotacc(n,i)+Erot(n,i)
            Edivacc(n,i)=Edivacc(n,i)+Ediv(n,i)
            Etotacc(n,i)=Etotacc(n,i)+Etot(n,i)
            Ztotacc(n,i)=Ztotacc(n,i)+Ztot(n,i)
         enddo
         Eglbacc(i)=Eglbacc(i)+Eglb(i)
         Zglbacc(i)=Zglbacc(i)+Zglb(i)
         Urmsacc(i)=Urmsacc(i)+sqrt(2.0*Eglb(i))
         Ecntacc(i)=Ecntacc(i)+Ecnt(i)/Enrm(i)
      enddo
      nacc=nacc+1


      if (iflag.eq.0) then 
         do i=1,4
            Urmsacc(i)=0.d0
            Eglbacc(i)=0.d0
            Zglbacc(i)=0.d0
            Ecntacc(i)=0.d0
            do n=0,nn
               Erotacc(n,i)=0.d0
               Edivacc(n,i)=0.d0
               Etotacc(n,i)=0.d0
               Ztotacc(n,i)=0.d0
            enddo
         enddo
         nacc=0
      endif

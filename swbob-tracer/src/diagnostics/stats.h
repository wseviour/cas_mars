#ifndef STATS_H_
#define STATS_H_

      REAL_TYPE l1num(plev)
      REAL_TYPE l2num(plev)
      REAL_TYPE linfnum(plev)

      REAL_TYPE l1den(plev)
      REAL_TYPE l2den(plev)
      REAL_TYPE linfden(plev)

      REAL_TYPE l1(plev)
      REAL_TYPE l2(plev)
      REAL_TYPE linf(plev)

      common /l_norms/ l1num,l2num,linfnum,
     $                 l1den,l2den,linfden,
     $                 l1,   l2,   linf

      REAL_TYPE phiave(plev)
      REAL_TYPE phitave(plev)
      REAL_TYPE phi0ave(plev)

      REAL_TYPE phimean(plev)
      REAL_TYPE phivar(plev)

      REAL_TYPE dphisq(plev)
      REAL_TYPE dphitsq(plev)
      REAL_TYPE dphi0sq(plev)

      common /phi_moments/ phiave,phitave,phi0ave,
     $                     phimean,phivar,
     $                     dphisq,dphitsq,dphi0sq
      
#endif

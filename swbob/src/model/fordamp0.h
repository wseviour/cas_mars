#ifndef FORDAMP0_H_
#define FORDAMP0_H_


      REAL_TYPE kv
      common /fdrag/ kv

      REAL_TYPE latj,dlat,polvor,kt,usb
      common /frlax/ latj,dlat,polvor,kt,usb

      REAL_TYPE latt,dlatt,wavet,ampt,freqt
      REAL_TYPE lath,dlath,waveh,amph,freqh
      REAL_TYPE latm,dlatm,wavem,ampm
      common /ftopo/ latt,dlatt,wavet,ampt,freqt
      common /fheat/ lath,dlath,waveh,amph,freqh
      common /fmech/ latm,dlatm,wavem,ampm

      REAL_TYPE eps0,gamma,fcor0,fcor1
      common /fspec0/ eps0,gamma,fcor0,fcor1

      INT_TYPE nf0,nf1,irl
      common /fspec1/ nf0,nf1,irl




#endif



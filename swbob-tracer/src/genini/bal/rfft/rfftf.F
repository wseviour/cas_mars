#include <type.h>
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c                            Copyright (C) 1995   
c 
c              University Corporation for Atmospheric Research
c 
c                           All rights reserved.
c 
c      No part of this  code may  be reproduced,  stored in  a retrieval
c      system,  translated,  transcribed,  transmitted or distributed in
c      any form or by any means, manual, electric, electronic, chemical,
c      electo-magnetic, mechanical, optical, photocopying, recording, or
c      otherwise, without the prior  explicit written permission  of the
c      author(s)  or their  designated proxies.  In no event  shall this
c      copyright notice be removed or altered in any way.
c 
c      This code is provided "as is",  without any warranty of any kind,
c      either  expressed or implied,  including but not limited to,  any
c      implied warranty of  merchantibility or fitness  for any purpose.
c      In no event will any party who distributed the code be liable for
c      damages or for any claim(s) by any other party, including but not
c      limited to,  any lost  profits,  lost monies,  lost data  or data
c      rendered inaccurate,  losses sustained  by third parties,  or any
c      other special, incidental or consequential damages arising out of
c      the use or inability to use the program,  even if the possibility
c      of such damages has  been advised against.  The entire risk as to
c      the quality,  the performance, and the fitness of the program for
c      any particular purpose lies with the party using the code.
c 
c      *****************************************************************
c      ANY USE OF  THIS CODE CONSTITUES  ACCEPTANCE OF  THE TERMS OF THE
c                               ABOVE STATEMENTS
c      *****************************************************************
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     SUBROUTINE RFFTF(N,R,WSAVE)
C
C     SUBROUTINE RFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
C     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
C     BELOW AT OUTPUT PARAMETER R.
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
C
C     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C             TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C             IN THE PROGRAM THAT CALLS RFFTF. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
C
C
C     OUTPUT PARAMETERS
C
C     R       R(1) = THE SUM FROM I=1 TO I=N OF R(I)
C
C             IF N IS EVEN SET L =N/2   , IF N IS ODD SET L = (N+1)/2
C
C               THEN FOR K = 2,...,L
C
C                  R(2*K-2) = THE SUM FROM I = 1 TO I = N OF
C
C                       R(I)*COS((K-1)*(I-1)*2*PI/N)
C
C                  R(2*K-1) = THE SUM FROM I = 1 TO I = N OF
C
C                      -R(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C             IF N IS EVEN
C
C                  R(N) = THE SUM FROM I = 1 TO I = N OF
C
C                       (-1)**(I-1)*R(I)
C
C      *****  NOTE
C                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
C                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
C                  SEQUENCE BY N.
C
C     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
C             CALLS OF RFFTF OR RFFTB.
C
      SUBROUTINE RFFTF (N,R,WSAVE)
      IMPLICIT NONE

      REAL_TYPE  R(*),WSAVE(*)
      INT_TYPE N
      REAL_TYPE  RN

      INT_TYPE I
C
      IF (N .GT. 1) THEN
         RN=1.0D0/N
         DO I=1,N
            R(I)=R(I)*RN
         END DO
         CALL RFFTF1(N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      END IF

      RETURN
      END

#include <type.h>
      SUBROUTINE RFFTF1 (N,C,CH,WA,IFAC)
      REAL_TYPE  CH(*),C(*),WA(*)
      INT_TYPE   N,IFAC(*)

      INT_TYPE NA,NF,L1,L2,IW,KH,K1,IP,IDO,IDL1
      INT_TYPE IX2,IX3,IX4
      INT_TYPE I

      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO K1=1,NF
         KH = NF-K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .EQ. 4) THEN
            IX2 = IW+IDO
            IX3 = IX2+IDO
            IF (NA .EQ. 0) THEN
               CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
            ELSE
               CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
            END IF
         ELSE IF (IP .EQ. 2) THEN
            IF (NA .EQ. 0) THEN
               CALL RADF2 (IDO,L1,C,CH,WA(IW))
            ELSE
               CALL RADF2 (IDO,L1,CH,C,WA(IW))
            END IF
         ELSE IF (IP .EQ. 3) THEN
            IX2 = IW+IDO
            IF (NA .EQ. 0) THEN
               CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
            ELSE
               CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
            END IF
         ELSE IF (IP .EQ. 5) THEN
            IX2 = IW+IDO
            IX3 = IX2+IDO
            IX4 = IX3+IDO
            IF (NA .EQ. 0) THEN
               CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            ELSE
               CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            ENDIF
         ELSE 
            IF (IDO .EQ. 1) NA = 1-NA
            IF (NA .EQ. 0) THEN
               CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
               NA = 1
            ELSE
               CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
               NA = 0
            END IF
         END IF
         L2 = L1
      END DO

      IF (NA .EQ. 0) THEN
         DO I=1,N
            C(I) = CH(I)
         END DO
      END IF

      END

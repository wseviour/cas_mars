#include <type.h>
      SUBROUTINE RFFTB1 (N,C,CH,WA,IFAC)
      IMPLICIT NONE

      REAL_TYPE  CH(*)      ,C(*)       ,WA(*)      
      INT_TYPE IFAC(*)
      INT_TYPE N

      INT_TYPE NA,NF,L1,L2,IW,K1,IP,IDO,IDL1
      INT_TYPE IX2,IX3,IX4
      INT_TYPE I

      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .EQ. 4) THEN
            IX2 = IW+IDO
            IX3 = IX2+IDO
            IF (NA .EQ. 0) THEN
               CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
            ELSE
               CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
            END IF
            NA = 1-NA
         ELSE IF (IP .EQ. 2) THEN
            IF (NA .EQ. 0) THEN
               CALL RADB2 (IDO,L1,C,CH,WA(IW))
            ELSE
               CALL RADB2 (IDO,L1,CH,C,WA(IW))
            END IF
            NA = 1-NA
         ELSE IF (IP .EQ. 3) THEN
            IX2 = IW+IDO
            IF (NA .EQ. 0) THEN
               CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
            ELSE
               CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
            ENDIF
            NA = 1-NA
         ELSE IF (IP .EQ. 5) THEN
            IX2 = IW+IDO
            IX3 = IX2+IDO
            IX4 = IX3+IDO
            IF (NA .EQ. 0) THEN
               CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            ELSE
               CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            END IF
            NA = 1-NA
         ELSE
            IF (NA .EQ. 0) THEN
               CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
            ELSE
               CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
            END IF
            IF (IDO .EQ. 1) NA = 1-NA
         END IF
         L1 = L2
         IW = IW+(IP-1)*IDO
      END DO

      IF (NA.EQ.1) THEN
         DO I=1,N
            C(I) = CH(I)
         END DO
      END IF

      RETURN
      END

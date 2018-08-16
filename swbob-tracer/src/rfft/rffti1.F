#include <type.h>
      SUBROUTINE RFFTI1 (N,WA,IFAC)
      IMPLICIT NONE

      REAL_TYPE   WA(*)      
      INT_TYPE    IFAC(*)
      INT_TYPE    N

      INT_TYPE    NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/

      INT_TYPE NL,NF,NFM1,NTRY,NQ,NR,IB,IS
      INT_TYPE I,J,K1,II,L1,L2,LD,IP,IPM,IDO

      REAL_TYPE TPI,ARG,ARGH,ARGLD,FI

      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF ((NTRY .NE. 2).OR.(NF .EQ. 1)) THEN
         IF (NL .NE. 1) GO TO 104
      ELSE
         DO I=2,NF
            IB = NF-I+2
            IFAC(IB+2) = IFAC(IB+1)
         END DO
         IFAC(3) = 2
      END IF

      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0D0*DATAN(1.0D0)
      ARGH = TPI/REAL(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1

      IF (NFM1 .NE. 0) THEN
         DO K1=1,NFM1
            IP = IFAC(K1+2)
            LD = 0
            L2 = L1*IP
            IDO = N/L2
            IPM = IP-1
            DO J=1,IPM
               LD = LD+L1
               I = IS
               ARGLD = REAL(LD)*ARGH
               FI = 0.0d0
               DO II=3,IDO,2
                  I = I+2
                  FI = FI+1.0d0
                  ARG = FI*ARGLD
                  WA(I-1) = DCOS(ARG)
                  WA(I  ) = DSIN(ARG)
               END DO
               IS = IS+IDO
            END DO
            L1 = L2
         END DO
      END IF

      END

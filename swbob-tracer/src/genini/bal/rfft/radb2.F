#include <type.h>
      SUBROUTINE RADB2 (IDO,L1,CC,CH,WA1)
      IMPLICIT NONE

      INT_TYPE   IDO,L1
      REAL_TYPE  CC(IDO,2,L1), 
     1           CH(IDO,L1,2),
     2           WA1(*)

      INT_TYPE I,K,IC,IDP2
      REAL_TYPE TI2,TR2

      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE

      IF (IDO-2) 107,105,102

  102 IDP2 = IDO+2

      DO K=1,L1
         DO I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            CH(I  ,K,1) = CC(I  ,1,K)-CC(IC  ,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            TI2 = CC(I  ,1,K)+CC(IC  ,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I  ,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
         END DO
      END DO

      IF (MOD(IDO,2) .EQ. 1) RETURN

  105 DO K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
      END DO

  107 RETURN
      END


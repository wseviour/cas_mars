#include <type.h>
      SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)
      IMPLICIT NONE

      INT_TYPE   IDO,L1

      REAL_TYPE  CH(IDO,2,L1),
     1           CC(IDO,L1,2),
     2           WA1(*)

      INT_TYPE I,K,IC,IDP2
      REAL_TYPE TI2,TR2

      DO 101 K=1,L1
         CH(1  ,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE

      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO K=1,L1
         DO I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I  ,K,2)
            TI2 = WA1(I-2)*CC(I  ,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I  ,1,K) = CC(I  ,K,1)+TI2
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC  ,2,K) = TI2-CC(I,K,1)
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
         END DO
      END DO

      IF (MOD(IDO,2) .EQ. 1) RETURN

  105 DO K=1,L1
         CH(1  ,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) =  CC(IDO,K,1)
      END DO

  107 RETURN
      END

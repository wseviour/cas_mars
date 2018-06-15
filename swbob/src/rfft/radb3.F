#include <type.h>
      SUBROUTINE RADB3 (IDO,L1,CC,CH,WA1,WA2)
      IMPLICIT NONE

      INT_TYPE  IDO,L1
      REAL_TYPE CC(IDO,3,L1),
     1          CH(IDO,L1,3),
     2          WA1(*),
     3          WA2(*)

      REAL_TYPE TAUR,TAUI
      DATA TAUR,TAUI /-.500000000000000D0,.866025403784439D0/

      REAL_TYPE TR2,TI2
      REAL_TYPE CR2,CI2
      REAL_TYPE CR3,CI3
      REAL_TYPE DR2,DI2
      REAL_TYPE DR3,DI3

      INT_TYPE I,K,IC,IDP2

      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE

      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO K=1,L1
         DO I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TI2 = CC(I  ,3,K)-CC(IC  ,2,K)
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I  ,3,K)+CC(IC  ,2,K))
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CI2 = CC(I  ,1,K)+TAUR*TI2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            CH(I  ,K,1) = CC(I  ,1,K)+TI2
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I  ,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I  ,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
         END DO
      END DO

      RETURN
      END

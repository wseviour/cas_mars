#include <type.h>
      SUBROUTINE PGAQD(NLATG,THETA,WTS,WORK,LWORK,IERROR)
      IMPLICIT NONE

C ============================================================================
C
C     SUBROUTINE GAQD COMPUTES GAUSSIAN POINTS (IN RADIANS) AND WEIGHTS
C     ON THE SPHERE IN THE INTERVAL (0,PI).  (THESE CAN BE USED IN
C     GAUSSIAN QUADRATURE FOR ESTIMATING INTEGRALS ON THE SPHERE)
C
C     4/8/99: SUBROUTINE PGAQD.
C     THIS VERSION RUNS SIGNIFICANTLY FASTER THAN GAQD. 
C     THIS IS ESPECIALLY USEFUL FOR *REALLY* HIGH RESOLUTIONS (E.G. >T1000)
C
C                                                          ... RICH LOFT
C ============================================================================

      INT_TYPE  LWORK,NLATG
      REAL_TYPE WORK(LWORK),THETA(NLATG),WTS(NLATG)
      INT_TYPE  IERROR

C ============================================================================
C     LOCAL VARIABLES
C ============================================================================

      INT_TYPE  I1,I2,I3
      REAL_TYPE X

C ============================================================================
C     CHECK WORK SPACE LENGTH
C ============================================================================

      IERROR = 1
      IF (LWORK.LT.NLATG*7) RETURN    ! formerly 4*(N*(N+1))+2 which WAY is too big!

      IERROR = 2
      IF (NLATG.LE.0) RETURN

      IERROR = 0

      IF (NLATG.GT.2) THEN

C ============================================================================
C        PARTITION WORK SPACE FOR REAL_TYPE EIGENVALUE(VECTOR COMPUTATION)
C ============================================================================

         I1 = 1
         I2 = I1+2*NLATG
         I3 = I2+2*NLATG

         CALL GAQD1(NLATG,THETA,WTS,WORK(I1),WORK(I2),WORK(I3),IERROR)

         IF (IERROR.NE.0) THEN
            IERROR = 3
            RETURN
         END IF

         RETURN
      ELSE IF (NLATG.EQ.1) THEN

         WTS(1) = 2.0D0
         THETA(1) = DACOS(0.0D0)

      ELSE IF (NLATG.EQ.2) THEN

C ============================================================================
C        COMPUTE WEIGHTS AND POINTS ANALYTICALLY WHEN N=2
C ============================================================================

         WTS(1) = 1.0D0
         WTS(2) = 1.0D0
         X = DSQRT(1.0D0/3.0D0)
         THETA(1) = DACOS(X)
         THETA(2) = DACOS(-X)
         RETURN

      END IF

      END

      SUBROUTINE GAQD1(NLATG,THETA,WTS,W,E,WRK,IER)
      IMPLICIT NONE

C ============================================================================
C
C     SET SYMMETRIC TRIDIAGNONAL MATRIX SUBDIAGONAL AND DIAGONAL
C     COEFFICIENTS FOR MATRIX COMING FROM COEFFICIENTS IN THE
C     RECURSION FORMULA FOR LEGENDRE POLYNOMIALS
C     A(N)*P(N-1)+B(N)*P(N)+C(N)*P(N+1) = 0.
C
C ============================================================================

      INT_TYPE  NLATG
      REAL_TYPE THETA(NLATG)
      REAL_TYPE WTS(NLATG)
      REAL_TYPE W(NLATG)
      REAL_TYPE E(NLATG)
      REAL_TYPE WRK(NLATG*3)
      INT_TYPE  IER

C ============================================================================
C     LOCAL VARIABLES
C ============================================================================

      INT_TYPE  MATZ,INDX
      INT_TYPE  I,J,NLATG2
      REAL_TYPE TEMP

      WRK(1)=0.D0
      WRK(NLATG+1) = 0.D0
      W(1)=0.D0
      E(1) = 0.D0

      DO 100 J=2,NLATG
         WRK(J)= (DBLE(J)-1.D0)/DSQRT((2.D0*DBLE(J)-1.D0)*(2.D0*DBLE(J)-3.D0))
         WRK(J+NLATG)=0.D0
         E(J) = WRK(J)
         W(J) = 0.D0
  100 CONTINUE

C ============================================================================
C     COMPUTE EIGENVALUES  OF MATRIX
C ============================================================================

      MATZ = 1
      INDX = 2*NLATG+1

      CALL DRST(NLATG,W,E,MATZ,WRK(INDX),IER)

      IF (IER.NE.0) RETURN

C ============================================================================
C     COMPUTE GAUSSIAN WEIGHTS AND POINTS
C ============================================================================

      DO 101 J=1,NLATG
         THETA(J) = DACOS(W(J))

C ============================================================================
C        SET GAUSSIAN WEIGHTS AS 1ST COMPONENTS OF EIGENVECTORS SQUARED
C ============================================================================

         INDX   = 2*NLATG + J
         WTS(J) = 2.0D0*WRK(INDX)**2
  101 CONTINUE

C ============================================================================
C     REVERSE ORDER OF GAUSSIAN POINTS TO BE
C     MONOTONIC INCREASING IN RADIANS
C ============================================================================

      NLATG2 = NLATG/2

      DO 102 I=1,NLATG2
         TEMP = THETA(I)
         THETA(I) = THETA(NLATG-I+1)
         THETA(NLATG-I+1) = TEMP
  102 CONTINUE

      RETURN
      END

      SUBROUTINE DRST(NLATG,W,E,MATZ,Z,IERR)
      IMPLICIT NONE

      INT_TYPE  NLATG
      REAL_TYPE W(NLATG),E(NLATG),Z(NLATG)
      INT_TYPE  IERR

C ============================================================================
C
C     DRST IS A REAL_TYPE MODIFICATION OF RST OFF EISPACK
C     TO BE USED  TO COMPUTE GAUSSIAN POINTS AND WEIGHTS
C
C ============================================================================

      INTEGER I,MATZ

C ============================================================================
C
C     FIND BOTH EIGENVALUES AND EIGENVECTORS 
C     ...INITIALIZE Z(I,J) = IDENTITY (ONLY NEED FIRST ROW)
C
C ============================================================================

      Z(1)=1.0D0      
      DO I = 2, NLATG
         Z(I) = 0.0D0
      END DO

      CALL  DINTQL(NLATG,W,E,Z,IERR)

      RETURN
      END

      SUBROUTINE DINTQL(NLATG,D,E,Z,IERR)
      IMPLICIT NONE

C ============================================================================
C     DINTQL IS A REAL_TYPE MODIFICATION OF INTQL2 OFF
C     EISPACK TO BE USED BY GAQD IN SPHEREPACK FOR COMPUTING
C     GAUSSIAN WEIGHTS AND POINTS
C
C ============================================================================

      INT_TYPE NLATG
      REAL_TYPE D(NLATG)
      REAL_TYPE E(NLATG)
      REAL_TYPE Z(NLATG)
      INT_TYPE  IERR

      INTEGER I,J,K,L,M,II,MML
      REAL_TYPE B,C,F,G,P,R,S,TST1,TST2,DPYTHA

      IERR = 0
      IF (NLATG .EQ. 1) GO TO 1001

      DO I = 2, NLATG
         E(I-1) = E(I)
      END DO
      E(NLATG) = 0.0D0

      DO 240 L = 1, NLATG
	 J = 0
C ============================================================================
C        LOOK FOR SMALL SUB-DIAGONAL ELEMENT ...
C ============================================================================

  105    DO 110 M = L, NLATG-1  
	    TST1 = DABS(D(M)) + DABS(D(M+1))
	    TST2 = TST1 + DABS(E(M))
	    IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
	 IF (M .EQ. L) GO TO 240
	 IF (J .EQ. 30) GO TO 1000
	 J = J + 1

C ============================================================================
C        FORM SHIFT ...
C ============================================================================

	 G = (D(L+1) - P) / (2.0D0 * E(L))
	 R = DPYTHA(G,1.0D0)
	 G = D(M) - P + E(L) / (G + SIGN(R,G))
	 S = 1.0D0
	 C = 1.0D0
	 P = 0.0D0
	 MML = M - L

C ============================================================================
C        FOR I=M-1 STEP -1 UNTIL L DO ...
C ============================================================================

         DO 200 I=M-1,L,-1

	    F = S * E(I)
	    B = C * E(I)
	    R = DPYTHA(F,G)
	    E(I+1) = R
	    IF (R .EQ. 0.0D0) GO TO 210
	    S = F / R
	    C = G / R
	    G = D(I+1) - P
	    R = (D(I) - G) * S + 2.0D0 * C * B
	    P = S * R
	    D(I+1) = G + P
	    G = C * R - B

C ============================================================================
C           FORM VECTOR ...
C ============================================================================

            F = Z(I+1)
	    Z(I+1) = S * Z(I) + C * F
	    Z(I) = C * Z(I) - S * F

  200    CONTINUE

	 D(L) = D(L) - P
	 E(L) = G
	 E(M) = 0.0D0
	 GO TO 105
C ============================================================================
C        RECOVER FROM UNDERFLOW ..........
C ============================================================================
  210    D(I+1) = D(I+1) - P
	 E(M) = 0.0D0
	 GO TO 105
  240 CONTINUE

C ============================================================================
C     ORDER EIGENVALUES AND EIGENVECTORS ...
C ============================================================================

      DO 300 II = 2, NLATG
	 I = II - 1
	 K = I
	 P = D(I)

	 DO 260 J = II, NLATG
	    IF (D(J) .GE. P) GO TO 260
	    K = J
	    P = D(J)
  260    CONTINUE

	 IF (K .EQ. I) GO TO 300

	 D(K) = D(I)
	 D(I) = P

         P = Z(I)
	 Z(I) = Z(K)
         Z(K) = P

  300 CONTINUE

      GO TO 1001

C ============================================================================
C     SET ERROR -- NO CONVERGENCE TO AN
C     EIGENVALUE AFTER 30 ITERATIONS ..........
C ============================================================================

 1000 IERR = L
 1001 RETURN
      END

      REAL_TYPE FUNCTION DPYTHA(A,B)
      REAL_TYPE A,B

C ============================================================================
C     DPYTHA IS A REAL_TYPE MODIFICATION OF PYTHAG OFF EISPACK
C     FOR USE BY DIMTQL
C
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
C ============================================================================

      REAL_TYPE P,R,S,T,U

      P = DABS(A)
      IF (DABS(B).GE.DABS(A)) P = DABS(B)
C     P = AMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DABS(A)/P)**2
      IF (DABS(B).LT.DABS(A)) R = (DABS(B)/P)**2
C     R = (AMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
	 T = 4.0D0 + R
	 IF (T .EQ. 4.0D0) GO TO 20
	 S = R/T
	 U = 1.0D0 + 2.0D0*S
	 P = U*P
	 R = (S/U)**2 * R
      GO TO 10
   20 DPYTHA = P

      RETURN
      END



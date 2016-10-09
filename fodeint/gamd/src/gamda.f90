!! THIS MODULE IS PART OF THE CODE GAMD
!! THE CODE GAMD NUMERICALLY SOLVES A (POSSIBLY STIFF)
!! SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS 
!! OR A LINEARLY IMPLICIT DAE
!!
!! Copyright (C)1997-2007   
!! Authors: FRANCESCA MAZZIA (mazzia@dm.uniba.it) 
!!          FELICE IAVERNARO (felix@dm.uniba.it) 
!!
!!
!!This program is free software; you can redistribute it and/or
!!modify it under the terms of the GNU General Public License
!!as published by the Free Software Foundation; either version 2
!!of the License, or (at your option) any later version.
!!
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details.
!!
!!Licensed under The GNU General Public License, Version 2 or later.
!!    http://www.gnu.org/licenses/info/GPLv2orLater.html
!!
!!You should have received a copy of the GNU General Public License
!!along with this program; if not, write to the Free Software
!!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
!!USA.
!!
!!

MODULE PRECISION
      IMPLICIT NONE
      INTEGER,PARAMETER     :: PREC = KIND(1.0d0)
      REAL(PREC), PARAMETER :: EPS  = EPSILON(1.d0)
END MODULE PRECISION
MODULE LINALGGAMD
!!-----------------------------------------------------------------------
!!     ADDITIONAL LINEAR ALGEBRA ROUTINES REQUIRED BY GAMD
!!-----------------------------------------------------------------------
!!     VERSION OF MAY 16, 2003
!!-----------------------------------------------------------------------
!!
       USE PRECISION
       IMPLICIT NONE
!! PARAMETER
       INTEGER :: MLLU, MULU, MDIAG, MBDIAG, MBB, MDIFF
  CONTAINS
    !!
    !!     SUBROUTINE DEC
    !!
SUBROUTINE DEC (N, NDIM, A, IP, IER)
  !! VERSION REAL DOUBLE PRECISION
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,NDIM
  INTEGER, INTENT(IN OUT) :: IP(N), IER
  INTEGER :: NM1,K,KP1,M,I,J
  REAL(PREC), INTENT(IN OUT):: A(NDIM,N)
  REAL(PREC) :: T
  !!-----------------------------------------------------------------------
  !!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
  !!  INPUT..
  !!     N = ORDER OF MATRIX.
  !!     NDIM = DECLARED DIMENSION OF ARRAY  A .
  !!     A = MATRIX TO BE TRIANGULARIZED.
  !!  OUTPUT..
  !!     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
  !!     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
  !!     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
  !!     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
  !!     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
  !!           SINGULAR AT STAGE K.
  !!  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
  !!  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
  !!  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
  !!
  !!  REFERENCE..
  !!     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
  !!     C.A.C.M. 15 (1972), P. 274.
  !!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN

!!----------------------- END OF SUBROUTINE DEC -------------------------
END SUBROUTINE DEC
                  !!
                  !!     SUBROUTINE SOL
                  !!
SUBROUTINE SOL (N, NDIM, A, B, IP)
  !! VERSION REAL DOUBLE PRECISION
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,NDIM
  INTEGER, INTENT(IN) :: IP(N)
  INTEGER :: K,NM1,KP1,M,I,KB,KM1
  REAL(PREC), INTENT(IN) :: A(NDIM,N)
  REAL(PREC), INTENT(IN OUT) :: B(N)
  REAL(PREC) :: T
  !!-----------------------------------------------------------------------
  !!  SOLUTION OF LINEAR SYSTEM, A*X = B .
  !!  INPUT..
  !!    N = ORDER OF MATRIX.
  !!    NDIM = DECLARED DIMENSION OF ARRAY  A .
  !!    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
  !!    B = RIGHT HAND SIDE VECTOR.
  !!    IP = PIVOT VECTOR OBTAINED FROM DEC.
  !!  DO NOT USE IF DEC HAS SET IER .NE. 0.
  !!  OUTPUT..
  !!    B = SOLUTION VECTOR, X .
  !!-----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
!!----------------------- END OF SUBROUTINE SOL -------------------------
END SUBROUTINE SOL
!!
!!     SUBROUTINE DECB
!!
SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, NDIM, ML, MU
  INTEGER, INTENT(OUT) :: IP(N), IER
  REAL(PREC), INTENT(IN OUT) :: A(NDIM,N)
  REAL(PREC) :: T
  INTEGER :: MD, MD1,NM1,M, JU, I, J, K, KP1, MDL, MM, JK, IJK
  !!-----------------------------------------------------------------------
  !!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
  !!  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
  !!  INPUT..
  !!     N       ORDER OF THE ORIGINAL MATRIX A.
  !!     NDIM    DECLARED DIMENSION OF ARRAY  A.
  !!     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS
  !!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
  !!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
  !!                ML+1 THROUGH 2*ML+MU+1 OF  A.
  !!     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!  OUTPUT..
  !!     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
  !!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
  !!     IP      INDEX VECTOR OF PIVOT INDICES.
  !!     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
  !!     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
  !!                SINGULAR AT STAGE K.
  !!  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
  !!  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
  !!  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
  !!
  !!  REFERENCE..
  !!     THIS IS A MODIFICATION OF
  !!     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
  !!     C.A.C.M. 15 (1972), P. 274.
  !!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J)
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
  !!----------------------- END OF SUBROUTINE DECB ------------------------
END SUBROUTINE DECB
                   !!
                   !!     SUBROUTINE SOLB
                   !!
SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N, NDIM, ML, MU
  INTEGER, INTENT(IN) :: IP(N)
  REAL(PREC), INTENT(IN) :: A(NDIM,N)
  REAL(PREC), INTENT(IN OUT) :: B(N)
  REAL(PREC) :: T
  INTEGER :: MD, MD1,NM1,M, JU, I, J, K, KB,MDM,KMD, LM, MDL, IMD

  !!-----------------------------------------------------------------------
  !!  SOLUTION OF LINEAR SYSTEM, A*X = B .
  !!  INPUT..
  !!    N      ORDER OF MATRIX A.
  !!    NDIM   DECLARED DIMENSION OF ARRAY  A .
  !!    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
  !!    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
  !!    B      RIGHT HAND SIDE VECTOR.
  !!    IP     PIVOT VECTOR OBTAINED FROM DECB.
  !!  DO NOT USE IF DECB HAS SET IER .NE. 0.
  !!  OUTPUT..
  !!    B      SOLUTION VECTOR, X .
  !!-----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
!!----------------------- END OF SUBROUTINE SOLB ------------------------
END SUBROUTINE SOLB
 
 
 
 
 
FUNCTION  MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,Y) 
 
  !!----------------------------------------------------------------------- 
  !!  MATRIX VECTOR MULTIPLICATION WITH BANDED MATRIX FMAS  
  !!  INPUT.. 
  !!    R      DECLARED DIMENSION OF ARRAY FMAS. 
  !!    LDMAS  LEADING DIMENSION OF ARRAY FMAS. 
  !!    MLMAS  LOWER BANDWIDTH OF FMAS (DIAGONAL IS NOT COUNTED). 
  !!    MUMAS  UPPER BANDWIDTH OF FMAS (DIAGONAL IS NOT COUNTED). 
  !!  OUTPUT.. 
  !!    MATMULB = FMAS*Y 
  !!----------------------------------------------------------------------- 
 
  USE PRECISION 
  IMPLICIT NONE 
 
  !!   INPUT VARIABLES 
  !!------------------------------------ 
  INTEGER, INTENT(IN) :: R, LDMAS, MLMAS, MUMAS 
  REAL(PREC), INTENT(IN) ::   FMAS(LDMAS,R), Y(R) 
  !! 
  !!   OUTPUT VARIABLES 
  !!------------------------------------ 
  REAL(PREC) :: MATMULB(R) 
  !! 
  !!   LOCAL VARIABLES 
  !!------------------------------------ 
  INTEGER  :: I, J 
  REAL(PREC) ::  S1 
  
  !! 
  !!   EXECUTABLE STATEMENTS 
  !!--------------------------------- 
  !! 
  !!--------- ONE STEP OF THE ITERATION PROCEDURE 
 
      MATMULB(1:R)=0.D0 
      DO I=1,R 
         S1=0.0D0 
         DO J=MAX(1,I-MLMAS),MIN(R,I+MUMAS) 
            S1=S1+FMAS(I-J+MBDIAG,J)*Y(J) 
         END DO 
         MATMULB(I)=MATMULB(I)+S1 
      END DO 
END FUNCTION MATMULB 


END MODULE LINALGGAMD
MODULE SUBGAMD
!!-----------------------------------------------------------------------
!!     ADDITIONAL ROUTINES REQUIRED BY GAMD
!!-----------------------------------------------------------------------
!!     VERSION OF AUGUST 29, 2003
!!-----------------------------------------------------------------------
!!
       USE PRECISION
       USE LINALGGAMD
       IMPLICIT NONE

!!  ORDER 3

      REAL(PREC), PARAMETER ::    &
     &            B311 = 5d0/12d0, &
     &            B312 = 8d0/12d0, &
     &            B313 = -1d0/12d0

       REAL(PREC), PARAMETER ::  Cp31 = 1d0/24d0
       REAL(PREC), PARAMETER :: &
     &            L321 =  4.012395124208693d-01,&
     &            L331 =  8.819910099032529d-04,&
     &            L341 =  1.728116022258560d-04,&
     &            L332 =  3.680857287181573d-01,&
     &            L342 =  1.635381132422046d-03,&
     &            L343 =  3.688541178419062d-01
!!--------- ORDER 5
!!--------- B3-B5 OF DIMENSION  4
       REAL(PREC), PARAMETER :: &
     &            B3511 = 49d0/720d0, &
     &            B3512 = -83d0/360d0,&
     &            B3513 = 17d0/60d0,  &
     &            B3514 = -53d0/360d0,&
     &            B3515 = 19d0/720d0, &
     &            B3521 = 19d0/720d0, &
     &            B3522 = -23d0/360d0,&
     &            B3523 = 1d0/30d0,   &
     & B3524 = 7d0/360d0,             &
     & B3525 = -11d0/720d0,           &
     & B3531 = -11d0/720d0,           &
     & B3532 = 37d0/360d0,            &
     & B3533 = -13d0/60d0,            &
     & B3534 = 67d0/360d0,            &
     & B3535 = -41d0/720d0,           &
     & B3541 = 19d0/720d0,            &
     & B3542 = -53d0/360d0,           &
     & B3543 =  17d0/60d0,            &
     & B3544 = -83d0/360d0,           &
     & B3545 =  49d0/720d0

!![ 49/720, -83/360,  17/60, -53/360,  19/720]
!![ 19/720, -23/360,   1/30,   7/360, -11/720]
!![-11/720,  37/360, -13/60,  67/360, -41/720]
!![ 19/720, -53/360,  17/60, -83/360,  49/720]

!!--------- A5, B5, B53, B56 Cp5 := MATRICES DEFINING  GAM5

       REAL(PREC), PARAMETER :: &
     & B511 = 251d0/720d0,      &
     & B512 = 323d0/360d0,      &
     & B513 = - 11d0/30d0,      &
     & B514 = 53d0/360d0,       &
     & B515 = -19d0/720d0,      &
     & B521 = -19d0/720d0,      &
     & B522 =  173d0/360d0,     &
     & B523 = 19d0/30d0,        &
     & B524 = -37d0/360d0,      &
     & B525 = 11d0/720d0

       REAL(PREC), PARAMETER:: &
     & B5711 = 1997d0/60480d0, &
     & B5712 = -113d0/630d0,   &
     & B5713 = 1619d0/4032d0,  &
     & B5714 = -715d0/1512d0,  &
     & B5715 = 1241d0/4032d0,  &
     & B5716 = -263d0/2520d0,  &
     & B5717 = 863d0/60480d0,  &
     & B5721 = -733d0/60480d0, &
     & B5722 = 41d0/630d0,     &
     & B5723 = -193d0/1344d0,  &
     & B5724 = 251d0/1512d0,   &
     & B5725 = -425d0/4032d0,  &
     & B5726 = 29d0/840d0,     &
     & B5727 = -271d0/60480d0, &
     & B5731 = -271d0/60480d0, &
     & B5732 = 97d0/5040d0,    &
     & B5733 = -13d0/448d0,    &
     & B5734 = 5d0/378d0,      &
     & B5735 = 37d0/4032d0,    &
     & B5736 = -19d0/1680d0,   &
     & B5737 = 191d0/60480d0,  &
     & B5741 = 191d0/60480d0,  &
     & B5742 = -67d0/2520d0,   &
     & B5743 = 115d0/1344d0,   &
     & B5744 = -211d0/1512d0,  &
     & B5745 = 499d0/4032d0,   &
     & B5746 = -2d0/35d0,      &
     & B5747 = 653d0/60480d0

!!--------- B5 - B7
!!
!![1997/60480,  -113/630, 1619/4032, -715/1512, 1241/4032, -263/2520,  863/60480]
!![-733/60480,    41/630, -193/1344,  251/1512, -425/4032,    29/840, -271/60480]
!![-271/60480,   97/5040,   -13/448,     5/378,   37/4032,  -19/1680,  191/60480]
!![ 191/60480,  -67/2520,  115/1344, -211/1512,  499/4032,     -2/35,  653/60480]
!![-271/60480,    29/840, -425/4032,  251/1512, -193/1344,    41/630, -733/60480]
!![ 863/60480, -263/2520, 1241/4032, -715/1512, 1619/4032,  -113/630, 1997/60480]
!!

       REAL(PREC), PARAMETER :: &
     &  L521 =  3.668340831928216D-01, &
     & L531 =  2.477905683677308D-03,  &
     & L541 = -1.919925047010838D-03,  &
     & L551 =  2.218385581234200D-03,  &
     & L561 = -5.442189351609260D-03,  &
     & L532 =  3.216639533696728D-01,  &
     & L542 =  1.231925763308414D-03,  &
     & L552 =  7.841944627374794D-03,  &
     & L562 =  1.002485104590053D-03,  &
     & L543 =  3.375100828961925D-01,  &
     & L553 = -2.614300734741796D-04,  &
     & L563 =  1.066631182323580D-03,  &
     & L554 =  3.523137378783708D-01,  &
     & L564 = -3.596681121610224D-04,  &
     & L565 =  3.617716171655064D-01,  &
     & CP51 =   3D0/160D0,             &
     & CP52 = -11D0/1440D0

!!--------- A7, B7, B57, B58 Cp7 := MATRICES DEFINING GAM7



        REAL(PREC), PARAMETER :: &
     & B711 = 19087d0/60480d0,   &
     & B712 = 2713d0/2520d0,     &
     & B713 = -15487d0/20160d0,  &
     & B714 = 586d0/945d0,       &
     & B715 = -6737d0/20160d0,   &
     & B716 = 263d0/2520d0,      &
     & B717 = -863d0/60480d0,    &
     & B721 = -863d0/60480d0,    &
     & B722 = 349d0/840d0,       &
     & B723 = 5221d0/6720d0,     &
     & B724 = -254d0/945d0,      &
     & B725 = 811d0/6720d0,      &
     & B726 = -29d0/840d0,       &
     & B727 = 271d0/60480d0,     &
     & B731 = 271d0/60480d0,     &
     & B732 = -23d0/504d0,       &
     & B733 = 10273d0/20160d0,   &
     & B734 = 586d0/945d0,       &
     & B735 = -2257d0/20160d0,   &
     & B736 = 67d0/2520d0,       &
     & B737 = -191d0/60480d0

!!
!!--------- THE LAST THREE ROWS ARE THE REVERSE OF THE FIRST THREE
!!
!! B79 =
!!
!![ 75203/3628800, -280187/1814400, 129781/259200, -238937/259200, 27289/25920, -197687/259200,  88531/259200, -156437/1814400,  33953/3628800]
!![-17827/3628800,   66043/1814400, -30389/259200,   55513/259200, -6281/25920,   44983/259200, -19859/259200,   34453/1814400,  -7297/3628800]
!![  8963/3628800,  -32987/1814400,  15061/259200,  -27257/259200,  3049/25920,  -21527/259200,   9331/259200,  -15797/1814400,   3233/3628800]
!![  3233/3628800,  -10067/1814400,   3601/259200,   -4337/259200,     23/3240,    1393/259200,  -2129/259200,    7123/1814400,  -2497/3628800]
!![ -2497/3628800,   12853/1814400,  -7859/259200,   18583/259200, -2681/25920,   24313/259200, -13589/259200,   30043/1814400,  -8227/3628800]
!![  3233/3628800,  -15797/1814400,   9331/259200,  -21527/259200,  3049/25920,  -27257/259200,  15061/259200,  -32987/1814400,   8963/3628800]
!![ -7297/3628800,   34453/1814400, -19859/259200,   44983/259200, -6281/25920,   55513/259200, -30389/259200,   66043/1814400, -17827/3628800]
!![ 33953/3628800, -156437/1814400,  88531/259200, -197687/259200, 27289/25920, -238937/259200, 129781/259200, -280187/1814400,  75203/3628800]
!!
       REAL(PREC), PARAMETER ::     &
     & B7911 =   75203d0/3628800d0, &
     & B7912 = -280187d0/1814400d0, &
     & B7913 =  129781d0/259200d0,  &
     & B7914 = -238937d0/259200d0,  &
     & B7915 =   27289d0/25920d0,   &
     & B7916 = -197687d0/259200d0,  &
     & B7917 =   88531d0/259200d0,  &
     & B7918 = -156437d0/1814400d0, &
     & B7919 =   33953d0/3628800d0

       REAL(PREC), PARAMETER ::    &
     & B7921 = -17827d0/3628800d0, &
     & B7922 =  66043d0/1814400d0, &
     & B7923 = -30389d0/259200d0,  &
     & B7924 =  55513d0/259200d0,  &
     & B7925 =  -6281d0/25920d0,   &
     & B7926 =  44983d0/259200d0,  &
     & B7927 = -19859d0/259200d0,  &
     & B7928 =  34453d0/1814400d0, &
     & B7929 =  -7297d0/3628800d0

       REAL(PREC), PARAMETER :: &
     & B7931 =   8963d0/3628800d0, &
     & B7932 = -32987d0/1814400d0, &
     & B7933 =  15061d0/259200d0,  &
     & B7934 = -27257d0/259200d0,  &
     & B7935 =   3049d0/25920d0,   &
     & B7936 = -21527d0/259200d0,  &
     & B7937 =   9331d0/259200d0,  &
     & B7938 = -15797d0/1814400d0, &
     & B7939 =   3233d0/3628800d0

       REAL(PREC), PARAMETER :: &
     & B7941 =   3233d0/3628800d0, &
     & B7942 = -10067d0/1814400d0, &
     & B7943 =   3601d0/259200d0,  &
     & B7944 =  -4337d0/259200d0,  &
     & B7945 =     23d0/3240d0,    &
     & B7946 =   1393d0/259200d0,  &
     & B7947 =  -2129d0/259200d0,  &
     & B7948 =   7123d0/1814400d0, &
     & B7949 =  -2497d0/3628800d0


      REAL(PREC),  PARAMETER :: &
     & B7951 =  -2497d0/3628800d0,&
     & B7952 =  12853d0/1814400d0,&
     & B7953 =  -7859d0/259200d0, &
     & B7954 =  18583d0/259200d0, &
     & B7955 =  -2681d0/25920d0,  &
     & B7956 =  24313d0/259200d0, &
     & B7957 = -13589d0/259200d0, &
     & B7958 =  30043d0/1814400d0,&
     & B7959 =  -8227d0/3628800d0

      REAL(PREC),  PARAMETER :: &
     & L721 = 3.023839891568610D-01, &
     & L731 = 3.201698610574002D-05, &
     & L741 = 4.193101163680004D-04, &
     & L751 = 1.686924996069667D-04, &
     & L761 = 4.806043527549464D-05, &
     & L771 = 3.598347048026785D-06, &
     & L781 = 7.892534649789167D-04, &
     & L732 = 2.559868364091398D-01, &
     & L742 = 1.336896192287030D-04, &
     & L752 = 3.080994719931695D-03, &
     & L762 = 1.457177183563680D-04, &
     & L772 = 9.259360509484074D-04, &
     & L782 = 2.397658879381223D-04, &
     & L743 = 2.639734712170458D-01, &
     & L753 = 1.734338929611258D-04, &
     & L763 = 6.704398263264620D-03, &
     & L773 = 4.559927214651730D-05, &
     & L783 = 6.396418554053151D-05, &
     & L754 = 2.817729090368562D-01, &
     & L764 = 2.877761776030408D-04, &
     & L774 = 1.810919475521773D-04, &
     & L784 = 1.009049833235848D-03, &
     & L765 = 2.993040718034231D-01, &
     & L775 = 2.009850887505898D-03, &
     & L785 = 1.748065618845750D-03, &
     & L776 = 3.150349043479135D-01, &
     & L786 = 3.243816792609449D-05, &
     & L787 = 3.271307059448932D-01

      REAL(PREC),  PARAMETER :: &
     & CP71 = 103D0/9061D0,     &
     & CP72 = -13D0/4480D0,     &
     & CP73 =  67D0/42431D0

!!--------- A8, B8, B86, B810 Cp8 := MATRICES DEFINING GAM9

       REAL(PREC), PARAMETER ::   &
     & B911 = 1070017D0/3628800D0,&
     & B912 = 2233547D0/1814400D0,&
     & B913 = -2302297D0/1814400D0,&
     & B914 = 2797679D0/1814400D0, &
     & B915 = -31457D0/22680D0,    &
     & B916 = 1573169D0/1814400D0, &
     & B917 = -645607D0/1814400D0, &
     & B918 = 156437D0/1814400D0,  &
     & B919 = -33953D0/3628800D0,  &
     & B921 = -33953D0/3628800D0,  &
     & B922 = 687797D0/1814400D0,  &
     & B923 =  1622393D0/1814400D0,&
     & B924 = -876271D0/1814400D0, &
     & B925 =   8233D0/22680D0,    &
     & B926 =    -377521D0/1814400D0, &
     & B927 =   147143D0/1814400D0,   &
     & B928 =  -34453D0/1814400D0,    &
     & B929 =   7297D0/3628800D0,     &
     & B931 = 7297D0/3628800D0,       &
     & B932 =  -49813D0/1814400D0,    &
     & B933 =  819143D0/1814400D0,    &
     & B934 =  1315919D0/1814400D0,   &
     & B935 = -5207D0/22680D0,        &
     & B936 =  198929D0/1814400D0,    &
     & B937 =  -71047D0/1814400D0,    &
     & B938 =  15797D0/1814400D0,     &
     & B939 = -3233D0/3628800D0

       REAL(PREC), PARAMETER ::  &
     & B941 = -3233D0/3628800D0,      &
     & B942 = 18197D0/1814400D0,      &
     & B943 =  -108007D0/1814400D0,   &
     & B944 =  954929D0/1814400D0,    &
     & B945 = 13903D0/22680D0,        &
     & B946 = -212881D0/1814400D0,    &
     & B947 =  63143D0/1814400D0,     &
     & B948 =  -12853D0/1814400D0,    &
     & B949 =  2497D0/3628800D0


!!   B910=B9-B10;
!!               prime 4 righe : la 5 e' uguale alla 4;
!!     le ultime 4 uguali alle prime 4 negate (simmetriche);
!!    le prime  5 colonne  : le altre sono simmetriche e negate;
!!
!![ 8183/1036800, -8183/115200,   8183/28800, -57281/86400,  57281/57600, -57281/57600,  57281/86400,  -8183/28800,  8183/115200, -8183/1036800]
!![  -425/290304,    425/32256,    -425/8064,     425/3456,    -425/2304,     425/2304,    -425/3456,     425/8064,   -425/32256,    425/290304]
!![      7/12800,    -63/12800,      63/3200,    -147/3200,     441/6400,    -441/6400,     147/3200,     -63/3200,     63/12800,      -7/12800]
!![-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600]
!![-2497/7257600,  2497/806400, -2497/201600,   2497/86400,  -2497/57600,   2497/57600,  -2497/86400,  2497/201600, -2497/806400,  2497/7257600]
!![ 2497/7257600, -2497/806400,  2497/201600,  -2497/86400,   2497/57600,  -2497/57600,   2497/86400, -2497/201600,  2497/806400, -2497/7257600]
!![     -7/12800,     63/12800,     -63/3200,     147/3200,    -441/6400,     441/6400,    -147/3200,      63/3200,    -63/12800,       7/12800]
!![   425/290304,   -425/32256,     425/8064,    -425/3456,     425/2304,    -425/2304,     425/3456,    -425/8064,    425/32256,   -425/290304]
!![-8183/1036800,  8183/115200,  -8183/28800,  57281/86400, -57281/57600,  57281/57600, -57281/86400,   8183/28800, -8183/115200,  8183/1036800]

       REAL(PREC), PARAMETER :: &
     & B91011 =   8183d0/1036800d0,  &
     & B91012 =  -8183d0/115200d0,  &
     & B91013 =   8183d0/28800d0,   &
     & B91014 = -57281d0/86400d0,   &
     & B91015 =  57281d0/57600d0

       REAL(PREC), PARAMETER ::   &
     & B91021 =   -425d0/290304d0,&
     & B91022 =    425d0/32256d0, &
     & B91023 =   -425d0/8064d0,  &
     & B91024 =    425d0/3456d0,  &
     & B91025 =   -425d0/2304d0

       REAL(PREC), PARAMETER :: &
     & B91031 =     7d0/12800d0,&
     & B91032 =   -63d0/12800d0,&
     & B91033 =    63d0/3200d0, &
     & B91034 =  -147d0/3200d0, &
     & B91035 =   441d0/6400d0

       REAL(PREC), PARAMETER ::    &
     & B91041 = -2497d0/7257600d0, &
     & B91042 =  2497d0/806400d0,  &
     & B91043 = -2497d0/201600d0,  &
     & B91044 =  2497d0/86400d0,   &
     & B91045 = -2497d0/57600d0

      REAL(PREC),  PARAMETER :: &
     & Cp91 =  7.892554012345216d-03, &
     & Cp92 = -1.463982583774219d-03, &
     & Cp93 =  5.468749999999983d-04, &
     & Cp94 = -3.440531305114634d-04

!!--------- THE OTHERS ARE THE SAME WITH CHANGED SIGN


       REAL(PREC), PARAMETER :: &
     & L921   =   2.590721934790442d-01,     &
     & L932   =   2.077575545359853d-01,     &
     & L943   =   2.032874698558627d-01,     &
     & L954   =   2.036384888660128d-01,     &
     & L965   =   2.039599505779785d-01,     &
     & L976   =   2.034044409161703d-01,     &
     & L987   =   2.017245408702437d-01,     &
     & L998   =   1.986549276295617d-01


CONTAINS

  SUBROUTINE DECLU(R,JF0,H,LDJAC,LU,LDLU,IPIV,FMAS,LDMAS,MLMAS,MUMAS,ORD,IER,IJOB)
    USE PRECISION
    USE  LINALGGAMD
    IMPLICIT NONE
    !!
    !!   INPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(IN):: R, LDJAC, LDLU, LDMAS, MLMAS, MUMAS, ORD, IJOB
    REAL(PREC),INTENT(IN) ::  JF0(LDJAC,R), FMAS(LDMAS,R),  H
    !!
    !!   OUTPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(OUT):: IER, IPIV(R)
    REAL(PREC),INTENT(OUT) :: LU(LDLU,R)
    !!
    !!   LOCAL VARIABLES
    !!------------------------------------
    INTEGER :: I, J, IB
    REAL(PREC) :: FAC


    REAL(PREC), PARAMETER :: L31  =  6.411501944628007d-01,  &
                           & L51  =  6.743555662880509D-01,  &
                           & L71  =  7.109158294404152D-01,  &
                           & L91  =  7.440547954061898d-01 

   
    !!
    !!   EXECUTABLE STATEMENTS
    !!---------------------------------
    !!
    SELECT CASE(ORD)

    CASE(1)
       FAC = -L31*H
    CASE(2)
       FAC = -L51*H
    CASE(3)    
       FAC = -L71*H
    CASE(4)
       FAC = -L91*H
    END SELECT

    SELECT CASE(IJOB)
    CASE(1)
       !! -------- ODE: JACOBIAN A FULL MATRIX
       DO J=1,R
          LU(1:R,J)= FAC*JF0(1:R,J)
          LU(J,J)=LU(J,J)+1d0
       END DO
       CALL DEC (R,LDLU,LU,IPIV,IER)
       RETURN
    CASE(2)
       !! -------- ODE: JACOBIAN A BAND MATRIX
       DO J=1,R
          LU(MLLU+1:MLLU+MDIAG,J)= FAC*JF0(1:MDIAG,J)
          LU(MDIAG,J)=LU(MDIAG,J)+1d0
       END DO
       CALL DECB (R,LDLU,LU,MLLU,MULU,IPIV,IER)
    CASE(3)
       !! -------- DAE: FMAS A BAND MATRIX, JACOBIAN A FULL MATRIX
       DO J=1,R
          LU(1:R,J)= FAC*JF0(1:R,J)
          DO I=MAX(1,J-MUMAS),MIN(R,J+MLMAS)
             LU(I,J)=LU(I,J)+FMAS(I-J+MBDIAG,J)
          END DO
       END DO
       CALL DEC (R,LDLU,LU,IPIV,IER)
    CASE(4)
       !! -------- DAE: FMAS A BAND MATRIX, JACOBIAN A BAND MATRIX
       DO J=1,R
          LU(MLLU+1:MLLU+MDIAG,J)= FAC*JF0(1:MDIAG,J)
          DO I=1,MBB
            IB=I+MDIFF
            LU(IB,J)=LU(IB,J)+FMAS(I,J)
          END DO
       END DO
       CALL DECB (R,LDLU,LU,MLLU,MULU,IPIV,IER)
    CASE(5)
       !! -------- DAE: FMAS A FULL MATRIX, JACOBIAN A FULL MATRIX
          LU(1:R,1:R)= FMAS(1:R,1:R) + FAC*JF0(1:R,1:R)
       CALL DEC (R,LDLU,LU,IPIV,IER)
    END SELECT
    RETURN

  END SUBROUTINE DECLU

!!
!!  SUBROUTINE SOLLU
!!

SUBROUTINE SOLLU(R,LU,LDLU,F,IPIV,IJOB)
  USE PRECISION
  USE LINALGGAMD
  IMPLICIT NONE

  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN):: R, LDLU, IPIV(R), IJOB
  REAL(PREC), INTENT(IN) ::  LU(LDLU,R)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: F(R)
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  SELECT CASE(IJOB)
  CASE(1,3,5)
     !! -------- JACOBIAN A FULL MATRIX
     CALL SOL (R,LDLU,LU(1,1),F(1),IPIV(1))
  CASE(2,4)
     !! -------- JACOBIAN A BAND MATRIX
     CALL SOLB (R,LDLU,LU(1,1),MLLU,MULU,F(1),IPIV(1))

  END SELECT

  RETURN

END SUBROUTINE SOLLU


!!
!!  SUBROUTINE NEWTGS
!!

SUBROUTINE NEWTGS(R,DBLK,LU,LDLU,FMAS,LDMAS,MLMAS,MUMAS,IPIV,F,DN,IJOB)
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN):: R, LDLU, LDMAS, MLMAS, MUMAS, IPIV(R), IJOB, DBLK
  REAL(PREC), INTENT(IN)::  F(R,DBLK), LU(LDLU,R), FMAS(LDMAS,R)
  !!
  !!   OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(OUT)::  DN(R,DBLK)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER     ::  J
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!

  SELECT CASE(IJOB)
  CASE(1,2)  !! ODE

     DN(1:R,1) = -F(1:R,1)
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1,1),IPIV(1),IJOB)
     DO J=2,DBLK
        DN(1:R,J) =  -F(1:R,J)+DN(1:R,J-1)
        CALL  SOLLU(R,LU(1,1),LDLU,DN(1,J),IPIV(1),IJOB)
     END DO

  CASE(3,4)  !! DAE: FMAS A BAND MATRIX

     DN(1:R,1) = -F(1:R,1)
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1,1),IPIV(1),IJOB)
     DO J=2,DBLK
        DN(1:R,J) =  -F(1:R,J) + MATMULB(R,FMAS(1,1),LDMAS,MLMAS,MUMAS,DN(1:R,J-1))
        CALL  SOLLU(R,LU(1,1),LDLU,DN(1,J),IPIV(1),IJOB)
     END DO

  CASE(5)    !! DAE: FMAS A FULL MATRIX

     DN(1:R,1) = -F(1:R,1)
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1,1),IPIV(1),IJOB)
     DO J=2,DBLK
        DN(1:R,J) =  -F(1:R,J) + MATMUL(FMAS(1:R,1:R),DN(1:R,J-1))
        CALL  SOLLU(R,LU(1,1),LDLU,DN(1,J),IPIV(1),IJOB)
     END DO
 


  END SELECT 

  RETURN

END SUBROUTINE NEWTGS

!!
!!   FINE SUBROUTINE NEWTGS
!!
SUBROUTINE INTERP(R,TP,YP,T1,F1,NT1,DBLKOLD,DBLK,T0,Y0)
  USE PRECISION
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: R, DBLK, DBLKOLD, NT1
  REAL(PREC), INTENT(IN) :: F1(R,*), T1(*), Y0(R), T0
  REAL(PREC), INTENT(IN OUT) :: YP(R,*), TP(*)
  INTEGER :: I,J, N, IT1, NT2, NTi
  
  IF (DBLK < DBLKOLD) THEN
     NTi = MAX(5,NT1)
  ELSE
     NTi = MAX(3,NT1)
  END IF
  NT2 = NTi+1
  N = DBLKOLD+1


  DO IT1=2,DBLK+1
     YP(1:R,IT1) = F1(1:R,NTi)
     DO J=NT2,N
        YP(1:R,IT1) = YP(1:R,IT1)*(T1(IT1-1)-TP(J)) + F1(1:R,J)
     END DO
  END DO

 

  YP(1:R,1) = Y0(1:R)

  TP(1) = T0
  TP(2:DBLK+1) = T1(1:DBLK)

  RETURN
END SUBROUTINE INTERP

!!
!!    SUBROUTINE DIFFDIV
!!

SUBROUTINE DIFFDIV(TP,YP,R,DBLK,NT1)
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!-----------------------------------
  INTEGER, INTENT(IN) :: R,  DBLK
  !!
  !!   INPUT/OUTPUT VARIABLES
  !------------------------------------
  INTEGER, INTENT(IN OUT):: NT1
  REAL(PREC), INTENT(IN OUT):: TP(DBLK+1), YP(R,DBLK+1)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  ::  J, K, N
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!

  N = DBLK+1
  NT1 = 1

!write(*,*) 'TP=',TP
!write(*,*) 'YP=',YP(1,:)

  DO J=N-1,NT1,-1
     DO K=1,J
        YP(1:R,K)= ( YP(1:R,K)- YP(1:R,K+1) )/( TP(K)-TP(K+N-J))
     END DO
  END DO
!write(*,*) 'YP=',YP(1,:)
!read(*,*)

  RETURN
END SUBROUTINE DIFFDIV
!!
!!     FUNCTION  CONTR
!!
FUNCTION CONTR(I,R,T,TP,FF,DBLK,NT1)
  !! ----------------------------------------------------------
  !!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT. IT PROVIDES AN
  !!     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT T.
  !!     IT GIVES THE VALUE OF THE INTERPOLATION POLYNOMIAL, DEFINED FOR
  !!     THE LAST SUCCESSFULLY COMPUTED STEP.
  !! ----------------------------------------------------------
  USE PRECISION
  IMPLICIT NONE
  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) ::  R, I, DBLK, NT1
  REAL(PREC) :: CONTR
  REAL(PREC), INTENT(IN) ::  T, TP(DBLK+1), FF(R,DBLK+1)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER :: J, N, NT2, NTi
  REAL(PREC) :: YP
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  N = DBLK+1
  NTi = MAX(1,NT1)
  NT2=NTi+1
  YP = FF(I,NTi)
  DO J=NT2,N
     YP = YP*(T-TP(J)) + FF(I,J)
  END DO
  CONTR = YP
  RETURN
END FUNCTION  CONTR


    !!
    !!  SUBROUTINE TERMNOT3  (ORDER 3)
    !!
SUBROUTINE  TERMNOT3(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    &
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            &
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR) 

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1)
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0,TETAK0,RPAR(1)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,5), FP(R,5), &
       &                   DN(R), F1(R,5),TP(5)
  INTEGER, INTENT(IN OUT) :: NFCN 
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)
 
  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE
  TER = .FALSE.
 
 
  SELECT CASE(IJOB) 
  CASE(1,2)  !! ODE 
 
      DN(1:R)=(YP(1:R,2)-YP(1:R,1))- & 
       &   H*(B311*FP(1:R,1)+B312*FP(1:R,2)+B313*FP(1:R,3)) 
 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,2)=YP(J,2)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
 
      DO J=1,R 
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2) & 
              &               +B312*FP(J,3)+ B313*FP(J,4) 
         DN(J)=(YP(J,3)-YP(J,2))-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,3)=YP(J,3)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(3),YP(1,3),F1(1,2),IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      DO J=1,R 
         SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3)) & 
              & +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5) 
         DN(J)=(YP(J,4)-YP(J,3))-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,4)=YP(J,4)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
      NFCN = NFCN + 1 
 
 
      DO J=1,R 
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3)) & 
              & +L343*(F1(J,3)-FP(J,4))  & 
              & +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5) 
         DN(J)=(YP(J,5)-YP(J,4))-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,5)=YP(J,5)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
  
 
      ERRVJ = sqrt(ERRVJ/R) 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
    
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
          TER = .TRUE. 
         RETURN 
      END IF 
 
 
  CASE(3,4)  !! DAE: FMAS A BAND MATRIX 
 
      DN(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1))- & 
       &   H*(B311*FP(1:R,1)+B312*FP(1:R,2)+B313*FP(1:R,3)) 
 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,2)=YP(J,2)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
      ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2)) 
      DO J=1,R 
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2) & 
              &               +B312*FP(J,3)+ B313*FP(J,4) 
         DN(J)=ZP(J)-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,3)=YP(J,3)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(3),YP(1,3),F1(1,2),IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3)) 
      DO J=1,R 
         SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3)) & 
              & +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5) 
         DN(J)=ZP(J)-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,4)=YP(J,4)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
	  NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
       
 
 
      ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4)) 
      DO J=1,R 
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3)) & 
              & +L343*(F1(J,3)-FP(J,4))  & 
              & +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5) 
         DN(J)=ZP(J)-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,5)=YP(J,5)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
  
 
      ERRVJ = sqrt(ERRVJ/R) 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
    
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
          TER = .TRUE. 
         RETURN 
      END IF 
 
 
  CASE(5)    !! DAE: FMAS A FULL MATRIX 
 
      DN(1:R)=MATMUL(FMAS(1:R,1:R),(YP(1:R,2)-YP(1:R,1)))- & 
       &   H*(B311*FP(1:R,1)+B312*FP(1:R,2)+B313*FP(1:R,3)) 
 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,2)=YP(J,2)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2)) 
      DO J=1,R 
         SUM = L321*(F1(J,1)-FP(J,2))+B311*FP(J,2) & 
              &               +B312*FP(J,3)+ B313*FP(J,4) 
         DN(J)=ZP(J)-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,3)=YP(J,3)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(3),YP(1,3),F1(1,2),IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3)) 
      DO J=1,R 
         SUM = L331*(F1(J,1)-FP(J,2))+L332*(F1(J,2)-FP(J,3)) & 
              & +B311*FP(J,3)+B312*FP(J,4)+B313*FP(J,5) 
         DN(J)=ZP(J)-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,4)=YP(J,4)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
      ERRVJ = sqrt(ERRVJ/R) 
 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
	  NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
       
 
 
      ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4)) 
      DO J=1,R 
         SUM = L341*(F1(J,1)-FP(J,2))+L342*(F1(J,2)-FP(J,3)) & 
              & +L343*(F1(J,3)-FP(J,4))  & 
              & +B313*FP(J,3)+B312*FP(J,4)+B311*FP(J,5) 
         DN(J)=ZP(J)-H*SUM 
      END DO 
      CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
      ERRVJ = 0D0 
      DO J=1,R 
         YP(J,5)=YP(J,5)-DN(J) 
         SUM = (DN(J)/SCAL(J)) 
         ERRVJ =  ERRVJ + SUM*SUM 
      END DO 
  
 
      ERRVJ = sqrt(ERRVJ/R) 
      ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
    
 
      IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
         TER = .TRUE. 
         RETURN 
      END IF 
 
      IERR = 0 
      CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
      NFCN = NFCN + 1 
      IF (IERR .NE.0) THEN 
          TER = .TRUE. 
         RETURN 
      END IF 
 
  END SELECT 

  FP(1:R,2:5) = F1(1:R,1:4)

  RETURN
END SUBROUTINE TERMNOT3
    !!
    !!  SUBROUTINE TERMNOT5  (ORDER 5)
    !!        
SUBROUTINE  TERMNOT5(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    & 
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            & 
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR) 

  USE PRECISION
  IMPLICIT NONE

  !!   INPUT VARIABLES
  !!------------------------------------ 
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1) 
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0, TETAK0, RPAR(1) 
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,7), FP(R,7), &
       &                   DN(R), F1(R,7),TP(7)
  INTEGER, INTENT(IN OUT) :: NFCN 
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)

  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE

 

  TER = .FALSE. 
  SELECT CASE(IJOB) 
  CASE(1,2)  !! ODE 

     DO J=1,R
        SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)&
             &    +B515*FP(J,5)
        DN(J)=YP(J,2)-YP(J,1)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF
 
     DO J=1,R
        SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)&
             &               +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5)
        DN(J)=YP(J,3)-YP(J,2)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 

  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


     DO J=1,R
        SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))&
             &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)+B524*FP(J,5)+B525*FP(J,6)
        DN(J)=YP(J,4)-YP(J,3)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))&
             &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)&
             &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7)
        DN(J)=YP(J,5)-YP(J,4)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF
  
     DO J=1,R
        SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))&
             &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))&
             &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7)
        DN(J)=YP(J,6)-YP(J,5)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )

  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))&
             &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))&
             &+L565*(F1(J,5)-FP(J,6))&
             &+B515*FP(J,3)+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7)
        DN(J)=YP(J,7)-YP(J,6)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6),IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF 

  CASE(3,4)  !! DAE: FMAS A BAND MATRIX 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1)) 
     DO J=1,R 
        SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)& 
             &    +B515*FP(J,5) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,2)=YP(J,2)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR = 0 
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
         TER = .TRUE. 
        RETURN 
     END IF 
  
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2)) 
     DO J=1,R 
        SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)& 
             &               +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,3)=YP(J,3)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
 
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3)) 
     DO J=1,R 
        SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))& 
             &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)+B524*FP(J,5)+B525*FP(J,6) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,4)=YP(J,4)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4)) 
     DO J=1,R 
        SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))& 
             &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)& 
             &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,5)=YP(J,5)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,6)-YP(1:R,5)) 
     DO J=1,R 
        SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))& 
             &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))& 
             &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,6)=YP(J,6)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,7)-YP(1:R,6)) 
     DO J=1,R 
        SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))& 
             &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))& 
             &+L565*(F1(J,5)-FP(J,6))& 
             &+B515*FP(J,3)+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,7)=YP(J,7)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR = 0 
     CALL FCN(R,TP(7),YP(1,7),F1(1,6),IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
 
 

  CASE(5)    !! DAE: FMAS A FULL MATRIX 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,2)-YP(1:R,1)) 
     DO J=1,R 
        SUM = B511*FP(J,1)+B512*FP(J,2)+B513*FP(J,3)+B514*FP(J,4)& 
             &    +B515*FP(J,5) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,2)=YP(J,2)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR = 0 
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
         TER = .TRUE. 
        RETURN 
     END IF 
  
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2)) 
     DO J=1,R 
        SUM = L521*(F1(J,1)-FP(J,2))+B521*FP(J,1)+B522*FP(J,2)& 
             &               +B523*FP(J,3)+B524*FP(J,4)+B525*FP(J,5) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,3)=YP(J,3)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
 
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3)) 
     DO J=1,R 
        SUM = L531*(F1(J,1)-FP(J,2))+L532*(F1(J,2)-FP(J,3))& 
             &+B521*FP(J,2)+B522*FP(J,3)+B523*FP(J,4)+B524*FP(J,5)+B525*FP(J,6) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,4)=YP(J,4)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4)) 
     DO J=1,R 
        SUM = L541*(F1(J,1)-FP(J,2))+L542*(F1(J,2)-FP(J,3))& 
             &+L543*(F1(J,3)-FP(J,4))+B521*FP(J,3)+B522*FP(J,4)& 
             &+B523*FP(J,5)+B524*FP(J,6)+B525*FP(J,7) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,5)=YP(J,5)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,6)-YP(1:R,5))   
     DO J=1,R 
        SUM =  L551*(F1(J,1)-FP(J,2))+L552*(F1(J,2)-FP(J,3))& 
             &+L553*(F1(J,3)-FP(J,4))+L554*(F1(J,4)-FP(J,5))& 
             &+B525*FP(J,3)+B524*FP(J,4)+B523*FP(J,5)+B522*FP(J,6)+B521*FP(J,7) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,6)=YP(J,6)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
 
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,7)-YP(1:R,6)) 
     DO J=1,R 
        SUM =  L561*(F1(J,1)-FP(J,2))+L562*(F1(J,2)-FP(J,3))& 
             &+L563*(F1(J,3)-FP(J,4))+L564*(F1(J,4)-FP(J,5))& 
             &+L565*(F1(J,5)-FP(J,6))& 
             &+B515*FP(J,3)+B514*FP(J,4)+B513*FP(J,5)+B512*FP(J,6)+B511*FP(J,7) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,7)=YP(J,7)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR = 0 
     CALL FCN(R,TP(7),YP(1,7),F1(1,6),IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
 
END SELECT

  FP(1:R,2:7) = F1(1:R,1:6)

  RETURN
END SUBROUTINE TERMNOT5

    !!
    !!  SUBROUTINE TERMNOT7  (ORDER 7)
    !!
SUBROUTINE  TERMNOT7(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    & 
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            & 
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR) 
 
  USE PRECISION 
  IMPLICIT NONE 
 
  !!   INPUT VARIABLES 
  !!------------------------------------ 
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1) 
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0, TETAK0, RPAR(1) 
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,9), FP(R,9), &
       &                   DN(R), F1(R,9),TP(9)
  INTEGER, INTENT(IN OUT) :: NFCN 
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)
 
  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE
 

  TER = .FALSE. 
 
  SELECT CASE(IJOB) 
  CASE(1,2)  !! ODE 

     DO J=1,R
        SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)&
             &          +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7)
        DN(J) =YP(J,2)-YP(J,1)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2) - DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
 
     IERR = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+&
             &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7)
        DN(J) = YP(J,3)-YP(J,2)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF
 
     IERR = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
         TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2) -FP(J,3))&
             & +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)&
             & +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7)
        DN(J) =YP(J,4)-YP(J,3)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF
    
     IERR  = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))&
             & +L743*(F1(J,3)-FP(J,4))&
             & +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)&
             & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8)
        DN(J) =YP(J,5)-YP(J,4)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
    
     IERR  = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))&
             &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))&
             & +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)&
             & +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9)
        DN(J)=YP(J,6)-YP(J,5)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
    
     IERR  = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))&
             &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))&
             &+L765*(F1(J,5)-FP(J,6))&
             & +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)&
             & +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9)
        DN(J) = YP(J,7)-YP(J,6)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF
   
     IERR  = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), IERR,  RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))&
        &+ L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))&
        &+ L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))&
        & +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)&
        & +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9)
        DN(J) = YP(J,8)-YP(J,7)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF


     DO J=1,R
        SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))&
        &+ L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))&
        &+ L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))&
        &+ L787*(F1(J,7)-FP(J,8))&
        & +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)&
        & +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9)
        DN(J) = YP(J,9)-YP(J,8)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
 
     IERR  = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF 
 
 
  CASE(3,4)  !! DAE: FMAS A BAND MATRIX 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1)) 
     DO J=1,R 
        SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)& 
             &          +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7) 
        DN(J) =ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,2)=YP(J,2) - DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR = 0 
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
         TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2)) 
     DO J=1,R 
        SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+& 
             &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,3)=YP(J,3)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR = 0 
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
         TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3)) 
     DO J=1,R 
        SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2) -FP(J,3))& 
             & +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)& 
             & +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7) 
        DN(J) =ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,4)=YP(J,4)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR  = 0 
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4)) 
     DO J=1,R 
        SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))& 
             & +L743*(F1(J,3)-FP(J,4))& 
             & +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)& 
             & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8) 
        DN(J) =ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,5)=YP(J,5)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR  = 0 
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,6)-YP(1:R,5)) 
     DO J=1,R 
        SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))& 
             &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))& 
             & +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)& 
             & +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,6)=YP(J,6)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR  = 0 
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,7)-YP(1:R,6)) 
     DO J=1,R 
        SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))& 
             &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))& 
             &+L765*(F1(J,5)-FP(J,6))& 
             & +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)& 
             & +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,7)=YP(J,7)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR  = 0 
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), IERR,  RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,8)-YP(1:R,7)) 
     DO J=1,R 
        SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))& 
        &+ L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))& 
        &+ L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))& 
        & +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)& 
        & +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,8)=YP(J,8)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,9)-YP(1:R,8)) 
     DO J=1,R 
        SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))& 
        &+ L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))& 
        &+ L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))& 
        &+ L787*(F1(J,7)-FP(J,8))& 
        & +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)& 
        & +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,9)=YP(J,9)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR  = 0 
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
  CASE(5)    !! DAE: FMAS A FULL MATRIX 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,2)-YP(1:R,1)) 
     DO J=1,R 
        SUM= B711*FP(J,1)+B712*FP(J,2)+B713*FP(J,3)& 
             &          +B714*FP(J,4)+B715*FP(J,5)+B716*FP(J,6)+B717*FP(J,7) 
        DN(J) =ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,2)=YP(J,2) - DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR = 0 
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
         TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2)) 
     DO J=1,R 
        SUM = L721*(F1(J,1)-FP(J,2))+B721*FP(J,1)+B722*FP(J,2)+& 
             &B723*FP(J,3)+B724*FP(J,4)+B725*FP(J,5)+B726*FP(J,6)+B727*FP(J,7) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,3)=YP(J,3)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR = 0 
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
         TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3)) 
     DO J=1,R 
        SUM = +L731*(F1(J,1)-FP(J,2))+L732*(F1(J,2) -FP(J,3))& 
             & +B731*FP(J,1)+B732*FP(J,2)+B733*FP(J,3)& 
             & +B734*FP(J,4)+B735*FP(J,5)+B736*FP(J,6)+B737*FP(J,7) 
        DN(J) =ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,4)=YP(J,4)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR  = 0 
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4)) 
     DO J=1,R 
        SUM =L741*(F1(J,1)-FP(J,2))+L742*(F1(J,2)-FP(J,3))& 
             & +L743*(F1(J,3)-FP(J,4))& 
             & +B731*FP(J,2)+B732*FP(J,3)+B733*FP(J,4)& 
             & +B734*FP(J,5)+B735*FP(J,6)+B736*FP(J,7)+B737*FP(J,8) 
        DN(J) =ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,5)=YP(J,5)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR  = 0 
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,6)-YP(1:R,5)) 
     DO J=1,R 
        SUM = L751*(F1(J,1)-FP(J,2))+L752*(F1(J,2)-FP(J,3))& 
             &+L753*(F1(J,3)-FP(J,4))+L754*(F1(J,4)-FP(J,5))& 
             & +B731*FP(J,3)+B732*FP(J,4)+B733*FP(J,5)& 
             & +B734*FP(J,6)+B735*FP(J,7)+B736*FP(J,8)+B737*FP(J,9) 
        DN(J)=ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,6)=YP(J,6)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR  = 0 
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,7)-YP(1:R,6)) 
     DO J=1,R 
        SUM = L761*(F1(J,1)-FP(J,2))+L762*(F1(J,2)-FP(J,3))& 
             &+L763*(F1(J,3)-FP(J,4))+L764*(F1(J,4)-FP(J,5))& 
             &+L765*(F1(J,5)-FP(J,6))& 
             & +B737*FP(J,3)+B736*FP(J,4)+B735*FP(J,5)& 
             & +B734*FP(J,6)+B733*FP(J,7)+B732*FP(J,8)+B731*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,7)=YP(J,7)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR  = 0 
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), IERR,  RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,8)-YP(1:R,7)) 
     DO J=1,R 
        SUM = L771*(F1(J,1)-FP(J,2))+L772*(F1(J,2)-FP(J,3))& 
        &+ L773*(F1(J,3)-FP(J,4))+L774*(F1(J,4)-FP(J,5))& 
        &+ L775*(F1(J,5)-FP(J,6))+L776*(F1(J,6)-FP(J,7))& 
        & +B727*FP(J,3)+B726*FP(J,4)+B725*FP(J,5)& 
        & +B724*FP(J,6)+B723*FP(J,7)+B722*FP(J,8)+B721*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,8)=YP(J,8)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,9)-YP(1:R,8)) 
     DO J=1,R 
        SUM = L781*(F1(J,1)-FP(J,2))+L782*(F1(J,2)-FP(J,3))& 
        &+ L783*(F1(J,3)-FP(J,4))+L784*(F1(J,4)-FP(J,5))& 
        &+ L785*(F1(J,5)-FP(J,6))+L786*(F1(J,6)-FP(J,7))& 
        &+ L787*(F1(J,7)-FP(J,8))& 
        & +B717*FP(J,3)+B716*FP(J,4)+B715*FP(J,5)& 
        & +B714*FP(J,6)+B713*FP(J,7)+B712*FP(J,8)+B711*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,9)=YP(J,9)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
  
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR  = 0 
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
 
 
 
  END SELECT

  FP(1:R,2:9) = F1(1:R,1:8)



  RETURN


END SUBROUTINE TERMNOT7
    !!
    !!  SUBROUTINE TERMNOT9  (ORDER 9)
    !!
SUBROUTINE  TERMNOT9(R,FCN,H,IT,DN, F1,FP,YP,TP,NFCN,                    & 
                    &  ERRNEWT,ERRNEWT0,TETAK0,LU, LDLU,IPIV,            & 
                    &  FMAS,LDMAS,MLMAS,MUMAS, SCAL,IJOB,TER,RPAR,IPAR) 
 
  USE PRECISION 
  IMPLICIT NONE 
 
  !!   INPUT VARIABLES 
  !!------------------------------------ 
  INTEGER, INTENT(IN) :: R,  IJOB, IPIV(R), LDLU, LDMAS, MLMAS, MUMAS, IT, IPAR(1) 
  REAL(PREC), INTENT(IN) ::  H, SCAL(R), LU(LDLU,R), FMAS(LDMAS,R), ERRNEWT0, TETAK0, RPAR(1) 
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) :: ERRNEWT, YP(R,10), FP(R,10), &
       &                   DN(R), F1(R,10),TP(10)
  INTEGER, INTENT(IN OUT) :: NFCN 
  LOGICAL, INTENT(OUT):: TER
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: J, IERR
  REAL(PREC) ::  ERRVJ, SUM, FAC, ZP(R)
 
  EXTERNAL FCN
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!
  !!--------- ONE STEP OF THE ITERATION PROCEDURE

 
  TER = .FALSE.
  SELECT CASE(IJOB) 
  CASE(1,2)  !! ODE 

     DO J=1,R
        SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)&
             &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9)
        DN(J) = YP(J,2)-YP(J,1)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,2)=YP(J,2)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN
        TER = .TRUE.
        RETURN
     END IF

     IERR  = 0
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)&
             &   +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)&
             &   +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9)
        DN(J) = YP(J,3)-YP(J,2)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,3)=YP(J,3)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
 
     IERR  = 0
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L932*(F1(J,2)-FP(J,3))&
             &+B931*FP(J,1)+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)&
             &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9)
             DN(J) = YP(J,4)-YP(J,3)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,4)=YP(J,4)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
   
     IERR = 0
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)&
             & +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)&
             & +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9)
        DN(J) =YP(J,5)-YP(J,4)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,5)=YP(J,5)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF
     
     IERR = 0
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF
   
     DO J=1,R
        SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)&
             & +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)&
             & +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10)
        DN(J) =YP(J,6)-YP(J,5)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,6)=YP(J,6)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
    
     IERR = 0
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF
   
     DO J=1,R
           SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)&
             & +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)&
             & +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10)
        DN(J) =YP(J,7)-YP(J,6)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,7)=YP(J,7)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
   
     IERR = 0
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)&
             & +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)&
             & +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10)
        DN(J) =YP(J,8)-YP(J,7)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,8)=YP(J,8)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN
        TER = .TRUE.
        RETURN
     END IF
   
     IERR = 0
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)&
             & +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)&
             & +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10)
        DN(J) =YP(J,9)-YP(J,8)-H*SUM
     END DO
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,9)=YP(J,9)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
     
     IERR = 0
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), IERR, RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF

     DO J=1,R
        SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)&
             & +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)&
             & +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10)
        DN(J) =YP(J,10)-YP(J,9)-H*SUM
     END DO
   
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB)
     ERRVJ = 0D0
     DO J=1,R
        YP(J,10)=YP(J,10)-DN(J)
        SUM = (DN(J)/SCAL(J))
        ERRVJ =  ERRVJ + SUM*SUM
     END DO
     ERRVJ = sqrt(ERRVJ/R)
     ERRNEWT = MAX( ERRNEWT, ERRVJ )
   
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN
        TER = .TRUE.
        RETURN
     END IF
     
     IERR = 0
     CALL FCN(R,TP(10),YP(1,10),F1(1,9),IERR,  RPAR,IPAR)
     NFCN = NFCN + 1
     IF (IERR .NE.0) THEN
        TER = .TRUE.
        RETURN
     END IF
  
  
  CASE(3,4)  !! DAE: FMAS A BAND MATRIX 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,2)-YP(1:R,1)) 
     DO J=1,R 
        SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)& 
             &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,2)=YP(J,2)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,3)-YP(1:R,2)) 
     DO J=1,R 
        SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)& 
             &   +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)& 
             &   +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,3)=YP(J,3)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR  = 0 
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,4)-YP(1:R,3)) 
     DO J=1,R 
        SUM = L932*(F1(J,2)-FP(J,3))& 
             &+B931*FP(J,1)+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)& 
             &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9) 
             DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,4)=YP(J,4)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR = 0 
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,5)-YP(1:R,4)) 
     DO J=1,R 
        SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)& 
             & +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)& 
             & +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,5)=YP(J,5)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
      
     IERR = 0 
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,6)-YP(1:R,5))   
     DO J=1,R 
        SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)& 
             & +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)& 
             & +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,6)=YP(J,6)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR = 0 
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,7)-YP(1:R,6))    
     DO J=1,R 
           SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)& 
             & +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)& 
             & +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,7)=YP(J,7)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR = 0 
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,8)-YP(1:R,7)) 
     DO J=1,R 
        SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)& 
             & +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)& 
             & +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,8)=YP(J,8)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR = 0 
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,9)-YP(1:R,8)) 
     DO J=1,R 
        SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)& 
             & +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)& 
             & +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,9)=YP(J,9)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
      
     IERR = 0 
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMULB(R,FMAS,LDMAS,MLMAS,MUMAS,YP(1:R,10)-YP(1:R,9)) 
     DO J=1,R 
        SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)& 
             & +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)& 
             & +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
    
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,10)=YP(J,10)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
      
     IERR = 0 
     CALL FCN(R,TP(10),YP(1,10),F1(1,9),IERR,  RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
   
  CASE(5)    !! DAE: FMAS A FULL MATRIX 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,2)-YP(1:R,1)) 
     DO J=1,R 
        SUM = B911*FP(J,1)+B912*FP(J,2)+B913*FP(J,3)+B914*FP(J,4)& 
             &+B915*FP(J,5)+B916*FP(J,6)+B917*FP(J,7)+B918*FP(J,8)+B919*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,2)=YP(J,2)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     IERR  = 0 
     CALL FCN(R,TP(2),YP(1,2),F1(1,1), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,3)-YP(1:R,2)) 
     DO J=1,R 
        SUM = L921*(F1(J,1)-FP(J,2))+B921*FP(J,1)+B922*FP(J,2)& 
             &   +B923*FP(J,3)+B924*FP(J,4)+B925*FP(J,5)& 
             &   +B926*FP(J,6)+B927*FP(J,7)+B928*FP(J,8)+B929*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,3)=YP(J,3)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     IERR  = 0 
     CALL FCN(R,TP(3),YP(1,3),F1(1,2), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,4)-YP(1:R,3)) 
     DO J=1,R 
        SUM = L932*(F1(J,2)-FP(J,3))& 
             &+B931*FP(J,1)+B932*FP(J,2)+B933*FP(J,3)+B934*FP(J,4)+B935*FP(J,5)& 
             &+B936*FP(J,6)+B937*FP(J,7)+B938*FP(J,8)+B939*FP(J,9) 
             DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,4)=YP(J,4)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR = 0 
     CALL FCN(R,TP(4),YP(1,4),F1(1,3), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,5)-YP(1:R,4)) 
     DO J=1,R 
        SUM = L943*(F1(J,3)-FP(J,4))+B941*FP(J,1)+B942*FP(J,2)& 
             & +B943*FP(J,3)+B944*FP(J,4)+B945*FP(J,5)+B946*FP(J,6)& 
             & +B947*FP(J,7)+B948*FP(J,8)+B949*FP(J,9) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,5)=YP(J,5)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
      
     IERR = 0 
     CALL FCN(R,TP(5),YP(1,5),F1(1,4), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
  
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,6)-YP(1:R,5))   
     DO J=1,R 
        SUM = L954*(F1(J,4)-FP(J,5))+B941*FP(J,2)+B942*FP(J,3)& 
             & +B943*FP(J,4)+B944*FP(J,5)+B945*FP(J,6)+B946*FP(J,7)& 
             & +B947*FP(J,8)+B948*FP(J,9)+B949*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,6)=YP(J,6)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
     
     IERR = 0 
     CALL FCN(R,TP(6),YP(1,6),F1(1,5), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,7)-YP(1:R,6))    
     DO J=1,R 
           SUM = L965*(F1(J,5)-FP(J,6))+B949*FP(J,2)+B948*FP(J,3)& 
             & +B947*FP(J,4)+B946*FP(J,5)+B945*FP(J,6)+B944*FP(J,7)& 
             & +B943*FP(J,8)+B942*FP(J,9)+B941*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,7)=YP(J,7)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR = 0 
     CALL FCN(R,TP(7),YP(1,7),F1(1,6), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,8)-YP(1:R,7)) 
     DO J=1,R 
        SUM = L976*(F1(J,6)-FP(J,7))+B939*FP(J,2)+B938*FP(J,3)& 
             & +B937*FP(J,4)+B936*FP(J,5)+B935*FP(J,6)+B934*FP(J,7)& 
             & +B933*FP(J,8)+B932*FP(J,9)+B931*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,8)=YP(J,8)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 ) ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
    
     IERR = 0 
     CALL FCN(R,TP(8),YP(1,8),F1(1,7), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,9)-YP(1:R,8)) 
     DO J=1,R 
        SUM = L987*(F1(J,7)-FP(J,8))+B929*FP(J,2)+B928*FP(J,3)& 
             & +B927*FP(J,4)+B926*FP(J,5)+B925*FP(J,6)+B924*FP(J,7)& 
             & +B923*FP(J,8)+B922*FP(J,9)+B921*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,9)=YP(J,9)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
     
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
      
     IERR = 0 
     CALL FCN(R,TP(9),YP(1,9),F1(1,8), IERR, RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
 
     ZP(1:R)=MATMUL(FMAS(1:R,1:R),YP(1:R,10)-YP(1:R,9)) 
     DO J=1,R 
        SUM = L998*(F1(J,8)-FP(J,9))+B919*FP(J,2)+B918*FP(J,3)& 
             & +B917*FP(J,4)+B916*FP(J,5)+B915*FP(J,6)+B914*FP(J,7)& 
             & +B913*FP(J,8)+B912*FP(J,9)+B911*FP(J,10) 
        DN(J) = ZP(J)-H*SUM 
     END DO 
    
     CALL  SOLLU(R,LU(1,1),LDLU,DN(1),IPIV(1),IJOB) 
     ERRVJ = 0D0 
     DO J=1,R 
        YP(J,10)=YP(J,10)-DN(J) 
        SUM = (DN(J)/SCAL(J)) 
        ERRVJ =  ERRVJ + SUM*SUM 
     END DO 
     ERRVJ = sqrt(ERRVJ/R) 
     ERRNEWT = MAX( ERRNEWT, ERRVJ ) 
    
     IF ((it.GE.1).AND.(errnewt/errnewt0 .GT. tetak0 )  ) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
      
     IERR = 0 
     CALL FCN(R,TP(10),YP(1,10),F1(1,9),IERR,  RPAR,IPAR) 
     NFCN = NFCN + 1 
     IF (IERR .NE.0) THEN 
        TER = .TRUE. 
        RETURN 
     END IF 
   
 
 
  END SELECT 
 

  FP(1:R,2:10) = F1(1:R,1:9)

  RETURN

  RETURN
END SUBROUTINE TERMNOT9


    !!
    !!  SUBROUTINE ESTERR
    !! ERRORS ESTIMATION.    ERRSAME: THE CURRENT ORDER
    !!                         ERRUP: GREATER ORDER (THAN ERRSAME)
    !!                       ERRDOWN: LOWER ORDER
    !!         

SUBROUTINE  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP,            &
     &     R, H, ORD, DBLK, LU, LDLU, FMAS, LDMAS, MLMAS, MUMAS, &
     &     IPIV, F, DN,SCAL,ORDMAX,ORDMIN,IJOB)


  !!
  !!   INPUT VARIABLES
  !!------------------------------------
  INTEGER, INTENT(IN) :: R, ORD,ORDMIN, ORDMAX,DBLK,IJOB,IPIV(R),LDLU, LDMAS, MLMAS, MUMAS
  REAL(PREC), INTENT(IN) :: H, SCAL(R),  LU(LDLU,R), FMAS(LDMAS,R)
  !!
  !!   INPUT/OUTPUT VARIABLES
  !!------------------------------------
  REAL(PREC), INTENT(IN OUT) ::  ERRV(DBLK), ERRSAME, ERRUP, ERRDOWN,   &
       &                  FP(R,11), F(R,11), DN(R,11)
  !!
  !!   LOCAL VARIABLES
  !!------------------------------------
  INTEGER  :: I, J
  REAL(PREC) ::  ERRVJ,   &
       &              FP1,FP2,FP3,FP4,FP5,FP6,FP7,FP8,FP9,FP10
  !!
  !!   EXECUTABLE STATEMENTS
  !!---------------------------------
  !!--------- ERRSAME ESTIMATION
  SELECT CASE(ORD)
  CASE(1)

     DO I=1,R
        FP1= FP(I,1)
        FP2= FP(I,2)
        FP3= FP(I,3)
        FP4= FP(I,4)
        FP5= FP(I,5)
        F(I,1) = H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
        F(I,2) = H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
        F(I,3) = H*(B3531*FP1+B3532*FP2+B3533*FP3+B3534*FP4+B3535*FP5)
        F(I,4) = H*(B3541*FP1+B3542*FP2+B3543*FP3+B3544*FP4+B3545*FP5)
       END DO
   
  CASE(2)
     DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
        F(J,1) = H*(B5711*FP1+B5712*FP2+B5713*FP3&
             &            +B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
        F(J,2) = H*(B5721*FP1+B5722*FP2+B5723*FP3&
             &            +B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
        F(J,3) = H*(B5731*FP1+B5732*FP2+B5733*FP3&
             &            +B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
        F(J,4) = H*(B5741*FP1+B5742*FP2+B5743*FP3&
             &            +B5744*FP4+B5745*FP5+B5746*FP6+B5747*FP7)
        F(J,5) = H*(B5727*FP1+B5726*FP2+B5725*FP3&
             &            +B5724*FP4+B5723*FP5+B5722*FP6+B5721*FP7)
        F(J,6) = H*(B5717*FP1+B5716*FP2+B5715*FP3&
             &            +B5714*FP4+B5713*FP5+B5712*FP6+B5711*FP7)
     END DO

  CASE(3)
     DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
        FP8= FP(J,8)
        FP9= FP(J,9)
        F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4&
             & +B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)
        F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4&
             & +B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)
        F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4&
             & +B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)
        F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4&
             & +B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)
        F(J,5) = H*(B7951*FP1+B7952*FP2+B7953*FP3+B7954*FP4&
             & +B7955*FP5+B7956*FP6+B7957*FP7+B7958*FP8+B7959*FP9)
        F(J,6) = H*(B7939*FP1+B7938*FP2+B7937*FP3+B7936*FP4&
             & +B7935*FP5+B7934*FP6+B7933*FP7+B7932*FP8+B7931*FP9)
        F(J,7) = H*(B7929*FP1+B7928*FP2+B7927*FP3+B7926*FP4&
             & +B7925*FP5+B7924*FP6+B7923*FP7+B7922*FP8+B7921*FP9)
        F(J,8) = H*(B7919*FP1+B7918*FP2+B7917*FP3+B7916*FP4&
             & +B7915*FP5+B7914*FP6+B7913*FP7+B7912*FP8+B7911*FP9)
     END DO

  CASE(4)
     DO J=1,R
        FP1= FP(J,1)
        FP2= FP(J,2)
        FP3= FP(J,3)
        FP4= FP(J,4)
        FP5= FP(J,5)
        FP6= FP(J,6)
        FP7= FP(J,7)
        FP8= FP(J,8)
        FP9= FP(J,9)
        FP10= FP(J,10)
        F(J,1) = H*(B91011*(FP1-FP10)+B91012*(FP2-FP9)+B91013*(FP3-FP8)&
             & +B91014*(FP4-FP7)+B91015*(FP5-FP6) )
        F(J,2) = H*(B91021*(FP1-FP10)+B91022*(FP2-FP9)+B91023*(FP3-FP8)&
             & +B91024*(FP4-FP7)+B91025*(FP5-FP6) )
        F(J,3) = H*(B91031*(FP1-FP10)+B91032*(FP2-FP9)+B91033*(FP3-FP8)&
             & +B91034*(FP4-FP7)+B91035*(FP5-FP6) )
        F(J,4) = H*(B91041*(FP1-FP10)+B91042*(FP2-FP9)+B91043*(FP3-FP8)&
             & +B91044*(FP4-FP7)+B91045*(FP5-FP6) )
        F(J,5) =  F(J,4)
        F(J,6) = -F(J,4)
        F(J,7) = -F(J,3)
        F(J,8) = -F(J,2)
        F(J,9) = -F(J,1)
     END DO

  END SELECT

  !!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRSAME

  CALL NEWTGS(R,DBLK,LU(1,1),LDLU,FMAS(1,1),LDMAS,MLMAS,MUMAS,&
     &      IPIV(1),F(1,1),DN(1,1),IJOB)


  !!--------- COMPUTE  ERRSAME AND ERRV (VECTOR ERROR)

  ERRSAME = 0D0
  DO J=1,DBLK
     ERRV(J) = 0D0
     DO I=1,R
        FP1 = (DN(I,J)/SCAL(I) )
        ERRV(J) =  ERRV(J)+ FP1*FP1
     END DO
     ERRV(J) = SQRT(ERRV(J)/R)
     ERRSAME = MAX( ERRSAME, ERRV(J) )
  END DO
  ERRSAME = MAX(ERRSAME, 1d-15)
  ERRDOWN = 0D0
  ERRUP   = 0D0
  IF ( ERRSAME .LE. 1d0) THEN

     IF (ORD .LT. ORDMAX) THEN
        !!--------- ERRUP ESTIMATION
        SELECT CASE(ORD)
        CASE(1)
           DO I=1,R
              FP1 =  F(I,1)/CP31
              FP2 =  F(I,2)/CP31
              FP3 =  F(I,3)/CP31
              FP4 = -F(I,4)/CP31
              F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP51
              F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP52
              F(I,3) = -(-FP2 + 2d0*FP3 - FP4)*CP52
              F(I,4) = -(-FP2 + 2d0*FP3 - FP4)*CP51
           END DO

        CASE(2)
           DO I = 1, R
              FP1 =  F(I,1)/CP51
              FP2 =  F(I,2)/CP52
              FP3 =  F(I,3)/CP52
              FP4 =  F(I,4)/CP52
              FP5 = -F(I,5)/CP52
              FP6 = -F(I,6)/CP51
              F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP71
              F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP72
              F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP73
              F(I,4) = -(-FP3 + 2d0*FP4 - FP5)*CP73
              F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP72
              F(I,6) = -(-FP4 + 2d0*FP5 - FP6)*CP71

           END DO
        CASE(3)
     
           DO I = 1, R
              FP1 =  F(I,1)/CP71
              FP2 =  F(I,2)/CP72
              FP3 =  F(I,3)/CP73
              FP4 =  F(I,4)/CP73
              FP5 =  F(I,5)/CP73
              FP6 = -F(I,6)/CP73
              FP7 = -F(I,7)/CP72
              FP8 = -F(I,8)/CP71
              F(I,1) =  (-FP1 + 2d0*FP2 - FP3)*CP91
              F(I,2) =  (-FP1 + 2d0*FP2 - FP3)*CP92
              F(I,3) =  (-FP2 + 2d0*FP3 - FP4)*CP93
              F(I,4) =  (-FP3 + 2d0*FP4 - FP5)*CP94
              F(I,5) = -(-FP4 + 2d0*FP5 - FP6)*CP94
              F(I,6) = -(-FP5 + 2d0*FP6 - FP7)*CP93
              F(I,7) = -(-FP6 + 2d0*FP7 - FP8)*CP92
              F(I,8) = -(-FP6 + 2d0*FP7 - FP8)*CP91
           END DO

        END SELECT

        !!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRUP


        CALL NEWTGS(R,DBLK,LU(1,1),LDLU,FMAS(1,1),LDMAS,MLMAS,MUMAS,IPIV(1),F(1,1),DN(1,1),IJOB)


        !!-------- COMPUTE ERRUP
        ERRUP = 0D0
        DO J=1,DBLK
           ERRVJ = 0D0
           DO I=1,R
              FP1 = (DN(I,J)/SCAL(I) )
              ERRVJ =  ERRVJ + FP1*FP1
           END DO
           ERRVJ = SQRT(ERRVJ/R)
           ERRUP = MAX( ERRUP, ERRVJ )
        END DO
        ERRUP = max(ERRUP, 1d-15)

     END IF

     IF (ORD .GT. ORDMIN) THEN
        !!--------- ERRDOWN ESTIMATION
        SELECT CASE(ORD)
        CASE(2)

           DO J=1,R
              FP1= FP(J,1)
              FP2= FP(J,2)
              FP3= FP(J,3)
              FP4= FP(J,4)
              FP5= FP(J,5)
              FP6= FP(J,6)
              FP7= FP(J,7)
              F(J,1) = H*(B3511*FP1+B3512*FP2+B3513*FP3+B3514*FP4+B3515*FP5)
              F(J,2) = H*(B3521*FP1+B3522*FP2+B3523*FP3+B3524*FP4+B3525*FP5)
              F(J,3) = H*(B3521*FP2+B3522*FP3+B3523*FP4+B3524*FP5+B3525*FP6)
              F(J,4) = H*(B3521*FP3+B3522*FP4+B3523*FP5+B3524*FP6+B3525*FP7)
              F(J,5) = H*(B3531*FP3+B3532*FP4+B3533*FP5+B3534*FP6+B3535*FP7)
              F(J,6) = H*(B3541*FP3+B3542*FP4+B3543*FP5+B3544*FP6+B3545*FP7)
           END DO
        CASE(3)

           DO J=1,R
              FP1= FP(J,1)
              FP2= FP(J,2)
              FP3= FP(J,3)
              FP4= FP(J,4)
              FP5= FP(J,5)
              FP6= FP(J,6)
              FP7= FP(J,7)
              FP8= FP(J,8)
              FP9= FP(J,9)
              F(J,1) = H*(B5711*FP1+B5712*FP2+B5713*FP3   &
                   &            +B5714*FP4+B5715*FP5+B5716*FP6+B5717*FP7)
              F(J,2) = H*(B5721*FP1+B5722*FP2+B5723*FP3   &
                   &            +B5724*FP4+B5725*FP5+B5726*FP6+B5727*FP7)
              F(J,3) = H*(B5731*FP1+B5732*FP2+B5733*FP3   &
                   &            +B5734*FP4+B5735*FP5+B5736*FP6+B5737*FP7)
              F(J,4) = H*(B5731*FP2+B5732*FP3+B5733*FP4   &
                   &            +B5734*FP5+B5735*FP6+B5736*FP7+B5737*FP8)
              F(J,5) = H*(B5731*FP3+B5732*FP4+B5733*FP5   &
                   &            +B5734*FP6+B5735*FP7+B5736*FP8+B5737*FP9)
              F(J,6) = H*(B5741*FP3+B5742*FP4+B5743*FP5   &
                   &            +B5744*FP6+B5745*FP7+B5746*FP8+B5747*FP9)
              F(J,7) = H*(B5727*FP3+B5726*FP4+B5725*FP5   &
                   &            +B5724*FP6+B5723*FP7+B5722*FP8+B5721*FP9)
              F(J,8) = H*(B5717*FP3+B5716*FP4+B5715*FP5   &
                   &            +B5714*FP6+B5713*FP7+B5712*FP8+B5711*FP9)

           END DO
        CASE(4)


           DO J=1,R
              FP1= FP(J,1)
              FP2= FP(J,2)
              FP3= FP(J,3)
              FP4= FP(J,4)
              FP5= FP(J,5)
              FP6= FP(J,6)
              FP7= FP(J,7)
              FP8= FP(J,8)
              FP9= FP(J,9)
              FP10= FP(J,10)
              F(J,1) = H*(B7911*FP1+B7912*FP2+B7913*FP3+B7914*FP4&
                   & +B7915*FP5+B7916*FP6+B7917*FP7+B7918*FP8+B7919*FP9)

              F(J,2) = H*(B7921*FP1+B7922*FP2+B7923*FP3+B7924*FP4&
                   & +B7925*FP5+B7926*FP6+B7927*FP7+B7928*FP8+B7929*FP9)

              F(J,3) = H*(B7931*FP1+B7932*FP2+B7933*FP3+B7934*FP4&
                   & +B7935*FP5+B7936*FP6+B7937*FP7+B7938*FP8+B7939*FP9)

              F(J,4) = H*(B7941*FP1+B7942*FP2+B7943*FP3+B7944*FP4&
                   & +B7945*FP5+B7946*FP6+B7947*FP7+B7948*FP8+B7949*FP9)

              F(J,5) = -H*(B7941*FP2+B7942*FP3+B7943*FP4+B7944*FP5&
                   & +B7945*FP6+B7946*FP7+B7947*FP8+B7948*FP9+B7949*FP10)


              F(J,6) = -H*(B7951*FP2+B7952*FP3+B7953*FP4+B7954*FP5&
                   & +B7955*FP6+B7956*FP7+B7957*FP8+B7958*FP9+B7959*FP10)

              F(J,7) = -H*(B7939*FP2+B7938*FP3+B7937*FP4+B7936*FP5&
                   & +B7935*FP6+B7934*FP7+B7933*FP8+B7932*FP9+B7931*FP10)

              F(J,8) = -H*(B7929*FP2+B7928*FP3+B7927*FP4+B7926*FP5&
                   & +B7925*FP6+B7924*FP7+B7923*FP8+B7922*FP9+B7921*FP10)

              F(J,9) = -H*(B7919*FP2+B7918*FP3+B7917*FP4+B7916*FP5&
                   & +B7915*FP6+B7914*FP7+B7913*FP8+B7912*FP9+B7911*FP10)

           END DO

        END SELECT
        !!--------- A SINGLE SPLITTING-NEWTON ITERATION FOR ERRDOWN

        CALL NEWTGS(R,DBLK,LU(1,1),LDLU,FMAS(1,1),LDMAS,MLMAS,MUMAS,IPIV(1),F(1,1),DN(1,1),IJOB)


        !!--------- COMPUTE ERRDOWN
        ERRDOWN = 0D0
        DO J=1,DBLK
           ERRVJ = 0D0
           DO I=1,R
              FP1 = (DN(I,J)/SCAL(I) )
              ERRVJ =  ERRVJ + FP1*FP1
           END DO
           ERRVJ = SQRT(ERRVJ/R)
           ERRDOWN = MAX(ERRDOWN, ERRVJ )
        END DO
        ERRDOWN = max(ERRDOWN, 1d-15)
     END IF

  END IF
  RETURN
END SUBROUTINE ESTERR


END MODULE SUBGAMD

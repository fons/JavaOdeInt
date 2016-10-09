C*************************************************************************
C     2/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 1 NAND GATE PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...   
      PARAMETER (ND=14,LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)   
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),
     +     ERROR(ND),RPAR(2),IPAR(2),T(0:16),
     +     MASBND(4),MBND(4),YPRIME(ND),jcount(15)
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'NAND GATE',
     +     IND=1)
C...
C... COMMON PARAMETERS
      double precision RGS, RGD, RBS, RBD, CGS, CGD,
     *                 CBD, CBS, C9,VT0, BETA, CGAMMA, DELTA, PHI,
     *                 CURIS, VTH, VDD,VBB
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB
C...
C...
      EXTERNAL PDERV,RESID
C...
      open(6,file='nandi.out')
      rewind(6)
      open(8,file='nandaci.out')     
      rewind(8) 
C
      write(6,2050)IND,PROBLEM,SOLVER
 2050 FORMAT(1X,'COMPUTATIONAL STATISTICS OF THE INDEX',I2,1X,A,
     &     ' PROBLEM USING ', A/)
      WRITE(6,2051)ND 
 2051 FORMAT(1X,'NUMBER OF EQUATIONS :',I3)
      WRITE(6,*)' '      
C...  
C... LOOP FOR DIFFERENT TOLERANCES
C...      
      NTOLMN=1
      NTOLMX=10
      NTOLDF=2
      NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
      TOLST=0.1D0**NTOLMN
      TOLFC=0.1D0**(1.D0/NTOLDF)      
      DO 30 NTOL=1,NRLOOP
         SUMM=0.0D+0
         SUMH=0.0D+0
C...  DIMENSION OF THE SYSTEM
         NEQN=14
C...  ENDPOINT OF INTEGRATION
         NDISC = 15
         T(0)=0.0D0         
         DO i = 1,NDISC + 1
            T(I) = DBLE(I)*5.D0
         ENDDO
C...  PROBLEM CONSTANT
         RGS    = .4D+1
         RGD    = .4D+1
         RBS    = .1D+2
         RBD    = .1D+2
         CGS    = .6D-4
         CGD    = .6D-4
         CBD    = 2.4D-5
         CBS    = 2.4D-5
         C9     = .5D-4
         DELTA  = 0.2D-1
         CURIS  = 1.D-14
         VTH    = 25.85D0
         VDD    = 5.d0
         VBB    = -2.5d0
C...  
C...  REQUIRED TOLERANCE
         RTOL =TOLST
         ATOL =RTOL                
         ITOL=2
         WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
         WRITE(6,*) 'RTOL =' ,RTOL
         WRITE(6,*) 'ATOL =' ,ATOL
C...  INITIAL VALUES
         CALL INIT(NEQN,T(0),Y,YPRIME)
C...  SET DEFAULT VALUES 
         MF=22
         INDEX=1
         LOUT=6
         MAXDER=7
         H0=1.d-2*RTOL
         H = H0
C... MAXIMAL NUMBER OF STEPS 
         iwork(14)=100000
         XOUT=T(NDISC+1)
         iwork(1)=14
         iwork(2)=0
         iwork(3)=0
         ierr = 0         
         work(1) = 0.d0
         XOUT=T(1)         
         mbnd(1) = neqn
         mbnd(2) = neqn
         mbnd(3) = neqn
         mbnd(4) = neqn
         do k=5,12
            jcount(k) =0
         enddo
         
         it1=mclock()
         TT = T(0)
C...
         DO 40 I=0,NDISC            
            H = H0
            XOUT = T(I+1)       
 220     CONTINUE
         
C...  CALL OF THE SUBROUTINE 
C...
         CALL MEBDFI(NEQN,TT,H,Y,YPRIME,XOUT,T(NDISC+1),MF,INDEX,LOUT,
     +        LWORK,WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,
     +        ATOL,RPAR,IPAR,PDERV, RESID,IERR)
C...
c         WRITE(*,*) 'INDEX =' , INDEX, 'IERR =',IERR
         IF (INDEX .EQ. 1) THEN
            INDEX = 2
            GOTO 220
         ELSEIF (INDEX .NE. 0) THEN
            WRITE(*,*) 'MEBDFI RETURN INDEX = ',INDEX
            WRITE(*,*)
c            stop
C...  GO TO THE NEXT TOLERANCE
            GOTO 25
         ENDIF 
C     
         index = 1
         TT = XOUT
         do k=5,12
            jcount(k) = jcount(k) + iwork(k)
         enddo
 40      CONTINUE
C
         it2=mclock()
         TIME=(it2-it1)/100.0D+0         
C...
C... PRINT AND CALCULATE THE ERROR AT THE END POINT
C...         
         CALL SOLN(NEQN,TRUE)
         CALL SCD(NEQN,Y,TRUE,AERR12,RERR12,US2)                 
C...  
C...  STATISTICS OF THE SIMULATION         
C...
         do k=5,12
            iwork(k) = jcount(k)
         enddo
         HUSED = WORK(2)
         NQUSED = IWORK(4)
         NSTEP  = IWORK(5)
         NFAIL  = IWORK(6)
         NFE    = IWORK(7)
         NJE    = IWORK(8)
         NDEC   = IWORK(9)
         NBSOL  = IWORK(10)
         NPSET  = IWORK(11)
         NCOSET = IWORK(12)
         MAXORD = IWORK(13)
C...
C...  THE CURRENT RUN IS COMPLETE, SO PRINT THE COMPUTATIONAL STAT-
C...  ISTICS FOR MEBDF AND GO ON TO THE NEXT RUN
C...
         WRITE(LOUT,8)HUSED,NQUSED,MAXORD,NSTEP,NFAIL,NFE,NJE,NDEC,
     +        NBSOL,US2,TIME
 8       FORMAT(1H /
     1 ' LAST STEP SIZE                       ',              D20.9,/,
     2 ' LAST ORDER OF THE METHOD             ',                I10,/,
     3 ' MAXIMUM ORDER USED SO FAR            ',                I10,/,
     4 ' TOTAL NUMBER OF STEPS TAKEN          ',                I10,/, 
     5 ' TOTAL NUMBER OF FAILED STEPS         ',                I10,/,
     6 ' NUMBER OF FUNCTION EVALUATIONS       ',                I10,/,
     7 ' NUMBER OF JACOBIAN EVALUATIONS       ',                I10,/,
     8 ' NUMBER OF FACTORIZATION              ',                I10,/,
     9 ' NUMBER OF BACKSOLVES                 ',                I10,/,
     + ' NUMBER OF CORRECT DIGITS             ',              F10.2,/,   
     + ' CPU TIME                             ',                F10.4) 
         WRITE(LOUT,*) '---------------------------------------------'    
         WRITE(8,*) US2,TIME
C...         
C...  NEW TOLERANCE
C...         
 25      TOLST=TOLST*TOLFC
 30   CONTINUE
      close(6)
      close(8)
      STOP
      END
C----------------------------------------------------------------------
      SUBROUTINE RESID(NEQN,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)           
      DOUBLE PRECISION T,Y(NEQN),DELTA(NEQN),YPRIME(NEQN),RPAR(*)
      integer nnodes,i,j,nj,ipar(*)      
      double precision am(14,14),fy(14),dum

      CALL   MAS(NEQN,AM,NEQN,Y,IPAR,RPAR,IERR)
      CALL   F(NEQN,T,Y,FY,IPAR,RPAR,IERR)

      if(ierr.eq.-1)return

      do 20 i=1,14
         dum = -fy(i)
         do 10 j=1,14
            dum = dum+AM(i,j)*yprime(j)
   10    continue
         DELTA(I) = dum
   20 continue
C
      RETURN
      END 
C-----------------------------------------------------------------------
      SUBROUTINE PDERV(T,Y,PD,NEQN,YPRIME,M,IPAR,RPAR,IERR,CON)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION Y(NEQN),PD(*),YPRIME(NEQN)
c
c     dummy subroutine
c
      RETURN
      END
C------------------------------------------------------------------------      
      subroutine init(neqn,t,y,yprime)
      integer neqn,i
      double precision t,y(neqn),yprime(neqn) 
c
      double precision VBB

      VBB = -2.5d0

      y(1)  = 5d0
      y(2)  = 5d0
      y(3)  = VBB
      y(4)  = VBB
      y(5)  = 5d0
      y(6)  = 3.62385d0
      y(7)  = 5d0
      y(8)  = VBB
      y(9)  = VBB
      y(10) = 3.62385d0
      y(11) = 0d0
      y(12) = 3.62385d0
      y(13) = VBB
      y(14) = VBB

      do 10 i=1,14
         yprime(i) = 0d0
   10 continue
      return
      end
c-------------------------------------------------------------------------
      subroutine soln(neqn,y)
      double precision y(neqn)
c
c
c computed at Cray C90, using Cray double precision:
C Solving NAND gate using PSIDE
C
C User input:
C
C give relative error tolerance: 1d-16
C give absolute error tolerance: 1d-16
C
C
C Integration characteristics:
C
C    number of integration steps       22083
C    number of accepted steps          21506
C    number of f evaluations          308562
C    number of Jacobian evaluations      337
C    number of LU decompositions       10532
C
C CPU-time used:                         451.71 sec

      y(  1) =  0.4971088699385777d+1
      y(  2) =  0.4999752103929311d+1
      y(  3) = -0.2499998781491227d+1
      y(  4) = -0.2499999999999975d+1
      y(  5) =  0.4970837023296724d+1
      y(  6) = -0.2091214032073855d+0
      y(  7) =  0.4970593243278363d+1
      y(  8) = -0.2500077409198803d+1
      y(  9) = -0.2499998781491227d+1
      y( 10) = -0.2090289583878100d+0
      y( 11) = -0.2399999999966269d-3
      y( 12) = -0.2091214032073855d+0
      y( 13) = -0.2499999999999991d+1
      y( 14) = -0.2500077409198803d+1
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE F(N,T,Y,DF,IPAR,RPAR,IERR)
C ---------------------------------------------------------------------
C
C Right-hand side f(X,t) for the network equation
C             C(Y) * Y' - f(Y,t) = 0
C describing the nand gate
C
C ---------------------------------------------------------------------
C
C Input parameters:
C          N......number of node potentials (14)
C          T......time point t
C          Y......node potentials at time point t
C Output parameter:
C          F......right-hand side f(Y,t)
C
C External reference:
C          IDS: Drain-source current
C          IBS: Nonlinear current characteristic for diode between
C               bulk and source
C          IBD: Nonlinear current characteristic for diode between
C               bulk and drain
C          PULSE: Input signal in pulse form
C
C ---------------------------------------------------------------------

      INTEGER N,ierr
      double precision T,Y(N),DF(N) ,IDS,IBS,IBD, V1,V2,V1D,V2D
      EXTERNAL IDS, IBS, IBD, PULSE

      double precision RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                 VT0, BETA, CGAMMA, DELTA, PHI,
     *                 CURIS, VTH, VDD,VBB
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB
      CALL PULSE(T,V1,V1D,0.D0,5.D0,5.D0,5.D0,5.D0,5.D0,20.D0)
      CALL PULSE(T,V2,V2D,0.D0,5.D0,15.D0,5.D0,15.D0,5.D0,40.D0)


      DF(1)=-(Y(1)-Y(5))/RGS - IDS(1,Y(2)-Y(1),Y(5)-Y(1),Y(3)-Y(5),
     *       Y(5)-Y(2),Y(4)-VDD,ierr)
      DF(2)=-(Y(2)-VDD)/RGD + IDS(1,Y(2)-Y(1),Y(5)-Y(1),Y(3)-Y(5),
     *       Y(5)-Y(2),Y(4)-VDD,ierr)
      DF(3)=-(Y(3)-VBB)/RBS + IBS(Y(3)-Y(5))
      DF(4)=-(Y(4)-VBB)/RBD + IBD(Y(4)-VDD)
      DF(5)=-(Y(5)-Y(1))/RGS- IBS(Y(3)-Y(5)) - (Y(5)-Y(7))/RGD -
     *       IBD(Y(9)-Y(5))
      DF(6)=CGS*V1D-(Y(6)-Y(10))/RGS-
     *       IDS(2,Y(7)-Y(6),V1-Y(6),Y(8)-Y(10),V1-Y(7),Y(9)-Y(5),ierr)
      DF(7)=CGD*V1D-(Y(7)-Y(5))/RGD+
     *       IDS(2,Y(7)-Y(6),V1-Y(6),Y(8)-Y(10),V1-Y(7),Y(9)-Y(5),ierr)
      DF(8)=-(Y(8)-VBB)/RBS + IBS(Y(8)-Y(10))
      DF(9)=-(Y(9)-VBB)/RBD + IBD(Y(9)-Y(5))
      DF(10)=-(Y(10)-Y(6))/RGS - IBS(Y(8)-Y(10)) -
     *         (Y(10)-Y(12))/RGD - IBD(Y(14)-Y(10))
      DF(11)=CGS*V2D-Y(11)/RGS-IDS(2,Y(12)-Y(11),V2-Y(11),Y(13),
     *       V2-Y(12),Y(14)-Y(10),ierr)
      DF(12)=CGD*V2D-(Y(12)-Y(10))/RGD+
     *       IDS(2,Y(12)-Y(11),V2-Y(11),Y(13),V2-Y(12),Y(14)-Y(10),ierr)
      DF(13)=-(Y(13)-VBB)/RBS + IBS(Y(13))
      DF(14)=-(Y(14)-VBB)/RBD + IBD(Y(14)-Y(10))

      if(ierr.eq.-1)return
      RETURN
      END
C---------------------------------------------------------------------------
C
      double precision FUNCTION IDS (NED,VDS, VGS, VBS, VGD, VBD, ierr)
C ---------------------------------------------------------------------------
C
C Function evaluating the drain-current due to the model of
C Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   NED  Integer parameter for MOSFET-type
C   VDS  Voltage between drain and source
C   VGS  Voltage between gate and source
C   VBS  Voltage between bulk and source
C   VGD  Voltage between gate and drain
C   VBD  Voltage between bulk and drain
C
C External reference:
C   GDSP, GDSM Drain function for VDS > 0 gevalp. VDS < 0
C
C ---------------------------------------------------------------------------

      INTEGER NED,ierr
      double precision VDS, VGS, VBS, VGD, VBD,GDSP, GDSM
      EXTERNAL GDSP, GDSM

      IF ( VDS .GT. 0.D0 ) THEN
       IDS = GDSP (NED,VDS, VGS, VBS,ierr)
      ELSE IF ( VDS .EQ. 0.D0) THEN
       IDS = 0.D0
      ELSE IF ( VDS .LT. 0.D0) THE N
       IDS = GDSM (NED,VDS, VGD, VBD,ierr)
      END IF

      if(ierr.eq.-1)return

      RETURN
      END
C ---------------------------------------------------------------------------
      double precision FUNCTION GDSP (NED,VDS, VGS, VBS, ierr)
      integer NED,ierr
      double precision  VDS, VGS, VBS,VTE

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB
C

      IF(NED.EQ.1) THEN
C --- Depletion-type
      VT0=-2.43D0
      CGAMMA=.2D0
      PHI=1.28D0
      BETA=5.35D-4
      ELSE
C --- Enhancement-type
      VT0=.2D0
      CGAMMA=0.035D0
      PHI=1.01D0
      BETA=1.748D-3
      END IF

      if(phi-vbs.lt.0d0.or.phi.lt.0d0)then
         ierr=-1
         return
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBS) - DSQRT(PHI) )

      IF ( VGS-VTE .LE. 0.D0) THEN
       GDSP = 0.D0
      ELSE IF ( 0.D0 .LT. VGS-VTE .AND. VGS-VTE .LE. VDS ) THEN
       GDSP = - BETA * (VGS - VTE)**2.D0 * (1.D0 + DELTA*VDS)
      ELSE IF ( 0.D0 .LT. VDS .AND. VDS .LT. VGS-VTE ) THEN
       GDSP = - BETA * VDS * (2.D0*(VGS - VTE) - VDS) *
     *          (1.D0 + DELTA*VDS)
      END IF

      RETURN
      END
C ---------------------------------------------------------------------------
      double precision FUNCTION GDSM (NED,VDS, VGD, VBD, ierr)
      integer NED,ierr
      double precision VDS, VGD, VBD,VTE


      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  VT0, BETA, CGAMMA, DELTA, PHI,
     *                  CURIS, VTH, VDD,VBB

      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB

      IF(NED.EQ.1) THEN
C --- Depletion-type
      VT0=-2.43D0
      CGAMMA=.2D0
      PHI=1.28D0
      BETA=5.35D-4
      ELSE
C --- Enhancement-type
      VT0=.2D0
      CGAMMA=0.035D0
      PHI=1.01D0
      BETA=1.748D-4
      END IF

      if(phi-vbd.lt.0d0.or.phi.lt.0d0)then
         ierr=-1
         return
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBD) - DSQRT(PHI) )

      IF ( VGD-VTE .LE. 0.D0) THEN
       GDSM = 0.D0
      ELSE IF ( 0.D0 .LT. VGD-VTE .AND. VGD-VTE .LE. -VDS ) THEN
       GDSM = BETA * (VGD - VTE)**2d0 * (1.D0 - DELTA*VDS)
      ELSE IF ( 0.D0 .LT. -VDS .AND. -VDS .LT. VGD-VTE ) THEN
       GDSM = - BETA * VDS * (2d0 *(VGD - VTE) + VDS) *
     *          (1.D0 - DELTA*VDS)
      END IF

      RETURN
      END
C ---------------------------------------------------------------------------

      double precision FUNCTION IBS (VBS)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C source due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and source
C
C ---------------------------------------------------------------------------

      double precision VBS
      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB

C

C
C     IBS = GBS (VBS)
C

      IF ( VBS .LE. 0.D0 ) THEN
       IBS = - CURIS * ( DEXP( VBS/VTH ) - 1.D0 )
      ELSE
       IBS = 0.D0
      END IF

      RETURN
      END
C ---------------------------------------------------------------------------

      double precision FUNCTION IBD (VBD)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C drain  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and drain
C
C ---------------------------------------------------------------------------

      double precision VBD
      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  VT0, BETA, CGAMMA, DELTA, PHI,
     *                  CURIS, VTH, VDD,VBB

      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB


C

C
C     IBD = GBD (VBD)
C
      IF ( VBD .LE. 0.D0 ) THEN
       IBD = - CURIS * ( DEXP( VBD/VTH ) - 1.D0 )
      ELSE
       IBD = 0.D0
      END IF
      RETURN
      END
C ---------------------------------------------------------------------------

      SUBROUTINE PULSE(X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD)
C ---------------------------------------------------------------------------
C
C Evaluating input signal at time point X
C
C Structure of input signal:
C
C                -----------------------                       HIGH
C               /                       \
C              /                         \
C             /                           \
C            /                             \
C           /                               \
C          /                                 \
C         /                                   \
C        /                                     \
C  ------                                       ---------      LOW
C
C |DELAY|   T1  |         T2           |   T3  |
C |          P     E     R     I     O     D            |
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   X                      Time-point at which input signal is evaluated
C   LOW                    Low-level of input signal
C   HIGH                   High-level of input signal
C   DELAY,T1,T2,T3, PERIOD Parameters to specify signal structure
C
C Output parameter:
C   VIN    Voltage of input signal at time point X
C   VIND   Derivative of VIN at time point X
C
C ---------------------------------------------------------------------------

      double precision X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD,TIME

      TIME = DMOD(X,PERIOD)

      IF (TIME.GT.(DELAY+T1+T2+T3)) THEN
      VIN = LOW
      VIND= 0.D0
      ELSE IF (TIME.GT.(DELAY+T1+T2)) THEN
      VIN = ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW
      VIND= -((HIGH-LOW)/T3)
      ELSE IF (TIME.GT.(DELAY+T1)) THEN
      VIN = HIGH
      VIND= 0.D0
      ELSE IF (TIME.GT.DELAY) THEN
      VIN = ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW
      VIND= ((HIGH-LOW)/T1)
      ELSE
      VIN = LOW
      VIND=0.D0
      END IF

      RETURN
      END
C --------------------------------------------------------------------------
      SUBROUTINE MAS(N,AM,M,Y,IPAR,RPAR,IERR)
C      SUBROUTINE CAP(N,Y,AM)
C ---------------------------------------------------------------------
C
C Voltage-dependent capacitance matrix C(Y) for the network equation
C             C(Y) * Y' - f(Y,t) = 0
C describing the nand gate
C
C ---------------------------------------------------------------------
C
C Input parameters:
C          N......number of node potentials (14)
C          Y......value of node potentials
C Output parameter:
C          AM.....voltage-dependent capacitance matrix
C
C External reference:
C          CBDBS: Voltage-dependent capacitance CBS(V) and CBD(V)
C
C ---------------------------------------------------------------------
      double precision CBDBS
      INTEGER N
      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *     VT0, BETA, CGAMMA, DELTA, PHI,
     *     CURIS, VTH, VDD,VBB
      double precision Y(N), AM(N,N)
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *     VT0, BETA, CGAMMA, DELTA, PHI,
     *     CURIS, VTH, VDD,VBB
      EXTERNAL CBDBS
      integer I,J
      
      DO 10 I=1,N
         DO 20 J=1,N
            AM(I,J)=0d0
 20      CONTINUE
 10   CONTINUE
      
      AM(1,1)= CGS
      AM(1,5)=-CGS
      AM(2,2)= CGD
      AM(2,5)=-CGD
      AM(3,3)= CBDBS(Y(3)-Y(5))
      AM(3,5)=-CBDBS(Y(3)-Y(5))
      AM(4,4)= CBDBS(Y(4)-VDD)
      AM(5,1)=-CGS
      AM(5,2)=-CGD
      AM(5,3)=-CBDBS(Y(3)-Y(5))
      AM(5,5)= CGS+CGD-AM(5,3)+
     *     CBDBS(Y(9)-Y(5)) +C9
      AM(5,9)=-CBDBS(Y(9)-Y(5))
      AM(6,6)= CGS
      AM(7,7)= CGD
      AM(8,8)= CBDBS(Y(8)-Y(10))
      AM(8,10)=-CBDBS(Y(8)-Y(10))
      AM(9,5)=-CBDBS(Y(9)-Y(5))
      AM(9,9)=CBDBS(Y(9)-Y(5))
      AM(10,8)=-CBDBS(Y(8)-Y(10))
      AM(10,10)=-AM(8,10)+CBDBS(Y(14)-Y(10))+C9
      AM(10,14)= -CBDBS(Y(14)-Y(10))
      AM(11,11)=CGS
      AM(12,12)=CGD
      AM(13,13)=CBDBS(Y(13))
      AM(14,10)=-CBDBS(Y(14)-Y(10))
      AM(14,14)= CBDBS(Y(14)-Y(10))

      RETURN
      END
C ---------------------------------------------------------------------------

      double precision FUNCTION CBDBS (V)
C ---------------------------------------------------------------------------
C
C Function evaluating the voltage-dependent capacitance between bulk and
C drain gevalp. source  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   V    Voltage between bulk and drain gevalp. source
C
C ---------------------------------------------------------------------------
      double precision V,PHIB

      double precision  RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *                  VT0, BETA, CGAMMA, DELTA, PHI,
     *                  CURIS, VTH, VDD,VBB

      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS, C9,
     *               VT0, BETA, CGAMMA, DELTA, PHI,
     *               CURIS, VTH, VDD,VBB

      PHIB=0.87D0

      IF ( V .LE. 0.D0 ) THEN
       CBDBS = CBD/DSQRT(1.D0-V/PHIB)
      ELSE
       CBDBS = CBD*(1.D0+V/(2.D0*PHIB))
      END IF

      RETURN
      END
C-----------------------------------------------------------------------------

      SUBROUTINE SCD(NEQN,Y,TRUE,AERR12,RERR12,US2)
      DOUBLE PRECISION Y(NEQN),TRUE(NEQN),AERR12,RERR12
      DOUBLE PRECISION AERRMX,RERRMX,AERR,RERR,US2
      LOGICAL SKIPI
      CHARACTER*1 SKIPC
C
      WRITE(*,*)
      WRITE(*,100) 
 100  FORMAT(' NUMERICAL SOLUTION')
      WRITE(*,*)
      WRITE(*,101)
 101  FORMAT(47X,'SCD')
      WRITE(*,102)
 102  FORMAT(5X,'SOLUTION COMPONENT',18X,'-------------    IGNORE')
      WRITE(*,103)
 103  FORMAT(41X,'ABS        REL ')
       
      NUMA = 0
      NUMR = 0
      AERR12 = 0D0
      RERR12 = 0D0
      AREEMX = 0D0
      RERRMX = 0D0
C            
      DO I=1,NEQN
         SKIPI =  I .NE. 5

         IF (SKIPI) THEN
            SKIPC='Y'
         ELSE
            SKIPC =' '
         ENDIF

         AERR = DABS(TRUE(I)-Y(I))
         IF (.NOT. SKIPI) THEN
            AERRMX = DMAX1(AERRMX,AERR)
            AERR12 = AERR12 + AERR**2
            NUMA = NUMA +1
         ENDIF
C
         IF (TRUE(I) .NE. 0.D0) THEN
            RERR = DABS((TRUE(I)-Y(I))/TRUE(I))
            IF (.NOT. SKIPI) THEN
               RERRMX = DMAX1(RERRMX,RERR)
               RERR12 = RERR12 + RERR**2
               NUMR = NUMR +1
            ENDIF
         ENDIF
C          
         IF (AERR .EQ. 0D0)THEN
            WRITE(*,111) I,Y(I), SKIPC
 111        format(I3,2X,e24.16,30X,a)
         ELSEIF (TRUE(I) .EQ. 0D0) THEN
            WRITE(*,121) I,Y(I),-LOG10(AERR), SKIPC
 121        FORMAT(I3,2X,e24.16,6X,F10.2,6X,a)
         ELSE
            write(*,122) i,y(i),-log10(aerr),-log10(rerr), SKIPC
 122        FORMAT(I3,2X,E24.16,6X,F10.2,1X,F10.2,6X,A)
         ENDIF

      ENDDO
      WRITE(*,*)
C
      AERR12 = DSQRT(AERR12)
      RERR12 = DSQRT(RERR12)
C
      WRITE(*,112) NUMA, NUMR
 112  format('USED COMPONENTS FOR SCD',12X,I10,1X,I10)
      WRITE(*,113) -log10(AERRMX),-log10(rerrmx) 
 113  FORMAT('SCD OF Y (MAXIMUM NORM)',12X,F10.2,1X,F10.2)
      WRITE(*,114) -log10(AERR12),-log(rerr12)
 114  FORMAT('SCD OF Y (EUCLIDIAN NORM)',10X,F10.2,1X,F10.2)           
            
      US2 = -LOG10(RERRMX)
      WRITE(*,115) US2
 115  FORMAT('USING RELATIVE ERROR YIELDS SCD : ', 12X,F10.2)
C
      RETURN
      END
C-----------------------------------------------------------------------

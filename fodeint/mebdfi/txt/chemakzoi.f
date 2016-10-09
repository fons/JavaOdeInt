C*************************************************************************
C     01/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 1 CHEMICAL AKZO  PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...        
      PARAMETER (ND=6,LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)      
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),
     +     ERROR(ND),RPAR(2),IPAR(2),MBND(4),YPRIME(ND)             
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'CHEMICAL AKZO',
     +     IND=1)
C...
      EXTERNAL PDERV, RESID
C...
      open(6,file='akzoi.out')
      rewind(6) 
      open(8,file='akzoaci.out')
      rewind(8)

      write(6,2050)IND,PROBLEM,SOLVER
 2050 FORMAT(1X,'COMPUTATIONAL STATISTICS OF THE INDEX',I2,1X,A,
     &     ' PROBLEM USING ', A/)
      WRITE(6,2051)ND 
 2051 FORMAT(1X,'NUMBER OF EQUATIONS :',I2)
      WRITE(6,*)' '
C...  
C... LOOP FOR DIFFERENT TOLERANCES
C...      
      NTOLMN=2
      NTOLMX=10
      NTOLDF=2
      NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
      TOLST=0.1D0**NTOLMN
      TOLFC=0.1D0**(1.D0/NTOLDF)        
      DO 30 NTOL=1,NRLOOP         
C...  DIMENSION OF THE SYSTEM
         NEQN=6
C...  ENDPOINT OF INTEGRATION    
         X=0.0D0             
         XEND=180.0D0         
C...  REQUIRED TOLERANCE
         RTOL=TOLST
         ATOL=RTOL
         ITOL=2
         WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
         WRITE(6,*) 'RTOL =' ,RTOL
         WRITE(6,*) 'ATOL =' ,ATOL
C...  INITIAL VALUES          
         CALL INIT(NEQN,X,Y,YPRIME)         
C...  SET DEFAULT VALUES 
         MF=21
         INDEX=1
         LOUT=6
         MAXDER=7
C...  MAXIMAL NUMBER OF STEPS         
         iwork(14)=100000          
         H = RTOL
         XOUT=180.0D0
         IWORK(1)= NEQN
         IWORK(2)=0
         IWORK(3)=0
         MBND(1) = NEQN
         MBND(2) = NEQN
         MBND(3) = NEQN
         MBND(4) = NEQN
         it1=mclock()        
      
 220     CONTINUE
C...  CALL OF THE SUBROUTINE           
         CALL MEBDFI(NEQN,X,H,Y,YPRIME,XOUT,XEND,MF,INDEX,LOUT, 
     +        LWORK,WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,ATOL,
     +        RPAR,IPAR,PDERV,RESID,IERR)
c         WRITE(*,*) 'INDEX =', INDEX, 'IERR =', IERR
         IF (INDEX .EQ. 1) THEN         
            INDEX = 0 
            GOTO 220
         ELSEIF (INDEX .NE. 0) THEN
            WRITE(*,*) 'MEBDFI RETURN INDEX = ', INDEX
c            stop
C...  GO TO THE NEXT TOLERANCE
            GOTO 25
         ENDIF
C
         it2=mclock()
         TIME=(it2-it1)/100.0D+0
C...
C...  CALCULATE AND PRINT THE ERROR AT THE END POINT
C...                       
         CALL SOLN(NEQN,TRUE)
         CALL SCD(NEQN,Y,TRUE,AERR12,RERR12,US2)                  
C...  
C...  STATISTICS OF THE SIMULATION         
C...
         HUSED  = WORK(2)
         NQUSED = IWORK(4)
         NSTEP  = IWORK(5)
         NFAIL  = IWORK(6)
         NRE    = IWORK(7)
         NJE    = IWORK(8)
         NDEC   = IWORK(9)
         NBSOL  = IWORK(10)
         NPSET  = IWORK(11)
         NCOSET = IWORK(12)
         MAXORD = IWORK(13)
C...
C...  THE CURRENT RUN IS COMPLETE, SO PRINT THE COMPUTATIONAL STAT-
C...  ISTICS FOR MEBDF AND GO ON TO THE NEXT RUN
         WRITE(LOUT,8)HUSED,NQUSED,MAXORD,NSTEP,NFAIL,NRE,NJE,NDEC,
     +        NBSOL,US2,TIME
 8       FORMAT(1H /
     1 ' LAST STEP SIZE                       ',              D13.6,/,
     2 ' LAST ORDER OF THE METHOD             ',                I10,/,
     3 ' MAXIMUM ORDER USED SO FAR            ',                I10,/,
     4 ' TOTAL NUMBER OF STEPS TAKEN          ',                I10,/, 
     5 ' TOTAL NUMBER OF FAILED STEPS         ',                I10,/,
     6 ' NUMBER OF RESIDUAL EVALUATIONS       ',                I10,/,
     7 ' NUMBER OF JACOBIAN EVALUATIONS       ',                I10,/,
     8 ' NUMBER OF FACTORIZATION              ',                I10,/,
     9 ' NUMBER OF BACKSOLVES                 ',                I10,/,
     + ' NUMBER OF CORRECT DIGITS           ',              F13.2,/,   
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
C-----------------------------------------------------------------------------
      SUBROUTINE RESID(N,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
      DOUBLE PRECISION T, Y(N), DELTA(N), YPRIME(N), RPAR(2)
      INTEGER N, IPAR(2)
C
      CALL F(N,T,Y,DELTA,IPAR,RPAR,IERR)
C
      DO I=1, N-1
         DELTA(I) = YPRIME(I) - DELTA(I)
      ENDDO
      DELTA(N) = -DELTA(N)
C
      RETURN
      END
C-----------------------------------------------------------------------------
      SUBROUTINE F(NEQN,X,Y,DF,IPAR,RPAR,IERR)
      DOUBLE PRECISION  X,Y(NEQN),DF(NEQN),RPAR(*)    
      INTEGER NEQN, IPAR(*)
      double precision k1,k2,k3,k4,kbig,kla,po2,hen,ks      
      parameter (
     +     k1   = 18.7d0,
     +     k2   = 0.58d0,
     +     k3   = 0.09d0,
     +     k4   = 0.42d0,
     +     kbig = 34.4d0,
     +     kla  = 3.3d0,
     +     ks   = 115.83d0,
     +     po2  = 0.9d0,
     +     hen  = 737d0
     +     )
      double precision r1,r2,r3,r4,r5,fin
C...      
      if (y(2) .lt. 0d0) then
         ierr = -1
         return
      endif
      
      r1  = k1*(y(1)**4)*sqrt(y(2))
      r2  = k2*y(3)*y(4)
      r3  = k2/kbig*y(1)*y(5)
      r4  = k3*y(1)*(y(4)**2)
      r5  = k4*(y(6)**2)*sqrt(y(2))
      fin = kla*(po2/hen-y(2))
      
      df(1) =   -2d0*r1 +r2 -r3     -r4
      df(2) = -0.5d0*r1             -r4     -0.5d0*r5 + fin
      df(3) =        r1 -r2 +r3
      df(4) =           -r2 +r3 -2d0*r4
      df(5) =            r2 -r3         +r5
      df(6) = ks*y(1)*y(4)-y(6)
C
      RETURN
      END
C-----------------------------------------------------------------------       
      SUBROUTINE PDERV(T,Y,DFDY,NEQN,YPRIME,M,CON,IPAR,RPAR,IERR)
      double precision T,Y(NEQN),DFDY(NEQN,NEQN),YPRIME(NEQN)      
      double precision k1,k2,k3,k4,kbig,kla,ks,RPAR(*),con
      parameter (
     +     k1   =18.7d0,
     +     k2   =0.58d0,
     +     k3   =0.09d0,
     +     k4   =0.42d0,
     +     kbig =34.4d0,
     +     kla  =3.3d0,
     +     ks   =115.83d0
     +     )
      integer i,j,ipar(*)
      double precision r11,r12,r23,r24,r31,r35,r41,r44,r52,r56,fin2
C      
      if (y(2) .lt. 0d0) then
         ierr = -1
         return
      endif
      r11  = 4d0*k1*(y(1)**3)*sqrt(y(2))
      r12  = 0.5d0*k1*(y(1)**4)/sqrt(y(2))
      r23  = k2*y(4)
      r24  = k2*y(3)
      r31  = (k2/kbig)*y(5)
      r35  = (k2/kbig)*y(1)
      r41  = k3*y(4)**2
      r44  = 2d0*k3*y(1)*y(4)
      r52  = 0.5d0*k4*(y(6)**2)/sqrt(y(2))
      r56  = 2d0*k4*y(6)*sqrt(y(2))
      fin2 = -kla
      dfdy(1,1) = -2d0*r11-r31-r41
      dfdy(1,2) = -2d0*r12
      dfdy(1,3) = r23
      dfdy(1,4) = r24-r44
      dfdy(1,5) = -r35
      dfdy(2,1) = -0.5d0*r11-r41
      dfdy(2,2) = -0.5d0*r12-0.5d0*r52+fin2
      dfdy(2,4) = -r44
      dfdy(2,6) = -0.5d0*r56
      dfdy(3,1) = r11+r31
      dfdy(3,2) = r12
      dfdy(3,3) = -r23
      dfdy(3,4) = -r24
      dfdy(3,5) = r35
      dfdy(4,1) = r31-2d0*r41
      dfdy(4,3) = -r23
      dfdy(4,4) = -r24-2d0*r44
      dfdy(4,5) = r35
      dfdy(5,1) = -r31
      dfdy(5,2) = r52
      dfdy(5,3) = r23
      dfdy(5,4) = r24
      dfdy(5,5) = -r35
      dfdy(5,6) = r56
      dfdy(6,1) = ks*y(4)
      dfdy(6,4) = ks*y(1)
      dfdy(6,6) = -1d0

      do i=1,neqn
         do j=1,neqn
            dfdy(i,j) =-dfdy(i,j)            
         enddo
      enddo
c compute pd = -df/dy + con*M
      do i=1,neqn-1
         dfdy(i,i) = 1.0d0/con+dfdy(i,i)
      enddo
      
      RETURN
      END
C-----------------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime)
      integer neqn
      double precision t,y(neqn),yprime(neqn)            
      double precision k1,k2,k3,k4,kbig,kla,ks
      parameter (ks   =115.83d0)

      y(1) = 0.444d0
      y(2) = 0.00123d0
      y(3) = 0d0
      y(4) = 0.007d0
      y(5) = 0d0
      y(6) = ks*y(1)*y(4)
c
      call f(neqn,0d0,y,yprime,ipar,rpar,ierr)
      
      RETURN
      END
C-----------------------------------------------------------------------------
      SUBROUTINE SOLN(N,TRUE)
      DOUBLE PRECISION TRUE(N)
      INTEGER N 
C
      TRUE(1)= 0.1150794920661702d0
      TRUE(2)= 0.1203831471567715d-2
      TRUE(3)= 0.1611562887407974d0
      TRUE(4)= 0.3656156421249283d-3
      TRUE(5)= 0.1708010885264404d-1
      TRUE(6)= 0.4873531310307455d-2
C
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
         SKIPI = .FALSE.

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

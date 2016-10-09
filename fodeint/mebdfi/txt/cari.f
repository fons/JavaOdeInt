C*************************************************************************
C     01/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 3 CAR AXIS PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFV
      PARAMETER (ND=10,LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),MBND(4),
     +     ERROR(ND),RPAR(2),IPAR(2),YPRIME(ND)
        CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'CAR AXIS',
     +     IND=3)
C...
      EXTERNAL PDERV, RESID
C        
      open(6,file='cari.out')
      rewind(6)
      open(8,file='caraci.out')
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
      NTOLMN=2
      NTOLMX=10
      NTOLDF=2
      NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
      TOLST=0.1D0**NTOLMN
      TOLFC=0.1D0**(1.D0/NTOLDF)        
      DO 30 NTOL=1,NRLOOP        
C...  DIMENSION OF THE SYSTEM
         N=10
C...  ENDPOINT OF INTEGRATION
         X=0.0D0        
         XEND=3.0D0
C...  REQUIRED TOLERANCE
         RTOL=TOLST
         ATOL=RTOL
         ITOL=2
         WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
         WRITE(6,*) 'RTOL =' ,RTOL
         WRITE(6,*) 'ATOL =' ,ATOL
C...  INITIAL VALUES         
         CALL INIT(N,X,Y,YPRIME)
C...  SET DEFAULT VALUES 
         MF=22
         INDEX=1
         LOUT=6                          
         MAXDER=7
C...  MAXIMAL NUMBER OF STEPS
         IWORK(14)=100000
         H=RTOL        
         XOUT=3.0D0
         iwork(1)=4
         iwork(2)=4
         iwork(3)=2      
         mbnd(1) = n
         mbnd(2) = n
         mbnd(3) = n
         mbnd(4) = n
C
        it1=mclock()
C     
 220    continue
C...  CALL OF THE SUBROUTINE  
        CALL MEBDFI(N,X,H,Y,YPRIME,XOUT,XEND,MF,INDEX,LOUT,LWORK, 
     +       WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,ATOL,
     +       RPAR,IPAR,PDERV,RESID,IERR)
c        WRITE(*,*) ' INDEX = ', INDEX, ' IERR = ', IERR
        IF (INDEX .EQ. 1) THEN
           INDEX = 0 
           GOTO 220
        ELSEIF (INDEX .NE. 0) THEN
           WRITE(*,*) 'MEBDFI RETURN INDEX = ',INDEX
c            stop
C...  GO TO THE NEXT TOLERANCE
           GOTO 25
        ENDIF 
C     
        it2=mclock()
        TIME=(it2-it1)/100.0D+0
C...
C... PRINT AND CALCULATE THE ERROR AT THE END POINT
C...                
        CALL SOLN(N,X,TRUE)
        CALL SCD(N,Y,TRUE,AERR12,RERR12,US2)                              
C...  
C...  STATISTICS OF THE SIMULATION         
C...      
        HUSED  =WORK(2)
        NQUSED =IWORK(4)
        NSTEP  =IWORK(5)
        NFAIL  =IWORK(6)
        NRE    =IWORK(7)         
        NJE    =IWORK(8)         
        NDEC   =IWORK(9)
        NBSOL  =IWORK(10)         
        MAXORD =IWORK(13) 
C...
C...  THE CURRENT RUN IS COMPLETE, SO PRINT THE COMPUTATIONAL STAT-
C...  ISTICS FOR MEBDF AND GO ON TO THE NEXT RUN
C...        
        WRITE(LOUT,8)HUSED,NQUSED,MAXORD,NSTEP,NFAIL,NRE,NJE,NDEC,
     +       NBSOL,US2, TIME
 8       FORMAT(1H /
     1 ' LAST STEP SIZE                       ',              F13.6,/,
     2 ' LAST ORDER OF THE METHOD             ',                I10,/,
     3 ' MAXIMUM ORDER USED SO FAR            ',                I10,/,
     4 ' TOTAL NUMBER OF STEPS TAKEN          ',                I10,/, 
     5 ' TOTAL NUMBER OF FAILED STEPS         ',                I10,/,
     6 ' NUMBER OF RESIDUAL EVALUATIONS       ',                I10,/,
     7 ' NUMBER OF JACOBIAN EVALUATIONS       ',                I10,/,
     8 ' NUMBER OF FACTORIZATION              ',                I10,/,
     9 ' NUMBER OF BACKSOLVES                 ',                I10,/,  
     + ' NUMBER OF CORRECT DIGITS             ',              F13.2,/,   
     + ' CPU TIME                             ',                F10.4) 
         WRITE(LOUT,*) '---------------------------------------------'      
         write(8,*) US2,TIME
C...
C...  NEW TOLERANCE 
C...
 25      TOLST=TOLST*TOLFC
 30   CONTINUE
      CLOSE(6)
      CLOSE(8)
      STOP
      END
c----------------------------------------------------------------------
      SUBROUTINE RESID(N,X,Y,DELTA,YPRIME,IPAR,RPAR,ier)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),DELTA(N),IPAR(*),RPAR(*),YPRIME(N)
      double precision m,eps, k
C
      M = 10d0
      eps = 1d-2
      k = m*eps*eps/2d0
      
      call feval(n,x,y,yprime,delta,ierr,rpar,ipar)
      do i=1,4
         delta(i) =yprime(i)-delta(i)
      enddo   
      do i=5,8
         delta(i) = k*yprime(i)- delta(i)
      enddo
      do i=9,10
         delta(i) = -delta(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine feval(neqn,t,y,yprime,f,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*), i
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*)      
      double precision M,L,L0,w,r,xb,yb,Ll,Lr,eps,g,
     +                 xl,yl,xr,yr,lam1,lam2

      eps = 1d-2
      M   = 10d0
      L   = 1d0
      L0  = 0.5d0
      r   = 0.1d0
      w   = 10d0
      g   = 1d0
      yb  = r*sin(w*t)
      xb  = sqrt(L*L-yb*yb)

      do 10 i=1,4
         f(i) = y(i+4)
   10 continue

      xl   = y(1)
      yl   = y(2)
      xr   = y(3)
      yr   = y(4)
      lam1 = y(9)
      lam2 = y(10)

      Ll = sqrt(xl**2+yl**2)
      Lr = sqrt((xr-xb)**2+(yr-yb)**2)

      f(5)  =(L0-Ll)*xl/Ll +lam1*xb+2d0*lam2*(xl-xr)
      f(6)  =(L0-Ll)*yl/Ll +lam1*yb+2d0*lam2*(yl-yr)-M*eps*eps*g/2d0
      f(7)  =(L0-Lr)*(xr-xb)/Lr -2d0*lam2*(xl-xr)
      f(8)  =(L0-Lr)*(yr-yb)/Lr -2d0*lam2*(yl-yr)-M*eps*eps*g/2d0

      f(9)  = xb*xl+yb*yl
      f(10) = (xl-xr)**2+(yl-yr)**2-L*L

      return
      end
c-----------------------------------------------------------------------      
        SUBROUTINE PDERV(X,Y,DFY,N,YPRIME,MEBAND,CON,IPAR,RPAR,IERR)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),DFY(MEBAND,N),IPAR(*),RPAR(*),yprime(*)
c
c     dummy subroutine
c        
        RETURN
        END
c-------------------------------------------------------------------------     
      subroutine init(neqn,t,y,yprime)
      integer neqn, i, ierr, ipar
      double precision t,y(neqn),yprime(neqn), rpar            
      double precision M,eps,L,L0,
     +                 xl,yl,xr,yr,xla,yla,xra,yra,lam1,lam2
            

      M   = 10d0
      eps = 1d-2
      L   = 1d0
      L0  = 0.5d0

      xr   = L
      xl   = 0
      yr   = L0
      yl   = yr
      xra  = -L0/L
      xla  = xra
      yra  = 0d0
      yla  = 0d0
      lam1 = 0d0
      lam2 = 0d0

      y(1)  =  xl
      y(2)  =  yl
      y(3)  =  xr
      y(4)  =  yr
      y(5)  =  xla
      y(6)  =  yla
      y(7)  =  xra
      y(8)  =  yra
      y(9)  =  lam1
      y(10) =  lam2

      call feval(neqn,0d0,y,y,yprime,ierr,rpar,ipar)
      do 10 i=1,4
         yprime(i) = y(i+4)
         yprime(i+4)= 2d0/(M*eps*eps)*yprime(i+4)
 10   continue
      do 20 i=9,10
         yprime(i)=0d0
 20   continue
      
      return
      end
c-----------------------------------------------------------------------
      subroutine soln(neqn,t,y)
      integer neqn
      double precision t,y(neqn)

      y(  1) =  0.4934557842755629d-001
      y(  2) =  0.4969894602303324d+000
      y(  3) =  0.1041742524885400d+001
      y(  4) =  0.3739110272652214d+000
      y(  5) = -0.7705836840321485d-001
      y(  6) =  0.7446866596327776d-002
      y(  7) =  0.1755681574942899d-001
      y(  8) =  0.7703410437794031d+000
      y(  9) = -0.4736886750784630d-002
      y( 10) = -0.1104680411345730d-002

      return
      end
c----------------------------------------------------------------------------
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
         SKIPI =  .FALSE.

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












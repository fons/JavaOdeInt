C*************************************************************************
C     01/03/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 1 TRANSISTOR  PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...
      PARAMETER (ND=8,LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)             
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),MBND(4),
     +     ERROR(ND),RPAR(2),IPAR(2),YPRIME(ND)             
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'TRANSISTOR',IND=1)      
C...
      EXTERNAL PDERV, RESID
C...     
      open(6,file='transi.out')
      rewind(6)
      open(8,file='transaci.out')
      rewind(8)
c      
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
         SUMM=0.0D+0
         SUMH=0.0D+0
C...  DIMENSION OF THE SYSTEM
         NEQN=8
C...  ENDPOINT OF INTEGRATION
         X=0.0D0                  
         XEND=0.2D+0
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
         MF=22
         INDEX=1
         LOUT=6
         MAXDER=7
C...  MAXIMAL NUMBER OF STEPS
         iwork(14)=200000 
         H=1.0D-2*RTOL
         XOUT=0.2D+0
         IWORK(1)=neqn
         IWORK(2)=0
         IWORK(3)=0
         mbnd(1) = neqn
         mbnd(2) = neqn
         mbnd(3) = neqn
         mbnd(4) = neqn                  
C...        
         it1=mclock()

 220     CONTINUE
C...  CALL OF THE SUBROUTINE          
        CALL MEBDFI(NEQN,X,H,Y,YPRIME,XOUT,XEND,MF,INDEX,LOUT,LWORK, 
     +        WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,ATOL,
     +        RPAR,IPAR,PDERV,RESID,IERR)
c        write(*,*) 'INDEX = ',INDEX,' IERR=',IERR                   
        if (index .eq. 1) then
            index = 0 
            goto 220
         elseif (index .ne. 0) THEN
            write(*,*) 'MEBDFI return index = ',index
c            stop
c...  go to the next tolerance
            goto 25
         ENDIF

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
         WRITE(*,8)HUSED,NQUSED,MAXORD,NSTEP,NFAIL,NRE,NJE,NDEC,
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
     + ' NUMBER OF CORRECT DIGITS             ',              F10.2,/,   
     + ' CPU TIME                             ',                F10.3) 
         WRITE(*,*) '---------------------------------------------'    
         WRITE(8,*) US2,TIME
C...         
C...  NEW TOLERANCE
C...
 25     TOLST=TOLST*TOLFC
 30   CONTINUE
      CLOSE(6)
      CLOSE(8)
      STOP
      END
C----------------------------------------------------------------------------
      SUBROUTINE RESID(NEQN,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),DELTA(NEQN),YPRIME(NEQN),IPAR(2),RPAR(2)
      double precision c1,c2,c3,c4,c5
      PARAMETER (C1=1D-6,C2=2D-6,C3=3D-6,C4=4D-6,C5=5D-6)
C
      CALL F(NEQN,T,Y,DELTA,IPAR,RPAR,IERR)
C
      DELTA(1) = -C1*YPRIME(1)+C1*YPRIME(2) -DELTA(1)
      DELTA(2) =  C1*YPRIME(1)-C1*YPRIME(2) -DELTA(2)
      DELTA(3) = -C2*YPRIME(3)              -DELTA(3)
      DELTA(4) = -C3*YPRIME(4)+C3*YPRIME(5) -DELTA(4)
      DELTA(5) =  C3*YPRIME(4)-C3*YPRIME(5) -DELTA(5)
      DELTA(6) = -C4*YPRIME(6)              -DELTA(6)
      DELTA(7) = -C5*YPRIME(7)+C5*YPRIME(8) -DELTA(7)
      DELTA(8) =  C5*YPRIME(7)-C5*YPRIME(8) -DELTA(8)
      
      RETURN
      END
C----------------------------------------------------------------------------
      subroutine f(neqn,t,y,dy,ipar,rpar,ierr)
      IMPLICIT double precision (A-H,O-Z)      
      double precision t,y(neqn),dy(neqn)      
      double precision ub,uf,alpha,beta,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,
     +     pi,uet,fac1,fac2
      parameter (ub=6d0,uf=0.026d0,alpha=0.99d0,beta=1d-6,
     +     r0=1000d0,r1=9000d0,r2=9000d0,r3=9000d0,
     +     r4=9000d0,r5=9000d0,r6=9000d0,r7=9000d0,
     +     r8=9000d0,r9=9000d0,pi=3.1415926535897931086244d0)
      
      uet   = 0.1d0*sin(200d0*pi*t)
      fac1  = beta*(exp((y(2)-y(3))/uf)-1d0)
      fac2  = beta*(exp((y(5)-y(6))/uf)-1d0)
      
      dy(1) = (y(1)-uet)/r0
      dy(2) = y(2)/r1+(y(2)-ub)/r2+(1d0-alpha)*fac1
      dy(3) = y(3)/r3-fac1
      dy(4) = (y(4)-ub)/r4+alpha*fac1
      dy(5) = y(5)/r5+(y(5)-ub)/r6+(1d0-alpha)*fac2
      dy(6) = y(6)/r7-fac2
      dy(7) = (y(7)-ub)/r8+alpha*fac2
      dy(8) = y(8)/r9
      
      return
      end
c-----------------------------------------------------------------------
      subroutine pderv(T,y,dfdy,n,yprime,mebnd,con,ipar,rpar,ierr)
      IMPLICIT double precision (A-H,O-Z)      
      double precision t,y(n),dfdy(mebnd,n),yprime(n),rpar(*)     
      integer n,i,ipar(*)
      double precision uf,alpha,beta,r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,
     +     fac1p,fac2p
      parameter (uf=0.026d0,alpha=0.99d0,beta=1d-6,
     +     r0=1000d0,r1=9000d0,r2=9000d0,r3=9000d0,
     +     r4=9000d0,r5=9000d0,r6=9000d0,r7=9000d0,
     +     r8=9000d0,r9=9000d0)
c
      fac1p = beta*dexp((y(2)-y(3))/uf)/uf
      fac2p = beta*dexp((y(5)-y(6))/uf)/uf
      do 10 i=1,8
         dfdy(1,i) = 0d0
         dfdy(3,i) = 0d0
         dfdy(4,i) = 0d0
   10 continue

      dfdy(1,3) = -(1d0-alpha)*fac1p
      dfdy(1,6) = -(1d0-alpha)*fac2p
      dfdy(2,1) = 1d0/r0
      dfdy(2,2) = 1d0/r1+1d0/r2+(1d0-alpha)*fac1p
      dfdy(2,3) = 1d0/r3+fac1p
      dfdy(2,4) = 1d0/r4
      dfdy(2,5) = 1d0/r5+1d0/r6+(1d0-alpha)*fac2p
      dfdy(2,6) = 1d0/r7+fac2p
      dfdy(2,7) = 1d0/r8
      dfdy(2,8) = 1d0/r9
      dfdy(3,2) = -fac1p
      dfdy(3,3) = -alpha*fac1p
      dfdy(3,5) = -fac2p
      dfdy(3,6) = -alpha*fac2p
      dfdy(4,2) = alpha*fac1p
      dfdy(4,5) = alpha*fac2p
      return
      end
C-----------------------------------------------------------------------
      SUBROUTINE INIT(NEQN,T,Y,YPRIME)
      DOUBLE PRECISION T,Y(NEQN),YPRIME(NEQN)
      double precision ub,r1,r2,r3,r5,r6,r7,c2,c4
      parameter (ub=6d0,r1=9000d0,r2=9000d0,r3=9000d0,
     +     r5=9000d0,r6=9000d0,r7=9000d0,c2=2d-6,c4=4d-6)
C
      Y(1) = 0d0
      Y(2) = ub/(r2/r1+1d0)
      Y(3) = Y(2)
      Y(4) = ub
      Y(5) = ub/(r6/r5+1d0)
      Y(6) = Y(5)
      Y(7) = Y(4)
      Y(8) = 0d0
C     
      yprime(3) = -y(2)/(c2*r3)
      yprime(6) = -y(5)/(c4*r7)
c-----------------------------------------------------------------------
c     the other initial values for yprime are determined numerically
c-----------------------------------------------------------------------
      yprime(1) = 51.338775d0
      yprime(2) = 51.338775d0
      yprime(4) = -24.9757667d0
      yprime(5) = -24.9757667d0
      yprime(7) = -10.00564453d0
      yprime(8) = -10.00564453d0
      RETURN
      END
C----------------------------------------------------------------------
      subroutine soln(neqn,y)
      double precision y(neqn)
C
      y(  1) = -0.5562145012262709d-002
      y(  2) =  0.3006522471903042d+001
      y(  3) =  0.2849958788608128d+001
      y(  4) =  0.2926422536206241d+001
      y(  5) =  0.2704617865010554d+001
      y(  6) =  0.2761837778393145d+001
      y(  7) =  0.4770927631616772d+001
      y(  8) =  0.1236995868091548d+001
C
      RETURN
      END
C-----------------------------------------------------------------------
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
 102  FORMAT(5X,'SOLUTION COMPONENT',18X,'-------------    IGNOR')
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

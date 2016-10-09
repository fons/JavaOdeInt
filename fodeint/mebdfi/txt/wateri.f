C*************************************************************************
C     01/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 2 WATER TUBE PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...   
      PARAMETER (ND=49,LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)   
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),
     +     ERROR(ND),RPAR(2),IPAR(2),ATOL(ND),RTOL(ND),T(0:2),
     +     MBND(4),YPRIME(ND)
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'WATER TUBE',
     +     IND=2)
C...
      EXTERNAL PDERV, RESID
C...
      open(6,file='wateri.out')
      rewind(6)
      open(8,file='wateraci.out')     
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
         NEQN=49
C...  ENDPOINT OF INTEGRATION
         T(0)=0.0D0         
         T(1)=17.0D0*3600D0
C...  REQUIRED TOLERANCE
         RTOL(1)=TOLST
         ATOL(1)=RTOL(1)
         DO I=2,49
            ATOL(I) = ATOL(1)
            RTOL(I) = RTOL(1)
         ENDDO
         DO I= 37,49
            ATOL(I) = 1D6*ATOL(1)
         ENDDO        
         ITOL=5
         WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
         WRITE(6,*) 'RTOL =' ,RTOL(1)
         WRITE(6,*) 'ATOL =' ,ATOL(1), 'and ', ATOL(37)
C...  INITIAL VALUES
         CALL INIT(NEQN,T(0),Y,YPRIME)
C...  SET DEFAULT VALUES 
         MF=22
         INDEX=1
         LOUT=6
         MAXDER=7
         H=RTOL(1)        
C... MAXIMAL NUMBER OF STEPS 
         iwork(14)=100000
         XOUT=17.0D0*3600D0
         iwork(1)=38
         iwork(2)=11
         iwork(3)=0
         ierr = 0         
         work(1) = 0.d0
         XOUT=T(1)         
         mbnd(1) = neqn
         mbnd(2) = neqn
         mbnd(3) = neqn
         mbnd(4) = neqn
         it1=mclock()
                  
 220     CONTINUE
         
C...  CALL OF THE SUBROUTINE  
         CALL MEBDFI(NEQN,T(0),H,Y,YPRIME,XOUT,T(1),MF,INDEX,LOUT,
     +        LWORK,WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,
     +        ATOL,RPAR,IPAR,PDERV,RESID,IERR)
c         write(*,*) 'INDEX =',INDEX,'IERR = ',IERR         
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
C...
         WRITE(LOUT,8)HUSED,NQUSED,MAXORD,NSTEP,NFAIL,NRE,NJE,NDEC,
     +        NBSOL,US2,TIME
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
c-----------------------------------------------------------------------
      SUBROUTINE RESID(N,X,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),DELTA(N),IPAR(*),RPAR(*),YPRIME(N)
      double precision nu,g,rho,rcrit,length,k,d,b,pi,a,c,v,x
c      
      parameter(nu    = 1.31d-6,
     +          g     = 9.8d0,
     +          rho   = 1.0d3,
     +          rcrit = 2.3d3,
     +          length= 1.0d3,
     +          k     = 2.0d-4,
     +          d     = 1.0d0,
     +          b     = 2.0d2,
     +          pi    = 3.141592653589793238462643383d0 )            

      a = pi*d**2/4d0
      c = b/(rho*g)
      v = rho*length/a
      
      call f(n,x,y,delta,ipar,rpar,ierr)
C
      do i=1,18
         delta(i) =v*yprime(i)-delta(i)         
      enddo   
      do i=19,36
         delta(i) = -delta(i)         
      enddo
      do i=37,38
         delta(i) = c*yprime(i)-delta(i)
      enddo
      do i=39,n
         delta(i) = -delta(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE F(NEQN,T,Y,DF,IPAR,RPAR,IERR)      
      DOUBLE PRECISION T,Y(NEQN),DF(NEQN),RPAR(*)
      integer nnodes,i,j,nj,ipar(*)
      parameter(nnodes=13)
      double precision phi(nnodes,nnodes),lambda(nnodes,nnodes),
     +                 p(nnodes),
     +                 fdba(nnodes,nnodes),netflo(nnodes),
     +                 ein(nnodes),eout(nnodes),rghres(nnodes,nnodes)
      double precision rtla,r,that,that2,a,mu
      double precision nu,g,rho,rcrit,length,k,d,b,pi
      parameter(nu    = 1.31d-6,
     +          g     = 9.8d0,
     +          rho   = 1.0d3,
     +          rcrit = 2.3d3,
     +          length= 1.0d3,
     +          k     = 2.0d-4,
     +          d     = 1.0d0,
     +          b     = 2.0d2,
     +          pi    = 3.141592653589793238462643383d0 )                  

      a = pi*d**2/4d0
      mu = nu*rho

      do 5 j=1,nnodes
         ein(j) = 0d0
         eout(j) = 0d0
    5 continue
      that=t/3600d0
      that2=that*that

      ein(1) = (1d0-cos(exp(-that)-1d0))/200d0
      ein(13) = (1d0-cos(exp(-that)-1d0))/80d0
      eout(10) = that2*(3d0*that2-92d0*that+720d0)/1d6

      do 8 j=1,nnodes
         do 7 i=1,nnodes
            phi(i,j) = 0d0
            lambda(i,j) = 1d0
    7    continue
    8 continue
c      write(*,*) 't=',t
c      do i=1,neqn
c         write(*,*) i,'y=',y(i)
c      enddo
c      write(*,*)

      phi( 1, 2) = y( 1)
      phi( 2, 3) = y( 2)
      phi( 2, 6) = y( 3)
      phi( 3, 4) = y( 4)
      phi( 3, 5) = y( 5)
      phi( 4, 5) = y( 6)
      phi( 5,10) = y( 7)
      phi( 6, 5) = y( 8)
      phi( 7, 4) = y( 9)
      phi( 7, 8) = y(10)
      phi( 8, 5) = y(11)
      phi( 8,10) = y(12)
      phi( 9, 8) = y(13)
      phi(11, 9) = y(14)
      phi(11,12) = y(15)
      phi(12, 7) = y(16)
      phi(12, 8) = y(17)
      phi(13,11) = y(18)      

      lambda( 1, 2) = y(19)
      lambda( 2, 3) = y(20)
      lambda( 2, 6) = y(21)
      lambda( 3, 4) = y(22)
      lambda( 3, 5) = y(23)
      lambda( 4, 5) = y(24)
      lambda( 5,10) = y(25)
      lambda( 6, 5) = y(26)
      lambda( 7, 4) = y(27)
      lambda( 7, 8) = y(28)
      lambda( 8, 5) = y(29)
      lambda( 8,10) = y(30)
      lambda( 9, 8) = y(31)
      lambda(11, 9) = y(32)
      lambda(11,12) = y(33)
      lambda(12, 7) = y(34)
      lambda(12, 8) = y(35)
      lambda(13,11) = y(36)

      p( 5) = y(37)
      p( 8) = y(38)
      p( 1) = y(39)
      p( 2) = y(40)
      p( 3) = y(41)
      p( 4) = y(42)
      p( 6) = y(43)
      p( 7) = y(44)
      p( 9) = y(45)
      p(10) = y(46)
      p(11) = y(47)
      p(12) = y(48)
      p(13) = y(49)

      do 20 j=1,nnodes         
         do 10 i=1,nnodes            
            if (lambda(i,j) .lt. 0d0) then
               ierr = -1
               return
            endif
            rtla=sqrt(lambda(i,j))
            r = abs(phi(i,j)*d/(nu*a))
            if (r.gt.rcrit) then
               rghres(i,j) = 1/rtla - 1.74d0 +
     +                   2d0*log10(2d0*k/d + 18.7d0/(r*rtla))
               fdba(i,j) = p(i) - p(j) -
     +                 lambda(i,j)*rho*length*phi(i,j)**2/(a**2*d)
            else
               rghres(i,j) = 1.d0/rtla - 1.74d0 +
     +                   2d0*log10(2d0*k/d + 18.7d0/(rcrit*rtla))
c              rghres(i,j)=lambda(i,j)-0.47519404529185289807d-1
               fdba(i,j) = p(i) - p(j) -
     +                     32d0*mu*length*phi(i,j)/(a*d**2)
            endif
   10    continue
   20 continue

      do 50 nj=1,nnodes
         netflo(nj) = ein(nj)-eout(nj)
         do 30 i=1,nnodes
            netflo(nj) = netflo(nj)+phi(i,nj)
   30    continue
         do 40 j=1,nnodes
            netflo(nj) = netflo(nj)-phi(nj,j)
   40    continue
   50 continue

      df( 1) = fdba( 1, 2)
      df( 2) = fdba( 2, 3)
      df( 3) = fdba( 2, 6)
      df( 4) = fdba( 3, 4)
      df( 5) = fdba( 3, 5)
      df( 6) = fdba( 4, 5)
      df( 7) = fdba( 5,10)
      df( 8) = fdba( 6, 5)
      df( 9) = fdba( 7, 4)
      df(10) = fdba( 7, 8)
      df(11) = fdba( 8, 5)
      df(12) = fdba( 8,10)
      df(13) = fdba( 9, 8)
      df(14) = fdba(11, 9)
      df(15) = fdba(11,12)
      df(16) = fdba(12, 7)
      df(17) = fdba(12, 8)
      df(18) = fdba(13,11)

      df(19) = rghres( 1, 2)
      df(20) = rghres( 2, 3)
      df(21) = rghres( 2, 6)
      df(22) = rghres( 3, 4)
      df(23) = rghres( 3, 5)
      df(24) = rghres( 4, 5)
      df(25) = rghres( 5,10)
      df(26) = rghres( 6, 5)
      df(27) = rghres( 7, 4)
      df(28) = rghres( 7, 8)
      df(29) = rghres( 8, 5)
      df(30) = rghres( 8,10)
      df(31) = rghres( 9, 8)
      df(32) = rghres(11, 9)
      df(33) = rghres(11,12)
      df(34) = rghres(12, 7)
      df(35) = rghres(12, 8)
      df(36) = rghres(13,11)

      df(37) = netflo( 5)
      df(38) = netflo( 8)
      df(39) = netflo( 1)
      df(40) = netflo( 2)
      df(41) = netflo( 3)
      df(42) = netflo( 4)
      df(43) = netflo( 6)
      df(44) = netflo( 7)
      df(45) = netflo( 9)
      df(46) = netflo(10)
      df(47) = netflo(11)
      df(48) = netflo(12)
      df(49) = netflo(13)
      
      RETURN
      END 
C-----------------------------------------------------------------------------
      SUBROUTINE PDERV(T,Y,PD,NEQN,YPRIME,M,CON,IPAR,RPAR,IERR)
      IMPLICIT double precision (A-H,O-Z)
      DIMENSION Y(NEQN),PD(*),YPRIME(NEQN)
c
c     dummy subroutine
c
      RETURN
      END
C----------------------------------------------------------------------------
      subroutine soln(neqn,y)
      double precision y(neqn)
c
      y(  1) =  0.2298488296477430d-002
      y(  2) =  0.1188984650746585d-002
      y(  3) =  0.1109503645730845d-002
      y(  4) =  0.1589620100314825d-003
      y(  5) =  0.1030022640715102d-002
      y(  6) =  0.8710606306836165d-003
      y(  7) =  0.3243571480903489d-002
      y(  8) =  0.1109503645730845d-002
      y(  9) =  0.7120986206521341d-003
      y( 10) =  0.6414613963833099d-003
      y( 11) =  0.9416978549524347d-003
      y( 12) =  0.3403428519096511d-002
      y( 13) =  0.2397639310739395d-002
      y( 14) =  0.2397639310739395d-002
      y( 15) =  0.3348581430454180d-002
      y( 16) =  0.1353560017035444d-002
      y( 17) =  0.1995021413418736d-002
      y( 18) =  0.5746220741193575d-002
      y( 19) =  0.4751940452918529d-001
      y( 20) =  0.4751940452918529d-001
      y( 21) =  0.4751940452918529d-001
      y( 22) =  0.4751940452918529d-001
      y( 23) =  0.4751940452918529d-001
      y( 24) =  0.4751940452918529d-001
      y( 25) =  0.4311196778792902d-001
      y( 26) =  0.4751940452918529d-001
      y( 27) =  0.4751940452918529d-001
      y( 28) =  0.4751940452918529d-001
      y( 29) =  0.4751940452918529d-001
      y( 30) =  0.4249217433601160d-001
      y( 31) =  0.4732336439609648d-001
      y( 32) =  0.4732336439609648d-001
      y( 33) =  0.4270002118868241d-001
      y( 34) =  0.4751940452918529d-001
      y( 35) =  0.4751940452918529d-001
      y( 36) =  0.3651427026675656d-001
      y( 37) =  0.1111268591478108d+006
      y( 38) =  0.1111270045592387d+006
      y( 39) =  0.1111271078730254d+006
      y( 40) =  0.1111269851929858d+006
      y( 41) =  0.1111269255355337d+006
      y( 42) =  0.1111269322658045d+006
      y( 43) =  0.1111269221703983d+006
      y( 44) =  0.1111270121140691d+006
      y( 45) =  0.1111274419515807d+006
      y( 46) =  0.1111255158881087d+006
      y( 47) =  0.1111278793439227d+006
      y( 48) =  0.1111270995171642d+006
      y( 49) =  0.1111298338971779d+006
      return
      end
c-------------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime)
      integer neqn,i
      double precision t,y(neqn),yprime(neqn) 
c
      do 10 i=1,neqn
         y(i) = 0d0         
         yprime(i) = 0d0
 10   continue
      do 20 i=19,36
         y(i) = 0.47519404529185289807d-1
 20   continue
      do 30 i=37,49
         y(i) = 109800d0
 30   continue
      return
      end
c-----------------------------------------------------------------------
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
c         IF (AERR .EQ. 0D0)THEN
c            WRITE(*,111) I,Y(I), SKIPC
c 111        format(I3,2X,e24.16,30X,a)
c         ELSEIF (TRUE(I) .EQ. 0D0) THEN
c            WRITE(*,121) I,Y(I),-LOG10(AERR), SKIPC
c 121        FORMAT(I3,2X,e24.16,6X,F10.2,6X,a)
c         ELSE
c            write(*,122) i,y(i),-log10(aerr),-log10(rerr), SKIPC
c 122        FORMAT(I3,2X,E24.16,6X,F10.2,1X,F10.2,6X,A)
c         ENDIF

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

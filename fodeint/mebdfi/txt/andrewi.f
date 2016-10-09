C*************************************************************************
C     02/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 3 ANDREW SQUEEZING  PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...
      PARAMETER (ND=27, LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)   
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),
     +     ERROR(ND),RPAR(2),IPAR(2),YPRIME(ND),MBND(4)             
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'ANDREW SQUEEZING',
     +     IND=3)
C...
      EXTERNAL PDERV, RESID
C...
      open(6,file='andrewi.out')
      rewind(6)
      open(8,file='andrewaci.out')     
      rewind(8)
C
      write(6,2050)IND,PROBLEM,SOLVER
 2050 FORMAT(1X,'COMPUTATIONAL STATISTICS OF THE INDEX',I2,1X,A,
     &     ' PROBLEM USING ', A/)
      WRITE(6,2051)ND 
 2051 FORMAT(1X,'NUMBER OF EQUATIONS :',I5)
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
         NEQN=27  
C...  ENDPOINT OF INTEGRATION
         X=0.0D0         
         XEND=3.0D-2       
C... REQUIRED TOLERANCE
         RTOL = TOLST
         ATOL=RTOL
         ITOL=2
         WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
         WRITE(6,*) 'RTOL =' ,RTOL
         WRITE(6,*) 'ATOL =' ,ATOL
C... INITIAL VALUES        
         CALL INITIAL(NEQN,X,Y,YPRIME)
C...  SET DEFAULT VALUES 
         MF=21
         INDEX=1
         LOUT=6
         MAXDER=7
C...  MAXIMAL NUMBER OF STEPS  
         IWORK(14)=100000
         H=RTOL 
         XOUT=3.0d-2
         IWORK(1)=7
         IWORK(2)=7
         IWORK(3)=13                  
         MBND(1) = NEQN
         MBND(2) = NEQN
         MBND(3) = NEQN
         MBND(4) = NEQN
         it1=mclock()

 220     CONTINUE
C...  CALL OF THE SUBROUTINE           
C... 
         CALL MEBDFI(NEQN,X,H,Y,YPRIME,XOUT,XEND,MF,INDEX,LOUT,LWORK, 
     +        WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,ATOL,
     +        RPAR,IPAR,PDERV,RESID,IERR)
c         WRITE(*,*) 'INDEX = ',INDEX,' IERR =',IERR
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
         CALL SOLN(NEQN,X,TRUE)
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
      STOP
      END
c-----------------------------------------------------------------------
      SUBROUTINE RESID(N,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
      DOUBLE PRECISION T, Y(N), DELTA(N), YPRIME(N),RPAR(2)
      INTEGER N, IPAR(2)
C
      CALL F(N,T,Y,DELTA,IPAR,RPAR, IERR)
C
      DO J = 1,14
         DELTA(J) = YPRIME(J)-DELTA(J)
      ENDDO
      DO I=15,N
         DELTA(I) = -DELTA(I)
      ENDDO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE F(NEQN,T,Y,DY,IPAR,RPAR,IERR)
      DOUBLE PRECISION T, Y(NEQN), DY(NEQN)         
      integer i,j,neqn
      double precision m1,m2,m3,m4,m5,m6,m7,xa,ya,xb,yb,xc,yc,c0,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,l0,
     +     ss,sa,sb,sc,sd,ta,tb,u,ua,ub,zf,zt,fa,mom,
     +     sibe,sith,siga,siph,side,siom,siep,
     +     cobe,coth,coga,coph,code,coom,coep,
     +     sibeth,siphde,siomep,cobeth,cophde,coomep,
     +     bep,thp,php,dep,omp,epp,
     +     m(7,7),ff(7),gp(6,7),g(6),xd,yd,lang,force,fx,fy
      parameter (m1=.04325d0,m2=.00365d0,m3=.02373d0,m4=.00706d0,
     +     m5=.07050d0,m6=.00706d0,m7=.05498d0,
     +     xa=-.06934d0,ya=-.00227d0,
     +     xb=-0.03635d0,yb=.03273d0,
     +     xc=.014d0,yc=.072d0,c0=4530d0,
     +     i1=2.194d-6,i2=4.410d-7,i3=5.255d-6,i4=5.667d-7,
     +     i5=1.169d-5,i6=5.667d-7,i7=1.912d-5,
     +     d=28d-3,da=115d-4,e=2d-2,ea=1421d-5,
     +     rr=7d-3,ra=92d-5,l0=7785d-5,
     +     ss=35d-3,sa=1874d-5,sb=1043d-5,sc=18d-3,sd=2d-2,
     +     ta=2308d-5,tb=916d-5,u=4d-2,ua=1228d-5,ub=449d-5,
     +     zf=2d-2,zt=4d-2,fa=1421d-5,mom=33d-3)
      
      sibe = dsin(y(1))
      sith = dsin(y(2))
      siga = dsin(y(3))
      siph = dsin(y(4))
      side = dsin(y(5))
      siom = dsin(y(6))
      siep = dsin(y(7))
c     
      cobe = dcos(y(1))
      coth = dcos(y(2))
      coga = dcos(y(3))
      coph = dcos(y(4))
      code = dcos(y(5))
      coom = dcos(y(6))
      coep = dcos(y(7))
c     
      sibeth = dsin(y(1)+y(2))
      siphde = dsin(y(4)+y(5))
      siomep = dsin(y(6)+y(7))
c     
      cobeth = dcos(y(1)+y(2))
      cophde = dcos(y(4)+y(5))
      coomep = dcos(y(6)+y(7))
c
      bep = y(8)
      thp = y(9)
      php = y(11)
      dep = y(12)
      omp = y(13)
      epp = y(14)
c
      do 20 j = 1,7
         do 10 i = 1,7
            m(i,j) = 0d0
 10      continue
 20   continue
c
      m(1,1) = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 + i2
      m(2,1) = m2*(da**2-da*rr*coth) + i2
      m(2,2) = m2*da**2 + i2
      m(3,3) = m3*(sa**2+sb**2) + i3
      m(4,4) = m4*(e-ea)**2 + i4
      m(5,4) = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
      m(5,5) = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)
     &     + i4 + i5
      m(6,6) = m6*(zf-fa)**2 + i6
      m(7,6) = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
      m(7,7) = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)
     &     + i6 + i7
      
      do 40 j=2,7
         do 30 i=1,j-1
            m(i,j) = m(j,i)
 30      continue
 40   continue
c     
      xd = sd*coga + sc*siga + xb
      yd = sd*siga - sc*coga + yb
      lang  = dsqrt ((xd-xc)**2 + (yd-yc)**2)
      force = - c0 * (lang - l0)/lang
      fx = force * (xd-xc)
      fy = force * (yd-yc)
      ff(1) = mom - m2*da*rr*thp*(thp+2*bep)*sith
      ff(2) = m2*da*rr*bep**2*sith
      ff(3) = fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga)
      ff(4) = m4*zt*(e-ea)*dep**2*coph
      ff(5) = - m4*zt*(e-ea)*php*(php+2*dep)*coph
      ff(6) = - m6*u*(zf-fa)*epp**2*coom
      ff(7) = m6*u*(zf-fa)*omp*(omp+2*epp)*coom
c
      do 60 j=1,7
         do 50 i=1,6
            gp(i,j) = 0d0
 50      continue
 60   continue
c
      gp(1,1) = - rr*sibe + d*sibeth
      gp(1,2) = d*sibeth
      gp(1,3) = - ss*coga
      gp(2,1) = rr*cobe - d*cobeth
      gp(2,2) = - d*cobeth
      gp(2,3) = - ss*siga
      gp(3,1) = - rr*sibe + d*sibeth
      gp(3,2) = d*sibeth
      gp(3,4) = - e*cophde
      gp(3,5) = - e*cophde + zt*side
      gp(4,1) = rr*cobe - d*cobeth
      gp(4,2) = - d*cobeth
      gp(4,4) = - e*siphde
      gp(4,5) = - e*siphde - zt*code
      gp(5,1) = - rr*sibe + d*sibeth
      gp(5,2) = d*sibeth
      gp(5,6) = zf*siomep
      gp(5,7) = zf*siomep - u*coep
      gp(6,1) = rr*cobe - d*cobeth
      gp(6,2) = - d*cobeth
      gp(6,6) = - zf*coomep
      gp(6,7) = - zf*coomep - u*siep
c
      g(1) = rr*cobe - d*cobeth - ss*siga - xb
      g(2) = rr*sibe - d*sibeth + ss*coga - yb
      g(3) = rr*cobe - d*cobeth - e*siphde - zt*code - xa
      g(4) = rr*sibe - d*sibeth + e*cophde - zt*side - ya
      g(5) = rr*cobe - d*cobeth - zf*coomep - u*siep - xa
      g(6) = rr*sibe - d*sibeth - zf*siomep + u*coep - ya
c
      do 70 i=1,14
         dy(i) = y(i+7)
   70 continue

      do 100 i=15,21
         dy(i) = -ff(i-14)
         do 80 j=1,7
            dy(i) = dy(i)+m(i-14,j)*y(j+14)
   80    continue
         do 90 j=1,6
            dy(i) = dy(i)+gp(j,i-14)*y(j+21)
   90    continue
  100 continue
      do 110 i=22,27
         dy(i) = g(i-21)
  110 continue

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE PDERV(T,Y,DFDY,NEQN,YPRIME,MN,CON,IPAR,RPAR,IERR)
      DOUBLE PRECISION T, Y(NEQN), DFDY(MN,NEQN),CON
      INTEGER NEQN  
      
c-----------------------------------------------------------------------
c     the Jacobian computed here is an approximation, see p. 540 of
c     Hairer & Wanner `solving ordinary differential equations II'
c-----------------------------------------------------------------------     
      double precision m1,m2,m3,m4,m5,m6,m7,
     +     i1,i2,i3,i4,i5,i6,i7,d,da,e,ea,rr,ra,
     +     ss,sa,sb,ta,tb,u,ua,ub,zf,zt,fa,
     +     sibe,siga,siph,side,siom,siep,
     +     cobe,coth,coga,code,coep,
     +     sibeth,siphde,siomep,cobeth,cophde,coomep,
     +     m(7,7),gp(6,7)
      parameter (m1=.04325d0,m2=.00365d0,m3=.02373d0,m4=.00706d0,
     +     m5=.07050d0,m6=.00706d0,m7=.05498d0,
     +     i1=2.194d-6,i2=4.410d-7,i3=5.255d-6,i4=5.667d-7,
     +     i5=1.169d-5,i6=5.667d-7,i7=1.912d-5,
     +     d=28d-3,da=115d-4,e=2d-2,ea=1421d-5,
     +     rr=7d-3,ra=92d-5,
     +     ss=35d-3,sa=1874d-5,sb=1043d-5,
     +     ta=2308d-5,tb=916d-5,u=4d-2,ua=1228d-5,ub=449d-5,
     +     zf=2d-2,zt=4d-2,fa=1421d-5)

      sibe = dsin(y(1))
      siga = dsin(y(3))
      siph = dsin(y(4))
      side = dsin(y(5))
      siom = dsin(y(6))
      siep = dsin(y(7))
c
      cobe = dcos(y(1))
      coth = dcos(y(2))
      coga = dcos(y(3))
      code = dcos(y(5))
      coep = dcos(y(7))
c
      sibeth = dsin(y(1)+y(2))
      siphde = dsin(y(4)+y(5))
      siomep = dsin(y(6)+y(7))
c
      cobeth = dcos(y(1)+y(2))
      cophde = dcos(y(4)+y(5))
      coomep = dcos(y(6)+y(7))
c
      do 51 j = 1,7
         do 52 i = 1,7
            m(i,j) = 0d0
 52      continue
 51      continue
c
      m(1,1) = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 +i2
      m(2,1) = m2*(da**2-da*rr*coth) + i2
      m(2,2) = m2*da**2 + i2
      m(3,3) = m3*(sa**2+sb**2) + i3
      m(4,4) = m4*(e-ea)**2 + i4
      m(5,4) = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
      m(5,5) = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)
     +         + i4 + i5
      m(6,6) = m6*(zf-fa)**2 + i6
      m(7,6) = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
      m(7,7) = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)
     +         + i6 + i7

      do 40 j=2,7
         do 30 i=1,j-1
            m(i,j) = m(j,i)
   30    continue
   40 continue
c
      do 60 j=1,7
         do 50 i=1,6
            gp(i,j) = 0d0
   50    continue
   60 continue
c
      gp(1,1) = - rr*sibe + d*sibeth
      gp(1,2) = d*sibeth
      gp(1,3) = - ss*coga
      gp(2,1) = rr*cobe - d*cobeth
      gp(2,2) = - d*cobeth
      gp(2,3) = - ss*siga
      gp(3,1) = - rr*sibe + d*sibeth
      gp(3,2) = d*sibeth
      gp(3,4) = - e*cophde
      gp(3,5) = - e*cophde + zt*side
      gp(4,1) = rr*cobe - d*cobeth
      gp(4,2) = - d*cobeth
      gp(4,4) = - e*siphde
      gp(4,5) = - e*siphde - zt*code
      gp(5,1) = - rr*sibe + d*sibeth
      gp(5,2) = d*sibeth
      gp(5,6) = zf*siomep
      gp(5,7) = zf*siomep - u*coep
      gp(6,1) = rr*cobe - d*cobeth
      gp(6,2) = - d*cobeth
      gp(6,6) = - zf*coomep
      gp(6,7) = - zf*coomep - u*siep
c    
      do 80 j=1,neqn
         do 70 i=1,neqn
            dfdy(i,j) = 0d0
 70      continue
 80   continue
      do 90 i=1,14
         dfdy(i,i+7) = 1d0
 90   continue
      do 110 i=1,7
         do 100 j=1,7
            dfdy(14+j,14+i) = m(j,i)
 100     continue
 110  continue
      do 130 i=1,6
         do 120 j=1,7
            dfdy(14+j,21+i) = gp(i,j)
 120     continue
 130  continue
      do 150 i=1,7
         do 140 j=1,6
            dfdy(21+j,i) = gp(j,i)
 140     continue
 150  continue
C
      do i=1,neqn
         do  j=1,neqn
            dfdy(i,j) = -dfdy(i,j)     
         enddo
      enddo
c compute pd = -df/dy + con*M
      do j=1,14              
         dfdy(j,j) = 1.0d0/con+dfdy(j,j)                        
      enddo
c      
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE INITIAL(NEQN,T,Y, YPRIME)
      DOUBLE PRECISION Y(NEQN),T, YPRIME(NEQN)
C
      Y(1 ) = -0.0617138900142764496358948458001D0
      Y(2 ) =  0D0
      Y(3 ) =  0.455279819163070380255912382449D0
      Y(4 ) =  0.222668390165885884674473185609D0
      Y(5 ) =  0.487364979543842550225598953530D0
      Y(6 ) = -0.222668390165885884674473185609D0
      Y(7 ) =  1.23054744454982119249735015568D0
      Y(8 ) =  0D0
      Y(9 ) =  0D0
      Y(10) =  0D0
      Y(11) =  0D0
      Y(12) =  0D0
      Y(13) =  0D0
      Y(14) =  0D0
      Y(15) =  14222.4439199541138705911625887D0
      Y(16) = -10666.8329399655854029433719415D0
      Y(17) =  0D0
      Y(18) =  0D0
      Y(19) =  0D0
      Y(20) =  0D0
      Y(21) =  0D0
      Y(22) =  98.5668703962410896057654982170D0
      Y(23) = -6.12268834425566265503114393122D0
      Y(24) =  0D0
      Y(25) =  0D0
      Y(26) =  0D0
      Y(27) =  0D0
C
      do 10 k=1,14
         yprime(k) = y(k+7)
   10 continue
      do 20 k=15,27
         yprime(k) = 0d0
   20 continue
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SOLN(NEQN,T,Y)
      DOUBLE PRECISION T,Y(NEQN)
C
      Y(  1) =  0.1581077119629904D+2
      Y(  2) = -0.1575637105984298D+2
      Y(  3) =  0.4082224013073101D-1
      Y(  4) = -0.5347301163226948D+0
      Y(  5) =  0.5244099658805304D+0
      Y(  6) =  0.5347301163226948D+0
      Y(  7) =  0.1048080741042263D+1
      Y(  8) =  0.1139920302151208D+4
      Y(  9) = -0.1424379294994111D+4
      Y( 10) =  0.1103291221937134D+2
      Y( 11) =  0.1929337464421385D+2
      Y( 12) =  0.5735699284790808D+0
      Y( 13) = -0.1929337464421385D+2
      Y( 14) =  0.3231791658026955D+0
      Y( 15) = -0.2463176316945196D+5
      Y( 16) =  0.5185037701610329D+5
      Y( 17) =  0.3241025686413781D+6
      Y( 18) =  0.5667493645176213D+6
      Y( 19) =  0.1674362929479361D+5
      Y( 20) = -0.5667493645176222D+6
      Y( 21) =  0.9826520791458422D+4
      Y( 22) =  0.1991753333731910D+3
      Y( 23) = -0.2975531228015052D+2
      Y( 24) =  0.2306654119098399D+2
      Y( 25) =  0.3145271365475927D+2
      Y( 26) =  0.2264249232082739D+2
      Y( 27) =  0.1161740700019673D+2
C
      RETURN
      END
C--------------------------------------------------------------------
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
         SKIPI =  I .GT. 7

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
      

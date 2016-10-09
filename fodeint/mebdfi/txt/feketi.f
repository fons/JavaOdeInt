C*************************************************************************
C     02/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 2 FEKETE PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...        
      PARAMETER (ND=160, LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),MBND(4),
     +     ERROR(ND),RPAR(2),IPAR(2), DY(ND),YPRIME(ND)          
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'FEKETE',
     +     IND=2)  
      Parameter(nart=20, neqn=8*nart)     
C...
      external PDERV, RESID
C
      open(6,file='feketi.out')
      rewind(6)
      open(8,file='feketaci.out')     
      rewind(8)
C
      write(6,2050)IND,PROBLEM,SOLVER
 2050 FORMAT(1X,'COMPUTATIONAL STATISTICS OF THE INDEX',I2,1X,A,
     &     ' PROBLEM USING ', A/)
      WRITE(6,2051)ND 
 2051 FORMAT(1X,'NUMBER OF EQUATIONS :',I4)
      WRITE(6,*)' '      
C...      
C... LOOP FOR DIFFERENT TOLERANCES
C...    
      NTOLMN=2
      NTOLMX=5
      NTOLDF=2
      NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
      TOLST=0.1D0**NTOLMN
      TOLFC=0.1D0**(1.D0/NTOLDF)       
      DO 30 NTOL=1,NRLOOP         
C...  DIMENSION OF THE SYSTEM
         N = 160
C...  ENDPOINT OF INTEGRATION
         X=0.0D0                 
         XEND=1.0d+3
C...  REQUIRED TOLERANCE
         rtol=tolst
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
         iwork(14)=200000            
         H=RTOL  
         XOUT=XEND
         IWORK(1)=6*nart
         IWORK(2)=2*nart
         IWORK(3)=0         
         MBND(1) = NEQN
         MBND(2) = NEQN
         MBND(3) = NEQN
         MBND(4) = NEQN
         it1=mclock()         
         
 220     CONTINUE
C...  CALL OF THE SUBROUTINE           
         CALL MEBDFI(NEQN,X,H,Y,YPRIME,XOUT,XEND,MF,INDEX,LOUT,LWORK, 
     +        WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,ATOL,
     +        RPAR,IPAR,PDERV,RESID,IERR)
c         WRITE(*,*) 'INDEX =',INDEX, 'IERR =',IERR
         IF (INDEX .EQ. 1) THEN
            INDEX = 0 
            GOTO 220
         ELSEIF (INDEX .NE. 0) THEN
            write(*,*) 'MEBDFI return index = ',index
C            stop
C...  GO TO THE NEXT TOLERANCE
            GOTO 25
         ENDIF

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
      CLOSE(6)
      CLOSE(8)
      STOP
      END
c----------------------------------------------------------------------
      SUBROUTINE RESID(NEQN,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),DELTA(NEQN),YPRIME(NEQN),IPAR(2),RPAR(2)
C
      CALL F(NEQN,T,Y,DELTA,IPAR,RPAR,IERR)
      NART = NEQN/8
      DO I=1,6*NART
         DELTA(I) = YPRIME(I) -DELTA(I)
      ENDDO
      DO I=6*NART+1, 8*NART
         DELTA(I) = DELTA(I)
      ENDDO

      RETURN
      END
C----------------------------------------------------------------------
      subroutine f(neqn,t,y,dy,ipar,rpar,ierr)
      integer i,j,k,nart,maxn,IPAR(*)
      parameter(maxn=150)
      double precision t,y(neqn),dy(neqn),RPAR(*)
      double precision p(maxn,3),q(maxn,3),lam(maxn),mu(maxn),
     +     pp(maxn,3),qp(maxn,3),phi(maxn),gpq(maxn),
     +     fa(maxn,maxn,3),rn,alpha
      
      alpha=.5d0
      nart=neqn/8
      
      do 20 i=1,nart
         do 10 k=1,3
            p(i,k)=y(3*(i-1)+k)
            q(i,k)=y(3*nart+3*(i-1)+k)
 10      continue
         lam(i)=y(6*nart+i)
         mu(i)=y(7*nart+i)
 20   continue
      do 70 i=1,nart
         do 60 j=1,nart
            if(i.eq.j)then
               do 30 k=1,3
                  fa(i,j,k)=0d0
 30            continue
            else
               rn=0d0
               do 40 k=1,3
                  rn=rn+(p(i,k)-p(j,k))**2
 40            continue
               do 50 k=1,3
                  fa(i,j,k)=(p(i,k)-p(j,k))/rn
 50            continue
            endif
 60      continue
 70   continue
      do 100 i=1,nart
         do 90 k=1,3
            pp(i,k)=q(i,k)+2*mu(i)*p(i,k)
            qp(i,k)=-alpha*q(i,k)+2*lam(i)*p(i,k)
            do 80 j=1,nart
               qp(i,k)=qp(i,k)+fa(i,j,k)
 80         continue
 90      continue
 100  continue
      do 120 i=1,nart
         phi(i)=-1d0
         gpq(i)=0d0
         do 110 k=1,3
            phi(i)=phi(i)+p(i,k)**2
            gpq(i)=gpq(i)+2*p(i,k)*q(i,k)
 110     continue
 120  continue
      do 140 i=1,nart
         do 130 k=1,3
            dy(3*(i-1)+k)=pp(i,k)
            dy(3*nart+3*(i-1)+k)=qp(i,k)
 130     continue
         dy(6*nart+i)=phi(i)
         dy(7*nart+i)=gpq(i)
 140  continue
      
      return
      end
c-----------------------------------------------------------------------
      subroutine pderv(T,Y,DFY,NEQN,YPRIME,NN,IPAR,RPAR,IERR,CON)      
      double precision t,y(neqn),dfy(neqn,neqn)      
      integer i,j,k,l,m,nart,maxn,neqn
      parameter(maxn=150)
      double precision p(maxn,3),q(maxn,3),lam(maxn),mu(maxn),
     +     rn(maxn,maxn),alpha
      
      alpha=.5d0
      nart=neqn/8
      
      do 20 i=1,nart
         do 10 k=1,3
            p(i,k)=y(3*(i-1)+k)
            q(i,k)=y(3*nart+3*(i-1)+k)
 10      continue
         lam(i)=y(6*nart+i)
         mu(i)=y(7*nart+i)
 20   continue
      
      do 50 j=1,nart
         do 40 i=1,nart
            rn(i,j)=0d0
            do 30 k=1,3
               rn(i,j)=rn(i,j)+(p(i,k)-p(j,k))**2
 30         continue
 40      continue
 50   continue
           
cJ_pp
      do 90 i=1,nart
         do 80 k=1,3
            dfy(3*(i-1)+k,3*(i-1)+k)=2d0*mu(i)
 80      continue
 90   continue
cJ_pq
      do 110 i=1,nart
         do 100 k=1,3
            dfy(3*(i-1)+k,3*nart+3*(i-1)+k)=1d0
 100     continue
 110  continue
c J_pm
      do 130 i=1,nart
         do 120 k=1,3
            dfy(3*(i-1)+k,7*nart+i)=2d0*p(i,k)
 120     continue
 130  continue
      
c J_qp l=i,m=k
      do 160 i=1,nart
         do 150 k=1,3
            dfy(3*nart+3*(i-1)+k,3*(i-1)+k)=2d0*lam(i)
            do 140 j=1,nart
               if(j.ne.i)then
                  dfy(3*nart+3*(i-1)+k,3*(i-1)+k)=
     +                 dfy(3*nart+3*(i-1)+k,3*(i-1)+k)+
     +                 (rn(i,j)-2d0*(p(i,k)-p(j,k))**2)/rn(i,j)**2
               endif
 140        continue
 150     continue
 160  continue
c J_qp l=i,m<>k
      do 200 i=1,nart
         do 190 k=1,3
            do 180 m=1,3
               if(m.ne.k)then
                  do 170 j=1,nart
                     if(j.ne.i)then
                        dfy(3*nart+3*(i-1)+k,3*(i-1)+m)=
     +                       dfy(3*nart+3*(i-1)+k,3*(i-1)+m)-
     +                       2d0*(p(i,k)-p(j,k))*(p(i,m)-p(j,m))/
     +                       rn(i,j)**2
                     endif
 170              continue
               endif
 180        continue
 190     continue
 200  continue
c J_qp l<>i,m=k
      do 230 i=1,nart
         do 220 l=1,nart
            if(l.ne.i)then
               do 210 k=1,3
                  dfy(3*nart+3*(i-1)+k,3*(l-1)+k)=
     +                 (-rn(i,l)+2d0*(p(i,k)-p(l,k))**2)/rn(i,l)**2
 210           continue
            endif
 220     continue
 230  continue
c J_qp l<>i,m<>k
      do 270 i=1,nart
         do 260 l=1,nart
            if(l.ne.i)then
               do 250 k=1,3
                  do 240 m=1,3
                     if(m.ne.k)then
                        dfy(3*nart+3*(i-1)+k,3*(l-1)+m)=
     +                       2d0*(p(i,k)-p(l,k))*(p(i,m)-p(l,m))/
     +                       rn(i,l)**2
                     endif
 240              continue
 250           continue
            endif
 260     continue
 270  continue

c J_qq
      do 290 i=1,nart
         do 280 k=1,3
            dfy(3*nart+3*(i-1)+k,3*nart+3*(i-1)+k)=-alpha
 280     continue
 290  continue
c J_ql
      do 310 i=1,nart
         do 300 k=1,3
            dfy(3*nart+3*(i-1)+k,6*nart+i)=2d0*p(i,k)
 300     continue
 310  continue

c J_lp
      do 330 i=1,nart
         do 320 k=1,3
            dfy(6*nart+i,3*(i-1)+k)=2d0*p(i,k)
 320     continue
 330  continue

c J_mp
      do 350 i=1,nart
         do 340 k=1,3
            dfy(7*nart+i,3*(i-1)+k)=2d0*q(i,k)
 340     continue
 350  continue
c J_mq
      do 370 i=1,nart
         do 360 k=1,3
            dfy(7*nart+i,3*nart+3*(i-1)+k)=2d0*p(i,k)
 360     continue
 370  continue
c compute pd = -df/dy + con*M
      do i=1,neqn
         do j=1,neqn
            dfy(i,j) =-dfy(i,j)            
         enddo
      enddo
      do i=1,6*nart
         dfy(i,i) = 1.0d0/con+dfy(i,i)
      enddo
      
      return
      end
C-----------------------------------------------------------------------
      SUBROUTINE INIT(N,T,Y,YPRIME)
      DOUBLE PRECISION Y(N),T,PI,RPAR(2),DY(160),YPRIME(N)
      INTEGER IPAR(2), N, NART
      double precision alpha,beta      
C     
      pi=3.141592653589793238462643383d0
      NART = N/8
      do 10 i=1,3
         alpha=2*pi*dble(i)/dble(3)+pi/dble(13)
         beta=3*pi/dble(8)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
 10   continue
      do 20 i=4,10
         alpha=2*pi*dble(i-3)/dble(7)+pi/dble(29)
         beta=pi/dble(8)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
 20   continue
      do 30 i=11,16
         alpha=2*pi*dble(i-10)/dble(6)+pi/dble(7)
         beta=-2*pi/dble(15)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
 30   continue
      do 40 i=17,20
         alpha=2*pi*dble(i-17)/dble(4)+pi/dble(17)
         beta=-3*pi/dble(10)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
 40   continue
      
      do 50 i=3*nart+1,6*nart
         y(i)=0d0
 50   continue
      do 60 i=6*nart+1,8*nart
         y(i)=0d0
 60   continue
      call f(n,0d0,y,yprime,ipar,rpar,ierr)
      do 80 i=1,nart
         do 70 j=1,3
            y(6*nart+i)=y(6*nart+i)+y(3*(i-1)+j)*
     +           yprime(3*nart+3*(i-1)+j)
 70      continue
         y(6*nart+i)=-y(6*nart+i)/2d0
 80   continue
      call f(n,0d0,y,yprime,ipar,rpar,ierr)
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE SOLN(NEQN,T,Y)
      DOUBLE PRECISION T, Y(NEQN)      
c
c computed at medusa
c problem           fekete
c solver            RADAU5
c rtol              0.100E-11
c atol              0.100E-11
c h0                0.100E-11
c ymax
c # scd                  0.55
c # steps                1249
c # steps accepted       1249
c # f-eval               9076
c # Jac-eval              589
c # LU-decomp             770
c CPU-time              22.87
c
      y(  1) =  -0.4070263380333202d+00
      y(  2) =   0.3463758772791802d+00
      y(  3) =   0.8451942450030429d+00
      y(  4) =   0.7752934752521549d-01
      y(  5) =  -0.2628662719972299d+00
      y(  6) =   0.9617122871829146d+00
      y(  7) =   0.7100577833343567d+00
      y(  8) =   0.1212948055586120d+00
      y(  9) =   0.6936177005172217d+00
      y( 10) =   0.2348267744557627d+00
      y( 11) =   0.7449277976923311d+00
      y( 12) =   0.6244509285956391d+00
      y( 13) =  -0.4341114738782885d+00
      y( 14) =   0.8785430442262876d+00
      y( 15) =   0.1992720444237660d+00
      y( 16) =  -0.9515059600312596d+00
      y( 17) =   0.2203508762787005d+00
      y( 18) =   0.2146669498274008d+00
      y( 19) =  -0.6385191643609878d+00
      y( 20) =  -0.4310833259390688d+00
      y( 21) =   0.6375425027722121d+00
      y( 22) =  -0.1464175087914336d+00
      y( 23) =  -0.9380871635228862d+00
      y( 24) =   0.3139337298744690d+00
      y( 25) =   0.5666974065069942d+00
      y( 26) =  -0.6739221885076542d+00
      y( 27) =   0.4740073135462156d+00
      y( 28) =   0.9843259538440293d+00
      y( 29) =  -0.1696995357819996d+00
      y( 30) =  -0.4800504290609090d-01
      y( 31) =   0.1464175087914331d+00
      y( 32) =   0.9380871635228875d+00
      y( 33) =  -0.3139337298744656d+00
      y( 34) =  -0.7092757549979014d+00
      y( 35) =   0.5264062637139616d+00
      y( 36) =  -0.4688542938854929d+00
      y( 37) =  -0.8665731819284478d+00
      y( 38) =  -0.4813878059756024d+00
      y( 39) =  -0.1315929352982178d+00
      y( 40) =  -0.2347897778700538d+00
      y( 41) =  -0.8594340408013130d+00
      y( 42) =  -0.4541441287957579d+00
      y( 43) =   0.5530976940074118d+00
      y( 44) =  -0.7674370265615124d+00
      y( 45) =  -0.3242273140037833d+00
      y( 46) =   0.7711050969896927d+00
      y( 47) =   0.6357041816577034d+00
      y( 48) =   0.3573685519777001d-01
      y( 49) =   0.7103951209379591d+00
      y( 50) =   0.2403570431280519d+00
      y( 51) =  -0.6614886725910596d+00
      y( 52) =  -0.3038208738735660d-01
      y( 53) =   0.4501923293640461d+00
      y( 54) =  -0.8924145871442046d+00
      y( 55) =  -0.5772996158107093d+00
      y( 56) =  -0.1766763414971813d+00
      y( 57) =  -0.7971892020969544d+00
      y( 58) =   0.2414481766969039d+00
      y( 59) =  -0.3416456818373135d+00
      y( 60) =  -0.9082846503446250d+00
      y( 61) =   0.2409619682166627d-15
      y( 62) =  -0.1139818460497816d-15
      y( 63) =   0.1627536276556335d-15
      y( 64) =   0.1745651819597609d-15
      y( 65) =  -0.1914278710633076d-15
      y( 66) =  -0.6639600671806291d-16
      y( 67) =   0.1708576733899083d-15
      y( 68) =  -0.2277602521390053d-15
      y( 69) =  -0.1350782790950654d-15
      y( 70) =   0.2411941341109454d-15
      y( 71) =  -0.1438238671800488d-15
      y( 72) =   0.8087033550666644d-16
      y( 73) =   0.1618239105233347d-15
      y( 74) =   0.1837556152070701d-16
      y( 75) =   0.2715177369929503d-15
      y( 76) =   0.7930078658689191d-16
      y( 77) =   0.7482020588342764d-16
      y( 78) =   0.2746974939098084d-15
      y( 79) =   0.8849338913035911d-16
      y( 80) =  -0.5940734725324115d-16
      y( 81) =   0.4845984056889910d-16
      y( 82) =  -0.3728835248155620d-16
      y( 83) =  -0.4600332954062859d-16
      y( 84) =  -0.1548568884846698d-15
      y( 85) =   0.2507541692375411d-16
      y( 86) =  -0.1560155223230823d-15
      y( 87) =  -0.2517946296860555d-15
      y( 88) =  -0.3739779361502470d-16
      y( 89) =  -0.1381663620885020d-15
      y( 90) =  -0.2784051540342329d-15
      y( 91) =   0.6624397102887671d-16
      y( 92) =   0.4226207488883120d-16
      y( 93) =   0.1571821772296610d-15
      y( 94) =  -0.4112243677286995d-16
      y( 95) =   0.1939960344265876d-15
      y( 96) =   0.2800184977692136d-15
      y( 97) =  -0.9189023375328813d-16
      y( 98) =   0.1392943179389155d-15
      y( 99) =   0.9556003995587458d-16
      y(100) =  -0.2234188557495892d-15
      y(101) =   0.1276804778190781d-15
      y(102) =  -0.1261196211463950d-15
      y(103) =  -0.1887754149742397d-15
      y(104) =  -0.2140788698695373d-16
      y(105) =  -0.2713591291421657d-15
      y(106) =   0.1107887633060814d-15
      y(107) =  -0.1318443715631340d-15
      y(108) =  -0.4521275683078691d-16
      y(109) =  -0.1277688851278605d-15
      y(110) =   0.4850914012115388d-16
      y(111) =  -0.1195891666741192d-15
      y(112) =  -0.1569641653843750d-15
      y(113) =   0.1856239009452638d-15
      y(114) =   0.9898466095646496d-16
      y(115) =  -0.2068030800303723d-15
      y(116) =   0.2451470336752085d-15
      y(117) =   0.9542986459336358d-16
      y(118) =  -0.2456074075580993d-15
      y(119) =   0.1532475480661800d-15
      y(120) =  -0.1229326332276474d-15
      y(121) =  -0.4750000000000000d+01
      y(122) =  -0.4750000000000001d+01
      y(123) =  -0.4750000000000000d+01
      y(124) =  -0.4750000000000000d+01
      y(125) =  -0.4750000000000000d+01
      y(126) =  -0.4750000000000000d+01
      y(127) =  -0.4750000000000000d+01
      y(128) =  -0.4750000000000000d+01
      y(129) =  -0.4750000000000000d+01
      y(130) =  -0.4750000000000000d+01
      y(131) =  -0.4750000000000001d+01
      y(132) =  -0.4750000000000001d+01
      y(133) =  -0.4750000000000000d+01
      y(134) =  -0.4750000000000000d+01
      y(135) =  -0.4750000000000000d+01
      y(136) =  -0.4750000000000000d+01
      y(137) =  -0.4749999999999999d+01
      y(138) =  -0.4750000000000000d+01
      y(139) =  -0.4750000000000000d+01
      y(140) =  -0.4750000000000000d+01
      y(141) =  -0.3537526598492654d-19
      y(142) =   0.2338193888161182d-18
      y(143) =  -0.3267771993164953d-18
      y(144) =   0.2915679914072042d-18
      y(145) =   0.1965183195887647d-18
      y(146) =  -0.6224992924096233d-19
      y(147) =  -0.1715878416756298d-18
      y(148) =  -0.2704741705248803d-18
      y(149) =   0.3008700893194513d-18
      y(150) =  -0.2703121624910402d-18
      y(151) =   0.4243755291982164d-18
      y(152) =   0.2862063003949612d-18
      y(153) =   0.1222125408406218d-19
      y(154) =  -0.4958862706817728d-18
      y(155) =  -0.7070673036251212d-18
      y(156) =  -0.4454983024194383d-18
      y(157) =  -0.1125384872521777d-18
      y(158) =   0.1512898724592511d-18
      y(159) =  -0.6163704221424137d-19
      y(160) =   0.6255426995473074d-19
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
c      WRITE(*,100) 
c 100  FORMAT(' NUMERICAL SOLUTION')
c      WRITE(*,*)
c      WRITE(*,101)
c 101  FORMAT(47X,'SCD')
c      WRITE(*,102)
c 102  FORMAT(5X,'SOLUTION COMPONENT',18X,'-------------    IGNOR')
c      WRITE(*,103)
c 103  FORMAT(41X,'ABS        REL ')
       
      NUMA = 0
      NUMR = 0
      AERR12 = 0D0
      RERR12 = 0D0
      AREEMX = 0D0
      RERRMX = 0D0
C            
      DO I=1,NEQN
         SKIPI =  I .GT. 60

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

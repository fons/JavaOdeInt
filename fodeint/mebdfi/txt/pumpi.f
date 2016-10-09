C*************************************************************************
C     16/03/2000
C*************************************************************************
C     DRIVER FOR MEBDFI ON THE INDEX 2  PUMP PROBLEM
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...
C...  PARAMETERS FOR MEBDFI
C...   
      PARAMETER (ND=9,LWORK=(38+3*ND)*ND+3,LIWORK=ND+14)   
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),TRUE(ND),
     +     ERROR(ND),RPAR(2),IPAR(2),ATOL(ND),RTOL(ND),T(0:41),
     +     MASBND(4),MBND(4),YPRIME(ND),jcount(10)
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFI',PROBLEM = 'PUMP',
     +     IND=2)
C...
      EXTERNAL PDERV, RESID
C...
      open(6,file='pumpi.out')
      rewind(6)
      open(8,file='pumpaci.out')     
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
C...  DIMENSION OF THE SYSTEM
         NEQN=9
C...  ENDPOINT OF INTEGRATION
         NDISC = 39
         T(0)=0.0D0         
         T(1)   = 50d-9
         T(2)   = 60d-9
         T(3)   = 110d-9
         do 10 i=4,39
            T(i) = t(mod(i,4)) + dble(int(i/4))*120d-9
 10      continue
         T(40)   = 1200d-9
C...  REQUIRED TOLERANCE
         RTOL(1)=TOLST
         ATOL(1)=RTOL(1)
         DO I=2,9
            ATOL(I) = ATOL(1)
            RTOL(I) = RTOL(1)
         ENDDO
         DO I= 1,5
            ATOL(I) = 1D-6*ATOL(1)
            rtol(i) = rtol(1)
         ENDDO 
         ATOL(9) = 1000
         RTOL(9) = 1000
         ITOL=5
         WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
         WRITE(6,*) 'RTOL =' ,RTOL(1)
         WRITE(6,*) 'ATOL =' ,ATOL(1), 'and ', ATOL(9)
C...  INITIAL VALUES
         CALL INIT(NEQN,T(0),Y,YPRIME)
C...  SET DEFAULT VALUES 
         MF=22
         INDEX=1
         LOUT=6
         MAXDER=7
         H0 = 1.0d-6*TOLST
         write(*,*) 'HSTART =',H0
C... MAXIMAL NUMBER OF STEPS 
         iwork(14)=100000         
         iwork(1)=8
         iwork(2)=1
         iwork(3)=0
         ierr = 0         
         work(1) = 0.d0
         XOUT=T(1)         
         mbnd(1) = neqn
         mbnd(2) = neqn
         mbnd(3) = neqn
         mbnd(4) = neqn
         do k =5,12
            jcount(k) = 0
         enddo
         
         it1=mclock()         
         TT = T(0)   
         do 40 i=0,ndisc            
            h = h0            
            XOUT = T(I+1)
 220        continue         
         
C...  CALL OF THE SUBROUTINE  
         CALL MEBDFI(NEQN,TT,H,Y,YPRIME,XOUT,T(ndisc+1),MF,INDEX,LOUT,
     +        LWORK,WORK,LIWORK,IWORK,MBND,MAXDER,ITOL,RTOL,
     +        ATOL,RPAR,IPAR,PDERV, RESID,IERR)
c         write(*,*) 'index, ierr=',index,ierr
         IF (INDEX .EQ. 1) THEN
            INDEX = 2
            GOTO 220
         ELSEIF (INDEX .NE. 0) THEN
            WRITE(*,*) 'MEBDFI RETURN INDEX = ',INDEX
c            stop
C...  GO TO THE NEXT TOLERANCE
            GOTO 25
         ENDIF         
         index = 1
         TT = XOUT
 
         do k= 5,12
            jcount(k) = jcount(k) + iwork(k)
         enddo
C
 40      continue
         
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
         do k= 5,12
            iwork(k) = jcount(k) 
         enddo
C
         HUSED = WORK(2)
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
     1 ' LAST STEP SIZE                       ',              D20.9,/,
     2 ' LAST ORDER OF THE METHOD             ',                I10,/,
     3 ' MAXIMUM ORDER USED SO FAR            ',                I10,/,
     4 ' TOTAL NUMBER OF STEPS TAKEN          ',                I10,/, 
     5 ' TOTAL NUMBER OF FAILED STEPS         ',                I10,/,
     6 ' NUMBER OF RESIDUAL EVALUATIONS       ',                I10,/,
     7 ' NUMBER OF JACOBIAN EVALUATIONS       ',                I10,/,
     8 ' NUMBER OF FACTORIZATION              ',                I10,/,
     9 ' NUMBER OF BACKSOLVES                 ',                I10,/,
     + ' NUMBER OF CORRECT DIGITS            ',              F13.2,/,   
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
      SUBROUTINE RESID(NEQN,T,Y,DELTA,YPRIME,IPAR,RPAR,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),DELTA(NEQN),YPRIME(NEQN),IPAR(2),RPAR(2)
C
      CALL F(NEQN,T,Y,DELTA,IPAR,RPAR,IERR)
            
      DELTA(1) = YPRIME(1)             -DELTA(1)
      DELTA(2) = YPRIME(2) + YPRIME(3) -DELTA(2)
      DELTA(3) = YPRIME(4) + YPRIME(5) -DELTA(3)
      DO I=4,NEQN
         DELTA(I) = -DELTA(I)
      ENDDO
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE F(NEQN,T,Y,DF,IPAR,RPAR,IERR)      
      DOUBLE PRECISION T,Y(NEQN),DF(NEQN),RPAR(*)
      integer nnodes,i,j,nj,ipar(*)
      external vin
      double precision qgate, qsrc, qdrain, vin
      double precision capd, caps
      parameter (capd = 0.40d-12, caps = 1.60d-12)

      DF(1) = -y(9)
      DF(2) = 0d0
      DF(3) = 0d0
      DF(4) = -y(6) + vin(t)
      DF(5) = y(1) - qgate  (y(6), y(6)-y(7), y(6)-y(8))
      DF(6) = y(2) - caps*y(7)
      DF(7) = y(3) - qsrc (y(6), y(6)-y(7), y(6)-y(8))
      DF(8) = y(4) - capd*y(8)
      DF(9) = y(5) - qdrain(y(6), y(6)-y(7), y(6)-y(8))
      
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
      subroutine soln(neqn,y)
      double precision y(neqn)

c     computed with PSIDE, tol=2d-8

      integer i

      y(1) = 0.1262800429876759d-12
      do 10 i=2,8
         y(i) = 0d0
   10 continue
      y(9) = 0.1522565920949845d-03
      return
      end
c-------------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime)
      integer neqn,i
      double precision t,y(neqn),yprime(neqn) 
c
      double precision qgate,qsrc,qdrain

      y(1) = qgate  (0d0,0d0,0d0)
      y(2) = 0d0
      y(3) = qsrc   (0d0,0d0,0d0)
      y(4) = 0d0
      y(5) = qdrain (0d0,0d0,0d0)
      y(6) = 0d0
      y(7) = 0d0
      y(8) = 0d0
      y(9) = 0d0

      yprime(1) = 0d0
      yprime(2) = 0d0
      yprime(3) = 0d0
      yprime(4) = 0d0
      yprime(5) = 0d0
      yprime(6) = 0d0
      yprime(7) = 0d0
      yprime(8) = 0d0
      yprime(9) = 0d0
      return
      end
c-----------------------------------------------------------------------
c auxiliary functions required by subroutine feval
c-----------------------------------------------------------------------
C     Copyright by Georg Denk, 1995
C
C     Reference: Michael G"unther, Georg Denk, Uwe Feldmann:
C     How models for MOS transistors reflect charge distribution
C     effects.
C
C     Submitted to the journal Mathematical Modelling of Systems,
C     also: Preprint 1745, May 1995, Technische Hochschule Darmstadt,
C     Fachbereich Mathematik.
C
C------------------ Begincharge.f
      double precision function qgate (vgb,vgs,vgd)
C
C ... computes gate charge
C
      double precision vgb, vgs, vgd
C
      double precision vt0, gamma, phi, cox
      parameter (vt0 = 0.20d+0, gamma = 0.350d-1, phi = 1.010d+0,
     ,           cox = 4.0d-12)
C
      double precision ugs, ugd, ugb, ubs, vfb, vte, ugst, ugdt
C
      if ((vgs - vgd) .le. 0.d0) then
         ugs = vgd
         ugd = vgs
      else
         ugs = vgs
         ugd = vgd
      endif
C
      ugb = vgb
      ubs = ugs - ugb
C
      vfb = vt0 - gamma * dsqrt(phi) - phi
      vte = vt0 + gamma *( dsqrt(phi - ubs) - dsqrt(phi) )
C
      if ( ugb .le. vfb) then
         qgate = cox * (ugb - vfb)
      else if ( ugb .gt. vfb .and. ugs .le. vte ) then
         qgate = cox * gamma *
     *         ( dsqrt((gamma/2.d0)**2.d0 + ugb - vfb) -
     *           gamma/2.d0 )
      else if ( ugb .gt .vfb .and. ugs .gt. vte ) then
         ugst = ugs - vte
         if (ugd.gt.vte) then
            ugdt = ugd - vte
         else
            ugdt = 0.d0
         endif
         qgate = cox * ( (2.d0/3.d0) * (ugdt + ugst -
     *           ((ugdt * ugst)/(ugdt+ugst)) ) +
     *           gamma * dsqrt(phi-ubs) )
      endif
C
      return
      end
C
      double precision function qbulk (vgb,vgs,vgd)
C
C ... computes bulk charge
C
      double precision vgb, vgs, vgd
C
      double precision vt0, gamma, phi, cox
      parameter (vt0 = 0.20d+0, gamma = 0.350d-1, phi = 1.010d+0,
     ,           cox = 4.0d-12)
C
      double precision ugs, ugd, ugb, ubs, vfb, vte
C
      if ((vgs - vgd) .le. 0.d0) then
         ugs = vgd
         ugd = vgs
      else
         ugs = vgs
         ugd = vgd
      endif
C
      ugb = vgb
      ubs = ugs - ugb
C
      vfb = vt0 - gamma * dsqrt(phi) - phi
      vte = vt0 + gamma *( dsqrt(phi - ubs) - dsqrt(phi) )
C
      if (ugb.le.vfb) then
         qbulk = - cox * (ugb - vfb)
      else if ((ugb.gt.vfb).and.(ugs.le.vte)) then
         qbulk = - cox * gamma *
     *         ( dsqrt((gamma/2.d0)**2.d0 + ugb - vfb) -
     *         gamma/2.d0 )
      else if ( ugb .gt .vfb .and. ugs .gt. vte ) then
         qbulk = - cox * gamma * dsqrt(phi-ubs)
      endif
C
      return
      end
C
      double precision function qsrc (vgb,vgs,vgd)
C
C ... computes gate charge
C
      double precision vgb, vgs, vgd
C
      double precision vt0, gamma, phi, cox
      parameter (vt0 = 0.20d+0, gamma = 0.350d-1, phi = 1.010d+0,
     ,           cox = 4.0d-12)
C
      double precision ugs, ugd, ugb, ubs, vfb, vte, ugst, ugdt
C
      if ((vgs - vgd) .le. 0.d0) then
         ugs = vgd
         ugd = vgs
      else
         ugs = vgs
         ugd = vgd
      endif
C
      ugb = vgb
      ubs = ugs - ugb
C
      vfb = vt0 - gamma * dsqrt(phi) - phi
      vte = vt0 + gamma *( dsqrt(phi - ubs) - dsqrt(phi) )
C
      if (ugb.le.vfb) then
         qsrc = 0.d0
      else if ((ugb.gt.vfb).and.(ugs.le.vte)) then
         qsrc = 0.d0
      else if ( ugb .gt .vfb .and. ugs .gt. vte ) then
         ugst = ugs - vte
         if (ugd.ge.vte) then
            ugdt = ugd - vte
         else
            ugdt = 0.d0
         endif
         qsrc = - cox * (1.d0/3.d0) * (ugdt + ugst -
     *           ((ugdt * ugst)/(ugdt+ugst)) )
      endif
C
      return
      end
C
      double precision function qdrain (vgb,vgs,vgd)
C
C ... computes drain charge
C
      double precision vgb, vgs, vgd
C
      double precision vt0, gamma, phi, cox
      parameter (vt0 = 0.20d+0, gamma = 0.350d-1, phi = 1.010d+0,
     ,           cox = 4.0d-12)
C
      double precision ugs, ugd, ugb, ubs, vfb, vte, ugst, ugdt
C
      if ((vgs - vgd) .le. 0.d0) then
         ugs = vgd
         ugd = vgs
      else
         ugs = vgs
         ugd = vgd
      endif
C
      ugb = vgb
      ubs = ugs - ugb
C
      vfb = vt0 - gamma * dsqrt(phi) - phi
      vte = vt0 + gamma *( dsqrt(phi - ubs) - dsqrt(phi) )
C
      if (ugb.le.vfb) then
         qdrain = 0.d0
      else if ((ugb.gt.vfb).and.(ugs.le.vte)) then
         qdrain = 0.d0
      else if ( ugb .gt .vfb .and. ugs .gt. vte ) then
         ugst = ugs - vte
         if (ugd.ge.vte) then
            ugdt = ugd - vte
         else
            ugdt = 0.d0
         endif
         qdrain = - cox * (1.d0/3.d0) * (ugdt + ugst -
     *           ((ugdt * ugst)/(ugdt+ugst)) )
      endif
C
      return
      end
C
      double precision function vin(t)
C
C ... computes pulsed input voltage
C
      double precision t
C
      double precision vhigh, deltat, t1, t2, t3, dummy
      parameter (vhigh = 20.0d+0)
C
      deltat = 120.0d-9
      t1 = 50.0d-9
      t2 = 60.0d-9
      t3 = 110.0d-9
C
      dummy = dmod(t, deltat)
      if (dummy .lt. t1) then
C ... 0 <= t' < 50ns
         vin = 0.0d+0
      else
         if (dummy .lt. t2) then
C ... 50ns <= t' < 60ns
            vin = (dummy - t1) * 0.10d+9 * vhigh
         else
            if (dummy .lt. t3) then
C ... 60ns <= t' < 110ns
               vin = vhigh
            else
C ... 110ns <= t' < 120ns
               vin = (deltat - dummy) * 0.10d+9 * vhigh
            endif
         endif
      endif
C
      return
      end
C---------------------------------------------------------------------------
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
         SKIPI =  I .EQ. 9

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

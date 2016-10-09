C*************************************************************************
C     02/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFDAE ON THE VANDERPOL PROBLEM
C************************************************************************* 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ND=2,LWORK=(33+3*ND)*ND+3,LIWORK=ND+14)
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),MBND(4),TRUE(ND),
     +     ERROR(ND),MASBND(4)
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFDAE',PROBLEM = 'VANDERPOL')
C     
      EXTERNAL F,PDERV,MAS
C       
      open(6,file='vdp.res')
      rewind(6)
C
      write(6,2050)IND,PROBLEM,SOLVER
 2050 FORMAT(1X,'COMPUTATIONAL STATISTICS OF THE INDEX',I2,1X,A,
     &     ' PROBLEM USING ', A/)
      WRITE(6,2051)ND 
 2051 FORMAT(1X,'NUMBER OF EQUATIONS :',I2)
      WRITE(6,*)' '
C...     
C...  LOOP FOR DIFFERENT TOLERANCES  
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
         N=2
C... ENDPOINT OF INTEGRATION
         X=0.0D0
         XOUT=1.0D0
         XEND=11.0D0
C...  REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=2
        WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
        WRITE(6,*) 'RTOL =' ,RTOL
        WRITE(6,*) 'ATOL =' ,ATOL
        WRITE(6,*)
C...  INITIAL VALUES        
        Y(1)=2.0D0
        Y(2)=0.0D0
C...  SET DEFAULT VALUES 
        MF=21 
        INDEX=1
        LOUT=6
        MAXDER=7
C...  MAXIMAL NUMBER OF STEPS
        iwork(14)=200000                
        H=1.0D-6
        MASBND(1)=0
        it1=mclock()
C...
C... CALL OF THE SUBROUTINE  
C...        
        DO 20 I=1,11
 220       CONTINUE        
           CALL MEBDF(N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK, 
     +          WORK,LIWORK,IWORK,MBND,masbnd,MAXDER,ITOL,RTOL,ATOL,
     +          Rpar,ipar,f,pderv,mas,ierr)
           IF (INDEX.EQ.1) THEN
              INDEX=0
              GOTO 220
           ELSE
              XOUT=XOUT+1.D0
           END IF
C     
           it2=mclock()
           TIME=(it2-it1)/100.0D+0 
C...  
C...  CALCULATE AND PRINT THE ERROR AT THE END POINT
C...
           CALL SOLN(N,T,TRUE,I)
           sum=0.0D+0
           do 270 k=1,n
              error(k)=dabs(true(k)-y(k))
              sum=sum+((error(k)/(rtol*dabs(y(k))+atol))**2)
 270       continue
           sum=dsqrt(sum/dble(n))
           sumh=sumh+sum
           summ=max(summ,sum)
 20     continue
        summ=summ*atol
        sumh=sumh*atol/(11.0D+0)
        us=-log10(summ)
        uv=-log10(sumh)
C...  
C...  STATISTICS OF THE SIMULATION  
C...
        HUSED  = WORK(2)
        NQUSED = IWORK(4)
        NSTEP  = IWORK(5)
        NFAIL  = IWORK(6)
        NFE    = IWORK(7)
        NJE    = IWORK(8)
        NDEC   = IWORK(9)
        NBSOL  = IWORK(10)
        MAXORD = IWORK(13)                  
C...  
C...  THE CURRENT RUN IS COMPLETE, SO PRINT THE COMPUTATIONAL STAT-
C...  ISTICS FOR MEBDF AND GO ON TO THE NEXT RUN
        WRITE(LOUT,8)HUSED,NQUSED,MAXORD,NSTEP,NFAIL,NFE,NJE,NDEC,
     +       NBSOL,US,TIME
 8      FORMAT(1H /
     1 ' LAST STEP SIZE                       ',              D13.6,/,
     2 ' LAST ORDER OF THE METHOD             ',                I10,/,
     3 ' MAXIMUM ORDER USED SO FAR            ',                I10,/,
     4 ' TOTAL NUMBER OF STEPS TAKEN          ',                I10,/, 
     5 ' TOTAL NUMBER OF FAILED STEPS         ',                I10,/,
     6 ' NUMBER OF FUNCTION EVALUATIONS       ',                I10,/,
     7 ' NUMBER OF JACOBIAN EVALUATIONS       ',                I10,/,
     8 ' NUMBER OF FACTORIZATION              ',                I10,/,
     9 ' NUMBER OF BACKSOLVES                 ',                I10,/,
     + ' NUMBER OF CORRECT DIGITS             ',              F10.3,/,   
     + ' CPU TIME                             ',                F10.4) 
        WRITE(LOUT,*) '---------------------------------------------'         
C...      
C...  NEW TOLERANCE
C...        
 25     TOLST=TOLST*TOLFC
 30   CONTINUE
      close(6)
      STOP
      END
C-------------------------------------------------------------------------
      SUBROUTINE F(N,X,Y,DF,IPAR,RPAR,ierr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),DF(N),IPAR(*),RPAR(*)
C
      EPS=1.D-6
      DF(1)=Y(2)
      PROD=1.D0-Y(1)*Y(1)
      DF(2)=(PROD*Y(2)-Y(1))/EPS
C
      RETURN
      END 
C--------------------------------------------------------------------------
      SUBROUTINE PDERV(X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),DFY(MEBAND,N),IPAR(*),RPAR(*)
C
      EPS=1.D-6
      DFY(1,1)=0.D0
      DFY(1,2)=1.D0
      DFY(2,1)=(-2.0D0*Y(1)*Y(2)-1.0D0)/EPS
      DFY(2,2)=(1.0D0-Y(1)**2)/EPS
C
      RETURN
      END
C-------------------------------------------------------------------------
      subroutine mas(n,am,m,IPAR,RPAR,ierr)
      dimension am(2,2),IPAR(8),RPAR(*)
C
C     DUMMY ROUTINE
C
      return
      end
C-------------------------------------------------------------------------
      SUBROUTINE SOLN(N,T,TRUE,I)
      DOUBLE PRECISION TRUE(N),T
      INTEGER I
C
      IF (I.EQ.1) THEN
         TRUE(1)=-0.1863646254808130E+01
         TRUE(2)= 0.7535430865435460E+00
      ELSEIF (I.EQ.2) THEN
         TRUE(1)=  0.1706167732170456E+01
         TRUE(2)= -0.8928097010248257E+00
      ELSEIF (I.EQ.3) THEN
         TRUE(1)=  -0.1510606936744095E+01
         TRUE(2)=   0.1178380000730945E+01
      ELSEIF (I.EQ.4) THEN
         TRUE(1)=   0.1194414677695905E+01
         TRUE(2)=  -0.2799585996540082E+01
      ELSEIF (I.EQ.5) THEN
         TRUE(1)=   0.1890428596416747E+01
         TRUE(2)=  -0.7345118680166940E+00
      ELSEIF (I.EQ.6) THEN
         TRUE(1)=-0.1737716306805883E+01
         TRUE(2)=   0.8604008653025923E+00
      ELSEIF (I.EQ.7) THEN
         TRUE(1)=   0.1551614645548223E+01
         TRUE(2)=  -0.1102382892533321E+01
      ELSEIF (I.EQ.8) THEN
         TRUE(1)=-0.1278631984330405E+01
         TRUE(2)=   0.2013890883009155E+01
      ELSEIF (I.EQ.9) THEN
         TRUE(1)=  -0.1916552949489830E+01
         TRUE(2)=   0.7169573003463228E+00
      ELSEIF (I.EQ.10) THEN
         TRUE(1)=   0.1768163792391936E+01
         TRUE(2)=-0.8315276407898496E+00
      ELSEIF (I.EQ.11) THEN
         TRUE(1)=-0.1590150544829062E+01
         TRUE(2)= 0.1040279389212485E+01
      ENDIF
C
        RETURN
        END
C-----------------------------------------------------------------------

















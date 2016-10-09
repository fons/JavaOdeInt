C*************************************************************************
C     2/04/2000
C*************************************************************************
C     DRIVER FOR MEBDFDAE ON THE INDEX 0 PENDULUM  PROBLEM
C*************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ND=5,LWORK=(41+2*ND)*ND+3,LIWORK=ND+14)
      DIMENSION Y(ND),WORK(LWORK),IWORK(LIWORK),MBND(4),TRUE(ND),
     +     ERROR(ND),MASBND(4)
      CHARACTER*(*) SOLVER, PROBLEM
      PARAMETER (SOLVER  = ' MEBDFDAE',PROBLEM = 'PENDULUM',IND=0)
C     
      EXTERNAL F,PDERV,MAS
C
      open(6,file='pendul0.res')
      rewind(6)
C
      write(6,2050)IND,PROBLEM,SOLVER
 2050 FORMAT(1X,'COMPUTATIONAL STATISTICS OF THE INDEX',I2,1X,A,
     &     ' PROBLEM USING ', A/)
      WRITE(6,2051)ND 
 2051 FORMAT(1X,'NUMBER OF EQUATIONS :',I2)
      WRITE(6,*)' '
C
C...  LOOP FOR DIFFERENT TOLERANCES   
      NTOLMN=2
      NTOLMX=10
      NTOLDF=2
      NRLOOP=(NTOLMX-NTOLMN)*NTOLDF+1
      TOLST=0.1D0**NTOLMN
      TOLFC=0.1D0**(1.D0/NTOLDF)       
      DO 30 NTOL=1,NRLOOP
         summ=0.0D+0
         sumh=0.0D+0
C...  DIMENSION OF THE SYSTEM
         N=5
C...  ENDPOINT OF INTEGRATION 
         X=0.0D0
         XEND=1.0D0
C...  REQUIRED TOLERANCE
        RTOL=TOLST
        ATOL=RTOL
        ITOL=2
        WRITE(6,*) 'RESULT WITH THE FOLLOWING TOL :'
        WRITE(6,*) 'RTOL =' ,RTOL
        WRITE(6,*) 'ATOL =' ,ATOL
        WRITE(6,*)
C...  INITIAL VALUES
         CALL INIT(N,X,Y)
C...  SET DEFAULT VALUES 
         MF=21 
         INDEX=1
         LOUT=6
         MAXDER=7
C...  MAXIMAL NUMBER OF STEPS
         iwork(14)=100000
         work(1)=0.0D+0
         H=1.0D-6
         XOUT=1.0D0
         iwork(1)=n
         iwork(2)=0
         iwork(3)=0
         masbnd(1)=0
         masbnd(2)=n        
         masbnd(3)=1
         masbnd(4)=n
         MBND(4)=n
C
        it1=mclock()
C...
C... CALL OF THE SUBROUTINE          
 220    CONTINUE        
        CALL MEBDF(N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK, 
     +       WORK,LIWORK,IWORK,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,
     +       rpar,ipar,f,pderv,mas,ierr)
        IF (INDEX.EQ.1) THEN
           INDEX=0
           GOTO 220
        ELSEIF (INDEX .NE. 0) THEN
            WRITE(*,*) 'MEBDFDAE RETURN INDEX = ',INDEX
C            stop
C...  GO TO THE NEXT TOLERANCE
            GOTO 25           
        END IF
C
        it2=mclock()
        TIME=(it2-it1)/100.0D+0 
C...  
C...  CALCULATE AND PRINT THE ERROR AT THE END POINT
C...
        cons1=1.0D+0-y(1)**2-y(2)**2
        cons2=y(1)*y(3)+y(2)*y(4)
        cons3=y(3)**2+y(4)**2-y(5)-y(2)
        write(6,1901)
 1901   format(1x,'the following are conservations laws which hold')
        write(6,*) cons1,cons2,cons3
C
        CALL SOLN(N,X,TRUE)
        sum=0.0D+0        
        do 270 k=1,n
           error(k)=dabs(true(k)-y(k))
           sum=sum+((error(k)/(rtol*dabs(y(k))+atol))**2)
 270    continue
        sum=dsqrt(sum/dble(n))                        
        summ=sum*atol
        us=-log10(summ)
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
        WRITE(LOUT,8)hused,nqused,maxord,NSTEP,NFAIL,NFE,NJE,NDEC,
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
      SUBROUTINE F(N,X,Y,DF,ipar,rpar,ierr)
C...  RIGHT-HAND SIDE OF PENDULUM PROBLEM
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DIMENSION Y(N),DF(N),ipar(*),rpar(*)
C
      Df(1)=Y(3)
      df(2)=y(4)
      df(3)=-y(1)*y(5)
      df(4)=-y(2)*y(5)-1.0D+0
      df(5) = -2.0D+0*(y(1)*y(3)+y(2)*y(4))*y(5)-3.0D+0*y(4)
C
      RETURN
      END 
C-------------------------------------------------------------------------
      SUBROUTINE PDERV(X,Y,DFY,N,MEBAND,ipar,rpar,ierr)
C...  JACOBIAN OF PENDULUM  PROBLEM
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      DIMENSION Y(N),DFY(MEBAND,N),ipar(*),rpar(*)
C
      dfy(1,3)=1.0D+0
      dfy(2,4)=1.0D+0
      dfy(3,1)=-y(5)
      dfy(3,5)=-y(1)
      dfy(4,2)=-y(5)
      dfy(4,5)=-y(2)
      dfy(5,1)=-2.0D+0*y(3)
      dfy(5,2)=-2.0D+0*y(4)
      dfy(5,3)=-2.0D+0*y(1)
      dfy(5,4)=-2.0D+0*y(2)-3.0D+0
C
      RETURN
      END
C-------------------------------------------------------------------------
      subroutine mas(m,am,n,ipar,rpar,ierr)
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      dimension am(5,5),ipar(*),rpar(*)
C
      am(1,1)=1.0D+0
      am(1,2)=0.0D+0
      am(2,1)=0.0D+0
      am(2,2)=1.0D+0
      am(3,3)=1.0D+0
      am(4,4)=1.0D+0
      am(5,5)=1.0D+0
C
      RETURN
      END
C---------------------------------------------------------------------
      SUBROUTINE INIT(N,T,Y)
      DOUBLE PRECISION Y(N),T
      INTEGER N
C
      Y(1)=1.0D+0
      Y(2)=0.0D+0
      Y(3)=0.0D+0
      Y(4)=1.0D+0
      Y(5)=1.0D+0
C
      RETURN
      END
C-------------------------------------------------------------------------
      SUBROUTINE SOLN(N,T,TRUE)
      DOUBLE PRECISION TRUE(N),T
      INTEGER N
C
      TRUE(1)=0.867348640604551746d+0
      TRUE(2)=0.497701050476373918D+0
      TRUE(3)=-0.337480180329866458D-1
      TRUE(4)=0.588130114660715686d-1
      TRUE(5)=-0.493103151528274442D+0
C
      RETURN
      END
C-----------------------------------------------------------------------

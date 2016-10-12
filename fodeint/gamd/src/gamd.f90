!! THE CODE GAMD NUMERICALLY SOLVES A (POSSIBLY STIFF)
!! SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS 
!! OR A LINEARLY IMPLICIT DAE
!!
!! Copyright (C)1997-2007   
!! Authors: FRANCESCA MAZZIA (mazzia@dm.uniba.it)
!!          FELICE IAVERNARO (felix@dm.uniba.it) 
!!
!!
!!This program is free software; you can redistribute it and/or
!!modify it under the terms of the GNU General Public License
!!as published by the Free Software Foundation; either version 2
!!of the License, or (at your option) any later version.
!!
!!This program is distributed in the hope that it will be useful,
!!but WITHOUT ANY WARRANTY; without even the implied warranty of
!!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!GNU General Public License for more details.
!!
!!Licensed under The GNU General Public License, Version 2 or later.
!!    http://www.gnu.org/licenses/info/GPLv2orLater.html
!!
!!You should have received a copy of the GNU General Public License
!!along with this program; if not, write to the Free Software
!!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
!!USA.
!!
!!

SUBROUTINE   GAMD(R,FCN,T0,Y0,TEND,H,            &
       &                  RTOL,ATOL,ITOL,        &
       &                  JAC ,IJAC,MLJAC,MUJAC, &
       &                  MAS ,IMAS,MLMAS,MUMAS, &
       &                  SOLOUT,IOUT,           &
       &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    
    !!
    !!     PURPOSE: THE CODE GAMD NUMERICALLY SOLVES A (POSSIBLY STIFF)
    !!              SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS
    !!              IN THE FORM  Y'=F(T,Y), OR A LINEARLY IMPLICIT DAE 
    !!              MY'=F(T,Y) WITH CONSTANT MASS (FULL OR BANDED) MATRIX M, 
    !!              WITH A GIVEN INITIAL CONDITION. IT IS INTENDED AS AN 
    !!              EXTENSION OF THE CODE GAM TO DAEs (SEE THE REVISION 
    !!              HISTORY FOR DETAILS)
    !!              
    !!     AUTHORS: F. IAVERNARO AND F. MAZZIA
    !!              UNIVERSITA' DEGLI STUDI DI BARI,
    !!              DIPARTIMENTO DI MATEMATICA
    !!              VIA ORABONA 4, 70125 BARI, ITALY
    !!              E-MAIL:  MAZZIA@DM.UNIBA.IT
    !!                       FELIX@DM.UNIBA.IT
    !!
    !!     METHODS: THE METHODS USED ARE IN THE CLASS OF BOUNDARY VALUE
    !!              METHODS (BVMs), NAMELY THE GENERALIZED ADAMS METHODS
    !!              (GAMs) OF ORDER 3-5-7-9 WITH STEP SIZE CONTROL.
    !!
    !!  REFERENCES: L.BRUGNANO, D.TRIGIANTE,  Solving Differential Problems
    !!              by Multistep Initial and Boundary Value Methods,
    !!              Gordon & Breach,1998.
    !!
    !!              F.IAVERNARO, F.MAZZIA,  Block-Boundary Value Methods
    !!              for the solution of Ordinary Differential Equations,
    !!              Siam J. Sci. Comput. 21 (1) (1999) 323--339.
    !!
    !!              F.IAVERNARO, F.MAZZIA,  Solving Ordinary Differential
    !!              Equations by Generalized Adams Methods: properties and
    !!              implementation techniques,
    !!              Appl. Num. Math. 28 (2-4) (1998) 107-126.
    !!
    !! DESCRIPTION: THE CODE GAMD CONSISTS OF TWO FILES:
    !!              - gamd.f90   CONTAINS THE MAIN SUBROUTINES THAT IMPLEMENT THE
    !!                           INTEGRATION PROCEDURE;
    !!
    !!              - gamda.f90  CONTAINS THE ADDITIONAL LINEAR ALGEBRA ROUTINES
    !!                           REQUIRED BY gamd.f90 PLUS SOME OTHER SUBROUTINES 
    !!                           PROPER OF THE USED METHODS;
    !!
    !!    COMMENTS: THE PHILOSOFY AND THE STYLE USED IN WRITING THE CODE ARE VERY
    !!              SIMILAR TO THOSE CHARACTERIZING THE CODE  RADAU5.
    !!              INDEED THE AUTHORS IMPORTED FROM RADAU5 SOME SUBROUTINES,
    !!              COMMENTS AND IMPLEMENTATION TECHNIQUES LEAVING UNCHANGED
    !!              THE NAME AND THE MEANING OF  A NUMBER OF VARIABLES.
    !!              THE AUTHORS ARE VERY GRATEFUL TO ANYONE USES THE CODE AND
    !!              WOULD APPRECIATE ANY CRITICISM AND REMARKS ON HOW IT PERFORMS.
    !!
    !! REVISION HISTORY (YYYY/MM/DD):                   
    !!                      
    !!                   1997/20/08
    !!                      - FIRST VERSION OF GAM
    !!   
    !!                   1999/11/25
    !!                      - CORRECTED OUTPUT LAST STEPSIZE
    !!                      - CORRECTED INPUT  H  
    !!     
    !!                   2003/23/08  
    !!                      - FIRST VERSION OF GAMD
    !!                      - REWRITTEN IN FORTRAN 90
    !!                      - EXTENDED TO LINEARLY IMPLICIT DAES  MY'=F(T,Y)
    !!                        WITH CONSTANT MASS MATRIX M
    !!                      - ADDED IERR IN FCN
    !!                     
    !!                   2006/24/01
    !!                      - CORRECTED THE VALUE CALJAC TO AVOID
    !!                        THE COMPUTATION OF JACOBIAN MORE THEN ONE
    !!                        TIME PER STEPS
    !!                      - CHANGED THE DEFINITION OF TP AND T1 IN ALL
    !!                        THE SUBROUTINE TP(DBLK+1) --> TP(*) 
    !!
    !!                   2006/15/02
    !!                      - CORRECTED ALL THE RUN TIM ERROR FOR
    !!                        SALFORD FTN95 FORTRAN
    !!
    !!                   2007/15/03
    !!                       - CORRECTED INSTRUCTION FOR WORK INDEX 5:7
    !!                   2007/24/05
    !!                       - ADDED THE GNU General Public License
    !! -------------------------------------------------------------------------
    !!     INPUT PARAMETERS
    !! -------------------------------------------------------------------------
    !!     R           DIMENSION OF THE SYSTEM
    !!
    !!     FCN         NAME (EXTERNAL) OF THE SUBROUTINE COMPUTING THE
    !!                 VALUE OF F(T,Y):
    !!                    SUBROUTINE FCN(R,T,Y,F,IERR,RPAR,IPAR)
    !!                    DOUBLE PRECISION X,Y(R),F(R)
    !!                    INTEGER IERR
    !!                    F(1)=...   ETC.
    !!                  IERR = -1 to prevent OVERFLOW
    !!                 (RPAR, IPAR    SEE BELOW)
    !!
    !!     T0          INITIAL T-VALUE
    !!
    !!     Y0          INITIAL VALUES FOR Y
    !!
    !!     TEND        FINAL T-VALUE (TEND-T MUST BE POSITIVE)
    !!
    !!     H           INITIAL STEP SIZE GUESS;
    !!                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,
    !!                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD.
    !!                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS
    !!                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6).
    !!
    !!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY
    !!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH R.
    !!
    !!     ITOL        SWITCH FOR RTOL AND ATOL:
    !!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS.
    !!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF
    !!                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL
    !!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS.
    !!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW
    !!                     RTOL(I)*ABS(Y(I))+ATOL(I).
    !!
    !!     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES
    !!                 THE PARTIAL DERIVATIVES OF F(T,Y) WITH RESPECT TO Y
    !!                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY
    !!                 A DUMMY SUBROUTINE IN THE CASE IJAC=0).
    !!                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM
    !!                    SUBROUTINE JAC(R,T,Y,DFY,LDFY,RPAR,IPAR)
    !!                    DOUBLE PRECISION T,Y(R),DFY(LDFY,R)
    !!                    DFY(1,1)= ...
    !!                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS
    !!                 FURNISHED BY THE CALLING PROGRAM.
    !!                 IF (MLJAC.EQ.R) THE JACOBIAN IS SUPPOSED TO
    !!                    BE FULL AND THE PARTIAL DERIVATIVES ARE
    !!                    STORED IN DFY AS
    !!                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J)
    !!                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND
    !!                    THE PARTIAL DERIVATIVES ARE STORED
    !!                    DIAGONAL-WISE AS
    !!                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J).
    !!
    !!     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN:
    !!                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE
    !!                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED.
    !!                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC.
    !!
    !!     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN:
    !!                    MLJAC=R: JACOBIAN IS A FULL MATRIX. THE LINEAR
    !!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
    !!                       0<=MLJAC<R: MLJAC IS THE LOWER BANDWITH OF JACOBIAN
    !!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW
    !!                       THE MAIN DIAGONAL).
    !!
    !!     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON-
    !!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
    !!                 NEED NOT BE DEFINED IF MLJAC=R.
    !!
    !!     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS-
    !!                 MATRIX M.
    !!                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY
    !!                 MATRIX AND NEEDS NOT TO BE DEFINED;
    !!                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
    !!                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM
    !!                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR)
    !!                    DOUBLE PRECISION AM(LMAS,N)
    !!                    AM(1,1)= ....
    !!                  IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED
    !!                  AS FULL MATRIX LIKE
    !!                         AM(I,J) = M(I,J)
    !!                  ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED
    !!                  DIAGONAL-WISE AS
    !!                         AM(I-J+MUMAS+1,J) = M(I,J).
    !!
    !!     IMAS       GIVES INFORMATION ON THE MASS-MATRIX:
    !!                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY
    !!                       MATRIX, MAS IS NEVER CALLED.
    !!                    IMAS=1: MASS-MATRIX  IS SUPPLIED.
    !!
    !!    MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX:
    !!                       MLMAS=N: THE FULL MATRIX CASE. THE LINEAR ALGEBRA 
    !!                                IS DONE BY FULL-MATRIX GAUSS-ELIMINATION.
    !!                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE MATRIX
    !!                                (>= NUMBER OF NON-ZERO DIAGONALS BELOW
    !!                                THE MAIN DIAGONAL).
    !!                 MLMAS IS SUPPOSED TO BE .LE. MLJAC.
    !!
    !!     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON-
    !!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL).
    !!                 NEED NOT BE DEFINED IF MLMAS=N.
    !!                 MUMAS IS SUPPOSED TO BE .LE. MUJAC.
    !!
    !!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
    !!                 NUMERICAL SOLUTION DURING INTEGRATION.
    !!                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
    !!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0.
    !!                 IT MUST HAVE THE FORM
    !!                    SUBROUTINE SOLOUT(R,TP,YP,FF,NT,DBLK,ORD,RPAR,IPAR,IRTRN)
    !!                    INTEGER R, DBLK, ORD, IPAR(*), IRTRN, NT
    !!                    DOUBLE PRECISION TP(*), YP(R,*), RPAR(*), FF(R,*)
    !!                    ....
    !!                 SOLOUT FURNISHES THE SOLUTION "YP" AT THE
    !!                    GRID-POINTS "TP(*)".
    !!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
    !!                    IS SET <0, GAM  RETURNS TO THE CALLING PROGRAM.
    !!
    !!                 CONTINUOUS OUTPUT:
    !!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
    !!                 FOR THE INTERVAL [TP(1),TP(DBLK+1)] IS AVAILABLE THROUGH
    !!                 THE FUNCTION
    !!                        >>>   CONTR(I,R,T,TP,FF,DBLK,NT)   <<<
    !!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
    !!                 COMPONENT OF THE SOLUTION AT THE POINT T. THE VALUE
    !!                 T SHOULD LIE IN THE INTERVAL [TP(1),TP(DBLK+1)] ON
    !!                 WHICH THE SOLUTION IS COMPUTED AT CURRENT STEP.
    !!                 DO NOT CHANGE THE ENTRIES OF FF(R,*) and NT, IF THE
    !!                 DENSE OUTPUT FUNCTION IS USED.
    !!
    !!     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT:
    !!                    IOUT=0: SUBROUTINE IS NEVER CALLED
    !!                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT.
    !!
    !!     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK".
    !!                 WORK(1), WORK(2),.., WORK(LWORK) SERVE AS PARAMETERS
    !!                 FOR THE CODE. FOR STANDARD USE OF THE CODE
    !!                 WORK(1),..,WORK(LWORK) MUST BE SET TO ZERO BEFORE
    !!                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE.
    !!
    !!     LWORK       DECLARED LENGTH OF ARRAY "WORK". IN THIS VERSION SET
	!!                 LWORK=21
    !!
    !!     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK".
    !!                 IWORK(1),IWORK(2),...,IWORK(LIWORK) SERVE 
    !!                 AS PARAMETERS FOR THE CODE. FOR STANDARD USE, 
    !!                 SET IWORK(1),..,IWORK(9) TO ZERO BEFORE CALLING.
    !!                 - IWORK(10),...,IWORK(24) SERVE AS WORKING AREA.
    !!                 - IWORK(25),IWORK(26),IWORK(27) CONTAIN THE DIMENSION OF
	!!                   THE INDEX 1, INDEX2, INDEX3, VARIABLES RESPECTIVELY.
	!!                   THEY MUST BE PASSED IN GAMD AS INPUT VARIABLES WHEN 
	!!                   SOLVING A DAE; THEY MAY BE SET EQUAL TO ZERO WHEN 
	!!                   SOLVING AN ODE.  
    !!
    !!
    !!     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". IN THIS VERSION SET
	!!                 LIWORK=27
    !!
    !!     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
    !!                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING
    !!                 PROGRAM AND THE FCN, JAC SUBROUTINES.
    !!
    !!     SEE THE EXAMPLE BELOW FOR A PRACTICAL EXPLANATION ON THE USE OF
    !!     SOME OF THE LISTED VARIABLES AND SUBROUTINES.
    !!
    !! -------------------------------------------------------------------------
    !!     EXAMPLE PROBLEM.
    !!---------------------------------------------------------------------------
    !!
    !! the following is a simple example problem, with the coding
    !! needed for its solution by GAMD.  The problem is from chemical
    !! kinetics, and consists of the following three rate equations..
    !!     dy1/dt = -.04*y1 + 1.e4*y2*y3
    !!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
    !!     dy3/dt = 3.e7*y2**2
    !! on the interval from t = 0.0 to t = 4.e10, with initial conditions
    !! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
    !!
    !! the following coding, written in fortran 77, solves this problem with GAMD,
    !! printing results at t = .4, 4., ..., 4.e10.  it uses
    !! itol = 1 and atol much smaller for y2 than for y1 or y3 because
    !! y2 has much smaller values.
    !! At the end of the run, statistical quantities of interest are
    !! printed (see optional outputs in the full description below). 
    !! In case you want to run it, please copy and paste it into a new file, 
	!! remove the first column and save it using the fortran 77 extension .f, 
	!! for example "sample.f". When the source code of the problem is written 
	!! in fortran 77, we suggest to create the object file of the problem separately
	!! and then link it with the code files.   
!        implicit none
!        integer     neq 
!        parameter ( neq=3 )
!        integer     lwork, liwork
!        parameter ( lwork=21, liwork=27)
!
!
!        integer  i, iwork(liwork), itol, iout, nsteps, naccept
!        integer  mljac, mujac, ijac, mlmas, mumas, imas, ipar, idid
!
!        double precision atol(neq), rtol(neq),work(lwork),t,tout,y(neq)
!        double precision h, rpar
!        external feval, jeval, meval, solout
!       
!
!    
!        y(1) = 1.0d0
!        y(2) = 0.0d0
!        y(3) = 0.0d0
!        t = 0.0d0
!        tout = 4.0d10
!        iout = 1
!        ijac = 0
!	  imas = 0 
!        mljac = neq
!        mujac = neq
!        h = 1d-6
!    
!        do i = 1,21
!           work(i) = 0.0d0
!        enddo
!        do i = 1,27
!           iwork(i) = 0.0d0
!        enddo
!    
!        itol = 1
!        rtol(1) = 1.0d-5
!        rtol(2) = 1.0d-5
!        rtol(3) = 1.0d-5
!        atol(1) = 1.0d-5
!        atol(2) = 1.0d-8
!        atol(3) = 1.0d-5
!    
!        rpar = 0.4d0
!
!        call GAMD(neq, feval, t, y, tout, h, rtol, atol, itol,
!     +            jeval, ijac, mljac, mujac,
!     +            meval, imas, mlmas, mumas,
!     +            solout,iout,
!     +            work,lwork,iwork,liwork,rpar,ipar,idid)
!
!
!
!   
!50      format(7h at t =,e12.4,6h   y =,3e14.6)
!        write(6,*)
!        write(6,50) t, y(1), y(2), y(3)
!    
!        nsteps = 0
!        do i=12,23
!           nsteps = nsteps + iwork(i)
!        end do
!        naccept = iwork(12)+iwork(13)+iwork(14)+iwork(15)
!        write(6,41) nsteps,naccept,iwork(10),iwork(11),iwork(24)
!41      format(  ' # steps  =         ',i8,/,
!     +       ' # accept =         ',i8,/,
!     +       ' # f-eval =         ',i8,/,
!     +       ' # Jac-eval =       ',i8,/,
!     +       ' # LU-decomp =      ',i8/)
!    
!        stop
!        end
!    
!        subroutine feval(neq, t, y, ydot, rpar, ipar)
!        double precision t, y(3), ydot(3), rpar
!        integer ipar
!        ydot(1) = -.04d0*y(1) + 1.0d4*y(2)*y(3)
!        ydot(3) = 3.0d7*y(2)*y(2)
!        ydot(2) = -ydot(1) - ydot(3)
!        return
!        end
!
!c-----------------------------------------------------------------------
!
!        subroutine jeval(neqn,t,y,jac,ldim,rpar,ipar)
!        double precision t,y(neqn),jac(ldim,neqn),rpar
!        integer neqn,ldim,ipar
!c
!c       dummy subroutine
!c
!
!        return
!        end
!
!c-----------------------------------------------------------------------
!
!        subroutine meval(ldim,neqn,t,y,yprime,dfddy,ierr,rpar,ipar)
!        integer ldim,neqn,ierr,ipar(*)
!        double precision t,y(neqn),yprime(neqn),dfddy(ldim,neqn),rpar(*)
!c
!c       dummy subroutine
!c
!
!        return
!        end
!
!c-----------------------------------------------------------------------
!
!        subroutine solout(r,tp,yp,ff,nt1,dblk,ord,rpar,ipar,irtrn)
!        use subgamd
!	  implicit none
!	  
!        integer r,dblk,ord,ipar(*),irtrn,nt1
!        double precision tp(*),yp(r,*),rpar(*),ff(r,*),y(3),t
!           t = rpar(1)
!           if ( (tp(1).le.t).and.(t.lt.tp(dblk+1)) ) then
!              y(1) = contr(1,r,t,tp,ff,dblk,nt1)
!              y(2) = contr(2,r,t,tp,ff,dblk,nt1)
!              y(3) = contr(3,r,t,tp,ff,dblk,nt1)
!              write(6,50) t, y(1), y(2), y(3)
!              rpar(1) = rpar(1)*10d0
!           endif
!50      format(7h at t =,e12.4,6h   y =,3e14.6)
!        return
!        end
    !!
    !! -------------------------------------------------------------------------
    !!     SOPHISTICATED SETTING OF PARAMETERS
    !! -------------------------------------------------------------------------
    !!              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK
    !!              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),...
    !!              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO.
    !!              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES:
    !!
    !!    IWORK(1)  NOT USED
    !!
    !!    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS.
    !!              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000.
    !!
    !!    IWORK(3)  ORDMIN, 3 <= ORDMIN <= 9,
    !!
    !!    IWORK(4)  ORDMAX, ORDMIN <= ORDMAX <= 9
    !!
    !!    IWORK(5)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATIONS FOR THE
    !!              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP FOR ORDER 3.
    !!              THE DEFAULT VALUE (FOR IWORK(5)=0) IS 10.
    !!
    !!    IWORK(6)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
    !!              ORDER 5, THE DEFAULT VALUE (FOR IWORK(6)=0) IS 18.
    !!
    !!    IWORK(7)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
    !!              ORDER 7, THE DEFAULT VALUE (FOR IWORK(7)=0) IS 26.
    !!
    !!    IWORK(8)  THE MAXIMUM NUMBER OF SPLITTING-NEWTON ITERATION FOR
    !!              ORDER 9, THE DEFAULT VALUE (FOR IWORK(5)=0) IS 36.
	!!   
	!!    IWORK(9)  NOT YET USED 
    !!
    !!    IWORK(10:24) USED AS OUTPUT PARAMETERS (SEE BELOW)   
    !!
    !!       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR
    !!       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1.
    !!       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT
    !!       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER. 
    !!       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE
    !!       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2.
    !!
    !!    IWORK(25)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR 
    !!               ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM.
    !!               DEFAULT IWORK(25)=N.
    !!
    !!    IWORK(26)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(26)=0.
    !!
    !!    IWORK(27)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(27)=0.
    !!
    !!
    !!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16.
    !!
    !!    WORK(2)   HMAX  MAXIMAL STEP SIZE, DEFAULT TEND-T0.
    !!
    !!    WORK(3)   THET DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
    !!
    !!    WORK(4)   FACNEWT:  stopping criterion for splitting-Newton method
    !!                 for small values of min(abs(y_i)) and min(abs(f_j)).
    !!
    !!    WORK(5)   TETAK0(1) stopping criterium for the order 3 
	!!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(1)
    !!
    !!    WORK(6)   TETAK0(2) stopping criterium for the order 5 
	!!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(2)
    !!
    !!    WORK(7)   TETAK0(3) stopping criterium for the order 7 
	!!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(3)
    !!
    !!    WORK(8)   TETAK0(4) stopping criterium for the order 9
	!!              splitting-Newton method:
    !!              the iterates must be decreasing by a factor tetak0(4)
    !!
    !!    WORK(9)   CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER  5 METHOD
    !!              USED IN THE ORDER VARIATION STRATEGY
    !!              (DEFAULT WORK(6) = 2.4D0)
    !!
    !!    WORK(10)  CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER  7 METHOD
    !!              USED IN THE ORDER VARIATION STRATEGY
    !!              (DEFAULT WORK(6) = 4.0D0)
    !!
    !!    WORK(11)  CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER  9 METHOD
    !!              USED IN THE ORDER VARIATION STRATEGY
    !!              (DEFAULT WORK(6) = 7.2D0)
    !!
    !!    WORK(12)-WORK(13)   FACL-FACR: PARAMETERS FOR STEP SIZE SELECTION
    !!               THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
    !!               FACL  <=  HNEW/HOLD <= FACR
    !!               (DEFAULT WORK(9) = 0.12, WORK(10) = 10 )
    !!
    !!    WORK(14)  SFDOWN:SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE LOWER ORDER METHOD
    !!                  (DEFAULT WORK(11) = 20D0)
    !!
    !!    WORK(15)  SFUP:SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE UPPER ORDER METHOD
    !!                  (DEFAULT WORK(12) = 40D0)
    !!
    !!    WORK(16)  SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE CURRENT ORDER METHOD
    !!                  (DEFAULT WORK(13) = 18D0)
    !!
    !!    WORK(17)  SF: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!                  USED FOR THE CURRENT ORDER METHOD WHEN IS
    !!                  FAILED THE ERROR CONTROL TEST (DEFAULT WORK(14) = 15D0)
    !!
    !!
    !!
    !!    WORK(18)  FACNEWT stopping criterion for splitting-Newton method ORDER 3
    !!                  (DEFAULT WORK(15) = 1.0D-3)
    !!
    !!    WORK(19)  FACNEWT stopping criterion for splitting-Newton method ORDER 5
    !!                  (DEFAULT WORK(16) = 9.0D-2)
    !!
    !!    WORK(20)  FACNEWT stopping criterion for splitting-Newton method ORDER 7
    !!                  (DEFAULT WORK(17) = 9.0D-1)
    !!
    !!    WORK(21)  FACNEWT stopping criterion for splitting-Newton method ORDER 9
    !!                  (DEFAULT WORK(18) = 9.9D-1)
    !!
    !! -------------------------------------------------------------------------
    !!     OUTPUT PARAMETERS
    !! -------------------------------------------------------------------------
    !!     T0          T-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED
    !!                 (AFTER SUCCESSFUL RETURN T0=TEND).
    !!
    !!     Y(N)        NUMERICAL SOLUTION AT T0
    !!
    !!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP
    !!
    !!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
    !!                   IDID= 1  COMPUTATION SUCCESSFUL,
    !!                   IDID=-1  INPUT IS NOT CONSISTENT,
    !!                   IDID=-2  LARGER NMAX IS NEEDED,
    !!                   IDID=-3  STEP SIZE BECOMES TOO SMALL,
    !!                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR.
    !!                   IDID=-5  GAM CANNOT HANDLES IERR=-1.
    !!
    !!   IWORK(10)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL
    !!                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)
    !!   IWORK(11)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY
    !!                      OR NUMERICALLY)
    !!   IWORK(12)  NSTEP(1)  NUMBER OF COMPUTED STEPS   ORD 3
    !!   IWORK(13)  NSTEP(2)  NUMBER OF COMPUTED STEPS   ORD 5
    !!   IWORK(14)  NSTEP(3)  NUMBER OF COMPUTED STEPS   ORD 7
    !!   IWORK(15)  NSTEP(4)  NUMBER OF COMPUTED STEPS   ORD 9
    !!   IWORK(16)  NNEWT(1)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 3
    !!   IWORK(17)  NNEWT(2)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 5
    !!   IWORK(18)  NNEWT(3)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 7
    !!   IWORK(19)  NNEWT(4)  NUMBER OF REJECTED STEPS (DUE TO NEWTON CONVERGENCE) 9
    !!   IWORK(20)  NERR(1)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 3
    !!   IWORK(21)  NERR(2)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 5
    !!   IWORK(22)  NERR(3)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 7
    !!   IWORK(23)  NERR(4)   NUMBER OF REJECTED STEPS (DUE TO ERROR TEST) 9
    !!   IWORK(24)  NDEC      NUMBER OF LU-DECOMPOSITIONS
    !!-----------------------------------------------------------------------
    !!     DECLARATIONS
    !! -------------------------------------------------------------------------
    USE PRECISION
    IMPLICIT NONE
    !!
    !!   INPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(IN) :: R, ITOL, IJAC, IMAS, IOUT, IPAR(1), LIWORK, LWORK

    REAL(PREC), INTENT(IN) :: TEND, ATOL(1), RTOL(1), RPAR(1)
    !!
    !!   INPUT/OUTPUT VARIABLES
    !!------------------------------------
    INTEGER, INTENT(IN OUT)    ::  IWORK(LIWORK), IDID,  &
         &  MLJAC, MUJAC, MLMAS, MUMAS
    REAL(PREC), INTENT(IN OUT) :: T0, Y0(R), H, WORK(LWORK)

    !!
    !!   LOCAL VARIABLES
    !!------------------------------------
    REAL(PREC) :: FACNORD(4), HMAX, THET, FACNEWT, TETAK0(4), CS(4), FACL, FACR,       &
         &                  SFDOWN, SFUP, SFSAME, SF, UROUND
    INTEGER :: ORDMIN, ORDMAX, ITINT(4), ITMAX, IJOB, NMAX, LDJAC, LDLU, LDMAS
    INTEGER ::  NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4)
    INTEGER ::  IEYP, IEFP, IEDN,IEF,IEF1, IEJF0,IELU, ISTORE, &
         &   I,IEIPIV, IESC, NIND1, NIND2, NIND3
    LOGICAL ::  ARRET, JBAND, IMPLCT

    !!
    !!   EXTERNAL FUNCTIONS
    !!------------------------------------
    !!EXTERNAL FCN, JAC, MAS, SOLOUT
    INTERFACE
       !!-----------------------------------------------------------------------
       SUBROUTINE fcn(neqn,t,y,dy,ierr,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn, ipar(*)
         INTEGER, INTENT(IN OUT):: ierr
         REAL(PREC), INTENT(IN) :: t,y(neqn),rpar(*)
         REAL(PREC), INTENT(OUT) :: dy(neqn)
       END SUBROUTINE fcn
       !!-----------------------------------------------------------------------
       SUBROUTINE jac(neqn,t,y,jacob,ldim,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn,ldim,ipar(*)
         REAL(PREC), INTENT(IN) :: t,y(neqn),rpar(*)
         REAL(PREC), INTENT(OUT) :: jacob(ldim,neqn)
       END SUBROUTINE jac
       !!-----------------------------------------------------------------------
       SUBROUTINE mas(neqn,am,ldim,rpar,ipar)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: neqn, ldim, ipar(*)
         REAL(PREC), INTENT(IN) :: rpar(*)
         REAL(PREC), INTENT(OUT) :: am(ldim,neqn)
        END SUBROUTINE mas
       !!-----------------------------------------------------------------------
       SUBROUTINE SOLOUT(R,TP,YP,F1,NT1,DBLK,ORD,RPAR,IPAR,IRTRN)
         USE PRECISION
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: R, DBLK, ORD, IPAR(*), IRTRN,NT1
         REAL(PREC), INTENT(IN) :: TP(*),YP(R,*),RPAR(*),F1(R,*)
       END SUBROUTINE solout
    END INTERFACE
    !! -------------------------------------------------------------------------
    !!     SETTING THE PARAMETERS
    !! -------------------------------------------------------------------------

 
    ARRET   = .FALSE.
    !! -------- NMAX := THE MAXIMAL NUMBER OF STEPS -----
    IF (IWORK(2).EQ.0) THEN
       NMAX=100000
    ELSE
       NMAX=IWORK(2)
       IF (NMAX.LE.0) THEN
          WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- ORDMIN  :=  MINIMAL ORDER
    IF (IWORK(3).EQ.0) THEN
       ORDMIN = 1
    ELSE
       ORDMIN=(IWORK(3)-1)/2
        IF ((ORDMIN.LE.0).OR.(ORDMIN.GT.4)) THEN
          WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- ORDMAX :=  MAXIMAL ORDER
    IF (IWORK(4).EQ.0) THEN
       ORDMAX = 4
    ELSE
       ORDMAX=(IWORK(4)-1)/2
       IF ((ORDMAX.LE.0).OR.(ORDMAX.GT.4).OR.(ORDMAX.LT.ORDMIN)) THEN
          WRITE(6,*)' CURIOUS INPUT IWORK(4)=',IWORK(4)
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- ITINT(1) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 3
    IF (IWORK(5).EQ.0) THEN
       ITINT(1)=12
    ELSE
       ITINT(1)=IWORK(5)
       IF (ITINT(1).LE.0) THEN
          WRITE(6,*)' CURIOUS INPUT IWORK(5)=',IWORK(5)
          ARRET=.TRUE.
       END IF
    END IF
    ITMAX = ITINT(1)
    !! -------- ITINT(2) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 5
    IF (IWORK(6).EQ.0) THEN
       ITINT(2)=18
    ELSE
       ITINT(2)=IWORK(6)
       IF (ITINT(2).LT.0) THEN
          WRITE(6,*)' CURIOUS INPUT IWORK(6)=',IWORK(6)
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- ITINT(3) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 7
    IF (IWORK(7).EQ.0) THEN
       ITINT(3)= 26
    ELSE
       ITINT(3)=IWORK(7)
       IF (ITINT(3).LT.0) THEN
          WRITE(6,*)' CURIOUS INPUT IWORK(7)=',IWORK(7)
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS
      NIND1=IWORK(25)
      NIND2=IWORK(26)
      NIND3=IWORK(27)
  
      IF (NIND1.EQ.0) NIND1=R
      IF (NIND1+NIND2+NIND3.NE.R) THEN
         WRITE(6,*)' CURIOUS INPUT FOR IWORK(25,26,27)=',NIND1,NIND2,NIND3
         ARRET=.TRUE.
      END IF

    !! -------- ITINT(4) :=  NUMBER OF SPLITTING-NEWTON ITERATIONS ORD 9
    IF (IWORK(8).EQ.0) THEN
       ITINT(4)= 36
    ELSE
       ITINT(4)=IWORK(8)
       IF (ITINT(4).LT.0) THEN
          WRITE(6,*)' CURIOUS INPUT IWORK(8)=',IWORK(8)
          ARRET=.TRUE.
       END IF
    END IF

    !! -------- UROUND :=  SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0
    IF (WORK(1).EQ.0.0D0) THEN
       UROUND=2.30D-16
    ELSE
       UROUND=WORK(1)
       IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN
          WRITE(6,*)'COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1)
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- HMAX := MAXIMAL STEP SIZE
    IF (WORK(2).EQ.0.D0) THEN
       HMAX=TEND-T0
    ELSE
       HMAX=WORK(2)
       IF (HMAX.GT.TEND-T0) THEN
          HMAX=TEND-T0
       END IF
    END IF
    !! -------- THET:  DECIDE WHETHER THE JACOBIAN SHOULD BE RECOMPUTED
    IF (WORK(3).EQ.0.D0) THEN
       THET = 0.005
    ELSE
       THET=WORK(3)
       IF (THET .GT. 1d0) THEN
          WRITE(6,*)' CURIOUS INPUT WORK(3)=',WORK(3)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNEWT: STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!--------           FOR SMALL VALUES OF min(abs(y_i)) and min(abs(f_j))
    IF (WORK(4).EQ.0.D0) THEN
       FACNEWT= 1d-2 
       FACNEWT=MAX(FACNEWT,EPS/RTOL(1) )
    ELSE
       FACNEWT=WORK(4)
       FACNEWT=MAX(FACNEWT,EPS/RTOL(1) )
       IF (FACNEWT.GE.1.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(4) ',WORK(4)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(1): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!---------             ORDER 3
    IF (WORK(15).EQ.0.D0) THEN
	   SELECT CASE (IMAS)
	   CASE(0)    ! ODE CASE
          FACNORD(1) = 1D-3
	   CASE(1)    ! DAE CASE
          FACNORD(1) = 5d-3
	   END SELECT
       FACNORD(1) = MAX(FACNORD(1) ,EPS/RTOL(1) )
    ELSE
       FACNORD(1) = WORK(15)
       FACNORD(1) = MAX(FACNORD(1) ,EPS/RTOL(1) )
       IF (FACNEWT.GE.1.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(15) ',WORK(15)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(2): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!--------              ORDER 5
    IF (WORK(16).EQ.0.D0) THEN
	   SELECT CASE (IMAS)
	   CASE(0)    ! ODE CASE
          FACNORD(2) = 5d-2
	   CASE(1)    ! DAE CASE
          FACNORD(2) = 5D-2
	   END SELECT
       FACNORD(2) = MAX(FACNORD(2) ,EPS/RTOL(1) )
    ELSE
       FACNORD(2) = WORK(16)
       FACNORD(2) = MAX(FACNORD(2) ,EPS/RTOL(1) )
       IF (FACNORD(2).GE.1.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(16) ',WORK(16)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(3): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!---------             ORDER 7
    IF (WORK(17).EQ.0.D0) THEN
	   SELECT CASE (IMAS)
	   CASE(0)    ! ODE CASE
          FACNORD(3) = 1.0d-1
	   CASE(1)    ! DAE CASE
          FACNORD(3) = 5.0d-2
	   END SELECT
       FACNORD(3) = MAX(FACNORD(3), EPS/RTOL(1) )
    ELSE
       FACNORD(3) = WORK(17)
       FACNORD(3) = MAX(FACNORD(3), EPS/RTOL(1) )
       IF (FACNORD(3).GE.1.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(17) ',WORK(17)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACNORD(4): STOPPING CRITERION FOR SPLITTING-NEWTON METHOD
    !!---------             ORDER 9
    IF (WORK(18).EQ.0.D0) THEN
	   SELECT CASE (IMAS)
	   CASE(0)    ! ODE CASE
          FACNORD(4) = 1.0d-1
	   CASE(1)    ! DAE CASE
          FACNORD(4) = 5.0d-2
	   END SELECT
       FACNORD(4) = MAX(FACNORD(4),EPS/RTOL(1) )
    ELSE
       FACNORD(4) = WORK(18)
       FACNORD(4) = MAX(FACNORD(4),EPS/RTOL(1) )
       IF (FACNORD(4).GE.1.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(18) ',WORK(18)
          ARRET=.TRUE.
       END IF
    END IF

    !!--------- TETAK0(1:4): STOPPING CRITERIUM FOR THE SPLITTING-NEWTON METHOD
    !!---------              THE ERROR IN THE ITERATES MUST BE DECREASING BY A FACTOR
    !!---------              TETAK0(I), i=1,..,4 ACCORDING TO THE SELECTED ORDER
	IF (IMAS==0) THEN 
	   TETAK0(1:4) = 0.9D0                         ! ODE CASE
	ELSE
       TETAK0(1:4) = (/ 0.9D0, 1.5d0, 20D0, 20D0 /)    ! DAE CASE
	END IF
    DO I=1,4
       IF (WORK(4+I).NE.0.D0) THEN
          TETAK0(I) = WORK(4+I)
	   END IF
    END DO
    IF (MINVAL(TETAK0).LE.0.0D0) THEN
       WRITE(6,*)'WRONG INPUT FOR WORK(5:8) ',WORK(5:8)
       ARRET=.TRUE.
    END IF

    CS(1) = 1.0D0
    !!--------- CS(2): EMPIRICAL COMPUTATIONAL COST FOR ORDER 5
    !!---------        USED IN THE ORDER VARIATION STRATEGY.
    IF (WORK(9).EQ.0.D0) THEN
       CS(2) = 2.4D0
    !   IF ((NIND2.NE.0).OR. (NIND3.NE.0)) CS(2) = 2.8d0;
    ELSE
       CS(2) = WORK(9)
       IF (CS(2).LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(9) ',WORK(9)
          ARRET=.TRUE.
       END IF
    END IF
    !!---------  CS(3): EMPIRICAL COMPUTATIONAL COST FOR ORDER 7
    !!---------         USED IN THE ORDER VARIATION STRATEGY.
    IF (WORK(10).EQ.0.D0) THEN
       CS(3) = 4.0D0
    !   IF  ((NIND2.NE.0).OR. (NIND3.NE.0)) CS(2) = 4.5d0;
    ELSE
       CS(3) = WORK(10)
       IF (CS(3).LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(10) ',WORK(10)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- CS(4): EMPIRICAL COMPUTATIONAL COST FOR ORDER 9
    !!---------        USED IN THE ORDER VARIATION STRATEGY.
    IF (WORK(11).EQ.0.D0) THEN
       CS(4) =7.2D0
    !  IF ((NIND2.NE.0).OR. (NIND3.NE.0)) CS(2) = 8d0;
    ELSE
       CS(4) = WORK(11)
       IF (CS(4).LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(11) ',WORK(11)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACL: PARAMETER FOR STEP SIZE SELECTION
    !!---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
    !!---------       FACL <= HNEW/HOLD
    IF (WORK(12).EQ.0.D0) THEN
       FACL = 0.12D0
    ELSE
       FACL = WORK(12)
       IF (FACL.LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(12) ',WORK(12)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- FACR: PARAMETER FOR STEP SIZE SELECTION
    !!---------       THE NEW STEPSIZE IS CHOSEN SUBJECT TO THE RESTRICTION
    !!---------       HNEW/HOLD <= FACR
    IF (WORK(13).EQ.0.D0) THEN
       FACR = 10D0
    ELSE
       FACR = WORK(13)
       IF (FACR.LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(13) ',WORK(13)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SFDOWN: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------         USED FOR THE LOWER ORDER METHOD
    IF (WORK(14).EQ.0.D0) THEN
       SFDOWN = 20.0D0
    ELSE
       SFDOWN = WORK(14)
       IF (SFDOWN.LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(14) ',WORK(14)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SFUP:  SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------        USED FOR THE UPPER ORDER METHOD
    IF (WORK(15).EQ.0.D0) THEN
       SFUP = 40.0D0
    ELSE
       SFUP = WORK(15)
       IF (SFUP.LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(15) ',WORK(15)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SFSAME: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------         USED FOR THE CURRENT ORDER METHOD
    IF (WORK(16).EQ.0.D0) THEN
       SFSAME = 18.0D0
    ELSE
       SFSAME = WORK(16)
       IF (SFSAME.LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(16) ',WORK(16)
          ARRET=.TRUE.
       END IF
    END IF
    !!--------- SF: SAFETY FACTOR IN STEP SIZE PREDICTION
    !!---------     USED FOR THE CURRENT ORDER METHOD WHEN
    !!---------     THE ERROR CONTROL TEST fails
    IF (WORK(17).EQ.0.D0) THEN
       SF = 15.0D0
    ELSE
       SF = WORK(17)
       IF (SF.LE.0.0D0) THEN
          WRITE(6,*)'WRONG INPUT FOR WORK(17) ',WORK(17)
          ARRET=.TRUE.
       END IF
    END IF
    !! -------- CHECK  THE TOLERANCES
    IF (ITOL.EQ.0) THEN
       IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE. EPS) THEN
  
          WRITE (6,*) ' TOLERANCES ARE TOO SMALL'
          ARRET=.TRUE.
       END IF
    ELSE
       DO I=1,R
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE. EPS) THEN
             WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL'
             ARRET=.TRUE.
          END IF
       END DO
    END IF
    !! -------------------------------------------------------------------------
    !!     COMPUTATION OF ARRAY ENTRIES
    !! -------------------------------------------------------------------------
    !! -------- BANDED OR NOT, IMPLICIT ?
    JBAND=MLJAC.LT.R
    IMPLCT=IMAS.NE.0
    !! -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS
    IF (JBAND) THEN
       LDJAC = MLJAC+MUJAC+1
       LDLU  = MLJAC+LDJAC
    ELSE
       LDJAC = R
       LDLU  = R
    END IF
!    IF (JBAND) THEN
!       IJOB=2
!    ELSE
!       IJOB=1
!    END IF
    !! --------  MASS MATRIX
      IF (IMPLCT) THEN
          IF (MLMAS.NE.R) THEN
              LDMAS=MLMAS+MUMAS+1
              IF (JBAND) THEN
                 IJOB=4
              ELSE
                 IJOB=3
              END IF
          ELSE
              MUMAS=R
              LDMAS=R
              IJOB=5
          END IF
    !! ------   BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC"
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN
             WRITE (6,*) 'BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC"'
            ARRET=.TRUE.
          END IF
      ELSE
          LDMAS=0
          IF (JBAND) THEN
             IJOB=2
          ELSE
             IJOB=1
          END IF
      END IF
      LDMAS=MAX(1,LDMAS)


!! -------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK_IN
!! -------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK
!!      IEYP  = 19
!!      IEFP  = IEYP  + 10*R
!!      IEDN  = IEFP  + 10*R
!!      IEF   = IEDN  + R
!!      IEF1  = IEF   + 10*R
!!      IESC  = IEF1  + 10*R
!!      IEJF0 = IESC  + R
!!      IELU  = IEJF0 + R*LDJAC
!!--------- TOTAL STORAGE REQUIREMENT
!!      ISTORE = IELU + R*LDLU - 1
      ISTORE = 21
      IF(ISTORE.GT.LWORK)THEN
         WRITE(6,*)' INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=',ISTORE
         ARRET=.TRUE.
      END IF
!! -------- ENTRY POINTS FOR INTEGER WORKSPACE
!!      IEIPIV=25
!! -------- TOTAL REQUIREMENT
!!      ISTORE=IEIPIV+R-1
      ISTORE = 27
      IF (ISTORE.GT.LIWORK) THEN
         WRITE(6,*)' INSUFF. STORAGE FOR IWORK, MIN. LIWORK=',ISTORE
         ARRET=.TRUE.
      END IF
    !! -------- WHEN A FAIL HAS OCCURED, GAM RETURNs WITH IDID=-1
    IF (ARRET) THEN
       IDID=-1
       RETURN
    END IF

 
    !! -------------------------------------------------------------------------
    !!     CALL TO CORE INTEGRATOR
    !! -------------------------------------------------------------------------
    CALL ETRO(R, FCN, T0, Y0(1), TEND, HMAX, H, RTOL, ATOL, ITOL,                       &
         &   JAC, IJAC, MLJAC, MUJAC, MAS, MLMAS, MUMAS, SOLOUT, IOUT, IDID, NMAX,   &
         &   UROUND, THET, FACNEWT, FACNORD, TETAK0, CS, FACL, FACR, SFDOWN,         &
         &   SFUP, SFSAME, SF, ORDMIN, ORDMAX, ITINT, ITMAX,                         &
         &   IMPLCT, JBAND, IJOB, LDJAC, LDLU, LDMAS, NIND1, NIND2, NIND3,           &
         &   NFCN, NJAC, NSTEP, NNEWT, NERR, NDEC, RPAR, IPAR)

    IWORK(10)= NFCN
    IWORK(11)= NJAC
    IWORK(12)= NSTEP(1)
    IWORK(13)= NSTEP(2)
    IWORK(14)= NSTEP(3)
    IWORK(15)= NSTEP(4)
    IWORK(16)= NNEWT(1)
    IWORK(17)= NNEWT(2)
    IWORK(18)= NNEWT(3)
    IWORK(19)= NNEWT(4)
    IWORK(20)= NERR(1)
    IWORK(21)= NERR(2)
    IWORK(22)= NERR(3)
    IWORK(23)= NERR(4)
    IWORK(24)= NDEC

    RETURN
  CONTAINS
    !!END SUBROUTINE GAM
    !!
    !!--------- END OF SUBROUTINE GAM
    !!
    !! -------------------------------------------------------------------------
    !!     SUBROUTINE  ETRO (Extended trapezoidal Rules of Odd order,
    !!                       that is GAMs)
    !! -------------------------------------------------------------------------
    SUBROUTINE  ETRO(R,FCN,T0,Y0,TEND,HMAX,H,RTOL,ATOL,ITOL,                  &
         &   JAC,IJAC,MLJAC,MUJAC,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID,           &
         &   NMAX,UROUND,THET,FACNEWT,FACNORD,TETAK0,CS,FACL,FACR,SFDOWN,     &
         &   SFUP,SFSAME,SF, ORDMIN,ORDMAX,ITINT,ITMAX,                       &
         &   IMPLCT,JBAND,IJOB,LDJAC,LDLU,LDMAS,NIND1,NIND2,NIND3,            &
         &   NFCN,NJAC,NSTEP,NNEWT,NERR,NDEC,RPAR,IPAR)
      !! -------------------------------------------------------------------------
      !!     CORE INTEGRATOR FOR GAM
      !!     PARAMETERS SAME AS IN GAM WITH ADDED WORKSPACE
      !! -------------------------------------------------------------------------
      !!     DECLARATIONS
      !! ----------------------------------------------------------

      USE LINALGGAMD
      USE SUBGAMD
      IMPLICIT NONE
      !!
      !!   COMMON
      !!------------------------------------
      !!COMMON/LINAL/MLLU,MULU,MDIAG
      !!
      !!   INPUT VARIABLES
      !!------------------------------------
      INTEGER, INTENT(IN) ::  R, ORDMIN, ORDMAX,  ITOL, IJAC, MLJAC, MUJAC, MLMAS, MUMAS, IOUT, &
           &         IPAR(1), ITINT(4),  IJOB, NMAX, LDJAC, LDLU, LDMAS, NIND1, NIND2, NIND3

      REAL(PREC), INTENT(IN)  :: TEND, ATOL(1), RTOL(1), RPAR(1), FACNORD(4),&
           &                  HMAX, THET, FACNEWT, TETAK0(4), CS(4), FACL, FACR,       &
           &                  SFDOWN, SFUP, SFSAME, SF, UROUND


      LOGICAL, INTENT(IN) :: IMPLCT
      !!
      !!   OUTPUT VARIABLES
      !!------------------------------------
      INTEGER, INTENT(OUT) :: NDEC, NFCN, NJAC, NSTEP(4), NNEWT(4), NERR(4)
      !!
      !!   INPUT/OUTPUT VARIABLES
      !!----
      REAL(PREC), INTENT(IN OUT) :: T0, Y0(R), H  

      INTEGER, INTENT(IN OUT) ::  IDID, ITMAX 
      !!
      !!   LOCAL VARIABLES
      !!------------------------------------
      REAL(PREC), ALLOCATABLE :: SCAL(:), YP(:,:), FP(:,:),             &
           &       F(:,:),DN(:), F1(:,:), JF0(:,:), LU(:,:), FMAS(:,:)

      INTEGER, ALLOCATABLE :: IPIV(:)



      INTEGER :: I, J, NSING,                                      &
           &          FAILNI, FAILEI, ORDOLD, ORD, ORD2, ORDN,     &
           &          IT, DBL(4), DBLK, DBLKOLD, ORDDEC,           &
           &          NSTEPS, IRTRN, NT1, IER

      REAL(PREC) :: ERRV(10), TP(11), T1(11), YSAFE, DELT,         &
           &                  THETA, TETAK, TETAKOLD,THETAPREC,    &
           &                  HOLD, HDEC,  ERRNEWT, ERRNEWT0,ESP,  &
           &                  ERRUP, ERRSAME, ERRDOWN, RR, RRN, TH, THN, FACN
      REAL(PREC) ::    DYTH,QNEWT, HHFAC, HACC, R0, HEXTRAP, CPORD(4), ERRSAMEOLD

      INTEGER ::  IFAC,  IERR
      LOGICAL :: JBAND, CALJAC, NEWJAC, JVAI, TER, EXTRAP
      LOGICAL :: VAR, INDEX1, INDEX2, INDEX3
      !!
      !!   EXTERNAL FUNCTIONS
      !!------------------------------------
      EXTERNAL FCN,JAC,MAS,SOLOUT

      !! -------- CONSTANTS


      ALLOCATE( SCAL(R), YP(R,11), FP(R,11),                            &
           &       F(R,11), DN(R), F1(R,11), JF0(LDJAC,R), LU(LDLU,R),  & 
           &       FMAS(LDMAS,R), IPIV(R)  )

      
      MLLU=MLJAC
      MULU=MUJAC
      MDIAG=MLLU + MULU +1
      !! ------- CHECK THE INDEX OF THE PROBLEM ----- 
      INDEX1=NIND1.NE.0
      INDEX2=NIND2.NE.0
      INDEX3=NIND3.NE.0
      !! ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ----------
   
      IF (IMPLCT) THEN
        CALL MAS(R,FMAS,LDMAS,RPAR,IPAR)
        MBDIAG=MUMAS+1
        MBB=MLMAS+MUMAS+1
        MDIFF=MLJAC+MUJAC-MUMAS
      END IF
   
   
      !!--------- DBL(1:4) := SIZE OF THE COEFFICIENT MATRICES DEFINING THE GAMs

      DBL(1) = 4
      DBL(2) = 6
      DBL(3) = 8
      DBL(4) = 9

      ORD  = ORDMIN
      NSTEPS = 0
      NFCN = 0
      NJAC = 0
      NDEC = 0
      NSTEP(1:4)=0
      NNEWT(1:4)=0
      NERR(1:4)=0
     

      !!--------- STARTING VALUES FOR NEWTON ITERATION
      DBLK = DBL(ORD)
      DBLKOLD = DBLK
      H    = MIN( H, ABS(TEND-T0)/DBLK )
      IERR = 0 
      CALL FCN(R,T0,Y0,FP(1,1),IERR, RPAR,IPAR) 
      NFCN = NFCN + 1
      IF (IERR.NE.0) THEN
         IDID=-5
         WRITE(*,*) 'GAM: ERROR: ',&
    &              'cannot handle IERR = -1 in (T0,Y0)'
          RETURN
      ENDIF
    	
	  !! SCALING FACTOR OF THE STEPSIZE FOR THE HIGHER INDEX PROBLEMS
	  !! USED TO DEFINE THE SCALING FACTOR FOR THE COMPUTATION OF THE ERROR

		  CPORD(1) = 1d0
		  CPORD(2) = 1d-1
		  CPORD(3) = 1d-2
		  CPORD(4) = 1d-3

      
	  
      !! -------- NUMBER OF FAILURES IN THE SPLITTING-NEWTON SCHEME
      FAILNI = 0
      !! -------- NUMBER OF FAILURES DUE TO THE ERROR TEST
      FAILEI = 0
      NSING  = 0
      ORDOLD = 2*ORD
      ORDDEC = 2*ORD
      HOLD   = 2*H  
      HACC   = 2*H
      HDEC   = 2*H
      CALJAC = .TRUE.
      EXTRAP = .FALSE.
      ITMAX  = ITINT(ORD)

      CALL FCN(R,T0,Y0,DN(1),IERR,RPAR,IPAR)

	
      !!--------- MAIN LOOP (ADVANCING IN TIME)
      LOOP_TIME : DO
         !! 100   CONTINUE
         !!--------- (EVENTUALLY) COMPUTE THE JACOBIAN MATRIX NUMERICALLY
         NEWJAC = .FALSE.
         ERRSAME = 0d0
         IF (CALJAC) THEN
            IF (IJAC.EQ.0) THEN
               DO I=1,R
                  YSAFE=Y0(I)
                  DELT=SQRT(EPS*MAX(1.D-5,ABS(YSAFE)))
                  Y0(I)=YSAFE+DELT
                  IERR = 0
                  CALL FCN(R,T0,Y0,DN,IERR,RPAR,IPAR)
                  IF (IERR.NE.0) THEN
                    IDID=-5
                     WRITE(*,*) 'GAM: ERROR: ',&
      &              'cannot handle IERR = -1 in numerical Jacobians'
                    RETURN
                  ENDIF
                  IF (JBAND) THEN
                     DO J=MAX(1,I-MUJAC),MIN(R,I+MLJAC)
                        JF0(J-I+MUJAC+1,I) = (DN(J)-FP(J,1))/DELT
                     END DO
                  ELSE                     
                    JF0(1:R,I)=(DN(1:R)-FP(1:R,1))/DELT                   
                  END IF
                  Y0(I)=YSAFE
               END DO
            ELSE
               !! -------- COMPUTE JACOBIAN MATRIX ANALYTICALLY
               
               CALL JAC(R,T0,Y0(1),JF0(1,1),LDJAC,RPAR,IPAR)
              
            END IF
            NJAC = NJAC + 1
            NEWJAC = .TRUE.
         END IF




         !!--------- DEFINE SCAL

         THN = 1d0
         J    = 0
         IF (ITOL.EQ.0) THEN
            DO I=1,R
               RRN = ABS(Y0(I))
               SCAL(I)=ATOL(1)+RTOL(1)*RRN
               IF (RRN .LT. THN) THEN
                  J = I
                  THN = RRN
               ENDIF
            END DO 
           

         ELSE
            DO I=1,R
               RRN = ABS( Y0(I) )
               SCAL(I)=ATOL(I)+RTOL(I)*RRN
               IF (RRN .LT. THN) THEN
                  J = I
                  THN = RRN
               ENDIF
            END DO
           
         END IF



         !!--------- DEFINE  FACN
         FACN = FACNORD(ORD)
         IF (THN .LT. 1d-1) THEN
            IF (ABS(FP(J,1)) .LT. 1d-5) THEN
               FACN = MIN(FACNEWT, FACNORD(ORD) )
            END IF
         END IF





         !!--------- DEFINE TP AND YP
         IER = 1
         DO WHILE (IER .NE. 0) 
          IF (EXTRAP) THEN
            T1(1) = T0+H
            DO I=2,DBLK+1
               T1(I) = T1(I-1)+H
            END DO
            CALL INTERP(R,TP(1),YP(1,1),T1(1),F1(1,1),NT1,DBLKOLD,DBLK,T0,Y0(1))
          ELSE
            TP(1) = T0
            YP(1:R,1) = Y0(1:R)
            DO I=2,DBLK+1
               YP(1:R,I) = Y0(1:R)
               TP(I) = TP(I-1)+H
            END DO
          END IF

          IER = 0
          DO J = 2, DBLK+1
             IERR = 0
             CALL FCN(R,TP(J),YP(1,J),FP(1,J),IERR, RPAR,IPAR)
             IER = MIN( IER, IERR)
          END DO
          NFCN = NFCN + DBLK
          IF (IER .NE. 0) THEN
            H = MIN(HACC,H/2d0)
          END IF
         END DO         
         HEXTRAP  = H 

         !!--------- FACTORIZE THE ITERATION MATRIX
         IF ((ORDDEC.NE.ORD) .OR. (HDEC.NE.H).OR.(NEWJAC) ) THEN            
            IER = 1
            DO WHILE ( IER .NE. 0)
               CALL DECLU(R,JF0(1,1),H,LDJAC,LU(1,1),LDLU,IPIV,FMAS,LDMAS,MLMAS,MUMAS,ORD,IER,IJOB)
               NDEC = NDEC + 1
               IF (IER.NE.0) THEN
                  NSING = NSING + 1
                  IF (NSING.GT.5) THEN
                     WRITE(6,*) 'MATRIX IS REPEATEDLY SINGULAR, IER= ',IER
                     WRITE(6,900) T0
                     IDID=-4
                     GOTO 800
                  ELSE
                     H = MIN(HACC,H/2D0)
                  END IF
               END IF
            END DO
            HDEC = H
            ORDDEC = ORD 
         END IF
         IF (H .NE. HEXTRAP) THEN 

         !!--------- DEFINE AGAIN  TP AND YP
            IER = 1
            DO WHILE (IER .NE. 0) 
            IF (EXTRAP) THEN
              T1(1) = T0+H
              DO I=2,DBLK+1
                 T1(I) = T1(I-1)+H
              END DO
              CALL INTERP(R,TP(1),YP(1,1),T1(1),F1(1,1),NT1,DBLKOLD,DBLK,T0,Y0(1))
            ELSE
              TP(1) = T0
              YP(1:R,1) = Y0(1:R)
              DO I=2,DBLK+1
                YP(1:R,I) = Y0(1:R)
                TP(I) = TP(I-1)+H
              END DO
            END IF
            IER = 0
            DO J = 2, DBLK+1
               IERR = 0
               CALL FCN(R,TP(J),YP(1,J),FP(1,J),IERR, RPAR,IPAR)
               IER = MIN(IERR, IER)
            END DO
            NFCN = NFCN + DBLK
 !           IF (IER .NE. 0) THEN
 !             H = MIN(HACC,H/2d0)
 !           END IF

           END DO
         END IF



		 !! DAE:  SCALE THE TOLERANCES ACCORDING TO THE INDEX OF THE VARIABLES
          IF (INDEX2) THEN
             DO I=NIND1+1,NIND1+NIND2
                SCAL(I)=SCAL(I)/(MIN(1D0,H*CPORD(ORD)))
             END DO
          END IF
          IF (INDEX3) THEN
             DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3
                SCAL(I)=SCAL(I)/(MIN(1D0,H*H*CPORD(ORD)**2))
             END DO
          END IF
         

         !!---------- COMPUTE THE NUMERICAL SOLUTION AT TIMES T1(1)...T1(DBLK)
         !!---------- DEFINE VARIABLES NEEDED IN THE ITERATION
         ERRNEWT  = FACN+1d0
         ERRNEWT0 = FACN+1d0
         !!----------
         TETAK    = 1.0D0
         THETA    = 1.0D0
         TETAKOLD = 1.0D0
         THETAPREC = THETA
         ITMAX    = ITINT(ORD)
         IT  = 0

         !!
         !!-------- SPLITTING NEWTON LOOP
         !!
         NEWT_LOOP : DO
            !! 300    CONTINUE
            ERRNEWT0 = ERRNEWT
            ERRNEWT  = 0D0

            !!--------- COMPUTE ONE ITERATION FOR THE SELECTED ORDER
            ORD_NEWT : SELECT CASE(ORD)
            CASE(1)

               CALL      TERMNOT3(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,   &
                    &  ERRNEWT,ERRNEWT0,TETAK0(1),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            CASE(2)
               CALL      TERMNOT5(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,        &
                    &  ERRNEWT,ERRNEWT0,TETAK0(2),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            CASE(3)
               CALL      TERMNOT7(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,        &
                    &  ERRNEWT,ERRNEWT0,TETAK0(3),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            CASE(4)
               CALL      TERMNOT9(R,FCN,H,IT,DN(1), F(1,1),FP(1,1),YP(1,1),TP(1),NFCN,        &
                    &  ERRNEWT,ERRNEWT0,TETAK0(4),LU(1,1), LDLU,IPIV(1),            &
                    &  FMAS(1,1),LDMAS,MLMAS,MUMAS, SCAL(1),IJOB,TER,RPAR,IPAR)

            END SELECT ORD_NEWT

            IF (TER)  THEN
               ERRNEWT = (FACN + 1)
               EXIT NEWT_LOOP
            END IF


            !!--------- COMPUTE TETAK, ETAK

            TETAKOLD = TETAK
            TETAK    = ERRNEWT/ERRNEWT0
       
            IF (IT.LT.2) THEN
               THETA=TETAK0(ORD)/2
            ELSE IF (IT .EQ. 2) THEN
               THETA = TETAK
            ELSE IF (IT .GT. 2) THEN
               THETA = SQRT(TETAK*TETAKOLD)
            END IF  
            IT = IT+1

            JVAI = (IT .LE. ITMAX).AND.(ERRNEWT.GT.FACN) .AND. &
                &((THETA.LT.TETAK0(ORD)).OR.(IT.LE.2)) .AND. &
                &  (ERRNEWT .GT.0d0)
          !   JVAI = (IT .LE. ITMAX).AND.(ERRNEWT.GT.FACN)
            IF (.NOT. (JVAI)) EXIT  NEWT_LOOP
            !!
            !!--------- END OF NEWTON LOOP
            !!
         END DO NEWT_LOOP
         !! 999   CONTINUE
       


            IF (ERRNEWT.GT.FACN) THEN
            !!--------- THE ITERATION DOES NOT CONVERGE
            FAILNI = FAILNI + 1
            NNEWT(ORD) = NNEWT(ORD)+1
            !!--------- CHOICE OF THE NEW STEPSIZE
            HOLD = H
            H=MIN(HACC,H/2d0)
            DBLKOLD = DBLK
            EXTRAP = .FALSE.
            IF (FAILNI .EQ. 1) THEN
               CALJAC = .NOT. NEWJAC
            ELSE
               CALJAC = .FALSE.
            END IF
            !!-------- RETURN TO THE MAIN LOOP
         ELSE
	
            !!--------- THE ITERATION CONVERGES
            !!--------- ERROR ESTIMATION

			ERRSAMEOLD = ERRSAME

            CALL  ESTERR(ERRV, ERRSAME, ERRUP, ERRDOWN, FP, &
                 &     R, H, ORD, DBLK, LU(1,1), LDLU, FMAS(1,1), LDMAS, MLMAS, MUMAS, &
                 &     IPIV(1), F(1,1), F1(1,1), SCAL(1), ORDMAX,ORDMIN,IJOB)
             
            
		    IF (  FAILEI > 5) THEN
			    IF ( ABS(ERRSAME-ERRSAMEOLD) < 1d-4) THEN
				   ERRSAME = 0.8d0
                END IF
			END IF


                
			
            IF ( ERRSAME .GT. 0.8d0  ) THEN
               FAILEI = FAILEI + 1
               NERR(ORD) = NERR(ORD) + 1
               IF (FAILEI .EQ. 1) THEN 
                 CALJAC = (THETA .GT. THET) .and. (.NOT. NEWJAC)
               ELSE
                 CALJAC=.FALSE.
               END IF
               !!           IF (IT .LE. DBLK+1 ) CALJAC = .FALSE.       
               !!--------- NEW STEPSIZE SELECTION
               ORD2 = 2*ORD
               ESP = 1D0/(ORD2+1D0)
               RRN=MAX(FACL,MIN(FACR,(SF*ERRSAME)**ESP))
               HOLD = H
               H = H/RRN
               DBLKOLD = DBLK
               F1(1:R,1:DBLKOLD+1) = YP(1:R,1:DBLKOLD+1)
               CALL DIFFDIV(TP,F1,R,DBLK,NT1)
               EXTRAP = .TRUE.
               
               !!--------- RETURN TO THE MAIN LOOP
            ELSE
			   
               !!--------- THE STEPSIZE IS ACCEPTED
			   
               NSTEP(ORD) = NSTEP(ORD)+1
               T0 = TP(DBLK+1)
               Y0(1:R) = YP(1:R,DBLK+1)
               FP(1:R,1) = FP(1:R,DBLK+1)
               !!--------- NEW STEPSIZE SELECTION
               ORD2 = 2*ORD
               ESP = 1D0/(ORD2+1d0)               
               RRN=MAX(FACL,MIN(FACR,(SFSAME*ERRSAME)**ESP))
               THN=DBL(ORD)/(CS(ORD)*RRN)
               ORDN = ORD

               IF  (ORD.LT.ORDMAX) THEN
                     ESP = 1D0/(ORD2+3D0)
                     RR=MAX(FACL,MIN(FACR,(SFUP*ERRUP)**ESP))
                  
                     TH=DBL(ORD+1)/(CS(ORD+1)*RR )
                     IF (TH .GT. THN ) THEN
                        ORDN = ORD + 1
                        RRN  = RR
                        THN  = TH
                     END IF
               END IF

               IF ( ORD.GT.ORDMIN)  THEN
                  ESP = 1D0/(ORD2-1d0)
                  RR=MAX(FACL,MIN(FACR,(SFDOWN*ERRDOWN)**ESP))    
                  TH=DBL(ORD-1)/(CS(ORD-1)*RR )
                  IF ( (TH .GT. THN )) THEN
                     ORDN = ORD - 1
                     RRN  = RR
                  END IF
               END IF
                
               HOLD = H
               HACC = H
            
               IF (ORDN.GT.ORD) THEN
                  H = MIN(H/RRN,HOLD)                
               ELSE
                   H = H/RRN
               END IF

              
               ORDOLD = ORD
               ORD = ORDN
               DBLKOLD = DBLK
               DBLK = DBL(ORD)

               CALJAC = (THETA .GT. THET)
 !               write(6,*) 'THETA ', THETA,  TETAK0(ORD)/2, THET
 
               IF ((FAILNI.NE.0).OR.(FAILEI.NE.0)) THEN
                  H = MIN( H, HOLD)
               END IF
               IF  (.NOT. CALJAC) THEN
                  IF ((H/HOLD.LE.1.1D0 ).AND.(H/HOLD.GE.0.9D0))THEN 
                     H = HOLD
                  END IF
               END IF
               H = MIN( H, MIN(HMAX, (TEND-T0)/DBLK) )
               
               
              
               F1(1:R,1:DBLKOLD+1) = YP(1:R,1:DBLKOLD+1)
               CALL DIFFDIV(TP,F1,R,DBLKOLD,NT1)
               EXTRAP = .TRUE.

               IF (IOUT.NE.0) THEN
                  !!--------- CALL SOLOUT
                  CALL SOLOUT(R,TP(1),YP(1,1),F1(1,1),NT1,DBLKOLD,ORDOLD,RPAR,IPAR,IRTRN)
                  IF (IRTRN.LT.0) GOTO 800
               END IF
               IF (NSTEPS .EQ. 0) THEN
                  FAILNI = 0
                  FAILEI = 0
               ELSE
                  FAILNI = MAX(FAILNI-1,0)
                  FAILEI = MAX(FAILEI-1,0)
                 
               END IF
               NSING  = 0
            END IF
            !!--------- END IF ERRSAME > 1
         END IF
         !!--------- END IF ERRNEWT > 1
         NSTEPS = NSTEPS + 1
         IF (0.1d0*ABS(T0-TEND)/DBLK .GT. ABS(T0)*EPS ) THEN
		   
            IF (0.1d0*ABS(H) .LE. ABS(T0+EPS)*EPS) THEN
               WRITE(6,*) ' STEPSIZE TOO SMALL, H=',H
               WRITE(6,900) T0
               IDID=-3
               !!  GOTO 800
               EXIT LOOP_TIME
            END IF
            IF (NSTEPS.GT.NMAX) THEN
               WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'
               WRITE(6,900) T0
               IDID=-2
               !!  GOTO 800
               EXIT LOOP_TIME
            END IF
            CYCLE LOOP_TIME
            !! GOTO 100
            !!
            !!---------- END WHILE T0 < T
            !!
         ELSE
            H    = HACC
            IDID = 1
            EXIT LOOP_TIME
         END IF
    END DO LOOP_TIME   
    DEALLOCATE( SCAL, YP, FP, F,DN, F1, JF0, LU, FMAS, IPIV  )

900      FORMAT(' EXIT OF GAM AT T=',E18.4)
800      RETURN
     END SUBROUTINE ETRO

END SUBROUTINE GAMD

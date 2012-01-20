      SUBROUTINE PERIOD (IVPSOL,N,M,T,X,P,EPS,FM,NRW,RW,NIW,IW,
     1                    ICALL,MYP,HES)
C*    Begin Prologue PERIOD
      EXTERNAL IVPSOL
      DOUBLE PRECISION EPS,P,T(M),X(N,M),FM(N,2),RW(NRW)
      DOUBLE PRECISION MYP(*),HES(N,N)
      INTEGER IW(NIW),ICALL(8)
C
C---------------------------------------------------------------------
C
C*  Title
C
C    (Period)ic Solution of Nonlinear Ordinary Differential Equations
C
C*  Written by        P. Deuflhard, R. Winzen
C*  Purpose           Solution of nonlinear two-point boundary value
C                     problems with period solutions of an unknown
C                     period length 
C*  Method            Local Nonlinear two-point Boundary Value
C                     Problems solver (Multiple shooting approach)
C*  Category          I1b3 - Differential and integral equations
C                            Eigenvalue problems
C*  Keywords          Nonlinear boundary value problems, Multiple
C                     shooting, Newton methods, Period solutions
C*  Version           0.5 (Test Version)
C*  Revision          July 1984
C*  Latest Change     January 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Lutz Weimann
C                     ZIB, Division Scientific Computing, 
C                          Department Numerical Analysis and Modelling
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: weimann@zib.de
C
C*    References:
C
C     /1/ P.Deuflhard:
C         Computation of Periodic Solutions of Nonlinear ODE's
C         Univ. Heidelberg, SFB 123, Tech. Rep. 261 (1984)
C
C     /2/ P.Deuflhard, G.Bader:
C         Multiple Shooting Techniques Revisited.
C         Univ. Heidelberg, SFB 123, Tech. Rep. 163 (1982)
C
C     /3/ P.Deuflhard
C         Lectures on Numerical Analysis
C         Univ. Heidelberg, 1982-1984
C
C     /4/ R.Winzen
C         University of Heidelberg, Inst. F. Angewandte Mathematik
C         Diplom Arbeit (1984)
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time.
C    In any case you should not deliver this code without a special
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from aquisition or application of this code.
C
C* Software status
C    This code is not under special care of ZIB and belongs to ZIB
C    software class 3.
C
C     ------------------------------------------------------------
C
C  SUBROUTINES GIVEN BY USER
C---------------------------
C    FCN (T,Z,DZ)          RIGHT-HAND SIDE OF SYSTEM OF FIRST-ORDER
C                          DIFFERENTIAL EQUATIONS
C    DFDY(T,Z,DF)          FUNKTIONAL MATRIX OF RIGHT-HAND SIDE
C
C    REMARK: BOTH ROUTINES MUST HAVE EXACTLY THOSE NAMES
C
C  EXTERNAL SUBROUTINE (TO BE SUPPLIED BY THE USER)
C----------------------
C
C    IVPSOL (N,FCN,T,Y,TEND,EPS,HMAX,H,HS,USCAL,NRW,RW,NIW,IW,IPAR)
C                          INTEGRATOR
C
C    REMARK: IF ICALL(4)=1,2 PERIOD REQUIRES IVPSOL TO SOLVE N*(N+1)
C            ODE'S FOR INTEGRATION OF THE VARIATIONAL EQUATION
C
C  INPUT PARAMETERS (* MARKS INOUT PARAMETERS)
C------------------
C    N             NUMBER OF FIRST-ORDER DIFFERENTIAL EQUATIONS
C    M             NUMBER OF NODES
C                  M.EQ.2    SINGLE SHOOTING
C                  M.GT.2    MULTIPLE SHOOTING
C  * T(M)          SHOOTING NODES IN UNIT-INTERVALL
C  * X(N,M)        INITIAL DATA
C  * P             PERIOD, OR ITS ESTIMATE
C    EPS           REQUIRED RELATIVE PRECISION OF SOLUTION
C    ICALL(1)      CLASSIFICATION OF RIGHT-HAND SIDE (RHS)
C                   0  NON-AUTONOMOUS RHS
C                   1      AUTONOMOUS RHS
C    ICALL(2)      BVP CLASSIFICATION BY USER
C                   0        LINEAR BVP
C                   1        NONLINEAR BVP
C                            GOOD INITIAL DATA AVAILABLE
C                   2        HIGHLY NONLINEAR BVP
C                            ONLY BAD INITIAL DATA AVAILABLE
C    ICALL(3)      MAXIMUM PERMITTED NUMBER OF ITERATIONS
C    ICALL(4)      DIFFERENTIATION PARAMETER FOR JACOBIAN UPDATING
C                   0  DIFFERENCE APPROXIMATION OF WRONSKIAN MATRICES
C                   1  COMPUTATION OF WRONSKIAN MATRICES BY INTEGRATION
C                      OF VARIATIONAL EQUATION (WITHOUT BLOCK-COLUMN)
C                   2  IN ADDITION TO 1 THE FIRST BLOCK-COLUMN OF JACO-
C                      BIAN IS USED (COMP. /4/) ( ONLY INTERESTING FOR
C                      STIFF ODE'S)
C    ICALL(5)      RANK1 UPDATING PARAMETER
C                   0  NO RANK1 UPDATES OF WRONSKIAN MATRICES
C                   1  RANK1 UPDATES ALLOWED IF POSSIBLE
C    ICALL(6)      REFINEMENT PARAMETER
C                   0   NO ITERATIVE REFINEMENT
C                   1 WITH ITERATIVE REFINEMENT
C    ICALL(7)      STIFFNESS PARAMETER
C                   0 IVPSOL IS SUITED TO NON-STIFF PROBLEMS
C                   1 IVPSOL IS SUITED TO     STIFF PROBLEMS
C  * ICALL(8)      PRINT PARAMETER
C                  -1        NO PRINT
C                   0        INITIAL DATA
C                            ITERATIVE VALUES OF LEVEL FUNCTIONS
C                            SOLUTION DATA (OR FINAL DATA, RESPECTIVELY)
C                  +1        ADDITIONALLY
C                            ITERATES T(J),X(I,J), I=1,N,J=1,M
C
C    NRW           DIMENSION OF REAL WORKSPACE RW(NRW)
C                  -----------------------------------------------------
C                  NRW .GE.   2*N*N*(M+2) + 14*M*N + 9*N + 2*M + 4
C                              + N*(N*N+14*N+16)+2  FOR INTEGRATOR
C                  -----------------------------------------------------
C    RW(NRW)       REAL WORKSPACE
C
C    NIW           DIMENSION OF INTEGER WORKSPACE IW(NIW)
C                  --------------------------------------
C                  NIW .GE.   N+1
C                           + N FOR INTEGRATOR
C                  --------------------------------------
C    IW(NIW)       INTEGER WORKSPACE
C
C  OUTPUT PARAMETERS
C-------------------
C    T(M),X(N,M)   SOLUTION DATA (OR FINAL DATA, RESPECTIVELY)
C    P             PERIOD
C    FM(N,2)       FLOQUET MULTIPLIERS
C
C    ICALL(8)      .GT.0     NUMBER OF ITERATIONS PERFORMED
C                            TO OBTAIN THE SOLUTION
C                  .LT.0     PERIOD TERMINATION
C                  -1        ITERATION STOPS AT STATIONARY POINT
C                  -2        ITERATION STOPS AFTER ITMAX ITERATION STEPS
C                            (AS INDICATED BY INPUT PARAMETER ICALL(3))
C                  -3        INTEGRATOR FAILED
C                            TO COMPLETE THE TRAJECTORY
C                  -4        GAUSS-NEWTON METHOD
C                            FAILED TO CONVERGE
C                  -5        ITERATIVE REFINEMENT
C                            FAILED TO CONVERGE
C                  -6        RELIABLE RELATIVE
C                            ACCURACY GREATER THAN 1.D-2
C                  -7        DISCRETE BVP APPEARS TO BE ILL-CONDITIONED
C                            FOR GIVEN SET OF NODES AND
C                            PRESCRIBED RELATIVE ACCURACY
C                  -8        REAL OR INTEGER WORK-SPACE EXHAUSTED
C
C      ADDITIONALLY THE EIGENVALUES OF THE WRONSKIAN
C      AT THE SOLUTION ARE COMPUTED BY 'EISPACK' ROUTINES
C      LIT.: NUM. MATH. 12 349-308 (1968) MARTIN,WILKINSON
C            NUM. MATH. 13 293-304 (1969) PARLETT,REINSCH
C            NUM. MATH. 14 219-231 (1970) MARTIN,PETERS,WILKINSON
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION EPMACH,SMALL
      INTEGER IPAR(10)
      COMMON /MACHIN/ EPMACH, SMALL
      COMMON  /UNIT/  MOUT
C
C---------------------------------------------------------------------
C       MACHINE DEPENDENT CONSTANTS
C      -----------------------------
C  (ADAPTED TO IBM 3081 D, UNIVERSITY OF HEIDELBERG)
C
C  RELATIVE MACHINE PRECISION
      CALL ZIBCONST(EPMACH,SMALL)
C
C  OUTPUT UNIT FOR ITERATION MONITOR
      MOUT=6
C
C-----------------------------------------------------------------------
C  CHECK FOR SUFFICIENT REAL/INTEGER WORKSPACE
C----------------------------------------------
      NDIF1=N*(N*N+14*N+16)+2
      NDIF2=N
      MINRW= 2*N*N*(M+2) + 14*M*N + 9*N + 2*M + 4   + NDIF1
      MINIW= N + 1  + NDIF2
      IF(ICALL(8).GE.0) WRITE(MOUT,1000) MINRW,MINIW
      IF(MINRW.GT.NRW .OR. MINIW.GT.NIW) GOTO 900
C
C  WORKSPACE SPLITTING
C
      N1=N
      IF(ICALL(1).GT.0) N1=N+1
C
      M1=M-1
      NM=N*M
      NM1=N*M1
C
      NP=N+1
      NP2=NP+1
C
      N2=N*N*M1+1
      N3=N2+N*NP2
      N4=N3+N*NP
      N5=N4+N1*N1
      N6=N5+N*N1
      N7=N6+N*N
      N8=N7+NM
      N9=N8+NM
      N10=N9+NM
      N11=N10+NM
      N12=N11+NM
      N13=N12+NM
      N14=N13+NM1
      N15=N14+NM1
      N16=N15+NM1
      N17=N16+NM1
      N18=N17+NM1
      N19=N18+NM1
      N20=N19+N1
      N21=N20+N1
      N22=N21+N1
      N23=N22+N
      N24=N23+N
      N25=N24+N
      N26=N25+N
      N27=N26+N
      N28=N27+M
      N29=N28+M1
      N30=N29+N*NP
      N31=N30+N*NP*M1
      N32=N31+NDIF1
      N33=N32+N1
C
      CALL PRPER(IVPSOL,N,N1,NP,NP2,NDIF1,NDIF2,M,M1,T,X,P,EPS,FM,IPAR
     1 ,ICALL,IW(1)
     2 ,IW(NP2),RW(1),RW(N2),RW(N3),RW(N4),RW(N5),RW(N6),RW(N7),RW(N8)
     3 ,RW(N9),RW(N10),RW(N11),RW(N12),RW(N13),RW(N14),RW(N15),RW(N16)
     4 ,RW(N17),RW(N18),RW(N19),RW(N20),RW(N21),RW(N22),RW(N23),RW(N24)
     5 ,RW(N25),RW(N26),RW(N27),RW(N28),RW(N29),RW(N30),RW(N31)
     6 ,RW(N32),RW(N33),MYP,HES)
C
C  SOLUTION EXIT
      RETURN
C
C  FAIL EXIT  WORK-SPACE EXHAUSTED
900   IF(ICALL(8).GE.0.AND.MINRW.GT.NRW) WRITE(MOUT,1001)
      IF(ICALL(8).GE.0.AND.MINIW.GT.NIW) WRITE(MOUT,1002)
      ICALL(8)=-8
      RETURN
C
1000  FORMAT(30H0 MINIMAL REQUIRED WORK-SPACE:,/,
     1       17H0 REAL ARRAY  RW(,I5,3H)  ,
     2       20H  INTEGER ARRAY  IW(,I4,1H))
1001  FORMAT(42H0 ERROR: PERIOD  REAL WORK-SPACE EXHAUSTED,/)
1002  FORMAT(45H0 ERROR: PERIOD  INTEGER WORK-SPACE EXHAUSTED,/)
C
C  END DRIVER ROUTINE PERIOD
C
      END
C
      SUBROUTINE PRPER(IVPSOL,N,N1,NP,NP2,NDIF1,NDIF2,M,M1,T,X,P,EPS,FM
     1 ,IPAR,ICALL,PIVOT,IW,G,DY,Y,QE,E,BG,DX,DDX,DXQ,DXQA,XA,XW,XU,HH
     2 ,DHH,HHA,FP,FPA,D,T1,DX1,DE,U,DU,QU,T2,RF,HNS1,VH,USCAL,RW,V,XTG
     3 ,MYP,HES)
C
C---------------------------------------------------------------------
C  PERIODIC SOLUTION OF NONLINEAR ORDINARY DIFFERENTIAL EQUATIONS
C---------------------------------------------------------------------
C
      DOUBLE PRECISION T(M),X(N,M), Y(N,NP), DY(N,NP2), VH(N,NP)
     1 ,G(N,N,M1), BG(N,N), E(N,N1), QE(N1,N1),USCAL(N,NP,M1),FM(N,2)
     2 ,DX(N,M), DDX(N,M), DXQ(N,M),  DXQA(N,M), XA(N,M), XW(N,M)
     3 ,XU(N,M1), HH(N,M1), DHH(N,M1), HHA(N,M1), FP(N,M1), FPA(N,M1)
     4 ,D(N1), DE(N), U(N), DU(N), QU(N), T1(N1), T2(N), DX1(N1), RF(M)
     5 ,HNS1(M1),RW(NDIF1),V(N1),XTG(N,M),MYP(*),HES(N,N)
C
      DOUBLE PRECISION      COND  ,CONDH ,COND1 ,CONV  ,CONVA ,
     1 P     ,DP    ,PA    ,DPQ   ,DPQA  ,PW    ,DDP   ,
     2 DABS  ,DEL   ,DSQRT ,EPH     ,EPMACH,EPS   ,EPSMIN,FC    ,
     3 FCA   ,FCDNM ,FCH   ,FCMIN ,FCMINH,FCMIN2,FCNUM ,FCS   ,EIGHT ,
     4 H     ,HALF  ,HMAX  ,HSTART,ONE   ,REDH  ,RELDIF,HS    ,PRSCLP ,
     5 S     ,SENS1 ,SIGDEL,SIGMA ,SKAP  ,SMALL ,ST    ,SUM1  ,SUM2  ,
     6 SUMF  ,SUMX  ,SUMXA ,TEN   ,TFAIL ,TH    ,TJ    ,TJ1   ,TOL   ,
     7 TOLH  ,TOLMIN,TOLF  ,TOLJ  ,TWO   ,TENTH ,XTHR  ,ZERO  ,PTG
C
      INTEGER PIVOT(N1),ICALL(8),IPAR(10),IW(NDIF2)
C
      EXTERNAL IVPSOL,FCN,PRFVAR,PRSCAL
C
      COMMON /MACHIN/ EPMACH, SMALL
      COMMON  /UNIT/  MOUT
      COMMON  /DIM/  NPROB, NALL
C
      DATA  REDH/1.D-2/ , ZERO/0.D0/ , HALF/0.5D0/ , TENTH/1.D-1/,
     1      FCS/0.7D0/ , ONE/1.D0/ , TWO/2.D0/ , EIGHT/8.D0/ ,TEN/1.D1/
C
C---------------------------------------------------------------------
C      INTERNAL PARAMETERS
C    -----------------------
C  STANDARD VALUES FIXED BELOW
C
C  MAXIMUM PERMITTED SUB-CONDITION NUMBER OF MATRIX E
      COND=ONE/EPMACH
C
C  STARTING VALUE FOR PSEUDO-RANK OF SENSITIVITY MATRIX E
      IRANK=N
C
C  PRESCRIBED RELATIVE PRECISION IN PERIOD
C  ADAPTED TO SUBROUTINES DIFEX AND METAN
      EPSMIN=DSQRT(EPMACH)*REDH
      IF(ICALL(4).EQ.0.OR.ICALL(7).NE.0) EPSMIN=DSQRT(EPMACH*TEN)
      IF(EPS.LT.EPSMIN) EPS=EPSMIN
C
C  PRESCRIBED RELATIVE PRECISION FOR NUMERICAL INTEGRATION
      TOL=EPS*TENTH
      TOLMIN=EPSMIN*TENTH
C
C  PRESCRIBED RELATIVE DEVIATION FOR NUMERICAL DIFFERENTIATION
      RELDIF=DSQRT(TOLMIN)
C
C  STARTING VALUE OF RELAXATION FACTOR  (1.D-2 .LE. FC .LE. 1.D0)
C  FOR LINEAR OR MILDLY NONLINEAR PROBLEMS
      FC=1.D0
C
C  MINIMUM PERMITTED VALUE OF RELAXATION FACTOR
      FCMIN=1.D-2
C
C  FOR HIGHLY NONLINEAR PROBLEMS
      IF(ICALL(2).GT.1) FC=FCMIN
C
C  MAXIMUM PERMITTED NUMBER OF ITERATIVE REFINEMENTS SWEEPS
      NYMAX=0
      IF(ICALL(6).GT.0) NYMAX=M-1
C
C
C  DECISION PARAMETER FOR JACOBIAN RANK-1 UPDATES (SIGMA.GT.1.)
C  RANK-1 UPDATES INHIBITED, IF SIGMA.GT.1./FCMIN IS SET
      SIGMA=2.D0
C
C  THRESHOLD SCALING INITIAL VALUE
      XTHR=SMALL
C
C---------------------------------------------------------------------
C
C
C  INITIAL PREPARATIONS
C-----------------------
      NPROB=N
      NALL=N*NP
      IAUTO=ICALL(1)
      NONLIN=ICALL(2)
      ITMAX=ICALL(3)
      INUMV=ICALL(4)
      INUM=INUMV
      IBROY=ICALL(5)
      KPRINT=ICALL(8)
      TOLJ=DSQRT(TOL)
      TOLF=TOL
      FCMIN2=FCMIN*FCMIN
      FCMINH=DSQRT(FCMIN)
      ITER=0
      KOUNT=0
      INIT=0
      LEVEL=0
      IREPET=0
      IRKMAX=0
      IFLO=0
      FCA=FC
      HSTART=(T(2)-T(1))*P*REDH
      SENS1=ZERO
      COND1=ONE
C
      DO 5 J=1,M1
5     HNS1(J)=HSTART
C
      DO 10 I=1,N
10    X(I,M)=X(I,1)
C
      IF(KPRINT.LT.0) GOTO 2000
      WRITE(MOUT,1001)
      WRITE(MOUT,1002)
      DO 200 J=1,M
      TH=T(J)
      IF(IAUTO.EQ.0) TH=TH*P
200   WRITE(MOUT,1003) TH,(X(I,J),I=1,N)
      IF(IAUTO.GT.0) WRITE(MOUT,1007) P
      WRITE(MOUT,1004) N,M,EPS,ITMAX
      WRITE(MOUT,1001)
      IF(KPRINT.GT.0) GOTO 2000
      WRITE(MOUT,1005)
      WRITE(MOUT,1006)
      GOTO 2000
C
C-----------------------------------------------------------------------
C                  PRELIMINARY NEW ITERATE
C-----------------------------------------------------------------------
1000  INIT=1
      DO 1100 J=1,M1
      DO 1100 I=1,N
1100  X(I,J)=XA(I,J)+FC*DX(I,J)
      IF(IAUTO.GT.0) P=PA+FC*DP
      DO 1101 I=1,N
1101  X(I,M)=X(I,1)
      IF(ITER.GT.ITMAX) GOTO 9200
C
C  COMPUTATION OF THE TRAJECTORIES
C  (SOLUTION OF M1 INITIAL VALUE PROBLEMS)
C------------------------------------------
2000  J=1
      IPAR(1)=0
      IPAR(2)=0
      IF(INUMV.GT.0) IPAR(2)=INUMV
      IPAR(3)=0
      IPAR(4)=N
      IPAR(5)=0
      IPAR(6)=0
      KOUNT=KOUNT+1
      H=HSTART
2100  J1=J+1
      TJ=T(J)*P
      TJ1=T(J1)*P
      HMAX=DABS(TJ1-TJ)
      DO 2110 K=1,N
2110  T1(K)=X(K,J)
      CALL IVPSOL(N,FCN,TJ,T1,TJ1,TOLF,HMAX,H,HS,T2,NDIF1,RW,NDIF2,
     1            IW,IPAR,MYP)
      IF(H.NE.0) GOTO 2200
C
C  SINGULAR TRAJECTORY
      KFLAG=-J
      TFAIL=TJ
      IF(INIT.EQ.0) GOTO 9300
      IF(KPRINT.GE.0) WRITE(MOUT,2001)
      FC=FC*HALF
      IF(FC.LT.FCMIN) GOTO 7700
      GOTO 1000
C
C  CONTINUITY CONDITIONS
C------------------------
2200  DO 2210 K=1,N
      TH=T1(K)
      XU(K,J)=TH
2210  HH(K,J)=TH-X(K,J1)
      J=J1
      IF(J.LT.M) GOTO 2100
C
      IF(INIT.EQ.0) GOTO 5100
      LEVEL=1
C
C  COMPUTATION OF CONDENSED RIGHT-HAND SIDE U(N)
C------------------------------------------------
3000  IF(IRANK.GT.0) CALL PRRHSP (N,M1,1,HH,G,U,DE,T1,BG,MYP)
C
C  (BEST) LEAST SQUARES SOLUTION OF LINEAR (N,N1)-SYSTEM
C--------------------------------------------------------
3100  IF(IRANK.GT.0)  CALL PESOLC
     1      (E,N,N1,0,N,N1,DX1,U,IRANK,D,PIVOT,IREPET,QE,T1)
C
      IF(LEVEL.GT.0 .OR. IREPET.NE.0 .OR. IRANK.EQ.0) GOTO 3116
      DO 3115 I=1,IRANK
3115  QU(I)=U(I)
3116  CONTINUE
C
C  DESCALING OF SOLUTION DX1
      DO 3120 L=1,N
3120  DXQ(L,1)=DX1(L)*XW(L,1)
      IF(IAUTO.GT.0) DPQ=PW*DX1(N1)
C
C
C  RECURSIVE COMPUTATION OF DXQ(N,2),...,DXQ(N,M)
C-------------------------------------------------
      CALL PRRECU (N,M,M1,1,IAUTO,HH,G,FP,DXQ,DPQ,T1,T2,MYP)
C
C-----------------------------------------------------------------------
C         ITERATIVE REFINEMENT SWEEPS  NY=1,..,NYMAX
C-----------------------------------------------------------------------
      CALL PRSWEP (N,N1,M,M1,NY,NYMAX,EPS,EPH,HH,G,FP,DXQ,DPQ,DHH,
     &   DU,DE,T1,T2,BG,DX1,NE,IRANK,PIVOT,D,E,QE,SIGDEL,
     &   XW,PW,DDX,RF,LEVEL,RELDIF,TOL,TOLMIN,IAUTO,IREPET,KPRINT,
     &   IERR,MYP)
      GOTO (9500,9600,9700),IERR
      TOLF=TOL
C
C  EVALUATION OF SCALED STANDARD LEVEL FUNCTION SUMF
C----------------------------------------------------
      SUMF=ZERO
      DO 3400 J=1,M1
      J1=J+1
      DO 3400 I=1,N
3400  SUMF=SUMF + (HH(I,J)/XW(I,J1))**2
C
C-----------------------------------------------------------------------
C  PROJECTION FOR MOORE-PENROSE PSEUDOINVERSE OF JACOBIAN
C--------------------------------------------------------
      IF (IAUTO.EQ.0 .OR. IRANK.LT.N) GOTO 3590
      IF (LEVEL.EQ.1) GOTO 3550
      IPIVS=PIVOT(N1)
      DO 3510 I=1,N
      IF (PIVOT(I).NE.N1)
     &   XTG(PIVOT(I),1)=-V(I)*XW(PIVOT(I),1)
      IF (PIVOT(I).EQ.N1) PTG=-V(I)*PW
3510  CONTINUE
      IF (IPIVS.NE.N1) XTG(IPIVS,1)=XW(IPIVS,1)
      IF (IPIVS.EQ.N1) PTG=PW
      DO 3520 I=1,N
      DO 3520 J=1,M1
3520  HH(I,J)=0.D0
      CALL PRRECU (N,M,M1,1,IAUTO,HH,G,FP,XTG,PTG,T1,T2,MYP)
      CALL PRSWEP (N,N1,M,M1,NY,NYMAX,EPS,EPH,HH,G,FP,XTG,PTG,DHH,
     &   DU,DE,T1,T2,BG,DX1,NE,IRANK,PIVOT,D,E,QE,SIGDEL,
     &   XW,PW,DDX,RF,LEVEL,RELDIF,TOL,TOLMIN,IAUTO,IREPET,KPRINT,
     &   IERR,MYP)
      GOTO (9500,9600,9700),IERR
3550  CONTINUE
      S=PRSCLP(N,M,DXQ,XTG,DPQ,PTG,XW,PW,T,IAUTO)
      ST=PRSCLP(N,M,XTG,XTG,PTG,PTG,XW,PW,T,IAUTO)
      S=S/ST
      DO 3570 I=1,N
      DO 3570 J=1,M
3570  DXQ(I,J)=DXQ(I,J)-S*XTG(I,J)
      DPQ=DPQ-S*PTG
3590  CONTINUE
C
C  EVALUATION OF SCALED NATURAL LEVEL FUNCTION SUMX
C  AND SCALED MAXIMUM ERROR NORM CONV
C---------------------------------------------------
4000  SUMX=PRSCLP(N,M,DXQ,DXQ,DPQ,DPQ,XW,PW,T,IAUTO)
      CONV=ZERO
      DO 4010 J=1,M1
      DO 4010 I=1,N
      S=DABS(DXQ(I,J))/XW(I,J)
      IF(CONV.LT.S) CONV=S
4010  CONTINUE
      IF(IAUTO.EQ.0) GOTO 4020
      S=DABS(DPQ/PW)
      IF(CONV.LT.S) CONV=S
C
4020  IF(LEVEL.GT.0) GOTO 4500
C
C-----------------------------------------------------------------------
C          ORDINARY GAUSS-NEWTON CORRECTIONS DX(N,M)
C-----------------------------------------------------------------------
      DO 4110 J=1,M
      DO 4110 I=1,N
      XA(I,J)=X(I,J)
4110  DX(I,J)=DXQ(I,J)
      IF(IAUTO.EQ.0) GOTO 4120
      DP=DPQ
      PA=P
C
C  EVALUATION OF SUBCONDITION AND SENSITIVITY NUMBERS
4120  SUMXA=SUMX
      CONVA=CONV
      COND1=ONE
      SENS1=ZERO
      IF(IRANK.EQ.0) GOTO 4200
      SENS1=DABS(D(1))
      COND1=SENS1/DABS(D(IRANK))
C
C  A-PRIORI ESTIMATE OF RELAXATION FACTOR FC
C--------------------------------------------
4200  JRED=0
      IF(ITER.EQ.0 .OR. (NONLIN.EQ.0.AND.IAUTO.EQ.0)) GOTO 4400
      IF( (NEW.GT.0 .OR. IRANK.LT.N.AND.IRANKA.LT.N)
     1                              .AND.IREPET.EQ.0 ) GOTO 4350
C
C  FULL RANK CASE (INDEPENDENT OF PRECEDING RANK)
C  COMPUTATION OF THE DENOMINATOR OF A-PRIORI ESTIMATE
      DO 4201 J=1,M
      DO 4201 I=1,N
4201  DDX(I,J)=DX(I,J)-DXQA(I,J)
      IF(IAUTO.EQ.1) DDP=DP-DPQA
      FCDNM=PRSCLP(N,M,DDX,DDX,DDP,DDP,XW,PW,T,IAUTO)
C
C  COMPUTATION OF THE PROJECTED DENOMINATOR OF A-PRIORI ESTIMATE
      IF(IRANK.LT.N) GOTO 4350
      IF(IAUTO.EQ.0) GOTO 4300
      SUM1=PRSCLP(N,M,DXQA,XTG,DPQA,PTG,XW,PW,T,IAUTO)
      SUM2=PRSCLP(N,M,XTG,XTG,PTG,PTG,XW,PW,T,IAUTO)
      DEL=SUM1*SUM1/SUM2
      FCDNM=FCDNM-DEL
C
4300  FC=FCA/FCMIN
      IF(FCDNM.GT.FCNUM*FCMIN2) FC=DSQRT(FCNUM/FCDNM)*FCA
C
4350  IREPET=0
      IF(FC.LT.FCMIN) GOTO 7700
      IF(FC.GT.FCS) FC=ONE
C
4400  IF(KPRINT.LT.0) GOTO 1000
      IF(KPRINT.GT.0) WRITE(MOUT,1005)
      IF(KPRINT.GT.0) WRITE(MOUT,1006)
      WRITE(MOUT,4401) ITER,NY,SUMF,SUMXA,NEW,IRANK
      IF(KPRINT.GT.0) WRITE(MOUT,1005)
      GOTO 1000
C
C-----------------------------------------------------------------------
C             SIMPLIFIED GAUSS-NEWTON CORRECTIONS DXQ(N,M)
C-----------------------------------------------------------------------
C
C  RANK INDEPENDENT CONVERGENCE TEST
C
4500  IF(CONV.LE.EPS.AND.IRKMAX.EQ.N) GOTO 9000
C
C  NATURAL MONOTONICITY TEST
C
      IF(SUMX.LE.SUMXA) GOTO 5000
C
C  REDUCTION OF RELAXATION FACTOR FC
C------------------------------------
      IF(KPRINT.LT.0) GOTO 4600
      IF(KPRINT.EQ.0) GOTO 4610
      WRITE(MOUT,1005)
      WRITE(MOUT,1006)
4610  WRITE(MOUT,5001) ITER,NY,SUMF,SUMX,FC
      IF(KPRINT.GT.0) WRITE(MOUT,1005)
C
4600  JRED=JRED+1
      IF(NONLIN.EQ.0.AND.IAUTO.EQ.0) GOTO 9600
      TH=DSQRT(SUMX/SUMXA)
      TH=DSQRT(EIGHT*(TH+FC-ONE)/FC+ONE)-ONE
      FC=FC/TH
      IF(FC.LT.FCMIN .OR. NEW.GT.0.AND.JRED.GT.1) GOTO 7700
      GOTO 1000
C-----------------------------------------------------------------------
C           PREPARATIONS TO START THE FOLLOWING ITERATION STEP
C-----------------------------------------------------------------------
5000  ITER=ITER+1
      LEVEL=0
C
      IF(KPRINT.LT.0) GOTO 5100
      IF(KPRINT.EQ.0) GOTO 5010
      WRITE(MOUT,1005)
      WRITE(MOUT,1006)
5010  WRITE(MOUT,5001) ITER,NY,SUMF,SUMX,FC
      IF(KPRINT.EQ.0) GOTO 5100
      WRITE(MOUT,1005)
      DO 5020 J=1,M
      TH=T(J)
      IF(IAUTO.EQ.0) TH=TH*P
5020  WRITE(MOUT,1003) TH,(X(I,J),I=1,N)
      IF(IAUTO.GT.0) WRITE(MOUT,1007) P
C
C  SCALING OF VARIABLES X(N,M)
5100  CALL PRSCAL (N,M,M1,X,XU,XW,XTHR,MYP)
      IF(IAUTO.GT.0) PW=P
C
      IF(INIT.EQ.0) GOTO 6000
C
C  SAVING OF VALUES DXQ(N,M)
      DO 5200 J=1,M
      DO 5200 I=1,N
5200  DXQA(I,J)=DXQ(I,J)
      IF(IAUTO.GT.0) DPQA=DPQ
C
C  PRELIMINARY PSEUDO-RANK
      IRANKA=IRANK
      IF(IRANK.GE.0.AND.FC.GT.FCMINH) IRANK=N
C
C  A-POSTERIORI ESTIMATE OF RELAXATION FACTOR FC
C------------------------------------------------
      TH=FC-ONE
      DO 5400 J=1,M
      DO 5400 I=1,N
5400  DDX(I,J)=DXQ(I,J)+TH*DX(I,J)
      IF(IAUTO.EQ.1) DDP=DPQ+TH*DP
      FCNUM=PRSCLP(N,M,DX,DX,DP,DP,XW,PW,T,IAUTO)
      FCDNM=PRSCLP(N,M,DDX,DDX,DDP,DDP,XW,PW,T,IAUTO)
5410  FCH=DSQRT(FCNUM/FCDNM)*FC*FC*HALF
C
C  DECISION CRITERION FOR JACOBIAN UPDATING TECHNIQUE
C  INUM=0: NUMERICAL DIFFERENTIATION, INUM=1,2: INTEGRATION OF
C  VARIATIONAL EQUATION, INUM=3: RANK1 UPDATING
C--------------------------------------------------------------
      IF(IBROY.EQ.0) GOTO 5420
      INUM=3
      IF(FC.LT.FCA.AND.NEW.GT.0 .OR. FCH.LT.FC*SIGMA
     1     .OR. EPH*REDH.GT.EPS .OR. IRANK.GT.IRANKA)  INUM=INUMV
5420  FCA=FC
      IF(NONLIN.GT.0) FC=FCH
C
6000  IRANKA=IRANK
C
6400  IF(INUM.EQ.3) GOTO 6700
C
C  DIFFERENCE APPROXIMATION OF WRONSKIAN MATRICES G(1),...,G(M1)
C----------------------------------------------------------------
6500  NEW=0
      KFLAG=0
      IF(INUM.GT.0) GOTO 6600
      CALL PRDERG (N,M,M1,T,X,P,XU,XW,T2,TFAIL,T1,HSTART,G,IVPSOL
     1                ,TOL,RELDIF,KFLAG,IPAR,NDIF1,RW,NDIF2,IW,MYP)
      IF(KFLAG.LT.0) GOTO 9310
C
      KOUNT=KOUNT+N
C
      GOTO 7000
C
C  COMPUTATION OF WRONSKIAN MATRICES G(1),...,G(M1) BY
C   NUMERICAL INTEGRATION OF THE VARIATIONAL EQUATION
C------------------------------------------------------
C
C  ADAPTION OF TOLJ AND TOLF
6600  IF(INIT.EQ.0) GOTO 6610
      TOLJ=DSQRT(SUMX)
      IF(TOLJ.GT.REDH) TOLJ=REDH
      IF(TOLJ.LT.TOLMIN) TOLJ=TOLMIN
      TOLF=TOL*TENTH
6610  CALL PRVARG (N,NP,M,M1,INIT,INUMV,T,X,P,TFAIL,USCAL,VH,HNS1,G,Y
     1            ,IVPSOL,EPMACH,TOLJ,KFLAG,IPAR,NDIF1,RW,NDIF2,IW,MYP)
      IF(KFLAG.LT.0) GOTO 9310
C
      KOUNT=KOUNT+N
C
      GOTO 7000
C
C  RANK-1 UPDATES OF WRONSKIAN MATRICES G(1),...,G(M1)
C------------------------------------------------------
6700  NEW=NEW+1
      CALL PRRK1G (N,M,M1,IAUTO,XW,DX,DP,HH,FP,HHA,FPA,T1,G,FCA,MYP)
C
C
C  COMPUTATION OF SENSITIVITY MATRIX E
C--------------------------------------
C
7000  IF(IRANK.EQ.0) GOTO 7500
C
C  STORING ROW SCALING VECTOR
      DO 7100 I=1,N
7100  DE(I)=SMALL/XW(I,1)
C
C  SCALED MATRIX PRODUCT OF WRONSKIAN MATRICES
      CALL PRGMUL (N,M,M1,G,DE,E,T2,MYP)
C
      IF(IFLO.GT.0) GOTO 930
C
C  INTERNAL ROW AND COLUMN SCALING OF MATRIX E
      DO 7200 K=1,N
      S=XW(K,1)
      DO 7300 I=1,N
7300  E(I,K)=-E(I,K)*S
7200  E(K,K)=E(K,K)+SMALL
C
C  EXTENDED MATRIX E
      IF(IAUTO.EQ.0) GOTO 7500
      DO 7400 J=1,M1
      J1=J+1
      TJ1=T(J1)*P
      DO 7420 K=1,N
7420  T1(K)=XU(K,J)
      CALL FCN(TJ1,T1,T2,MYP)
      DO 7440 K=1,N
7440  FP(K,J)=T2(K)*(T(J1)-T(J))
7400  CONTINUE
      CALL PRRHSP(N,M1,1,FP,G,T2,DE,T1,BG,MYP)
      DO 7450 I=1,N
7450  E(I,N1)=-T2(I)*PW
C
C  MONITOR FOR ACTUALLY APPLIED MAXIMUM RANK
7500  IF(IRKMAX.LT.IRANK) IRKMAX=IRANK
C
C  SAVE VALUES OF FP(N,M1) AND HH(N,M1)
      IF(IREPET.NE.0) GOTO 7600
      DO 7510 I=1,N
      DO 7510 J=1,M1
7510  HHA(I,J)=HH(I,J)
      IF(IAUTO.EQ.0) GOTO 7600
      DO 7520 I=1,N
      DO 7520 J=1,M1
7520  FPA(I,J)=FP(I,J)
C
C  QR-DECOMPOSITION OF (N,N1)-MATRIX E
C--------------------------------------
7600  IF(IRANK.LE.0) GOTO 7610
      CONDH=COND
      CALL PEDECC (E,N,N1,0,N,N1,IRANK,CONDH,D,PIVOT,IREPET,QE,V)
7610  IF(IREPET) 3100,3000,3000
C
C  RESTORE FORMER VALUES
C------------------------
7700  IREPET=1
      LEVEL=0
      IF(IAUTO.EQ.0) GOTO 7720
      DO 7710 I=1,N
      DO 7710 J=1,M1
7710  FP(I,J)=FPA(I,J)
      P=PA
7720  DO 7730 I=1,N
      X(I,1)=XA(I,1)
      DO 7730 J=1,M1
      J1=J+1
      X(I,J1)=XA(I,J1)
      XU(I,J)=X(I,J1)+HHA(I,J)
7730  HH(I,J)=HHA(I,J)
      IF(KPRINT.GE.0) WRITE(MOUT,7701) ITER,FC,IRANK
      IF(ITER.EQ.0) FC=FCMIN
      INUM=INUMV
      IF(NEW.GT.0) GOTO 6500
C
C  PSEUDO-RANK REDUCTION
C------------------------
      IREPET=-1
      IF(IRANK.EQ.0) GOTO 9400
      DO 7750 I=1,IRANK
7750  U(I)=QU(I)
      IRANK=IRANK-1
      GOTO 7600
C
C
C
C------------------------ EXIT -----------------------------------------
C
C-----------------------------------------------------------------------
C                    SOLUTION EXIT
C-----------------------------------------------------------------------
C
9000  ITER=ITER+1
      DO 900 J=1,M1
      DO 900 I=1,N
900   X(I,J)=X(I,J)+DXQ(I,J)
      DO 910 I=1,N
910   X(I,M)=X(I,1)
      IF(IAUTO.GT.0) P=P+DPQ
C
C  COMPUTATION OF THE FLOQUET MULTIPLIERS
C-----------------------------------------
C
      IFLO=1
      INIT=0
      INUM=INUMV
      TOLJ=TOL
      DO 920 J=1,M1
920   HNS1(J)=HSTART
      GOTO 2000
C
930   DO 940 I=1,N
      DO 940 K=1,N
940   E(I,K)=E(I,K)/DE(K)
C
C  COMPUTATION OF MULTIPLIERS OF WRONSKIAN (SIMILARITY TRANSFORMED)
C    BY CALLING STANDARD SOFTWARE FOR COMPUTATION OF EIGENVALUES
C-------------------------------------------------------------------
      DO 941 I=1,N
      DO 941 K=1,N
941   HES(I,K)=E(I,K)

      CALL PRBALA(N,N,E,LOW,IGH,T2,MYP)
      CALL PRORTH(N,N,LOW,IGH,E,T2,MYP)
      CALL PRHQR(N,N,LOW,IGH,E,T1,T2,IERR,MYP)
C
      DO 950 I=1,N
      FM(I,1)=T1(I)
950   FM(I,2)=T2(I)
C
C  SPECIAL VALUES OF SOLUTION
C-----------------------------
      IF(IRANK.LT.N.AND.KPRINT.LT.0) GOTO 9100
      ICALL(8)=ITER
      IF(KPRINT.LT.0) RETURN
      IF(KPRINT.EQ.0) GOTO 9010
      WRITE(MOUT,1005)
      WRITE(MOUT,1006)
9010  WRITE(MOUT,5001) ITER,NY,SUMF,SUMX,FC
      WRITE(MOUT,1005)
      WRITE(MOUT,1001)
      IF(IRANK.LT.N) GOTO 9100
      WRITE(MOUT,9001) ITER,KOUNT
      WRITE(MOUT,9002) CONV
      WRITE(MOUT,9007) TOL
      IF(EPH.GT.CONV) CONV=EPH
      WRITE(MOUT,9003) CONV
9020  J1=1
      SMALL=ONE/SMALL
      IF(SENS1.GT.ONE) GOTO 9038
      SENS1=SENS1*SMALL
      WRITE(MOUT,9004) J1,IRANK,COND1,J1,IRANK,SENS1
      GOTO 9039
9038  WRITE(MOUT,9049) J1,IRANK,COND1,J1,IRANK,SENS1,SMALL
9039  WRITE(MOUT,1001)
      IF(ICALL(8).GT.0) WRITE(MOUT,9005)
      IF(ICALL(8).LT.0) WRITE(MOUT,9905)
      DO 9031 K=1,M
      TH=T(K)
      IF(IAUTO.EQ.0) TH=TH*P
9031  WRITE(MOUT,1003) TH,(X(I,K),I=1,N)
      IF(IAUTO.GT.0) WRITE(MOUT,1007) P
C
      IF(ICALL(8).LT.-1) RETURN
      WRITE(MOUT,10010)
      DO 10011 I=1,N
10011 WRITE(MOUT,10020) T1(I),T2(I)
      IF(IERR.NE.0) WRITE(MOUT,10030) IERR
      RETURN
C
C-----------------------------------------------------------------------
C                   FAIL EXIT MESSAGES
C-----------------------------------------------------------------------
C
C  RANK-DEFICIENCY: BEST LEAST SQUARES SOLUTION OF BVP OBTAINED
9100  ICALL(8)=-1
      IF(KPRINT.LT.0) GOTO 9900
      WRITE(MOUT,9101)
      WRITE(MOUT,9002) CONVA
      IF(EPH.GT.CONVA) CONVA=EPH
      WRITE(MOUT,9003) CONVA
      IF(ITER.EQ.0) GOTO 9900
      SKAP=ZERO
      IF(FCA.EQ.ONE.AND.FC.EQ.ONE.AND.IRANKA.EQ.IRANK)
     1 SKAP=DSQRT(SUMXA/FCNUM)
      IF(SKAP.GT.ZERO) WRITE(MOUT,9102) SKAP
      GOTO 9900
C
C  TERMINATION AFTER MORE THAN ITMAX ITERATIONS
9200  ICALL(8)=-2
      IF(KPRINT.GE.0) WRITE(MOUT,9201) ITMAX
      GOTO 9900
C
C  SINGULAR TRAJECTORY
9310  IF(INUM.EQ.0) WRITE(MOUT,9302)
      IF(INUM.GT.0) WRITE(MOUT,9303)
9300  ICALL(8)=-3
      IF(KPRINT.LT.0) GOTO 9900
      J1=-KFLAG
      WRITE(MOUT,9301) J1,TFAIL
      GOTO 9900
C
C  CONVERGENCE FAIL OF GAUSS-NEWTON METHOD
9400  ICALL(8)=-4
      IF(KPRINT.GE.0) WRITE(MOUT,9401)
      GOTO 9900
C
C  CONVERGENCE FAIL OF ITERATIVE REFINEMENT SWEEPS
9500  ICALL(8)=-5
      IF(KPRINT.LT.0) GOTO 9900
      WRITE(MOUT,9601)
      JN=JN-1
      IF(JN.GT.0) WRITE(MOUT,9602) JN
      GOTO 9900
C
C  INSUFFICIENT ERROR TOLERANCE FOR INTEGRATOR
9600  ICALL(8)=-6
      IF(KPRINT.LT.0) GOTO 9900
      TOLH=EPS/SIGDEL
      RELDIF=DSQRT(TOLH/SIGDEL)
      WRITE(MOUT,3203) TOLH,RELDIF
      WRITE(MOUT,9704) TOLH
      WRITE(MOUT,9703)
      GOTO 9900
C
C  ILL-CONDITIONED DISCRETE BOUNDARY VALUE PROBLEM
9700  ICALL(8)=-7
      IF(KPRINT.GE.0) WRITE(MOUT,9801)
C
C COMMON FAIL EXIT
9900  IF(KPRINT.LT.0) RETURN
      GOTO 9020
C
C-----------------------------------------------------------------------
1001  FORMAT(1H1)
1002  FORMAT(14H0 INITIAL DATA,//)
1003  FORMAT(1H ,D13.5,3(D20.10),5(/14X,3(D20.10)))
1004  FORMAT(4H0 N=,I2,3H M=,I2,/,31H  PRESCRIBED RELATIVE PRECISION,
     1   D10.2,/,45H  MAXIMUM PERMITTED NUMBER OF ITERATION STEPS,I3//)
C
C  ITERATION MONITOR
1005  FORMAT(2H0 ,72(1H*))
1006  FORMAT(2H0 ,4X,2HIT,4X,2HNY,7X,6HLEVELF,9X,6HLEVELX,
     1                               7X,7HREL.FC.,3X,3HNEW,4X,4HRANK)
1007  FORMAT(2H0 ,24X,8H PERIOD:,D20.10,/)
2001  FORMAT(44H0 SINGULAR TRAJECTORY, RELAXATION FACTOR OR ,
     1                                      19HPSEUDO-RANK REDUCED,/)
3201  FORMAT(22H0 ITERATIVE REFINEMENT,/)
3202  FORMAT(7H0 SWEEP,I3,10H STARTS AT,I3,/,10(/,5D12.3))
3203  FORMAT(31H0 SUGGESTED INTEGRATOR ACCURACY,D10.1,/
     1      ,40H0 SUGGESTED RELATIVE DEVIATION PARAMETER,D10.1/)
3204  FORMAT(36H0 ADAPTED IN THE NEXT ITERATION STEP/)
4401  FORMAT(2H0 ,2(4X,I2),2(5X,D10.3),15X,I2,6X,I2)
5001  FORMAT(2H0 ,2(4X,I2),2(5X,D10.3),6X,F5.3)
7701  FORMAT(2H0 ,4X,I2,31H NOT ACCEPTED RELAXATION FACTOR,
     1                                              11X,F5.3,12X,I2)
C
C  SOLUTION OUTPUT
9001  FORMAT(45H0 SOLUTION OF BOUNDARY VALUE PROBLEM OBTAINED,/,
     1 17H0 PERIOD REQUIRED,I3,21H ITERATION STEPS WITH,I4,
     2 23H TRAJECTORY EVALUATIONS,//)
9002  FORMAT(30H0   ACHIEVED RELATIVE ACCURACY,D10.3)
9003  FORMAT(30H0   RELIABLE RELATIVE ACCURACY,D10.3,/)
9004  FORMAT(17H0  SUBCONDITION (,I2,1H,,I2,2H) ,D10.3,/,
     1       17H0  SENSITIVITY  (,I2,1H,,I2,2H) ,D10.3,/)
9007  FORMAT(30H0 APPLIED INTEGRATOR TOLERANCE,D10.3)
9049  FORMAT(17H0  SUBCONDITION (,I2,1H,,I2,2H) ,D10.3,/,
     1       17H0  SENSITIVITY  (,I2,1H,,I2,2H) ,D10.3,2H *,1PD7.0/)
9005  FORMAT(15H0 SOLUTION DATA,/)
C
C  ERROR MESSAGES
9101  FORMAT(42H0 ITERATION TERMINATES AT STATIONARY POINT,/)
9102  FORMAT(30H0 INCOMPATIBILITY FACTOR KAPPA  ,D10.3,/)
9201  FORMAT(35H0 ITERATION TERMINATES AFTER ITMAX=,I3,
     1                                          17H  ITERATION STEPS)
9301  FORMAT(20H0 PERIOD TERMINATES ,/,
     1     13H0 SUBINTERVAL,I3,26H POSSIBLY INSERT NEW NODE ,D20.11/)
9302  FORMAT(52H0 SINGULAR TRAJECTORY BY DIFFERENCE APPROXIMATION OF,
     1                                     20H THE JACOBIAN MATRIX,/)
9303  FORMAT(52H0 SINGULAR TRAJECTORY BY INTEGRATION OF VARIATIONAL ,
     1                                                  8HEQUATION,/)
9401  FORMAT(39H0 GAUSS NEWTON METHOD FAILS TO CONVERGE,/)
9601  FORMAT(46H0 TERMINATION SINCE ITERATIVE REFINEMENT FAILS,
     1                    12H TO CONVERGE,/,19H  INSERT NEW NODES ,/)
9602  FORMAT(17H0 IN SUBINTERVAL ,I3)
9703  FORMAT(47H0 RELIABLE RELATIVE ACCURACY GREATER THAN 1.D-2,/)
9704  FORMAT(51H0 REDUCE RELATIVE ERROR TOLERANCE FOR INTEGRATOR TO,
     1            D10.1,/29H0 OR INCREASE NUMBER OF NODES/)
9801  FORMAT(49H0 ILL-CONDITIONED DISCRETE BVP   INSERT NEW NODES,/)
9905  FORMAT(12H0 FINAL DATA,/)
10010 FORMAT(1X,/,/,22H0 FLOQUET MULTIPLIERS:,/)
10020 FORMAT(10X,2D14.7)
10030 FORMAT(1X,/,23H0 ERROR CODE OF HQR IS:,I3,/)
C-----------------------------------------------------------------------
C
C  END SUBROUTINE PRPER
C
      END
      SUBROUTINE PRDERG (N,M,M1,T,X,P,XU,XW,XJ,TJ,T1,HSTART,
     1       G,IVPSOL,TOL,RELDIF,KFLAG,IPAR,ND1,RW,ND2,IW,MYP)
C----------------------------------------------------------------------
C  DIFFERENCE APPROXIMATION OF WRONSKIAN MATRICES G(1),..,G(M1)
C  ADAPTED FOR SUBROUTINE PERIOD
C----------------------------------------------------------------------
      DOUBLE PRECISION HMAX,H,HSAVE,HSTART,ONE,S,RELDIF,TJ,TJA,TJ1,TH
     1                ,TOL,DABS,ZERO,P,HS
      DOUBLE PRECISION T(M),X(N,M),G(N,N,M1),XW(N,M),XU(N,M1),XJ(N)
     1                ,T1(N),RW(ND1),MYP(*)
      INTEGER IPAR(10),IW(ND2)
      EXTERNAL FCN,IVPSOL
      DATA  ONE/1.D0/ , ZERO/0.D0/
C
      IPAR(1)=0
      IPAR(2)=0
      IPAR(3)=0
      IPAR(4)=N
      IPAR(5)=0
      IPAR(6)=0
      HSAVE=HSTART
      J=1
50    J1=J+1
      TJA=T(J)*P
      TJ1=T(J1)*P
      HMAX=DABS(TJ1-TJA)
      DO 500 I=1,N
C
      DO 503 K=1,N
503   XJ(K)=X(K,J)
      TH=XJ(I)
      S=XW(I,J)*RELDIF
      IF(TH.LT.ZERO) S=-S
      XJ(I)=TH+S
      S=ONE/S
      H=HSAVE
      TJ=TJA
      CALL IVPSOL (N,FCN,TJ,XJ,TJ1,TOL,HMAX,H,HS,T1,ND1,RW,ND2,IW,
     1             IPAR,MYP)
      IF(H.EQ.ZERO) GOTO 999
      DO 52 K=1,N
52    G(K,I,J)=S*(XJ(K)-XU(K,J))
C
500   CONTINUE
      HSAVE=H
      J=J1
      IF(J.LT.M) GOTO 50
C
      KFLAG=0
      RETURN
999   KFLAG=-J
C  ERROR RETURN
      RETURN
C
C  END SUBROUTINE PRDERG
C
      END
      SUBROUTINE PRVARG (N,NP,M,M1,INIT,INUMV,T,X,P,TJ,USCAL,VH,HNS1,G
     1                ,Y,IVPSOL,EPMACH,TOL,KFLAG,IPAR,ND1,RW,ND2,IW,MYP)
C-----------------------------------------------------------------------
C  INTEGRATION OF VARIATIONAL EQUATION
C  ADAPTED FOR SUBROUTINE PERIOD
C-----------------------------------------------------------------------
      DOUBLE PRECISION G(N,N,M1),Y(N,NP),X(N,M),T(M),VH(N,NP)
     1                 ,USCAL(N,NP,M1),HNS1(M1),EPMACH,HS,U
     2                 ,H,TOL,TJ,TJ1,P,HMAX,ONE,ZERO,RW(ND1),MYP(*)
      INTEGER IPAR(10),IW(ND2)
      EXTERNAL PRFVAR,IVPSOL
      DATA  ZERO/0.D0/, ONE/1.D0/
C
      NINT=N*NP
      IPAR(1)=INIT
      IPAR(2)=1
      IPAR(3)=1
      IPAR(4)=N
      IPAR(5)=0
      IF(INUMV.EQ.2) IPAR(5)=1
      IPAR(6)=0
      DO 100 J=1,M1
C
      J1=J+1
      DO 10 I=1,N
      Y(I,1)=X(I,J)
      DO 20 K=2,NP
20    Y(I,K)=ZERO
      L=I+1
10    Y(I,L)=ONE
C
C  RESTORE SCALING VECTOR
      IF(IPAR(1).EQ.0) GOTO 15
      DO 30 I=1,N
      DO 30 K=1,NP
30    VH(I,K)=USCAL(I,K,J)
C
15    TJ=T(J)*P
      TJ1=T(J1)*P
      H=HNS1(J)
      HMAX=DABS(TJ1-TJ)
      CALL IVPSOL (NINT,PRFVAR,TJ,Y,TJ1,TOL,HMAX,H,HS,VH,
     1             ND1,RW,ND2,IW,IPAR,MYP)
      IF(H.EQ.ZERO) GOTO 99
C
C  STORING SCALING VECTOR AND NEW INITIAL STEPSIZE
      DO 40 I=1,N
      U=VH(I,1)
      IF(U.LT.EPMACH) U=EPMACH
      USCAL(I,1,J)=U
      DO 40 K=2,NP
      U=VH(I,K)
      IF(U.LT.ONE) U=ONE
      USCAL(I,K,J)=U
40    CONTINUE
      HNS1(J)=HS
C
C  STORING WRONSKIAN
      DO 50 I=1,N
      DO 50 K=1,N
      KH=K+1
50    G(I,K,J)=Y(I,KH)
C
100   CONTINUE
      KFLAG=0
      RETURN
C
C  ERROR RETURN
99    KFLAG=-J
      RETURN
C
C  END SUBROUTINE PRVARG
C
      END
      SUBROUTINE PRRK1G (N,M,M1,IAUTO,XW,DX,DP,HH,FP,HHA,FPA,DXJ,G
     1                   ,FCA,MYP)
C-----------------------------------------------------------------------
C  RANK-1 UPDATES OF WRONSKIAN MATRICES G(1),...,G(M1)
C  ADAPTED FOR SUBROUTINE PERIOD
C-----------------------------------------------------------------------
      DOUBLE PRECISION DNM,FCA,FCH,ONE,S,T,ZERO,DP
      DOUBLE PRECISION G(N,N,M1),DX(N,M),XW(N,M),HH(N,M1),FP(N,M1),
     1                 HHA(N,M1),FPA(N,M1),DXJ(N),MYP(*)
      DATA  ZERO/0.D0/ , ONE/1.D0/
C
      FCH=FCA-ONE
C
      DO 100 J=1,M1
      DNM=ZERO
      DO 110 I=1,N
      T=DX(I,J)/XW(I,J)
      DXJ(I)=T/XW(I,J)
110   DNM=DNM+T*T
      DNM=DNM*FCA
      IF(DNM.EQ.ZERO) GOTO 100
      DO 120 K=1,N
      T=DXJ(K)/DNM
      DO 120 I=1,N
      S=G(I,K,J)
      IF(S.EQ.ZERO) GOTO 120
      S=S+T*(HH(I,J)+FCH*HHA(I,J))
      IF(IAUTO.GT.0) S=S+T*FCA*DP*(FPA(I,J)-FP(I,J))
      G(I,K,J)=S
120   CONTINUE
100   CONTINUE
      RETURN
C
C  END SUBROUTINE PRRK1G
C
      END
      SUBROUTINE PRGMUL (N,M,M1,G,DE,E,T1,MYP)
C-----------------------------------------------------------------------
C  SCALED MATRIX MULTPLICATION OF WRONSKIANS:  E=DE*G(M1)*.....*G(1)
C-----------------------------------------------------------------------
      DOUBLE PRECISION  G(N,N,M1),E(N,N),T1(N),DE(N),S,ZERO,MYP(*)
      DATA  ZERO/0.D0/
C
      DO 100 I=1,N
      DO 110 K=1,N
110   E(I,K)=ZERO
100   E(I,I)=DE(I)
C
      DO 200 JJ=1,M1
      J=M-JJ
      DO 200 I=1,N
      DO 210 K=1,N
      S=ZERO
      DO 211 L=1,N
211   S=S+E(I,L)*G(L,K,J)
210   T1(K)=S
      DO 220 K=1,N
220   E(I,K)=T1(K)
200   CONTINUE
C
      RETURN
C
C  END SUBROUTINE PRGMUL
C
      END
      SUBROUTINE PRRHSP (N,M1,JIN,HH,G,U,DE,V,BG,MYP)
C----------------------------------------------------------------------
C  COMPUTATION OF CONDENSED RIGHT-HAND SIDE U(N) OF PERIOD
C----------------------------------------------------------------------
      DOUBLE PRECISION S,ZERO
      DOUBLE PRECISION G(N,N,M1), HH(N,M1), BG(N,N)
     1                ,U(N),DE(N),V(N),MYP(*)
      DATA  ZERO/0.D0/
C
      IF(JIN.GT.M1) RETURN
C
      DO 110 I=1,N
      DO 111 K=1,N
111   BG(I,K)=ZERO
      BG(I,I)=DE(I)
110   U(I)=DE(I)*HH(I,M1)
      IF(M1.EQ.1.OR.JIN.EQ.M1) RETURN
C
      M2=M1-1
      DO 200 JJ=JIN,M2
      J=M2+JIN-JJ
      J1=J+1
      DO 200 I=1,N
      DO 210 K=1,N
      S=ZERO
      DO 211 L=1,N
211   S=S+BG(I,L)*G(L,K,J1)
210   V(K)=S
      S=U(I)
      DO 220 K=1,N
      S=S+V(K)*HH(K,J)
220   BG(I,K)=V(K)
200   U(I)=S
C
C  END SUBROUTINE PRRHSP
C
      RETURN
      END
      SUBROUTINE PRRECU (N,M,M1,JIN,IAUTO,HH,G,FP,DX,DP,U,V,MYP)
C---------------------------------------------------------------------
C  RECURSIVE SOLUTION OF M1 LINEAR (N,N)-SYSTEMS
C  ADAPTED FOR SUBROUTINE PERIOD
C---------------------------------------------------------------------
      DOUBLE PRECISION S,ZERO,DP,MYP(*)
      DOUBLE PRECISION G(N,N,M1),FP(N,M1),DX(N,M),HH(N,M1),U(N),V(N)
      DATA  ZERO/0.D0/
C
      DO 10 I=1,N
10    U(I)=DX(I,1)
C
      DO 100 J=1,M1
      J1=J+1
      DO 110 I=1,N
      S=ZERO
      IF(J.LT.JIN) GOTO 112
      S=HH(I,J)
112   IF(IAUTO.GT.0) S=S+FP(I,J)*DP
      DO 111 K=1,N
111   S=S+G(I,K,J)*U(K)
      V(I)=S
110   DX(I,J1)=S
      DO 120 I=1,N
120   U(I)=V(I)
100   CONTINUE
C
C  END SUBROUTINE PRRECU
C
      RETURN
      END
      SUBROUTINE PRSCAL (N,M,M1,X,XU,XW,XTHR,MYP)
C----------------------------------------------------------------------
C  PROVIDES SCALING XW(N,M) OF VARIABLES X(N,M)
C----------------------------------------------------------------------
C
      DOUBLE PRECISION X(N,M),XW(N,M), XU(N,M1),MYP(*)
      DOUBLE PRECISION DABS,EPMACH,HALF,ONE,RED,SMALL,XMAX,XTHR,ZERO
      COMMON /MACHIN/ EPMACH, SMALL
      DATA ZERO/0.D0/, HALF/0.5D0/, ONE/1.D0/, RED/1.D-2/
C
      DO 220 I=1,N
220   XW(I,1)=DABS(X(I,1))
C
C  ARITHMETIC MEAN FOR XW(N,2),...,XW(N,M)
      DO 221 J=1,M1
      J1=J+1
      DO 221 I=1,N
221   XW(I,J1)=(DABS(X(I,J1))+DABS(XU(I,J)))*HALF
C
C  THRESHOLD DETERMINATION
      DO 222 I=1,N
      XMAX=ZERO
      DO 223 J=1,M
      IF(XMAX.LT.XW(I,J)) XMAX=XW(I,J)
223   CONTINUE
      IF(XMAX.LT.XTHR*RED) XMAX=XTHR
      XMAX=XMAX*RED
      IF(XMAX.LT.SMALL) XMAX=ONE
      DO 224 J=1,M1
      IF(XMAX.GT.XW(I,J)) XW(I,J)=XMAX
224   CONTINUE
      XW(I,M)=XW(I,1)
222   CONTINUE
      XTHR=XMAX
C
C  END SUBROUTINE PRSCAL
C
      RETURN
      END
      SUBROUTINE PRFVAR (T,Y,DY,MYP)
C-----------------------------------------------------------------------
C  RIGHT HAND SIDE FOR CALLING IVPSOL IN SUBROUTINE PRVARG
C-----------------------------------------------------------------------
      DOUBLE PRECISION  Y(1),DY(1),S,T,ZERO,MYP(*)
      COMMON /DIM/ N, NINT
      DATA  ZERO/0.D0/
C
      CALL FCN(T,Y,DY,MYP)
      NH=N+1
      CALL DFDY(T,Y,DY(NH),MYP)
C
      DO 10 L=1,N
      DO 20 I=1,N
      S=ZERO
      DO 30 K=1,N
      K1=K*N+L
      K2=I*N+K
30    S=S+DY(K1)*Y(K2)
      IH=NINT+I
20    DY(IH)=S
      DO 40 K=1,N
      K1=K*N+L
      KH=NINT+K
40    DY(K1)=DY(KH)
10    CONTINUE
C
      RETURN
C
C  END SUBROUTINE PRFVAR
C
      END
      SUBROUTINE PRSWEP (N,N1,M,M1,NY,NYMAX,EPS,EPH,HH,G,FP,DXQ,DPQ,DHH,
     &   DU,DE,T1,T2,BG,DX1,NE,IRANK,PIVOT,D,E,QE,SIGDEL,
     &   XW,PW,DDX,RF,LEVEL,RELDIF,TOL,TOLMIN,IAUTO,IREPET,KPRINT,IERR
     &   ,MYP)
C  ITERATIVE REFINEMENT SWEEPS  NY=1,..,NYMAX
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION HH(N,M1),G(N,N,M1),FP(N,M1),DXQ(N,M),DHH(N,M1)
      DOUBLE PRECISION DU(N),DE(N),T1(N1),T2(N),BG(N,N),MYP(*)
      DOUBLE PRECISION DX1(N1),D(N1),E(N,N1),QE(N1,N1),XW(N,M),DDX(N,M),
     $RF(M)
      INTEGER PIVOT(N1)
      COMMON /UNIT/ MOUT
      DATA REDH/1.D-2/
      NY=0
      SIGDEL=10.D0
      SIGDLH=0.D0
      EPH=EPS
      IF(NYMAX.EQ.0 .OR. IRANK.LT.N) GOTO 9000
C
C  COMPUTATION OF REQUIRED CONTINUITY RESIDUALS DHH(N,M1)
      JN=1
      JIN=M
      DO 3150 I=1,N
3150  DU(I)=0.D0
      GOTO 3230
C
3200  DO 3220 J=JN,M1
      J1=J+1
      DO 3220 I=1,N
      S=HH(I,J)
      IF(IAUTO.EQ.0) S=S+FP(I,J)*DPQ
      DO 3221 K=1,N
3221  S=S+G(I,K,J)*DXQ(K,J)
3220  DHH(I,J)=S-DXQ(I,J1)
C
C  COMPUTATION OF CONDENSED RESIDUAL DU(NE)
      IF(IRANK.GT.0) CALL PRRHSP(N,M1,JIN,DHH,G,DU,DE,T1,BG,MYP)
C
3230  DO 3240 I=1,N
3240  DU(I)=DU(I)+DE(I)*(DXQ(I,1)-DXQ(I,M))
C
C  COMPUTATION OF CORRECTION DDX(N,1)
      IF (IRANK.GT.0) CALL PESOLC
     &   (E,N,N1,0,N,N1,DX1,DU,IRANK,D,PIVOT,IREPET,QE,T1)
C
C  DESCALING OF DDX(N,1), REFINEMENT OF DXQ(N,1)
C
      CORR=0.D0
      DO 3260 L=1,N
      S=DX1(L)
      IF (CORR.LT.DABS(S)) CORR=DABS(S)
      S=S*XW(L,1)
      DDX(L,1)=S
3260  DXQ(L,1)=DXQ(L,1)+S
      IF(IAUTO.EQ.0) GOTO 3262
      S=DX1(N1)
      IF (CORR.LT.DABS(S)) CORR=DABS(S)
      S=S*PW
      DDP=S
      DPQ=DPQ+S
3262  IF (CORR.LT.EPH) GOTO 3269
      EPH=CORR
      GOTO 9800
3269  RF(1)=CORR
C
C  RECURSIVE COMPUTATION OF DDX(N,2),...,DDX(N,M)
      CALL PRRECU(N,M,M1,JIN,IAUTO,DHH,G,FP,DDX,DDP,T1,T2,MYP)
C
C  REFINEMENT OF DXQ(N,2),...,DXQ(N,M)
      DO 3270 J=2,M
      CORR=0.D0
      DO 3271 I=1,N
      S=DDX(I,J)
      DXQ(I,J)=DXQ(I,J)+S
      S=DABS(S)/XW(I,J)
      IF(CORR.LT.S) CORR=S
3271  CONTINUE
3270  RF(J)=CORR
C
C  DETERMINATION OF SWEEP INDEX JN
      JA=JN
      DO 3280 J=1,M
      IF(RF(J).GT.EPH) GOTO 3290
3280  JN=J
C
3290  NY=NY+1
      IF(JN.LE.JA) GOTO 9600
      IF(JN.EQ.M) GOTO 3900
      JIN=JN
      IF(NY.GT.1 .OR. LEVEL.EQ.0) GOTO 3200
C
C  DETERMINATION AND ADAPTATION OF PARAMETERS TOL AND RELDIF
3900  IF(LEVEL.EQ.0 .OR. NY.GT.1) GOTO 3920
      DO 3910 J=1,M1
      S=0.D0
      IF (RF(J).NE.0.D0) S=RF(J+1)/RF(J)
      IF(SIGDLH.LT.S) SIGDLH=S
      RF(J)=S
3910  CONTINUE
      SIGDEL=DMAX1(SIGDLH,SIGDEL)
      TH=TOL*SIGDEL
      IF(TH.GT.REDH) GOTO 9700
      IF(TH.GT.EPH) EPH=TH
      TOLH=EPS/SIGDEL
      IF(TOLH.LT.TOLMIN) TOLH=TOLMIN
      TOL=TOLH
      RELDIF=DSQRT(TOL/SIGDEL)
3920  IF(JN.NE.M) GOTO 3200
9000  IERR=0
      RETURN
C  FAIL EXIT
9600  IERR=1
      IF (KPRINT.GE.0) WRITE(MOUT,60001)
      RETURN
9700  IERR=2
      IF (KPRINT.GE.0) WRITE(MOUT,60002)
      RETURN
9800  IERR=3
      IF (KPRINT.GE.0) WRITE(MOUT,60003)
      RETURN
60001 FORMAT(' ITERATIVE REFINEMENT FAILED TO CONVERGE')
60002 FORMAT(' RELIABLE RELATIVE ACCURACY NOT SUFFICIENT')
60003 FORMAT(' GAUSSIAN BLOCK ELIMINATION FAILED BY ILL-CONDITIONED',
     &   ' CONDENSED LINEAR SYSTEM')
      END
      FUNCTION PRSCLP(N,M,X1,X2,P1,P2,XW,PW,T,IAUTO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION X1(N,M),X2(N,M),XW(N,M),T(M)
      RL=T(M)-T(1)
      PRSCLP=0.D0
      DO 1000 J=1,M
      IF (J.EQ.1) S=T(2)-T(1)
      IF (J.GT.1 .AND. J.LT.M) S=T(J+1)-T(J-1)
      IF (J.EQ.M) S=T(M)-T(M-1)
      S=S/RL
      SUM=0.D0
      DO 1010 I=1,N
      S1=X1(I,J)/XW(I,J)
      S2=X2(I,J)/XW(I,J)
1010  SUM=SUM+S1*S2
1000  PRSCLP=PRSCLP+S*SUM
      PRSCLP=0.5D0*PRSCLP
      IF (IAUTO.EQ.1) PRSCLP=PRSCLP+(P1/PW)*(P2/PW)
      RETURN
      END
C
      SUBROUTINE PRDIFX (N,F,X,Y,XEND,EPS,HMAX,H,HS,
     1                      USCAL,NRW,RW,NIW,IW,IPAR,MYP)
C
C
C  EXPLICIT EXTRAPOLATION INTEGRATOR
C  FOR NON-STIFF SYSTEMS OF ORDINARY DIFFERENTIAL EQUATIONS
C  (BASED ON THE EXPLICIT MID-POINT DISCRETIZATION)
C  OR THEIR VARIATIONAL EQUATION
C
C  LATEST CHANGE: APR  06,'84  BY R. WINZEN
C
C  INTERNAL OPTION FOR POLYNOMIAL OR RATIONAL EXTRAPOLATION
C
C
C REFERENCES:
C
C /1/ W.B.GRAGG:
C     ON EXTRAPOLATION ALGORITHMS FOR ORDINARY INITIAL VALUE PROBLEMS
C     SIAM J. NUMER. ANAL. 2, 384-404 (1965)
C
C /2/ R.BULIRSCH, J.STOER:
C     NUMERICAL TREATMENT OF ORDINARY DIFFERENTIAL EQUATIONS BY
C     EXTRAPOLATION METHODS
C     NUMER. MATH. 8, 1-13 (1966)
C
C /3/ P.DEUFLHARD:
C     ORDER AND STEPSIZE CONTROL IN EXTRAPOLATION METHODS
C     UNIVERSITY OF HEIDELBERG, SFB 123: TECH. REP. 93 (1980)
C
C
C  EXTERNAL SUBROUTINE (TO BE SUPPLIED BY THE USER)
C
C    F (X,Y,DY)         RIGHT-HAND SIDE OF SYSTEM OF FIRST ORDER
C                       ORDINARY DIFFERENTIAL EQUATIONS
C      N                NUMBER OF ODE'S
C      X                ACTUAL POSITION
C      Y(N)             VALUES AT T
C      DY(N)            DERIVATIVES AT T
C
C
C  INPUT PARAMETERS (* MARKS INOUT PARAMETERS)
C
C    N                  NUMBER OF ODE'S
C  * X                  STARTING POINT OF INTEGRATION
C  * Y(N)               INITIAL VALUES Y(1),...,Y(N)
C    XEND               PRESCRIBED FINAL POINT OF INTEGRATION
C    EPS                PRESCRIBED RELATIVE PRECISION (.GT.0)
C    HMAX               MAXIMUM PERMITTED STEPSIZE
C  * H                  INITIAL STEPSIZE GUESS
C                       IF H.EQ.ZERO, THEN DIFEX INTERNALLY
C                       GENERATES AN INITIAL STEPSIZE GUESS H
C  * USCAL(N)           SCALING VECTOR (COMPARE PARAMETER IPAR(1))
C
C    NRW                DIMENSION OF REAL WORKSPACE
C                       ----------------------------
C                         NRW . GE .   N*15 +2*NH
C                       ----------------------------
C    RW(NRW)            REAL WORKSPACE
C
C    NIW                DIMENSION OF INTEGER WORKSPACE
C                       ----------------------------
C                         NIW . GE .   1 (DUMMY)
C                       ----------------------------
C    IW(NIW)            INTEGER WORKSPACE
C
C    IPAR(1)            PARAMETER MANAGES USER-SCALING
C                        0  NO USER SCALING
C                        1  USER SCALING
C    IPAR(2)            DUMMY
C    IPAR(3)            EQUATION STATUS
C                        0  NORMAL ODE
C                        1  VARIATIONAL EQUATION OF THAT ODE
C    IPAR(4)            DIMENSION IN REFERRENCE TO IPAR(3)
C                        N  IF IPAR(3)=0
C                        NH  OTHER CASE  ( THEN  N=NH*(NH+1) )
C    IPAR(5)            DUMMY
C  * IPAR(6)            PRINT PARAMETER
C                        0   NO OUTPUT
C                        1   INTEGRATION MONITOR
C                        2   ADDITIONALLY INTERMEDIATE SOLUTION POINTS
C                            T,Y(I),I=1,N
C
C  OUTPUT PARAMETERS
C
C    X                  ACTUAL FINAL POINT OF INTEGRATION
C    Y(N)               FINAL VALUES AT T
C    H                  STEPSIZE PROPOSAL FOR NEXT INTEGRATION STEP
C                       (H.EQ.0. ,IF DIFEX FAILS TO PROCEED)
C    HS                 INITIAL STEPSIZE PROPOSAL FOR NEXT INTEGRATION
C                       OF THE SAME EQUATION (NAMELY IN ITERATION CODES)
C    USCAL(N)           USED SCALING VECTOR
C                       USCAL(N)=DABS(Y(N))
C    IPAR(6)     .GE. 0:SUCCESSFUL INTEGRATION
C                       (KFLAG NOT ALTERED INTERNALLY)
C                .EQ.-1:MORE THAN NSTMAX BASIC INTEGRATION STEPS PER
C                       INTERVAL HAVE BEEN PERFORMED
C                .EQ.-2:MORE THAN JRMAX STEPSIZE REDUCTIONS
C                       OCCURRED PER BASIC INTEGRATION STEP
C                .EQ.-3:STEPSIZE PROPOSAL FOR NEXT BASIC INTEGRATION
C                       TOO SMALL
C                .EQ.-4 REAL- OR INTEGER-WORKSPACE EXHAUSTED
C
C    IPAR(7)            NUMBER OF INTEGRATION STEPS
C    IPAR(8)            NUMBER OF F-EVALUATIONS
C    IPAR(9)            DUMMY
C    IPAR(10)           DUMMY
C
C  REMARK: ALL DUMMY'S ARE NECCESSARY TO HAVE THE SAME PARAMETER LIST
C          AS SUBROUTINE METAN (STIFF SOLVER) HAS
C
C
      DOUBLE PRECISION Y(N),USCAL(N),RW(NRW),X,XEND,EPS,HMAX,H,HS
     1                ,EPMACH,SMALL,MYP(*)
      INTEGER IW(NIW),IPAR(10)
      EXTERNAL F
C
C---------------------------------------------------------------------
C       MACHINE DEPENDENT CONSTANTS
C      -----------------------------
C  (ADAPTED TO IBM 3081 D, UNIVERSITY OF HEIDELBERG)
C
C  RELATIVE MACHINE PRECISION
      EPMACH=2.D-16
C
C  SQRT(SMALLEST POSITIVE MACHINE NUMBER / EPMACH)
      SMALL=1.D-30
C
C  OUTPUT UNIT FOR ITERATION MONITOR
      LOUT=6
C
C-----------------------------------------------------------------------
C  CHECK FOR SUFFICIENT REAL/INTEGER WORKSPACE
C----------------------------------------------
      NH=IPAR(4)
      MINRW= N*15+2*NH
      IF(MINRW.GT.NRW) GOTO 900
C
C  WORKSPACE SPLITTING
C
      N1=N
      IF(IPAR(3).GT.0) N1=N1+NH
C
      N2=N+1
      N3=N2+N
      N4=N3+N1
      N5=N4+N1
      N6=N5+N
C
      CALL PRDIF (N,N1,F,X,Y,XEND,EPS,HMAX,H,HS,USCAL,EPMACH,SMALL,LOUT,
     1 IPAR,RW(1),RW(N2),RW(N3),RW(N4),RW(N5),RW(N6),MYP)
C
C  SOLUTION EXIT
      RETURN
C
C  FAIL EXIT  WORK-SPACE EXHAUSTED
900   IF(IPAR(6).GE.0.AND.MINRW.GT.NRW) WRITE(LOUT,1001)
      IPAR(6)=-4
      RETURN
C
1001  FORMAT(41H0 ERROR: DIFEX  REAL WORK-SPACE EXHAUSTED,/)
C
C  END DRIVER ROUTINE PRDIFX
C
      END
C
      SUBROUTINE PRDIF (N,N1,F,X,Y,XEND,EPS,HMAX,H,HS,USCAL,EPMACH,SMALL
     1 ,LOUT,IPAR,YL,YM,DY,DZ,S,DT,MYP)
C
      INTEGER NJ(10),INCR(10),NRED( 9),IPAR(10)
C
      DOUBLE PRECISION Y(N),YL(N),YM(N),DY(N1),DZ(N1),DT(N,10),
     1  S(N),USCAL(N),D(10,10),A(10),AL(10,10),MYP(*)
C
      DOUBLE PRECISION B,B1,C,DABS,DBLE,DIFF,DM,DMAX,DSQRT,EPH,EPMACH,
     1 EPS,ERR,FC,FCM,FCO,FIVE,FJ,FJ1,FMIN,FN,G,H,HALF,HMAX,HMAXU,HR,
     2 H1,OMJ,OMJO,ONE,ONE1,Q,QUART,RED,RO,SAFE,SK,SMALL,TA,TEN,
     3 U,V,W,X,XEND,XEPS,XN,YH,ZERO,HS,HS1,HS2,DMAX1
C
      DATA ZERO/0.D0/,FMIN/1.D-2/,RO/0.25D0/,QUART/0.25D0/,HALF/0.5D0/,
     1     SAFE/0.5D0/,ONE/1.D0/,ONE1/1.01D0/,TEN/1.D1/,FIVE/5.D0/
C
      EXTERNAL F
C
C  INTERNAL PARAMETERS
C----------------------
C
C  STEPSIZE SEQUENCE HA (DUE TO /3/ )
      DATA NJ/2,4,6,8,10,12,14,16,18,20/
      EPMACH=1.D-16
C
C  ASSOCIATED MAXIMUM COLUMN NUMBER (1.LE.KM.LE.9)
      KM=6
C
C  ASSOCIATED MAXIMUM ROW NUMBER (2.LE.JM.LE.10)
      JM=KM+1
C
C  POLYNOMIAL EXTRAPOLATION (AITKEN-NEVILLE ALGORITHM)
C  (RECOMMENDED STANDARD OPTION)
      IPOL=1
C  RATIONAL EXTRAPOLATION (STOER ALGORITHM)
C      IPOL=0
C
C  MAXIMUM PERMITTED NUMBER OF INTEGRATION STEPS PER INTERVAL
      NSTMAX=50000
C
C  MAXIMUM PERMITTED NUMBER OF STEPSIZE REDUCTIONS
      JRMAX=5
C
C  INITIAL PREPARATIONS
C-----------------------
      DO 1 I=1,N
      DO 1 J=1,10
1     DT(I,J)=ZERO
      IN=IPAR(1)
      KFLAG=IPAR(6)
      EPH=RO*EPS
      FJ1=DBLE(NJ(1))
      A(1)=FJ1+ONE
      DO 60 J=2,JM
      J1=J-1
      INCR(J1)=0
      NRED(J1)=0
      FJ=DBLE(NJ(J))
      V=A(J1)+FJ
      A(J)=V
      DO 61 K=1,J1
      W=FJ/DBLE(NJ(K))
61    D(J,K)=W*W
      IF(J.EQ.2) GOTO 60
      W=V-FJ1
      DO 62 K1=2,J1
      K=K1-1
      U=(A(K1)-V)/(W*DBLE(K+K1))
      U=EPH**U
62    AL(J1,K)=U
60    CONTINUE
      KOH=1
      JOH=2
65    IF(JOH.GE.JM) GOTO 66
      IF(A(JOH+1)*ONE1.GT.A(JOH)*AL(JOH,KOH)) GOTO 66
      KOH=JOH
      JOH=JOH+1
      GOTO 65
66    K=0
      KM=KOH
      JM=KM+1
      INCR(JM)=-1
      OMJO=ZERO
      IF(KFLAG.GE.1) WRITE(LOUT,1001) EPS,KM
      EPMACH=EPMACH*TEN
      HMAX=DABS(HMAX)
      NSTEP=0
      NFCN=0
      XEPS=(DABS(X)+DABS(XEND))*EPMACH
      FN=DBLE(N)
      H1=XEND-X
      HMAXU=HMAX
      HR=HMAX
      HS2=HMAX
      DMAX=FIVE
C
C  INITIAL SCALING
      DO 7 I=1,N
      IF(IN.GT.0) S(I)=DABS(USCAL(I))
7     USCAL(I)=ZERO
C
C  BASIC INTEGRATION STEP
C-------------------------
401   IF(DABS(H1).LE.XEPS) GOTO 403
      Q=H1/H
      IF(Q.LE.EPMACH) GOTO 403
      IF(KFLAG.GT.1) WRITE(LOUT,1009) NSTEP,NFCN,X,K,KOH
      IF(KFLAG.GT.1) WRITE(LOUT,1000) NSTEP,NFCN,X,(Y(I),I=1,N)
      IF(Q.GE.ONE1) GOTO 402
      HR=H
      H=H1
402   JRED=0
      DO 405 K=1,KM
405   INCR(K)=INCR(K)+1
      HMAX=DABS(H1)
      IF(HMAXU.LT.HMAX) HMAX=HMAXU
C
C  SCALING
C  (FOR REAL LIFE APPLICATIONS TO BE POSSIBLY ALTERED BY THE USER)
      IF(IN.GT.0.AND.NSTEP.EQ.0) GOTO 8
      DO 5 I=1,N
      U=DABS(Y(I))
      USCAL(I)=DMAX1(USCAL(I),U)
      IF(IN.GT.0) GOTO 5
      IF(U.LT.EPMACH) U=ONE
      S(I)=U
5     CONTINUE
C
C
C  EXPLICIT EULER STARTING STEP
8     CALL  F (X,Y,DZ,MYP)
      NFCN=NFCN+1
C
10    XN=X+H
      FCM=DABS(H)/HMAX
      IF(FCM.LT.FMIN) FCM=FMIN
C
C  DISCRETIZATION
C-----------------------------------------------------------------------
      DO 260 J=1,JM
      M=NJ(J)
      G=H/DBLE(M)
      B=G+G
      DO 210 I=1,N
      YL(I)=Y(I)
210   YM(I)=Y(I)+G*DZ(I)
      M=M-1
C  EXPLICIT MID-POINT RULE
      DO 220 K=1,M
      CALL  F (X+G*DBLE(K),YM,DY,MYP)
      NFCN=NFCN+1
      DO 220 I=1,N
      U=YL(I)+B*DY(I)
      YL(I)=YM(I)
      YM(I)=U
220   CONTINUE
C  FINAL STEP
      CALL  F (XN,YM,DY,MYP)
      NFCN=NFCN+1
      DM=ZERO
      DO 2200 I=1,N
      YH=YL(I)+G*DY(I)
      DIFF=YH-YM(I)
      YM(I)=(YM(I)+YH)*HALF
      DIFF=DABS(DIFF)
      SK=DABS(YM(I))
      IF(SK.LT.S(I)) SK=S(I)
      DIFF=DIFF/SK
      IF(DIFF.GT.DM) DM=DIFF
2200  CONTINUE
C
C  STABILITY CHECK
      IF(DM.LT.DMAX) GOTO 2209
C
C  EMERGENCY EXIT
      IF(KFLAG.GT.0) WRITE(LOUT,1006)
      GOTO 2601
C
C  PREVENTION OF POSSIBLE ORDER INCREASE
2209  IF(J.GT.2.OR.DM.LT.DMAX*HALF) GOTO 2207
      DO 2208 L=JOH,JM
      IF(INCR(L).GT.0) INCR(L)=0
      INCR(L)=INCR(L)-2
2208  CONTINUE
C
C  EXTRAPOLATION
C  ( IPOL=1: POLYNOMIAL, IPOL=0: RATIONAL)
2207  ERR=ZERO
      DO 234 I=1,N
      V=DT(I,1)
      C=YM(I)
      DT(I,1)=C
      IF(J.EQ.1) GOTO 234
      TA=C
      DO 231 K=2,J
      JK=J-K+1
      B1=D(J,JK)
      W=C-V
      IF(IPOL.EQ.0) GOTO 229
      U=W/(B1-ONE)
      C=B1*U
      GOTO 230
229   B1=B1*V
      B=B1-C
      U=V
      IF(B.EQ.ZERO) GOTO 230
      B=W/B
      U=C*B
      C=B1*B
230   V=DT(I,K)
      DT(I,K)=U
231   TA=U+TA
      YM(I)=TA
      TA=DABS(TA)
      IF(TA.LT.S(I)) TA=S(I)
      U=U/TA
      ERR=ERR+U*U
234   CONTINUE
      IF(J.EQ.1) GOTO 260
C ERROR (SCALED ROOT MEAN SQUARE)
      ERR=DSQRT(ERR/FN)
      KONV=0
      IF(ERR.LT.EPS) KONV=1
      ERR=ERR/EPH
C
C ORDER CONTROL
      K=J-1
      L=J+K
      FC=ERR**(ONE/DBLE(L))
      IF(FC.LT.FCM) FC=FCM
C  OPTIMAL ORDER DETERMINATION
      OMJ=FC*A(J)
      IF(J.GT.2.AND.OMJ*ONE1.GT.OMJO.OR.K.GT.JOH) GOTO 235
      KO=K
      JO=J
      OMJO=OMJ
      FCO=FC
235   IF(J.LT.KOH.AND.NSTEP.GT.0) GOTO 260
      IF(KONV.EQ.0) GOTO 236
      IF(KO.LT.K.OR.INCR(J).LT.0) GOTO 20
C  POSSIBLE INCREASE OF ORDER
      IF(NRED(KO).GT.0) NRED(KO)=NRED(KO)-1
      FC=FCO/AL(J,K)
      IF(FC.LT.FCM) FC=FCM
      J1=J+1
      IF(A(J1)*FC*ONE1.GT.OMJO) GOTO 20
      FCO=FC
      KO=JO
      JO=JO+1
      GOTO 20
C
C
C  CONVERGENCE MONITOR
236   RED=ONE/FCO
      JK=KM
      IF(JOH.LT.KM) JK=JOH
      IF(K.GE.JK) GOTO 239
      IF(KO.LT.KOH) RED=AL(KOH,KO)/FCO
237   IF(AL(JK,KO).LT.FCO) GOTO 239
260   CONTINUE
C-----------------------------------------------------------------------
C
C STEPSIZE REDUCTION (DUE TO EXTRAPOLATION TABLE)
239   RED=RED*SAFE
      H=H*RED
2392  IF(NSTEP.EQ.0) GOTO 2390
      NRED(KOH)=NRED(KOH)+1
      DO 2391 L=KOH,KM
2391  INCR(L)=-2-NRED(KOH)
2390  JRED=JRED+1
      IF(KFLAG.GT.0) WRITE(LOUT,1002) JRED,RED,KOH
      IF(JRED.GT.JRMAX) GOTO 32
      GOTO 10
C
C  STEPSIZE REDUCTION (DUE TO STABILITY)
2601  HMAX=G*FJ1*QUART
      RED=HMAX/DABS(H)
      H=HMAX
      IF(JRED.GT.0) GOTO 2390
      GOTO 2392
C
C  PREPARATIONS FOR NEXT BASIC INTEGRATION STEP
C-----------------------------------------------
20    X=XN
      IF(NSTEP.EQ.0) HS1=H
      IF(NSTEP.EQ.1) HS2=H
      H1=XEND-X
      DO  2606 I=1,N
2606  Y(I)=YM(I)
      NSTEP=NSTEP+1
      IF(NSTEP.GT.NSTMAX) GO TO 31
C
C STEPSIZE PREDICTION
      IF(FCO.NE.FCM) HR=H
      H=H/FCO
      KOH=KO
      JOH=KOH+1
      IF(DABS(H).GT.DABS(X)*EPMACH) GO TO 401
      GO TO 33
C
C  SOLUTION EXIT
C----------------
403   H=HR
      HS=(HS1+HS2)*HALF
      IPAR(6)=KFLAG
      IPAR(7)=NSTEP
      IPAR(8)=NFCN
      IF(KFLAG.GT.1) WRITE(LOUT,1009) NSTEP,NFCN,X,K,KOH
      IF(KFLAG.GT.1) WRITE(LOUT,1000) NSTEP,NFCN,X,(Y(I),I=1,N)
      HMAX=HMAXU
      RETURN
C
C  FAIL EXIT
C------------
31    IF(KFLAG.GE.1)WRITE(LOUT,1008) NSTMAX
      IPAR(6)=-1
      GOTO 39
32    IF(KFLAG.GE.1) WRITE(LOUT,1010)JRMAX
      IPAR(6)=-2
      GOTO 39
33    IF(KFLAG.GE.1) WRITE(LOUT,1004)
      IPAR(6)=-3
39    H=ZERO
      HMAX=HMAXU
      RETURN
C
1000  FORMAT(1H ,2I9,5D20.11,/,(1H ,38X,4D20.11))
1001  FORMAT(1H0,24HDIFEX     REL.PREC. EPS ,D10.3,8HMAX.COL.,I3
     *        ,//,5X,4HSTEP,3X,7HF-CALLS,8X,1HX,19X,7HY1(X)..,//)
1002  FORMAT(1H ,I3,17HREDUCTION FACTOR ,D10.3,I9,/)
1004  FORMAT(/,40H0  STEPSIZE REDUCTION FAILED TO SUCCEED, //)
1006  FORMAT(/,46H0  STABILITY CHECK ACTIVATED  STEPSIZE REDUCED,/)
1008  FORMAT(18H0MORE THAN NSTMAX=,I6,18H INTEGRATION STEPS,//)
1009  FORMAT(1H ,2I9,D20.11,I9,I6,/)
1010  FORMAT(17H0MORE THAN JRMAX=,I3,29H STEPSIZE REDUCTIONS PER STEP,/)
C
C
C
C  END PRDIFX
C
      END
C
C
      SUBROUTINE PRBALA(NM,N,A,LOW,IGH,SCALE,MYP)
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      DOUBLE PRECISION A(NM,N),SCALE(N)
      DOUBLE PRECISION C,F,G,R,S,B2,RADIX
      DOUBLE PRECISION DABS,MYP(*)
      LOGICAL NOCONV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT:
C
C        A CONTAINS THE BALANCED MATRIX;
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N;
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     PRBALA  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.
C                RADIX = 16.0D0 FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C     DATA RADIX/Z4210000000000000/
      DATA RADIX/1.6D+1/
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     :::::::::: IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ::::::::::
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     :::::::::: SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ::::::::::
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     :::::::::: FOR J=L STEP -1 UNTIL 1 DO -- ::::::::::
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0D0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     :::::::::: SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ::::::::::
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0D0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     :::::::::: NOW BALANCE THE SUBMATRIX IN ROWS K TO L ::::::::::
      DO 180 I = K, L
  180 SCALE(I) = 1.0D0
C     :::::::::: ITERATIVE LOOP FOR NORM REDUCTION ::::::::::
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0D0
         R = 0.0D0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + DABS(A(J,I))
            R = R + DABS(A(I,J))
  200    CONTINUE
C     :::::::::: GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ::::::::::
         IF (C .EQ. 0.0D0 .OR. R .EQ. 0.0D0) GO TO 270
         G = R / RADIX
         F = 1.0D0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     :::::::::: NOW BALANCE ::::::::::
  240    IF ((C + R) / F .GE. 0.95D0 * S) GO TO 270
         G = 1.0D0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
C
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
C     :::::::::: LAST CARD OF PRBALA ::::::::::
      END
C
C
      SUBROUTINE PRHQR(NM,N,LOW,IGH,H,WR,WI,IERR,MYP)
C
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITS,LOW,MP2,ENM2,IERR
      DOUBLE PRECISION H(NM,N),WR(N),WI(N)
      DOUBLE PRECISION P,Q,R,S,T,W,X,Y,ZZ,NORM,MACHEP
      DOUBLE PRECISION DSQRT,DABS,DSIGN,MYP(*)
      INTEGER MIN0
      LOGICAL NOTLAS
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  PRBALA.  IF  PRBALA  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N;
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG
C          FORM BY  ELMHES  OR  PRORTH, IF PERFORMED, IS STORED
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C     ON OUTPUT:
C
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED
C          BEFORE CALLING  PRHQR  IF SUBSEQUENT CALCULATION AND
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED;
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C     DATA MACHEP/Z3410000000000000/
      DATA MACHEP/2.22D-16/
C
      IERR = 0
      NORM = 0.0D0
      K = 1
C     :::::::::: STORE ROOTS ISOLATED BY PRBALA
C                AND COMPUTE MATRIX NORM ::::::::::
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + DABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0D0
   50 CONTINUE
C
      EN = IGH
      T = 0.0D0
C     :::::::::: SEARCH FOR NEXT EIGENVALUES ::::::::::
   60 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     :::::::::: LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ::::::::::
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = DABS(H(L-1,L-1)) + DABS(H(L,L))
         IF (S .EQ. 0.0D0) S = NORM
         IF (DABS(H(L,L-1)) .LE. MACHEP * S) GO TO 100
   80 CONTINUE
C     :::::::::: FORM SHIFT ::::::::::
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     :::::::::: FORM EXCEPTIONAL SHIFT ::::::::::
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
      X = 0.75D0 * S
      Y = X
      W = -0.4375D0 * S * S
  130 ITS = ITS + 1
C     :::::::::: LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ::::::::::
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         IF (DABS(H(M,M-1)) * (DABS(Q) + DABS(R)) .LE. MACHEP * DABS(P)
     X    * (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0D0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0D0
  160 CONTINUE
C     :::::::::: DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ::::::::::
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0D0
         IF (NOTLAS) R = H(K+2,K-1)
         X = DABS(P) + DABS(Q) + DABS(R)
         IF (X .EQ. 0.0D0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     :::::::::: ROW MODIFICATION ::::::::::
         DO 210 J = K, EN
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C
         J = MIN0(EN,K+3)
C     :::::::::: COLUMN MODIFICATION ::::::::::
         DO 230 I = L, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     :::::::::: ONE ROOT FOUND ::::::::::
  270 WR(EN) = X + T
      WI(EN) = 0.0D0
      EN = NA
      GO TO 60
C     :::::::::: TWO ROOTS FOUND ::::::::::
  280 P = (Y - X) / 2.0D0
      Q = P * P + W
      ZZ = DSQRT(DABS(Q))
      X = X + T
      IF (Q .LT. 0.0D0) GO TO 320
C     :::::::::: REAL PAIR ::::::::::
      ZZ = P + DSIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0D0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0D0
      WI(EN) = 0.0D0
      GO TO 330
C     :::::::::: COMPLEX PAIR ::::::::::
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = EN
 1001 RETURN
C     :::::::::: LAST CARD OF PRHQR ::::::::::
      END
C
C
      SUBROUTINE PRORTH(NM,N,LOW,IGH,A,ORT,MYP)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      DOUBLE PRECISION A(NM,N),ORT(IGH)
      DOUBLE PRECISION F,G,H,SCALE
      DOUBLE PRECISION DSQRT,DABS,DSIGN,MYP(*)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE PRORTH,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  PRBALA.  IF  PRBALA  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N;
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT:
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX;
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0D0
         ORT(M) = 0.0D0
         SCALE = 0.0D0
C     :::::::::: SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO 90 I = M, IGH
   90    SCALE = SCALE + DABS(A(I,M-1))
C
         IF (SCALE .EQ. 0.0D0) GO TO 180
         MP = M + IGH
C     :::::::::: FOR I=IGH STEP -1 UNTIL M DO -- ::::::::::
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
C
         G = -DSIGN(DSQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
C     :::::::::: FORM (I-(U*UT)/H) * A ::::::::::
         DO 130 J = M, N
            F = 0.0D0
C     :::::::::: FOR I=IGH STEP -1 UNTIL M DO -- ::::::::::
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
C
            F = F / H
C
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
C
  130    CONTINUE
C     :::::::::: FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ::::::::::
         DO 160 I = 1, IGH
            F = 0.0D0
C     :::::::::: FOR J=IGH STEP -1 UNTIL M DO -- ::::::::::
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
C
            F = F / H
C
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
C
  160    CONTINUE
C
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
C
  200 RETURN
C     :::::::::: LAST CARD OF PRORTH ::::::::::
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C procedure CALLPERIOD calculates periodic solution
C and its derivatives
C
C input: NN - number of diferential equations
C        M - number of nodes
C        IPR - output (print) settings (-1 - no output, 0 - standard, 1 - detailed)
C        X0 - guess for the intial point X(0)
C        P - guess for the period
C        NPAR - number of parameters
C        TRPAR - array of parameter values
C
C output: XT - values in nodes for scaled time (t/P)
C         Y - periodic solution in nodes
C         FM - eigenvalues of the Poincare map
C         HES - jacobian of the Poincare map (Fx)
C         POUT - period for the periodic solution
C         IFAIL - integer value which indicates if calculation of PERIOD is succeeded (TRUE if IFAIL>0)
C         ERRY - value which indicates if solution of PERIOD is constant (TRUE if ERRY>10^(-7))
C         FAL - matrix of derivarives with respect to parameters (Fp)
C         FAL - tensor of derivarives with respect to variables and parameters (Fxp)
C         FXX - tensor of second derivatives with respect to variables (Fxx)
C         FX - jacobian of the Poincare map calculated with finite differences
C         HDIF - step size for finite difference method
C
C
C revision:
C 2011-01-26 written by dka
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CALLPERIOD(NN,M,IPR,X0,P,NPAR,TRPAR,
     1             XT,Y,FM,HES,POUT,IFAIL,ERRY,FAL,
     2             FXAL,FXX,FX,H2)

      EXTERNAL PRDIFX,FCN,DFDY,PERIOD
      INTEGER NN,M,NPAR
      DOUBLE PRECISION ZZ(NN),Y(NN,M),XT(M),RW(6000),
     1          USCAL(NN),FM(NN,2),IPAR(10),RH(2000),
     2          IH(100),X0(NN)
      INTEGER IPR,M1,NRW,NIW,NRH,NIH,ICALL(8)
      INTEGER IW(100),I,I2,SIZEM
      DOUBLE PRECISION P,ERRY,POUT,ERRY1,ERRY2
      DOUBLE PRECISION TRPAR(NPAR),HES(NN,NN)
      DOUBLE PRECISION RED,H,HS,TOL,T,TEND,HMAX,EPS
      DOUBLE PRECISION FAL(NN,NPAR),FXAL(NN,NN,NPAR)
      DOUBLE PRECISION FXX(NN,NN,NN),FX(NN,NN)
      DOUBLE PRECISION HD1,XST(NN),EPSER,HDIF,H2

      HDIF=1.D-7
C
      M1=M-1

C
C  AMOUNT OF REAL- AND INTEGER-WORKSPACE
      NRW=6000
      NIW=100
C
C  SHOOTING NODES
      XT(1)=0.D0
      RED=1.D0/DBLE(M1)

      DO 10 I=2,M
10    XT(I)=XT(I-1)+RED

C
C  INITIAL VALUES OF TRAJECTORY

C
C  CALLING DIFEXP TO GENERATE INITIAL TRAJECTORY
      DO 2 I=1,NN
2     Y(I, 1) = X0(I)

      NRH=2000
      NIH=100
      IPAR(1)=0
      IPAR(2)=0
      IPAR(3)=0
      IPAR(4)=NN
      IPAR(5)=0
      IPAR(6)=0
      H=P/M1

      DO 3 I=1,NN
3     ZZ(I) = Y(I,1)

      TOL=1.D-7
      T=0

      DO 18 I=1,100
      IW(I)=0
18    CONTINUE

      DO 22 I=2,M
      TEND=XT(I)*P
      HMAX=TEND-T
      CALL PRDIFX (NN,FCN,T,ZZ,TEND,TOL,HMAX,H,HS,USCAL,NRH,RH,NIH,IH,
     2             IPAR,TRPAR)

      DO 4 I2=1,NN
4     Y(I2,I)=ZZ(I2)

22    CONTINUE

C  DESIRED RELATIVE ACCURACY FOR SOLUTION
      EPS=1.D-7
C  CLASSIFICATION OF RIGHT-HAND SIDE (AUTONOMOUS SYSTEM)
      ICALL(1)=1
C  PROBLEM IS HIGHLY NONLINEAR
      ICALL(2)=2
C  MAXIMUM PERMITTED NEWTON ITERATIONS
      ICALL(3)=200
C  SOLVING VARIATIONAL EQUATION
      ICALL(4)=2
C  RANK-1 UPDATES ALLOWED
      ICALL(5)=1
C  ITERATIVE REFINEMENT IS ACTIVATED
      ICALL(6)=1
C  THE SYSTEM IS NONSTIFF
      ICALL(7)=0
C  PRINT PARAMETER
      ICALL(8)=IPR

      CALL PERIOD(PRDIFX,NN,M,XT,Y,P,EPS,FM,6000,RW,100,IW,
     1            ICALL,TRPAR,HES)

      IFAIL=ICALL(8)
      POUT=P
      ERRY1=DABS(Y(1,1)-Y(1,2))+DABS(Y(1,1)-Y(1,3))
      ERRY2=DABS(Y(1,1)-Y(1,4))+DABS(Y(1,1)-Y(1,5))
      ERRY=ERRY1+ERRY2

      EPSER=1.D-7

C     IF PERIODIC SOLUTION IS FOUND
C     THEN CALCULATE DERIVATIVES
c      IF ((IFAIL.GT.0).AND.(ERRY.GT.EPSER)) THEN
      IF (IFAIL.GT.0) THEN

C     STEP SIZE H FOR DERIVATIVES CALCULATION
      HD1=HDIF

C     POINT X(0) ON THE POINCARE SECTION
      DO 5 I2=1,NN
5     XST(I2)=Y(I2,1)

C     SIZE OF OUTPUT ROW IF TILL 2ND ORDER DIRAVATIVAS WAS CALCULATED
      SIZEM=1+NN+(NN+NPAR)*NN+((NN+NPAR)*(NN+NPAR+1)*NN)/2

      CALL FTOTIDES(NN,NPAR,TRPAR,XST,POUT,SIZEM,
     1              FX,FAL,FXX,FXAL)



      ENDIF


      END

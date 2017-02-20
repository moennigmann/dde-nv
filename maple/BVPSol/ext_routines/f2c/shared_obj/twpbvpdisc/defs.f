CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The procedure defines the initial guess for solution of BVP.
C It uses the external procedure from PERIOD
C
C Parameters: NN - number of system variables
C             M -  number of nodes for the guess
C             X0 - initial point,
C             P - guess for period,
C             NPAR - numver of parameters,
C             TRPAR - parameters values,
C             XT - guess for x nodes (output),
C             Y - guess for y nodes (output)
C
C Revision: written by dka 2012-01-16
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
      
      SUBROUTINE DEFS(NN,M,X0,P,NPAR,TRPAR,
     1             XT,Y)

      EXTERNAL PRDIFX,FCN,DFDY
      INTEGER NN,M,NPAR
      DOUBLE PRECISION ZZ(NN),Y(NN,M),XT(M),RW(6000),
     1          USCAL(NN),IPAR(10),RH(6000),
     2          IH(100),X0(NN)
      INTEGER IPR,M1,NRW,NIW,NRH,NIH,ICALL(8)
      INTEGER IW(100),I,I2
      DOUBLE PRECISION P
      DOUBLE PRECISION TRPAR(NPAR)
      DOUBLE PRECISION RED,H,HS,TOL,T,TEND,HMAX

C  PRESICION OD A SOLUTION
      HDIF=1.D-7
C  NUMBER OF INTERVALS FOR SOLUTION GUESS
      M1=M-1

C
C  AMOUNT OF REAL- AND INTEGER-WORKSPACE
      NRW=6000
      NIW=100
C
C  NODES WHERE INITIAL VALUES WILL BE CALCULATED
      XT(1)=0.D0
      RED=1.D0/DBLE(M1)

      DO 10 I=2,M
10    XT(I)=XT(I-1)+RED


C  CALLING DIFEXP TO GENERATE INITIAL TRAJECTORY
      DO 2 I=1,NN
2     Y(I, 1) = X0(I)

C  PREPAIRING PARAMETERS SETTINGS FOR CALLING 'DIFEXP'
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

C  FOR EACH INTERVAL CALLING 'DIFEXP'
      DO 22 I=2,M
      TEND=XT(I)*P
      HMAX=TEND-T
      CALL PRDIFX (NN,FCN,T,ZZ,TEND,TOL,HMAX,H,HS,USCAL,NRH,RH,NIH,IH,
     2             IPAR,TRPAR)

      DO 4 I2=1,NN
4     Y(I2,I)=ZZ(I2)

22    CONTINUE

      RETURN
      END 

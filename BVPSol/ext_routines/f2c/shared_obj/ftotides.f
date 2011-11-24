      SUBROUTINE FTOTIDES(NN,NPAR,PAR,XST,P,SIZEM,
     1                      FX,FP,FXX,FXP)
      INTEGER NN,NN1,NPAR,NPAR1,NGEN
      DOUBLE PRECISION PAR(NPAR),P,XST(NN)
      INTEGER SIZEM, I1,J1, NOFF, I2,J2,K1
      DOUBLE PRECISION OUTM(SIZEM)
      DOUBLE PRECISION FP(NN,NPAR),FXP(NN,NN,NPAR)
      DOUBLE PRECISION FXX(NN,NN,NN),FX(NN,NN)
      DOUBLE PRECISION FPP(NN,NPAR,NPAR)
      DOUBLE PRECISION COMR,EPS

      NN1=NN
      NPAR1=NPAR
      call calltides(NN1,NPAR1,PAR,XST,P,OUTM)

c     PROVE IF INTEGRATION WITH TIDES RETURNS PERIODIC ORBIT
      EPS=1.D-8

      DO 2 I1=1,NN1
      COMR=ABS(OUTM(I1+1)-XST(I1))

      IF (COMR.GT.EPS) THEN
      PRINT *,'Integration with TIDES failed'
      PRINT *, COMR
      ENDIF

2     CONTINUE

C     EXTRACT FROM OUTM MATRIX FX
      NOFF=1+NN1
      I2=0

      DO 3 I1=1,NN1
      DO 4 J1=1,NN1
      I2=I2+1
      FX(J1,I1)=OUTM(NOFF+I2)
4     CONTINUE
3     CONTINUE

C     EXTRACT FROM OUTM MATRIX FP
      NOFF=NOFF+I2
      I2=0

      DO 5 I1=1,NPAR1
      DO 6 J1=1,NN1
      I2=I2+1
      FP(J1,I1)=OUTM(NOFF+I2)
6     CONTINUE
5     CONTINUE

C     EXTRACT FROM OUTM MATRIX FXX,FXP,FPP
      NOFF=NOFF+I2
      I2=0
      NGEN=NPAR1+NN1

      DO 7 I1=1,NGEN
      DO 8 J1=I1,NGEN
      DO 9 K1=1,NN1
      I2=I2+1
      IF ((I1.LE.NN1).AND.(J1.LE.NN1)) THEN
      FXX(K1,I1,J1)=OUTM(NOFF+I2)
      FXX(K1,J1,I1)=OUTM(NOFF+I2)
      END IF

      IF ((I1.LE.NN1).AND.(J1.GT.NN1)) THEN
      FXP(K1,I1,J1-NN1)=OUTM(NOFF+I2)
      END IF

      IF ((I1.GT.NN1).AND.(J1.GT.NN1)) THEN
      FPP(K1,I1-NN1,J1-NN1)=OUTM(NOFF+I2)
      FPP(K1,J1-NN1,I1-NN1)=OUTM(NOFF+I2)
      END IF

9     CONTINUE
8     CONTINUE
7     CONTINUE

      END
      SUBROUTINE PEDECC(A,NROW,NCOL,MCON,M,N,IRANK,COND,D,PIVOT,
     *KRED,AH,V)
C*    Begin Prologue DECCON
      INTEGER IRANK,MCON
      INTEGER M,N,NROW,NCOL,KRED
      INTEGER PIVOT(NCOL)
      DOUBLE PRECISION COND
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION D(NCOL),V(NCOL)
C     ------------------------------------------------------------
C
C*  Title
C
C*    Deccon - Constrained Least Squares QR-Decomposition
C
C*  Written by        P. Deuflhard,  L.Weimann
C*  Purpose           Solution of least squares problems, optionally
C                     with equality constraints.
C*  Method            Constrained Least Squares QR-Decomposition
C                     (see references below)
C*  Category          D9b1. -  Singular, overdetermined or
C                              underdetermined systems of linear 
C                              equations, generalized inverses. 
C                              Constrained Least Squares solution
C*  Keywords          Linear Least Square Problems, constrained, 
C                     QR-decomposition, pseudo inverse.
C*  Version           1.0
C*  Revision          August 1987
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
C     ===========
C
C       /1/ P.Deuflhard, V.Apostolescu:
C           An underrelaxed Gauss-Newton method for equality
C           constrained nonlinear least squares problems.
C           Lecture Notes Control Inform. Sci. vol. 7, p.
C           22-32 (1978)
C       /2/ P.Deuflhard, W.Sautter:
C           On rank-deficient pseudoinverses.
C           J. Lin. Alg. Appl. vol. 29, p. 91-111 (1980)
C    
C*    Related Programs:     SOLCON (PESOLC)
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
C    This code is under partial care of ZIB and belongs to ZIB software
C    class 2.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Constrained QR-decomposition of (M,N)-system  with
C     computation of pseudoinverse in case of rank-defeciency .
C     First MCON rows belong to equality constraints.
C
C     ------------------------------------------------------------
C
C*    Parameters:
C     ===========
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       A(NROW,NCOL)          Input matrix
C                             A(M,N)contains actual input
C       NROW                  declared number of rows of A and AH
C       NCOL                  declared number of columns of A and
C                             AH
C     * MCON                  number of equality constraints
C                             MCON.LE.N
C                             internally reduced if equality
C                             constraints
C                             are linearly dependent
C       M                     treated number of rows of matrix A
C       N                     treated number of columns of matrix A
C     * IRANK                 pseudo-rank of matrix A
C     * COND                  permitted upper bound of DABS(D(1)/D
C                             (IRANK))and of DABS(D(IRANK+1)/D(
C                             IRANK))
C                             ( sub - condition numbers of A )
C       KRED                  Type of operation
C                             >= 0  Householder triangularization
C                                   ( build up of pseudo -
C                                   inverse,if IRANK.LT.N )
C                             <  0  reduction of pseudo-rank of
C                                   matrix A , skipping Householder
C                                   triangularization, build-up of
C                                   new pseudo-inverse
C       V(N)                  real work array
C
C*    Output parameters
C     =================
C
C       A(M,N)                Output matrix updating product of
C                             Householder transformations and
C                             upper triangular matrix
C       MCON                  Pseudo-rank of constrained part of
C                             matrix A
C       IRANK                 Pseudo-rank of total matrix A
C       D(IRANK)              Diagonal elements of upper
C                             triangular matrix
C       PIVOT(N)              Index vector storing permutation of
C                             columns due to pivoting
C       COND                  Sub-condition number of A
C                             (in case of rank reduction:
C                             sub-condition number which led to
C                             rank reduction)
C       AH(N,N)               Updating matrix for part of pseudo
C                             inverse
C
C*    Subroutines called: none
C
C*    Machine constants used
C     ======================
C
C     EPMACH = relative machine precision
C     DOUBLE PRECISION EPMACH
C
C     ____________________________________________________________
C*    End Prologue
      INTRINSIC DABS,DSQRT
      DOUBLE PRECISION EPMACH,SMALL
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      INTEGER L1
      DOUBLE PRECISION S1
      INTEGER I,II,IRK1,ISUB,I1,J,JD,JJ,K,K1,LEVEL,MH
      DOUBLE PRECISION DD,H,HMAX,S,SH,SMALLS,T
C*    Begin
C     ______________________________________________________________
C     1 Initialization
      CALL ZIBCONST(EPMACH,SMALL)
      SMALLS = DSQRT(EPMACH*10.0D0)
      IF(IRANK.GT.N) IRANK = N
      IF(IRANK.GT.M) IRANK = M
C     ______________________________________________________________
C     1.1 Special case M=1 and N=1
      IF(M.EQ.1.AND.N.EQ.1)THEN
        PIVOT(1)=1
        D(1)=A(1,1)
        COND = ONE
        RETURN
      ENDIF
      IF(KRED.GE.0)THEN
C       ____________________________________________________________
C       1.1 Initialize pivot-array
        DO 11 J=1,N
          PIVOT(J)=J
11      CONTINUE
C       ____________________________________________________________
C       2. Constrained Householder triangularization
        JD = 1
        ISUB = 1
        MH = MCON
        IF(MH.EQ.0) MH = M
        K1 = 1
        LEVEL = 1
C       DO (Until)
2       CONTINUE
          K = K1
          IF(K.NE.N)THEN
            K1 = K+1
C           DO (Until)
20          CONTINUE
              IF(JD.NE.0)THEN
                DO 201 J=K,N
                  S = ZERO
                  DO 2011 L1=K,MH
                    S = S+A(L1,J)**2
2011              CONTINUE
                  D(J)=S
201             CONTINUE
              ENDIF
C             ______________________________________________________
C             2.1 Column pivoting
              S1 = D(K)
              JJ = K
              DO 21 L1=K,N
                IF(D(L1).GT.S1) THEN
                  S1=D(L1)
                  JJ = L1
                ENDIF
21            CONTINUE
              H = D(JJ)
              IF(JD.EQ.1) HMAX = H*SMALLS
              JD = 0
              IF(H.LT.HMAX) JD = 1
            IF(.NOT.(H.GE.HMAX)) GOTO 20
C           UNTIL ( expression - negated above)
            IF(JJ.NE.K)THEN
C             ______________________________________________________
C             2.2 Column interchange
              I = PIVOT(K)
              PIVOT(K)=PIVOT(JJ)
              PIVOT(JJ)=I
              D(JJ)=D(K)
              DO 221 L1=1,M
                S1=A(L1,JJ)
                A(L1,JJ)=A(L1,K)
                A(L1,K)=S1
221           CONTINUE
            ENDIF
          ENDIF
          H = ZERO
          DO 222 L1=K,MH
            H = H+A(L1,K)**2
222       CONTINUE
          T = DSQRT(H)
C         __________________________________________________________
C         2.3.0 A-priori test on pseudo-rank
          IF(ISUB.GT.0) DD = T/COND
          ISUB = 0
          IF(T.GT.DD.OR.K.LE.MCON)THEN
            IF(T.LE.DD)THEN
C             ______________________________________________________
C             2.3.1 Rank reduction
              MCON = K-1
              K1 = K
              MH = M
              JD = 1
              ISUB = 1
            ELSE
C             ______________________________________________________
C             2.4 Householder step
              S = A(K,K)
              IF(S.GT.ZERO) T =-T
              D(K)=T
              A(K,K)=S-T
              IF(K.NE.N)THEN
                T = ONE/(H-S*T)
                DO 24 J=K1,N
                  S = ZERO
                  DO 241 L1=K,MH
                    S = S+A(L1,K)*A(L1,J)
241               CONTINUE
                  S = S*T
                  S1 =-S
                  DO 242 L1=K,M
                    A(L1,J) = A(L1,J)+A(L1,K)*S1
242               CONTINUE
                  D(J)=D(J)-A(K,J)**2
24              CONTINUE
                IF(K.NE.IRANK)THEN
                  IF(K.EQ.MCON)THEN
                    MH = M
                    JD = 1
                    ISUB = 1
                  ENDIF
                ELSE
                  LEVEL = 3
                ENDIF
              ELSE
                LEVEL = 4
              ENDIF
            ENDIF
          ELSE
            IRANK = K-1
            IF(IRANK.EQ.0)THEN
              LEVEL = 4
            ELSE
              LEVEL = 3
            ENDIF
          ENDIF
        IF(.NOT.(LEVEL.NE.1)) GOTO  2
C       UNTIL ( expression - negated above)
      ELSE
        LEVEL = 3
      ENDIF
      IF(LEVEL.EQ.3)THEN
C     ______________________________________________________________
C     3 Rank-deficient pseudo-inverse
        IRK1 = IRANK+1
        DO 3 J=IRK1,N
          DO 31 II=1,IRANK
            I = IRK1-II
            S = A(I,J)
            IF(II.NE.1)THEN
              SH = ZERO
              DO 311 L1=I1,IRANK
                SH=SH+A(I,L1)*V(L1)
311           CONTINUE
              S = S-SH
            ENDIF
            I1 = I
            V(I)=S/D(I)
            AH(I,J)=V(I)
31        CONTINUE
          DO 32 I=IRK1,J
            S = ZERO
            DO 321 L1=1,I-1
              S = S+AH(L1,I)*V(L1)
321         CONTINUE
            IF(I.NE.J)THEN
              V(I)=-S/D(I)
              AH(I,J)=-V(I)
            ENDIF
32        CONTINUE
          D(J)=DSQRT(S+ONE)
3       CONTINUE
      ENDIF
C    ______________________________________________________________
C     9 Exit
      IF(K.EQ.IRANK)THEN
        T = D(IRANK)
        IF(T.NE.ZERO) COND = DABS(D(1)/T)
      ENDIF
      RETURN
      END
      SUBROUTINE PESOLC(A,NROW,NCOL,MCON,M,N,X,B,IRANK,D,PIVOT,
     *KRED,AH,V)
C*    Begin Prologue SOLCON
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION X(NCOL),B(NROW),D(NCOL),V(NCOL)
      INTEGER NROW,NCOL,MCON,M,N,IRANK,KRED
      INTEGER PIVOT(NCOL)
C     ____________________________________________________________
C
C*    Summary
C     =======
C
C     Best constrained linear least squares solution of (M,N)-
C     system . First MCON rows comprise MCON equality constraints.
C     To be used in connection with subroutine DECCON (PEDECC)
C     References:       See DECCON
C     Related Programs: DECCON
C
C*    Parameters:
C     ===========
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       A(M,N)            see output of DECCON
C       NROW              see output of DECCON
C       NCOL              see output of DECCON
C       M                 see output of DECCON
C       N                 see output of DECCON
C       MCON              see output of DECCON
C       IRANK             see output of DECCON
C       D(N)              see output of DECCON
C       PIVOT(N)          see output of DECCON
C       AH(N,N)           see output of DECCON
C       KRED              see output of DECCON
C     * B(M)              right-hand side of linear system, if
C                         KRED.GE.0
C                         right-hand side of upper linear system,
C                         if KRED.LT.0
C       V(N)              real work array
C
C*    Output parameters
C     =================
C
C       X(N)              best lsq-solution of linear system
C       B(M)              right-hand of upper trigular system
C                         ( transformed right-hand side of linear
C                         system )
C
C     ____________________________________________________________
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,L2
      DOUBLE PRECISION S1
      INTEGER I,II,I1,IRK1,J,JJ,J1,MH
      DOUBLE PRECISION S,SH
C*    Begin
C     ____________________________________________________________
C     1 Solution for pseudo-rank zero
      IF(IRANK.EQ.0)THEN
        S1 = ZERO
        DO 1 L1=1,N
          X(L1)=S1
1       CONTINUE
        RETURN
      ENDIF
      IF(KRED.GE.0.AND.(M.NE.1.OR.N.NE.1))THEN
C       __________________________________________________________
C       2 Constrained householder transformations of right-hand side
        MH = MCON
        IF(MH.EQ.0) MH = M
        DO 2 J=1,IRANK
          S = ZERO
          DO 21 L1=J,MH
            S = S+A(L1,J)*B(L1)
21        CONTINUE
          S = S/(D(J)*A(J,J))
          S1 = S
          DO 22 L1=J,M
            B(L1)=B(L1)+A(L1,J)*S1
22        CONTINUE
          IF(J.EQ.MCON) MH = M
2       CONTINUE
      ENDIF
C     ____________________________________________________________
C     3 Solution of upper triangular system
      IRK1 = IRANK+1
      DO 31 II=1,IRANK
        I = IRK1-II
        I1 = I+1
        S = B(I)
        IF(I1.LE.IRANK)THEN
          SH = ZERO
          DO 311 L1=I1,IRANK
            SH=SH+A(I,L1)*V(L1)
311       CONTINUE
          S = S-SH
        ENDIF
        V(I)=S/D(I)
31    CONTINUE
      IF(IRK1.LE.N)THEN
C       __________________________________________________________
C       3.2 Computation of the best constrained least squares-
C           solution
        DO 321 J=IRK1,N
          S = ZERO
          DO 3211 L1=1,J-1
            S = S+AH(L1,J)*V(L1)
3211      CONTINUE
          V(J)=-S/D(J)
321     CONTINUE
        DO 322 JJ=1,N
          J = N-JJ+1
          S = ZERO
          IF(JJ.NE.1)THEN
            DO 3221 L1=J1,N
              S=S+AH(J,L1)*V(L1)
3221        CONTINUE
          ENDIF
          IF(JJ.NE.1.AND.J.LE.IRANK)THEN
            V(J)=V(J)-S
          ELSE
            J1 = J
            V(J)=-(S+V(J))/D(J)
          ENDIF
322     CONTINUE
      ENDIF
C     ____________________________________________________________
C     4 Back-permutation of solution components
      DO 4 L1=1,N
        L2 = PIVOT(L1)
        X(L2) = V(L1)
4     CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C procedure CALLPERIOD calculates periodic solution
C and its derivatives with method TWPBVP
C
C input: NN - number of diferential equations
C        M - number of nodes
C        IPR - output (print) settings (-1 - no output, 0 - standard, 1 - detailed)
C        X0 - guess for the intial point X(0)
C        P - guess for the period
C        NPAR - number of parameters
C        TRPAR - array of parameter values
C        NPHASE - number of variable for which the phase condition is applied
C
C output: XT - values in nodes for scaled time (t/P)
C         Y - periodic solution in nodes in the case of calculation with PERIOD
C             or Y(1) starting condition in the case of calulation with TWPBVP
C         FM - eigenvalues of the Poincare map
C              (in the case of calulation with TWPBVP are not calculated)
C         POUT - period for the periodic solution
C         IFAIL - integer value which indicates if calculation of PERIOD is succeeded (TRUE if IFAIL=0)
C         ERRY - value which indicates if solution of PERIOD is constant (TRUE if ERRY>10^(-7))
C         FAL - matrix of derivarives with respect to parameters (Fp)
C         FAL - tensor of derivarives with respect to variables and parameters (Fxp)
C         FXX - tensor of second derivatives with respect to variables (Fxx)
C         FX - jacobian of the Poincare map
C         FPP - tensor of second derivatives with respect to variables (Fpp)
C         FXXP, FXPP - third order derivatives
C
C revision:
C 2012-01-07 written by dka
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CALLPERIOD(NN,M,IPR,X0,P,NPAR,TRPAR,
     1             Y,FM,POUT,IFAIL,ERRY,FAL,
     2             FXAL,FXX,FX,NPHASE,FPP,FXXP,FXPP)

      DOUBLE PRECISION P, PGUESS
      INTEGER NN, M, NUDIM, LTOL(2*NN), NLBC, NMSH, IPR, NPAR
      DOUBLE PRECISION X0(NN),XT(M),YU(NN,M),TRPAR(NPAR)
      DOUBLE PRECISION FM(NN,2),Y(NN,M),POUT,ERRY,ERRY1,ERRY2
      DOUBLE PRECISION FAL(NN,NPAR),FXAL(NN,NN,NPAR)
      DOUBLE PRECISION FXX(NN,NN,NN),FX(NN,NN)
      DOUBLE PRECISION FPP(NN,NPAR,NPAR)
      DOUBLE PRECISION FXXP(NN,NN,NN,NPAR)
      DOUBLE PRECISION FXPP(NN,NN,NPAR,NPAR)
      INTEGER NPHASE, offsetI, I2, SIZEM
      DOUBLE PRECISION XST(NN)

      DOUBLE PRECISION TOL(2*NN), FIXPNT(1)
      INTEGER LWRKFL, LWRKIN, NXXDIM, LISERIES
      PARAMETER(LWRKFL=2500000)
      PARAMETER(LWRKIN=500000)
      PARAMETER(NXXDIM=100000)
      PARAMETER (LISERIES=10000)
      DOUBLE PRECISION precis(3), xguess(NXXDIM), yguess(2*NN,NXXDIM)
      DOUBLE PRECISION U(2*NN,NXXDIM), XX(NXXDIM)
      INTEGER iseries(NXXDIM), NCOMP, NTOL, IWRK(LWRKIN), ind, I
      INTEGER NFXPNT, IPAR(1), NMAX
      DOUBLE PRECISION CKAPPA1, GAMMA1, CKAPPA, sigma, ckappa2
      LOGICAL   LINEAR, GIVEU, GIVMSH
      DOUBLE PRECISION ALEFT, ARIGHT, ETOL
      DOUBLE PRECISION WRK(LWRKFL)
      LOGICAL full, useC
      INTEGER IFAIL, IFLBVP, indnms, iset(6)

      EXTERNAL  TWPBVPC, FSUB, DFSUB, GSUB, DGSUB


C     TRPAR (real) and IPAR (integer) are arrays that are passed to
C     the F function (RHS of the proble), G function (bondary condition 
C     of the problem), and their Jacobians.
C     TRPAR is typically used as a parameter in the problem formulation, and
C     IPAR is normally used to select a problem from a set.

C     Dimension of the problem =2*n
C     n additional variable is substituted into BVP problem to consider period
C     as unknown variable and separate boudary conditions
      NUDIM=2*NN

C     If you do not like to use the conditioning in the mesh selection
C     useC  = .false.
      useC  =  .false.

C     WRK is the floating point workspace and IWRK
C     is the integer work space. We set these to zero here,
C     to ensure consistent behaviour with different compilers.
      DO ind=1,LWRKFL
         WRK(ind) = 0d0
      ENDDO
      DO ind=1,LWRKIN
         IWRK(ind) = 0
      ENDDO

C     ALEFT and ARIGHT are the values of x at the left
C     and right boundaries.
      ALEFT  = 0.D0
      ARIGHT = 1.D0

C     ETOL is the required error tolerance of the solution.
C     Decreasing ETOL will give a more accurate solution, but
C     more mesh points are likely to be required and the code
C     will take longer to finish.
C     ETOL is passed to this subroutine as an argument.

C     NTOL is the number of components that have to satisfy
C     the prescribed tolerance.
C     LTOL is the component that needs to be tested.
C     Most of the time one will set NTOL to the number of system components,
C     LTOL(i) to component i, and set the tolerance for each component TOL to be equal.
      ETOL=1.D-7
      NTOL = NUDIM
      TOL(1) = ETOL
      DO ind=1,ntol
        LTOL(ind)=ind
        TOL(ind) =TOL(1)
      ENDDO


C     logical variable full controls whether or not the solver prints out
C     verbose output.
      full=.TRUE.
      IF (IPR.LT.0) THEN
      full=.FALSE.
      ENDIF

C     The number of components in the system
      NCOMP = NUDIM

C     The number of boundary conditions at the left, the
C     number of the right being equal to NCOMP-NLBC
      NLBC=NN

C     NMSH is the number of initial points we supply to
C     the solver in the form of a guess, we set this
C     to zero to indicate that we have no initial guess
      NMSH  = M

C     Try with different guesses for periods
      PGUESS=P

C     Compute initial trajectory
      call DEFS(NN,M,X0,PGUESS,NPAR,TRPAR,XT,YU)

C     set initial guess
      DO 22 I=1,M
C     set for xguess equal intervals from 0 to 1 with step 1/(M-1),
C     which is output from procedure DEFS
      XX(I)=XT(I)
      xguess(I)=XT(I)

c     for first NN variables set values that are estimated with procedure DEFS
      DO 23 ind=1,NN
      U(ind,I)=YU(ind,I)
      yguess(ind,I)=YU(ind,I)
23    continue

C     for variable NN+1 set constant value of period guess
      U(NN+1,I)=PGUESS
      yguess(NN+1,I)=PGUESS

C     for the rest variables set constant values that are equal to initial data;
C     variable offsetI is introduced to omit variable which is involved into
C     the phase condition
      offsetI=0
      DO 24 ind=1,NN
      IF (ind.EQ.NPHASE) THEN
      offsetI=1
      ELSE
      U(ind+NN+1-offsetI,I)=YU(ind,1)
      yguess(ind+NN+1-offsetI,I)=YU(ind,1)
      ENDIF
24    continue
22    continue

C     values for the last node set equal to the values for the first node
      XX(M)=1.0D0
      xguess(M)=1.0D0

      DO 33 ind=1,NN
      U(ind,M)=YU(ind,1)
      yguess(ind,M)=YU(ind,1)
33    continue

C     fixed points are x values which are required by the
C     user to be part of the returned mesh.
C     NFXPNT is the number of fixed points, which we set to be zero
      NFXPNT= 0
      indnms=0

C     the problem is nonlinear so we specify .false. here
      LINEAR = .FALSE.

C     we do not supply any initial guesses, so again we
C     choose .false.
      GIVEU  = .TRUE.

C     No initial mesh (a set of x values) are given, so
C     this is .false. too
      GIVMSH = .TRUE.


      CALL TWPBVPC(NCOMP,NLBC,ALEFT,ARIGHT,NFXPNT,FIXPNT,NTOL,
     +            LTOL,TOL,LINEAR,GIVMSH,GIVEU,NMSH,NXXDIM,
     +            XX,NUDIM,U,NMAX,LWRKFL,WRK,LWRKIN,IWRK,
     +			  precis,
     +            FSUB,DFSUB,GSUB,DGSUB,
     +            ckappa1,gamma1,
     +            sigma,ckappa,
     +       ckappa2,TRPAR,ipar,iflbvp,liseries,iseries,indnms,
     +       full, useC, xguess, yguess, iset)


C     When returning from TWPBVPC, one should immediately
C     check IFLBVP, to see whether or not the solver was
C     succesful
c      IF (IFLBVP .LT. 0) THEN
c         WRITE(6,*) 'The code failed to converge!'
c         RETURN
c      END IF
      IFAIL=iflbvp

C     POUT IS CALCULATED PERIOD (which is (n+1)st variable)
      POUT=U(NN+1,1)

C     TO CHECK IF SOLUTION IS DIFFER FROM CONSTANT CALCULATING ERRY
      ERRY1=DABS(U(1,1)-U(1,2))+DABS(U(1,1)-U(1,3))
      ERRY2=DABS(U(1,1)-U(1,4))+DABS(U(1,1)-U(1,5))
      ERRY=ERRY1+ERRY2

C     IF PERIODIC SOLUTION IS FOUND AND IT IS NOT CONSTANT
C     THEN CALCULATE DERIVATIVES
      IF ((IFAIL.EQ.0).AND.(ERRY.GT.ETOL)) THEN

C     POINT X(0) ON THE POINCARE SECTION
      DO 5 I2=1,NN
5     XST(I2)=U(I2,1)

C     SIZE OF OUTPUT ROW IF TILL 3ND ORDER DIRAVATIVAS WAS CALCULATED
C    (BY CONSTRUCTION SUCH A ROW IS OUTPUT FROM TIDES)
      SIZEM=1+NN+(NN+NPAR)*NN
     1      +((NN+NPAR)*(NN+NPAR+1)*NN)/2
     2      +((NN+NPAR)*(NN+NPAR+1)*(NN+NPAR+2)*NN)/6

C     GET DERIVATIVES WITH TIDES
      CALL FTOTIDES(NN,NPAR,TRPAR,XST,POUT,SIZEM,
     1              FX,FAL,FXX,FXAL,FPP,FXXP,FXPP)

      ENDIF

C     MUNBER OF OUTPUT NODES
      M=NMSH
C     SOLUTION IN INITIAL POINT X(0)
      DO 6 I2=1,NN
6     Y(I2,1)=U(I2,1)

C     the solution x values are stored in XX, the Y
C     values are stored in U.
C     U(i,j) refers to component i of point j in the mesh.

      RETURN
      END

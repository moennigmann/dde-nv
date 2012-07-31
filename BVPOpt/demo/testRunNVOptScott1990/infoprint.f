CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C procedure printParConFun print information for critical parameters values
C           from constraint function template
C
C input: NPAR - number of parameters
C        PAR - array of parameter values
C
C revision:
C 2012-07-26 written by dka
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE printParConFun(NPAR,PAR)
      INTEGER NPAR
      DOUBLE PRECISION PAR(NPAR)
      INTEGER I1

      write(6,*) 'CONSTR PARS:', (PAR(I1),I1=1,NPAR)

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C procedure printParObjFun print information for parameters values
C           from objective function template
C
C input: NPAR - number of parameters
C        PAR - array of parameter values
C        OBJF - value of objective function
C
C revision:
C 2012-07-26 written by dka
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE printParObjFun(NPAR,PAR,OBJF)
      INTEGER NPAR
      DOUBLE PRECISION PAR(NPAR), OBJF
      INTEGER I1

      write(6,*) 'OBJ PARS:', (PAR(I1),I1=1,NPAR), 'OBJF:', OBJF 

      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C procedure printObjFun print objective function values
C
C input: OBJF - value of objective function
C
C revision:
C 2012-07-26 written by dka
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE printObjFun(OBJF)
      DOUBLE PRECISION OBJF

      write(6,*) 'OBJF:', OBJF 

      END
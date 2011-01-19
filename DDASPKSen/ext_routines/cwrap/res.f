C
C Copyright (C) 2002 by M. Monnigmann for 
C Lehrstuhl f. Prozesstechnik,
C RWTH Aachen. All rights reserved.
C
      SUBROUTINE RES (T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
      
      IMPLICIT NONE

      INTEGER NEQ
      PARAMETER (NEQ=2)

C
C     variables which are arguments of subroutine
C     note that cj, rpar, ipar, senpar are dummies for the simple example      
C
      INTEGER IRES, IPAR     
      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), CJ, DELTA(NEQ), RPAR,
     &                 SENPAR

C      
C     variables which are parameters of the ode system
C      
      DOUBLE PRECISION OMEGA, K

      K= 1.0
      OMEGA= 1.0

      DELTA(1)= YPRIME(1)+ K* Y(1)+ OMEGA* OMEGA* Y(2)
      DELTA(2)= YPRIME(2)- Y(1)
      
      END
      

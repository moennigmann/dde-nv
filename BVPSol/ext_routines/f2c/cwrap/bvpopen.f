      SUBROUTINE BVPOPEN
      
      INTEGER NOUT
      
C----------------------------------------------------------------------
C
C open file for PIRIOD output
C
C----------------------------------------------------------------------
      
      NOUT= 6
      
      OPEN (NOUT, FILE = './bvp_output.txt', 
     1      ACCESS = 'APPEND',STATUS = 'NEW')

      END 

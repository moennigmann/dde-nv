      SUBROUTINE OPENOPT
      IMPLICIT NONE
      INTEGER IOPTNS, INFORM, NOUT
      
C----------------------------------------------------------------------
C
C open options file and call NPFILE to set options for future
C call of npsol
C
C----------------------------------------------------------------------
      
C NOUT= unit to which output is directed

      NOUT= 6

C unit of file ./npsol_options.txt
      
      IOPTNS= 5
      OPEN(UNIT= IOPTNS, FILE= "./npsol_options.txt", STATUS= "old")

      CALL NPFILE(IOPTNS, INFORM)
      IF (INFORM.NE.0) THEN
         WRITE(NOUT, *) "ERROR IN CALL TO NPFILE:"
         WRITE(NOUT, *) "INFORM= ", INFORM
         STOP
      ENDIF

      CLOSE(UNIT= IOPTNS)

      END 

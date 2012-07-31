C-----------------------------------------------------------------------
C dummy routine called by cwrap_npsol
C forces output for fort.9 which is kept in memory buffer to be
C written into file fort.9
C-----------------------------------------------------------------------      
      SUBROUTINE CLOSE
      IMPLICIT NONE

      CLOSE(9)

      END
      

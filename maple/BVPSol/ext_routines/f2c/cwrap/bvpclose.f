C-----------------------------------------------------------------------
C dummy routine called by cwrap_npsol
C forces output for fort.9 which is kept in memory buffer to be
C written into file fort.9
C-----------------------------------------------------------------------      
      SUBROUTINE BVPCLOSE
      IMPLICIT NONE

      CLOSE(6)

      END
      

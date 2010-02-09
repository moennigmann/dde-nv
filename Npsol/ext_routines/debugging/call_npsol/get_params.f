C------------------------------------------------------------
C
C NPSOL/ext_routines/cwrap/get_params.f
C
C sets parameter values in common block
C
C notes: (i) not used by cwrap_npsol, but used in 
C   check* targets in ext_routines/debugging
C
C revision history:
C 021217 added stop command
C 021217 written by mmo
C
C------------------------------------------------------------
      SUBROUTINE GETPARAMS(NUMPARAMS, VALPARAMS)
      IMPLICIT NONE
      INTEGER NUMPARAMS
      INTEGER I1
      DOUBLEPRECISION VALPARAMS(NUMPARAMS)
      DOUBLEPRECISION PARAMS(200)
      
      COMMON/GLOBAL/PARAMS

C
C make sure common block is large enough
C
      IF (.NOT. NUMPARAMS.LE.200) THEN
         PRINT*, "Error in getparams, common block must not"
         PRINT*, "be larger than 200."
         STOP
      ENDIF
C     
C PARAMS cannot be input argument and in common block at the 
C same time
C
      DO 10 I1=1, NUMPARAMS, 1
         VALPARAMS(I1)= PARAMS(I1)
 10   CONTINUE

      END 

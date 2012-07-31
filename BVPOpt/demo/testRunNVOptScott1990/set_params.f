C------------------------------------------------------------
C
C NPSOL/ext_routines/cwrap/set_params.f
C
C sets parameter values in common block
C
C notes: (i) Some compilers warn about mismatch in common
C   block size. Here size of common block is fixed,
C   in problem-dependent files it will generally be fixed
C   to value passed to this subroutine in numparams.
C      
C revision history:
C 030324 increased max common block size to 1000
C 021217 added loop which assigns values passed in
C   VALPARAMS to common block variable PARAMS; (ii)
C   added stop command
C 021216 written by mmo
C
C------------------------------------------------------------
      SUBROUTINE SETPARAMS(NUMPARAMS, VALPARAMS)
      IMPLICIT NONE
      INTEGER NUMPARAMS
      INTEGER I1
      DOUBLEPRECISION VALPARAMS(NUMPARAMS)
      DOUBLEPRECISION PARAMS(1000)
      
      COMMON/GLOBAL/PARAMS

C
C make sure common block is large enough
C
      IF (.NOT. NUMPARAMS.LE.1000) THEN
         PRINT*, "Error in setparams, common block must not"
         PRINT*, "be larger than 1000."
         STOP
      ENDIF
C     
C PARAMS cannot be input argument and in common block at the 
C same time
C
      DO 10 I1=1, NUMPARAMS, 1
         PARAMS(I1)= VALPARAMS(I1)
 10   CONTINUE

      END 

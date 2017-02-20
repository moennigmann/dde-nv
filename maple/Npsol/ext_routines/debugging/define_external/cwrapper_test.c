#include <stdio.h>
#include <stdlib.h>
/***********************************************************************
  problem independent c subroutine cwrapper_npsol is called by maple
  and in turn calls fortran subroutine npsol_

  note that arguments NumParams and Params are not part of the 
  variables passed to npsol, but they are used to pass problem parameters
  in a common block
  **********************************************************************/
void cwrapper_test(int N, int M, double *B)
{
/*NOUT=&N; 
int counter = 0, counter2 = 0;  
 for( counter = 0; counter < N; counter++)
  {
    AOUT[counter]=A[counter];
  }

 for( counter2 = 0; counter2 < N*N; counter2++)
  {
    BOUT[counter2]=B[counter2];
  }*/
 
  return;
}

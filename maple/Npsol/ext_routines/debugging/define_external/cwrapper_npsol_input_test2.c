#include <stdio.h>
#include <stdlib.h>
/***********************************************************************
  problem independent c subroutine cwrapper_npsol is called by maple
  and in turn calls fortran subroutine npsol_

  note that arguments NumParams and Params are not part of the 
  variables passed to npsol, but they are used to pass problem parameters
  in a common block
  **********************************************************************/
void cwrapper_npsol(int N0, int N1, int N2, int N3, int N4,
		    int N5,int N6, int N7, int N8, int N9, int N10,
		    int N11, int N12, int N13, int N14, int N15, int N16,
                    int N17, int N18, int N19, int N20
		    )
{
 
}



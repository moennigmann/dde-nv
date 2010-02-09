/***********************************************************************
  problem independent c subroutine cwrapper_npsol is called by maple
  and in turn calls fortran subroutine npsol_

  note that arguments NumParams and Params are not part of the 
  variables passed to npsol, but they are used to pass problem parameters
  in a common block
  **********************************************************************/
void cwrapper_npsol(int N, int NCLIN, int NCNLN, int NROWA, int NROWJ,
		    int NROWR,
		    double *A, double *BL, double *BU,
		    int *INFORM, int *ITER , int *ISTATE,
		    double *C, double *CJAC, double *CLAMBDA,
		    double *OBJF, 
		    double *GRAD, double *R, double *XVEC,
		    int LENIW, int LENW,
		    int NumParams, double *Params 
		    )
{
  int *IW;
  double *W;
  
  /* allocate workspace for npsol */
  
  IW= (int *) malloc(LENIW* sizeof(int));
  W= (double *) malloc(LENW* sizeof(double));

  /*
    print variables passed from maple to c to screen
  */
}



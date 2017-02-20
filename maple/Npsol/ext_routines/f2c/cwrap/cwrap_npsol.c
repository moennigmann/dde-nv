/*------------------------------------------------------------

  NPSOL/ext_routines/cwrap/cwrap_npsol.c

  revision history:
  021216 added NumParams and Params which are used to pass parameter
    values to npsol in a common block; added call to set_params 
    which assigns values to parameters in common block; 
  01xxxx written by mmo

  ------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>

/***********************************************************************
  problem independent fortran subroutine npsol_
  as provided by Gill et al.
  **********************************************************************/
extern void npsol_(int *, int *, int *, int *, int *, int *,
		   double *, double *, double *,
		   void(*)(), void(*)(),
		   int *, int * , int *,
		   double *, double *, double *, double *, double *,
		   double *, double *,
		   int*, int *, double *, int *);

/***********************************************************************
  problem dependent fortran subroutine npsolconfun_
  which returns residuals of constraints,
  created automatically by maple
  **********************************************************************/
extern void npsolconfun_(int *, int *, int *, int *,
			  double *, double *, double *, double *,
			  int*
			  );

/***********************************************************************
  problem dependent fortran subroutine npsolobjfun_
  which returns value of profit function,
  created automatically by maple
  **********************************************************************/
extern void npsolobjfun_(int *, int *,
			  double *, double *, double *,
			  int*
			  );

/***********************************************************************
  problem independent fortran subroutine openopt_
  which reads npsol_options.txt and passes them to npsol_
  **********************************************************************/
extern void openopt_();


/***********************************************************************
  problem independent c subroutine cwrapper_npsol is called by maple
  and in turn calls fortran subroutine npsol_

  note that arguments NumParams and Params are not part of the 
  variables passed to npsol, but they are used to pass problem parameters
  in a common block
  **********************************************************************/
void cwrapper_npsol2(int N, int NCLIN, int NCNLN, int NROWA, int NROWJ,
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

  /* set parameter values in common block */
  setparams_(&NumParams, Params);

  /* open options file and transfer options to npsol */
  
  openopt_();

  
  npsol_(&N, &NCLIN, &NCNLN, &NROWA, &NROWJ,
	 &NROWR,
	 A, BL, BU,
	 npsolconfun_, npsolobjfun_, 
	 INFORM, ITER, ISTATE,
	 C, CJAC, CLAMBDA,
	 OBJF, GRAD, R, XVEC,
	 IW, &LENIW, W, &LENW
	 );

  /* close print file (unit 9) to make sure no output is lost */
  
  close_();
}

void cwrapper_npsol(int N, int NCLIN, int NCNLN, int NROWA, int NROWJ,
		    int NROWR,
		    double *A, double *B,
		    int *INFORM, int *ITER , int *ISTATE,
		    double *C, double *CJAC,
		    double *OBJF, 
		    double *GRAD, double *R, double *XVEC,
		    int *LE,
		    int NumParams, double *Params 
		    )
{
  int counter1, counter2, counter3;
  int leniw, lenw;
  int dimbl; 

  /* extract from LEN two arguments LENIW and LENW*/
  
  leniw = LE[0];
  lenw = LE[1]; 

  /* extract from B arguments BL, BU and CLAMBDA*/

 double *bl, *bu, *clambda;

  dimbl=N+NCLIN+NCNLN; 

  bl= (double *) malloc(dimbl* sizeof(double));
  bu= (double *) malloc(dimbl* sizeof(double));
  clambda= (double *) malloc(dimbl* sizeof(double));

 
 for( counter1 = 0; counter1 < dimbl; counter1++)
  {
    bl[counter1]=B[counter1];
  }
  
  for( counter2 = 0; counter2 < dimbl; counter2++)
  {
    bu[counter2]=B[counter2+dimbl];
  }
  
  for( counter3 = 0; counter3 < dimbl; counter3++)
  {
    clambda[counter3]=B[counter3+2*dimbl];
  }

  /* call cwrapper_npsol2 */
 cwrapper_npsol2(N, NCLIN, NCNLN, NROWA, NROWJ,
		    NROWR,
		    A, bl, bu,
		    INFORM, ITER, ISTATE,
		    C, CJAC, clambda,
		    OBJF, 
		    GRAD, R, XVEC,
		    leniw, lenw,
		    NumParams, Params
		    );
}






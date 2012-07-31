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
#include <math.h>

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
  //W= (double *) malloc(LENW* sizeof(double));
  W= (double *) calloc(LENW, sizeof(double));

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



/***********************************************************************
  Main procedure
  **********************************************************************/

int main(){
  
	int i1, i2, counter;
	  int mode, ncnln, nclin, n, nrowj, nrowa, nrowr, needc, nstate, numparams, leniw, lenw;
	  double *x, *c, *cjac, *par, *a, *bl, *bu, *clambda, *objf, *grad, *r, *pz;
	  int *inform, *iter, *istate;

	  mode= 0;
	  nstate= 0;
	  needc= 0;

	  /*

	     file system_dimension.h is created by NPSOL:-createDebugFiles()
	     contains values for ncnln, n, nrowj, numparams

	  */
	  ncnln= 17; //number of constraints
	  nclin= 0;
	  n    = 18; //number of variables
	  nrowj= 17; //=max(1, NCNLN)
	  nrowa =1; //=max(1, NCLIN)
	  nrowr =18; //=N
	  numparams= 0;


	  /*

	     allocate memory for variables,
	     constraint values and jacobian

	  */
	  x= (double *) malloc(n* sizeof(double));
	  par= (double *) malloc(numparams* sizeof(double));
	  c= (double *) malloc(ncnln* sizeof(double));
	  cjac= (double *) malloc(ncnln* n* sizeof(double));

	  /*

	     file variables.h is created by NPSOL:-createDebugFiles()
	     contains values for variables x[0]..x[n-1]
	     and parameters par[0]..par[numparams-1]

	  */

      	x[0]= (double) 4.283839e-01; /* mu0Star */
	x[1]= (double) 2.453232e-01; /* gStar */
	x[2]= (double) 4.670660e-01; /* mu0 */
	x[3]= (double) 2.136420e-01; /* g */
	x[4]= (double) 8.338739e-01; /* w[1] */
	x[5]= (double) -1.186127e-01; /* w[2] */
	x[6]= (double) -5.390597e-01; /* w[3] */
	x[7]= (double) -6.964741e-02; /* v[1] */
	x[8]= (double) -1.064988e+01; /* v[2] */
	x[9]= (double) 3.805392e-01; /* v[3] */
	x[10]= (double) 9.247846e-05; /* g1 */
	x[11]= (double) 7.029249e-01; /* u[1] */
	x[12]= (double) 3.367198e+01; /* u[2] */
	x[13]= (double) -1.825048e+00; /* u[3] */
	x[14]= (double) -8.405803e-05; /* kappa */
	x[15]= (double) -4.625613e+02; /* r[1] */
	x[16]= (double) 3.788442e+02; /* r[2] */
	x[17]= (double) 5.000000e-02; /* dist */


	  /*

	    print warning about option scale to screen

	  */

	  printf("#\n");
	  printf("# WARNING: If option scale was used when creating\n");
	  printf("#   instance of interface to NPSOL, residues and derivatives\n");
	  printf("#   are residues and derivaties of scaled system.\n");
	  printf("#\n");
	  printf("\n");

	  a=(double *) malloc(n*nrowa* sizeof(double));
	  for( counter = 0; counter < n*nrowa; counter++)
	  {
	    a[counter]= (double) 0.0;
	  }

	  bl=(double *) malloc((n+ncnln+nclin)* sizeof(double));
	  for( counter = 0; counter < n+ncnln+nclin; counter++)
	  {
		bl[counter]= (double) 0.0;
	  }
	  bl[0]=(double) 0.0;
	  bl[1]=(double) 0.0;
	  bl[2]=(double) 0.0;
	  bl[3]=(double) 0.0;
	  bl[4]=(double) -1e+30;
	  bl[5]=(double) -1e+30;
	  bl[6]=(double) -1e+30;
	  bl[7]=(double) -1e+30;
	  bl[8]=(double) -1e+30;
	  bl[9]=(double) -1e+30;
	  bl[10]=(double) -1e+30;
	  bl[11]=(double) -1e+30;
	  bl[12]=(double) -1e+30;
	  bl[13]=(double) -1e+30;
	  bl[14]=(double) -1e+30;
	  bl[15]=(double) -1e+30;
	  bl[16]=(double) -1e+30;
	  bl[17]=(double) 0.02828427124;


	  bu=(double *) malloc((n+ncnln+nclin)* sizeof(double));
	  for( counter = 0; counter < n+ncnln+nclin; counter++)
	  {
		 bu[counter]= (double) 0.0;
	  }
	  bu[0]=(double) 0.4717157288;
	  bu[1]=(double) 0.5;
	  bu[2]=(double) 0.5;
	  bu[3]=(double) 0.5;
	  bu[4]=(double) 1e+30;
	  bu[5]=(double) 1e+30;
	  bu[6]=(double) 1e+30;
	  bu[7]=(double) 1e+30;
	  bu[8]=(double) 1e+30;
	  bu[9]=(double) 1e+30;
	  bu[10]=(double) 1e+30;
	  bu[11]=(double) 1e+30;
	  bu[12]=(double) 1e+30;
	  bu[13]=(double) 1e+30;
	  bu[14]=(double) 1e+30;
	  bu[15]=(double) 1e+30;
	  bu[16]=(double) 1e+30;
	  bu[17]=(double) 1.0;

	  bu[34]=(double) 1e+30;


	  inform=(int *) malloc(sizeof(int));
	  iter=(int *) malloc(sizeof(int));

	  istate=(int *) malloc((n+ncnln+nclin)* sizeof(int));

	  for( counter = 0; counter < n+ncnln+nclin; counter++)
	  {
		  istate[counter]= (int) 0;
	  }

	  for( counter = 0; counter < ncnln; counter++)
	  {
	  	  c[counter]= (double) 0.0;
	  }

	  pz=(double *) malloc(122* sizeof(double));



	  for( counter = 0; counter < ncnln* n; counter++)
	  {
		  cjac[counter]= (double) 0.0;
	  }


	  clambda=(double *) malloc((n+ncnln+nclin)* sizeof(double));
	  for( counter = 0; counter < n+ncnln+nclin; counter++)
	  {
		  clambda[counter]= (double) 0.0;
	  }

	  objf= (double *) malloc(sizeof(double));

	  grad= (double *) malloc(n* sizeof(double));
	  for( counter = 0; counter < n; counter++)
	  {
		  grad[counter]= (double) 0.0;
	  }


	  r =  (double *) malloc(n*nrowr* sizeof(double));
	  for( counter = 0; counter < n*nrowr; counter++)
	  {
		  r[counter]= (double) 0.0;
	  }

  
  leniw=3* n+ nclin+ 2* ncnln;
  lenw=(2* n*n+ n* nclin+ 2* n* ncnln+ 20* n+ 11* nclin+ 21*ncnln)*5;
  
  cwrapper_npsol(n, nclin, ncnln, nrowa, nrowj,
		    nrowr,
		    a, bl, bu,
		    inform, iter, istate,
		    c, cjac, clambda,
		    objf, 
		    grad, r, x,
		    leniw, lenw,
		    numparams, par
		    );

  return 0; 
}

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


int main(){
  int i1, i2; 
  int mode, ncnln, nclin, n, nrowj, nrowa, nrowr, needc, nstate, numparams, leniw, lenw;
  double *x, *c, *cjac, *par, *a, *bl, *bu, *clambda, *objf, *grad, *r;
  int *inform, *iter, *istate;

  mode= 0;    
  nstate= 0;  
  needc= 0;   
  
  /* 

     file system_dimension.h is created by NPSOL:-CreateDebugFile()
     contains values for ncnln, n, nrowj, numparams

  */
  ncnln= 1;
  nclin= 1;
  n    = 2;
  nrowj= 1;
  nrowa =1;
  nrowr =2;
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

     file variables.h is created by NPSOL:-CreateDebugFile()
     contains values for variables x[0]..x[n-1]
     and parameters par[0]..par[numparams-1]

  */
  x[0]= (double) 0.000000e+00; /* x */
  x[1]= (double) 1.000000e+00; /* y */

  /*

    print warning about option scale to screen

  */
  printf("\n"); 
  printf("#\n");
  printf("# WARNING: If option scale was used when creating\n");
  printf("#   instance of interface to NPSOL, residues and derivatives\n");
  printf("#   are residues and derivaties of scaled system.\n");
  printf("#\n");
  printf("\n"); 

  a=(double *) malloc(n*nrowa* sizeof(double));
  a[0]=(double) -1.000000e+00;
  a[1]=(double) -1.000000e+00;

  bl=(double *) malloc((n+ncnln+nclin)* sizeof(double));
  bl[0]=(double) -100000000000000000000.000000;
  bl[1]=(double) -100000000000000000000.000000;
  bl[2]=(double) -2.000000e+00;
  bl[3]=(double)  0.000000e+00;

  bu=(double *) malloc((n+ncnln+nclin)* sizeof(double));
  bu[0]=(double) 100000000000000000000.000000;
  bu[1]=(double) 100000000000000000000.000000;
  bu[2]=(double) 100000000000000000000.000000;
  bu[3]=(double) 100000000000000000000.000000;

  inform=(int *) malloc(sizeof(int));
  iter=(int *) malloc(sizeof(int));

  istate=(int *) malloc((n+ncnln+nclin)* sizeof(int));
  istate[0]=(int) 0;
  istate[1]=(int) 0;
  istate[2]=(int) 0;
  istate[3]=(int) 0;

  c[0]=(double) 0.000000e+00;

  cjac[0]=(double) 0.000000e+00;
  cjac[1]=(double) 0.000000e+00;

  clambda=(double *) malloc((n+ncnln+nclin)* sizeof(double));
  clambda[0]=(double) 0.000000e+00;
  clambda[1]=(double) 0.000000e+00;
  clambda[2]=(double) 0.000000e+00;
  clambda[3]=(double) 0.000000e+00;

  objf= (double *) malloc(sizeof(double));

  grad= (double *) malloc(n* sizeof(double));
  grad[0]=(double) 0.000000e+00;
  grad[1]=(double) 0.000000e+00;

  r =  (double *) malloc(n*nrowr* sizeof(double));
  r[0]=(double) 0.000000e+00;
  r[1]=(double) 0.000000e+00;
  r[2]=(double) 0.000000e+00;
  r[3]=(double) 0.000000e+00;  

  leniw=9;
  lenw=86;
  
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

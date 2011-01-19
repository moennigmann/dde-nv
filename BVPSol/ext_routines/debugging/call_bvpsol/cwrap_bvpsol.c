/*------------------------------------------------------------

  BVPsol/ext_routines/cwrap/cwrap_bvpsol.c

  revision history:
  2011-01-14 written by dka

  ------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>


extern void callperiod_(int *,int *,int *,double *,double *,
                        int *,double *,double *,double *,
                        double *,double *,double *);

/***********************************************************************
  problem independent c subroutine cwrapper_bvpsol is called by maple
  and in turn calls fortran subroutine callperiod__
  **********************************************************************/
void cwrapper_bvpsol(int NN,int M,int IFAIL,double *X0,double *P,
                     int NPAR, double *TRPAR, double *XT,double *Y,
                     double *FM, double *HES,double *ERRY)
{
  
  callperiod_(&NN, &M, &IFAIL, X0, P,
	      &NPAR, TRPAR, XT, Y,
	      FM, HES, ERRY);
  
}


int main(){
  int n, m, ifail, npar, counter,counter1,counter2,counter3; 
  double per = 0,erry = 0;
  double *x0, *pars, *xt, *y, *fmult, *jac;
  
  n=4;
  m=30;
  npar=2;


  x0= (double *) malloc(n* sizeof(double));
  pars= (double *) malloc(npar* sizeof(double));
  xt= (double *) malloc(m* sizeof(double));
  y= (double *) malloc(n*m* sizeof(double));
  fmult= (double *) malloc(n*2* sizeof(double));
  jac= (double *) malloc(n*n* sizeof(double));



  ifail=0;
  
  per=14;
  
  pars[0]= (double) 0.19;
  pars[1]= (double) 0.0066651056;
  
  x0[0]= (double) 6.0;
  x0[1]= (double) 195.0;
  x0[2]= (double) 0.03;
  x0[3]= (double) 0.07;
  
  for( counter = 0; counter < m; counter++)
  {
    xt[counter]= (double ) 0;
  }

  for( counter1 = 0; counter1 < m*n; counter1++)
  {
    y[counter1]= (double ) 0;
  }


  for( counter2 = 0; counter2 < n*2; counter2++)
  {
    fmult[counter2]= (double ) 0;
  }

  for( counter3 = 0; counter3 < n*n; counter3++)
  {
    jac[counter3]= (double ) 0;
  }
  
  erry = 0;
  
  printf("\n"); 
  printf("#\n");
  printf("# Starting of PERIOD\n");
  printf("#\n");
  printf("\n"); 
  
  cwrapper_bvpsol(n, m, ifail, x0, &per,
  	          npar, pars, xt, y,
  	          fmult, jac, &erry);

  return 0; 
}

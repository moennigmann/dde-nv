/*---------------------------------------------------------------------- 

  NPSOL/ext_routines/debugging/check1.c 

  calls gnpsolconfun.c and prints constraints and their gradients
  to screen for debugging

  revision history:
  021217 added setting of parameters in common block used to
    pass parameter values to npsol
  021216 written by mmo

  ----------------------------------------------------------------------*/
#include <stdio.h> 
#include <stdlib.h>

int main(){
  int i1, i2; 
  int mode, ncnln, n, nrowj, needc, nstate, numparams;
  double *x, *c, *jac, *par; 

  mode= 0;    
  nstate= 0;  
  needc= 0;   
  
  /* 

     file system_dimension.h is created by NPSOL:-CreateDebugFile()
     contains values for ncnln, n, nrowj, numparams

  */
#include "system_dimension.h"

  /* 

     allocate memory for variables, 
     constraint values and jacobian

  */
  x= (double *) malloc(n* sizeof(double));
  par= (double *) malloc(numparams* sizeof(double)); 
  c= (double *) malloc(ncnln* sizeof(double));
  jac= (double *) malloc(ncnln* n* sizeof(double));

  /* 

     file variables.h is created by NPSOL:-CreateDebugFile()
     contains values for variables x[0]..x[n-1]
     and parameters par[0]..par[numparams-1]

  */
#include "variables.h"

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

  /* 

     print variables to screen

  */
  printf("Variables:= [\n");
  for (i1= 1; i1<= n-1; i1++) {
    printf("  x[%i]= %lf, \n", i1, x[i1-1]);
  }
  printf("  x[%i]= %lf \n];\n", n, x[n-1]);

  /* 

     print parameters to screen

  */
  printf("Parameters:= [\n");
  for (i1= 1; i1<= numparams-1; i1++) {
    printf("  par[%i]= %lf, \n", i1, par[i1-1]);
  }
  printf("  par[%i]= %lf \n];\n", numparams, par[numparams-1]);

  /*

    set parameter values in common block

  */
  printf("setting parameters ...\n");
  setparams_(&numparams, par);
  printf("...done\n");

  return 0; 
}

/* 

  for unknown reason libf2c.so expects a MAIN__

*/ 
int MAIN__(){return 0;}

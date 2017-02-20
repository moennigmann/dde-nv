/* 
  NPSOL/ext_routines/debugging/check2.c

  calls npscolconfun.c and prints constraints and their gradients
  to screen for debugging

  021216 written by mmo
*/
#include <stdio.h> 
#include <stdlib.h>
// extern int npsolconfun_(double *variables, double *residues)

int main(){
  int i1, i2; 
  int mode, ncnln, n, nrowj, needc, nstate, numparams; 
  double *x, *c, *par; 

  mode= 0;    // currently unused by maple/adifor-generated code
  nstate= 0;  // currently unused by maple/adifor-generated code
  needc= 0;   // currently unused by maple/adifor-generated code
  
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

  /* 
     file variables.h is created by NPSOL:-CreateDebugFile()
     contains values for variables x[0]..x[n-1] and
     for parameters par[0]..par[numparams-1]
  */
#include "variables.h"

  /*
    print brief explanation
  */
  printf("\n"); 
  printf("#\n");
  printf("# Using npsolconfun.f as generated for use with ADIFOR\n");
  printf("#\n");

  /*
    print warning about option scale to screen
  */
  printf("\n"); 
  printf("#\n");
  printf("# WARNING: If option scale was used when creating\n");
  printf("#   instance of interface to NPSOL, residues are\n");
  printf("#   residues of scaled system.\n");
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
  setparams_(&numparams, par);

  /*
    evaluate nonlinear constraints and jacobian
  */
  npsolconfun_(x, c);

  /* 
     print residues to screen
  */
  printf("\n"); 
  printf("Residues:= [\n");
  for (i1= 1; i1<= ncnln-1; i1++) {
    printf("  c[%i]= %lf, \n", i1, c[i1-1]);
  }
  printf("  c[%i]= %lf \n];\n", ncnln, c[ncnln-1]);

  return 0; 
}

/* 
  for unknown reason libf2c.so expects a MAIN__
*/ 
int MAIN__(){return 0;}

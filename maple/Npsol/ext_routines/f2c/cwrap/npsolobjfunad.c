/* 
  NPSOL/ext_routines/cwrap/npsolobjfunad.c

  wrapper for call to adifor-generated subroutine

  revision history:
  030418 replaced leading dimension of objgrad, nobjgrad, by number
    of variables w.r.t. which derivatives are needed, n
  021122 written by mmo
*/
void npsolobjfun_(int *mode, int *n, double x[*n], double *objf,
		  double objgrad[*n], int *nstate)
{
    double g_x[*n][*n];
    int i1, i2; 
    {
      /* debugging 
      printf("\n"); 
      printf("have entered npsolobjfunad\n");
	printf("*n= %i\n", *n);
	for (i1= 1; i1<= *n; i1++) {
	    printf("x(%i)= %lf \n", i1, x[i1-1]);
	}
      */

      /* setup seed matrix */
        for (i1= 1; i1<= *n; i1++){
	    for (i2= 1; i2<= *n; i2++){
		g_x[i1-1][i2-1]= 0.0;
	    }
	    g_x[i1-1][i1-1]= 1.0;
	}


      /* call adifor-generated code */
	gnpsolobjfun_(n, x, g_x, n, objf, objgrad, n);

      /* call adifor error handler */
        ehrpt_();

	return;
    }

}


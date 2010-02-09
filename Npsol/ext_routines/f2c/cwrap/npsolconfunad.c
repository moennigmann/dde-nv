/* 
  NPSOL/ext_routines/cwrap/npsolconfunad.c

  wrapper for call to adifor-generated subroutine

  notes:
  (i) fortran expects jac(ncnln, n), c stores transposed matrix,
  since in the c code we have to loop over elements, we must
  declare jac(n, ncnln)

  revision history:
  021121 removed transposing of jacobian, not necessary since c stores
    columns, fortran stores rows
  021120 written by mmo
*/
void npsolconfun_(int *mode, int *ncnln, int *n, int *nrowj, int *needc, double x[*n], double c[*ncnln], double jac[*ncnln* *n], int *nstate)
{
    double g_x[*n][*n];
    double g_c[*n* *ncnln];
    int i1, i2; 
    {
      /* debugging 
	printf("*n= %i\n", *n);
	printf("*ncnln= %i\n", *ncnln);
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
	gnpsolconfun_(n, x, g_x, n, c, g_c, n);

	/* call adifor error handler */
        ehrpt_();

	/* debugging 
	for (i1= 1; i1<= *ncnln; i1++) {
	    printf("c(%i)= %lf \n", i1, c[i1-1]);
	}
	*/

	/* debugging 
	for (i1= 1; i1<= *n* *ncnln; i1++) {
	    printf("g_c(%i)= %lf \n ", i1, g_c[i1-1]); 
	}
	*/

	for (i1= 1; i1<= *ncnln; i1++) {
	  for (i2= 1; i2<= *n; i2++) {
	    jac[(i2- 1)* *ncnln + i1- 1]= g_c[(i1- 1)* *n + i2- 1];
	  }
	}

	return; 

    }

}


#include <stdio.h>
#include <stdlib.h>
#include "olsen.h"

void calltides_(int *nvar1, int *npar1, double *p1, double *v1, double *xend, double *outvals) {

	int nvar;
	int npar;
    int nfun = 0;
	int nipt = 2;
	nvar=*nvar1;
	npar=*npar1;


	double v[nvar], p[npar], lt[nipt];
	double tolrel, tolabs;
	double** olsen_DataMatrix;
	FILE   *fd;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	int counter;
	for( counter = 0; counter < npar; counter++)
	{
	  p[counter]= p1[counter];
	}

/* --- INITIAL VALUES --- */
	for( counter = 0; counter < nvar; counter++)
	{
	  v[counter]= v1[counter];
	}

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0 ;
	lt[1]=*xend;

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-16 ;
	tolabs = 1.e-16 ;

/***********************************************************/
/***********************************************************/
/*        OUTPUT: file   &data matrix                     */
/***********************************************************/
/***********************************************************/

	fd = fopen("olsenResult.txt", "w");
	Array2DB_init(&olsen_DataMatrix, nipt, olsen_columns());


/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

	dp_tides(olsen, nvar, npar, nfun, v, p,
			lt, nipt, tolrel, tolabs, olsen_DataMatrix, fd);

///// TRANSLATE OUTPUT INTO FORTRAN

	int sizematr;
	sizematr=1+nvar+(nvar+npar)*nvar+((nvar+npar)*(nvar+npar+1)*nvar)/2;

	for( counter = 0; counter < sizematr; counter++)
	{
		outvals[counter]= olsen_DataMatrix[1][counter];
	}

	fclose(fd);

	return;
}



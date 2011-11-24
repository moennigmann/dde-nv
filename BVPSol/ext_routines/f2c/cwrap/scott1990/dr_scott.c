/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (1.20) June 8, 2011, 15:42

	Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
	Grupo de Mecanica Espacial
	University of Zaragoza
	SPAIN

	http://gme.unizar.es/software/tides
	Contact: <tides@unizar.es>

	This file is part of TIDES.

	TIDES is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	TIDES is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with TIDES.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "scott.h"

int main() {

	int nvar = 3;
	int npar = 2;
	int nfun = 0;
	int nipt = 2;
	double v[nvar], p[npar], lt[nipt];
	double tolrel, tolabs; 
	double** scott_DataMatrix;
	FILE   *fd;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 0.3 ; 
	p[1] = 0.4 ; 

/* --- INITIAL VALUES --- */
	v[0] = 1. ; 
	v[1] = 2. ; 
	v[2] = 3. ; 

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0 ; 
	lt[1] = 25.5 ; 

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-16 ;
	tolabs = 1.e-16 ;

/***********************************************************/
/***********************************************************/
/*        OUTPUT: file   &data matrix                     */
/***********************************************************/
/***********************************************************/

	fd = fopen("scottResult.txt", "w");
	Array2DB_init(&scott_DataMatrix, nipt, scott_columns());


/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

	dp_tides(scott, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, scott_DataMatrix, fd);


	fclose(fd); 

	return 0;
}



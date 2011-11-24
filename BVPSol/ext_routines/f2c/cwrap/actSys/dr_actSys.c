/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (1.20) March 7, 2011, 15:28

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
#include "actSys.h"

int main() {

	int nvar = 2;
	int npar = 2;
	int nfun = 0;
	int nipt = 2;
	double v[nvar], p[npar], lt[nipt];
	double tolrel, tolabs; 
	double** actSys_DataMatrix;
	FILE   *fd;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 50.636141 ; 
	p[1] = 109.77726 ; 

/* --- INITIAL VALUES --- */
	v[0] = 19.6878809 ; 
	v[1] = 60.115845 ; 

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0 ; 
	lt[1] = 5.59887 ; 

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-16 ;
	tolabs = 1.e-16 ;

/***********************************************************/
/***********************************************************/
/*        OUTPUT: file   &data matrix                     */
/***********************************************************/
/***********************************************************/

	fd = fopen("actSysResult.txt", "w");
	Array2DB_init(&actSys_DataMatrix, nipt, actSys_columns());


/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

	dp_tides(actSys, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, actSys_DataMatrix, fd);


	fclose(fd); 

	return 0;
}



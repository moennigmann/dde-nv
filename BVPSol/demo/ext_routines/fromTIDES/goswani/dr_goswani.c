/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (1.20) June 13, 2013, 11:46

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
#include "goswani.h"

int main() {

	int nvar = 4;
	int npar = 2;
	int nfun = 0;
	int nipt = 2;
	double v[nvar], p[npar], lt[nipt];
	double tolrel, tolabs; 
	double** goswani_DataMatrix;
	FILE   *fd;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 1. ; 
	p[1] = 2. ; 

/* --- INITIAL VALUES --- */
	v[0] = -0.3233 ; 
	v[1] = 0.2186 ; 
	v[2] = -0.3771 ; 
	v[3] = -1.0918 ; 

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0 ; 
	lt[1] = 0.73472 ; 

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-16 ;
	tolabs = 1.e-16 ;

/***********************************************************/
/***********************************************************/
/*        OUTPUT: file   &data matrix                     */
/***********************************************************/
/***********************************************************/

	fd = fopen("goswaniResult.txt", "w");
	Array2DB_init(&goswani_DataMatrix, nipt, goswani_columns());


/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

	dp_tides(goswani, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, goswani_DataMatrix, fd);


	fclose(fd); 

	return 0;
}



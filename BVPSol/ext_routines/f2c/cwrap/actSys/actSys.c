/****************************************************************************
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

#include "actSys.h"


long  actSys(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1])
{
	static int   VARIABLES        = 2;
	static int   PARAMETERS       = 2;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 16;
	static int   PARTIALS_VARS    = 4;
	static long  NUM_DERIVATIVES  = 35;
	static long  NUM_COLUMNS      = 70;

	static int   POS__PARTIALS[4] = {1,2,3,4};
	static int   POS_FUNCTIONS[1] = {0};

	static long  POS_ACCUM[36] = {0,1,3,5,7,9,12,16,20,24,27,31,35,38,42,45,49,55,61,67,73,81,89,95,103,109,113,119,125,131,139,145,149,155,161,165};
	static long  POS_COEFS[165] = {1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,3,3,1,1,1,2,2,1,1,1,1,2,2,1,1,1,1,2,2,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,3,3,1,1,1,2,2,1,1,1,1,2,2,1,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,3,3,1,1,1,2,2,1,1,1,2,1,1,2,1,1,3,3,1};
	static long  POS_PREVI[165] = {0,0,1,0,2,0,3,0,4,0,1,5,0,2,1,6,0,3,1,7,0,4,1,8,0,2,9,0,3,2,10,0,4,2,11,0,3,12,0,4,3,13,0,4,14,0,1,5,15,0,2,1,6,5,16,0,3,1,7,5,17,0,4,1,8,5,18,0,2,9,1,6,19,0,3,2,10,1,7,6,20,0,4,2,11,1,8,6,21,0,3,12,1,7,22,0,4,3,13,1,8,7,23,0,4,14,1,8,24,0,2,9,25,0,3,2,10,9,26,0,4,2,11,9,27,0,3,12,2,10,28,0,4,3,13,2,11,10,29,0,4,14,2,11,30,0,3,12,31,0,4,3,13,12,32,0,4,14,3,13,33,0,4,14,34};
	static long  POS_PREIV[165] = {0,1,0,2,0,3,0,4,0,5,1,0,6,1,2,0,7,1,3,0,8,1,4,0,9,2,0,10,2,3,0,11,2,4,0,12,3,0,13,3,4,0,14,4,0,15,5,1,0,16,5,6,1,2,0,17,5,7,1,3,0,18,5,8,1,4,0,19,6,1,9,2,0,20,6,7,1,10,2,3,0,21,6,8,1,11,2,4,0,22,7,1,12,3,0,23,7,8,1,13,3,4,0,24,8,1,14,4,0,25,9,2,0,26,9,10,2,3,0,27,9,11,2,4,0,28,10,2,12,3,0,29,10,11,2,13,3,4,0,30,11,2,14,4,0,31,12,3,0,32,12,13,3,4,0,33,13,3,14,4,0,34,14,4,0};

	static long  POS_ACCUM_S[36] = {0,1,2,3,4,5,7,9,11,13,15,17,19,21,23,25,28,31,34,37,40,44,48,51,55,58,61,64,67,70,74,77,80,83,86,89};
	static long  POS_COEFS_S[89] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1};
	static long  POS_PREVI_S[89] = {0,0,0,0,0,0,1,0,2,0,3,0,4,0,2,0,3,0,4,0,3,0,4,0,4,0,1,5,0,1,5,0,1,5,0,1,5,0,2,9,0,3,2,10,0,4,2,11,0,3,12,0,4,3,13,0,4,14,0,2,9,0,2,9,0,2,9,0,3,12,0,4,3,13,0,4,14,0,3,12,0,3,12,0,4,14,0,4,14};
	static long  POS_PREIV_S[89] = {0,1,2,3,4,5,1,6,1,7,1,8,1,9,2,10,2,11,2,12,3,13,3,14,4,15,5,1,16,6,2,17,7,3,18,8,4,19,6,1,20,6,7,1,21,6,8,1,22,7,1,23,7,8,1,24,8,1,25,9,2,26,10,3,27,11,4,28,10,2,29,10,11,2,30,11,2,31,12,3,32,13,4,33,13,3,34,14,4};


	if(ORDER < 0) return NUM_COLUMNS;

	static int  NOT_INITIALIZED = 1;
	if(NOT_INITIALIZED)
	{
		set_iterations();
		NOT_INITIALIZED = 0; 
	}
	set_max_order(ORDER);

	realNUM var[VARIABLES+1][NUM_DERIVATIVES][ORDER+1];
	realNUM par[PARAMETERS][NUM_DERIVATIVES][ORDER+1];
	realNUM link[LINKS][NUM_DERIVATIVES][ORDER+1];
	variables_init(var,v,t);
	parameters_init(par,p);
	links_init(link);
	derivatives_init(var,par,v);

	int i;
	for(i=0;  i<=ORDER; i++) {
		var_t(link[15],var[1], i);
		var_t(link[14],var[2], i);
		add_t_c("1.",var[1],link[0],i);
		mul_t_c("-0.1",par[0],link[1],i);
		mul_t_c("-0.1",var[2],link[2],i);
		mul_t(var[1],var[1],link[3],i);
		mul_t_c("0.1",link[3],link[4],i);
		mul_t(link[1],var[1],link[5],i);
		mul_t(link[2],par[0],link[6],i);
		add_t(link[0],link[4],link[7],i);
		divide_t(link[5],link[7],link[8],i);
		divide_t(link[6],link[7],link[9],i);
		add_t_c("-1.",link[9],link[10],i);
		add_t_c("-0.2",link[8],link[11],i);
		mul_t(link[10],var[1],link[12],i);
		mul_t(link[11],var[2],link[13],i);
		add_t_c("100.",link[13],link[14],i);
		add_t(link[12],par[1],link[15],i);
	}

	write_solution(cvfd,var,link);

	return NUM_COLUMNS;
}

long  actSys_columns()
{
	 return 71;
}

long  actSys_pos_der(char *der)
{
	static char* STR_DER[35] = {"0000","1000","0100","0010","0001","2000","1100","1010","1001","0200","0110","0101","0020","0011","0002","3000","2100","2010","2001","1200","1110","1101","1020","1011","1002","0300","0210","0201","0120","0111","0102","0030","0021","0012","0003"};
	long i;
	for(i=0; i < 35; i++)
		if(strcmp(der,STR_DER[i]) == 0) return i;
	return -1;
}

long  actSys_variable_column(int v, char *der)
{
	 return position_variable(v, actSys_pos_der, der);
}

long  actSys_function_column(int f, char *der)
{
	 return position_function(f, actSys_pos_der, der);
}



/****************************************************************************
	This file has been created by MathTIDES (1.20) February 22, 2011, 14:47

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

#include "olsen3D.h"


long  olsen3D(realNUM t, realNUM v[], realNUM p[], int ORDER, realNUM cvfd[][ORDER+1])
{
	static int   VARIABLES        = 4;
	static int   PARAMETERS       = 3;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 33;
	static int   PARTIALS_VARS    = 6;
	static long  NUM_DERIVATIVES  = 28;
	static long  NUM_COLUMNS      = 112;

	static int   POS__PARTIALS[6] = {1,2,3,4,5,6};
	static int   POS_FUNCTIONS[1] = {0};

	static long  POS_ACCUM[29] = {0,1,3,5,7,9,11,13,16,20,24,28,32,36,39,43,47,51,55,58,62,66,70,73,77,81,84,88,91};
	static long  POS_COEFS[91] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,1};
	static long  POS_PREVI[91] = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,1,7,0,2,1,8,0,3,1,9,0,4,1,10,0,5,1,11,0,6,1,12,0,2,13,0,3,2,14,0,4,2,15,0,5,2,16,0,6,2,17,0,3,18,0,4,3,19,0,5,3,20,0,6,3,21,0,4,22,0,5,4,23,0,6,4,24,0,5,25,0,6,5,26,0,6,27};
	static long  POS_PREIV[91] = {0,1,0,2,0,3,0,4,0,5,0,6,0,7,1,0,8,1,2,0,9,1,3,0,10,1,4,0,11,1,5,0,12,1,6,0,13,2,0,14,2,3,0,15,2,4,0,16,2,5,0,17,2,6,0,18,3,0,19,3,4,0,20,3,5,0,21,3,6,0,22,4,0,23,4,5,0,24,4,6,0,25,5,0,26,5,6,0,27,6,0};

	static long  POS_ACCUM_S[29] = {0,1,2,3,4,5,6,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49};
	static long  POS_COEFS_S[49] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	static long  POS_PREVI_S[49] = {0,0,0,0,0,0,0,0,1,0,2,0,3,0,4,0,5,0,6,0,2,0,3,0,4,0,5,0,6,0,3,0,4,0,5,0,6,0,4,0,5,0,6,0,5,0,6,0,6};
	static long  POS_PREIV_S[49] = {0,1,2,3,4,5,6,7,1,8,1,9,1,10,1,11,1,12,1,13,2,14,2,15,2,16,2,17,2,18,3,19,3,20,3,21,3,22,4,23,4,24,4,25,5,26,5,27,6};


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
		var_t(link[28],var[1], i);
		var_t(link[29],var[2], i);
		var_t(link[30],var[3], i);
		var_t(link[32],var[4], i);
		add_t_c("-8.",var[1],link[0],i);
		mul_t_c("-1.",par[2],link[1],i);
		mul_t_c("-1.",var[2],link[2],i);
		mul_t_c("-0.1",var[1],link[3],i);
		mul_t_c("-0.1",var[1],link[4],i);
		mul_t_c("0.3",var[1],link[5],i);
		mul_t(par[1],var[2],link[6],i);
		mul_t(var[3],var[3],link[7],i);
		add_t_c("-20.",link[6],link[8],i);
		mul_t_c("-500.",link[7],link[9],i);
		mul_t_c("500.",link[7],link[10],i);
		mul_t(link[0],link[1],link[11],i);
		mul_t(link[2],par[1],link[12],i);
		mul_t(link[3],var[2],link[13],i);
		mul_t(link[4],var[2],link[14],i);
		mul_t(link[5],var[2],link[15],i);
		add_t_c("0.00001",link[9],link[16],i);
		mul_t(link[8],var[3],link[17],i);
		mul_t(link[12],var[3],link[18],i);
		mul_t(link[13],par[0],link[19],i);
		mul_t(link[14],par[0],link[20],i);
		mul_t(link[15],par[0],link[21],i);
		add_t_c("-5.35",link[19],link[22],i);
		add_t_c("0.825",link[18],link[23],i);
		add_t(link[16],link[17],link[24],i);
		mul_t(link[19],var[4],link[25],i);
		mul_t(link[20],var[4],link[26],i);
		mul_t(link[21],var[4],link[27],i);
		add_t(link[11],link[26],link[28],i);
		add_t(link[23],link[25],link[29],i);
		add_t(link[24],link[27],link[30],i);
		mul_t(link[22],var[4],link[31],i);
		add_t(link[10],link[31],link[32],i);
	}

	write_solution(cvfd,var,link);

	return NUM_COLUMNS;
}

long  olsen3D_columns()
{
	 return 113;
}

long  olsen3D_pos_der(char *der)
{
	static char* STR_DER[28] = {"000000","100000","010000","001000","000100","000010","000001","200000","110000","101000","100100","100010","100001","020000","011000","010100","010010","010001","002000","001100","001010","001001","000200","000110","000101","000020","000011","000002"};
	long i;
	for(i=0; i < 28; i++)
		if(strcmp(der,STR_DER[i]) == 0) return i;
	return -1;
}

long  olsen3D_variable_column(int v, char *der)
{
	 return position_variable(v, olsen3D_pos_der, der);
}

long  olsen3D_function_column(int f, char *der)
{
	 return position_function(f, olsen3D_pos_der, der);
}



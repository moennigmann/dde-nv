/******************************************************************************
 revision history:
 040220 changed name of the adifor file for the res jacobian from
   gres to jres since gres denotes the adifor file for sensitivities -jge
 040205 changes to allow DAE system sensitivity parameters passing between
   maple and ddaspk: senpar to *senpar -jge
 020815 changed argument g_res_ in call to ddaspk_ to gres_ to be in
   acc. with sed command in CreateSharedObject which removes underscore
   in g_res; 
 020311 changes to allow for DAE system parameter passing between maple and 
   ddaspk: changed rpar to *rpar, call to fortran ddaspk changed accordingly  
 020309 changed T, rtol, atol to *T, *rtol, *atol in cwrapper_ddaspk and
   &T, &rtol, &atol to T, rtol, atol in ddaspk accordingly
 020301 going from maple6 to maple7: introduced pointers in cwrapper_ddapsk 
   according to changes in define_external in DDASPK:-CreateInstance,
   introduced dereferencings in call to ddaspk_ accordingly 
   (see repository revision 1.1.1.1 for version before these changes) 
 0111xx written by mkl, modeled after corresponding NPSOL file

 Copyright (C) 2002 by M. Monnigmann for 
 Lehrstuhl f. Prozesstechnik,
 RWTH Aachen. All rights reserved.

 *****************************************************************************/
#include <stdio.h>


/***********************************************************************
  problem independent fortran subroutine npsol_
  as provided by Gill et al.
  **********************************************************************/
extern void ddaspk_(void (*)(), int *, double *, double *, double *,
		    double *, double *, double *, double *, int *,
		    double *, int *, int *, int *,
		    double *, int *, void (*)(), void(*)(),
		    double *, void (*)()
		    );

extern void res_(double T, double y, double yprime, double *cj,
		 double delta, int ires, int ipar, double rpar,
		 double senpar
		 );

extern void jac_();

extern void psol_();

extern void gres_();

extern void jres_();
 		    
/***********************************************************************
  problem independent c subroutine cwrapper_ddaspk is called by maple
  and in turn calls fortran subroutine ddaspk_
  **********************************************************************/
void cwrapper_ddaspk(int neq, double *T, double *y, double *yprime,
		     double tout, double *info, double *rtol, double *atol,
		     int *idid, double *RWORK, int LRW, int *IWORK, int LIW,
		     double *rpar, int ipar, double *senpar
		     )
{
    /*printf("y= %f, %f \n", rpar[0], rpar[1]);*/
    
  ddaspk_(res_, &neq, T, y, yprime, &tout, info, rtol, atol, idid,
	  RWORK, &LRW, IWORK, &LIW, rpar, &ipar, jac_, psol_, senpar, gres_
	  );

}

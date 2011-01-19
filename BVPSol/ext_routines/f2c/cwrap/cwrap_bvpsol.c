/*------------------------------------------------------------

  BVPsol/ext_routines/cwrap/cwrap_bvpsol.c

  revision history:
  2011-01-14 written by dka

  ------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>


extern void callperiod_(int *,int *,int *,double *,double *,
                        int *,double *,double *,double *,
                        double *,double *,double *,int *,
                        double *,double *,double *,double *);
                        
extern void bvpopen_();                        

extern void bvpclose_();    

/***********************************************************************
  problem independent c subroutine cwrapper_bvpsol is called by maple
  and in turn calls fortran subroutine callperiod__
  **********************************************************************/
void cwrapper_bvpsol(int NN, int M, int IPRINT,double *X0, double P,
                     int NPAR, double *TRPAR, double *XT,double *Y,
                     double *FM, double *HES, double *POUT, int *IFAIL, 
                     double *ERRY,double *FAL,double *FXAL,double *FXX)
{
 
   bvpopen_();
   
   callperiod_(&NN, &M, &IPRINT, X0, &P,
	      &NPAR, TRPAR, XT, Y,
	      FM, HES, POUT, IFAIL, 
	      ERRY, FAL, FXAL, FXX);
	      
   bvpclose_();
	      
  
}

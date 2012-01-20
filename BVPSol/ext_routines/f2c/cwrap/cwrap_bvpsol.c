/*------------------------------------------------------------

  BVPsol/ext_routines/cwrap/cwrap_bvpsol.c

  revision history:
  2011-01-14 written by dka

  ------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>


extern void callperiod_(int *,int *,int *,double *,double *,
                        int *,double *,double *,
                        double *,double *,int *,
                        double *,double *,double *,double *, 
                        double *, int *,double *,double *, 
                        double *);
                        
extern void bvpopen_();                        

extern void bvpclose_();    

/***********************************************************************
  problem independent c subroutine cwrapper_bvpsol is called by maple
  and in turn calls fortran subroutine callperiod__
  **********************************************************************/
void cwrapper_bvpsol(int NN, int M, int IPRINT,double *X0, double P,
                     int NPAR, double *TRPAR, double *Y,
                     double *FM, double *POUT, int *IFAIL, 
                     double *ERRY,double *FAL,double *FXAL,double *FXX, 
                     double *FX, int NPHASE)
{
 
   double *FPP;
   double *FXXP;
   double *FXPP;
   
   FPP= (double *) malloc(NN*NPAR*NPAR*sizeof(double));
   FXXP= (double *) malloc(NN*NN*NN*NPAR*sizeof(double));
   FXPP= (double *) malloc(NN*NN*NPAR*NPAR*sizeof(double));
   
   bvpopen_();
   
   callperiod_(&NN, &M, &IPRINT, X0, &P,
	      &NPAR, TRPAR, Y,
	      FM, POUT, IFAIL, 
	      ERRY, FAL, FXAL, FXX, 
	      FX, &NPHASE, FPP, FXXP, FXPP);
	      
   bvpclose_();
	      
  
}

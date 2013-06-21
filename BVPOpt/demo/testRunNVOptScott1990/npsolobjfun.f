      subroutine npsolobjfun(mode,n,x,objf,objgrad,nstate)
      integer mode
      integer n
      doubleprecision x(18)
      doubleprecision objf
      doubleprecision objgrad(18)
      integer nstate

      doubleprecision ERRY
      doubleprecision FM(3,2)
      doubleprecision FP(3,2)
      doubleprecision FPP(3,2,2)
      doubleprecision FX(3,3)
      doubleprecision FXP(3,3,2)
      doubleprecision FXPP(3,3,2,2)
      doubleprecision FXX(3,3,3)
      doubleprecision FXXP(3,3,3,2)
      integer IFAIL
      integer NPHASE
      doubleprecision OBJPAR(2)
      doubleprecision P
      doubleprecision POUT
      doubleprecision X0(3)
      doubleprecision Y(3,1000)
      integer i1
      integer m1
      integer n1
      integer n2
      integer n3
      integer n_in_objfun
      doubleprecision z(1)

      common/global/par
      doubleprecision par(1)

        n_in_objfun = 18
        if ( .not. n_in_objfun .eq. n) then
          write(6,*) 'Error in confun: n differs from number of variable
     #s'
          call exit(0)
        endif
        X0(1) = 0.5D1
        X0(2) = 0.69D-1
        X0(3) = 0.15D1
        P = 0.25D2
        OBJPAR(1) = x(1)
        OBJPAR(2) = x(2)
        NPHASE = 2
        n1 = 3
        m1 = 1000
        n2 = 2
        n3 = -1
        call CALLPERIOD(n1,m1,n3,X0,P,n2,OBJPAR,Y,FM,POUT,IFAIL,ERRY,FP,
     #FXP,FXX,FX,NPHASE,FPP,FXXP,FXPP)
        z(1) = 0
        objf = 0.D0
        do 1000 i1 = 1,m1,1
          objf = objf+0.1D1*Y(2,i1)/DBLE(m1)
1000    continue
        objf = -objf
        objf = 10000*objf
        call printParObjFun(n2,OBJPAR,objf)
        if (IFAIL .ne. 0) then
          write(6,*) 'FAILED IN OBJECTIVE'
          call exit(0)
        endif
        return
      end

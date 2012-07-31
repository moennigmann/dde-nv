      subroutine npsolconfun(mode,ncnln,n,nrowj,needc,x,c,cjac,nstate)
      integer mode
      integer ncnln
      integer n
      integer nrowj
      integer needc
      doubleprecision x(18)
      doubleprecision c(17)
      doubleprecision cjac(17,18)
      integer nstate

      doubleprecision ERRY
      doubleprecision F0(3)
      doubleprecision F0P(3,2)
      doubleprecision FM(3,2)
      doubleprecision FP(3,2)
      doubleprecision FPP(3,2,2)
      doubleprecision FX(3,3)
      doubleprecision FXP(3,3,2)
      doubleprecision FXPP(3,3,2,2)
      doubleprecision FXX(3,3,3)
      doubleprecision FXXP(3,3,3,2)
      integer IFAIL
      doubleprecision J0(3,3)
      doubleprecision J0P(3,3,2)
      integer NPHASE
      doubleprecision P
      doubleprecision POUT
      doubleprecision T1
      doubleprecision TRPAR(2)
      doubleprecision X0(3)
      doubleprecision X0S(3)
      doubleprecision Y(3,50)
      integer i1
      integer m1
      integer n1
      integer n2
      integer n3
      integer n_in_confun
      integer ncnln_in_confun
      doubleprecision z(1)

      common/global/par
      doubleprecision par(1)

        n_in_confun = 18
        if ( .not. n_in_confun .eq. n) then
          write(6,*) 'Error in confun: n differs from number of variable
     #s'
          call exit(0)
        endif
        ncnln_in_confun = 17
        if ( .not. ncnln_in_confun .eq. ncnln) then
          write(6,*) 'Error in confun: ncnln differs from number of equa
     #tions'
          call exit(0)
        endif
        z(1) = 0
        P = 0.25D2
        X0(1) = 0.5D1
        X0(2) = 0.69D-1
        X0(3) = 0.15D1
        TRPAR(1) = x(3)
        TRPAR(2) = x(4)
        NPHASE = 2
        n1 = 3
        m1 = 50
        n2 = 2
        n3 = -1
        call printParConFun(n2,TRPAR)
        call CALLPERIOD(n1,m1,n3,X0,P,n2,TRPAR,Y,FM,POUT,IFAIL,ERRY,FP,F
     #XP,FXX,FX,NPHASE,FPP,FXXP,FXPP)
        do 1000 i1 = 1,3,1
          X0S(i1) = Y(i1,1)
1000    continue
        T1 = 0.D0
        call FCN(T1,X0S,F0,TRPAR)
        call DFDY(T1,X0S,J0,TRPAR)
        call funcF0P(X0S,F0P,TRPAR)
        call funcJ0P(X0S,J0P,TRPAR)
        c(1) = FX(1,1)*x(5)+FX(1,2)*x(6)+FX(1,3)*x(7)+x(5)
        c(2) = FX(2,1)*x(5)+FX(2,2)*x(6)+FX(2,3)*x(7)+x(6)
        c(3) = FX(3,1)*x(5)+FX(3,2)*x(6)+FX(3,3)*x(7)+x(7)
        c(4) = x(5)**2+x(6)**2+x(7)**2-1.D0
        c(5) = FX(1,1)*x(8)+FX(2,1)*x(9)+FX(3,1)*x(10)+x(8)+x(11)*x(5)
        c(6) = FX(1,2)*x(8)+FX(2,2)*x(9)+FX(3,2)*x(10)+x(9)+x(11)*x(6)
        c(7) = FX(1,3)*x(8)+FX(2,3)*x(9)+FX(3,3)*x(10)+x(10)+x(11)*x(7)
        c(8) = x(5)*x(8)+x(6)*x(9)+x(7)*x(10)-1.D0
        c(9) = FX(1,1)*x(12)+FX(2,1)*x(13)+FX(3,1)*x(14)-1.D0*x(12)+x(15
     #)+x(8)*(FXX(1,1,1)*x(5)+FXX(1,2,1)*x(6)+FXX(1,3,1)*x(7))+x(9)*(FXX
     #(2,1,1)*x(5)+FXX(2,2,1)*x(6)+FXX(2,3,1)*x(7))+x(10)*(FXX(3,1,1)*x(
     #5)+FXX(3,2,1)*x(6)+FXX(3,3,1)*x(7))
        c(10) = FX(1,2)*x(12)+FX(2,2)*x(13)+FX(3,2)*x(14)-1.D0*x(13)+x(8
     #)*(FXX(1,1,2)*x(5)+FXX(1,2,2)*x(6)+FXX(1,3,2)*x(7))+x(9)*(FXX(2,1,
     #2)*x(5)+FXX(2,2,2)*x(6)+FXX(2,3,2)*x(7))+x(10)*(FXX(3,1,2)*x(5)+FX
     #X(3,2,2)*x(6)+FXX(3,3,2)*x(7))
        c(11) = FX(1,3)*x(12)+FX(2,3)*x(13)+FX(3,3)*x(14)-1.D0*x(14)+x(8
     #)*(FXX(1,1,3)*x(5)+FXX(1,2,3)*x(6)+FXX(1,3,3)*x(7))+x(9)*(FXX(2,1,
     #3)*x(5)+FXX(2,2,3)*x(6)+FXX(2,3,3)*x(7))+x(10)*(FXX(3,1,3)*x(5)+FX
     #X(3,2,3)*x(6)+FXX(3,3,3)*x(7))
        c(12) = x(12)*F0(1)+x(13)*F0(2)+x(14)*F0(3)+(J0(1,1)*x(5)+J0(1,2
     #)*x(6)+J0(1,3)*x(7))*x(8)+(J0(2,1)*x(5)+J0(2,2)*x(6)+J0(2,3)*x(7))
     #*x(9)+(J0(3,1)*x(5)+J0(3,2)*x(6)+J0(3,3)*x(7))*x(10)
        c(13) = FP(1,1)*x(12)+FP(2,1)*x(13)+FP(3,1)*x(14)+x(8)*(x(5)*FXP
     #(1,1,1)+x(6)*FXP(1,2,1)+x(7)*FXP(1,3,1))+x(9)*(x(5)*FXP(2,1,1)+x(6
     #)*FXP(2,2,1)+x(7)*FXP(2,3,1))+x(10)*(x(5)*FXP(3,1,1)+x(6)*FXP(3,2,
     #1)+x(7)*FXP(3,3,1))-1.D0*x(16)
        c(14) = FP(1,2)*x(12)+FP(2,2)*x(13)+FP(3,2)*x(14)+x(8)*(x(5)*FXP
     #(1,1,2)+x(6)*FXP(1,2,2)+x(7)*FXP(1,3,2))+x(9)*(x(5)*FXP(2,1,2)+x(6
     #)*FXP(2,2,2)+x(7)*FXP(2,3,2))+x(10)*(x(5)*FXP(3,1,2)+x(6)*FXP(3,2,
     #2)+x(7)*FXP(3,3,2))-1.D0*x(17)
        c(15) = x(1)-1.D0*x(3)-1.D0*x(18)*x(16)/sqrt(x(16)**2+x(17)**2)
        c(16) = x(2)-1.D0*x(4)-1.D0*x(18)*x(17)/sqrt(x(16)**2+x(17)**2)
        c(17) = x(18)-0.2828427124D-1
        if (IFAIL .ne. 0) then
          write(6,*) 'FAILED IN CONSTRAINTS'
          call exit(0)
        endif
        return
      end

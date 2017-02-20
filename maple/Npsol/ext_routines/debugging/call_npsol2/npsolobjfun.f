      subroutine npsolobjfun(mode,n,x,objf,objgrad,nstate)
      integer mode
      integer n
      doubleprecision x(2)
      doubleprecision objf
      doubleprecision objgrad(2)
      integer nstate

      integer n_in_objfun
      doubleprecision z(1)

      common/global/par
      doubleprecision par(1)

        n_in_objfun = 2
        if ( .not. n_in_objfun .eq. n) then
          write(6,*) 'Error in confun: n differs from number of variable
     #s'
          call exit(0)
        endif
        z(1) = 0
        objf = (x(1)-2.D0)**2+(x(2)-1.D0)**2
        objgrad(1) = 2*x(1)-4
        objgrad(2) = 2*x(2)-2
        return
      end

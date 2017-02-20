      subroutine npsolconfun(mode,ncnln,n,nrowj,needc,x,c,cjac,nstate)
      integer mode
      integer ncnln
      integer n
      integer nrowj
      integer needc
      doubleprecision x(2)
      doubleprecision c(1)
      doubleprecision cjac(1,2)
      integer nstate

      integer n_in_confun
      integer ncnln_in_confun
      doubleprecision z(1)

      common/global/par
      doubleprecision par(1)

        n_in_confun = 2
        if ( .not. n_in_confun .eq. n) then
          write(6,*) 'Error in confun: n differs from number of variable
     #s'
          call exit(0)
        endif
        ncnln_in_confun = 1
        if ( .not. ncnln_in_confun .eq. ncnln) then
          write(6,*) 'Error in confun: ncnln differs from number of equa
     #tions'
          call exit(0)
        endif
        z(1) = 0
        c(1) = x(2)-1.D0*x(1)**2
        cjac(1,1) = -2*x(1)
        cjac(1,2) = 1
        return
      end

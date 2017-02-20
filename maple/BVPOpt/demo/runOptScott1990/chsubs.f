*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  chsubs.f
*
*     chcore   chfd     chfdlc   chfdls   chfgrd   chcJac   chfJac
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chcore( debug, done, first, epsa, epsr, fx,
     $                   inform, iter, itmax,
     $                   cdest, fdest, sdest, errbnd, f1,
     $                   f2, h, hopt, hphi )

      implicit           double precision (a-h,o-z)
      logical            debug, done, first

*     ==================================================================
*     chcore  implements algorithm  FD, the method described in
*     Gill, P.E., Murray, W., Saunders, M.A., and Wright, M. H.,
*     Computing Forward-Difference Intervals for Numerical Optimization,
*     Siam Journal on Scientific and Statistical Computing, vol. 4,
*     pp. 310-321, June 1983.
*     
*     The procedure is based on finding an interval (hphi) that
*     produces an acceptable estimate of the second derivative, and
*     then using that estimate to compute an interval that should
*     produce a reasonable forward-difference approximation.
*     
*     One-sided difference estimates are used to ensure feasibility with
*     respect to an upper or lower bound on x.  If x is close to an 
*     upper bound, the trial intervals will be negative.  The final 
*     interval is always positive.
*     
*     chcore has been designed to use a reverse communication
*     control structure, i.e., all evaluations of the function occur
*     outside this routine.  The calling routine repeatedly calls  
*     chcore  after computing the indicated function values.
*     
*     chcore  is similar to subroutine fdcore described in Report
*     SOL 83-6, Documentation of fdcore and fdcalc, by P.E. Gill,
*     W. Murray,  M.A. Saunders, and M.H. Wright, Department of
*     Operations Research,  Stanford University, Stanford, California
*     94305, June 1983.
*     
*     Systems Optimization Laboratory, Stanford University.
*     Department of Mathematics, University of California, San Diego.
*     Based on Fortran 66 Version 2.1 of  fdcore  written June 1983.
*     Fortran 77 Version written 25-May-1985.
*     This version of  chcore  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 

      logical            ce1big, ce2big, te2big, overfl
      save               cdsave, fdsave, hsave, oldh, rho, sdsave
      save               ce1big, ce2big, te2big

      parameter         (bndlo  =1.0d-3, bndup  =1.0d-1                )
      parameter         (zero   =0.0d+0, sixth  =1.6d-1, fourth =2.5d-1)
      parameter         (half   =5.0d-1, two    =2.0d+0, three  =3.0d+0)
      parameter         (four   =4.0d+0, ten    =1.0d+1                )

*     ------------------------------------------------------------------
*     Explanation of local variables...
*
*     bndlo, bndup, and rho control the logic of the routine.
*     bndlo  and  bndup are the lower and upper bounds that define an
*     acceptable value of the bound on the relative condition error in
*     the second derivative estimate.
*
*     The scalar rho is the factor by which the interval is multiplied
*     or divided, and also the multiple of the well-scaled interval
*     that is used as the initial trial interval.
*
*     All these values are discussed in the documentation.
*     ------------------------------------------------------------------

      iter  = iter + 1

*     Compute the forward-,  backward-,  central-  and second-order
*     difference estimates.

      fdest  = ddiv  ( f1 - fx,     h, overfl )
      fdest2 = ddiv  ( f2 - fx, two*h, overfl )

      oldcd = cdest
      cdest = ddiv  ( four*f1 - three*fx - f2, two*h, overfl )

      oldsd = sdest
      sdest = ddiv  ( fx      - two*f1   + f2, h*h  , overfl )

*     Compute  fdcerr  and  sdcerr,  bounds on the relative condition
*     errors in the first and second derivative estimates.

      afdmin = min( abs( fdest ), abs( fdest2 ) )
      fdcerr = ddiv  ( epsa, half*abs( h )*afdmin, overfl )
      sdcerr = ddiv  ( epsa, fourth*abs( sdest )*h*h, overfl )

      if (debug)
     $   write(iPrint, 9000) iter  , fx   , h,
     $                       f1    , fdest,
     $                       f2    , fdest2,
     $                       cdest , sdest,
     $                       fdcerr, sdcerr

*     ==================================================================
*     Select the correct case.
*     ==================================================================
      if (first) then
*        ---------------------------------------------------------------
*        First time through.
*        Check whether sdcerr lies in the acceptable range.
*        ------------------------------------------------------------
         first  = .false.
         done   = sdcerr .ge. bndlo  .and.  sdcerr .le. bndup
         te2big = sdcerr .lt. bndlo
         ce2big = sdcerr .gt. bndup
         ce1big = fdcerr .gt. bndup

         if (.not. ce1big) then
            hsave  = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest
         end if

         rho  = epsr**(-sixth)/four
         if (te2big) then

*           The truncation error may be too big  (same as saying
*           sdcerr is too small).  Decrease the trial interval.

            rho    = ten*rho
            oldh   = h
            h      = h / rho
         else if (ce2big) then

*           sdcerr is too large.  Increase the trial interval.

            oldh   = h
            h      = h*rho
         end if
      else if (ce2big) then
*        ---------------------------------------------------------------
*        During the last iteration,  the trial interval was
*        increased in order to decrease sdcerr.
*        ---------------------------------------------------------------
         if (ce1big  .and.  fdcerr .le. bndup) then
            ce1big = .false.
            hsave  = h
            fdsave = fdest
            cdsave = cdest
            sdsave = sdest
         end if

*        If sdcerr is small enough, accept h.  Otherwise,
*        increase h again.

         done   = sdcerr .le. bndup
         if (.not. done) then
            oldh   = h
            h      = h*rho
         end if
      else if (te2big) then
*        ---------------------------------------------------------------
*        During the last iteration,  the interval was decreased in order
*        to reduce the truncation error.
*        ---------------------------------------------------------------
         done   = sdcerr .gt. bndup
         if (done) then

*           sdcerr has jumped from being too small to being too
*           large.  Accept the previous value of h.

            h     = oldh
            sdest = oldsd
            cdest = oldcd
         else

*           Test whether fdcerr is sufficiently small.

            if (fdcerr .le. bndup) then
               ce1big = .false.
               hsave  = h
               fdsave = fdest
               cdsave = cdest
               sdsave = sdest
            end if

*           Check whether sdcerr is in range.

            done  = sdcerr .ge. bndlo

            if (.not. done) then

*              sdcerr is still too small, decrease h again.

               oldh = h
               h    = h / rho
            end if
         end if
      end if

*     ==================================================================
*     We have either finished or have a new estimate of h.
*     ==================================================================
      if (done) then

*        Sufficiently good second-derivative estimate found.
*        Compute the optimal interval.

         hphi   = abs( h )
         hopt   = two * sqrt( epsa ) / sqrt( abs( sdest ) )

*        err1 is the error bound on the forward-difference estimate
*        with the final value of h.  err2 is the difference of fdest
*        and the central-difference estimate with hphi.

         err1   = hopt*abs( sdest )
         err2   = abs( fdest - cdest )
         errbnd = max( err1, err2 )

*        Set inform = 4  if the forward- and central-difference
*        estimates are not close.

         inform = 0
         if (errbnd .gt. half*abs( fdest )) inform = 4
      else
*        ---------------------------------------------------------------
*        Check whether the maximum number of iterations has been
*        exceeded.  If not, exit.
*        ---------------------------------------------------------------
         done = iter .ge. itmax
         if (done) then
            if (ce1big) then

*              fdcerr was never small.  Probably a constant function.

               inform = 1
               hphi   = hopt
               fdest  = zero
               cdest  = zero
               sdest  = zero
               errbnd = zero
            else if (ce2big) then

*              fdcerr was small,  but sdcerr was never small.
*              Probably a linear or odd function.

               inform = 2
               hphi   = abs( hsave )
               hopt   = hphi
               fdest  = fdsave
               cdest  = cdsave
               sdest  = zero
               errbnd = two*epsa / hopt
            else

*              The only remaining case occurs when the second
*              derivative is changing too rapidly for an adequate
*              interval to be found (sdcerr remained small even
*              though h was decreased itmax times).

               inform = 3
               hphi   = abs( hsave )
               hopt   = hphi
               fdest  = fdsave
               cdest  = cdsave
               sdest  = sdsave
               errbnd = hopt*abs( sdest )/two + two*epsa/hopt
            end if
         end if
      end if

      if (debug) then
         write(iPrint, 9001) ce1big, ce2big, te2big
         if (done)
     $      write(iPrint, 9002) inform, hopt, errbnd
      end if

      return

 9000 format(/ ' //chcore//  Itn ', i3,
     $                             ' fx     h      ', 5x, 1p, 2d16.6
     $       / ' //chcore//  f1      fdest         ', 5x, 1p, 2d16.6
     $       / ' //chcore//  f2      fdest2        ', 5x, 1p, 2d16.6
     $       / ' //chcore//  cdest   sdest         ', 5x, 1p, 2d16.6
     $       / ' //chcore//  fdcerr  sdcerr        ', 5x, 1p, 2d16.6)
 9001 format(  ' //chcore//  ce1big  ce2big  te2big', 5x, 3l2       )
 9002 format(  ' //chcore//  inform  hopt    errbnd', i5, 1p, 2d16.6)

*     end of chcore
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chfd  ( inform, msglvl, lvlder,
     $                   n, ncnln, ldcJ, ldcJu,
     $                   bigbnd, epsrf, fdnorm, objf,
     $                   funobj, funcon, needc,
     $                   bl, bu, c, c1, c2, cJac, cJacu,
     $                   grad, gradu, hforwd, hcntrl, x, y )

      implicit           double precision(a-h,o-z)
      integer            needc(*)
      double precision   bl(n), bu(n)
      double precision   c(*), c1(*), c2(*),
     $                   cJac(ldcJ,*), cJacu(ldcJu,*)
      double precision   grad(n), gradu(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   x(n), y(n)
      external           funobj, funcon

*     ==================================================================
*     chfd    computes difference intervals for the missing gradients of
*     f(x) and c(x).  Intervals are computed using a procedure that
*     usually requires about two function evaluations if the function is
*     well scaled.  Central-difference gradients are obtained as a
*     by-product of the computation.
*
*     On entry...
*        objf and c contain the problem functions at the point x.
*        An element of cJac or grad not equal to rdummy signifies a
*        known gradient value.  Such values are not estimated by
*        differencing.  cJacu and gradu have dummy elements in the same
*        positions as  cJac and gradu.
*
*     On exit...                          
*        cJac and grad contain central-difference derivative estimates.
*        Elements of cJacu and gradu are unaltered except for those
*        corresponding to constant derivatives, which are given the same
*        values as cJac or grad.
*
*     Systems Optimization Laboratory, Department of Operations Research
*     Department of Mathematics, University of California, San Diego.
*     Original version written 28-July-1985.
*     This version of  chfd  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      logical            debug , done  , first , prtHdr, needed
      parameter         (rdummy =-11111.0d+0)
      parameter         (factor = 0.97d+0   )
      parameter         (zero   = 0.0d+0, half   =0.5d+0, one   =1.0d+0)
      parameter         (two    = 2.0d+0, four   =4.0d+0, ten   =1.0d+1)

      inform = 0
      needed = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or.  lvlder .eq. 1  .and.  ncnln .gt. 0
      if (.not. needed) return

      debug  = .false. 
      if (lfdset .eq. 0) then
         if (msglvl .gt. 0) then
            if (iPrint .gt. 0) write(iPrint, 1000)
         end if
         nstate = 0
         itmax  = 3
         mode   = 0

         nccnst = 0
         nfcnst = 0
         prtHdr = .true.

         fdnorm = zero

*        ===============================================================
*        For each column of the Jacobian augmented by the transpose of
*        the objective gradient, rows irow1 thru irow2 are searched for
*        missing elements.
*        ===============================================================
         irow1  = 1
         irow2  = ncnln + 1
         if (lvlder .eq. 1) irow2 = ncnln
         if (lvlder .eq. 2) irow1 = ncnln + 1

         biglow = - bigbnd
         bigupp =   bigbnd

         if (ncnln  .gt. 0)
     $      call iload ( ncnln, 0, needc, 1 )

         do 600, j = 1, n
            xj     = x(j)
            nColj  = 0
            sumsd  = zero
            sumeps = zero
            hfd    = zero
            hcd    = zero
            hmax   = zero
            hmin   = one / epspt3
            errmax = zero
            errmin = zero

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            do 500, i  = irow1, irow2
               if (i .le. ncnln) then
                  test = cJacu(i,j)
               else
                  test = gradu(j)
               end if

               if (test .eq. rdummy) then
*                 ======================================================
*                 Get the difference interval for this component.
*                 ======================================================
                  nColj = nColj + 1

                  if (i .le. ncnln) then
                     needc(i) = 1
                     fx       = c(i)
                     epsa     = epsrf*(one + abs( c(i) ))
                  else
                     fx       = objf
                     epsa     = epsrf*(one + abs( fx ))
                  end if

*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  iter   = 0
                  hopt   = two*(one + abs( xj ))*sqrt( epsrf )
                  h      = signh*ten*hopt
                  cdest  = zero
                  sdest  = zero
                  first  = .true.

*+                repeat
  400                x(j)  = xj + h
                     if (i .le. ncnln) then
                        call funcon( mode, ncnln, n, ldcJu,
     $                               needc, x, c1, cJacu, nstate )
                        if (mode .lt. 0) go to 9999
                        f1 = c1(i)
                     else
                        call funobj( mode, n, x, f1, gradu, nstate )
                        if (mode .lt. 0) go to 9999
                     end if

                     x(j)  = xj + h + h
                    if (i .le. ncnln) then
                       call funcon( mode, ncnln, n, ldcJu,
     $                              needc, x, c1, cJacu, nstate )
                        if (mode .lt. 0) go to 9999
                        f2 = c1(i)
                     else
                        call funobj( mode, n, x, f2, gradu, nstate )
                        if (mode .lt. 0) go to 9999
                     end if

                     call chcore( debug, done, first, epsa, epsrf, fx,
     $                            info, iter, itmax,
     $                            cdest, fdest, sdest, errbnd, f1,
     $                            f2, h, hopt, hphi )
      
*+                until     done
                  if (.not. done) go to 400

                  if (i .le. ncnln) then
                     cJac(i,j) = cdest
                     if (info .eq. 1  .or.  info .eq. 2) then
                        nccnst    =   nccnst + 1
                        ncdiff    =   ncdiff - 1
                        cJacu(i,j) = - rdummy
                     end if
                  else
                     grad(j)   = cdest
                     if (info .eq. 1  .or.  info .eq. 2) then
                        nfcnst    =   nfcnst + 1
                        nfdiff    =   nfdiff - 1
                        gradu(j)  = - rdummy
                     end if
                  end if

                  sumsd  = sumsd  + abs( sdest )
                  sumeps = sumeps +      epsa
                  if (hopt .gt. hmax) then
                     hmax   = hopt
                     errmax = errbnd
                  end if
                  if (hopt .lt. hmin) then
                     hmin   = hopt
                     errmin = errbnd
                  end if

                  if (info .eq. 0) hcd  = max ( hcd, hphi )
               end if
  500       continue

            if (nColj .gt. 0) then
               if (hmin .gt. hmax) then
                  hmin   = hmax
                  errmin = errmax
               end if

               if      (four*sumeps .lt. hmin*hmin*sumsd) then
                  hfd    = hmin
                  errmax = errmin
               else if (four*sumeps .gt. hmax*hmax*sumsd) then
                  hfd    = hmax
               else
                  hfd    = two*sqrt( sumeps / sumsd )
                  errmax = two*sqrt( sumeps * sumsd )
               end if

               if (hcd .eq. zero) hcd = ten*hfd

               if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
                  if (prtHdr) write(iPrint, 1100)
                  write(iPrint, 1200) j, xj, hfd, hcd, errmax
                  prtHdr = .false.
               end if
               fdnorm    = max (fdnorm, hfd)
               hforwd(j) = hfd / (one + abs(xj))
               hcntrl(j) = hcd / (one + abs(xj))
            end if
            x(j)      = xj
  600    continue

         if (nccnst + nfcnst .gt. 0) then

*           Check that the constants have been set properly by
*           evaluating the gradients at a strange (but feasible) point.

            d      =   one / n

            do 710, j = 1, n
               xj     =   x(j)
               stepbl = - one
               stepbu =   one
               if (bl(j) .gt. biglow)
     $            stepbl = max( stepbl, bl(j) - xj )
               if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $            stepbu = min( stepbu, bu(j) - xj )

               if (half*(stepbl + stepbu) .lt. zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if

               d = factor*d
  710       continue

            if (ncnln .gt. 0) then
               call iload ( ncnln, 1, needc, 1 )
               call funcon( mode, ncnln, n, ldcJu,
     $                      needc, y, c2, cJacu, nstate )
               if (mode .lt. 0) go to 9999
            end if

            call funobj( mode, n, y, objf2, gradu, nstate )
            if (mode .lt. 0) go to 9999

*           ------------------------------------------------------------
*           Loop over each of the components of  x.
*           ------------------------------------------------------------
            do 800, j = 1, n
               yj     = y(j)
               dx     = half*(x(j) - yj)
               y(j)   = yj + dx

               if (ncnln .gt. 0) then
                  nColj    = 0
                  do 720 i = 1, ncnln
                     if (cJacu(i,j) .eq. - rdummy) then
                        needc(i) = 1
                        nColj    = nColj + 1
                     else
                        needc(i) = 0
                     end if
  720             continue

                  if (nColj .gt. 0) then
                     call funcon( mode, ncnln, n, ldcJu,
     $                            needc, y, c1, cJacu, nstate )
                     if (mode .lt. 0) go to 9999

                     do 730 i = 1, ncnln
                        if (needc(i) .eq. 1) then
                           cJdiff = ( c1(i) -  c2(i) ) / dx
                           if (cJdiff .eq. cJac(i,j)) then
                              cJacu(i,j) = cJdiff
                           else
                              cJacu(i,j) = rdummy
                              nccnst    = nccnst - 1
                              ncdiff    = ncdiff + 1
                           end if
                        end if
  730                continue
                  end if
               end if

*              Now check the objective gradient component.

               if (gradu(j) .eq. - rdummy) then

                  call funobj( mode, n, y, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  gdiff = (f1 - objf2)/dx
                  if (gdiff .eq. grad(j)) then
                     gradu(j) = gdiff
                  else
                     gradu(j) = rdummy
                     nfdiff   = nfdiff + 1
                     nfcnst   = nfcnst - 1
                  end if
               end if

               y(j)  = yj
  800       continue

            if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
               if (lvlder .lt. 2  .and.  nccnst .gt. 0)
     $            write(iPrint, 1300) nccnst
               if (lvlder .ne. 1  .and.  nfcnst .gt. 0)
     $            write(iPrint, 1400) nfcnst
            end if

            if (ncnln .gt. 0) then
               if (ncdiff .eq. 0  .and.  lvlder .lt. 2) then
                  if (lvlder .eq. 0) lvlder = 2
                  if (lvlder .eq. 1) lvlder = 3
                  if (msglvl .gt. 0) then
                     if (iPrint .gt. 0) write(iPrint, 1500) lvlder
                  end if
               end if
            end if

            if (nfdiff .eq. 0  .and.  lvlder .ne. 1) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msglvl .gt. 0) then
                  if (iPrint .gt. 0) write(iPrint, 1600) lvlder
               end if
            end if
         end if
      else if (lfdset .eq. 2) then

*        The user has supplied hforwd and hcntrl.
*        Check for wild values.

         do 900, j = 1, n
            if (hforwd(j) .le. zero) then
               if (iPrint .gt. 0) 
     $            write(iPrint, 2000) j, hforwd(j), epspt5
               if (iSumm  .gt. 0) 
     $            write(iSumm , 2000) j, hforwd(j), epspt5
               hforwd(j) = epspt5
            end if
  900    continue
         do 910 j = 1, n
            if (hcntrl(j) .le. zero) then
               if (iPrint .gt. 0) 
     $            write(iPrint, 2100) j, hcntrl(j), epspt3
               if (iSumm  .gt. 0) 
     $            write(iSumm , 2100) j, hcntrl(j), epspt3
               hcntrl(j) = epspt3
            end if
  910    continue
      end if

      return

 9999 inform = mode
      return

 1000 format(//' Computation of the finite-difference intervals'
     $       / ' ----------------------------------------------' )
 1100 format(//'    j      x(j)   Forward dx(j)   Central dx(j) ',
     $         '     Error est.' /)
 1200 format(  i5, 1p, e10.2, e16.6, 2e16.6 )
 1300 format(/ i5,  '  constant constraint gradient elements assigned.')
 1400 format(/ i5,  '  constant  objective gradient elements assigned.')
 1500 format(//' All missing Jacobian elements are constants.  ',
     $         ' Derivative level increased to ', i4 )
 1600 format(//' All missing objective gradients are constants.  ',
     $         ' Derivative level increased to ', i4 )
 2000 format(' XXX  ', i4,'-th difference interval ',
     $          1p, e10.2, ' replaced by ', e10.2 )
 2100 format(' XXX  ', i4,'-th central-difference interval ',
     $          1p, e10.2, ' replaced by ', e10.2 )

*     end of chfd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chfdlc( inform, msglvl, lvlder, n, 
     $                   bigbnd, epsrf, fdnorm, objf,
     $                   funobj, bl, bu, 
     $                   grad, gradu, hforwd, hcntrl, x, y )

      implicit           double precision(a-h,o-z)
      double precision   bl(n), bu(n)
      double precision   grad(n), gradu(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   x(n), y(n)
      external           funobj

*     ==================================================================
*     chfdlc  computes difference intervals for the missing gradients of
*     f(x).  Difference intervals are computed using a procedure that 
*     usually requires about two evaluations if  f(x)  is well scaled.
*     Central-difference gradients are obtained as a by-product of the
*     algorithm.
*
*     On entry...
*       objf  contains the objective function at the point x.
*       An element of  grad  not equal to  rdummy  signifies a user-
*       specified value.  Such values are not estimated by differencing.
*       grad  has dummy elements in the same positions as gradu.
*
*     On exit...                          
*       grad  contain central-difference derivative estimates.
*       Elements of  gradu  are unaltered except for those corresponding 
*       to constant derivatives, which are given values values found
*       by differencing.
*
*     Systems Optimization Laboratory, Stanford University.
*     Department of Mathematics, University of California, San Diego.
*     Original version written  5-July-1990.
*     This version of  chfdlc  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      logical            debug , done  , first , prtHdr, needed
      parameter         (rdummy =-11111.0d+0)
      parameter         (factor = 0.97d+0   )
      parameter         (zero   = 0.0d+0, half =  0.5d+0, one = 1.0d+0)
      parameter         (two    = 2.0d+0, ten  = 10.0d+0)

      inform = 0
      needed = lvlder .eq. 0  .or.  lvlder .eq. 2
      if (.not. needed) return

      debug  = .false.
      if (lfdset .eq. 0) then
         if (msglvl .gt. 0  .and.  iPrint .gt. 0) write(iPrint, 1000)

         nstate = 0
         itmax  = 3
         mode   = 0

         nfcnst = 0
         prtHdr = .true.

         fdnorm = zero

         biglow = - bigbnd
         bigupp =   bigbnd

         do 600, j = 1, n
            xj     = x(j)

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            if (gradu(j) .eq. rdummy) then
*              =========================================================
*              This component needs a difference interval.
*              =========================================================
               epsa     = epsrf*(one + abs( objf ))
               iter   = 0
               hopt   = two*(one + abs( xj ))*sqrt( epsrf )
               h      = signh*ten*hopt
               cdest  = zero
               sdest  = zero
               first  = .true.

*+             repeat
  400             x(j)  = xj + h
                  call funobj( mode, n, x, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  x(j)  = xj + h + h
                  call funobj( mode, n, x, f2, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  call chcore( debug, done, first, epsa, epsrf, objf,
     $                         info, iter, itmax,
     $                         cdest, fdest, sdest, errbnd, f1,
     $                         f2, h, hopt, hphi )
*+             until     done
               if (.not. done) go to 400

               grad(j) = cdest
               if (info .eq. 1  .or.  info .eq. 2) then
                  nfcnst   =   nfcnst + 1
                  nfdiff   =   nfdiff - 1
                  gradu(j) = - rdummy
               end if

               if (info .eq. 0) then
                  hfd = hopt
                  hcd = hphi
               else 
                  hfd = two*(one + abs( xj ))*sqrt( epsrf )
                  hcd = ten*hfd
               end if

               if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
                  if (prtHdr) write(iPrint, 1100)
                  write(iPrint, 1200) j, xj, hfd, hcd, errbnd
                  prtHdr = .false.
               end if

               fdnorm    = max (fdnorm, hfd)
               hforwd(j) = hfd / (one + abs(xj))
               hcntrl(j) = hcd / (one + abs(xj))
            end if
            x(j)      = xj
  600    continue

         if (nfcnst .gt. 0) then

*           Check that the constant elements of  g  have been identified
*           correctly.  Re-evaluate the gradients at a strange feasible
*           point.

            d      =   one / n

            do 710, j = 1, n
               xj     =   x(j)
               stepbl = - one
               stepbu =   one
               if (bl(j) .gt. biglow)
     $            stepbl = max( stepbl, bl(j) - xj )
               if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $            stepbu = min( stepbu, bu(j) - xj )

               if (half*(stepbl + stepbu) .lt. zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if

               d = factor*d
  710       continue

            call funobj( mode, n, y, objf2, gradu, nstate )
            if (mode .lt. 0) go to 9999

*           ------------------------------------------------------------
*           Loop over each component of  x.
*           ------------------------------------------------------------
            do 800, j = 1, n
               yj     = y(j)
               dx     = half*(x(j) - yj)
               y(j)   = yj + dx

               if (gradu(j) .eq. - rdummy) then

                  call funobj( mode, n, y, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  gdiff = (f1 - objf2)/dx
                  if (gdiff .eq. grad(j)) then
                     gradu(j) = gdiff
                  else
                     gradu(j) = rdummy
                     nfdiff   = nfdiff + 1
                     nfcnst   = nfcnst - 1
                  end if
               end if

               y(j)  = yj
  800       continue

            if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
               if (nfcnst .gt. 0)
     $            write(iPrint, 1400) nfcnst
            end if

            if (nfdiff .eq. 0) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msglvl .gt. 0) then
                  if (iPrint .gt. 0) write(iPrint, 1600) lvlder
               end if
            end if
         end if
      else if (lfdset .eq. 2) then

*        The user has supplied  hforwd  and  hcntrl.
*        Check for wild values.

         do 900, j = 1, n
            if (hforwd(j) .le. zero) then
               if (iPrint .gt. 0) 
     $            write(iPrint, 2000) j, hforwd(j), epspt5
               if (iSumm  .gt. 0) 
     $            write(iSumm , 2000) j, hforwd(j), epspt5
               hforwd(j) = epspt5
            end if
  900    continue
         do 910, j = 1, n
            if (hcntrl(j) .le. zero) then
               if (iPrint .gt. 0) 
     $            write(iPrint, 2100) j, hcntrl(j), epspt3
               if (iSumm  .gt. 0) 
     $            write(iSumm , 2100) j, hcntrl(j), epspt3
               hcntrl(j) = epspt3
            end if
  910    continue
      end if

      return

 9999 inform = mode

 1000 format(//' Computation of the finite-difference intervals'
     $       / ' ----------------------------------------------' )
 1100 format(//'    j      x(j)   Forward dx(j)   Central dx(j) ',
     $         '     Error est.' /)
 1200 format(  i5, 1p, e10.2, e16.6, 2e16.6 )
 1400 format(/ I5,  '  constant  objective gradient elements assigned.')
 1600 format(//' All missing objective gradients are constants.  ',
     $         ' Derivative level increased to ', i4 )
 2000 format(' XXX  ', i4,'-th difference interval ',         1pe10.2,
     $       ' replaced by ', 1pe10.2 )
 2100 format(' XXX  ', i4,'-th central-difference interval ', 1pe10.2,
     $       ' replaced by ', 1pe10.2 )

*     end of chfdlc
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chfdls( inform, msglvl, lvlder,
     $                   m, n, ncnln, ldcJ, ldcJu, ldfJ, ldfJu,
     $                   bigbnd, epsrf, fdnorm,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, c2, cJac, cJacu,
     $                   f, f1, f2, fJac, fJacu, hforwd, hcntrl,
     $                   x, y )

      implicit           double precision(a-h,o-z)
      integer            needc(*)
      double precision   bl(n), bu(n)
      double precision   c(*), c1(*), c2(*),
     $                   cJac(ldcJ,*), cJacu(ldcJu,*)
      double precision   f(m), f1(m), f2(m),
     $                   fJac(ldfJ,*), fJacu(ldfJu,*)
      double precision   hforwd(*), hcntrl(*)
      double precision   x(n), y(n)
      external           funobj, funcon

*     ==================================================================
*     chfdls  computes difference intervals for missing Jacobian
*     elements of f(x) and c(x).  Intervals are computed using a
*     procedure that requires about two function evaluations if the
*     function is well scaled.  Central-difference gradients are
*     obtained as a by-product of the computation.
*
*     On entry...
*        f and c contain the problem functions at the point x.
*        An element of cJac or fJac not equal to rdummy signifies a known
*        gradient value.  Such values are not estimated by differencing.
*        cJacu and fJacu have dummy elements in the same positions as
*        cJac and fJac.
*
*     On exit...
*        cJac and fJac contain central-difference derivative estimates.
*        Elements of cJacu and fJacu are unaltered except for those
*        corresponding to constant derivatives, which are given the same
*        values as cJac or fJac.
*
*     Systems Optimization Laboratory, Stanford University.
*     Department of Mathematics, University of California, San Diego.
*     Original version written 10-May-1988.
*     This version of  chfdls  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      logical            debug , done  , first , prtHdr, needed
      parameter         (rdummy =-11111.0d+0)
      parameter         (factor = 0.97d+0    )
      parameter         (zero   = 0.0d+0, half = 0.5d+0, one =  1.0d+0)
      parameter         (two    = 2.0d+0, four = 4.0d+0, ten = 10.0d+0)

      inform = 0
      needed = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or.  lvlder .eq. 1  .and.  ncnln .gt. 0
      if (.not. needed) return

      debug  = .false.
      if (lfdset .eq. 0) then
         if (msglvl .gt. 0  .and.  iPrint .gt. 0) write(iPrint, 1000)

         nstate = 0
         itmax  = 3
         mode   = 0

         nccnst = 0
         nfcnst = 0
         prtHdr = .true.

         fdnorm = zero

*        ===============================================================
*        For each column of the matrix
*              ( cJac  )
*              ( fJac  ),
*        rows irow1 thru irow2 are searched for missing elements.
*        ===============================================================
         irow1  = 1
         irow2  = ncnln + m
         if (lvlder .eq. 1) irow2 = ncnln
         if (lvlder .eq. 2) irow1 = ncnln + 1

         biglow = - bigbnd
         bigupp =   bigbnd

         if (ncnln  .gt. 0)
     $      call iload ( ncnln, 0, needc, 1 )

         do 600, j = 1, n
            xj     = x(j)
            nColj  = 0
            sumsd  = zero
            sumeps = zero
            hfd    = zero
            hcd    = zero
            hmax   = zero
            hmin   = one / epspt3
            errmax = zero
            errmin = zero

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            do 500, irow = irow1, irow2
               if (irow .le. ncnln) then
                  i    = irow
                  test = cJacu(i,j)
               else
                  i    = irow - ncnln
                  test = fJacu(i,j)
               end if

               if (test .eq. rdummy) then
*                 ======================================================
*                 Get the difference interval for this component.
*                 ======================================================
                  nColj = nColj + 1

                  if (irow .le. ncnln) then
                     needc(i) = 1
                     fx       = c(i)
                  else
                     fx       = f(i)
                  end if

                  epsa     = epsrf*(one + abs( fx ))

*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  iter   = 0
                  hopt   = two*(one + abs( xj ))*sqrt( epsrf )
                  h      = signh*ten*hopt
                  cdest  = zero
                  sdest  = zero
                  first  = .true.

*+                repeat
  400                x(j)  = xj + h
                     if (irow .le. ncnln) then
                        call funcon( mode, ncnln, n, ldcJu,
     $                               needc, x, c1, cJacu, nstate )
                        if (mode .lt. 0) go to 9999
                        fforw = c1(i)
                     else
                        call funobj( mode, m    , n, ldfJu,
     $                                      x, f1, fJacu, nstate )
                        if (mode .lt. 0) go to 9999
                        fforw = f1(i)
                     end if

                     x(j)  = xj + h + h
                    if (irow .le. ncnln) then
                       call funcon( mode, ncnln, n, ldcJu,
     $                              needc, x, c1, cJacu, nstate )
                        if (mode .lt. 0) go to 9999
                        fback = c1(i)
                     else
                        call funobj( mode, m    , n, ldfJu,
     $                                      x, f1, fJacu, nstate )
                        if (mode .lt. 0) go to 9999
                        fback = f1(i)
                     end if

                     call chcore( debug, done, first, epsa, epsrf, fx,
     $                            info, iter, itmax,
     $                            cdest, fdest, sdest, errbnd, fforw,
     $                            fback, h, hopt, hphi )

*+                until     done
                  if (.not. done) go to 400

                  if (irow .le. ncnln) then
                     cJac(i,j) = cdest
                     if (info .eq. 1  .or.  info .eq. 2) then
                        nccnst     =   nccnst + 1
                        ncdiff     =   ncdiff - 1
                        cJacu(i,j) = - rdummy
                     end if
                  else
                     fJac(i,j) = cdest
                     if (info .eq. 1  .or.  info .eq. 2) then
                        nfcnst     =   nfcnst + 1
                        nfdiff     =   nfdiff - 1
                        fJacu(i,j) = - rdummy
                     end if
                  end if

                  sumsd  = sumsd  + abs( sdest )
                  sumeps = sumeps +      epsa
                  if (hopt .gt. hmax) then
                     hmax   = hopt
                     errmax = errbnd
                  end if
                  if (hopt .lt. hmin) then
                     hmin   = hopt
                     errmin = errbnd
                  end if

                  if (info .eq. 0) hcd  = max ( hcd, hphi )
               end if
  500       continue

            if (nColj .gt. 0) then
               if (hmin .gt. hmax) then
                  hmin   = hmax
                  errmin = errmax
               end if

               if      (four*sumeps .lt. hmin*hmin*sumsd) then
                  hfd    = hmin
                  errmax = errmin
               else if (four*sumeps .gt. hmax*hmax*sumsd) then
                  hfd    = hmax
               else
                  hfd    = two*sqrt( sumeps / sumsd )
                  errmax = two*sqrt( sumeps * sumsd )
               end if

               if (hcd .eq. zero) hcd = ten*hfd

               if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
                  if (prtHdr) write(iPrint, 1100)
                  write(iPrint, 1200) j, xj, hfd, hcd, errmax
                  prtHdr = .false.
               end if
               fdnorm    = max (fdnorm, hfd)
               hforwd(j) = hfd / (one + abs(xj))
               hcntrl(j) = hcd / (one + abs(xj))
            end if

            x(j)      = xj
  600    continue

         if (nccnst + nfcnst .gt. 0) then

*           Check that the constants have been set properly by
*           evaluating the gradients at a strange (but feasible) point.

            d      =   one / n

            do 710, j = 1, n
               xj     =   x(j)
               stepbl = - one
               stepbu =   one
               if (bl(j) .gt. biglow)
     $            stepbl = max( stepbl, bl(j) - xj )
               if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $            stepbu = min( stepbu, bu(j) - xj )

               if (half*(stepbl + stepbu) .lt. zero) then
                  y(j) = xj + d*stepbl
               else
                  y(j) = xj + d*stepbu
               end if

               d = factor*d
  710       continue

            if (ncnln .gt. 0) then
               call iload ( ncnln, 1, needc, 1 )
               call funcon( mode, ncnln, n, ldcJu,
     $                      needc, y, c2, cJacu, nstate )
               if (mode .lt. 0) go to 9999
            end if

            if (m     .gt. 0) then
               call funobj( mode, m, n, ldfJu, y, f2, fJacu, nstate )
               if (mode .lt. 0) go to 9999
            end if

*           ------------------------------------------------------------
*           Loop over each of the components of  x.
*           ------------------------------------------------------------
            do 800, j = 1, n
               yj     = y(j)
               dx     = half*(x(j) - yj)
               y(j)   = yj + dx

               if (ncnln .gt. 0) then
                  nColj     = 0
                  do 720, i = 1, ncnln
                     if (cJacu(i,j) .eq. - rdummy) then
                        needc(i) = 1
                        nColj    = nColj + 1
                     else
                        needc(i) = 0
                     end if
  720             continue

                  if (nColj .gt. 0) then
                     call funcon( mode, ncnln, n, ldcJu,
     $                            needc, y, c1, cJacu, nstate )
                     if (mode .lt. 0) go to 9999

                     do 730, i = 1, ncnln
                        if (needc(i) .eq. 1) then
                           cJdiff = ( c1(i) -  c2(i) ) / dx
                           if (cJdiff .eq. cJac(i,j)) then
                              cJacu(i,j) = cJdiff
                           else
                              cJacu(i,j) = rdummy
                              nccnst    = nccnst - 1
                              ncdiff    = ncdiff + 1
                           end if
                        end if
  730                continue
                  end if
               end if

*              Now check the objective Jacobian.

               nColj = 0
               do 740, i = 1, m
                  if (fJacu(i,j) .eq. - rdummy) then
                     nColj = nColj + 1
                  end if
  740          continue

               if (nColj .gt. 0) then
                  call funobj( mode, m, n, ldfJu,
     $                         y, f1, fJacu, nstate )
                  if (mode .lt. 0) go to 9999

                  do 750, i = 1, m
                     if (fJacu(i,j) .eq. - rdummy) then
                        fjdiff = ( f1(i) -  f2(i) ) / dx
                        if (fjdiff .eq. fJac(i,j)) then
                           fJacu(i,j) = fjdiff
                        else
                           fJacu(i,j) = rdummy
                           nfcnst     = nfcnst - 1
                           nfdiff     = nfdiff + 1
                        end if
                     end if
  750             continue
               end if

               y(j)  = yj
  800       continue

            if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
               if (lvlder .lt. 2  .and.  nccnst .gt. 0)
     $            write(iPrint, 1300) nccnst
               if (lvlder .ne. 1  .and.  nfcnst .gt. 0)
     $            write(iPrint, 1400) nfcnst
            end if

            if (ncnln .gt. 0) then
               if (ncdiff .eq. 0  .and.  lvlder .lt. 2) then
                  if (lvlder .eq. 0) lvlder = 2
                  if (lvlder .eq. 1) lvlder = 3
                  if (msglvl .gt. 0) then
                     if (iPrint .gt. 0) write(iPrint, 1500) lvlder
                  end if
               end if
            end if

            if (nfdiff .eq. 0  .and.  lvlder .ne. 1) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msglvl .gt. 0) then
                  if (iPrint .gt. 0) write(iPrint, 1600) lvlder
               end if
            end if
         end if
      else if (lfdset .eq. 2) then

*        The user has supplied hforwd and hcntrl.
*        Check for wild values.

         do 900, j = 1, n
            if (hforwd(j) .le. zero) then
               if (iPrint .gt. 0) 
     $            write(iPrint, 2000) j, hforwd(j), epspt5
               if (iSumm  .gt. 0) 
     $            write(iSumm , 2000) j, hforwd(j), epspt5
               hforwd(j) = epspt5
            end if
  900    continue
         do 910, j = 1, n
            if (hcntrl(j) .le. zero) then
               if (iPrint .gt. 0) 
     $            write(iPrint, 2100) j, hcntrl(j), epspt3
               if (iSumm  .gt. 0) 
     $            write(iSumm , 2100) j, hcntrl(j), epspt3
               hcntrl(j) = epspt3
            end if
  910    continue
      end if

      return

 9999 inform = mode
      return

 1000 format(//' Computation of the finite-difference intervals'
     $       / ' ----------------------------------------------' )
 1100 format(//'    j      x(j)   Forward dx(j)   Central dx(j) ',
     $         '     Error est.' /)
 1200 format(  i5, 1p, e10.2, e16.6, 2e16.6 )
 1300 format(/ i5,  '  constant constraint gradient elements assigned.')
 1400 format(/ i5,  '  constant objective  gradient elements assigned.')
 1500 format(//' All missing constraint Jacobian elements',
     $         ' are constants.   Derivative level increased to ', i4 )
 1600 format(//' All missing objective  Jacobian elements',
     $         ' are constants.   Derivative level increased to ', i4 )
 2000 format(' XXX  ', i4,'-th difference interval ',
     $       1p, e10.2, ' replaced by ', e10.2 )
 2100 format(' XXX  ', i4,'-th central-difference interval ',
     $       1p, e10.2, ' replaced by ', e10.2 )

*     end of  chfdls
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        
      subroutine chfgrd( inform, msglvl, n,
     $                   bigbnd, epsrf, oktol, fdchk, objf, xnorm,
     $                   funobj,
     $                   bl, bu, grad, gradu, dx, x, y )

      implicit           double precision(a-h,o-z)

      double precision   bl(n), bu(n), grad(n), gradu(n), dx(n)
      double precision   x(n), y(n)
      external           funobj

*     ==================================================================
*     chfgrd  checks if the gradients of the objective function have
*     been coded correctly.
*
*     On input,  the value of the objective function at the point x is
*     stored in objf.  The corresponding gradient is stored in gradu.
*     If any gradient component has not been specified,  it will have a
*     dummy value.  Missing values are not checked.
*
*     A cheap test is first undertaken by calculating the directional
*     derivative using two different methods.  If this proves 
*     satisfactory and no further information is desired, chfgrd is 
*     terminated.  Otherwise, the routine chcore is called to give 
*     optimal step-sizes and a forward-difference approximation to each 
*     component of the gradient for which a test is deemed necessary,
*     either by the program or the user.
*
*     Other inputs:
*
*           x         The n-dimensional point at which the
*                     gradient is to be verified.
*           epsrf     The positive bound on the relative error
*                     associated with computing the function at
*                     the point x.
*           oktol     The desired relative accuracy which the
*                     components of the gradient should satisfy.
* 
*     lvrfyc has the following meaning...
*
*       -1        do not perform any check.
*        0        do the cheap test only.
*        1 or 3   do both cheap and full test.
*
*     Systems Optimization Laboratory, Stanford University.
*     Department of Mathematics, University of California, San Diego.
*     Original version written  19-May-1985.
*     This version of  chfgrd  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5np/ lvrfyc, jverfy(4)

      logical            const , debug , done  , first , prtHdr
      logical            needed, okay
      character*4        key   , lbad  , lgood
      character*18       result(0:4)
      parameter         (rdummy =-11111.0d+0              )
      parameter         (zero   =0.0d+0, half  = 0.5d+0, point9 =0.9d+0)
      parameter         (one    =1.0d+0, two   = 2.0d+0, ten    =1.0d+1)
      parameter         (lbad   ='BAD?', lgood = '  OK')
      data               result
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      inform = 0
      needed = lvrfyc .eq. 0  .or.  lvrfyc .eq. 1  .or.  lvrfyc .eq. 3
      if (.not. needed) return
 
      if (msglvl .gt. 0  .and.  iPrint .gt. 0) write(iPrint, 1000)
      debug  = .false. 
      nstate = 0

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      h =     (one + xnorm)*fdchk

      if (     n .le. 100) then
         dxmult = 0.9d+0
      else if (n .le. 250) then
         dxmult = 0.99d+0
      else 
         dxmult = 0.999d+0
      end if

      dxj  = one / n
      do 110, j = 1, n
         dx(j)  =   dxj
         dxj    = - dxj*dxmult
  110 continue

*     ------------------------------------------------------------------
*     Do not perturb x(j) if the  j-th  element is missing.
*     Compute the directional derivative.
*     ------------------------------------------------------------------
      ncheck = 0
      do 120,   j = 1, n
         if (grad(j) .eq. rdummy) then
            dx(j) = zero
         else
            ncheck = ncheck + 1

            xj     =   x(j)                             
            stepbl = - one       
            stepbu =   one
            if (bl(j) .gt. biglow)
     $         stepbl = max( stepbl, bl(j) - xj )
            if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $         stepbu = min( stepbu, bu(j) - xj )

            if (half*(stepbl + stepbu) .lt. zero) then
               dx(j) = dx(j)*stepbl
            else
               dx(j) = dx(j)*stepbu
            end if
         end if
  120 continue

      if (ncheck .eq. 0) then
         if (msglvl .gt. 0  .and.  iPrint .gt. 0) write(iPrint, 3500)
         return
      end if
      gdx    = ddot  ( n, gradu, 1, dx, 1 )

*     ------------------------------------------------------------------
*     Make forward-difference approximation along  p.
*     ------------------------------------------------------------------
      call dcopy ( n,     x, 1, y, 1 )
      call daxpy ( n, h, dx, 1, y, 1 )
      
      mode   = 0
      call funobj( mode, n, y, objf1, gradu, nstate )
      if (mode .lt. 0) go to 9999

      gdiff =    (objf1 - objf) / h
      error = abs(gdiff - gdx ) / (abs(gdx) + one)

      okay  = error .le. oktol

      if (msglvl .gt. 0) then
         if ( okay ) then
            if (iPrint .gt. 0) write(iPrint, 1100)
         else
            if (iPrint .gt. 0) write(iPrint, 1200)
            if (iSumm  .gt. 0) write(iSumm , 1200)
         end if
         if (iPrint .gt. 0) write(iPrint, 1300) gdx, gdiff
      end if

      if (error .ge. point9) inform = 1

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      if (lvrfyc .eq. 1  .or.  lvrfyc .eq. 3) then
         prtHdr = .true.
         itmax  = 3
         ncheck = 0
         nwrong = 0
         ngood  = 0
         jmax   = 0
         emax   = zero
         j1     = jverfy(1)
         j2     = jverfy(2)

*        ---------------------------------------------------------------
*        Loop over each of the components of  x.
*        ---------------------------------------------------------------
         do 500, j = j1, j2

            if (grad(j) .ne. rdummy) then
*              ---------------------------------------------------------
*              Check this gradient component.
*              ---------------------------------------------------------
               ncheck = ncheck + 1
               gj     = grad(j)
               gsize  = abs( gj )
               xj     = x(j)

*              ---------------------------------------------------------
*              Find a finite-difference interval by iteration.
*              ---------------------------------------------------------
               iter   = 0
               epsa   = epsrf*(one + abs( objf ))
               cdest  = zero
               sdest  = zero
               first  = .true.

               stepbl = biglow
               stepbu = bigupp
               if (bl(j) .gt. biglow) stepbl = bl(j) - xj
               if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

               hopt   = two*(one + abs( xj ))*sqrt( epsrf )
               h      = ten*hopt
               if (half*(stepbl + stepbu) .lt. zero) h =  - h

*+             repeat
  400             x(j)  = xj + h
                  call funobj( mode, n, x, f1, gradu, nstate )
                  if (mode .lt. 0) go to 9999

                  x(j)  = xj + h + h
                  call funobj( mode, n, x, f2, gradu, nstate )
                  if (mode .lt. 0) go to 9999
                              
                  call chcore( debug, done, first, epsa, epsrf, objf,
     $                         info, iter, itmax,
     $                         cdest, fdest, sdest, errbnd, f1,
     $                         f2, h, hopt, hphi )

*+             until     done
               if (.not. done) go to 400

*              ---------------------------------------------------------
*              Exit for this variable.
*              ---------------------------------------------------------
               gdiff = cdest
               x(j)  = xj

               error = abs(gdiff - gj) / (gsize + one)
               if (error .ge. emax) then
                  emax  = error
                  jmax  = j
               end if

               okay =  error .le. oktol
               if (okay) then
                  key    = lgood
                  ngood  = ngood  + 1
               else
                  key    = lbad
                  nwrong = nwrong + 1
               end if

               if (msglvl .gt. 0  .and.  iPrint .gt. 0) then

*                 Zero components are not printed.

                  const = info    .eq. 1       .and. 
     $                    abs(gj) .lt. epspt8  .and.  okay
                  if (.not. const) then
                     if (prtHdr) then
                        write(iPrint, 3000)
                        prtHdr = .false.
                     end if

                     if ( okay ) then
                        write(iPrint, 3100) j, xj, hopt, gj, gdiff,
     $                                       key, iter
                     else
                        write(iPrint, 3110) j, xj, hopt, gj, gdiff,
     $                                       key, iter, result(info)
                     end if
                  end if
               end if
            end if
  500    continue

*        ===============================================================
*        Done.
*        ===============================================================
         inform = 0
         if (msglvl .gt. 0) then
            if (nwrong .eq. 0) then
               if (iPrint .gt. 0) write(iPrint, 3200) ngood , ncheck,
     $                                                j1    , j2
            else
               if (iPrint .gt. 0) write(iPrint, 3300) nwrong, ncheck,
     $                                                j1    , j2
               if (iSumm  .gt. 0) write(iSumm , 3300) nwrong, ncheck,
     $                                                j1    , j2
            end if
            write(iPrint, 3400) emax, jmax
         end if
         if (error .ge. point9) inform = 1
      end if

      call dcopy ( n, grad, 1, gradu, 1 )

      return

 9999 inform = mode
      return

 1000 format(/// ' Verification of the objective gradients.'
     $       /   ' ----------------------------------------' )
 1100 format(/   ' The objective gradients seem to be ok.')
 1200 format(/   ' XXX  The objective gradients seem to be incorrect.')
 1300 format(/   ' Directional derivative of the objective', 1p, e18.8/
     $           ' Difference approximation               ', 1p, e18.8 )
 3000 format(// 4x, 'j', 4x, 'x(j)', 5x, 'dx(j)', 11x,
     $           'g(j)', 9x, '  Difference approxn  Itns' /)
 3100 format(  i5, 1p, 2e10.2,      2e18.8, 2x, a4, i6          )
 3110 format(  i5, 1p, 2e10.2,      2e18.8, 2x, a4, i6, 2x, a18 )
 3200 format(/ i7, '  Objective gradients out of the', i6,
     $             '  set in cols', i6, '  through', i6,
     $             '  seem to be ok.')
 3300 format(/   ' XXX  There seem to be', i6,
     $           '  incorrect objective gradients out of the', i6,
     $           '  set in cols', i6, '  through', i6 )
 3400 format(/   ' The largest relative error was', 1p, e12.2,
     $           '   in element', i6 /)
 3500 format(/   ' No gradient elements assigned.' )

*     end of chfgrd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chcJac( inform, lvlder, msglvl,
     $                   ncset, n, ncnln, ldcJ, ldcJu,
     $                   bigbnd, epsrf, oktol, fdchk, xnorm,
     $                   funcon, needc,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   dx, err, x, y )

      implicit           double precision(a-h,o-z)
      integer            needc(*)
      double precision   bl(n), bu(n), c(*), c1(*), cJdx(*),
     $                   cJac(ldcJ,*), cJacu(ldcJu,*), err(*)
      double precision   dx(n), x(n), y(n)
      external           funcon

*     ==================================================================
*     chcJac  checks if the gradients of the constraints have been coded
*     correctly.
*
*     On input,  the values of the constraints at the point x are stored
*     in c.  Their corresponding gradients are stored in cJacu.  If any
*     Jacobian component has not been specified,  it will have a dummy
*     value.  Missing values are not checked.
*
*     A cheap test is first undertaken by calculating the directional
*     derivative using two different methods.  If this proves 
*     satisfactory and no further information is desired, chcJac is 
*     terminated.  Otherwise, chcore is called to give optimal stepsizes
*     and a central-difference approximation to each component of the 
*     Jacobian for which a test is deemed necessary, either by the 
*     program or the user.
*
*     lvrfyc has the following meaning...
*
*       -1        do not perform any check.
*        0        do the cheap test only.
*        2 or 3   do both cheap and full test.
*
*     Systems Optimization Laboratory, Stanford University.
*     Department of Mathematics, University of California, San Diego.
*     Original version written  19-May-1985.
*     This version of  chcJac  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol5np/ lvrfyc, jverfy(4)

      logical            const , debug , done  , first , prtHdr
      logical            needed, okay
      character*4        key   , lbad  , lgood
      character*18       result(0:4)
      parameter         (rdummy =-11111.0d+0              )
      parameter         (zero   =0.0d+0, half   =0.5d+0, point9 =0.9d+0)
      parameter         (one    =1.0d+0, two    =2.0d+0, ten    =1.0d+1)
      parameter         (lbad   ='BAD?', lgood  ='  OK')
      data               result
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'small derivative?'                   /

      inform = 0
      needed = ncnln  .gt. 0  .and.
     $         lvrfyc .eq. 0  .or.   lvrfyc .eq. 2  .or.  lvrfyc .eq. 3
      if (.not. needed) return

      if (msglvl .gt. 0  .and.  iPrint .gt. 0) write(iPrint, 1000)
      debug  = .false.
      nstate = 0

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      h = (one + xnorm)*fdchk

      if (     n .le. 100) then
         dxmult = 0.9d+0
      else if (n .le. 250) then
         dxmult = 0.99d+0
      else 
         dxmult = 0.999d+0
      end if

      dxj  = one / n
      do 110, j = 1, n
         dx(j) =   dxj
         dxj   = - dxj*dxmult
  110 continue

*     ------------------------------------------------------------------
*     Do not perturb  x(j)  if the  j-th  column contains any
*     unknown elements.  Compute the directional derivative for each
*     constraint gradient.
*     ------------------------------------------------------------------
      ncheck = 0
      do 140, j = 1, n
         do 130, i = 1, ncnln
            if (cJac(i,j) .eq. rdummy) then
               dx(j) = zero      
               go to 140
            end if           
  130    continue
         ncheck = ncheck + 1

         xj     =   x(j)
         stepbl = - one
         stepbu =   one
         if (bl(j) .gt. biglow)
     $      stepbl = max( stepbl, bl(j) - xj )
         if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $      stepbu = min( stepbu, bu(j) - xj )

         if (half*(stepbl + stepbu) .lt. zero) then
            dx(j) = dx(j)*stepbl
         else
            dx(j) = dx(j)*stepbu
         end if
  140 continue

      if (ncheck .eq. 0) then
         if (msglvl .gt. 0  .and.  iPrint .gt. 0) write(iPrint, 2300)
      else
*        ---------------------------------------------------------------
*        Compute  (Jacobian)*dx and check it against a 
*        forward-difference approximation along dx.
*        ---------------------------------------------------------------
         call dgemv ( 'Normal', ncnln, n, one, cJacu, ldcJu,
     $                dx, 1, zero, cJdx, 1 )

         call dcopy ( n,     x, 1, y, 1 )
         call daxpy ( n, h, dx, 1, y, 1 )

         call iload ( ncnln, 1, needc, 1 )

         mode   = 0
         call funcon( mode, ncnln, n, ldcJu,
     $                needc, y, c1, cJacu, nstate )
         if (mode .lt. 0) go to 9999

*        Set  err = (c1 - c)/h  - Jacobian*dx.  This should be small.

         do 170, i = 1, ncnln
            err(i) = (c1(i) - c(i)) / h  -  cJdx(i)
  170    continue            
         imax  = idamax( ncnln, err, 1 )
         emax  = abs(err(imax)) / (abs(cJdx(imax)) + one)

         if (msglvl .gt. 0) then
            if (emax .le. oktol) then
               if (iPrint .gt. 0) write(iPrint, 2000)
            else
               if (iPrint .gt. 0) write(iPrint, 2100)
               if (iSumm  .gt. 0) write(iSumm , 2100)
            end if
            if (iPrint .gt. 0) write(iPrint, 2200) emax, imax
         end if
         if (emax .ge. point9) inform = 1
      end if

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      if (lvrfyc .ge. 2) then
         if (lvlder .eq. 3) then

*           Recompute the Jacobian to find the non-constant elements.

            call f06qhf( 'General', ncnln, n, rdummy, rdummy, 
     $                   cJacu, ldcJu )

            call iload ( ncnln, 1, needc, 1 )
            nstate = 0
            mode   = 2
            call funcon( mode, ncnln, n, ldcJu,
     $                   needc, x, c1, cJacu, nstate )
            if (mode .lt. 0) go to 9999
         end if

         call iload ( ncnln, 0, needc, 1 )

         itmax  =   3
         ncheck =   0
         nwrong =   0
         ngood  =   0
         colmax = - one
         jcol   =   0
         irow   =   0
         mode   =   0
         j3     =   jverfy(3)
         j4     =   jverfy(4)

*        ---------------------------------------------------------------
*        Loop over each column.
*        ---------------------------------------------------------------
         do 600, j = j3, j4
            call dload ( ncnln, zero, err, 1 )
            nColj  = 0
            prtHdr = .true.
            xj     = x(j)

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            do 500, i = 1, ncnln
               epsaci   = epsrf*(one + abs( c(i) ))

               if (cJacu(i,j) .ne. rdummy) then
*                 ------------------------------------------------------
*                 Check this Jacobian element.
*                 ------------------------------------------------------
                  ncheck   = ncheck + 1
                  nColj    = nColj + 1
                  needc(i) = 1

                  cij    = cJac(i,j)
                  cJsize = abs( cij )
*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  iter   = 0
                  hopt   = two*(one + abs( xj ))*sqrt( epsrf )
                  h      = ten*hopt*signh
                  cdest  = zero
                  sdest  = zero
                  first  = .true.

*+                repeat
  400                x(j)  = xj + h
                     call funcon( mode, ncnln, n, ldcJu,
     $                            needc, x, c1, cJacu, nstate )
                     if (mode .lt. 0) go to 9999
                     f1    = c1(i)

                     x(j)  = xj + h + h
                     call funcon( mode, ncnln, n, ldcJu,
     $                            needc, x, c1, cJacu, nstate )
                     if (mode .lt. 0) go to 9999
                     f2    = c1(i)
                                                     
                     call chcore( debug,done,first,epsaci,epsrf,c(i),
     $                            info, iter, itmax,
     $                            cdest, fdest, sdest, errbnd, f1,
     $                            f2, h, hopt, hphi )
*+                until     done
                  if (.not. done) go to 400

*                 ------------------------------------------------------
*                 Exit for this element.
*                 ------------------------------------------------------
                  cJdiff   = cdest
                  err(i)   = abs(cJdiff - cij) / (cJsize + one)

                  okay     = err(i) .le. oktol
                  if ( okay ) then
                     key    = lgood
                     ngood  = ngood  + 1
                  else
                     key    = lbad
                     nwrong = nwrong + 1
                  end if

                  if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
                     const = info       .eq. 1       .and.
     $                       abs( cij ) .lt. epspt8  .and.  okay

                     if (.not. const) then
                        if (prtHdr) then
                           write(iPrint, 4000)
                           if ( okay ) then
                              write(iPrint, 4100)
     $                           j  , xj    , hopt, i,
     $                           cij, cJdiff, key , iter
                           else
                              write(iPrint, 4110)
     $                           j  , xj    , hopt, i,
     $                           cij, cJdiff, key , iter, result(info)
                           end if
                           prtHdr = .false.
                        else
                           if ( okay ) then
                              write(iPrint, 4200) 
     $                                        hopt, i,
     $                           cij, cJdiff, key , iter
                           else 
                              write(iPrint, 4210)
     $                                        hopt, i, 
     $                           cij , cJdiff, key , iter, result(info)
                           end if
                        end if
                     end if
                  end if
                  needc(i) = 0
               end if
  500       continue

*           ------------------------------------------------------------
*           Finished with this column.
*           ------------------------------------------------------------
            if (nColj .gt. 0) then
               imax = idamax( ncnln, err, 1 )
               emax = abs( err(imax) )

               if (emax .ge. colmax) then
                  irow   = imax
                  jcol   = j
                  colmax = emax
               end if
            end if
            x(j) = xj
  600    continue

         inform = 0
         if (colmax .ge. point9) inform = 1

         if (msglvl .gt. 0) then
            if (ncheck .eq. 0) then
               if (iPrint .gt. 0) write(iPrint, 4600) ncset
            else
               if (nwrong .eq. 0) then
                  if (iPrint .gt. 0) write(iPrint, 4300) ngood , ncheck,
     $                                                   j3, j4
               else
                  if (iPrint .gt. 0) write(iPrint, 4400) nwrong, ncheck,
     $                                                   j3, j4
                  if (iSumm  .gt. 0) write(iSumm , 4400) nwrong, ncheck,
     $                                                   j3, j4
                end if
               if (iPrint .gt. 0) write(iPrint, 4500) colmax, irow, jcol
            end if
         end if
      end if

*     Copy  ( constants + gradients + dummy values )  back into cJacu.
      
      call f06qff( 'General', ncnln, n, cJac, ldcJ, cJacu, ldcJu )

      return            

 9999 inform = mode
      return

 1000 format(/// ' Verification of the constraint gradients.'
     $       /   ' -----------------------------------------' )
 2000 format(/   ' The constraint Jacobian seems to be ok.')
 2100 format(/   ' XXX  The constraint Jacobian seems to be incorrect.')
 2200 format(/   ' The largest relative error was', 1p, e12.2,
     $           '  in row', i5 /)
 2300 format(/   ' Every column contains a constant or',
     $           ' missing element.')
 4000 format(// ' Column    x(j)     dx(j)    Row   ',
     $          ' Jacobian Value      Difference Approxn  Itns'    )
 4100 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6         )
 4110 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4200 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6         )
 4210 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4300 format(/ i7, '  constraint Jacobian elements out of the', i6,
     $             '  set in cols', i6, '  through', i6,
     $             '  seem to be ok.')
 4400 format(/   ' XXX  There seem to be', i6,
     $           '  incorrect Jacobian elements out of the', i6,
     $           '  set in cols', i6, '  through', i6 )
 4500 format(/ ' The largest relative error was', 1p, e12.2,
     $         '  in row', i5, ',  column', i5 /)
 4600 format(  ' All', i6, '   assigned Jacobian elements are',
     $         ' constant.' )

*     end of  chcJac
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine chfJac( inform, lvlder, msglvl,
     $                   nfset, m, n, ldfJ, ldfJu,
     $                   bigbnd, epsrf, oktol, fdchk, xnorm,
     $                   funobj,
     $                   bl, bu, f, f1, fJac, fJacu, fJdx,
     $                   dx, err, x, y )

      implicit           double precision(a-h,o-z)
      double precision   bl(n), bu(n), dx(n), x(n), y(n)
      double precision   f(m), f1(m), fJdx(m),
     $                   fJac(ldfJ,*), fJacu(ldfJu,*), err(m)
      external           funobj

*     ==================================================================
*     chfJac  checks if the objective Jacobian matrix has been coded
*     correctly.
*
*     On input,  the values of the objective vector at the point x are
*     stored in f.  Their corresponding gradients are stored in fJacu.
*     If any Jacobian component has not been specified,  it will have a
*     dummy value.  Missing values are not checked.
*
*     A cheap test is first undertaken by calculating the directional
*     derivative using two different methods.  If this proves
*     satisfactory and no further information is desired, chfJac is
*     terminated.  Otherwise, chcore is called to give optimal stepsizes
*     and a central-difference approximation to each component of the
*     Jacobian for which a test is deemed necessary, either by the
*     program or the user.
*
*     lvrfyc has the following meaning...
*
*       -1        do not perform any check.
*        0        do the cheap test only.
*        2 or 3   do both cheap and full test.
*
*     Systems Optimization Laboratory, Stanford University.
*     Department of Mathematics, University of California, San Diego.
*     Original version written  10-May-1988.
*     This version of  chfJac  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5np/ lvrfyc, jverfy(4)

      logical            const , debug , done  , first , prtHdr
      logical            needed, okay
      character*4        key   , lbad  , lgood
      character*18       result(0:4)
      parameter         (rdummy =-11111.0d+0              )
      parameter         (zero   =0.0d+0, half   =0.5d+0, point9 =0.9d+0)
      parameter         (one    =1.0d+0, two    =2.0d+0, ten    =1.0d+1)
      parameter         (lbad   ='BAD?', lgood  ='  OK')
      data               result
     $                 / '                 ', 'Constant?      ',
     $                   'Linear or odd?   ', 'Too nonlinear?',
     $                   'Small derivative?'                   /

      inform = 0
      needed = lvrfyc .eq. 0  .or.  lvrfyc .eq. 1  .or.  lvrfyc .eq. 3
      if (.not. needed) return

      if (iPrint .gt. 0) write(iPrint, 1000)
      debug  = .false. 
      nstate = 0

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Perform the cheap test.
*     ==================================================================
      h = (one + xnorm)*fdchk

      if (     n .le. 100) then
         dxmult = 0.9d+0
      else if (n .le. 250) then
         dxmult = 0.99d+0
      else 
         dxmult = 0.999d+0
      end if

      dxj  = one / n
      do 110, j = 1, n
         dx(j) =   dxj
         dxj   = - dxj*dxmult
  110 continue

*     ------------------------------------------------------------------
*     Do not perturb  x(j)  if the  j-th  column contains any
*     unknown elements.  Compute the directional derivative for each
*     objective gradient.
*     ------------------------------------------------------------------
      ncheck = 0
      do 140, j = 1, n
         do 130, i = 1, m
            if (fJac(i,j) .eq. rdummy) then
               dx(j) = zero
               go to 140
            end if
  130    continue
         ncheck = ncheck + 1

         xj     =   x(j)
         stepbl = - one
         stepbu =   one
         if (bl(j) .gt. biglow)
     $      stepbl = max( stepbl, bl(j) - xj )
         if (bu(j) .lt. bigupp  .and.  bu(j) .gt. bl(j))
     $      stepbu = min( stepbu, bu(j) - xj )

         if (half*(stepbl + stepbu) .lt. zero) then
            dx(j) = dx(j)*stepbl
         else
            dx(j) = dx(j)*stepbu
         end if
  140 continue

      if (ncheck .eq. 0) then
         if (iPrint .gt. 0) write(iPrint, 2300)
      else
*        ---------------------------------------------------------------
*        Compute  (Jacobian)*dx and check it against a 
*        forward-difference approximation along dx.
*        ---------------------------------------------------------------
         call dgemv ( 'Normal', m, n, one, fJacu, ldfJu,
     $                 dx, 1, zero, fJdx, 1 )

         call dcopy ( n,     x, 1, y, 1 )
         call daxpy ( n, h, dx, 1, y, 1 )

         mode   = 0
         call funobj( mode, m, n, ldfJu, y, f1, fJacu, nstate )
         if (mode .lt. 0) go to 9999

*        Set  err = (f1 - f)/h  - Jacobian*dx.  This should be small.

         do 170, i = 1, m
            err(i) = (f1(i) - f(i)) / h  -  fJdx(i)
  170    continue
         imax  = idamax( m, err, 1 )
         emax  = abs(err(imax)) / (abs(fJdx(imax)) + one)

         if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
            if (emax .le. oktol) then
               if (iPrint .gt. 0) write(iPrint, 2000)
            else
               if (iPrint .gt. 0) write(iPrint, 2100)
               if (iSumm  .gt. 0) write(iSumm , 2100)
            end if
            if (iPrint .gt. 0) write(iPrint, 2200) emax, imax
         end if
         if (emax .ge. point9) inform = 1
      end if

*     ==================================================================
*     Component-wise check.
*     ==================================================================
      if (lvrfyc .ge. 2) then
         if (lvlder .eq. 3) then

*           Recompute the Jacobian to find the non-constant elements.

            call f06qhf( 'General', m, n, rdummy, rdummy,
     $                   fJacu, ldfJu )

            nstate = 0
            mode   = 2
            call funobj( mode, m, n, ldfJu, x, f1, fJacu, nstate )
            if (mode .lt. 0) go to 9999
         end if

         itmax  =   3
         ncheck =   0
         nwrong =   0
         ngood  =   0
         colmax = - one
         jcol   =   0
         irow   =   0
         mode   =   0
         j1     =   jverfy(1)
         j2     =   jverfy(2)

*        ---------------------------------------------------------------
*        Loop over each column.
*        ---------------------------------------------------------------
         do 600, j = j1, j2
            call dload ( m, zero, err, 1 )
            prtHdr = .true.
            nColj  = 0 
            xj     = x(j)

            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj   
                                                 
            signh  = one
            if (half*(stepbl + stepbu) .lt. zero) signh =  - one

            do 500, i = 1, m
               epsafi   = epsrf*(one + abs( f(i) ))

               if (fJacu(i,j) .ne. rdummy) then
*                 ------------------------------------------------------
*                 Check this Jacobian element.
*                 ------------------------------------------------------
                  ncheck = ncheck + 1
                  nColj  = nColj  + 1
                  fij    = fJac(i,j)
                  fjsize = abs( fij )
*                 ------------------------------------------------------
*                 Find a finite-difference interval by iteration.
*                 ------------------------------------------------------
                  iter   = 0
                  hopt   = two*(one + abs( xj ))*sqrt( epsrf )
                  h      = ten*hopt*signh
                  cdest  = zero
                  sdest  = zero
                  first  = .true.

*+                repeat
  400                x(j)  = xj + h
                     call funobj( mode, m, n, ldfJu,
     $                            x, f1, fJacu, nstate )
                     if (mode .lt. 0) go to 9999
                     fforw = f1(i)

                     x(j)  = xj + h + h
                     call funobj( mode, m, n, ldfJu,
     $                            x, f1, fJacu, nstate )
                     if (mode .lt. 0) go to 9999
                     fback = f1(i)

                     call chcore( debug,done,first, epsafi, epsrf, f(i),
     $                            info, iter, itmax,
     $                            cdest, fdest, sdest, errbnd, fforw,
     $                            fback, h, hopt, hphi )
*+                until     done
                  if (.not. done) go to 400

*                 ------------------------------------------------------
*                 Exit for this element.
*                 ------------------------------------------------------
                  fjdiff   = cdest
                  err(i)   = abs(fjdiff - fij) / (fjsize + one)

                  okay     = err(i) .le. oktol
                  if ( okay ) then
                     key    = lgood
                     ngood  = ngood  + 1
                  else
                     key    = lbad
                     nwrong = nwrong + 1
                  end if

                  if (msglvl .gt. 0  .and.  iPrint .gt. 0) then
                     const = info       .eq. 1       .and. 
     $                       abs( fij ) .lt. epspt8  .and.  okay

                     if (.not. const) then
                        if (prtHdr) then
                           write(iPrint, 4000)
                           if ( okay ) then
                              write(iPrint, 4100)
     $                             j, xj    , hopt, i,
     $                           fij, fjdiff, key , iter
                           else
                              write(iPrint, 4110)
     $                             j, xj    , hopt, i,
     $                           fij, fjdiff, key , iter, result(info)
                           end if
                           prtHdr = .false.
                        else
                           if ( okay ) then
                              write(iPrint, 4200)
     $                                        hopt, i,
     $                           fij, fjdiff, key , iter
                           else
                              write(iPrint, 4210)
     $                                        hopt, i,
     $                           fij, fjdiff, key , iter, result(info)
                           end if
                        end if
                     end if
                  end if
               end if
  500       continue

*           ------------------------------------------------------------
*           Finished with this column.
*           ------------------------------------------------------------
            if (nColj .gt. 0) then
               imax = idamax( m, err, 1 )
               emax = abs( err(imax) )

               if (emax .ge. colmax) then
                  irow   = imax
                  jcol   = j
                  colmax = emax
               end if
            end if
            x(j) = xj
  600    continue

         inform = 0
         if (colmax .ge. point9) inform = 1

         if (msglvl .gt. 0) then
            if (ncheck .eq. 0) then
               if (iPrint .gt. 0) write(iPrint, 4600) nfset
            else
               if (nwrong .eq. 0) then
                  if (iPrint .gt. 0) write(iPrint, 4300) ngood , ncheck,
     $                                                   j1, j2
               else
                  if (iPrint .gt. 0) write(iPrint, 4400) nwrong, ncheck,
     $                                                   j1, j2
                  if (iSumm  .gt. 0) write(iSumm , 4400) nwrong, ncheck,
     $                                                   j1, j2
               end if
               if (iPrint .gt. 0) write(iPrint, 4500) colmax, irow, jcol
            end if
         end if
      end if

*     Copy  ( constants + gradients + dummy values )  back into fJacu.

      call f06qff( 'General', m, n, fJac, ldfJ, fJacu, ldfJu )

      return

 9999 inform = mode
      return

 1000 format(/// ' Verification of the objective gradients.'
     $       /   ' ----------------------------------------' )
 2000 format(/   ' The objective Jacobian seems to be ok.')
 2100 format(/   ' XXX  The objective Jacobian seems to be incorrect.')
 2200 format(/   ' The largest relative error was', 1p, e12.2,
     $           '  in row', i5 /)
 2300 format(/   ' Every column contains a constant or',
     $           ' missing element.')
 4000 format(// ' Column    x(j)     dx(j)    Row   ',
     $          ' Jacobian Value      Difference Approxn  Itns'    )
 4100 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6         )
 4110 format(/ i7,      1p, 2e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4200 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6         )
 4210 format(  7x, 10x, 1p,  e10.2, i5, 2e18.8, 2x, a4, i6, 2x, a18)
 4300 format(/ i7, '  Objective Jacobian elements out of the', i6,
     $             '  set in cols', i6, '  through', i6,
     $             '  seem to be ok.')
 4400 format(/   ' XXX  There seem to be', i6,
     $           '  incorrect objective Jacobian elements out of the',
     $             i6, '  set in cols', i6, '  through', i6 )
 4500 format(/ ' The largest relative error was', 1p, e12.2,
     $         '  in row', i5, ',  column', i5 /)
 4600 format(  ' All', i6, '   assigned Jacobian elements are',
     $         ' constant.' )

*     end of chfJac
      end

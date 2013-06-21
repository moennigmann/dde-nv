*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  rtcomsubs.f
*
*     cmalf1   cmalf    cmchk    cmperm   cmtsol   cmwrp 
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmalf1( firstv, negstp, bigalf, bigbnd, pnorm,
     $                   jadd1 , jadd2 , palfa1, palfa2,
     $                   istate, n, nctotl,
     $                   Anorm, Ap, Ax, bl, bu, featol, p, x )

      implicit           double precision(a-h,o-z)
      logical            firstv, negstp
      integer            istate(nctotl)
      double precision   Anorm(*), Ap(*), Ax(*)
      double precision   bl(nctotl), bu(nctotl), featol(nctotl),
     $                   p(n), x(n)

*     ==================================================================
*     cmalf1  finds steps palfa1, palfa2 such that
*        x + palfa1*p  reaches a linear constraint that is currently not
*                      in the working set but is satisfied.
*        x + palfa2*p  reaches a linear constraint that is currently not
*                      in the working set but is violated.
*     The constraints are perturbed by an amount featol, so that palfa1
*     is slightly larger than it should be,  and palfa2 is slightly
*     smaller than it should be.  This gives some leeway later when the
*     exact steps are computed by cmalf.
*
*     Constraints in the working set are ignored  (istate(j) .ge. 1).
*
*     If negstp is true, the search direction will be taken to be  - p.
*
*
*     Values of istate(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     The values  -2  and  -1  do not occur once a feasible point has
*     been found.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written  May 1980.
*     This version of cmalf1 dated 27-Oct-92.
*     ==================================================================
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            lastv
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      lastv  = .not. firstv
      jadd1  = 0
      jadd2  = 0
      palfa1 = bigalf

      palfa2 = zero
      if (firstv) palfa2 = bigalf

      do 200 j = 1, nctotl
         js = istate(j)
         if (js .le. 0) then
            if (j .le. n) then
               atx    = x(j)
               atp    = p(j)
               rownrm = one
            else
               i      = j - n
               atx    = Ax(i)
               atp    = Ap(i)
               rownrm = one  +  Anorm(i)
            end if
            if (negstp) atp = - atp

            if ( abs( atp ) .le. epspt9*rownrm*pnorm) then

*              This constraint appears to be constant along P.  It is
*              not used to compute the step.  Give the residual a value
*              that can be spotted in the debug output.

               res = - one
            else if (atp .le. zero  .and.  js .ne. -2) then
*              ---------------------------------------------------------
*              a'x  is decreasing and the lower bound is not violated.
*              ---------------------------------------------------------
*              First test for smaller PALFA1.

               absatp = - atp
               if (bl(j) .gt. (-bigbnd)) then
                  res    = atx - bl(j) + featol(j)
                  if (bigalf*absatp .gt. abs( res )) then
                     if (palfa1*absatp .gt. res)  then
                        palfa1 = res / absatp
                        jadd1  = j
                     end if
                  end if
               end if

               if (js .eq. -1) then

*                 The upper bound is violated.  Test for either larger
*                 or smaller PALFA2, depending on the value of FIRSTV.

                  res    = atx - bu(j) - featol(j)
                  if (bigalf*absatp .gt. abs( res )) then
                     if (firstv  .and.  palfa2*absatp .gt. res  .or.
     $                    lastv  .and.  palfa2*absatp .lt. res) then
                        palfa2 = res / absatp
                        jadd2  = j
                     end if
                  end if
               end if
            else if (atp .gt. zero  .and.  js .ne. -1) then
*              ---------------------------------------------------------
*              a'x  is increasing and the upper bound is not violated.
*              ---------------------------------------------------------
*              Test for smaller PALFA1.

               if (bu(j) .lt. bigbnd) then
                  res = bu(j) - atx + featol(j)
                  if (bigalf*atp .gt. abs( res )) then
                     if (palfa1*atp .gt. res) then
                        palfa1 = res / atp
                        jadd1  = j
                     end if
                  end if
               end if

               if (js .eq. -2) then

*                 The lower bound is violated.  Test for a new PALFA2.

                  res  = bl(j) - atx - featol(j)
                  if (bigalf*atp .gt. abs( res )) then
                     if (firstv  .and.  palfa2*atp .gt. res  .or.
     $                    lastv  .and.  palfa2*atp .lt. res) then
                        palfa2 = res / atp
                        jadd2  = j
                     end if
                  end if
               end if
            end if
         end if
  200 continue

*     end of cmalf1
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmalf ( firstv, hitlow, 
     $                   istate, inform, jadd, n, nctotl, numinf,
     $                   alfa, palfa, atphit, 
     $                   bigalf, bigbnd, pnorm,
     $                   Anorm, Ap, Ax, bl, bu, 
     $                   featol, p, x )

      implicit           double precision(a-h,o-z)
      integer            istate(nctotl)
      double precision   Anorm(*), Ap(*), Ax(*),
     $                   bl(nctotl), bu(nctotl), featol(nctotl),
     $                   p(n), x(n)
      logical            firstv, hitlow

*     ==================================================================
*     cmalf   finds a step alfa such that the point x + alfa*p reaches
*     one of the linear constraints (including bounds).  Two possible
*     steps are defined as follows...
*
*     alfa1   is the maximum step that can be taken without violating
*             one of the linear constraints that is currently satisfied.
*     alfa2   reaches a linear constraint that is currently violated.
*             Usually this will be the furthest such constraint along p,
*             but if firstv = .true. it will be the first one along p.
*             This is used only when the problem has been determined to 
*             be infeasible, and the sum of infeasibilities are being
*             minimized.  (alfa2  is not defined if numinf = 0.)
*
*     alfa will usually be the minimum of alfa1 and alfa2.
*     alfa could be negative (since we allow inactive constraints
*     to be violated by as much as featol).  In such cases, a
*     third possible step is computed, to find the nearest satisfied
*     constraint (perturbed by featol) along the direction  - p.
*     alfa  will be reset to this step if it is shorter.  This is the
*     only case for which the final step  alfa  does not move x exactly
*     onto a constraint (the one denoted by jadd).
*
*     Constraints in the working set are ignored  (istate(j) ge 1).
*
*     jadd    denotes which linear constraint is reached.
*
*     hitlow  indicates whether it is the lower or upper bound that
*             has restricted alfa.
*
*     Values of istate(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     The values -2 and -1 do not occur once a feasible point has been
*     found.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written  May 1980.
*     This version of  cmalf  dated  28-Oct-92.
*     ==================================================================
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            hlow1, hlow2, lastv, negstp, step2
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      inform = 0

*     ------------------------------------------------------------------
*     First pass -- find steps to perturbed constraints, so that
*     palfa1 will be slightly larger than the true step, and
*     palfa2 will be slightly smaller than it should be.
*     In degenerate cases, this strategy gives us some freedom in the
*     second pass.  The general idea follows that described by P.M.J.
*     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
*     ------------------------------------------------------------------

      negstp = .false.
      call cmalf1( firstv, negstp, bigalf, bigbnd, pnorm,
     $             jadd1, jadd2, palfa1, palfa2,
     $             istate, n, nctotl,
     $             Anorm, Ap, Ax, bl, bu, featol, p, x )

      jsave1 = jadd1
      jsave2 = jadd2

*     ------------------------------------------------------------------
*     Second pass -- recompute step-lengths without perturbation.
*     Amongst constraints that are less than the perturbed steps,
*     choose the one (of each type) that makes the largest angle
*     with the search direction.
*     ------------------------------------------------------------------
      alfa1  = bigalf
      alfa2  = zero
      if (firstv) alfa2 = bigalf

      apmax1 = zero
      apmax2 = zero
      atp1   = zero
      atp2   = zero
      hlow1  = .false.
      hlow2  = .false.
      lastv  = .not. firstv

      do 400 j = 1, nctotl
         js = istate(j)
         if (js .le. 0) then
            if (j  .le. n)  then
               atx    = x(j)
               atp    = p(j)
               rownrm = one
            else
               i      = j - n
               atx    = Ax(i)
               atp    = Ap(i)
               rownrm = Anorm(i) + one
            end if

            if ( abs( atp ) .le. epspt9*rownrm*pnorm) then

*              This constraint appears to be constant along p.  It is
*              not used to compute the step.  Give the residual a value
*              that can be spotted in the debug output.

               res = - one
            else if (atp .le. zero  .and.  js .ne. -2) then
*              ---------------------------------------------------------
*              a'x  is decreasing.
*              ---------------------------------------------------------
*              The lower bound is satisfied.  Test for smaller alfa1.

               absatp = - atp
               if (bl(j) .gt. (-bigbnd)) then
                  res    = atx - bl(j)
                  if (palfa1*absatp .ge. res  .or.  j .eq. jsave1) then
                     if (apmax1*rownrm*pnorm .lt. absatp) then
                        apmax1 = absatp / (rownrm*pnorm)
                        alfa1  = res / absatp
                        jadd1  = j
                        atp1   = atp
                        hlow1  = .true.
                     end if
                  end if
               end if

               if (js .eq. -1)  then

*                 The upper bound is violated.  Test for either a bigger
*                 or smaller alfa2,  depending on the value of firstv.

                  res    = atx - bu(j)
                  if (     (firstv  .and.  palfa2*absatp .ge. res
     $                 .or.  lastv  .and.  palfa2*absatp .le. res)
     $                 .or.  j .eq.  jsave2) then
                     if (apmax2*rownrm*pnorm .lt. absatp) then
                        apmax2 = absatp / (rownrm*pnorm)
                        if      (absatp .ge. one          ) then
                           alfa2 = res / absatp
                        else if (res    .lt. bigalf*absatp) then
                           alfa2 = res / absatp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2  = j
                        atp2   = atp
                        hlow2  = .false.
                     end if
                  end if
               end if
            else if (atp .gt. zero  .and.  js .ne.  -1)  then
*              ---------------------------------------------------------
*              a'x  is increasing and the upper bound is not violated.
*              ---------------------------------------------------------
*              Test for smaller alfa1.

               if (bu(j) .lt. bigbnd) then
                  res = bu(j) - atx
                  if (palfa1*atp .ge. res  .or.  j .eq. jsave1) then
                     if (apmax1*rownrm*pnorm .lt. atp) then
                        apmax1 = atp / (rownrm*pnorm)
                        alfa1  = res / atp
                        jadd1  = j
                        atp1   = atp
                        hlow1  = .false.
                     end if
                  end if
               end if

               if (js .eq. -2)  then

*                 The lower bound is violated.  Test for a new ALFA2.

                  res    = bl(j) - atx
                  if (     (firstv  .and.  palfa2*atp .ge. res
     $                 .or.  lastv  .and.  palfa2*atp .le. res)
     $                 .or.  j .eq.  jsave2) then
                     if (apmax2*rownrm*pnorm .lt. atp) then
                        apmax2 = atp / (rownrm*pnorm)
                        if      (atp .ge. one       ) then
                           alfa2 = res / atp
                        else if (res .lt. bigalf*atp) then
                           alfa2 = res / atp
                        else
                           alfa2 = bigalf
                        end if
                        jadd2  = j
                        atp2   = atp
                        hlow2  = .true.
                     end if
                  end if
               end if
            end if
         end if
  400 continue

*     ==================================================================
*     Determine alfa, the step to be taken.
*     ==================================================================
*     In the infeasible case, check whether to take the step alfa2
*     rather than alfa1...

      step2 = numinf .gt. 0  .and.  jadd2 .gt. 0

*     We do so if alfa2 is less than alfa1 or (if firstv is false)
*     lies in the range  (alfa1, palfa1)  and has a smaller value of
*     atp.

      step2 = step2 .and. (alfa2 .lt. alfa1   .or.   lastv  .and.
     $                     alfa2 .le. palfa1  .and.  apmax2 .ge. apmax1)

      if (step2) then
         alfa   = alfa2
         palfa  = palfa2
         jadd   = jadd2
         atphit = atp2
         hitlow = hlow2
      else
         alfa   = alfa1
         palfa  = palfa1
         jadd   = jadd1
         atphit = atp1
         hitlow = hlow1

*        If alfa1 is negative, the constraint to be added (jadd)
*        remains unchanged, but alfa may be shortened to the step
*        to the nearest perturbed satisfied constraint along  - p.

         negstp = alfa .lt. zero
         if (negstp) then
            call cmalf1( firstv, negstp, bigalf, bigbnd, pnorm,
     $                   jadd1, jadd2, palfa1, palfa2,
     $                   istate, n, nctotl,
     $                   Anorm, Ap, Ax, bl, bu, featol, p, x )

            alfa = - min( abs( alfa ), palfa1 )
         end if
      end if

*     Test for undefined or infinite step.

      if (jadd .eq. 0) then
         alfa   = bigalf
         palfa  = bigalf
      end if

      if (alfa .ge. bigalf) inform = 3

*     end of cmalf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmchk ( nerror, msglvl, lcrash, userkx,
     $                   liwork, lwork, litotl, lwtotl,
     $                   n, nclin, ncnln,
     $                   istate, kx, named, names,
     $                   bigbnd, bl, bu, clamda, x )

      implicit           double precision(a-h,o-z)
      character*16       names(*)
      logical            named, userkx
      integer            istate(n+nclin+ncnln), kx(n)
      double precision   bl(n+nclin+ncnln), bu(n+nclin+ncnln), x(n)
      double precision   clamda(n+nclin+ncnln)

*     ==================================================================
*     cmchk   checks the data input to various optimizers.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written 10-May-1980.
*     Fortran 77 version written  5-October-1984.
*     This version of cmchk dated  04-Oct-93.       
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 

      logical            ok
      parameter         (zero = 0.0d+0)
      character*5        id(3)
      data                id(1)   ,  id(2)   ,  id(3)
     $                 / 'varbl'  , 'lncon'  , 'nlcon'   /

      nerror = 0

*     ------------------------------------------------------------------
*     Check that there is enough workspace to solve the problem.
*     ------------------------------------------------------------------
      ok     = litotl .le. liwork  .and.  lwtotl .le. lwork
      if (.not. ok)  then
         nerror = nerror + 1
         if (iPrint .gt. 0) then
            write (iPrint, 1100) liwork, lwork, litotl, lwtotl
            write (iPrint, 1110)
         end if
      else if (msglvl .gt. 0)  then
         if (iPrint .gt. 0) then
            write (iPrint, 1100) liwork, lwork, litotl, lwtotl
         end if
      end if

      if (userkx) then
*        ---------------------------------------------------------------
*        Check for a valid KX.
*        ---------------------------------------------------------------
         ifail = 1
         call cmperm( kx, 1, n, ifail )
         if (ifail .ne. 0) then
            if (iPrint .gt. 0) write (iPrint, 1300)
            nerror = nerror + 1
         end if
      end if

*     ------------------------------------------------------------------
*     Check the bounds on all variables and constraints.
*     ------------------------------------------------------------------
      do 200, j = 1, n+nclin+ncnln
         b1     = bl(j)
         b2     = bu(j)
         ok     = b1 .lt. b2  .or. 
     $            b1 .eq. b2  .and.  abs(b1) .lt. bigbnd

         if (.not. ok)  then
            nerror = nerror + 1
            if (j .gt. n+nclin)  then
               k  = j - n - nclin
               l  = 3
            else if (j .gt. n)  then
               k  = j - n
               l  = 2
            else
               k = j
               l = 1
            end if
            if (iPrint .gt. 0) then
               if (named) then
                  if (b1 .eq. b2) then
                     write (iPrint, 1210) names(j), b1, bigbnd
                  else
                     write (iPrint, 1215) names(j), b1, b2
                  end if
               else 
                  if (b1 .eq. b2) then
                     write (iPrint, 1200) id(l), k, b1, bigbnd
                  else
                     write (iPrint, 1205) id(l), k, b1, b2
                  end if
               end if
            end if
         end if
  200 continue

*     ------------------------------------------------------------------
*     If warm start, check  istate and clamda.
*     ------------------------------------------------------------------
      if (lcrash .eq. 1) then         
         do 420, j = 1, n+nclin+ncnln
            is     = istate(j)
            ok     = is .ge. (- 2)   .and.   is .le. 4
            if (.not. ok)  then
               nerror = nerror + 1
               if (iPrint .gt. 0) write (iPrint, 1500) j, is
            end if
  420    continue

         if (nerror .eq. 0) then
            do 430, i = 1, ncnln
               j      = n + nclin + i
               is     = istate(j)
               cmul   = clamda(j)

               if      (is .eq. 0) then
                  cmul = zero

               else if (is .eq. 1) then
                  if (bl(j) .le. -bigbnd) is = 0
                  if (cmul  .lt.  zero  .or.  is .eq. 0) cmul = zero 

               else if (is .eq. 2) then
                  if (bu(j) .ge.  bigbnd) is = 0
                  if (cmul  .gt.  zero  .or.  is .eq. 0) cmul = zero 

               else if (is .eq. 3) then
                  if (bl(j) .lt.   bu(j)) is = 0
               end if

               istate(j) = is
               clamda(j) = cmul
  430       continue
         end if
      end if

      return

 1100 format(/ ' Workspace provided is     iw(', i8,
     $         '),  w(', i8, ').' /
     $         ' To solve problem we need  iw(', i8,
     $         '),  w(', i8, ').')
 1110 format(/ ' XXX  Not enough workspace to solve problem.')
 1200 format(/ ' XXX  The equal bounds on  ', a5, i3,
     $         '  are infinite.   Bounds =', g16.7,
     $         '  bigbnd =', g16.7)
 1205 format(/ ' XXX  The bounds on  ', a5, i3,
     $         '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1210 format(/ ' XXX  The equal bounds on  ', a8,
     $         '  are infinite.   Bounds =', g16.7,
     $         '  bigbnd =', g16.7)
 1215 format(/ ' XXX  The bounds on  ', a8,
     $         '  are inconsistent.   bl =', g16.7, '   bu =', g16.7)
 1300 format(/ ' XXX  kx has not been supplied as a valid',
     $         '  permutation.' )
 1500 format(/ ' XXX  Component', i5, '  of  istate  is out of',
     $         ' range...', i10)

*     end of cmchk
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmperm( kx, m1, m2, ifail )

      integer            ifail, m1, m2
      integer            kx(m2)

*     ==================================================================
*     CMPERM checks that elements M1 to M2 of KX contain a valid
*     permutation of the integers M1 to M2. The contents of KX are
*     unchanged on exit.
*
*     SOL version of NAG Library routine M01ZBF.
*     Written by N.N.Maclaren, University of Cambridge.
*     This version of CMPERM dated 18-June-1986.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 

      integer            i, ierr, j, k

*     Check the parameters.

      if (m2 .lt. 1  .or.  m1 .lt. 1  .or.  m1 .gt. m2) then
         ierr = 1
      else
         ierr = 0

*        Check that KX is within range.

         do 20 i = m1, m2
            j = kx(i)
            if ((j .lt. m1) .or. (j .gt. m2)) go to 100
            if (i .ne. j) kx(i) = -j
   20    continue

*        Check that no value is repeated.

         do 60 i = m1, m2
            k = - kx(i)
            if (k .ge. 0) then
               j     = i
   40          kx(j) = k
               j     = k
               k     = - kx(j)
               if (k .gt. 0) go to 40
               if (j .ne. i) go to 120
            end if
   60    continue
      end if

*     Return

   80 if (ierr .ne. 0) then
         ifail = ierr
      else
         ifail = 0
      end if
      return
  100 ierr = 2
      if (iPrint .gt. 0) write (iPrint, fmt=1200) i, j
      go to 140
  120 ierr = 3
      if (iPrint .gt. 0) write (iPrint, fmt=1300) j

*     Restore KX.

  140 do 160 i = m1, m2
         kx(i) = abs(kx(i))
  160 continue
      go to 80

 1200 format(/ ' XXX  kx(',I6,') contains an out-of-range value =', i16)
 1300 format(/ ' XXX  kx contains a duplicate value =',             i16)

*     end of cmperm
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmtsol( mode, ldt, n, t, y )

      implicit           double precision(a-h,o-z)
      integer            mode, ldt, n
      double precision   t(ldt,*), y(n)

*     ==================================================================
*     cmtsol  solves equations involving a reverse-triangular matrix  T
*     and a right-hand-side vector  y,  returning the solution in  y.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written February-1985.
*     ==================================================================
      parameter        ( zero = 0.0d+0 )

      n1 = n + 1
      if (mode .eq. 1) then

*        Mode = 1  ---  Solve  T * y(new) = y(old).

         do 100, j = 1, n
            jj     = n1 - j
            yj     = y(j)/t(j,jj)
            y(j)   = yj
            l      = jj - 1
            if (l .gt. 0  .and.  yj .ne. zero)
     $         call daxpy( l, (-yj), t(j+1,jj), 1, y(j+1), 1 )
  100    continue
      else

*        Mode = 2  ---  Solve  T' y(new) = y(old).

         do 500, j = 1, n
            jj     = n1 - j
            yj     = y(j)/t(jj,j)
            y(j)   = yj
            l      = jj - 1
            if (l .gt. 0  .and.  yj .ne. zero)
     $         call daxpy( l, (-yj), t(jj,j+1), ldt, y(j+1), 1 )
  500    continue
      end if

*     Reverse the solution vector.

      if (n .gt. 1) then
         l = n/2
         do 800, j = 1, l
            jj     = n1 - j
            yj     = y(j)
            y(j)   = y(jj)
            y(jj)  = yj
  800    continue
      end if

*     end of cmtsol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmwrp ( nfree, ldA,
     $                   n, nclin, nctotl,
     $                   nactiv, istate, kactiv, kx,
     $                   A, bl, bu, c, clamda, featol, r, rlamda, x )

      implicit           double precision(a-h,o-z)
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   A(ldA,*), bl(nctotl), bu(nctotl), c(*),
     $                   clamda(nctotl), featol(nctotl), r(nctotl)
      double precision   rlamda(n), x(n)

*     ==================================================================
*     cmwrp   creates the expanded Lagrange multiplier vector clamda.
*     and resets istate for the printed solution.
*
*     This version of cmwrp  is for reverse-triangular T.
*
*     Original Fortran 77 version written  05-May-93.
*     This version of  cmwrp   dated  05-May-93.
*     ==================================================================
      parameter         (zero  = 0.0d+0)

      nfixed = n     - nfree
      nplin  = n     + nclin
      nZ     = nfree - nactiv

*     Expand multipliers for bounds, linear and nonlinear constraints
*     into the  clamda  array.

      call dload ( nctotl, zero, clamda, 1 )
      do 150, k = 1, nactiv+nfixed
         if (k .le. nactiv) j = kactiv(k) + n
         if (k .gt. nactiv) j = kx(nZ+k)
         clamda(j) = rlamda(k)
  150 continue

*     Reset istate if necessary.

      do 500, j = 1, nctotl
         b1     = bl(j)
         b2     = bu(j)

         if (j .le. n) then
            rj  = x(j)
         else if (j .le. nplin) then
            i   = j - n
            rj  = ddot  ( n, A(i,1), ldA, x, 1 )
         else
            i   = j - nplin
            rj  = c(i)
         end if

         is     = istate(j)
         slk1   = rj - b1
         slk2   = b2 - rj
         tol    = featol(j)
         if (                  slk1 .lt. -tol) is = - 2
         if (                  slk2 .lt. -tol) is = - 1
         if (is .eq. 1  .and.  slk1 .gt.  tol) is =   0
         if (is .eq. 2  .and.  slk2 .gt.  tol) is =   0
         istate(j) = is
         r(j)      = rj
  500 continue

*     end of cmwrp 
      end

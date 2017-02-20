*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     File  cmsubs.f
*
*     cmalf1   cmalf    cmchk    cmperm   cmprt   +cmqmul  +cmr1md
*    +cmrswp   cmtsol
*   (+) Also in cmoptsubs.f
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
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      logical            lastv
      intrinsic          abs
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      if (cmdbg  .and.  icmdbg(3) .gt. 0) write (iPrint, 1100)
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

            if (cmdbg  .and.  icmdbg(3) .gt. 0)
     $         write (iPrint, 1200) j, js, featol(j), res,
     $                            atp, jadd1, palfa1, jadd2, palfa2
         end if
  200 continue

      return

 1100 format(/ '    j  js         featol        res             Ap',
     $         '     jadd1       palfa1     jadd2       palfa2' /)
 1200 format(i5, i4, 3g15.5, 2(i6, g17.7))

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
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      logical            hlow1, hlow2, lastv, negstp, step2
      intrinsic          abs, min
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
      if (cmdbg  .and.  icmdbg(3) .gt. 0) write (iPrint, 1000)
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

               if (js. eq. -1)  then

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

            if (cmdbg  .and.  icmdbg(3) .gt. 0)
     $      write (iPrint, 1200) j, js, featol(j), res, atp, jadd1,
     $                         alfa1, jadd2, alfa2
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

            if (cmdbg  .and.  icmdbg(1) .gt. 0)
     $         write (iPrint, 9000) alfa, palfa1

            alfa = - min( abs( alfa ), palfa1 )
         end if
      end if

*     Test for undefined or infinite step.

      if (jadd .eq. 0) then
         alfa   = bigalf
         palfa  = bigalf
      end if

      if (alfa .ge. bigalf) inform = 3
      if (cmdbg  .and.  icmdbg(1) .gt. 0  .and.  inform .gt. 0)
     $   write (iPrint, 9010) jadd, alfa
      return

 1000 format(/ ' cmalf  entered'
     $       / '    j  js         featol        res             Ap',
     $         '     jadd1        alfa1     jadd2        alfa2 '/)
 1200 format( i5, i4, 3g15.5, 2(i6, g17.7) )
 9000 format(/ ' //cmalf //  negative step',
     $       / ' //cmalf //           alfa          palfa'
     $       / ' //cmalf //', 2g15.4 )
 9010 format(/ ' //cmalf //  unbounded step.'
     $       / ' //cmalf //  jadd           alfa'
     $       / ' //cmalf //  ', i4, g15.4 )

*     end of cmalf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmchk ( nerror, msglvl, lcrash, userkx,
     $                   liwork, lwork, litotl, lwtotl,
     $                   n, nclin, ncnln,
     $                   istate, kx, named, names,
     $                   bigbnd, bl, bu, x )

      implicit           double precision(a-h,o-z)
      character*8        names(*)
      logical            named, userkx
      integer            istate(n+nclin+ncnln), kx(n)
      double precision   bl(n+nclin+ncnln), bu(n+nclin+ncnln), x(n)

*     ==================================================================
*     cmchk   checks the data input to various optimizers.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written 10-May-1980.
*     Fortran 77 version written  5-October-1984.
*     This version of cmchk dated  21-Mar-93.       
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            ok
      intrinsic          abs
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
*     If warm start, check  istate.
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

*     end of  cmchk
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
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      integer            i, ierr, j, k
      intrinsic          abs

*     Check the parameters.

      if (m2 .lt. 1  .or.  m1 .lt. 1  .or.  m1 .gt. m2) then
         ierr = 1
         if (cmdbg  .and.  icmdbg(3) .gt. 0)
     $      write (iPrint, fmt=1100) m1, m2
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

 1100 format(/ ' //cmperm//  Illegal parameter values,'
     $       / ' //cmperm//    m1    m1'
     $       / ' //cmperm//', 2i6 )
 1200 format(/ ' XXX  kx(',I6,') contains an out-of-range value =', i16)
 1300 format(/ ' XXX  kx contains a duplicate value =',             i16)

*     end of cmperm
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmprt ( msglvl, nfree, ldA,
     $                   n, nclin, nctotl, bigbnd,
     $                   named, names,
     $                   nactiv, istate, kactiv, kx,
     $                   A, bl, bu, c, clamda, rlamda, x )

      implicit           double precision(a-h,o-z)
      character*8        names(*)
      logical            named
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   A(ldA,*), bl(nctotl), bu(nctotl), c(*),
     $                   clamda(nctotl), rlamda(n), x(n)

*     ==================================================================
*     cmprt   creates the expanded Lagrange multiplier vector clamda.
*     If msglvl .eq 1 or msglvl .ge. 10,  cmprt prints  x,  A*x,
*     c(x),  their bounds, the multipliers, and the residuals (distance
*     to the nearer bound).
*     cmprt is called by lscore and npcore just before exiting.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written  October 1984.
*     This version of  cmprt  dated  10-June-1986.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      character*2        ls, lstate(7)
      character*5        id(3), id3
      character*8        id4
      external           ddot
      intrinsic          abs

      parameter        ( zero  = 0.0d+0 )
      data               id(1) / 'varbl' /
      data               id(2) / 'lncon' /
      data               id(3) / 'nlcon' /
      data               lstate(1) / '--' /, lstate(2) / '++' /
      data               lstate(3) / 'fr' /, lstate(4) / 'll' /
      data               lstate(5) / 'ul' /, lstate(6) / 'eq' /
      data               lstate(7) / 'tb' /

      nplin  = n     + nclin
      nz     = nfree - nactiv

*     Expand multipliers for bounds, linear and nonlinear constraints
*     into the  clamda  array.

      call dload ( nctotl, zero, clamda, 1 )
      nfixed = n - nfree
      do 150, k = 1, nactiv+nfixed
         if (k .le. nactiv) j = kactiv(k) + n
         if (k .gt. nactiv) j = kx(nz+k)
         clamda(j) = rlamda(k)
  150 continue

      if (iPrint .eq. 0
     $    .or.  (msglvl .lt. 10  .and.  msglvl .ne. 1)) return

      write (iPrint, 1100)
      id3 = id(1)

      do 500, j = 1, nctotl
         b1     = bl(j)
         b2     = bu(j)
         wlam   = clamda(j)
         is     = istate(j)
         ls     = lstate(is + 3)
         if (j .le. n) then

*           Section 1 -- the variables  x.
*           ------------------------------
            k      = j
            v      = x(j)

         else if (j .le. nplin) then

*           Section 2 -- the linear constraints  A*x.
*           -----------------------------------------
            if (j .eq. n + 1) then
               write (iPrint, 1200)
               id3 = id(2)
            end if

            k      = j - n
            v      = ddot  ( n, A(k,1), ldA, x, 1 )
         else

*           Section 3 -- the nonlinear constraints  c(x).
*           ---------------------------------------------

            if (j .eq. nplin + 1) then
               write (iPrint, 1300)
               id3 = id(3)
            end if

            k      = j - nplin
            v      = c(k)
         end if

*        Print a line for the j-th variable or constraint.
*        -------------------------------------------------
         res    = v - b1
         res2   = b2 - v
         if (abs(res) .gt. abs(res2)) res = res2
         ip     = 1
         if (b1 .le. ( - bigbnd )) ip = 2
         if (b2 .ge.     bigbnd  ) ip = ip + 2
         if (named) then

            id4 = names(j)
            if (ip .eq. 1) then
               write (iPrint, 2100) id4,    ls, v, b1, b2, wlam, res
            else if (ip .eq. 2) then
               write (iPrint, 2200) id4,    ls, v,     b2, wlam, res
            else if (ip .eq. 3) then
               write (iPrint, 2300) id4,    ls, v, b1,     wlam, res
            else
               write (iPrint, 2400) id4,    ls, v,         wlam, res
           end if

         else

            if (ip .eq. 1) then
               write (iPrint, 3100) id3, k, ls, v, b1, b2, wlam, res
            else if (ip .eq. 2) then
               write (iPrint, 3200) id3, k, ls, v,     b2, wlam, res
            else if (ip .eq. 3) then
               write (iPrint, 3300) id3, k, ls, v, b1,     wlam, res
            else
               write (iPrint, 3400) id3, k, ls, v,         wlam, res
           end if
         end if
  500 continue
      return

 1100 format(// ' Variable        State', 5x, ' Value',
     $   6x, ' Lower bound', 4x, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 1200 format(// ' Linear constr   State', 5x, ' Value',
     $   6x, ' Lower bound', 4x, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 1300 format(// ' Nonlnr constr   State', 5x, ' Value',
     $   6x, ' Lower bound', 4x, ' Upper bound',
     $   '  Lagr multiplier', '     Residual' /)
 2100 format(1x, a8, 10x, a2, 3g16.7, g16.7, g16.4)
 2200 format(1x, a8, 10x, a2, g16.7, 5x, ' None', 6x, g16.7,
     $   g16.7, g16.4)
 2300 format(1x, a8, 10x, a2, 2g16.7, 5x, ' None', 6x, g16.7, g16.4)
 2400 format(1x, a8, 10x, a2,  g16.7, 5x, ' None', 11x, ' None',
     $   6x, g16.7, g16.4)
 3100 format(1x, a5, i3, 10x, a2, 3g16.7, g16.7, g16.4)
 3200 format(1x, a5, i3, 10x, a2,  g16.7,
     $   5x, ' None', 6x, g16.7, g16.7, g16.4)
 3300 format(1x, a5, i3, 10x, a2, 2g16.7, 5x, ' None', 6x,
     $   g16.7, g16.4)
 3400 format(1x, a5, i3, 10x, a2,  g16.7,
     $   5x, ' None', 11x, ' None', 6x, g16.7, g16.4)

*     end of  cmprt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmqmul( mode, n, nz, nfree, ldQ, unitQ,
     $                   kx, v, Q, w )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kx(n)
      double precision   v(n), Q(ldQ,*), w(n)

*     ==================================================================
*     cmqmul  transforms the vector  v  in various ways using the
*     matrix  Q = ( Z  Y )  defined by the input parameters.
*
*        MODE               result
*        ----               ------
*
*          1                v = Z v
*          2                v = Y v
*          3                v = Q v
*
*     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
*     On output, v  is a full n-vector.
*
*
*          4                v = Z'v
*          5                v = Y'v
*          6                v = Q'v
*
*     On input,  v  is a full n-vector.
*     On output, v  is ordered as  ( v(free)  v(fixed) ).
*
*          7                v = Y'v
*          8                v = Q'v
*
*     On input,  v  is a full n-vector.
*     On output, v  is as in modes 5 and 6 except that v(fixed) is not
*     set.
*
*     Modes  1, 4, 7 and 8  do not involve  v(fixed).
*     Original F66 version  April 1983.
*     Fortran 77 version written  9-February-1985.
*     Level 2 BLAS added 10-June-1986.
*     This version of cmqmul dated 14-Sep-92.
*     ==================================================================
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      nfixed = n - nfree
      j1     = 1
      j2     = nfree
      if (mode .eq. 1  .or.  mode .eq. 4) j2 = nz
      if (mode .eq. 2  .or.  mode .eq. 5  .or.  mode .eq. 7) j1 = nz + 1
      lenv   = j2 - j1 + 1
      if (mode .le. 3) then
*        ===============================================================
*        Mode = 1, 2  or  3.
*        ===============================================================

         if (nfree .gt. 0) call dload ( nfree, zero, w, 1 )

*        Copy  v(fixed)  into the end of  wrk.

         if (mode .ge. 2  .and.  nfixed .gt. 0)
     $      call dcopy ( nfixed, v(nfree+1), 1, w(nfree+1), 1 )

*        Set  W  =  relevant part of  Q * V.

         if (lenv .gt. 0)  then
            if (unitQ) then
               call dcopy ( lenv, v(j1), 1, w(j1), 1 )
            else
               call dgemv ( 'n', nfree, j2-j1+1, one, Q(1,j1), ldQ,
     $                      v(j1), 1, one, w, 1 )
            end if
         end if

*        Expand  w  into  v  as a full n-vector.

         call dload ( n, zero, v, 1 )
         do 220, k = 1, nfree
            j      = kx(k)
            v(j)   = w(k)
  220    continue

*        Copy  w(fixed)  into the appropriate parts of  v.

         if (mode .gt. 1)  then
            do 320, l = 1, nfixed
               j       = kx(nfree+l)
               v(j)    = w(nfree+l)
  320       continue
         end if

      else
*        ===============================================================
*        Mode = 4, 5, 6, 7  or  8.
*        ===============================================================
*        Put the fixed components of  v  into the end of  w.

         if (mode .eq. 5  .or.  mode .eq. 6)  then
            do 420, l = 1, nfixed
               j          = kx(nfree+l)
               w(nfree+l) = v(j)
  420       continue
         end if

*        Put the free  components of  v  into the beginning of  w.

         if (nfree .gt. 0)  then
            do 520, k = 1, nfree
               j      = kx(k)     
               w(k) = v(j)
  520       continue

*           Set  v  =  relevant part of  Q' * w.

            if (lenv .gt. 0)  then
               if (unitQ) then
                  call dcopy ( lenv, w(j1), 1, v(j1), 1 )
               else
                  call dgemv ( 'T', nfree, j2-j1+1, one, Q(1,j1), ldQ,
     $                         w, 1, zero, v(j1), 1 )
               end if
            end if
         end if

*        Copy the fixed components of  w  into the end of  v.

         if (nfixed .gt. 0  .and.  (mode .eq. 5  .or.  mode .eq. 6))
     $      call dcopy ( nfixed, w(nfree+1), 1, v(nfree+1), 1 )
      end if

*     end of  cmqmul
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmr1md( n, nu, nrank, ldR, lenv, lenw,
     $                   R, u, v, w, c, s )

      implicit           double precision(a-h,o-z)
      integer            n, nu, nrank, ldR, lenv, lenw
      double precision   R(ldR,*), u(n,*), v(n), w(n)
      double precision   c(n), s(n)

*     ==================================================================
*     cmr1md  modifies the  nrank*n  upper-triangular matrix  R  so that
*     Q*(R + v*w')  is upper triangular,  where  Q  is orthogonal,
*     v  and  w  are vectors, and the modified  R  overwrites the old.
*     Q  is the product of two sweeps of plane rotations (not stored).
*     If required,  the rotations are applied to the nu columns of
*     the matrix  U.
*
*     The matrix v*w' is an (lenv) by (lenw) matrix.
*     The vector v is overwritten.
*                                           
*     Systems Optimization Laboratory, Stanford University.
*     Original version   October  1984.
*     Level-2 matrix routines added 22-Apr-1988.
*     This version of  cmr1md  dated 22-Apr-1988.
*     ==================================================================
      intrinsic          min

      j = min( lenv, nrank )
      if (nrank .gt. 0) then
*        ---------------------------------------------------------------
*        Reduce  v to beta*e( j )  using a backward sweep of rotations
*        in planes (j-1, j), (j-2, j), ..., (1, j).
*        ---------------------------------------------------------------
         call f06fqf( 'Fixed', 'Backwards', j-1, v(j), v, 1, c, s )

*        ---------------------------------------------------------------
*        Apply the sequence of rotations to U.
*        ---------------------------------------------------------------
         if (nu .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Backwards', j, nu,
     $                   1, j, c, s, u, n )

*        ---------------------------------------------------------------
*        Apply the sequence of rotations to R. This generates a spike in
*        the j-th row of R, which is stored in s.
*        ---------------------------------------------------------------
         call f06qwf( 'Left', n, 1, j, c, s, R, ldR )

*        ---------------------------------------------------------------
*        Form  beta*e(j)*w' + R.  This a spiked matrix, with a row
*        spike in row j.
*        ---------------------------------------------------------------
         call daxpy( min( j-1, lenw ), v(j), w   , 1, s     , 1     )
         call daxpy( lenw-j+1        , v(j), w(j), 1, R(j,j), ldR )

*        ---------------------------------------------------------------
*        Eliminate the spike using a forward sweep of rotations in
*        planes (1, j), (2, j), ..., (j-1, j).
*        ---------------------------------------------------------------
         call f06qsf( 'Left', n, 1, j, c, s, R, ldR )

*        ---------------------------------------------------------------
*        Apply the rotations to U.
*        ---------------------------------------------------------------
         if (nu .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Forwards', j, nu,
     $                   1, j, c, s, u, n )
      end if

*     end of  cmr1md
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmrswp( n, nU, nrank, ldR, i, j, R, U, c, s )

      implicit           double precision(a-h,o-z)
      double precision   R(ldR,*), U(n,*)
      double precision   c(n), s(n)

*     ==================================================================
*     CMRSWP  interchanges the  i-th  and  j-th  (i .lt. j)  columns of
*     an  nrank x n  upper-trapezoidal matrix  R   and restores the
*     resulting matrix to upper-trapezoidal form using two sweeps of 
*     plane rotations applied on the left.  R is overwritten.  
*
*     If nU .gt. 0,  the rotations are applied to the  nU  columns of
*     the matrix  U.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level-2 matrix routines added 13-May-1988.
*     This version of  CMRSWP  dated  26-Aug-1991.
*     ==================================================================

      parameter        ( zero = 0.0d+0 )
      intrinsic          min

*     Swap the elements of the i-th and j-th columns of R on, or above,
*     the main diagonal.

      call dswap( min(i,nrank), R(1,i), 1, R(1,j), 1 )
      lenRj = min(j, nrank)

      if (lenRj .gt. i) then
*        ---------------------------------------------------------------
*        Reduce elements  R(i+1,j), ..., R(lenRj,j)  to  beta*e(lenRj)  
*        using a backward sweep in planes
*        (lenRj-1,lenRj), (lenRj-2,lenRj), ..., (i+1,lenRj).
*        If required, apply the sequence of rotations to U.
*        ---------------------------------------------------------------
         call f06fqf( 'Fixed', 'Backwards', lenRj-i-1, R(lenRj,j),
     $                R(i+1,j), 1, c(i+1), s(i+1) )

         if (nU .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Backwards', n, nU, 
     $                   i+1, lenRj, c, s, U, n )

*        Put zeros into the j-th column of R in positions corresponding 
*        to the sub-diagonals of the i-th column.

         s(i) = R(lenRj,j)
         call dload ( lenRj-i, zero, R(i+1,j), 1 )

*        Apply the sequence of rotations to R.  This generates a spike 
*        in the (lenRj)-th row of R, which is stored in s.

         call f06qwf( 'Left', n, i+1, lenRj, c, s, R, ldR )

*        Eliminate the spike using a forward sweep in planes
*        (i,lenRj), (i+1,lenRj), ..., (lenRj-1,lenRj).
*        If necessary, apply the sequence of rotations to U.

         call f06qsf( 'Left', n, i, lenRj, c, s, R, ldR )

         if (nU .gt. 0)
     $      call f06qxf( 'Left', 'Bottom', 'Forwards', lenRj, nU,
     $                   i, lenRj, c, s, U, n )
      end if

*     end of cmrswp
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

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     file  npcomsubs.f
*
*     npalf    npchkd   npfd     npfeas   npmrt    npsrch
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npalf ( inform, n, nclin, ncnln,
     $                   alfa, alfmin, alfmax, bigbnd, dxnorm,
     $                   Anorm, Adx, Ax, bl, bu,
     $                   dslk, dx, slk, x )

      implicit           double precision (a-h,o-z)
      double precision   Anorm(*), Adx(*), Ax(*), bl(*), bu(*),
     $                   dslk(*), dx(n), slk(*), x(n)

*     ==================================================================
*     npalf   finds a step alfa such that the point x + alfa*p reaches
*     one of the slacks or linear constraints.  The step alfa is the
*     maximum step that can be taken without violating one of the slacks
*     or linear constraints that is currently satisfied.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written  June 1986.
*     This version of npalf dated  13-Jun-1987.
*     ==================================================================
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      parameter         (zero = 0.0d+0, one = 1.0d+0)

      alfa   = alfmax
      j      = 1

*+    while (j .le. n+nclin+ncnln .and. alfa .gt. alfmin) do
  100 if    (j .le. n+nclin+ncnln .and. alfa .gt. alfmin) then

         if      (j .le. n      ) then
            Axi    =  x(j)
            Adxi   = dx(j)
            rownrm = one
         else if (j .le. n+nclin) then

            i      = j - n
            Axi    = Ax(i)
            Adxi   = Adx(i)
            rownrm = Anorm(i) + one
         else

            i      = j - n - nclin
            Axi    = slk(i)
            Adxi   = dslk(i)
            rownrm = one
         end if

         res = - one
         if (     Adxi .le. - epspt9*rownrm*dxnorm) then

*           Constraint decreasing.

            Adxi = - Adxi
            if (bl(j) .gt. -bigbnd) res = Axi   - bl(j)
         else if (Adxi .gt.   epspt9*rownrm*dxnorm) then

*           Constraint increasing.

            if (bu(j) .lt.  bigbnd) res = bu(j) - Axi
         end if

         if (res .gt. zero  .and.  alfa*Adxi .gt. res)
     $      alfa  = res / Adxi

         j = j + 1
         go to 100
*+    end while
      end if

*     ==================================================================
*     Determine alfa, the bound on the step to be taken.
*     ==================================================================
      alfa   = max( alfa, alfmin )

      inform = 0
      if (alfa .ge. alfmax) inform = 1

*     end of npalf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npchkd( inform, msgNP, nstate, lvlder, nfun, ngrad,
     $                   ldJ, ldJu, n, ncnln,
     $                   funcon, funobj, needc,
     $                   bigbnd, epsrf, cdint, fdint,
     $                   fdchk, fdnorm, objf, xnorm,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   dx, grad, gradu, hforwd, hcntrl,
     $                   x, wrk1, wrk2 )

      implicit           double precision (a-h,o-z)
      integer            needc(*)
      double precision   c(*), c1(*), cJac(ldJ,*), cJacu(ldJu,*),
     $                   cJdx(*)
      double precision   bl(n), bu(n), dx(n), grad(n), gradu(n), x(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   wrk1(n), wrk2(*)
      external           funcon, funobj

*     ==================================================================
*     npchkd  performs the following...
*     (1)  Computes the objective and constraint values objf and c.
*     (2)  Evaluates the user-provided gradients in cJacu and gradu.
*     (3)  Counts the missing gradients.
*     (4)  Loads the known gradients into grad and cJac.
*     (5)  Checks that the known gradients are programmed correctly.
*     (6)  Computes the missing gradient elements.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written 4-September-1985.
*     This version of npchkd dated 26-Nov-1989.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)

      logical            centrl, needfd
      parameter        ( rdummy =-11111.0d+0)

      infog  = 0
      infocj = 0
      nfdiff = 0
      ncdiff = 0
      ncset  = n*ncnln

      if (ncnln .gt. 0) then
*        ===============================================================
*        Compute the constraints and Jacobian matrix.
*        ===============================================================
*        If some derivatives are missing, load the Jacobian with dummy
*        values.  Any elements left unaltered after the call to funcon
*        must be estimated.  Any missing Jacobian elements are stored
*        in  cJacu.

         needfd = lvlder .eq. 0  .or.  lvlder .eq. 1

         if (needfd)
     $      call f06qhf( 'General', ncnln, n, rdummy, rdummy,
     $                   cJacu, ldJu )

         call iload ( ncnln, 1, needc, 1 )

         mode   = 2
         call funcon( mode, ncnln, n, ldJu,
     $                needc, x, c, cJacu, nstate )
         if (mode .lt. 0) go to 9999

         call f06qff( 'General', ncnln, n, cJacu, ldJu, cJac, ldJ )

         if (needfd) then

*           Count the number of missing Jacobian elements.

            do 220, j = 1, n
               do 210, i = 1, ncnln
                  if (cJacu(i,j) .eq. rdummy) ncdiff = ncdiff + 1
  210          continue
  220       continue

            ncset = ncset - ncdiff
            if (nstate .eq. 1) then
               if (ncdiff .eq. 0) then
                  if (lvlder .eq. 0) lvlder = 2
                  if (lvlder .eq. 1) lvlder = 3
                  if (msgNP .gt. 0  .and.  iPrint .gt. 0)
     $               write(iPrint, 1000) lvlder
               else
                  if (msgNP .gt. 0  .and.  iPrint .gt. 0)
     $               write(iPrint, 1100) ncset, n*ncnln, ncdiff
               end if
            end if
         end if
      end if

*     ==================================================================
*     Repeat the procedure above for the objective function.
*     ==================================================================
      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2

      if (needfd)
     $   call dload ( n, rdummy, gradu, 1 )

      mode  = 2
      call funobj( mode, n, x, objf, gradu, nstate )
      if (mode .lt. 0) go to 9999

      call dcopy ( n, gradu, 1, grad, 1 )

      if (needfd) then

*        Count the number of missing gradient elements.

         do 300, j = 1, n
            if (gradu(j) .eq. rdummy) nfdiff = nfdiff + 1
  300    continue

         if (nstate .eq. 1) then
            if (nfdiff .eq. 0) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msgNP .gt. 0  .and.  iPrint .gt. 0)
     $            write(iPrint, 2000) lvlder
            else
               if (msgNP .gt. 0  .and.  iPrint .gt. 0)
     $            write(iPrint, 2100) n - nfdiff, n, nfdiff
            end if
         end if
      end if

      nfun  = nfun  + 1
      ngrad = ngrad + 1

*     ==================================================================
*     Check any gradients that have been provided.
*     ==================================================================
      if (lvrfyc .ge. 0) then
         if (ncset .gt. 0) then
            call chcJac( mode, lvlder, msgNP,
     $                   ncset, n, ncnln, ldJ, ldJu,
     $                   bigbnd, epsrf, epspt3, fdchk, xnorm,
     $                   funcon, needc,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   dx, wrk2, x, wrk1 )
            if (mode .lt. 0) go to 9999
            infocj = mode
         end if

         if (nfdiff .lt. n) then
            call chfgrd( mode, msgNP, n,
     $                   bigbnd, epsrf, epspt3, fdchk, objf, xnorm,
     $                   funobj,
     $                   bl, bu, grad, gradu, dx, x, wrk1 )
            if (mode .lt. 0) go to 9999
            infog   = mode
         end if
      end if

      needfd = ncdiff .gt. 0  .or.  nfdiff .gt. 0
      if (needfd) then
*        ===============================================================
*        Compute the missing gradient elements.
*        ===============================================================
*        chfd computes the finite-difference intervals or checks
*        preassigned (scalar) values.

         call chfd  ( mode, msgNP, lvlder,
     $                n, ncnln, ldJ, ldJu,
     $                bigbnd, epsrf, fdnorm, objf,
     $                funobj, funcon, needc,
     $                bl, bu, c, c1, cJdx, cJac, cJacu,
     $                grad, gradu, hforwd, hcntrl, x, dx )

         if (mode .lt. 0) go to 9999

         if (lfdset .gt. 0) then
            centrl = lvldif .eq. 2
            call npfd  ( centrl, mode,
     $                   ldJ, ldJu, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm, objf,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, cJdx, cJac, cJacu,
     $                   grad, gradu, hforwd, hcntrl, x )

            if (mode .lt. 0) go to 9999
         end if
      end if

      inform = infocj + infog

      return

*     The user requested termination.

 9999 inform = mode
      return

 1000 format(//' All Jacobian elements have been set.  ',
     $         ' Derivative level increased to ', i4 )
 1100 format(//' The user sets ', i6, '   out of', i6,
     $         '   Jacobian elements.'
     $       / ' Each iteration, ', i6,
     $         '   Jacobian elements will be estimated numerically.' )
 2000 format(//' All objective gradient elements have been set.  ',
     $         ' Derivative level increased to ', i4 )
 2100 format(//' The user sets ', i6, '   out of', i6,
     $         '   objective gradient elements.'
     $       / ' Each iteration, ', i6,
     $         '   gradient elements will be estimated numerically.' )

*     end of npchkd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfd  ( centrl, inform,
     $                   ldJ, ldJu, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm, objf,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, c2, cJac, cJacu,
     $                   grad, gradu, hforwd, hcntrl, x )

      implicit           double precision (a-h,o-z)
      logical            centrl
      integer            needc(*)

      double precision   bl(n), bu(n), c(*), c1(*), c2(*),
     $                   cJac(ldJ,*), cJacu(ldJu,*)
      double precision   grad(n), gradu(n), hforwd(n), hcntrl(n), x(n)
      external           funcon, funobj

*     ==================================================================
*     npfd   evaluates any missing gradients.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written 3-July-1986.
*     This version of npfd   dated 14-Jul-94.
*     ==================================================================
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      parameter         (rdummy=-11111.0d+0)
      parameter         (zero  = 0.0d+0, half  = 0.5d+0, one   = 1.0d+0)
      parameter         (three = 3.0d+0, four  = 4.0d+0                )

      inform = 0

*     ==================================================================
*     Use the pre-assigned difference intervals to approximate the
*     derivatives.
*     ==================================================================
*     Use either the same interval for each component (lfdset = 1),
*     or the intervals already in hforwd or hcntrl (lfdset = 0 or 2).

      nstate =   0
      mode   =   0

      biglow = - bigbnd
      bigupp =   bigbnd

      fdnorm =   zero

      do 340, j = 1, n
         xj     = x(j)
         nColj  = 0
         if (ncdiff .gt. 0) then
            do 310, i = 1, ncnln
               if (cJacu(i,j) .eq. rdummy) then
                  needc(i) = 1
                  nColj    = nColj + 1
               else
                  needc(i) = 0
               end if
  310       continue
         end if

         if (nColj .gt. 0  .or.  gradu(j) .eq. rdummy) then
            stepbl = biglow
            stepbu = bigupp
            if (bl(j) .gt. biglow) stepbl = bl(j) - xj
            if (bu(j) .lt. bigupp) stepbu = bu(j) - xj

            if (centrl) then
               if (lfdset .eq. 1) then
                  delta = cdint
               else
                  delta = hcntrl(j)
               end if
            else
               if (lfdset .eq. 1) then
                  delta = fdint
               else
                  delta = hforwd(j)
               end if
            end if

            delta  = delta*(one + abs(xj))
            fdnorm = max (fdnorm, delta)
            if (half*(stepbl + stepbu) .lt. zero) delta =  - delta

            x(j) = xj + delta
            if (nColj .gt. 0) then
               call funcon( mode, ncnln, n, ldJu,
     $                      needc, x, c1, cJacu, nstate )
               if (mode .lt. 0) go to 999
            end if

            if (gradu(j) .eq. rdummy) then
               call funobj( mode, n, x, objf1, gradu, nstate )
               if (mode .lt. 0) go to 999
            end if

            if (centrl) then
*              ---------------------------------------------------------
*              Central differences.
*              ---------------------------------------------------------
               x(j)  = xj + delta + delta

               if (nColj .gt. 0) then
                  call funcon( mode, ncnln, n, ldJu,
     $                         needc, x, c2, cJacu, nstate )
                  if (mode .lt. 0) go to 999

                  do 320, i = 1, ncnln
                     if (needc(i) .eq. 1)
     $                  cJac(i,j) = (four*c1(i) - three*c(i) - c2(i))
     $                                  / (delta + delta)
  320             continue
               end if

               if (gradu(j) .eq. rdummy) then
                  call funobj( mode, n, x, objf2, gradu, nstate )
                  if (mode .lt. 0) go to 999

                  grad(j) = (four*objf1 - three*objf - objf2)
     $                                  / (delta + delta)

               end if
            else
*              ---------------------------------------------------------
*              Forward Differences.
*              ---------------------------------------------------------
               if (nColj .gt. 0) then
                  do 330, i = 1, ncnln
                     if (needc(i) .eq. 1)
     $                  cJac(i,j) = (c1(i) -  c(i))/  delta
  330             continue
               end if

               if (gradu(j) .eq. rdummy)
     $            grad(j) = (objf1 - objf) /  delta

            end if
         end if
         x(j) = xj

  340 continue

      return

  999 inform = mode
      return

*     end of npfd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfeas( n, nclin, ncnln, istate,
     $                   bigbnd, cvnorm, errmax, jmax, nviol,
     $                   Ax, bl, bu, c, featol, x, work )

      implicit           double precision (a-h,o-z)
      integer            istate(n+nclin+ncnln)
      double precision   Ax(*), bl(n+nclin+ncnln), bu(n+nclin+ncnln)
      double precision   c(*), featol(n+nclin+ncnln), x(n)
      double precision   work(n+nclin+ncnln)

*     ==================================================================
*     npfeas  computes the following...
*     (1)  The number of constraints that are violated by more
*          than  featol  and the 2-norm of the constraint violations.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version      April    1984.
*     This version of  npfeas  dated  16-October-1985.
*     ==================================================================
      parameter        ( zero = 0.0d+0 )

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Compute nviol, the number of constraints violated by more than
*     featol,  and cvnorm,  the 2-norm of the constraint
*     violations and residuals of the constraints in the qp working set.
*     ==================================================================
      nviol  = 0

      do 200, j = 1, n+nclin+ncnln
         feasj  = featol(j)
         res    = zero

         if (j .le. n + nclin) then

*           Bound or general linear constraint.

            if (j .le. n) then
               con =  x(j)
            else
               con = Ax(j-n)
            end if

            tolj   = feasj
         else

*           Nonlinear constraint.

            con    = c(j-n-nclin)
            tolj   = zero
         end if

*        Check for constraint violations.

         if (bl(j) .gt. biglow) then
            res    = bl(j) - con
            if (res .gt.   feasj ) nviol = nviol + 1
            if (res .gt.    tolj ) go to 190
         end if

         if (bu(j) .lt. bigupp) then
            res    = bu(j) - con
            if (res .lt. (-feasj)) nviol = nviol + 1
            if (res .lt.  (-tolj)) go to 190
         end if

*        This constraint is satisfied,  but count the residual as a
*        violation if the constraint is in the working set.

         is     = istate(j)

         if (is .eq. 0) then
            res = zero
         else if (is .eq. 1  .or.  is .le. -2) then
            res = bl(j) - con
         else if (is .ge. 2  .or.  is .eq. -1) then
            res = bu(j) - con
         end if

         if (abs( res ) .gt. feasj) nviol = nviol + 1

*        Set the array of violations.

  190    work(j) = res
  200 continue

      jmax   = idamax( n+nclin+ncnln, work, 1 )
      errmax = abs ( work(jmax) )
      cvnorm = dnrm2 ( n+nclin+ncnln, work, 1 )

*     end of npfeas
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npmrt ( feasqp, n, nclin, ncnln,
     $                   objalf, grdalf, qpcurv,
     $                   istate,
     $                   cJdx, cmul, cs,
     $                   dlam, Pen, violn,
     $                   work1, work2 )

      implicit           double precision (a-h,o-z)

      logical            feasqp

      integer            istate(*)

      double precision   cJdx(*), cmul(*), cs(*),
     $                   dlam(*), Pen(*), violn(*)
      double precision   work1(*), work2(*)

*     ==================================================================
*     npmrt   computes the value and directional derivative of the
*     augmented Lagrangian merit function.  The penalty parameters
*     Pen(j) are boosted if the directional derivative of the resulting
*     augmented Lagrangian function is not sufficiently negative.  If
*     Pen needs to be increased,  the perturbation with minimum two-norm
*     is found that gives a directional derivative equal to  - p'Hp.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written  27-May-1985.
*     This version of  NPMRT  dated 14-November-1985.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      logical            incrun
      common    /sol6np/ PenMax, PenNrm, PenDmp, PenScl, incrun

      logical            boost , overfl
      parameter        ( zero   = 0.0d+0, half = 0.5d+0, one = 1.0d+0 )
      parameter        ( two    = 2.0d+0                              )

      if (ncnln .eq. 0) return

      rtmin  = wmach(6)

      objalf = objalf - ddot  ( ncnln, cmul, 1, cs, 1 )
      grdalf = grdalf - ddot  ( ncnln, dlam, 1, cs, 1 )

      call dcopy ( ncnln, cs, 1, work1, 1 )

      if (.not. feasqp) then
         nplin  = n + nclin

         do 100 i = 1, ncnln
            if (istate(nplin+i) .lt. 0  .or.  violn(i) .ne. zero)
     $         work1(i) = - cJdx(i)
  100    continue
      end if

      grdalf = grdalf + ddot  ( ncnln, work1, 1, cmul, 1 )

      if (feasqp) then

*        Find the quantities that define  PenMin, the vector of minimum
*        two-norm such that the directional derivative is one half of
*        approximate curvature   - (dx)'H(dx).

         do 350, i = 1, ncnln
            if (abs( cs(i) ) .le. rtmin) then
               work2(i) = zero
            else
               work2(i) = cs(i)**2
            end if
  350    continue

         qnorm  = dnrm2 ( ncnln, work2, 1 )
         tscl   = ddiv  ( grdalf + half*qpcurv, qnorm, overfl )
         if (abs( tscl ) .le. PenMax  .and.  .not. overfl) then
*           ------------------------------------------------------------
*           Bounded  PenMin  found.  The final value of  Pen(J)  will
*           never be less than  PenMin(j).  If the  QP  was feasible,  a
*           trial value  PenNew  is computed that is equal to the
*           geometric mean of the previous  Pen  and a damped value of
*           PenMin.  The new  Pen  is defined as  PenNew  if it is less
*           than half the previous  Pen  and greater than  PenMin.
*           ------------------------------------------------------------
            PenScl  = one
            do 400, i = 1, ncnln
               PenMin = max(  (work2(i)/qnorm)*tscl, zero )
               Peni   = Pen(i)

               PenNew = sqrt( Peni*(PenDmp + PenMin) )
               if (PenNew .lt. half*Peni  ) Peni = PenNew
               if (Peni   .lt.      PenMin) Peni = PenMin
               Pen(i) = Peni
  400       continue

            Pen1   = PenNrm
            PenNrm = dnrm2 ( ncnln, Pen, 1 )

*           ------------------------------------------------------------
*           If  incrun = true,  there has been a run of iterations in
*           which the norm of  Pen  has not decreased.  Conversely,
*           incrun = false  implies that there has been a run of
*           iterations in which the norm of Pen has not increased.  If
*           incrun changes during this iteration the damping parameter
*           PenDmp is increased by a factor of two.  This ensures that
*           Pen(j) will oscillate only a finite number of times.
*           ------------------------------------------------------------
            boost  = .false.
            if (      incrun  .and.  PenNrm .lt. Pen1) boost = .true.
            if (.not. incrun  .and.  PenNrm .gt. Pen1) boost = .true.
            if (boost) then
               PenDmp = two*PenDmp
               incrun = .not. incrun
            end if
         end if
      else

*        The  QP  was infeasible.  Do not alter the penalty parameters,
*        but compute the scale factor so that the constraint violations
*        are reduced.

         call ddscl ( ncnln, Pen, 1, work1, 1 )
         pterm2 = ddot  ( ncnln, work1, 1, cs, 1 )

         PenScl  = PenMax
         tscl   = ddiv  ( grdalf, pterm2, overfl )
         if (tscl .gt. PenScl  .and.  tscl .le. PenMax/(one+PenNrm)
     $                        .and.  .not. overfl)
     $      PenScl = tscl

         call dcopy ( ncnln, cs, 1, work1, 1 )
      end if

*     ------------------------------------------------------------------
*     Compute the new value and directional derivative of the
*     merit function.
*     ------------------------------------------------------------------
      call ddscl ( ncnln, Pen, 1, work1, 1 )

      pterm  = ddot  ( ncnln, work1, 1, cs, 1 )
      objalf = objalf + half*PenScl*pterm

      if (feasqp)
     $  pterm2 = pterm

      grdalf = grdalf -      PenScl*pterm2

*     end of npmrt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsrch( needfd, inform, n, ncnln,
     $                   ldJ, ldJu, nfun, ngrad,
     $                   needc, funcon, funobj,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, gL1, gL,
     $                   objf, objalf, curvQP, xnorm,
     $                   c, c2, cJac, cJacu, cJdx, cJdx2,
     $                   cmul1, cmul, cs1, cs, dx, dlam, dslk,
     $                   grad, gradu, qpmul, Pen,
     $                   slk1, slk, x1, x, work )

      implicit           double precision (a-h,o-z)
      logical            needfd
      integer            needc(*)
      double precision   dx(n), grad(n), gradu(n), x1(n), x(n)
      double precision   c(*), c2(*), cJac(ldJ,*), cJacu(ldJu,*),
     $                   cJdx(*), cJdx2(*), cmul1(*), cmul(*)
      double precision   cs1(*), cs(*), dlam(*), dslk(*), qpmul(*),
     $                   Pen(*), slk1(*), slk(*), work(*)
      external           funobj, funcon

*     ==================================================================
*     npsrch finds the steplength alfa that gives sufficient decrease in
*     the augmented Lagrangian merit function.
*
*     On exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
*     with an associated merit function value  objalf  which is lower
*     than that at the base point. If  inform = 4, 5, 6, 7 or 8,  alfa
*     is zero and  objalf will be the merit value at the base point.
*
*     Original version written  27-May-1985.
*     12-Jun-1986: Level 2 BLAS added 
*     13-Oct-1994: MXS added backtrack if functions return mode = -1.
*     This version of npsrch dated  17-Sep-95.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      logical            incrun
      common    /sol6np/ PenMax, PenNrm, PenDmp, PenScl, incrun

      logical            debug , done  , first , imprvd
      parameter        ( zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0 )
      parameter        ( two    =2.0d+0                                )
      parameter        ( tolg   =1.0d-1                                )
      parameter        ( rmu    =1.0d-4                                )

      epsmch = wmach(3)

      if (.not. needfd  .and.  ncnln .gt. 0)
     $   cs1Jdx = ddot( ncnln, cs1, 1, cJdx, 1 )

*     ------------------------------------------------------------------
*     Set the input parameters and tolerances for srchc and srchq.
*
*     tolrx   is the tolerance on relative changes in dx resulting from
*             changes in alfa.
*
*     tolax   is the tolerance on absolute changes in dx resulting from
*             changes in alfa.
*
*     tolabs  is the tolerance on absolute changes in alfa.
*
*     tolrel  is the tolerance on relative changes in alfa.
*
*     toltny  is the magnitude of the smallest allowable value of alfa.
*             if  M(tolabs) - M(0) .gt. epsaf,  the line search tries
*             steps in the range  toltny .le. alfa .le. tolabs.
*     ------------------------------------------------------------------
      nstate = 0
      debug  = .false.

      if (needfd) then
         maxf = 15
      else
         maxf = 10
      end if

      epsaf  = epsrf*(one + abs( objalf ))
      tolax  = epspt8
      tolrx  = epspt8

      if (tolrx*xnorm + tolax .lt. dxnorm*alfmax) then
         tolabs = (tolrx*xnorm + tolax) /  dxnorm
      else
         tolabs = alfmax
      end if
      tolrel = max( tolrx , epsmch )

      t      = zero
      do 10, j = 1, n
         s = abs( dx(j) )
         q = abs( x(j) )*tolrx + tolax
         if (s .gt. t*q) t = s / q
   10 continue

      if (t*tolabs .gt. one) then
         toltny = one / t
      else
         toltny = tolabs
      end if

      oldf   = objalf
      oldg   = grdalf
      
   50 alfbst = zero
      fbest  = zero
      gbest  = (one - rmu)*oldg
      targtg = (rmu - eta)*oldg
      g0     = gbest

      if (ncnln .gt. 0) call iload ( ncnln, 1, needc, 1 )

      if (needfd) then
         mode = 0
      else
         mode = 2
      end if

      first  = .true.

*     ------------------------------------------------------------------
*     Commence main loop, entering srchc or srchq two or more times.
*     first = true for the first entry,  false for subsequent entries.
*     done  = true indicates termination, in which case the value of
*     inform gives the result of the search.
*     inform = 1 if the search is successful and alfa < alfmax.
*            = 2 if the search is successful and alfa = alfmax.
*            = 3 if a better point was found but too many functions
*                were needed (not sufficient decrease).
*            = 4 if alfmax < tolabs (too small to do a search).
*            = 5 if alfa < alfsml (srchq only -- maybe want to switch
*                to central differences to get a better direction).
*            = 6 if the search found that there is no useful step.
*                The interval of uncertainty is less than 2*tolabs.
*                The minimizer is very close to alfa = zero
*                or the gradients are not sufficiently accurate.
*            = 7 if there were too many function calls.
*            = 8 if the input parameters were bad
*                (alfmax le toltny  or  oldg ge 0).
*     ------------------------------------------------------------------
*+    repeat
  100    if (needfd) then
            call srchq ( first , debug , done  , imprvd, inform,
     $                   maxf  , numf  , iPrint, 
     $                   alfmax, alfsml, epsaf , 
     $                   g0    , targtg, ftry  ,         
     $                   tolabs, tolrel, toltny,
     $                   alfa  , alfbst, fbest  )
         else
            call srchc ( first , debug , done  , imprvd, inform,
     $                   maxf  , numf  , iPrint,
     $                   alfmax,         epsaf , 
     $                   g0    , targtg, ftry  , gtry  , 
     $                   tolabs, tolrel, toltny,
     $                   alfa  , alfbst, fbest , gbest )
         end if

         if (imprvd) then
            objf   = tobj
            objalf = tobjM

            if (ncnln .gt. 0)
     $         call dcopy ( ncnln, c2, 1, c, 1 )

            if (.not. needfd) then
               call dcopy ( n, gradu, 1, grad, 1 )
               gdx = tgdx
               gL  = tgL

               if (ncnln .gt. 0) then
                  call dcopy ( ncnln, cJdx2, 1, cJdx, 1 )
                  call f06qff( 'General', ncnln, n, cJacu, ldJu,
     $                         cJac, ldJ )
               end if
            end if
         end if

*        ---------------------------------------------------------------
*        If done = false,  the problem functions must be computed for
*        the next entry to srchc or srchq.
*        If done = true,   this is the last time through.
*        ---------------------------------------------------------------
         if (.not. done) then
            call dcopy ( n,       x1, 1, x, 1 )
            call daxpy ( n, alfa, dx, 1, x, 1 )

            if (ncnln .gt. 0) then

*              Compute the new estimates of the multipliers and slacks.
*              If the step length is greater than one,  the multipliers
*              are fixed as the QP-multipliers.

               if (alfa .le. one) then
                  call dcopy ( ncnln,       cmul1, 1, cmul, 1 )
                  call daxpy ( ncnln, alfa, dlam , 1, cmul, 1 )
               end if
               call dcopy ( ncnln,       slk1, 1, slk, 1 )
               call daxpy ( ncnln, alfa, dslk, 1, slk, 1 )

*              ---------------------------------------------------------
*              Compute the new constraint vector and Jacobian.
*              ---------------------------------------------------------
               call funcon( mode, ncnln, n, ldJu,
     $                      needc, x, c2, cJacu, nstate )
               if (mode .lt. 0) go to 999

               call dcopy ( ncnln,         c2 , 1, cs, 1 )
               call daxpy ( ncnln, (-one), slk, 1, cs, 1 )

               call dcopy ( ncnln, cs , 1,  work, 1 )
               call ddscl ( ncnln, Pen, 1,  work, 1 )

               fterm  =             ddot( ncnln, cmul, 1, cs, 1 ) -
     $                  half*PenScl*ddot( ncnln, work, 1, cs, 1 )
            end if

*           ------------------------------------------------------------
*           Compute the value and gradient of the objective function.
*           ------------------------------------------------------------
            call funobj( mode, n, x, tobj, gradu, nstate )
            if (mode .lt. 0) go to 999

            if (ncnln .gt. 0) then
               tobjM = tobj  - fterm
            else
               tobjM = tobj
            end if

            ftry  = tobjM - oldf - rmu*oldg*alfa 

*           ------------------------------------------------------------
*           Compute auxiliary gradient information.
*           ------------------------------------------------------------
            if (.not. needfd) then
               gtry = ddot( n, gradu, 1, dx, 1 )
               tgdx = gtry
               tgL  = gtry
               if (ncnln .gt. 0) then

*                 Compute the Jacobian times the search direction.

                  call dgemv ( 'n', ncnln, n, one, cJacu, ldJu, dx, 1,
     $                         zero, cJdx2, 1 )

                  call dcopy ( ncnln,         cJdx2, 1, work, 1 )
                  call daxpy ( ncnln, (-one), dslk , 1, work, 1 )

                  gtry   = gtry - ddot( ncnln, cmul, 1, work, 1 )
                  if (alfa .le. one)
     $               gtry   = gtry - ddot( ncnln, dlam, 1, cs      , 1 )

                  call ddscl ( ncnln, Pen , 1, work, 1 )
                  gtry = gtry + PenScl*ddot( ncnln, work , 1, cs   , 1 )

                  tgL  = tgdx -        ddot( ncnln, cJdx2, 1, qpmul, 1 )

*                 ------------------------------------------------------
*                 If alfbnd .le. alfa .lt. alfmax and the norm of the
*                 quasi-Newton update is bounded, set alfmax to be alfa.
*                 This will cause the line search to stop if the merit
*                 function is decreasing at the boundary.
*                 ------------------------------------------------------
                  if (alfbnd .le. alfa  .and.  alfa .lt. alfmax) then
                     csJdx  = ddot   ( ncnln, cs, 1, cJdx2, 1 )
                     curvL  = tgL - gL1
                     curvc  = abs( csJdx - cs1Jdx )
                     PenBFS = max( curvQP*tolg - curvL, zero )
                     if (PenBFS .le. curvc*PenMax) then
                        alfmax = alfa
                     else
                        alfbnd = min( two*alfa, alfmax )
                     end if
                  end if
               end if

               gtry = gtry  - rmu*oldg

            end if
         end if
*+    until (      done)
      if    (.not. done) go to 100

      nfun = nfun + numf
      if (.not. needfd) ngrad = ngrad + numf
      alfa = alfbst

      if (.not. imprvd) then
         call dcopy ( n,       x1, 1, x, 1 )
         call daxpy ( n, alfa, dx, 1, x, 1 )
         if (ncnln .gt. 0) then
            if (alfa .le. one) then
               call dcopy ( ncnln,       cmul1, 1, cmul, 1 )
               call daxpy ( ncnln, alfa, dlam , 1, cmul, 1 )
            end if
            call dcopy ( ncnln,         slk1 , 1, slk, 1 )
            call daxpy ( ncnln,   alfa, dslk , 1, slk, 1 )
            call dcopy ( ncnln,         c    , 1, cs , 1 )
            call daxpy ( ncnln, (-one), slk  , 1, cs , 1 )
         end if
      end if

      return

*     ------------------------------------------------------------------
*     One of the function routines set mode < 0.
*     mode = -1 means "undefined function at step alfa".
*     We start the whole search again with a smaller alfmax and alfa.
*     (The search routines will say if this is too small.)
*     ------------------------------------------------------------------
  999 if (mode .eq. -1) then
         alfmax = 0.1d+0 * alfa
         alfa   = alfmax
         go to 50
      end if

*     The user wants to stop.  Who am I to object?

      inform = mode
      return

*     end of npsrch
      end


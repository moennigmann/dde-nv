*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  nlssolsubs.f
*
*     nlssol   nlchkd   nlcore   nldflt   nlfd    nlfile    nlloc
*     nlJtJ    nlkey    nlnkey   nlobjf   nloptn  nlopti    nloptr
*     nlsrch
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlssol( m, n, nclin, ncnln, ldA, ldcJu, ldfJu, ldR,
     $                   A, bl, bu,
     $                   funcon, funobj,
     $                   inform, iter, istate,
     $                   c, cJacu, y, f, fJacu, clamda, objf, R, x,
     $                   iw, leniw, w, lenw )

      implicit           double precision(a-h,o-z)
      external           funcon, funobj
      integer            istate(n+nclin+ncnln)
      integer            iw(leniw) 
      double precision   A(ldA,*), bl(n+nclin+ncnln),
     $                   bu(n+nclin+ncnln)
      double precision   c(*), cJacu(ldcJu,*), clamda(n+nclin+ncnln)
      double precision   y(m), f(m), fJacu(ldfJu,*)
      double precision   R(ldR,*), x(n)
      double precision   w(lenw)

*     ==================================================================
*     nlssol  solves the constrained least-squares problem
*
*               minimize           (1/2)(y - f(x))'(y - f(x))
*
*                                       (    x  )
*               subject to    bl  .le.  (  A*x  )  .le.  bu
*                                       (  c(x) )
*
*     where  y  is a given constant m-vector,  f(x)  is an m-vector of 
*     smooth functions,  A  is a constant matrix and  c(x)  is a vector
*     of smooth nonlinear functions.  The feasible region is defined 
*     by a mixture of linear and nonlinear equality or inequality
*     constraints on  x.
*
*     The dimensions of the problem are...
*
*     m        the number of data points (the dimension of y). 
*
*     n        the number of variables (dimension of  x),
*
*     nclin    the number of linear constraints (rows of the matrix  A),
*
*     ncnln    the number of nonlinear constraints (dimension of  c(x)),
*
*
*     nlssol  uses a sequential quadratic programming algorithm, with a
*     positive-definite quasi-Newton approximation to the transformed
*     Hessian  Q'HQ  of the Lagrangian function (which will be stored in
*     the array  R).
* 
*     Version 4.04, June        9, 1989. Matches NPSOL Version 4.04.
*     Version 4.06, November    5, 1991. Matches NPSOL Version 4.06.
*     Version 5.01, July       12, 1994. Debug printing eliminated.
*                                        Data vector y added.
*     Version 5.02, September  15, 1995. Printing revamped.
*
*     This version of  NLSSOL  dated 15-Sep-95.
*     Copyright  1988--1995  Optimates.
*     This software is not in the public domain. Its use is governed by
*     a license agreement with Optimates.  It is illegal to 
*     make copies except as authorized by the license agreement.
*     ==================================================================

*     Common blocks.

*     +Include lsparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   itmax1, itmax2, lcrash, lformH, lprob , msgLS ,
     $                   nn    , nnclin, nprob , ipadls(21)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolOpt, tolrnk, rpadls(22)

      equivalence       (iprmls(1), itmax1 ), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, lvlder, 
     $                   lverfy, msgNP , nlnf  , nlnj  , nlnx  , nncnln,
     $                   nsave , nload , ksave , ipadnp(15)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , Hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), itmxnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     +Include nlparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnl(mxparm), ipsvnl
      double precision   rprmnl(mxparm), rpsvnl

      common    /nlpar1/ ipsvnl(mxparm),
     $                   ltypeH, nreset, ipadnl(28)

      common    /nlpar2/ rpsvnl(mxparm), rpadnl(30)

      equivalence       (iprmnl(1), ltypeH), (rprmnl(1), rpadnl(1))

      save      /nlpar1/, /nlpar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (itmxnp, nmajor), (itmax2, nminor), (msgLS , msgQP )

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol3cm/ lennam, ldT   , ncolT , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize , dTmax , dTmin
      common    /sol6cm/ Rcndbd, Rfrobn, dRmax , dRmin

      logical            unitQ
      common    /sol1sv/ nactiv, nfree , nZ   , unitQ
      save      /sol1sv/

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)
      logical            incrun
      common    /sol6np/ Penmax, Pennrm, Pendmp, scale, incrun

      parameter         (lennl = 20)
      common    /sol1nl/ locnl(lennl)

*     Local variables.

      character*16       names(1)
      logical            cold  , linobj, named , overfl, rowerr, vertex
      logical            needfd
      parameter         (zero   =0.0d+0, point3 =3.3d-1, point8 =0.8d+0)
      parameter         (point9 =0.9d+0, one    =1.0d+0)
      parameter         (hundrd =1.0d+2, growth =1.0d+2                )

      character*40       title
      data               title
     $                 / 'NLSSOL ---  Version 5.0-2      Sep  1995' /
*                         1234567890123456789012345678901234567890

*     Set the machine-dependent constants.

      call mchpar()

      epsmch = wmach( 3)
      rteps  = wmach( 4)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      Penmax = one/epsmch
      rootn  = sqrt(dble(n))

*     Default names will be provided for variables during printing.

      named  = .false.
      inform = 0
      iter   = 0

*     Set the default values for the parameters.

      call nldflt( n, nclin, ncnln, title )

      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or. (lvlder .eq. 1  .and.  ncnln .gt. 0)
      cold   = lcrash .eq. 0
      lvldif = 0
      if (needfd) lvldif = 1

      nplin  = n     + nclin
      nctotl = nplin + ncnln

*     Assign the dimensions of arrays in the parameter list of nlcore.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter minact and minfxd
*     accordingly.

      minact = 0
      minfxd = 0

      mxfree = n - minfxd
      maxact = max( 1, min( n, nclin ) )
      maxnZ  = n - ( minfxd + minact )

      if (nclin + ncnln .eq. 0) then
         ldQ   = 1
         ldT   = 1
         ncolT = 1
      else
         ldQ   = max( 1, mxfree )
         ldT   = max( maxnZ, maxact )
         ncolT = mxfree
      end if

      lennam = 1

      ldAqp = max( nclin+ncnln, 1 )

*     nlloc  defines the arrays that contain the locations of various
*     work arrays within  w  and  iw.

      litotl = 0
      lwtotl = 0
      call nploc(    n, nclin, ncnln, nctotl, litotl, lwtotl)
      call nlloc( m, n,                       litotl, lwtotl)

      lkactv = locls( 1)
      lanorm = locls( 2)
      lcJdx  = locls( 3)
      lres   = locls( 5)
      lres0  = locls( 6)
      lgq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk1  = locls(14)

      lkx    = locnp( 1)
      liperm = locnp( 2)
      lAqp   = locnp( 3)
      ldx    = locnp( 7)
      lfeatl = locnp(10)
      lwrk2  = locnp(12)

      lcmul  = locnp(16)
      lwrk3  = locnp(21)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)
      lcJac  = locnp(27)
      lgrad  = locnp(28)

      lf     = locnl( 1)
      lfJac  = locnl( 2)
      lfJdx  = locnl( 3)
      lyf    = locnl( 4) 

      ldcJ   = max ( ncnln, 1 )
      ldfJ   = max ( m    , 1 )

*     Allocate addresses not allocated in nploc or nlloc.

      lAx    = lwtotl + 1
      lwtotl = lAx    + nclin - 1
      lAx    = min( lAx, lwtotl )

*     Check input parameters and storage limits.

      call cmchk ( nerror, msgNP, lcrash, .false.,
     $             leniw, lenw, litotl, lwtotl,
     $             n, nclin, ncnln,
     $             istate, iw, named, names,
     $             bigbnd, bl, bu, clamda, x )

      if (nerror .gt. 0) then
         inform = 9
         go to 800
      end if

      tolrnk = one / hcndbd
      Rcndbd = sqrt( hcndbd )

      if (tolfea .gt. zero)
     $   call dload ( nplin, tolfea, w(lfeatl), 1 )

      if (ncnln .gt. 0  .and.  ctol .gt. zero)
     $   call dload ( ncnln, ctol, w(lfeatl+nplin), 1 )

      if (lfdset .eq. 0) then
         fdchk = sqrt( epsrf )
      else if (lfdset .eq. 1) then
         fdchk = fdint
      else
         fdchk = w(lhfrwd)
      end if

      nfun   = 0
      ngrad  = 0
      nstate = 1

      xnorm  = dnrm2 ( n, x, 1 )
      xnorm0 = xnorm

*     ------------------------------------------------------------------
*     If required,  compute the problem functions.
*     If the constraints are nonlinear,  the first call of funcon
*     sets up any constant elements in the Jacobian matrix.  A copy of
*     the Jacobian (with constant elements set) is placed in  cJacu.
*     ------------------------------------------------------------------
      if (lverfy .ge. 10) then
         lvrfyc = lverfy - 10

         call nlchkd( info, msgNP, nstate, lvlder, nfun, ngrad,
     $                ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                funcon, funobj, iw(lneedc),
     $                bigbnd, epsrf, cdint, fdint,
     $                fdchk, fdnorm, xnorm,
     $                bl, bu, c, w(lwrk3), w(lcJac), cJacu, w(lcJdx),
     $                f, w(lf), w(lfJac), fJacu, w(lfJdx),
     $                w(ldx), w(lhfrwd), w(lhctrl),
     $                x, w(lwrk1), w(lwrk2), w(lyf) )

         if (info .ne. 0) then
            if (info .gt. 0) inform = 7
            if (info .lt. 0) inform = info
            go to 800
         end if
         nstate = 0
      end if

      if (lcrash .lt. 2) then
*        ===============================================================
*        Cold or warm start.  Use  lscore  to obtain a point that
*        satisfies the linear constraints.
*        ===============================================================
         if (nclin .gt. 0) then
            ianrmj = lanorm
            do 110, j = 1, nclin
               w(ianrmj) = dnrm2 ( n, A(j,1), ldA )
               ianrmj    = ianrmj + 1
  110       continue
            call dcond ( nclin, w(lanorm), 1, Asize, amin )
         end if

         call dcond ( nplin, w(lfeatl), 1, feamax, feamin )
         call dcopy ( nplin, w(lfeatl), 1, w(lwtinf), 1 )
         call dscal ( nplin, (one/feamin), w(lwtinf), 1 )

*        ===============================================================
*        The input values of x and (optionally)  istate are used by
*        lscrsh  to define an initial working set.
*        ===============================================================
         vertex = .false.
         call lscrsh( cold, vertex,
     $                nclin, nplin, nactiv, nartif,
     $                nfree, n, ldA,
     $                istate, iw(lkactv),
     $                bigbnd, tolact,
     $                A, w(lAx), bl, bu, x, w(lwrk1), w(lwrk2) )

         unitQ  = .true.
         nres   = 0
         ngq    = 0
C-->     condmx = max( one/epspt5, hundrd )
C-->     condmx = max( one/epspt3, hundrd )
         condmx = max( one/epspt5, hundrd )

         ikx    = lkx
         do 120,  i = 1, n
            iw(ikx) = i
            ikx     = ikx + 1
  120    continue

         if (cold) then
            nrank  = 0
         else
            nrank  = nlnx
            call dload ( nlnx, zero, w(lres0), 1 )
         end if

*        ---------------------------------------------------------------
*        Re-order kx so that the free variables come first.
*        If a warm start is required, nrank will be nonzero and the
*        factor R will be updated.
*        ---------------------------------------------------------------
         call lsbnds( unitQ,
     $                inform, nZ, nfree, nrank, nres, ngq,
     $                n, ldQ, ldA, ldR, ldT,
     $                istate, iw(lkx), condmx,
     $                A, R, w(lT), w(lres0), w(lgq), w(lQ),
     $                w(lwrk1), w(lwrk2), w(lrlam) )

*        ---------------------------------------------------------------
*        Factorize the initial working set.
*        ---------------------------------------------------------------
         if (nactiv .gt. 0) then
            nact1  = nactiv
            nactiv = 0

            call lsadds( unitQ, vertex,
     $                   inform, 1, nact1, nactiv, nartif, nZ, nfree,
     $                   nrank, nrejtd, nres, ngq,
     $                   n, ldQ, ldA, ldR, ldT,
     $                   istate, iw(lkactv), iw(lkx), condmx,
     $                   A, R, w(lT), w(lres0), w(lgq), w(lQ),
     $                   w(lwrk1), w(lwrk2), w(lrlam) )
         end if

         ssq1 = zero

         linobj = .false.
         call lssetx( linobj, rowerr, unitQ,
     $                nclin, nactiv, nfree, nrank, nZ,
     $                n, nplin, ldQ, ldA, ldR, ldT,
     $                istate, iw(lkactv), iw(lkx),
     $                jmax, errmax, ctx, xnorm,
     $                A, w(lAx), bl, bu, w(lgq), w(lres), w(lres0),
     $                w(lfeatl), R, w(lT), x, w(lQ), w(lwrk1),w(lwrk2) )

*        ---------------------------------------------------------------
*        Call  lscore  to find a feasible  x.
*        ---------------------------------------------------------------
*        Use  work2  as the multiplier vector.

         jinf   = 0
         lclam  = lwrk2

         call lscore( 'FP problem', named, names, linobj, unitQ,
     $                nLPerr, itns, jinf, nclin, nplin,
     $                nactiv, nfree, nrank, nZ, nZ1,
     $                n, ldA, ldR,
     $                istate, iw(lkactv), iw(lkx),
     $                ctx, obj, ssq1, suminf, numinf, xnorm,
     $                bl, bu, A, w(lclam), w(lAx),
     $                w(lfeatl), R, x, w )

         if (nLPerr .gt. 0) then
            inform = 2
            go to 800
         else if (msgQP .gt. 0) then 
            if (iPrint .gt. 0) write(iPrint, 7000)
            if (iSumm  .gt. 0) write(iSumm , 7000)
         end if
      end if

*     ==================================================================
*     Check the gradients at this first feasible x.
*     ==================================================================
      if (lverfy .ge. 10  .and. xnorm .eq. xnorm0) then
*        Relax, we already have everything at this x.
      else
         lvrfyc = lverfy
         if (lverfy .ge. 10) lvrfyc = -1

         call nlchkd( info, msgNP, nstate, lvlder, nfun, ngrad,
     $                ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                funcon, funobj, iw(lneedc),
     $                bigbnd, epsrf, cdint, fdint,
     $                fdchk, fdnorm, xnorm,
     $                bl, bu, c, w(lwrk3), w(lcJac), cJacu, w(lcJdx),
     $                f, w(lf), w(lfJac), fJacu, w(lfJdx),
     $                w(ldx), w(lhfrwd), w(lhctrl),
     $                x, w(lwrk1), w(lwrk2), w(lyf) )

         if (info .ne. 0) then
            if (info .gt. 0) inform = 7
            if (info .lt. 0) inform = info
            go to 800
         end if
      end if

      mode = 2
      call nlobjf( mode, n, m, y, f, w(lfJac), ldfJ, 
     $             objf, w(lgq), w(lyf) )
      call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $             iw(lkx), w(lgq), w(lQ), w(lwrk1) )

      if ( cold ) then
*        ---------------------------------------------------------------
*        Cold start.
*        ---------------------------------------------------------------
         if (ltypeH .eq. 0) then

*           Initialize  R'R  as J'J.

            call nlJtJ ( n, ldQ, nfree, unitQ,
     $                   iw(lkx),
     $                   m, w(lfJac), ldfJ,
     $                   R, ldR, w(lQ),
     $                   w(lwrk1) )
            Rfrobn = f06qgf( 'Frobenius norm', 'Upper', n, n, R, ldR )
         else if (ltypeH .eq. 1) then

*           Initialize  R'R  as the identity.

            call f06qhf( 'Upper-triangular', n, n, zero, one, R, ldR )
            Rfrobn = rootn
         end if         
         if (ncnln .gt. 0) call dload ( ncnln, zero, w(lcmul), 1 )
      else
*        ---------------------------------------------------------------
*        Warm start.
*        Set the multipliers for the nonlinear constraints.
*        Check the condition of the initial factor R.
*        ---------------------------------------------------------------
         if (ncnln .gt. 0)
     $      call dcopy ( ncnln, clamda(nplin+1), 1, w(lcmul), 1 )
         Rfrobn = f06qgf( 'Frobenius norm', 'Upper', n, n, R, ldR )
      end if

*     ------------------------------------------------------------------
*     Check the condition of the initial factor R.
*     Compute the diagonal condition estimator of R.
*     ------------------------------------------------------------------
      call dcond ( n, R, ldR+1, dRmax, dRmin )
      cond   = ddiv  ( dRmax, dRmin, overfl )

      if (      cond   .gt. Rcndbd
     $    .or.  Rfrobn .gt. rootn*growth*dRmax) then
*        ---------------------------------------------------------------
*        Refactorize the Hessian and bound the condition estimator.
*        ---------------------------------------------------------------
         call nprset( unitQ,
     $                n, nfree, nZ, ldQ, ldR,
     $                iw(liperm), iw(lkx),
     $                w(lgq), R, w(lQ), w(lwrk1), w(lres0) )
      end if

*     ==================================================================
*     Solve the problem.
*     ==================================================================
      if (ncnln .eq. 0) then
*        ---------------------------------------------------------------
*        The problem has only linear constraints and bounds.
*        ---------------------------------------------------------------
         call nlcore( named, names, unitQ, inform, iter,
     $                m, n, nclin, ncnln, nactiv, nfree, nZ,
     $                ldcJ, ldcJu, ldfJ, ldfJu, ldA, ldR,
     $                nfun, ngrad, istate, iw(lkactv), iw(lkx),
     $                objf, fdnorm, xnorm, funcon, funobj,
     $                A, w(lAx), bl, bu, c, w(lcJac), cJacu, clamda,
     $                y, f, w(lfJac), fJacu, w(lfeatl), w(lgrad), R, x,
     $                iw, w, lenw )
      else
*        ---------------------------------------------------------------
*        The problem has some nonlinear constraints.
*        ---------------------------------------------------------------
         if (nclin .gt. 0)
     $      call f06qff( 'General', nclin, n, A, ldA, w(lAqp), ldAqp )

*        Try and add some nonlinear constraint indices to kactiv.

         call npcrsh( cold, n, nclin, ncnln,
     $                nctotl, nactiv, nfree, nZ,
     $                istate, iw(lkactv), bigbnd, tolact,
     $                bl, bu, c )

         call nlcore( named, names, unitQ, inform, iter,
     $                m, n, nclin, ncnln, nactiv, nfree, nZ,
     $                ldcJ, ldcJu, ldfJ, ldfJu, ldAqp, ldR,
     $                nfun, ngrad, istate, iw(lkactv), iw(lkx),
     $                objf, fdnorm, xnorm, funcon, funobj,
     $                w(lAqp), w(lAx), bl, bu, c, w(lcJac),cJacu,clamda,
     $                y, f, w(lfJac), fJacu, w(lfeatl), w(lgrad), R, x,
     $                iw, w, lenw )
      end if

*     ------------------------------------------------------------------
*     If required, overwrite R with the factor of the Hessian.
*     ------------------------------------------------------------------
      if (lformH .gt. 0) then
         call lsfrmH( 'Hessian', unitQ, 
     $                nfree, n, n, ldQ, ldR,
     $                iw(lkx), R, w(lQ), w(lwrk1), w(lwrk2) )
      end if

*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 if (msgNP .gt.   0) then
         if (iPrint .gt. 0) then
            if (inform .lt.   0) write(iPrint, 3000)
            if (inform .eq.   0) write(iPrint, 4000)
            if (inform .eq.   1) write(iPrint, 4100)
            if (inform .eq.   2) write(iPrint, 4200)
            if (inform .eq.   3) write(iPrint, 4300)
            if (inform .eq.   4) write(iPrint, 4400)
            if (inform .eq.   5) write(iPrint, 4500)
            if (inform .eq.   6) write(iPrint, 4600)
            if (inform .eq.   7) write(iPrint, 4700)
            if (inform .eq.   9) write(iPrint, 4900) nerror
         end if

         if (iSumm  .gt. 0) then
            if (inform .lt.   0) write(iSumm , 3000)
            if (inform .eq.   0) write(iSumm , 4000)
            if (inform .eq.   1) write(iSumm , 4100)
            if (inform .eq.   2) write(iSumm , 4200)
            if (inform .eq.   3) write(iSumm , 4300)
            if (inform .eq.   4) write(iSumm , 4400)
            if (inform .eq.   5) write(iSumm , 4500)
            if (inform .eq.   6) write(iSumm , 4600)
            if (inform .eq.   7) write(iSumm , 4700)
            if (inform .eq.   9) write(iSumm , 4900) nerror
         end if

         if (inform .ge. 0  .and.  inform .lt. 7) then
            if (nLPerr .eq. 0) then
               if (iPrint .gt. 0) write(iPrint, 5000) objf
               if (iSumm  .gt. 0) write(iSumm , 5000) objf
            else
               if (nLPerr .eq. 3) then
                  if (iPrint .gt. 0) write(iPrint, 5010) suminf
                  if (iSumm  .gt. 0) write(iSumm , 5010) suminf
               else
                  if (iPrint .gt. 0) write(iPrint, 5020) suminf
                  if (iSumm  .gt. 0) write(iSumm , 5020) suminf
               end if
            end if
         end if
      end if

*     Recover the optional parameters set by the user.

      call icopy ( mxparm, ipsvls, 1, iprmls, 1 )
      call dcopy ( mxparm, rpsvls, 1, rprmls, 1 )
      call icopy ( mxparm, ipsvnp, 1, iprmnp, 1 )
      call dcopy ( mxparm, rpsvnp, 1, rprmnp, 1 )
      call icopy ( mxparm, ipsvnl, 1, iprmnl, 1 )
      call dcopy ( mxparm, rpsvnl, 1, rprmnl, 1 )

      if (inform .lt. 9) then
         if (ncnln .gt. 0) then
            call f06qff( 'General', ncnln, n, w(lcJac), ldcJ, 
     $                   cJacu, ldcJu)
         end if
         call f06qff( 'General', m    , n, w(lfJac), ldfJ, 
     $                fJacu, ldfJu)
      end if

      return

 3000 format(/ ' Exit NLSSOL - User requested termination.'          )
 4000 format(/ ' Exit NLSSOL - Optimal solution found.'              )
 4100 format(/ ' Exit NLSSOL - Optimal solution found, ',
     $         ' but the requested accuracy could not be achieved.'  )
 4200 format(/ ' Exit NLSSOL - No feasible point for the linear',
     $         ' constraints.')
 4300 format(/ ' Exit NLSSOL - No feasible point for the nonlinear',
     $         ' constraints.')
 4400 format(/ ' Exit NLSSOL - Too many major iterations.             ')
 4500 format(/ ' Exit NLSSOL - Problem is unbounded (or badly scaled).')
 4600 format(/ ' Exit NLSSOL - Current point cannot be improved upon. ')
 4700 format(/ ' Exit NLSSOL - Large errors found in the derivatives. ')

 4900 format(/ ' Exit NLSSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.')
 5000 format(/ ' Final nonlinear objective value =', G16.7 )
 5010 format(/ ' Minimum sum of infeasibilities =',  G16.7 )
 5020 format(/ ' Final sum of infeasibilities =',    G16.7 )

 7000 format(/ ' The linear constraints are feasible.')

*     end of  nlssol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlchkd( inform, msgNP, nstate, lvlder, nfun, ngrad,
     $                   ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                   funcon, funobj, needc,
     $                   bigbnd, epsrf, cdint, fdint,
     $                   fdchk, fdnorm, xnorm,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   f, f1, fJac, fJacu, fJdx, dx, hforwd, hcntrl,
     $                   x, wrk1, wrk2, wrk4 )

      implicit           double precision(a-h,o-z)
      integer            needc(*)
      double precision   c(*), c1(*), cJac(ldcJ,*), cJacu(ldcJu,*),
     $                   cJdx(*)
      double precision   f(m), f1(m), fJac(ldfJ,*), fJacu(ldfJu,*),
     $                   fJdx(m)
      double precision   bl(n), bu(n), dx(n), x(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   wrk1(n), wrk2(*), wrk4(m)
      external           funcon, funobj

*     ==================================================================
*     nlchkd  performs the following...
*     (1)  Computes the objective and constraint values f and c.
*     (2)  Evaluates the user-provided gradients in cJacu and fJacu.
*     (3)  Counts the missing gradients.
*     (4)  Loads the known gradients into cJac and fJac.
*     (5)  Checks that the known gradients are programmed correctly.
*     (6)  Computes the missing gradient elements.
*
*     Original version based on npchkd, written 4-September-1985.
*     This version of nlchkd dated  18-May-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)
                                    
      logical            centrl, needfd
      parameter         (rdummy =-11111.0d+0)

      infocj = 0
      infofj = 0
      nfdiff = 0
      ncdiff = 0
      ncset  = n*ncnln
      nfset  = n*m

      if (ncnln .gt. 0) then
*        ===============================================================
*        Compute the constraints and Jacobian matrix.
*        ===============================================================
*        If some derivatives are missing, load the Jacobian with dummy
*        values.  Any elements left unaltered after the call to funcon
*        must be estimated.  A record of the missing Jacobian elements
*        is stored in  cJacu.

         needfd = lvlder .eq. 0  .or.  lvlder .eq. 1

         if (needfd)
     $      call f06qhf( 'General', ncnln, n, rdummy, rdummy,
     $                   cJacu, ldcJu )

         call iload ( ncnln, 1, needc, 1 )

         mode   = 2
         call funcon( mode, ncnln, n, ldcJu,
     $                needc, x, c, cJacu, nstate )
         if (mode .lt. 0) go to 9999

         call f06qff( 'General', ncnln, n, cJacu, ldcJu, cJac, ldcJ )

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
     $     call f06qhf( 'General', m, n, rdummy, rdummy,
     $                   fJacu, ldfJu )

      mode   = 2
      call funobj( mode, m, n, ldfJu, x, f, fJacu, nstate )
      if (mode .lt. 0) go to 9999

      call f06qff( 'general', m, n, fJacu, ldfJu, fJac, ldfJ )

      if (needfd) then
         do 420, j = 1, n
            do 410, i = 1, m
               if (fJacu(i,j) .eq. rdummy) nfdiff = nfdiff + 1
  410       continue
  420    continue

         nfset = nfset - nfdiff
         if (nstate .eq. 1) then
            if (nfdiff .eq. 0) then
               if (lvlder .eq. 0) lvlder = 1
               if (lvlder .eq. 2) lvlder = 3
               if (msgNP .gt. 0  .and.  iPrint .gt. 0)
     $            write(iPrint, 2000) lvlder
            else
               if (msgNP .gt. 0  .and.  iPrint .gt. 0)
     $            write(iPrint, 2100) nfset, n*m, nfdiff
            end if
         end if
      end if

      nfun  = nfun  + 1
      ngrad = ngrad + 1

*     ==================================================================
*     Check whatever gradient elements have been provided.
*     ==================================================================
      if (lvrfyc .ge. 0) then
         if (ncset .gt. 0) then
            call chcJac( mode, lvlder, msgNP,
     $                   ncset, n, ncnln, ldcJ, ldcJu,
     $                   bigbnd, epsrf, epspt3, fdchk, xnorm,
     $                   funcon, needc,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   dx, wrk2, x, wrk1 )
            if (mode .lt. 0) go to 9999
            infocj = mode
         end if

         if (nfset .gt. 0) then
            call chfJac( mode, lvlder, msgNP,
     $                   nfset, m, n, ldfJ, ldfJu,
     $                   bigbnd, epsrf, epspt3, fdchk, xnorm,
     $                   funobj,
     $                   bl, bu, f, f1, fJac, fJacu, fJdx,
     $                   dx, wrk4, x, wrk1 )
            if (mode .lt. 0) go to 9999
            infofj = mode
         end if
      end if

      needfd = ncdiff .gt. 0  .or.  nfdiff .gt. 0
      if (needfd) then
*        ===============================================================
*        Compute the missing gradient elements.
*        ===============================================================
         call chfdls( mode, msgNP, lvlder,
     $                m, n, ncnln, ldcJ, ldcJu, ldfJ, ldfJu,
     $                bigbnd, epsrf, fdnorm,
     $                funcon, funobj, needc,
     $                bl, bu, c, c1, cJdx, cJac, cJacu,
     $                f, f1, fJdx, fJac, fJacu, hforwd, hcntrl,
     $                x, dx ) 
         if (mode .lt. 0) go to 9999

         if (lfdset .gt. 0) then
            centrl = .false.
            call nlfd  ( centrl, mode,
     $                   ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, cJdx, cJac, cJacu,
     $                   f, f1, fJdx, fJac, fJacu,
     $                   hforwd, hcntrl, x )

            if (mode .lt. 0) go to 9999
         end if
      end if

      inform = infocj + infofj

      return

*     The user requested termination.

 9999 inform = mode
      return

 1000 format(//' All constraint Jacobian elements have been set.  ',
     $         ' Derivative level increased to ', i4 )
 1100 format(//' The user sets ', i6, '   out of', i6,
     $         '   constraint Jacobian elements.'
     $       / ' Each iteration, ', i6, '   constraint Jacobian',
     $         ' elements will be estimated numerically.' )
 2000 format(//' All objective Jacobian elements have been set.  ',
     $         ' Derivative level increased to ', i4 )
 2100 format(//' The user sets ', i6, '   out of', i6,
     $         '   objective Jacobian elements.'
     $       / ' Each iteration, ', i6, '   objective Jacobian',
     $         ' elements will be estimated numerically.' )

*     end of  nlchkd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlcore( named, names, unitQ, inform, majIts,
     $                   m, n, nclin, ncnln, nactiv, nfree, nZ,
     $                   ldcJ, ldcJu, ldfJ, ldfJu, ldAqp, ldR,
     $                   nfun, ngrad, istate, kactiv, kx,
     $                   objf, fdnorm, xnorm, funcon, funobj,
     $                   Aqp, Ax, bl, bu, c, cJac, cJacu, clamda,
     $                   y, f, fJac, fJacu, featol, grad, R, x,
     $                   iw, w, lenw )

      implicit           double precision(a-h,o-z)

      logical            named
      integer            istate(*), kactiv(n), kx(n)
      integer            iw(*)
      double precision   Aqp(ldAqp,*), Ax(*),
     $                   c(*), cJac(ldcJ,*), cJacu(ldcJu,*)
      double precision   y(m), f(m), fJac(ldfJ,*), fJacu(ldfJu,*),
     $                   grad(n)
      double precision   R(ldR,*), x(n)
      double precision   bl(n+nclin+ncnln), bu(n+nclin+ncnln),
     $                   clamda(n+nclin+ncnln), featol(n+nclin+ncnln)
      double precision   w(lenw)
      external           funobj, funcon

      character*16       names(*)

*     ==================================================================
*     nlcore  is the core routine for  nlssol,  a sequential quadratic
*     programming (SQP) method for minimizing a sum of squares subject
*     to linear and nonlinear constraints.
*
*     Original version based on npcore,  11-May-1988.
*     This version of nlcore dated  12-Jul-94
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol3cm/ lennam, ldT   , ncolT , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize , dTmax , dTmin
      common    /sol6cm/ Rcndbd, Rfrobn, dRmax, dRmin

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      logical            incrun
      common    /sol6np/ Penmax, Pennrm, Pendmp, scale, incrun

      parameter         (lennl = 20)
      common    /sol1nl/ locnl(lennl)

*     +Include lsparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   itmax1, itmax2, lcrash, lformH, lprob , msgLS ,
     $                   nn    , nnclin, nprob , ipadls(21)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolOpt, tolrnk, rpadls(22)

      equivalence       (iprmls(1), itmax1 ), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, lvlder, 
     $                   lverfy, msgNP , nlnf  , nlnj  , nlnx  , nncnln,
     $                   nsave , nload , ksave , ipadnp(15)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , Hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), itmxnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     +Include nlparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnl(mxparm), ipsvnl
      double precision   rprmnl(mxparm), rpsvnl

      common    /nlpar1/ ipsvnl(mxparm),
     $                   ltypeH, nreset, ipadnl(28)

      common    /nlpar2/ rpsvnl(mxparm), rpadnl(30)

      equivalence       (iprmnl(1), ltypeH), (rprmnl(1), rpadnl(1))

      save      /nlpar1/, /nlpar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (itmxnp, nmajor), (itmax2, nminor), (msgLS , msgQP )

      logical            goodgq, newgq
      logical            centrl, convrg, convpt, done  , error , feasqp
      logical            infeas, needfd, optiml, overfl, reset , unitQ
      logical            KTcond(2)
      
      character*5        MjrMsg
      parameter         (zero   = 0.0d+0, one = 1.0d+0) 
      parameter         (growth = 1.0d+2              ) 

*     Specify machine-dependent parameters.

      flmax  = wmach(7)
      rtmax  = wmach(8)

      lanorm = locls( 2)
      lRpq   = locls( 5)
      lqrwrk = locls( 6)
      lHpq   = locls( 8)
      lgq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk1  = locls(14)
      lqptol = locls(15)

      liperm = locnp( 2)
      lAqp   = locnp( 3)
      lAdx   = locnp( 4)
      lbl    = locnp( 5)
      lbu    = locnp( 6)
      ldx    = locnp( 7)
      lgq1   = locnp( 8)
      lx1    = locnp(11)
      lwrk2  = locnp(12)
      lcs1   = locnp(13)
      lcs2   = locnp(14)
      lc1mul = locnp(15)
      lcmul  = locnp(16)
      lcJdx1 = locnp(17)
      ldlam  = locnp(18)
      ldslk  = locnp(19)
      lPen   = locnp(20)
      lwrk3  = locnp(21)
      lslk1  = locnp(22)
      lslk   = locnp(23)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)

      lcJac1 = lAqp   + nclin
      lcJdx  = lAdx   + nclin
      lvioln = lwrk3

      lf     = locnl( 1)
      lfJdx  = locnl( 3)
      lyf    = locnl( 4) 

*     Initialize

      MjrMsg = '     '
      nQPinf = 0

      majIt0 = majIts
      nplin  = n     + nclin
      nctotl = nplin + ncnln
      ncqp   = nclin + ncnln
      nl     = min( nplin + 1, nctotl )

      ldcJ1 = max( ncqp , 1 )

      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or. (lvlder .eq. 1  .and.  ncnln .gt. 0)

      alfa   = zero
      alfdx  = zero
      rtftol = sqrt( ftol )
      rootn  = sqrt( dble(n) )

*     ------------------------------------------------------------------
*     Information from the feasibility phase will be used to generate a
*     hot start for the first QP subproblem.
*     ------------------------------------------------------------------
      call dcopy ( nctotl, featol, 1, w(lqptol), 1 )

      nstate = 0

      objalf = objf
      if (ncnln .gt. 0) then
         objalf = objalf - ddot  ( ncnln, w(lcmul), 1, c, 1 )

         incrun = .true.
         Pennrm = zero
         Pendmp = one
         scale  = one
         call dload ( ncnln, zero, w(lPen), 1 )
      end if

      newgq  = .false.

**    ==================================================================
*+    repeat                             (until converged or error exit)

  100    minIts = 0

**       ===============================================================
*+       repeat                         (Until a good gradient is found)

  110       centrl = lvldif .eq. 2

            if (newgq) then
               if (needfd) then
*                 ------------------------------------------------------
*                 Compute any missing gradient elements and the
*                 transformed gradient of the objective.
*                 ------------------------------------------------------
                  call nlfd  ( centrl, mode,
     $                         ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                         bigbnd, cdint, fdint, fdnorm,
     $                         funcon, funobj, iw(lneedc),
     $                         bl, bu, c, w(lwrk2), w(lwrk3),cJac,cJacu,
     $                         f, w(lf), w(lfJdx), fJac, fJacu,
     $                         w(lhfrwd), w(lhctrl), x )
                  inform = mode
                  if (mode .lt. 0) go to 800
               end if

               mode = 2
               call nlobjf( mode, n, m, y, f, fJac, ldfJ, 
     $                      dumobj, w(lgq), w(lyf) )
               call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $                      kx, w(lgq), w(lQ), w(lwrk1) )
               newgq  = .false.
            end if

*           ============================================================
*           (1) Solve an inequality quadratic program (IQP) for the
*               search direction and multiplier estimates.
*           (2) For each nonlinear inequality constraint,  compute
*               the slack variable for which the merit function is
*               minimized.
*           (3) Compute the search direction for the slack variables
*               and multipliers.
*
*           Note that the array violn is wrk3.
*           ============================================================
            call npiqp ( feasqp, unitQ, nQPerr, majIts, Mnr,
     $                   n, nclin, ncnln, ldcJ, ldAqp, ldR,
     $                   linact, nlnact, nactiv, nfree, nZ, numinf,
     $                   istate, kactiv, kx,
     $                   dxnorm, gdx, curvQP,
     $                   Aqp, w(lAdx), w(lanorm), Ax, bl, bu,
     $                   c, cJac, clamda, w(lcmul), w(lcs1),
     $                   w(ldlam), w(ldslk), w(ldx), w(lbl), w(lbu),
     $                   w(lqptol), R, w(lPen), w(lslk), w(lvioln), x,
     $                   w(lwtinf), iw, w )

            minIts = minIts + Mnr

            if (feasqp) then
               nQPinf = 0
            else
               nQPinf = nQPinf + 1
               MjrMsg(2:2) = 'infeasible subproblem'
            end if

*           ============================================================
*           Compute quantities needed for the convergence test.
*           ============================================================
*           Compute the norms of the reduced gradient and the
*           gradient with respect to the free variables.

            gznorm = zero
            if (nZ .gt. 0)
     $         gznorm = dnrm2 ( nZ   , w(lgq), 1 )
            gfnorm = gznorm
            if (nfree .gt. 0  .and.  nactiv .gt. 0)
     $         gfnorm = dnrm2 ( nfree, w(lgq), 1 )

*           If the forward-difference estimate of the transformed
*           gradient of the Lagrangian function is small,  switch to
*           central differences, recompute the derivatives and re-solve
*           the QP.

            goodgq = .true.
            if (needfd  .and.  .not. centrl) then
               glnorm = dnrm2 ( n, w(lHpq), 1 )
               if (ncnln .eq. 0) then
                  cnorm = zero
               else
                  cnorm = dnrm2 ( ncnln, c, 1 )
               end if

               gltest = (one + abs(objf) + abs(cnorm))*epsrf/fdnorm
               if (glnorm .le. gltest) then
                  goodgq      = .false.
                  MjrMsg(3:3) = 'central differences'
                  lvldif      = 2
                  newgq       = .true.
                  if (msgNP .ge. 5  .and.  iPrint .gt. 0) then
                     if (minIts .gt. 0) write(iPrint, 3000) minIts 
                  end if
               end if
            end if
*+       until     (goodgq)
         if (.not.  goodgq ) go to 110

*        ===============================================================
*        (1) Compute the number of constraints that are violated by more
*            than featol.
*        (2) Compute the 2-norm of the residuals of the constraints in
*            the QP working set.
*        ===============================================================
         call npfeas( n, nclin, ncnln, istate,
     $                bigbnd, cvnorm, errmax, jmax, nviol,
     $                Ax, bl, bu, c, featol, x, w(lwrk2) )

*        Define small quantities that reflect the magnitude of objf and
*        the norm of grad(free).

         objsiz = one + abs( objf )
         xsize  = one +  xnorm
         gtest  = max( objsiz, gfnorm )
         dinky  = rtftol * gtest

         if (nactiv .eq. 0) then
            condT = zero
         else if (nactiv .eq. 1) then
            condT = dTmin
         else
            condT = ddiv  ( dTmax, dTmin, overfl )
         end if

         call dcond ( n, R, ldR+1, dRmax, dRmin )

         condH = ddiv  ( dRmax, dRmin, overfl )
         if (condH .lt. rtmax) then
            condH = condH*condH
         else
            condH = flmax
         end if

         if (nZ .eq. 0) then
            condHz = one
         else if (nZ .eq. n) then
            condHz = condH
         else
            call dcond ( nZ, R, ldR+1, drzmax, drzmin )
            condHz = ddiv  ( drzmax, drzmin, overfl )
            if (condHz .lt. rtmax) then
               condHz = condHz*condHz
            else
               condHz = flmax
            end if
         end if

*        ---------------------------------------------------------------
*        Test for convergence.
*        The point test convpt checks for a K-T point at the initial
*        point or after a large change in x.
*        ---------------------------------------------------------------
         convpt    = dxnorm .le. epspt8*gtest  .and.  nviol  .eq. 0
     $                                         .and.  nQPerr .le. 1
         KTcond(1) = gznorm .lt. dinky
         KTcond(2) = nviol  .eq. 0
         optiml    = KTcond(1)  .and.  KTcond(2)

         convrg    = majIts .gt. 0  .and.  alfdx .le. rtftol*xsize

         infeas    =       convrg         .and.  .not. feasqp
     $               .or.  nQPinf .gt. 7

         done      = convpt  .or.  (convrg  .and. optiml)
     $                       .or.   infeas

         objalf = objf
         grdalf = gdx
         gL1    = gdx
         if (ncnln .gt. 0) then
            gL1 = gL1 - ddot( ncnln, w(lcJdx), 1, clamda(nl), 1 )

*           Compute the value and directional derivative of the
*           augmented Lagrangian merit function.
*           The penalty parameters may be increased or decreased.

            call npmrt ( feasqp, n, nclin, ncnln,
     $                   objalf, grdalf, curvQP,
     $                   istate,
     $                   w(lcJdx), w(lcmul), w(lcs1),
     $                   w(ldlam), w(lPen), w(lvioln),
     $                   w(lwrk1), w(lwrk2) )
         end if

*        ===============================================================
*        Print the details of this iteration.
*        ===============================================================
         call npprt ( KTcond, convrg, MjrMsg, msgNP, msgQP,
     $                ldR, ldT, n, nclin, ncnln,
     $                nctotl, nactiv, linact, nlnact, nZ, nfree,
     $                majIt0, majIts, minIts, istate, alfa, nfun,
     $                condHz, condT, objalf, objf, gznorm, cvnorm,
     $                Ax, c, R, w(lT), w(lvioln), x, w(lwrk1) )

         alfa  = zero
         error = majIts .ge. nmajor

         if (.not. (done  .or.  error)) then
            majIts = majIts + 1

*           Make copies of information needed for the BFGS update.

            call dcopy ( n, x     , 1, w(lx1) , 1 )
            call dcopy ( n, w(lgq), 1, w(lgq1), 1 )

            if (ncnln .gt. 0) then
               call dcopy ( ncnln, w(lcJdx), 1, w(lcJdx1), 1 )
               call dcopy ( ncnln, w(lcmul), 1, w(lc1mul), 1 )
               call dcopy ( ncnln, w(lslk) , 1, w(lslk1) , 1 )
            end if

*           ============================================================
*           Compute parameters for the line search.
*           ============================================================
*           alfmin is the smallest allowable step predicted by the QP
*           subproblem.

            alfmin = one
            if (.not. feasqp) alfmin = zero

*           ------------------------------------------------------------
*           alfmax is the largest feasible steplength subject to a user-
*           defined limit alflim on the change in X.
*           ------------------------------------------------------------
            if (ncnln .gt. 0  .and.  needfd) then
               alfmax = one
            else
               alfmax = ddiv  ( bigdx, dxnorm, overfl )
               call npalf ( info, n, nclin, ncnln,
     $                      alfa, alfmin, alfmax, bigbnd, dxnorm,
     $                      w(lanorm), w(lAdx), Ax, bl, bu,
     $                      w(ldslk), w(ldx), w(lslk), x )
               alfmax = alfa
               if (alfmax .lt. one + epspt3  .and.  feasqp)
     $            alfmax = one
            end if

*           ------------------------------------------------------------
*           alfbnd is a tentative upper bound on the steplength.  If the
*           merit function is decreasing at alfbnd and certain
*           conditions hold,  alfbnd will be increased in multiples of
*           two (subject to not being greater than alfmax).
*           ------------------------------------------------------------
            if (ncnln .eq. 0) then
               alfbnd = alfmax
            else
               alfbnd = min( one, alfmax )
            end if

*           ------------------------------------------------------------
*           alfsml trips the computation of central differences.  If a
*           trial steplength falls below alfsml, the linesearch is
*           terminated.
*           ------------------------------------------------------------
            alfsml = zero
            if (needfd  .and. .not. centrl) then
               alfsml = ddiv  ( fdnorm, dxnorm, overfl )
               alfsml = min   ( alfsml, alfmax )
            end if

*           ============================================================
*           Compute the steplength using safeguarded interpolation.
*           ============================================================
            alflim = ddiv ( (one+xnorm)*dxlim, dxnorm, overfl )
            alfa   = min  ( alflim, one )

            call nlsrch( needfd, nlserr, m, n, ncnln,
     $                   ldcJ, ldcJu, ldfJ, ldfJu, nfun, ngrad,
     $                   iw(lneedc), funcon, funobj,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, gL1, gL2,
     $                   objf, objalf, curvQP, xnorm,
     $                   c, cJac, cJacu, w(lcJdx),
     $                   w(lc1mul), w(lcmul), w(lcs1), w(lcs2),
     $                   w(ldx), w(ldlam), w(ldslk), y, f, fJac, fJacu,
     $                   grad, w(lgq), clamda(nl), w(lPen),
     $                   w(lslk1), w(lslk), w(lx1), x, w, lenw )

*           ------------------------------------------------------------
*           nlsrch  sets nlserr to the following values...
*
*           nlserr will be negative if the user set mode lt 0.
*
*           values of nlserr occurring with a nonzero value of alfa.
*           1 -- if the search was successful and alfa lt alfmax.
*           2 -- if the search was successful and alfa  = alfmax.
*           3 -- if the search ended after mfsrch iterations.
*
*           values of nlserr occurring with a zero value of alfa....
*           4 -- if alfmax was too small.
*           6 -- if no improved point could be found.
*           7 -- if the input value of gdx is non-negative.
*           ------------------------------------------------------------
            if (nlserr .lt. 0) then
               inform = nlserr
               go to 800
            end if

            if (alfa .gt. alflim) MjrMsg(4:4) = 'l'

            error  = nlserr .ge. 4
            if (error) then
*              ---------------------------------------------------------
*              The line search failed to find a better point.
*              If exact gradients or central differences are being used,
*              or the KT conditions are satisfied, stop.  Otherwise,
*              switch to central differences and solve the QP again.
*              ---------------------------------------------------------
               if (needfd  .and.  .not. centrl) then
                  if (.not. optiml) then
                     error       = .false.
                     MjrMsg(3:3) = 'central differences'
                     lvldif      = 2
                     newgq       = .true.
                     if (msgNP .ge. 5  .and.  iPrint .gt. 0) then
                        write(iPrint, 3000) minIts
                     end if
                  end if
               end if
            else
               if (needfd) then
*                 ======================================================
*                 Compute any missing gradients.
*                 ======================================================
                  mode  = 1
                  ngrad = ngrad + 1

                  if (ncnln .gt. 0) then
                     call iload ( ncnln, 1, iw(lneedc), 1 )

                     call funcon( mode, ncnln, n, ldcJu, iw(lneedc),
     $                            x, w(lwrk1), cJacu, nstate )
                     inform = mode
                     if (mode .lt. 0) go to 800

                     call f06qff( 'General', ncnln, n, cJacu, ldcJu,
     $                            cJac, ldcJ )
                  end if

                  call funobj( mode, m, n, ldfJu,
     $                         x, w(lf), fJacu, nstate )
                  inform = mode
                  if (mode .lt. 0) go to 800

                  call f06qff( 'General', m, n, fJacu, ldfJu,
     $                         fJac, ldfJ )

                  call nlfd  ( centrl, mode,
     $                         ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                         bigbnd, cdint, fdint, fdnorm,
     $                         funcon, funobj, iw(lneedc),
     $                         bl, bu, c, w(lwrk2), w(lwrk3),cJac,cJacu,
     $                         f, w(lf), w(lfJdx), fJac, fJacu,
     $                         w(lhfrwd), w(lhctrl), x )

                  inform = mode
                  if (mode .lt. 0) go to 800

                  mode = 2
                  call nlobjf( mode, n, m, y, f, fJac, ldfJ,
     $                         objf, grad, w(lyf) )
                  gdx = ddot( n, grad, 1, w(ldx), 1 )
                  gL2 = gdx
                  if (ncnln .gt. 0) then
                     call dgemv ( 'No trans', ncnln, n, one, cJac, ldcJ,
     $                            w(ldx), 1, zero, w(lcJdx), 1 )
                     gL2 = gL2 -
     $                         ddot( ncnln, w(lcJdx), 1, clamda(nl), 1 )
                  end if
               end if

               call dcopy ( n, grad, 1, w(lgq), 1 )
               call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $                      kx, w(lgq), w(lQ), w(lwrk1) )

               xnorm  = dnrm2 ( n, x, 1 )

               if (ncnln .gt. 0  .and.  alfa .ge. one)
     $            call dcopy ( ncnln, clamda(nl), 1, w(lcmul), 1 )

               if (nclin .gt. 0)
     $            call daxpy ( nclin, alfa, w(lAdx), 1, Ax, 1 )
               alfdx = alfa * dxnorm

               reset = nlnact .eq. 0  .and.  mod(majIts,nreset) .eq. 0
               if (.not. reset) then
*                 ======================================================
*                 Update the factors of the approximate Hessian of the
*                 Lagrangian function.
*                 ======================================================
                  call npupdt( MjrMsg, unitQ,
     $                         n, ncnln, nfree, nZ,
     $                         ldcJ1, ldcJ, ldQ, ldR, kx,
     $                         alfa, gL1, gL2, curvQP,
     $                         w(lcJac1), cJac, w(lcJdx1), w(lcJdx),
     $                         w(lcs1), w(lcs2), w(lgq1), w(lgq),
     $                         w(lHpq), w(lRpq), clamda(nl), R,
     $                         w(lwrk3), w(lQ), w(lwrk2), w(lwrk1) )

                  call dcond ( n, R, ldR+1, dRmax, dRmin )
                  cond   = ddiv  ( dRmax, dRmin, overfl )
                  reset  =       cond   .gt. Rcndbd
     $                     .or.  Rfrobn .gt. rootn*growth*dRmax
               end if

               if (reset) then
*                 ------------------------------------------------------
*                 Reset R'R  as J'J.
*                 ------------------------------------------------------
                  call nlJtJ ( n, ldQ, nfree, unitQ, kx, m, fJac, ldfJ,
     $                         R, ldR, w(lQ), w(lwrk1) )
                  Rfrobn = f06qgf( 'Frobenius norm', 'Upper', 
     $                             n, n, R, ldR )
                  call dcond ( n, R, ldR+1, dRmax, dRmin )
                  cond   = ddiv  ( dRmax, dRmin, overfl )
                  reset  =       cond   .gt. Rcndbd
     $                     .or.  Rfrobn .gt. rootn*growth*dRmax
               end if
               
               if (reset) then
*                 ------------------------------------------------------
*                 Reset the range-space partition of Q'HQ.
*                 ------------------------------------------------------
                  MjrMsg(5:5) = 'refactorize Hessian'

                  call nprset( unitQ,
     $                         n, nfree, nZ, ldQ, ldR,
     $                         iw(liperm), kx,
     $                         w(lgq), R, w(lQ), w(lwrk1), w(lqrwrk) )
               end if
            end if
         end if

*+    until     (done  .or.  error)
      if (.not. (done  .or.  error) ) go to 100

*     ======================end of main loop============================

      if (done) then
         if (convrg  .and.  optiml) then
            inform = 0
         else if (convpt) then
            inform = 1
         else if (infeas) then
            inform = 3
         end if
      else if (error) then
         if (majIts .ge. nmajor) then
            inform = 4
         else if (optiml) then
            inform = 1
         else
            inform = 6
         end if
      end if

*     ------------------------------------------------------------------
*     Set  clamda.  Print the full solution.
*     ------------------------------------------------------------------
  800 if (msgNP .gt. 0  .and.  iPrint .gt. 0) then
         write(iPrint, 2100) inform, majIts, nfun, ngrad
      end if

      call cmwrp ( nfree, ldAqp,
     $             n, nclin, nctotl,
     $             nactiv, istate, kactiv, kx,
     $             Aqp, bl, bu, c, clamda, featol,
     $             w(lwrk1), w(lrlam), x )
      call cmprnt( msgNP, n, nclin, nctotl, bigbnd,
     $             named, names, istate,
     $             bl, bu, clamda, featol, w(lwrk1) )

      return

 2100 format(/ ' Exit  NP phase.  Inform = ', i2, ' MajIts = ', i5,
     $         '   nfun = ', i5, '   ngrad = ', i5 )
 3000 format(  ' Minor itn', i6, '.  Central-differences computed. ',
     $         ' QP re-solved.' )

*     end of nlcore
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nldflt( n, nclin, ncnln, title )

      implicit           double precision(a-h,o-z)

      character*(*)      title

*     ==================================================================
*     nldflt  loads the default values of parameters not set in the
*     options file.
*
*     Original Fortran 77 version written 10-September-1985.
*     This version of nldflt dated 19-May-95.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)

      logical            newOpt, listOp
      common    /sol7np/ newOpt, listOp, ncalls
      save      /sol7np/

*     +Include lsparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   itmax1, itmax2, lcrash, lformH, lprob , msgLS ,
     $                   nn    , nnclin, nprob , ipadls(21)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolOpt, tolrnk, rpadls(22)

      equivalence       (iprmls(1), itmax1 ), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, lvlder, 
     $                   lverfy, msgNP , nlnf  , nlnj  , nlnx  , nncnln,
     $                   nsave , nload , ksave , ipadnp(15)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , Hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), itmxnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     +Include nlparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnl(mxparm), ipsvnl
      double precision   rprmnl(mxparm), rpsvnl

      common    /nlpar1/ ipsvnl(mxparm),
     $                   ltypeH, nreset, ipadnl(28)

      common    /nlpar2/ rpsvnl(mxparm), rpadnl(30)

      equivalence       (iprmnl(1), ltypeH), (rprmnl(1), rpadnl(1))

      save      /nlpar1/, /nlpar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (itmxnp, nmajor), (itmax2, nminor), (msgLS , msgQP )

      parameter         (zero   =  0.0d+0, one    =   1.0d+0) 
      parameter         (point3 =  3.3d-1, point8 =   0.8d+0) 
      parameter         (point9 =  0.9d+0, two    =   2.0d+0) 
      parameter         (tenp6  =  1.0d+6, hundrd = 100.0d+0) 
      parameter         (rdummy = -11111.0d+0, idummy = -11111 ) 
      parameter         (gigant =  1.0d+20*.99999d+0) 
      parameter         (wrktol =  1.0d-2) 

      character*4        icrsh(0:2)
      character*3        cHess(0:1)
      data                icrsh(0),  icrsh(1),  icrsh(2)
     $                 / 'Cold'   , 'Warm'   , 'Hot '    /
      data                cHess(0),  cHess(1)
     $                 / ' no',      'yes'   /

      epsmch = wmach( 3)
      condbd = max ( one/(hundrd*epsmch*dble(n)), tenp6 )

      nplin  = n     + nclin
      nctotl = nplin + ncnln

*     Make a dummy call nlnkey to ensure that the defaults are set.

      call nlnkey()
      newopt = .true.

*     Save any optional parameters set by the user.  nldflt may have
*     to change some components of iprmls, rprmls, iprmnp, rprmnp, 
*     iprmnl and rprmnl to their default values.
                     
      call icopy ( mxparm, iprmls, 1, ipsvls, 1 )
      call dcopy ( mxparm, rprmls, 1, rpsvls, 1 )
      call icopy ( mxparm, iprmnp, 1, ipsvnp, 1 )
      call dcopy ( mxparm, rprmnp, 1, rpsvnp, 1 )
      call icopy ( mxparm, iprmnl, 1, ipsvnl, 1 )
      call dcopy ( mxparm, rprmnl, 1, rpsvnl, 1 )

      if (          iPrint .lt. 0     )   call mcout ( iPrint, iSumry )
      if (          iSumm  .lt. 0     )   call mcout ( iPrntr, iSumm  )
      if (          iSumm  .eq. iPrint)   iPrint  = 0

      if (          lcrash .lt. 0
     $    .or.      lcrash .gt. 2     )   lcrash  =  0
      if (          lvlder .lt. 0
     $    .or.      lvlder .gt. 3     )   lvlder  =  3
      if (          lformH .lt. 0
     $    .or.      lformH .gt. 1     )   lformH  =  0
      if (          nmajor .lt. 0     )   nmajor  = max(50, 3*nplin+
     $                                                     10*ncnln )
      if (          nminor .lt. 1     )   nminor  = max(50, 3*nctotl)
      if (          msgNP  .eq. idummy)   msgNP   = 10
      if (          msgQP  .eq. idummy)   msgQP   =  0
                                          nlnf    =  n
                                          nlnj    =  n
                                          nlnx    =  n
      if (          jvrfy2 .le. 0
     $    .or.      jvrfy2 .gt. n     )   jvrfy2  =  n
      if (          jvrfy1 .le. 0
     $    .or.      jvrfy1 .gt. jvrfy2)   jvrfy1  =  1
      if (          jvrfy4 .le. 0
     $    .or.      jvrfy4 .gt. n     )   jvrfy4  =  n
      if (          jvrfy3 .le. 0
     $    .or.      jvrfy3 .gt. jvrfy4)   jvrfy3  =  1
      if (         (lverfy .lt. -1
     $    .or.      lverfy .gt. 13) .or.
     $             (lverfy .ge.  4 
     $    .and.     lverfy .le.  9)   )   lverfy  =  0
      if (          tolact .lt. zero
     $    .or.      tolact .ge. one   )   tolact  =  wrktol
      if (          tolfea .lt. epsmch
     $    .or.      tolfea .ge. one   )   tolfea  =  epspt5
      if (          tolOpt .lt. epsmch
     $    .or.      tolOpt .ge. one   )   tolOpt  =  epspt8
      if (          epsrf  .lt. epsmch
     $    .or.      epsrf  .ge. one   )   epsrf   =  epspt9
                                          lfdset  =  0
      if (          fdint  .lt. zero  )   lfdset  =  2
      if (          fdint  .eq. rdummy)   lfdset  =  0
      if (          fdint  .ge. epsmch
     $    .and.     fdint  .lt. one   )   lfdset  =  1
      if (          lfdset .eq. 1
     $    .and.    (cdint  .lt. epsmch
     $    .or.      cdint  .ge. one  ))   cdint   = epsrf**point3
      if (          bigbnd .le. zero  )   bigbnd  = gigant
      if (          bigdx  .le. zero  )   bigdx   = max( gigant,bigbnd )
      if (          dxlim  .le. zero  )   dxlim   = two
      if (          eta    .lt. zero
     $    .or.      eta    .ge. one   )   eta     = point9
      if (          ftol   .lt. epsrf
     $    .or.      ftol   .ge. one   )   ftol    = epsrf**point8

      if (          hcndbd .lt. one   )   hcndbd  = condbd

                                          dctol   = epspt5
      if (          lvlder .lt. 2     )   dctol   = epspt3
      if (          ctol   .lt. epsmch
     $    .or.      ctol   .ge. one   )   ctol    = dctol
      if (          ltypeH .lt. 0     )   ltypeH  = 0
      if (          nreset .le. 0     )   nreset  = 2

      itmax1    = max( 50, 3*(n + nclin + ncnln) )
      jverfy(1) = jvrfy1
      jverfy(2) = jvrfy2
      jverfy(3) = jvrfy3
      jverfy(4) = jvrfy4

      if (msgNP .gt. 0) then
*        ----------------
*        Print the title.
*        ----------------
         lent = len( title )
         nspace = (81 - lent)/2 + 1
         
         if (iPrint .gt. 0) then
            write(iPrint, '(///// (80a1) )')
     $            (' ', j=1, nspace), (title(j:j), j=1,lent)
            write(iPrint, '(80a1 //)')
     $            (' ', j=1, nspace), ('='       , j=1,lent)
         end if

         if (iSumm  .gt. 0) then
            write(iSumm , '(///// (80a1) )')
     $            (' ', j=1, nspace), (title(j:j), j=1,lent)
            write(iSumm , '(80a1 //)')
     $            (' ', j=1, nspace), ('='       , j=1,lent)
         end if

         if (iPrint .gt. 0  .and.  lcrash .le. 1) then
            write(iPrint, 2000)
            write(iPrint, 2100) nclin , tolfea, icrsh(lcrash) ,
     $                          n     , bigbnd, tolact,
     $                          dxlim , bigdx , cHess(lformH)
            write(iPrint, 2200) ncnln , ftol  , epsrf ,
     $                          nlnj  , ctol  , epsmch,
     $                          nlnf  , eta   , iPrint,
     $                          lvlder, lverfy, iSumm
            write(iPrint, 2300) nmajor, msgNP,
     $                          nminor, msgQP

            if (lvlder .lt. 3) then
               if      (lfdset .eq. 0) then
                  write(iPrint, 2400)
               else if (lfdset .eq. 1) then
                  write(iPrint, 2401) fdint, cdint
               else if (lfdset .eq. 2) then
                  write(iPrint, 2402)
               end if
            end if

            if (ltypeH .eq. 0) then
               write(iPrint, 2500) nreset
            else
               write(iPrint, 2501) nreset
            end if
         end if
      end if

      return

 2000 format(
     $//' Parameters'
     $/ ' ----------' )
 2100 format(
     $/ ' Linear constraints.....',     i10,   2x,
     $  ' Linear feasibility.....', 1p, e10.2, 2x,
     $  1x, a4, ' start.............'
     $/ ' Variables..............',     i10,   2x,
     $  ' Infinite bound size....', 1p, e10.2, 2x,
     $  ' Crash tolerance........',     e10.2
     $/ ' Step limit.............', 1p, e10.2, 2x,
     $  ' Infinite step size.....',     e10.2, 2x,
     $  ' Hessian................', 7x, a3  )
 2200 format(
     $/ ' Nonlinear constraints..',     i10,   2x,
     $  ' Optimality tolerance...', 1p, e10.2, 2x,
     $  ' Function precision.....',     e10.2
     $/ ' Nonlinear Jacobian vars',     i10,   2x,
     $  ' Nonlinear feasibility..', 1p, e10.2, 2x,
     $  ' eps (machine precision)', 1p, e10.2
     $/ ' Nonlinear objectiv vars',     i10,   2x,
     $  ' Linesearch tolerance...', 1p, e10.2, 2x,
     $  ' Print file.............',     i10
     $/ ' Derivative level.......',     i10,   2x,
     $  ' Verify level...........',     i10,   2x,
     $  ' Summary file...........',     i10)
 2300 format(
     $/ ' Major iterations limit.',     i10,   2x,
     $  ' Major print level......',     i10
     $/ ' Minor iterations limit.',     i10,   2x,
     $  ' Minor print level......',     i10 )
 2400 format(/ ' Difference intervals to be computed.' )
 2401 format(/ ' Difference interval....', 1p, e10.2, 2x,
     $         ' Central diffce interval',     e10.2 )
 2402 format(/ ' User-supplied difference intervals.' )
 2500 format(/ ' J''J initial Hessian....',12x,
     $         ' Reset frequency........', i10 )
 2501 format(/ ' Unit initial Hessian...', 12x,
     $         ' Reset frequency........', i10 )

*     end of nldflt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlfd  ( centrl, inform,
     $                   ldcJ, ldcJu, ldfJ, ldfJu, m, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, c2, cJac, cJacu,
     $                   f, f1, f2, fJac, fJacu,
     $                   hforwd, hcntrl, x )

      implicit           double precision(a-h,o-z)
      logical            centrl
      integer            needc(*)

      double precision   bl(n), bu(n), c(*), c1(*), c2(*),
     $                   cJac(ldcJ,*), cJacu(ldcJu,*)
      double precision   f(m), f1(m), f2(m),
     $                   fJac(ldfJ,*), fJacu(ldfJu,*)
      double precision   hforwd(n), hcntrl(n), x(n)
      external           funcon, funobj

*     ==================================================================
*     nlfd   evaluates any missing gradients.
*
*     Original version based on npfd written 3-July-1986.
*     This version of nlfd   dated 12-Jul-94.
*     ==================================================================
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      parameter         (rdummy=-11111.0d+0)
      parameter         (zero  = 0.0d+0, half  = 0.5d+0, one   = 1.0d+0)
      parameter         (three = 3.0d+0, four  = 4.0d+0                )

      inform = 0

*     ==================================================================
*     Use pre-assigned difference intervals to approximate derivatives.
*     ==================================================================
*     Use either the same interval for each component (lfdset = 1),
*     or the intervals already in hforwd or hcntrl (lfdset = 0 or 2).

      nstate =   0
      mode   =   0

      biglow = - bigbnd
      bigupp =   bigbnd

      fdnorm =   zero

      do 400, j  = 1, n

         xj     = x(j)

         ncmiss = 0
         if (ncdiff .gt. 0) then
            do 310, i = 1, ncnln
               if (cJacu(i,j) .eq. rdummy) then
                  needc(i) = 1
                  ncmiss   = ncmiss + 1
               else
                  needc(i) = 0
               end if
  310       continue
         end if

         nfmiss = 0
         if (nfdiff .gt. 0) then
            do 320, i = 1, m
               if (fJacu(i,j) .eq. rdummy) nfmiss   = nfmiss + 1
  320       continue
         end if

         if (ncmiss .gt. 0  .or.  nfmiss .gt. 0) then
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

            if (ncmiss .gt. 0) then
               call funcon( mode, ncnln, n, ldcJu,
     $                      needc, x, c1, cJacu, nstate )
               if (mode .lt. 0) go to 999
            end if

            if (nfmiss .gt. 0) then
               call funobj( mode, m, n, ldfJu, x, f1, fJacu, nstate )
               if (mode .lt. 0) go to 999
            end if

            if (centrl) then
*              ---------------------------------------------------------
*              Central differences.
*              ---------------------------------------------------------
               x(j)  = xj + delta + delta

               if (ncmiss .gt. 0) then
                  call funcon( mode, ncnln, n, ldcJu,
     $                         needc, x, c2, cJacu, nstate )
                  if (mode .lt. 0) go to 999

                  do 330, i = 1, ncnln
                     if (needc(i) .eq. 1)
     $                  cJac(i,j) = (four*c1(i) - three*c(i) - c2(i))
     $                                  / (delta + delta)
  330             continue
               end if

               if (nfmiss .gt. 0) then
                  call funobj( mode, m, n, ldfJu, x, f2, fJacu, nstate )
                  if (mode .lt. 0) go to 999

                  do 340, i = 1, m
                     if (fJacu(i,j) .eq. rdummy)
     $                  fJac(i,j) = (four*f1(i) - three*f(i) - f2(i))
     $                                  / (delta + delta)
  340             continue
               end if
            else
*              ---------------------------------------------------------
*              Forward Differences.
*              ---------------------------------------------------------
               if (ncmiss .gt. 0) then
                  do 350, i = 1, ncnln
                     if (needc(i) .eq. 1)
     $                  cJac(i,j) = (c1(i) -  c(i))/  delta
  350             continue
               end if

               if (nfmiss .gt. 0) then
                  do 360, i = 1, m
                     if (fJacu(i,j) .eq. rdummy)
     $                  fJac(i,j) = (f1(i) -  f(i))/  delta
  360             continue
               end if
            end if
         end if
         x(j) = xj
  400 continue

      return

  999 inform = mode

*     end of nlfd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlfile( iOptns, inform )
      integer            iOptns, inform

*     ==================================================================
*     nlfile  reads the options file from unit  iOptns  and loads the
*     options into the relevant elements of  iprmnp  and  rprmnp.
*
*     If  iOptns .lt. 0  or  iOptns .gt. 99  then no file is read,
*     otherwise the file associated with unit  iOptns  is read.
*
*     output:
*
*         inform = 0  if a complete  options  file was found
*                     (starting with  begin  and ending with  end);
*                  1  if  iOptns .lt. 0  or  iOptns .gt. 99;
*                  2  if  begin  was found, but end-of-file
*                     occurred before  end  was found;
*                  3  if end-of-file occurred before  begin  or
*                     endrun  were found;
*                  4  if  endrun  was found before  begin.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      logical            newOpt, listOp
      common    /sol7np/ newOpt, listOp, ncalls
      save      /sol7np/
      external           nlkey
*     ------------------------------------------------------------------
*     Update ncalls, the number of calls of nloptn and nlfile since the
*     start of this problem.
*     On the very first call, the default parameters are set.

      call nlnkey()
      call opfile( iOptns, iPrint, iSumm, 
     $             listOp, newOpt, inform, nlkey )
      newOpt = .false.

*     end of nlfile
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlloc ( m, n, litotl, lwtotl)

      implicit           double precision(a-h,o-z)

*     ==================================================================
*     nlloc  allocates additional addresses for nlssol.
*
*     Original version   11-May-1988.
*     This version of  nlloc  dated 12-Jul-94.
*     ==================================================================
      parameter         (lennl = 20)
      common    /sol1nl/ locnl(lennl)

      miniw     = litotl + 1
      minw      = lwtotl + 1

      lf        = minw
      lfJac     = lf     + m
      lfJdx     = lfJac  + m*n
      lyf       = lfJdx  + m
      minw      = lyf    + m

      locnl( 1) = lf
      locnl( 2) = lfJac
      locnl( 3) = lfJdx
      locnl( 4) = lyf

      litotl    = miniw - 1
      lwtotl    = minw  - 1

*     end of  nlloc
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlJtJ ( n, ldQ, nfree, unitQ, kx,
     $                   m, fJac, ldfJ, R, ldR, Q, work )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kx(n)
      double precision   fJac(ldfJ,*)
      double precision   R(ldR,*), Q(ldQ,*)
      double precision   work(n)

*     ==================================================================
*     nlJtJ  loads the Cholesky factor of Q'HQ, where  H = J'J, the
*     approximate Hessian.
*
*     Original version written 9-May-1989.
*     This version of  nlJtJ  dated  14-Jan-1991.
*     ==================================================================
      external           dcopy, dgeqr, dgeapq, dgemv
      parameter         (zero   =0.0d+0, one    =1.0d+0)

      nfixed = n - nfree
      mr     = min(m,n)

*     ------------------------------------------------------------------
*     Form the QR factorization of fJac.
*     Note that dgeqr requires m .ge. n.
*     ------------------------------------------------------------------
      call dgeqr ( m, mr, fJac, ldfJ, work, info )

      if (m .lt. n) then
         info = 0
         call dgeapq( 'Transpose', 'Separate', m, min(m,n),
     $                fJac, ldfJ, work,
     $                n-m, fJac(1,m+1), ldfJ, work(m+1), info )
      end if

      if (nfixed .gt. 0) then
         do 200, i = 1, mr
            do 100, l = 1, nfixed
               j = kx(nfree+l)
               if (i .le. j) then
                  R(i,nfree+l) = fJac(i,j)
               else
                  R(i,nfree+l) = zero
               end if
  100       continue
  200    continue
      end if
         
      if (nfree .gt. 0) then
         do 400, i = 1, mr
            do 300, l = 1, nfree
               j = kx(l)
               if (i .le. j) then
                  R(i,l) = fJac(i,j)
               else
                  R(i,l) = zero
               end if    
  300       continue
  400    continue
      end if

      if ((nfree .gt. 0) .and. (.not. unitQ)) then
         do 500, i = 1, mr
            call dgemv ( 'Transpose', nfree, nfree,
     $                   one, Q, ldQ, R(i,1), ldR,
     $                   zero, work, 1 )
            call dcopy ( nfree, work, 1, R(i,1), ldR )
  500    continue
      end if

      call dgeqr ( mr, mr, R, ldR, work, info )

      if (mr .lt. n) then
         info = 0
         call dgeapq( 'Transpose', 'Separate', mr, mr, R, ldR, work,
     $                n-mr, R(1,mr+1), ldR, work(mr+1), info )
         call f06qhf( 'General', n-mr, n, zero, zero, R(mr+1,1), ldR )
      end if

*     end of nlJtJ
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlkey ( iPrint, iSumm, listOp, buffer, key )

      implicit           double precision (a-h,o-z)
      character*(*)      buffer
      logical            listOp

*     ==================================================================
*     nlkey   decodes the option contained in  buffer  in order to set
*     a parameter value in the relevant element of the parameter arrays.
*
*     Input:
*        iPrint   the print   file for error messages
*        iSumm    the summary file for error messages.
*     Output:
*        key    The first keyword contained in buffer.
*
*        nlkey  calls opnumb and the subprograms
*               lookup, scannrl tokens, upcase
*        (now called oplook, opscan, optokn, opuppr)
*        supplied by Informatics General, Inc., Palo Alto, California.
*
*     This version of nlkey  dated 14-Sep-95.
*     ==================================================================
*     +Include lsparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   itmax1, itmax2, lcrash, lformH, lprob , msgLS ,
     $                   nn    , nnclin, nprob , ipadls(21)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolOpt, tolrnk, rpadls(22)

      equivalence       (iprmls(1), itmax1 ), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, lvlder, 
     $                   lverfy, msgNP , nlnf  , nlnj  , nlnx  , nncnln,
     $                   nsave , nload , ksave , ipadnp(15)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , Hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), itmxnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     +Include nlparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnl(mxparm), ipsvnl
      double precision   rprmnl(mxparm), rpsvnl

      common    /nlpar1/ ipsvnl(mxparm),
     $                   ltypeH, nreset, ipadnl(28)

      common    /nlpar2/ rpsvnl(mxparm), rpadnl(30)

      equivalence       (iprmnl(1), ltypeH), (rprmnl(1), rpadnl(1))

      save      /nlpar1/, /nlpar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (itmxnp, nmajor), (itmax2, nminor), (msgLS , msgQP )

      external           opnumb
      logical            more  , number, opnumb, sorted

      parameter         (     maxkey = 44,  maxtie = 21,   maxtok = 10)
      character*16       keys(maxkey), ties(maxtie), token(maxtok)
      character*16       key, key2, key3, value

      parameter         (idummy = -11111,  rdummy = -11111.0d+0,
     $                   sorted = .true.,  zero   =  0.0d+0    )

      data   keys
     $ / 'BEGIN           ',
     $   'CENTRAL         ', 'COLD            ', 'CONDITION       ',
     $   'CONSTRAINTS     ', 'CRASH           ', 'DEFAULTS        ',
     $   'DERIVATIVE      ', 'DIFFERENCE      ', 'END             ',
     $   'FEASIBILITY     ', 'FUNCTION        ', 'HESSIAN         ',
     $   'HOT             ', 'INFINITE        ', 'IPRMLS          ',
     $   'ITERATIONS      ', 'ITERS:ITERATIONS', 'ITNS :ITERATIONS',
     $   'JTJ             ',
     $   'LINE            ', 'LINEAR          ', 'LINESEARCH:LINE ',
     $   'LIST            ',
     $   'LOWER           ', 'MAJOR           ', 'MINOR           ',
     $   'NOLIST          ', 'NONLINEAR       ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'RESET           ',
     $   'ROW             ', 'RPRMLS          ', 'START           ',
     $   'STEP            ', 'STOP            ', 'SUMMARY         ',
     $   'UNIT            ', 'UPPER           ', 'VARIABLES       ',
     $   'VERIFY          ', 'WARM            '/

      data   ties
     $ / 'BOUND           ', 'CONSTRAINTS     ',
     $   'FEASIBILITY     ', 'FILE            ', 'GRADIENTS       ',
     $   'INITIAL         ',
     $   'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'JACOBIAN        ', 'LEVEL           ',
     $   'NO              ',
     $   'NO.      :NUMBER',
     $   'NUMBER          ', 'OBJECTIVE       ', 'PRINT           ',
     $   'SEARCH          ', 'STEP            ', 'TOLERANCE       ',
     $   'VARIABLES       ', 'YES             '/
*-----------------------------------------------------------------------

*     Eliminate comments and empty lines.
*     A '*' appearing anywhere in buffer terminates the string.

      i      = index( buffer, '*' )
      if (i .eq. 0) then
         lenbuf = len( buffer )
      else
         lenbuf = i - 1
      end if
      if (lenbuf .le. 0) then
         key = '*'
         go to 900
      end if

*     ------------------------------------------------------------------
*     Extract up to maxtok tokens from the record.
*     ntoken returns how many were actually found.
*     key, key2, key3 are the first tokens if any, otherwise blank.
*     ------------------------------------------------------------------
      ntoken = maxtok
      call optokn( buffer(1:lenbuf), ntoken, token )
      key    = token(1)
      key2   = token(2)
      key3   = token(3)

*     Certain keywords require no action.

      if (key .eq. ' '     .or.  key .eq. 'BEGIN' ) go to 900
      if (key .eq. 'LIST'  .or.  key .eq. 'NOLIST') go to 900
      if (key .eq. 'END'                          ) go to 900

*     Most keywords will have an associated integer or real value,
*     so look for it no matter what the keyword.

      i      = 1
      number = .false.

   50 if (i .lt. ntoken  .and.  .not. number) then
         i      = i + 1
         value  = token(i)
         number = opnumb( value )
         go to 50
      end if

      if (number) then
         read (value, '(bn, e16.0)') rvalue
      else
         rvalue = zero
      end if

*     convert the keywords to their most fundamental form
*     (upper case, no abbreviations).
*     sorted says whether the dictionaries are in alphabetic order.
*     loci   says where the keywords are in the dictionaries.
*     loci = 0 signals that the keyword wasn't there.

      call oplook( maxkey, keys, sorted, key , loc1 )
      call oplook( maxtie, ties, sorted, key2, loc2 )

*     ------------------------------------------------------------------
*     Decide what to do about each keyword.
*     The second keyword (if any) might be needed to break ties.
*     Some seemingly redundant testing of more is used
*     to avoid compiler limits on the number of consecutive else ifs.
*     ------------------------------------------------------------------
      more   = .true.
      if (more) then
         more   = .false.
         if      (key .eq. 'CENTRAL     ') then
            cdint  = rvalue
         else if (key .eq. 'COLD        ') then
            lcrash = 0
         else if (key .eq. 'CONDITION   ') then
            hcndbd = rvalue
         else if (key .eq. 'CONSTRAINTS ') then
            nnclin = rvalue
         else if (key .eq. 'CRASH       ') then
            tolact = rvalue
         else if (key .eq. 'DEFAULTS    ') then
            call mcout ( iPrint, iSumm )
            listOp = .true.
            do 20,     i = 1, mxparm
               iprmls(i) = idummy
               rprmls(i) = rdummy
               iprmnp(i) = idummy
               rprmnp(i) = rdummy
               iprmnl(i) = idummy
               rprmnl(i) = rdummy
   20       continue
         else if (key .eq. 'DERIVATIVE  ') then
            lvlder = rvalue
         else if (key .eq. 'DIFFERENCE  ') then
            fdint  = rvalue
         else if (key .eq. 'FEASIBILITY ') then
            tolfea = rvalue
            ctol   = rvalue
         else if (key .eq. 'FUNCTION    ') then
            epsrf  = rvalue
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'HESSIAN     ') then
            lformH = 1
            if   (key2.eq. 'NO          ') lformH = 0
         else if (key .eq. 'HOT         ') then
            lcrash = 2
         else if (key .eq. 'INFINITE    ') then
              if (key2.eq. 'BOUND       ') bigbnd = rvalue * 0.99999d+0
              if (key2.eq. 'STEP        ') bigdx  = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'IPRMLS      ') then
*           allow things like  IPRMLS 21 = 100  to set IPRMLS(21) = 100
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
               read (key3, '(bn, i16)') iprmls(ivalue)
            else
               if (iPrint .gt. 0)          write(iPrint, 2400) ivalue
               if (iSumm  .gt. 0)          write(iSumm , 2400) ivalue
            end if
         else if (key .eq. 'ITERATIONS  ') then
            nmajor = rvalue
         else if (key .eq. 'JTJ         ') then
            ltypeH = 0
         else if (key .eq. 'LINEAR      ') then
            if (key2  .eq. 'CONSTRAINTS ') nnclin = rvalue
            if (key2  .eq. 'FEASIBILITY ') tolfea = rvalue
            if (key2  .eq. 'SEARCH      ') eta    = rvalue
            if (loc2 .eq.  0             ) then
               if (iPrint .gt. 0)          write(iPrint, 2320) key2
               if (iSumm  .gt. 0)          write(iSumm , 2320) key2
            end if
         else if (key .eq. 'LINESEARCH  ') then
            eta    = rvalue
         else if (key .eq. 'LOWER       ') then
            bndlow = rvalue
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'MAJOR       ') then
              if (key2.eq. 'ITERATIONS  ') nmajor = rvalue
              if (key2.eq. 'PRINT       ') msgNP  = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'MINOR       ') then
              if (key2.eq. 'ITERATIONS  ') nminor = rvalue
              if (key2.eq. 'PRINT       ') msgQP  = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'NONLINEAR   ') then
              if (key2.eq. 'CONSTRAINTS ') nncnln = rvalue
              if (key2.eq. 'FEASIBILITY ') ctol   = rvalue
              if (key2.eq. 'JACOBIAN    ') nlnj   = rvalue
              if (key2.eq. 'OBJECTIVE   ') nlnf   = rvalue
              if (key2.eq. 'VARIABLES   ') nlnx   = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'OPTIMALITY  ') then
            ftol   = rvalue
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'PRINT       ') then
              if (key2.eq. 'FILE        ') iPrint = rvalue
              if (key2.eq. 'LEVEL       ') msgNP  = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'PROBLEM     ') then
              if (key2.eq. 'NUMBER      ') nprob  = rvalue
         else if (key .eq. 'RESET       ') then
            nreset = rvalue
         else if (key .eq. 'ROW         ') then
              if (key2.eq. 'TOLERANCE   ') ctol   = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'RPRMLS      ') then
*           Allow things like  rprmls 21 = 2  to set rprmls(21) = 2.0
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
               read (key3, '(bn, e16.0)') rprmls(ivalue)
            else
               if (iPrint .gt. 0) write(iPrint, 2400) ivalue
               if (iSumm  .gt. 0) write(iSumm , 2400) ivalue
            end if
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if (key .eq. 'START       ') then
              if (key2.eq. 'CONSTRAINTS ') jvrfy3 = rvalue
              if (key2.eq. 'OBJECTIVE   ') jvrfy1 = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'STEP        ') then
            dxlim  = rvalue
         else if (key .eq. 'STOP        ') then
              if (key2.eq. 'CONSTRAINTS ') jvrfy4 = rvalue
              if (key2.eq. 'OBJECTIVE   ') jvrfy2 = rvalue
              if (loc2.eq.  0            ) then
                 if (iPrint .gt. 0)        write(iPrint, 2320) key2
                 if (iSumm  .gt. 0)        write(iSumm , 2320) key2
              end if
         else if (key .eq. 'SUMMARY     ') then
            iSumm  = rvalue
         else if (key .eq. 'UNIT        ') then
            ltypeH = 1
         else if (key .eq. 'UPPER       ') then
            bndupp = rvalue
         else if (key .eq. 'VARIABLES   ') then
            nn     = rvalue
         else if (key .eq. 'VERIFY      ') then
              if (key2.eq. 'OBJECTIVE   ') lverfy =  1
              if (key2.eq. 'CONSTRAINTS ') lverfy =  2
              if (key2.eq. 'NO          ') lverfy = -1
              if (key2.eq. 'YES         ') lverfy =  3
              if (key2.eq. 'GRADIENTS   ') lverfy =  3
              if (key2.eq. 'LEVEL       ') lverfy =  rvalue
              if (loc2.eq.  0            ) lverfy =  3
         else if (key .eq. 'WARM        ') then
            lcrash = 1
         else
            if (iPrint .gt. 0) write(iPrint, 2300) key
            if (iSumm  .gt. 0) write(iSumm , 2300) key
         end if
      end if

  900 return

 2300 format(' XXX  Keyword not recognized:         ', a)
 2320 format(' XXX  Second keyword not recognized:  ', a)
 2400 format(' XXX  The parm subscript is out of range:', i10)

*     end of nlkey
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlnkey( )

      implicit           double precision (a-h,o-z)

*     ==================================================================
*     nlnkey  counts  consecutive calls of nloptn or nlfile.
*
*     Original version written  11-Sep-95,
*     This version of  nlnkey  dated  14-Sep-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      logical            newOpt, listOp
      common    /sol7np/ newOpt, listOp, ncalls
      save      /sol7np/

*     +Include lsparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   itmax1, itmax2, lcrash, lformH, lprob , msgLS ,
     $                   nn    , nnclin, nprob , ipadls(21)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolOpt, tolrnk, rpadls(22)

      equivalence       (iprmls(1), itmax1 ), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4, lvlder, 
     $                   lverfy, msgNP , nlnf  , nlnj  , nlnx  , nncnln,
     $                   nsave , nload , ksave , ipadnp(15)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , Hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), itmxnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     +Include nlparm-Sep-95++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnl(mxparm), ipsvnl
      double precision   rprmnl(mxparm), rpsvnl

      common    /nlpar1/ ipsvnl(mxparm),
     $                   ltypeH, nreset, ipadnl(28)

      common    /nlpar2/ rpsvnl(mxparm), rpadnl(30)

      equivalence       (iprmnl(1), ltypeH), (rprmnl(1), rpadnl(1))

      save      /nlpar1/, /nlpar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (itmxnp, nmajor), (itmax2, nminor), (msgLS , msgQP )

      parameter         (rdummy = -11111.0d+0, idummy = -11111)

      logical             first
      save                first
      data                first /.true./

      if ( first ) then
         nCalls = 0
         first  = .false.
         newOpt = .true.
         listOp = .true.

         call mcout ( iPrint, iSumm )
         do 10,     i = 1, mxparm
            rprmls(i) = rdummy
            iprmls(i) = idummy
            rprmnp(i) = rdummy
            iprmnp(i) = idummy
            rprmnl(i) = rdummy
            iprmnl(i) = idummy
   10    continue
         first  = .false.
      end if

      if ( newOpt ) then
         nCalls = 1
      else
         nCalls = nCalls + 1
      end if
  
*     end of nlnkey
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlobjf( mode, n, m, y, f, fJac, ldfJ,
     $                   objf, grad, yf )

      implicit           double precision(a-h,o-z)
      double precision   fJac(ldfJ,*)
      double precision   f(m), y(m), yf(m)
      double precision   grad(n)

*     ==================================================================
*     nlobjf  loads the objective and gradient of the function
*                      (1/2)(y - f(x))'(y - f(x)).
*
*     Original version written 12-Jul-94.
*     This version of  nlobjf  dated 12-Jul-94.
*     ==================================================================
      parameter         (zero   =0.0d+0, half = 0.5d+0, one = 1.0d+0)

      call dcopy ( m,         y, 1, yf, 1 )
      call daxpy ( m, (-one), f, 1, yf, 1 ) 
      if (mode .ne. 1) then
         objf = half*(dnrm2( m, yf, 1 ))**2
      end if

      if (mode .ne. 0) then
         call dgemv ( 'Transpose', m, n, (-one), fJac, ldfJ,
     $                yf, 1, zero, grad, 1 )
      end if

*     end of nlobjf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nloptn( string )
      character*(*)      string

*     ==================================================================
*     nloptn  loads the option supplied in string into the relevant
*     element of iprmlc, rprmlc, iprmnp or rprmnp.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      logical            newOpt, listOp
      common    /sol7np/ newOpt, listOp, ncalls
      save      /sol7np/

      character*16       key
      character*72       buffer
*     ------------------------------------------------------------------
      buffer = string

*     If this is the first call of nlnkey, set newOpt and default values
*     of the optional parameters. The default is to list the options.
*     Increment ncalls, the number of calls of nloptn and nlfile for
*     this optimization.

      call nlnkey()

*     Call  nlkey  to decode the option and set the parameter value.
*     If required, print a heading at the start of a new run.
*     Note that the following call to nlkey may reset iPrint and iSumm.

      call nlkey ( iPrint, iSumm, listOp, buffer, key )
      if (key .eq.  'LIST'  ) listOp = .true.
      if (key .eq.  'NOLIST') listOp = .false.

      if ( listOp ) then 
         if ( newOpt ) then
            if (iPrint .gt. 0) then
               write ( iPrint, '(// a / a /)' )
     $                         ' Optional Parameters',
     $                         ' -------------------'
            end if
            newOpt = .false.
         end if
         if (iPrint .gt. 0) write ( iPrint, '( 6x, a )'    ) buffer
      end if

*     end of nloptn
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlopti( string, ivalue )

      implicit           double precision (a-h,o-z)
      character*(*)      string
      integer            ivalue

*     ==================================================================
*     nlopti decodes the option contained in  string // ivalue.
*
*     14 Sep 1995: first version.
*     ==================================================================
      character*16       key
      character*72       buff72

      write(key, '(i16)') ivalue
      lenbuf = len(string)
      buff72 = string
      buff72(lenbuf+1:lenbuf+16) = key
      call nloptn( buff72 )

*     end of nlopti
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nloptr( string, rvalue )

      implicit           double precision (a-h,o-z)
      character*(*)      string
      double precision   rvalue

*     ==================================================================
*     nloptr decodes the option contained in  string // rvalue.
*
*     14 Sep 1995: first version.
*     ==================================================================
      character*16       key
      character*72       buff72

      write(key, '(1p, e16.8)') rvalue
      lenbuf = len(string)
      buff72 = string
      buff72(lenbuf+1:lenbuf+16) = key
      call nloptn( buff72 )

*     end of nloptr
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nlsrch( needfd, inform, m, n, ncnln,
     $                   ldcJ, ldcJu, ldfJ, ldfJu, nfun, ngrad,
     $                   needc, funcon, funobj,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, gL1, gL,
     $                   objf, objalf, curvQP, xnorm,
     $                   c, cJac, cJacu, cJdx, cmul1, cmul, cs1, cs,
     $                   dx, dlam, dslk, y, f, fJac, fJacu, 
     $                   grad, gradu, qpmul, Pen,
     $                   slk1, slk, x1, x, w, lenw )

      implicit           double precision (a-h,o-z)
      logical            needfd
      integer            needc(*)
      double precision   dx(n), x1(n), x(n)
      double precision   c(*), cJac(ldcJ,*), cJacu(ldcJu,*), cJdx(*),
     $                   cmul1(*), cmul(*), cs1(*), cs(*)
      double precision   y(m), f(m), fJac(ldfJ,*), fJacu(ldfJu,*)
      double precision   grad(n), gradu(n)
      double precision   dlam(*), dslk(*), qpmul(*),
     $                   Pen(*), slk1(*), slk(*)
      double precision   w(lenw)
      external           funcon, funobj

*     ==================================================================
*     nlsrch finds the steplength alfa that gives sufficient decrease in
*     the augmented Lagrangian merit function.
*
*     On exit, if inform = 1, 2 or 3,  alfa will be a nonzero steplength
*     with an associated merit function value  objalf  which is lower
*     than that at the base point. If  inform = 4, 5, 6, 7 or 8  alfa 
*     is zero and  objalf  will be the merit value at the base point.
*
*     Original version written  27-May-1985.
*     Level 2 BLAS added 12-June-1986.
*     This version of nlsrch dated 12-Jul-94.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)
      logical            incrun
      common    /sol6np/ Penmax, Pennrm, Pendmp, scale, incrun

      parameter         (lennl = 20)
      common    /sol1nl/ locnl(lennl)

      logical            debug, done, first, imprvd
      parameter         (zero = 0.0d+0, half = 0.5d+0)
      parameter         (two  = 2.0d+0, one  = 1.0d+0) 
      parameter         (tolg = 1.0d-1, rmu  = 1.0d-4) 

      epsmch = wmach(3)

      lc     = locls(14)                               
      lwork  = locnp(12)
      lcJdx  = locnp(21)                      
      lf     = locnl( 1)
      lyf    = locnl( 4) 

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
*             if  m(tolabs) - m(0) .gt. epsaf,  the linesearch tries
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
      alfbst = zero
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
      gtry   = zero
      tobj   = zero

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
     $         call dcopy ( ncnln, w(lc), 1, c, 1 )

            call dcopy ( m, w(lf), 1, f, 1 )

            if (.not. needfd) then
               gdx   = tgdx
               gL    = tgL
               call dcopy ( n, gradu, 1, grad, 1 )

               if (ncnln .gt. 0) then
                  call dcopy ( ncnln, w(lcJdx), 1, cJdx, 1 )
                  call f06qff( 'General', ncnln, n, cJacu, ldcJu,
     $                         cJac, ldcJ )
               end if

               call f06qff( 'General', m, n, fJacu, ldfJu,
     $                      fJac, ldfJ )
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

*              Compute new estimates of the multipliers and slacks.
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
               call funcon( mode, ncnln, n, ldcJu,
     $                      needc, x, w(lc), cJacu, nstate )
               if (mode .lt. 0) go to 999

               call dcopy ( ncnln,         w(lc), 1, cs, 1 )
               call daxpy ( ncnln, (-one), slk  , 1, cs, 1 )

               call dcopy ( ncnln, cs , 1, w(lwork), 1 )
               call ddscl ( ncnln, Pen, 1, w(lwork), 1 )

               fterm  =      ddot( ncnln, cmul    , 1, cs, 1 ) -
     $                  half*ddot( ncnln, w(lwork), 1, cs, 1 )
            end if

*           ------------------------------------------------------------
*           Compute the value and (if required) the Jacobian matrix of
*           the objective function.
*           ------------------------------------------------------------
            call funobj( mode, m, n, ldfJu, x, w(lf), fJacu, nstate )
            if (mode .lt. 0) go to 999

            call nlobjf( mode, n, m, y, w(lf), fJacu, ldfJu,
     $                   tobj, gradu, w(lyf) )

            if (ncnln .gt. 0) then
               tobjM = tobj  - fterm
            else
               tobjM = tobj
            end if

            ftry  = tobjM - oldf - rmu*oldg*alfa 

            if (.not. needfd) then
*              ---------------------------------------------------------
*              Compute auxiliary gradient information.
*              ---------------------------------------------------------
               gtry  = ddot( n, gradu, 1, dx, 1 )
               tgdx  = gtry
               tgL   = gtry

               if (ncnln .gt. 0) then

*                 Compute the Jacobian times the search direction.

                  call dgemv ( 'No', ncnln, n, one, cJacu, ldcJu, dx, 1,
     $                         zero, w(lcJdx), 1 )

                  call dcopy ( ncnln,         w(lcJdx), 1, w(lwork), 1 )
                  call daxpy ( ncnln, (-one), dslk    , 1, w(lwork), 1 )

                  gtry = gtry - ddot( ncnln, cmul, 1, w(lwork), 1 )
                  if (alfa .le. one)
     $               gtry = gtry - ddot( ncnln, dlam, 1, cs     , 1 )

                  call ddscl ( ncnln, Pen , 1, w(lwork), 1 )
                  gtry = gtry  + ddot( ncnln, w(lwork), 1, cs   , 1 )
                  tgL  = tgdx  - ddot( ncnln, w(lcJdx), 1, qpmul, 1 )

*                 ------------------------------------------------------
*                 If alfbnd .le. alfa .lt. alfmax and the norm of the
*                 quasi-Newton update is bounded, set alfmax to be alfa.
*                 This will cause the line search to stop if the merit
*                 function is decreasing at the boundary.
*                 ------------------------------------------------------
                  if (alfbnd .le. alfa  .and.  alfa .lt. alfmax) then
                     csJdx  = ddot   ( ncnln, cs, 1, w(lcJdx), 1 )
                     curvL  = tgL  - gL1
                     curvc  = abs( csJdx - cs1Jdx )
                     PenBFS = max( curvQP*tolg - curvL, zero )
                     if (PenBFS .le. curvc*Penmax) then
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

*     The user wants to stop.

  999 inform = mode
      return

*     end of nlsrch
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     file  npsolsubs.f
*
*     npsol    npalf    npchkd   npcore   npcrsh   npdflt   npfd     
*     npfeas   npfile   npgetr   npiqp    npkey    nploc    npmrt    
*     npoptn   npprt    nprset   npsavr   npsetx   npsrch   npupdt
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsol ( n, nclin, ncnln, ldA, ldcJu, ldR,
     $                   A, bl, bu,
     $                   funcon, funobj,
     $                   inform, iter, istate,
     $                   c, cJacu, clamda, objf, gradu, R, x,
     $                   iw, leniw, w, lenw )

      implicit           double precision (a-h,o-z)
      external           funcon, funobj
      integer            istate(n+nclin+ncnln)
      integer            iw(leniw)
      double precision   A(ldA,*), bl(n+nclin+ncnln),
     $                   bu(n+nclin+ncnln)
      double precision   c(*), cJacu(ldcJu,*), clamda(n+nclin+ncnln)
      double precision   gradu(n), R(ldR,*), x(n)
      double precision   w(lenw)

*     ==================================================================
*     npsol   solves the nonlinear program
*
*            minimize                   f(x)
*
*                                    (    x  )
*            subject to    bl  .le.  (  A*x  )  .le.  bu
*                                    (  c(x) )
*
*     where  f(x)  is a smooth scalar function,  A  is a constant matrix
*     and  c(x)  is a vector of smooth nonlinear functions.  The 
*     feasible region is defined by a mixture of linear and nonlinear 
*     equality or inequality constraints on  x.
*
*     The dimensions of the problem are...
*
*     n        the number of variables (dimension of  x),
*
*     nclin    the number of linear constraints (rows of the matrix  A),
*
*     ncnln    the number of nonlinear constraints (dimension of  c(x)),
*
*
*     npsol   uses a sequential quadratic programming algorithm, with a
*     positive-definite quasi-Newton approximation to the transformed
*     Hessian  Q'HQ  of the Lagrangian function (which will be stored in
*     the array  R).
*
*
*     Complete documentation for  NPSOL  is contained in Report
*     SOL 86-2, Users guide for NPSOL (Version 4.0), by P.E. Gill,
*     W. Murray, M.A. Saunders and M.H. Wright, Department of Operations
*     Research,  Stanford University, Stanford, California 94305.
*
*     Systems Optimization Laboratory, Stanford University.
*     Version 1.1,  April     12, 1983. The less said about this one...
*     Version 2.0,  April     30, 1984.
*     Version 3.0,  March     20, 1985. First Fortran 77 version
*     Version 3.2,  August    20, 1985.
*     Version 4.0,  April     16, 1986. First version with differences
*     Version 4.01, June      30, 1986. Level 2 BLAS + F77 linesearch
*     Version 4.02, August     5, 1986. Reset SSBFGS. One call to chfd
*     Version 4.03, June      14, 1987. Step limit
*     Version 4.04, June      28, 1989. Vectorizable BLAS
*     Version 4.05, November  28, 1989. Load and save files added
*     Version 4.06, November   5, 1991. srchq and srchc updated
*                   October   29, 1992. Summary/print file option.
*
*     Copyright  1983/1993  Stanford University.
*
*     This material may be reproduced by or for the U.S. Government 
*     pursuant to the copyright license under DAR Clause 7-104.9(a)
*     (1979 Mar).
*
*     This material is based upon work partially supported by the 
*     National Science Foundation under Grants MCS-7926009 and 
*     ECS-8312142; the Department of Energy Contract AM03-76SF00326,
*     PA No. DE-AT03-76ER72018; the Army Research Office Contract
*     DAA29-84-K-0156; and the Office of Naval Research Grant
*     N00014-75-C-0267.
*     ==================================================================
*     ------------------------------------------------------------------
*     Common blocks.
*     ------------------------------------------------------------------
*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     +Include npparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformH, lvlder, lverfy, msgNP , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgQP )

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
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
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            cmdbg, lsdbg, npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

      intrinsic          abs   , max   , min   , mod   , sqrt  , real

*     Local variables.

      external           ddiv  , dnrm2
      character*8        names(1)
      logical            cold  , linobj, named , overfl, rowerr, vertex
      logical            needfd
      parameter         (zero   =0.0d+0, point3 =3.3d-1, point8 =0.8d+0)
      parameter         (point9 =0.9d+0, one    =1.0d+0                )
      parameter         (hundrd =1.0d+2, growth =1.0d+2                )

      character*40       title
      data               title
     $                 / 'NPSOL  ---  Version 4.06-2     Nov  1992' /
      
*     Set the machine-dependent constants.

      call mchpar()

      epsmch = wmach( 3)
      rteps  = wmach( 4)
      nout   = wmach(11)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      rhomax = one/epsmch
      rootn  = sqrt(real(n))

*     Default names will be provided for variables during printing.
      
      named  = .false.
      inform = 0

*     Set the default values for the parameters.

      call npdflt( n, nclin, ncnln, leniw, lenw, title )

      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or. (lvlder .eq. 1  .and.  ncnln .gt. 0)
      cold   = lcrash .eq. 0
      lvldif = 0
      if (needfd) lvldif = 1

      nplin  = n     + nclin
      nctotl = nplin + ncnln

*     Assign the dimensions of arrays in the parameter list of NPCORE.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter MINACT and MINFXD
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
      if (ncnln .eq. 0  .and.  nclin .gt. 0) ldAqp = ldA

*     nploc  defines the arrays that contain the locations of various
*     work arrays within  w  and  iw.

      litotl = 0
      lwtotl = 0
      
      call nploc( n, nclin, ncnln, nctotl, litotl, lwtotl)

*     Allocate certain addresses that are not allocated in NPLOC.

      lax    = lwtotl + 1
      lwtotl = lax    + nclin - 1
      lAx    = min( lAx, lwtotl )

*     Check input parameters and storage limits.
      call cmchk ( nerror, msgNP, lcrash, .false.,
     $             leniw, lenw, litotl, lwtotl,
     $             n, nclin, ncnln,
     $             istate, iw, named, names,
     $             bigbnd, bl, bu, x )

      if (nerror .gt. 0) then
         inform = 9
         go to 800
      end if

      lkactv = locls( 1)
      lAnorm = locls( 2)
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
      lrho   = locnp(20)
      lwrk3  = locnp(21)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)
      lcJac  = locnp(27)
      lgrad  = locnp(28)

      ldcJ   = max ( ncnln, 1 )

      tolrnk = one / hcndbd
      Rcndbd = sqrt( hcndbd )

*     ==================================================================
*     If a unit number for a load file has been set, read initial values
*     from an old run.  These values override existing settings. 
*     ==================================================================
      if (nload .gt. 0) then
         call npgetr( nerror, unitQ, n, nclin, ncnln, ldR, ldQ,
     $                nfree0, iter, istate, iw(lkx),
     $                w(lhfrwd), w(lhctrl),
     $                w(lcmul), R, w(lrho), x, w(lQ) )

         if (nerror .gt. 0) then
            inform = 9
            go to 800
         end if
      end if

*     ==================================================================
*     Load the arrays of feasibility tolerances.
*     ==================================================================
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

*     ------------------------------------------------------------------
*     If required,  compute the problem functions.
*     If the constraints are nonlinear,  the first call of funcon
*     sets up any constant elements in the Jacobian matrix.  A copy of
*     the Jacobian (with constant elements set) is placed in  cJacu.
*     ------------------------------------------------------------------
      if (lverfy .ge. 10) then
         xnorm  = dnrm2 ( n, x, 1 )
         lvrfyc = lverfy - 10
         call npchkd( info, msgNP, nstate, lvlder, nfun, ngrad,
     $                ldcJ, ldcJu, n, ncnln,
     $                funcon, funobj, iw(lneedc),
     $                bigbnd, epsrf, cdint, fdint,
     $                fdchk, fdnorm, objf, xnorm,
     $                bl, bu, c, w(lwrk3), w(lcJac), cJacu, w(lcJdx),
     $                w(ldx), w(lgrad), gradu, w(lhfrwd), w(lhctrl),
     $                x, w(lwrk1), w(lwrk2), w, lenw )

         if (info .ne. 0) then
            if (info .gt. 0) inform = 7
            if (info .lt. 0) inform = info
            go to 800
         end if
         nstate = 0
      end if
      
      call icopy ( ldbg, ilsdbg, 1, icmdbg, 1 )

      if (nclin .gt. 0) then
         ianrmj = lAnorm
         do 110, j = 1, nclin
            w(ianrmj) = dnrm2 ( n, A(j,1), ldA )
            ianrmj    = ianrmj + 1
  110    continue
         call dcond ( nclin, w(lAnorm), 1, Asize, Amin )
      end if

      call dcond ( nplin, w(lfeatl), 1, feamax, feamin )
      call dcopy ( nplin, w(lfeatl), 1, w(lwtinf), 1 )
      call dscal ( nplin, (one/feamin), w(lwtinf), 1 )

*     ==================================================================
*     The input values of x and (optionally)  istate are used by
*     lscrsh  to define an initial working set.
*     ==================================================================
      vertex = .false.
      call lscrsh( cold, vertex,
     $             nclin, nplin, nactiv, nartif,
     $             nfree, n, ldA,
     $             istate, iw(lkactv),
     $             bigbnd, tolact,
     $             A, w(lax), bl, bu, x, w(lwrk1), w(lwrk2) )

      nres   = 0
      ngq    = 0
C-->  condmx = max( one/epspt5, hundrd )
C-->  condmx = max( one/epspt3, hundrd )
      condmx = max( one/epspt5, hundrd )

      if (lcrash .le. 1) then
*        ===============================================================
*        Cold or warm start. The upper-triangular matrix R is the factor
*        of an approximate Lagrangian Hessian.
*        ===============================================================
         unitQ  = .true.
         iter   = 0

         ikx    = lkx
         do 120, i = 1, n
            iw(ikx) = i
            ikx     = ikx + 1
  120    continue

         if (cold) then
            call f06qhf( 'Upper-triangular', n, n, zero, one, R, ldR )
            Rfrobn = rootn

            nrank  = 0
            if (ncnln .gt. 0) call dload ( ncnln, (zero), w(lcmul), 1 )
         else

*           R will be updated while finding a feasible x.

            nrank  = nlnx
            call dload ( nlnx, (zero), w(lres0), 1 )
            if (ncnln .gt. 0)
     $         call dcopy ( ncnln, clamda(nplin+1), 1, w(lcmul), 1 )

         end if

         incrun = .true.
         rhonrm =  zero
         rhodmp =  one
         scale  =  one
         
         call dload ( ncnln, (zero), w(lrho), 1 )

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

      else
*        ===============================================================
*        Hot start.
*        Stop if the computed and input values of nfree don't match.
*        ===============================================================
         if (nfree0 .ne. nfree) then
            nerror = 1
            inform = 9
            go to 800
         end if
      end if

*     ------------------------------------------------------------------
*     Factorize the linear constraints in the initial working set.
*     ------------------------------------------------------------------
      if (nactiv .gt. 0) then
         nact1  = nactiv
         nactiv = 0

         call lsadds( unitQ, vertex,
     $                inform, 1, nact1, nactiv, nartif, nZ, nfree,
     $                nrank, nrejtd, nres, ngq,
     $                n, ldQ, ldA, ldR, ldT,
     $                istate, iw(lkactv), iw(lkx), condmx,
     $                A, R, w(lT), w(lres0), w(lgq), w(lQ),
     $                w(lwrk1), w(lwrk2), w(lrlam) )
      end if

      if (lcrash .le. 1) then
*        ===============================================================
*        Cold or warm start.  Move  x  on to the linear constraints and
*        find a feasible point.
*        ===============================================================
         ssq1   =  zero
         linobj = .false.
         call lssetx( linobj, rowerr, unitQ,
     $                nclin, nactiv, nfree, nrank, nZ,
     $                n, nplin, ldQ, ldA, ldR, ldT,
     $                istate, iw(lkactv), iw(lkx),
     $                jmax, errmax, ctx, xnorm,
     $                A, w(lAx), bl, bu, w(lgq), w(lres), w(lres0),
     $                w(lfeatl), R, w(lT), x, w(lQ),w(lwrk1),w(lwrk2) )

*        ---------------------------------------------------------------
*        Call  lscore  to find a feasible  x.
*        ---------------------------------------------------------------
*        Use  work2  as the multiplier vector.

         jinf   = 0
         lclam  = lwrk2

         idbgsv = idbg
         if (idbg .gt. 0) then
            idbg = nminor + 1
         end if

         call lscore( 'FP problem', named, names, linobj, unitQ,
     $                nLPerr, itns, jinf, nclin, nplin,
     $                nactiv, nfree, nrank, nZ, nZr,
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

         idbg  = idbgsv
         call icopy ( ldbg, inpdbg, 1, icmdbg, 1 )
      else
*        ---------------------------------------------------------------
*        Hot start.
*        The point  x  is preassigned.  Compute the 2-norm of  x.
*        Initialize  Ax  for the linear constraints.
*        ---------------------------------------------------------------
         nrank  = nlnx
         xnorm  = dnrm2 ( n, x, 1 )
         if (nclin .gt. 0)
     $      call dgemv ( 'N', nclin, n, one, A, ldA,
     $                   x, 1, zero, w(lAx), 1 )
      end if

      if (lcrash .gt. 0) then

*        Check for a bad R.

         Rfrobn = f06qgf( 'Frobenius norm', 'Upper', n, n, R, ldR )
         call dcond ( n, R, ldR+1, dRmax, dRmin )
         cond   = ddiv  ( dRmax, dRmin, overfl )

         if (      cond   .gt. Rcndbd
     $       .or.  Rfrobn .gt. rootn*growth*dRmax) then
*           ------------------------------------------------------------
*           Refactorize the Hessian and bound the condition estimator.
*           ------------------------------------------------------------
            if (iPrint .gt. 0) write(iPrint, 9000)
            if (iSumm  .gt. 0) write(iSumm , 9000)

            call nprset( unitQ,
     $                   n, nfree, nZ, ldQ, ldR,
     $                   iw(liperm), iw(lkx),
     $                   w(lgq), R, w(lQ), w(lwrk1), w(lres0) )
         end if
      end if
      
*     ==================================================================
*     Now we can check the gradients at a feasible x.
*     ==================================================================
      lvrfyc = lverfy
      if (lverfy .ge. 10) lvrfyc = -1

      call npchkd( info, msgNP, nstate, lvlder, nfun, ngrad,
     $             ldcJ, ldcJu, n, ncnln,
     $             funcon, funobj, iw(lneedc),
     $             bigbnd, epsrf, cdint, fdint,
     $             fdchk, fdnorm, objf, xnorm,
     $             bl, bu, c, w(lwrk3), w(lcJac), cJacu, w(lcJdx),
     $             w(ldx), w(lgrad), gradu, w(lhfrwd), w(lhctrl),
     $             x, w(lwrk1), w(lwrk2), w, lenw )

      if (info .ne. 0) then
         if (info .gt. 0) inform = 7
         if (info .lt. 0) inform = info
         go to 800
      end if
      
      call dcopy ( n, w(lgrad), 1, w(lgq), 1 )
      call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $             iw(lkx), w(lgq), w(lQ), w(lwrk1) )

*     ==================================================================
*     Solve the problem.
*     ==================================================================
      if (ncnln .eq. 0) then
*        ---------------------------------------------------------------
*        The problem has only linear constraints and bounds.
*        ---------------------------------------------------------------
         call npcore( named, names, unitQ, inform, iter,
     $                n, nclin, ncnln, nctotl, nactiv, nfree, nZ,
     $                ldcJ, ldcJu, ldAqp, ldR,
     $                nfun, ngrad, istate, iw(lkactv), iw(lkx),
     $                objf, fdnorm, xnorm, funcon, funobj,
     $                A, w(lAx), bl, bu, c, w(lcJac), cJacu, clamda,
     $                w(lfeatl), w(lgrad), gradu, R, x, iw, w, lenw )
      else
*        ---------------------------------------------------------------
*        The problem has some nonlinear constraints.
*        ---------------------------------------------------------------
         if (nclin .gt. 0)
     $      call f06qff( 'General', nclin, n, A, ldA, w(lAqp), ldAqp )

*        Try and add some nonlinear constraint indices to KACTIV.
*
         call npcrsh( cold, n, nclin, ncnln,
     $                nctotl, nactiv, nfree, nZ,
     $                istate, iw(lkactv), bigbnd, tolact,
     $                bl, bu, c )

         call npcore( named, names, unitQ, inform, iter,
     $                n, nclin, ncnln, nctotl, nactiv, nfree, nZ,
     $                ldcJ, ldcJu, ldAqp, ldR,
     $                nfun, ngrad, istate, iw(lkactv),iw(lkx),
     $                objf, fdnorm, xnorm, funcon, funobj,
     $                w(lAqp), w(lAx), bl, bu, c, w(lcJac),cJacu,clamda,
     $                w(lfeatl), w(lgrad), gradu, R, x, iw, w, lenw )

      end if

*     ==================================================================
*     If a unit number for a save file has been set, save the details of
*     this run.
*     ==================================================================
      if (nsave .gt. 0  .and.  ksave .gt. nmajor) then
         call npsavr( unitQ, n, nclin, ncnln, ldR, ldQ,
     $                nfree, nsave, iter, istate, iw(lkx),
     $                w(lhfrwd), w(lhctrl),
     $                w(lcmul), R, w(lrho), x, w(lQ) )
      end if

*     ------------------------------------------------------------------
*     If required, form the triangular factor of the Hessian.
*     ------------------------------------------------------------------
*     First,  form the square matrix  R  such that  H = R'R.
*     Compute the  QR  factorization of  R.

      if (lformH .gt. 0) then
         lv     = lwrk2
         do 400 j = 1, n
            if (j .gt. 1)
     $         call dload ( j-1, zero, w(lv), 1 )

            lvj = lv + j - 1
            call dcopy ( n-j+1, R(j,j), ldR, w(lvj), 1     )
            call cmqmul( 3, n, nZ, nfree, ldQ, unitQ,
     $                   iw(lkx), w(lv), w(lQ), w(lwrk1) )
            call dcopy ( n    , w(lv) , 1    , R(j,1), ldR )
  400    continue
         
         call dgeqr ( n, n, R, ldR, w(lwrk1), info )
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
*     Load the final constraint and objective gradients.

      call icopy ( mxparm, ipsvls, 1, iprmls, 1 )
      call dcopy ( mxparm, rpsvls, 1, rprmls, 1 )
      call icopy ( mxparm, ipsvnp, 1, iprmnp, 1 )
      call dcopy ( mxparm, rpsvnp, 1, rprmnp, 1 )

      if (inform .lt. 9) then
         if (ncnln .gt. 0) then
            call f06qff( 'General', ncnln, n, w(lcJac), ldcJ, 
     $                   cJacu, ldcJu)
         end if
         call dcopy ( n, w(lgrad), 1, gradu, 1 )
      end if

      return

 3000 format(/ ' Exit NPSOL - User requested termination.'          )
 4000 format(/ ' Exit NPSOL - Optimal solution found.'              )
 4100 format(/ ' Exit NPSOL - Optimal solution found, ',
     $         ' but the requested accuracy could not be achieved.' )
 4200 format(/ ' Exit NPSOL - No feasible point for the linear',
     $         ' constraints.')
 4300 format(/ ' Exit NPSOL - No feasible point for the nonlinear',
     $         ' constraints.')
 4400 format(/ ' Exit NPSOL - Too many major iterations.             ')
 4600 format(/ ' Exit NPSOL - Current point cannot be improved upon. ')
 4700 format(/ ' Exit NPSOL - Large errors found in the derivatives. ')

 4900 format(/ ' Exit NPSOL - ', i10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.')
 5000 format(/ ' Final nonlinear objective value =', g16.7 )
 5010 format(/ ' Minimum sum of infeasibilities =',  g16.7 )
 5020 format(/ ' Final sum of infeasibilities =',    g16.7 )

 7000 format(/ ' The linear constraints are feasible.')
 9000 format(  ' XXX  Bad initial Hessian,   R  refactorized.' )

*     end of npsol
      end

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
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg
      integer            lcmdbg
      parameter         (lcmdbg = 5)
      common    /cmdebg/ icmdbg(lcmdbg), cmdbg

      intrinsic          abs, max, min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      if (cmdbg  .and.  icmdbg(3) .gt. 0) write(iPrint, 1000)

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

         if (cmdbg  .and.  icmdbg(3) .gt. 0)
     $      write(iPrint, 1200) j, res, Adxi, alfa

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

      if (cmdbg  .and.  icmdbg(1) .gt. 0  .and.  inform .gt. 0)
     $   write(iPrint, 9010) alfa

      return

 1000 format(/ ' npalf  entered'
     $       / '    j            res             Ap           alfa '/)
 1200 format( i5, 3g15.5 )
 9010 format(/ ' //npalf //  No finite step.'
     $       / ' //npalf //             alfa'
     $       / ' //npalf //  ', g15.4 )

*     end of npalf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npchkd( inform, msgNP, nstate, lvlder, nfun, ngrad,
     $                   ldcJ, ldcJu, n, ncnln,
     $                   funcon, funobj, needc,
     $                   bigbnd, epsrf, cdint, fdint,
     $                   fdchk, fdnorm, objf, xnorm,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   dx, grad, gradu, hforwd, hcntrl,
     $                   x, wrk1, wrk2, w, lenw )

      implicit           double precision (a-h,o-z)
      integer            needc(*)
      double precision   c(*), c1(*), cJac(ldcJ,*), cJacu(ldcJu,*),
     $                   cJdx(*)
      double precision   bl(n), bu(n), dx(n), grad(n), gradu(n), x(n)
      double precision   hforwd(*), hcntrl(*)
      double precision   wrk1(n), wrk2(*), w(lenw)
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
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
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
*        must be estimated.  A record of the missing jacobian elements
*        is stored in  cJacu.

         needfd = lvlder .eq. 0  .or.  lvlder .eq. 1

         if (needfd)
     $      call f06qhf( 'General', ncnln, n, rdummy, rdummy,
     $                   cJacu, ldcJu )

         call iload ( ncnln, (1), needc, 1 )

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
*     Check whatever gradient elements have been provided.
*     ==================================================================
      if (lvrfyc .ge. 0) then
         if (ncset .gt. 0) then
            call chcJac( mode, lvlder, msgNP,
     $                   ncset, n, ncnln, ldcJ, ldcJu,
     $                   bigbnd, epsrf, epspt3, fdchk, xnorm,
     $                   funcon, needc,
     $                   bl, bu, c, c1, cJac, cJacu, cJdx,
     $                   dx, wrk2, x, wrk1, w, lenw )
            if (mode .lt. 0) go to 9999
            infocj = mode
         end if

         if (nfdiff .lt. n) then
            call chfgrd( mode, msgNP, n,
     $                   bigbnd, epsrf, epspt3, fdchk, objf, xnorm,
     $                   funobj,
     $                   bl, bu, grad, gradu, dx, x, wrk1, w, lenw )
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
     $                n, ncnln, ldcJ, ldcJu,
     $                bigbnd, epsrf, fdnorm, objf,
     $                funobj, funcon, needc,
     $                bl, bu, c, c1, cJdx, cJac, cJacu,
     $                grad, gradu, hforwd, hcntrl, x,
     $                dx, w, lenw )

         if (mode .lt. 0) go to 9999

         if (lfdset .gt. 0) then
            centrl = lvldif .eq. 2
            call npfd  ( centrl, mode,
     $                   ldcJ, ldcJu, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm, objf,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, cJdx, cJac, cJacu,
     $                   grad, gradu, hforwd, hcntrl, x,
     $                   w, lenw )

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

      subroutine npcore( named, names, unitQ, inform, majIts,
     $                   n, nclin, ncnln, nctotl, nactiv, nfree, nZ,
     $                   ldcJ, ldcJu, ldAqp, ldR,
     $                   nfun, ngrad, istate, kactiv, kx,
     $                   objf, fdnorm, xnorm, funcon, funobj,
     $                   Aqp, Ax, bl, bu, c, cJac, cJacu, clamda,
     $                   featol, grad, gradu, R, x, iw, w, lenw )

      implicit           double precision (a-h,o-z)
      logical            named
      integer            istate(*), kactiv(n), kx(n)
      integer            iw(*)
      double precision   Aqp(ldAqp,*), Ax(*), bl(nctotl), bu(nctotl),
     $                   c(*), cJac(ldcJ,*), cJacu(ldcJu,*)
      double precision   clamda(nctotl), featol(nctotl), grad(n),
     $                   gradu(n), R(ldR,*), x(n)
      double precision   w(lenw)
      external           funcon, funobj

      double precision   Asize, dTmax, dTmin
      character*8        names(*)
               
*     ==================================================================
*     npcore  is the core routine for  npsol,  a sequential quadratic
*     programming (SQP) method for nonlinearly constrained optimization.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version      February-1982.
*     This version of npcore dated 23-Dec-92.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT   , ncolT , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize , dTmax , dTmin
      common    /sol6cm/ Rcndbd, Rfrobn, dRmax, dRmin

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      parameter         (ldbg = 5)
      logical            cmdbg, npdbg
      common    /npdebg/ inpdbg(ldbg), npdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     +Include npparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformH, lvlder, lverfy, msgNP , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgQP )

      logical            goodgq, newgq
      logical            centrl, convrg, convpt, done  , error , feasqp
      logical            infeas, needfd, optiml, overfl, unitQ
      logical            KTcond(2)
      intrinsic          abs   , max   , min   , mod   , real  , sqrt
      external           ddiv  , ddot  , dnrm2

      character*5        MjrMsg
      parameter        ( zero = 0.0d+0, one = 1.0d+0                )
      parameter        ( growth=1.0d+2                              )

*     specify machine-dependent parameters.

      flmax  = wmach(7)
      rtmax  = wmach(8)

      lAnorm = locls( 2)
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
      lrho   = locnp(20)
      lwrk3  = locnp(21)
      lslk1  = locnp(22)
      lslk   = locnp(23)
      lneedc = locnp(24)
      lhfrwd = locnp(25)
      lhctrl = locnp(26)

      lcJac1 = lAqp   + nclin
      lcJdx  = lAdx   + nclin
      lvioln = lwrk3

*     Initialize

      MjrMsg = '     '
      nQPinf = 0

      majit0 = majIts
      nplin  = n     + nclin
      ncqp   = nclin + ncnln
      nl     = min( nplin + 1, nctotl )

      ldcJ1 = max( ncqp , 1 )

      needfd = lvlder .eq. 0  .or.  lvlder .eq. 2
     $                        .or. (lvlder .eq. 1  .and.  ncnln .gt. 0)

      alfa   = zero
      alfdx  = zero
      rtftol = sqrt( ftol )
      rootn  = sqrt( real(n) )

*     If debug printing is required,  turn off any extensive printing
*     until iteration  idbg.

      msgsv1 = msgNP
      msgsv2 = msgQP
      if (idbg .le. nmajor  .and.  idbg .gt. 0) then
         msgNP = 0
         if (msgsv1 .ge. 5) msgNP = 5
         msgQP = 0
         if (msgsv2 .ge. 5) msgQP = 5
      end if

*     ------------------------------------------------------------------
*     Information from the feasibility phase will be used to generate a
*     hot start for the first QP subproblem.
*     ------------------------------------------------------------------
      call dcopy ( nctotl, featol, 1, w(lqptol), 1 )

      nstate = 0

      objalf = objf
      if (ncnln .gt. 0) then
         objalf = objalf - ddot  ( ncnln, w(lcmul), 1, c, 1 )
      end if

      newgq  = .false.

**    ==================================================================
*+    repeat                             (until converged or error exit)

*        ===============================================================
*        See if we want to save the details of this iteration.
*        ===============================================================
  100    if (mod(majIts,ksave) .eq. 0 .and. majIts .ne. majit0) then
            call npsavr( unitQ, n, nclin, ncnln, ldR, ldQ,
     $                   nfree, nsave, majIts, istate, kx,
     $                   w(lhfrwd), w(lhctrl),
     $                   w(lcmul), R, w(lrho), x, w(lQ) )
         end if

         minIts = 0

**       ===============================================================
*+       repeat                         (Until a good gradient is found)

  110       centrl = lvldif .eq. 2

            if (newgq) then
               if (needfd) then
*                 ------------------------------------------------------
*                 Compute any missing gradient elements and the
*                 transformed gradient of the objective.
*                 ------------------------------------------------------
                  call npfd  ( centrl, mode,
     $                         ldcJ, ldcJu, n, ncnln,
     $                         bigbnd, cdint, fdint, fdnorm, objf,
     $                         funcon, funobj, iw(lneedc),
     $                         bl, bu, c, w(lwrk2), w(lwrk3),cJac,cJacu,
     $                         grad, gradu, w(lhfrwd), w(lhctrl), x,
     $                         w, lenw )
                  inform = mode
                  if (mode .lt. 0) go to 800
               end if

               call dcopy ( n, grad, 1, w(lgq), 1 )
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
     $                   dxnorm, gdx, qpcurv,
     $                   Aqp, w(lAdx), w(lAnorm), Ax, bl, bu,
     $                   c, cJac, clamda, w(lcmul), w(lcs1),
     $                   w(ldlam), w(ldslk), w(ldx), w(lbl), w(lbu),
     $                   w(lqptol), R, w(lrho), w(lslk), w(lvioln), x,
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
*           Compute the norms of the projected gradient and the
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
            call dcond ( nZ, R, ldR+1, dRzmax, dRzmin )
            condHz = ddiv  ( dRzmax, dRzmin, overfl )
            if (condHz .lt. rtmax) then
               condHz = condHz*condHz
            else
               condHz = flmax
            end if
         end if

*        ---------------------------------------------------------------
*        Test for convergence.
*        The point test  convpt  checks for a K-T point at the initial
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
         glf1   = gdx
         if (ncnln .gt. 0) then
            glf1   = glf1
     $                 - ddot( ncnln, w(lcJdx), 1, clamda(nl), 1 )

*           Compute the value and directional derivative of the
*           augmented Lagrangian merit function.
*           The penalty parameters may be increased or decreased.

            call npmrt ( feasqp, n, nclin, ncnln,
     $                   objalf, grdalf, qpcurv,
     $                   istate,
     $                   w(lcJdx), w(lcmul), w(lcs1),
     $                   w(ldlam), w(lrho), w(lvioln),
     $                   w(lwrk1), w(lwrk2) )
         end if

*        ===============================================================
*        Print the details of this iteration.
*        ===============================================================
         call npprt ( KTcond, convrg, MjrMsg, msgNP, msgQP,
     $                ldR, ldT, n, nclin, ncnln,
     $                nctotl, nactiv, linact, nlnact, nZ, nfree,
     $                majit0, majIts, minIts, istate, alfa, nfun,
     $                condHz, condH, condT, objalf, objf,
     $                gfnorm, gznorm, cvnorm,
     $                Ax, c, R, w(lT), w(lvioln), x, w(lwrk1) )

         alfa  = zero
         error = majIts .ge. nmajor

         if (.not. (done  .or.  error)) then
            majIts = majIts + 1

            if (majIts .eq. idbg  .and.  iPrint .gt. 0) then
               npdbg = .true.
               cmdbg =  npdbg
               msgNP =  msgsv1
               msgQP =  msgsv2
            end if

*           Make copies of information needed for the BFGS update.

            call dcopy ( n, x     , 1, w(lx1) , 1 )
            call dcopy ( n, w(lgq), 1, w(lgq1), 1 )

            if (ncnln .gt. 0) then
               call dcopy ( ncnln, w(lcJdx), 1, w(lcJdx1), 1 )
               call dcopy ( ncnln, w(lcmul), 1, w(lc1mul), 1 )
               call dcopy ( ncnln, w(lslk) , 1, w(lslk1) , 1 )
            end if

*           ============================================================
*           Compute the parameters for the linesearch.
*           ============================================================
*           alfmin is the smallest allowable step predicted by the QP
*           subproblem.

            alfmin = one
            if (.not. feasqp) alfmin = zero

*           ------------------------------------------------------------
*           alfmax is the largest feasible steplength subject to a user-
*           defined limit alflim on the change in x.
*           ------------------------------------------------------------
            if (ncnln .gt. 0  .and.  needfd) then
               alfmax = one
            else
               alfmax = ddiv  ( bigdx, dxnorm, overfl )
               call npalf ( info, n, nclin, ncnln,
     $                      alfa, alfmin, alfmax, bigbnd, dxnorm,
     $                      w(lAnorm), w(lAdx), Ax, bl, bu,
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

            call npsrch( needfd, nlserr, n, ncnln,
     $                   ldcJ, ldcJu, nfun, ngrad,
     $                   iw(lneedc), funcon, funobj,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, glf1, glf2,
     $                   objf, objalf, qpcurv, xnorm,
     $                   c, w(lwrk1), cJac, cJacu, w(lcJdx), w(lwrk3),
     $                   w(lc1mul), w(lcmul), w(lcs1),
     $                   w(lcs2), w(ldx), w(ldlam), w(ldslk), grad,
     $                   gradu, clamda(nl), w(lrho),
     $                   w(lslk1), w(lslk), w(lx1), x,
     $                   w(lwrk2), w, lenw )

*           ------------------------------------------------------------
*           npsrch  sets nlserr to the following values...
*
*           < 0  if the user wants to stop.
*             1  if the search is successful and alfa < alfmax.
*             2  if the search is successful and alfa = alfmax.
*             3  if a better point was found but too many functions
*                were needed (not sufficient decrease).
*
*           Values of nlserr occurring with a nonzero value of alfa.
*             4  if alfmax < tolabs (too small to do a search).
*             5  if alfa  < alfsml (srchq only -- maybe want to switch
*                to central differences to get a better direction).
*             6  if the search found that there is no useful step.
*                The interval of uncertainty is less than 2*tolabs.
*                The minimizer is very close to alfa = zero
*                or the gradients are not sufficiently accurate.
*             7  if there were too many function calls.
*             8  if the input parameters were bad
*                (alfmax le toltny  or  uphill).
*           ------------------------------------------------------------
            if (nlserr .lt. 0) then
               inform = nlserr
               go to 800
            end if

            if (alfa .gt. alflim) MjrMsg(4:4) = 'l'

            error  = nlserr .ge. 4
            if (error) then
*              ---------------------------------------------------------
*              The linesearch failed to find a better point.
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
*                 Compute the missing gradients.
*                 ======================================================
                  mode  = 1
                  ngrad = ngrad + 1

                  if (ncnln .gt. 0) then
                     call iload ( ncnln, (1), iw(lneedc), 1 )

                     call funcon( mode, ncnln, n, ldcJu, iw(lneedc),
     $                            x, w(lwrk1), cJacu, nstate )
                     inform = mode
                     if (mode .lt. 0) go to 800

                     call f06qff( 'General', ncnln, n, cJacu, ldcJu,
     $                            cJac, ldcJ )
                  end if

                  call funobj( mode, n, x, obj, gradu, nstate )
                  inform = mode
                  if (mode .lt. 0) go to 800

                  call dcopy ( n, gradu, 1, grad, 1 )

                  call npfd  ( centrl, mode,
     $                         ldcJ, ldcJu, n, ncnln,
     $                         bigbnd, cdint, fdint, fdnorm, objf,
     $                         funcon, funobj, iw(lneedc),
     $                         bl, bu, c, w(lwrk2), w(lwrk3),cJac,cJacu,
     $                         grad, gradu, w(lhfrwd), w(lhctrl), x,
     $                         w, lenw )

                  inform = mode
                  if (mode .lt. 0) go to 800

                  gdx  =  ddot( n, grad, 1, w(ldx), 1 )
                  glf2 =  gdx
                  if (ncnln .gt. 0) then
                     call dgemv ( 'N', ncnln, n, one, cJac, ldcJ,
     $                            w(ldx), 1, zero, w(lcJdx), 1 )
                     glf2 = glf2 -
     $                      ddot( ncnln, w(lcJdx), 1, clamda(nl), 1 )
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
               alfdx   = alfa * dxnorm

*              =========================================================
*              Update the factors of the approximate Hessian of the
*              Lagrangian function.
*              =========================================================
               call npupdt( MjrMsg, unitQ,
     $                      n, ncnln, nfree, nZ,
     $                      ldcJ1, ldcJ, ldQ, ldR, kx,
     $                      alfa, glf1, glf2, qpcurv,
     $                      w(lcJac1), cJac, w(lcJdx1), w(lcJdx),
     $                      w(lcs1), w(lcs2), w(lgq1), w(lgq),
     $                      w(lHpq), w(lRpq), clamda(nl), r,
     $                      w(lwrk3), w(lQ), w(lwrk2), w(lwrk1) )

               call dcond ( n, R, ldR+1, dRmax, dRmin )
               cond   = ddiv  ( dRmax, dRmin, overfl )

               if (      cond   .gt. Rcndbd
     $             .or.  Rfrobn .gt. rootn*growth*dRmax) then
*                 ------------------------------------------------------
*                 Reset the condition estimator and range-space
*                 partition of Q'HQ.
*                 ------------------------------------------------------
                  if (npdbg  .and.  inpdbg(1) .gt. 0) then
                     write(iPrint, 9000) Rfrobn, dRmax , dRmin,
     $                                   cond  , Rcndbd
                  end if

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
  800 msgNP = msgsv1
      msgQP = msgsv2
      if (msgNP .gt. 0  .and.  iPrint .gt. 0) then
         write(iPrint, 2100) inform, majIts, nfun, ngrad
      end if

      call cmprt ( msgNP, nfree, ldAqp,
     $             n, nclin, nctotl, bigbnd,
     $             named, names,
     $             nactiv, istate, kactiv, kx,
     $             Aqp, bl, bu, c, clamda, w(lrlam), x )
      if (ncnln .gt. 0)
     $call dcopy ( ncnln, w(lcmul), 1, clamda(n+nclin+1), 1 )

      return

 2100 format(/ ' Exit  NP phase.  Inform = ', i2, '  Majits = ', i4,
     $         '   nfun = ', i4, '   ngrad = ', i4 )
 3000 format(  ' Minor itn', i6, '.  Central-differences computed. ',
     $         ' QP re-solved.' )
 9000 format(/ ' //npcore//        Rfrobn         dRmax         dRmin'
     $       / ' //npcore//', 1p, 3e14.2
     $       / ' //npcore//          Cond        Rcndbd'
     $       / ' //npcore//', 1p, 2e14.2 )

*     end of npcore
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npcrsh( cold, n, nclin, ncnln,
     $                   nctotl, nactiv, nfree, nZ,
     $                   istate, kactiv, bigbnd, tolact,
     $                   bl, bu, c )

      implicit           double precision (a-h,o-z)
      logical            cold
      integer            istate(nctotl), kactiv(n)
      double precision   c(*), bl(nctotl), bu(nctotl)

*     ==================================================================
*     npcrsh  adds indices of nonlinear constraints to the initial
*     working set.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version   14-February 1985.
*     This version of  npcrsh  dated 14-November-1985.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      intrinsic          abs, min
      parameter        ( one = 1.0d+0 )

      nfixed = n      - nfree
      linact = nactiv
      nplin  = n      + nclin

*     If a cold start is being made, initialize the status of the QP
*     working set.  First,  if  BL(j) = BU(j),  set ISTATE(j)=3.

      if (cold) then
         do 130, j = nplin+1, nctotl
            istate(j) = 0
            if (bl(j) .eq. bu(j)) istate(j) = 3
  130    continue
      end if

*     Increment NACTIV and KACTIV.
*     Ensure that the number of bounds and general constraints in the
*     QP  working set does not exceed N.

      do 200, j = nplin+1, nctotl
         if (nfixed + nactiv .eq. n) istate(j) = 0
         if (istate(j) .gt. 0) then
            nactiv = nactiv + 1
            kactiv(nactiv) = j - n
         end if
  200 continue

      if (cold) then

*        ---------------------------------------------------------------
*        If a cold start is required, an attempt is made to add as many
*        nonlinear constraints as possible to the working set.
*        ---------------------------------------------------------------
*        The following loop finds the most violated constraint.  If
*        there is room in kactiv, it will be added to the working set
*        and the process will be repeated.

         is     =   1
         biglow = - bigbnd
         bigupp =   bigbnd
         toobig =   tolact + tolact

*        while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
  500    if    (is .gt. 0  .and.  nfixed + nactiv .lt. n) then
            is   = 0
            cmin = tolact

            do 520, i = 1, ncnln
               j      = nplin + i
               if (istate(j) .eq. 0) then
                  b1     = bl(j)
                  b2     = bu(j)
                  resl   = toobig
                  resu   = toobig
                  if (b1 .gt. biglow)
     $            resl   = abs( c(i) - b1 ) / (one + abs( b1 ))
                  if (b2 .lt. bigupp)
     $            resu   = abs( c(i) - b2 ) / (one + abs( b2 ))
                  res    = min( resl, resu )
                  if (res .lt. cmin) then
                     cmin = res
                     imin = i
                     is   = 1
                     if (resl .gt. resu) is = 2
                  end if
               end if
  520       continue

            if (is .gt. 0) then
               nactiv         = nactiv + 1
               kactiv(nactiv) = nclin  + imin
               j              = nplin  + imin
               istate(j)      = is
            end if
            go to 500
*        end while
         end if
      end if

*     ------------------------------------------------------------------
*     An initial working set has now been selected.
*     ------------------------------------------------------------------
      nlnact = nactiv - linact
      nZ     = nfree  - nactiv
      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1000) nfixed, linact, nlnact

      return

 1000 format(/ ' //npcrsh//  Working set selected....'
     $       / ' //npcrsh// nfixed linact nlnact     '
     $       / ' //npcrsh//', 3i7 )

*     end of npcrsh
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npdflt( n, nclin, ncnln, leniw, lenw, title )

      implicit           double precision (a-h,o-z)

      character*(*)      title

*     ==================================================================
*     npdflt  loads the default values of parameters not set in the
*     options file.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written 10-September-1985.
*     This version of npdflt dated 09-Nov-92.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
                 
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      common    /sol5np/ lvrfyc, jverfy(4)

      logical            cmdbg, lsdbg, npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

      logical            newopt
      common    /sol7np/ newopt
      save      /sol7np/

*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     +Include npparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformH, lvlder, lverfy, msgNP , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgQP )

      intrinsic          abs    , len    , max   , mod   , real
      parameter        ( zero   =  0.0d+0, one    =  1.0d+0 )
      parameter        ( point3 =  3.3d-1, point8 =  0.8d+0 )
      parameter        ( point9 =  0.9d+0, two    =  2.0d+0 )
      parameter        ( tenp6  =  1.0d+6, hundrd = 10.0d+1 )
      parameter        ( rdummy = -11111., idummy = -11111  )
      parameter        ( gigant =  1.0d+20*.99999           )
      parameter        ( wrktol =  1.0d-2                   )

      character*4        icrsh(0:2)
      character*16       key
      data                icrsh(0),  icrsh(1),  icrsh(2)
     $                 / 'cold'   , 'warm'   , 'hot '    /

      epsmch = wmach( 3)
      nout   = wmach(11)

      condbd = max ( one/(hundrd*epsmch*real(n)), tenp6 )

      nplin  = n     + nclin
      nctotl = nplin + ncnln

*     Make a dummy call npkey to ensure that the defaults are set.

      call npkey ( nout, '*', key )
      newopt = .true.

*     Save the optional parameters set by the user.  The values in
*     iprmls, rprmls, iprmnp and rprmnp may be changed to their
*     default values.

      call icopy ( mxparm, iprmls, 1, ipsvls, 1 )
      call dcopy ( mxparm, rprmls, 1, rpsvls, 1 )
      call icopy ( mxparm, iprmnp, 1, ipsvnp, 1 )
      call dcopy ( mxparm, rprmnp, 1, rpsvnp, 1 )

      if (          iPrnt  .lt. 0     )   iPrnt   = nout
      if (          iSumry .lt. 0     )   iSumry  = 6
      if (          iSumry .eq. iPrnt )   iSumry  = 0
                                          iPrint  = iPrnt
                                          iSumm   = iSumry
      if (          lcrash .lt. 0
     $    .or.      lcrash .gt. 2     )   lcrash  =  0
      if (          lvlder .lt. 0
     $    .or.      lvlder .gt. 3     )   lvlder  =  3
      if (          lformH .lt. 0
     $    .or.      lformH .gt. 1     )   lformH  =  0

      if (          nmajor .lt. 0     )   nmajor  = max(50, 3*nplin+
     $                                                     10*ncnln )
      if (          nminor .lt. 1     )   nminor  = max(50, 3*nctotl)
      if (          mjrdbg .lt. 0     )   mjrdbg  =  0
      if (          mnrdbg .lt. 0     )   mnrdbg  =  0
      if (          idbg   .lt. 0
     $    .or.      idbg   .gt. nmajor)   idbg    =  0
      if (          mjrdbg .eq. 0
     $    .and.     mnrdbg .eq. 0     )   idbg    = nmajor + 1
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
      if (          lverfy .eq. idummy
     $    .or.      lverfy .gt. 13    )   lverfy  =  0

      if (          ksave  .le. 0     )   ksave   =  nmajor + 1
      if (          nsave  .lt. 0     )   nsave   =  0
      if (          nsave  .eq. 0     )   ksave   =  nmajor + 1
      if (          nload  .lt. 0     )   nload   =  0
      if (          lcrash .le. 1     )   nload   =  0
      if (          nload  .eq. 0 
     $    .and.     lcrash .eq. 2     )   lcrash  =  0

      if (          tolact .lt. zero
     $    .or.      tolact .ge. one   )   tolact  =  wrktol
      if (          tolfea .lt. epsmch
     $    .or.      tolfea .ge. one   )   tolfea  =  epspt5
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

      itmax1    = max( 50, 3*(n + nclin + ncnln) )
      jverfy(1) = jvrfy1
      jverfy(2) = jvrfy2
      jverfy(3) = jvrfy3
      jverfy(4) = jvrfy4

      npdbg = idbg .eq. 0  .and.  iPrint .gt. 0
      cmdbg = npdbg
      lsdbg = npdbg

      k     = 1
      msg1  = mjrdbg
      msg2  = mnrdbg
      do 200, i = 1, ldbg
         inpdbg(i) = mod( msg1/k, 10 )
         icmdbg(i) = inpdbg(i)
         ilsdbg(i) = mod( msg2/k, 10 )
         k = k*10
  200 continue

      if (msgNP .gt. 0) then

*        Print the title.  If no hot start is specified,  the parameters
*        are final and can be printed.

         lenT = len( title )
         nspace = (81 - lenT)/2 + 1
         if (iPrint .gt. 0) then
            write(iPrint, '(///// (80a1) )')
     $            (' ', j=1, nspace), (title(j:j), j=1,lenT)
            write(iPrint, '(80a1 //)')
     $            (' ', j=1, nspace), ('='       , j=1,lenT)
         end if

         if (iSumm .gt. 0) then
            write(iSumm, '(///// (80a1) )')
     $            (' ', j=1, nspace), (title(j:j), j=1,lenT)
            write(iSumm, '(80a1 //)')
     $            (' ', j=1, nspace), ('='       , j=1,lenT)
         end if

         if (iPrint .gt. 0  .and.  lcrash .le. 1) then
            write(iPrint, 2000)
            write(iPrint, 2100) nclin , tolfea, icrsh(lcrash) ,
     $                          n     , bigbnd, tolact,
     $                          dxlim , bigdx
            write(iPrint, 2200) ncnln , ftol  , epsrf ,
     $                          nlnj  , ctol  , epsmch,
     $                          nlnf  , eta   , iPrnt ,
     $                          lvlder, lverfy, iSumry
            write(iPrint, 2300) nmajor, msgNP ,
     $                          nminor, msgQP ,
     $                          nload , nsave , ksave
            
            if (lvlder .lt. 3) then
               if      (lfdset .eq. 0) then
                  write(iPrint, 2400)
               else if (lfdset .eq. 1) then
                  write(iPrint, 2401) fdint, cdint
               else if (lfdset .eq. 2) then
                  write(iPrint, 2402)
               end if
            end if
         end if
      end if

      return

 2000 format(
     $//' Parameters'
     $/ ' ----------' )
 2100 format(
     $/ ' Linear constraints.....',     i10,   6x,
     $  ' Linear feasibility.....', 1p, e10.2, 6x,
     $  1x, a4, ' start.............'
     $/ ' Variables..............',     i10,   6x,
     $  ' Infinite bound size....', 1p, e10.2, 6x,
     $  ' Crash tolerance........',     e10.2
     $/ ' Step limit.............', 1p, e10.2, 6x,
     $  ' Infinite step size.....',     e10.2  )
 2200 format(
     $/ ' Nonlinear constraints..',     i10,   6x,
     $  ' Optimality tolerance...', 1p, e10.2, 6x,
     $  ' Function precision.....',     e10.2
     $/ ' Nonlinear Jacobian vars',     i10,   6x,
     $  ' Nonlinear feasibility..', 1p, e10.2, 6x,
     $  ' eps (machine precision)', 1p, e10.2
     $/ ' Nonlinear objectiv vars',     i10,   6x,
     $  ' Linesearch tolerance...', 1p, e10.2, 6x,
     $  ' Print file.............',     i10
     $/ ' Derivative level.......',     i10,   6x,
     $  ' Verify level...........',     i10,   6x,
     $  ' Summary file...........',     i10)
 2300 format(
     $/ ' Major iterations limit.',     i10, 6x,
     $  ' Major print level......',     i10
     $/ ' Minor iterations limit.',     i10, 6x,
     $  ' Minor print level......',     i10
     $/ ' RUN loaded from file...',     i10, 6x,
     $  ' RUN to be saved on file',     i10, 6x,
     $  ' Save frequency.........',     i10)

 2400 format(/ ' Difference intervals to be computed.' )
 2401 format(/ ' Difference interval....', 1p, e10.2, 6x,
     $         ' Central diffce interval',     e10.2 )
 2402 format(/ ' User-supplied difference intervals.' )

*     end of npdflt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfd  ( centrl, inform,
     $                   ldcJ, ldcJu, n, ncnln,
     $                   bigbnd, cdint, fdint, fdnorm, objf,
     $                   funcon, funobj, needc,
     $                   bl, bu, c, c1, c2, cJac, cJacu,
     $                   grad, gradu, hforwd, hcntrl, x,
     $                   w, lenw )

      implicit           double precision (a-h,o-z)
      logical            centrl
      integer            needc(*)

      double precision   bl(n), bu(n), c(*), c1(*), c2(*),
     $                   cJac(ldcJ,*), cJacu(ldcJu,*)
      double precision   grad(n), gradu(n), hforwd(n), hcntrl(n), x(n)
      double precision   w(lenw)
      external           funcon, funobj

*     ==================================================================
*     npfd   evaluates any missing gradients.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written 3-July-1986.
*     This version of npfd   dated 14-Sep-92.
*     ==================================================================
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      intrinsic          abs   , max

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
               call funcon( mode, ncnln, n, ldcJu,
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
                  call funcon( mode, ncnln, n, ldcJu,
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
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      external           idamax, dnrm2
      intrinsic          abs
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

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1000) errmax, jmax

      cvnorm = dnrm2 ( n+nclin+ncnln, work, 1 )

      return

 1000 format(/ ' //npfeas//  The maximum violation is ', 1p, e14.2,
     $                     ' in constraint', i5 )

*     end of npfeas
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npfile( ioptns, inform )
      integer            ioptns, inform

*     ==================================================================
*     npfile  reads the options file from unit  ioptns  and loads the
*     options into the relevant elements of  iprmnp  and  rprmnp.
*
*     If  ioptns .lt. 0  or  ioptns .gt. 99  then no file is read,
*     otherwise the file associated with unit  ioptns  is read.
*
*     Output:
*
*         inform = 0  if a complete  options  file was found
*                     (starting with  begin  and ending with  end);
*                  1  if  ioptns .lt. 0  or  ioptns .gt. 99;
*                  2  if  begin  was found, but end-of-file
*                     occurred before  end  was found;
*                  3  if end-of-file occurred before  begin  or
*                     endrun  were found;
*                  4  if  endrun  was found before  begin.
*     ==================================================================
      logical             newopt
      common     /sol7np/ newopt
      save       /sol7np/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar, npkey
      logical             first
      save                first , nout
      data                first /.true./

*     If first time in, set  nout.
*     newopt is true first time into npfile or npoptn
*     and just after a call to npsol.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout   = wmach(11)
      end if

      call opfile( ioptns, nout, inform, npkey )

      return

*     end of npfile
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npgetr( inform, unitQ, n, nclin, ncnln, ldR, ldQ,
     $                   nfree, iter, istate, kx,
     $                   hforwd, hcntrl,
     $                   cmul, R, rho, x, Q )

      implicit           double precision (a-h,o-z)
      logical            unitQ
      integer            istate(n+nclin+ncnln), kx(n)
      double precision   R(ldR,*), x(n), Q(ldQ,*)
      double precision   hforwd(*), hcntrl(*)
      double precision   cmul(*), rho(*)

*     ==================================================================
*     npgetr  loads details of a previous run from unit nload.
*
*     Original version   24-Nov-89.
*     This version of  npgetr  dated 14-Sep-92.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformH, lvlder, lverfy, msgNP , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgQP )

      character*4        icrsh(0:2)
      data                icrsh(0),  icrsh(1),  icrsh(2)
     $                 / 'COLD'   , 'WARM'   , 'HOT '    /

      if (nload .le. 0) return

      if (iPrint .gt. 0) write(iPrint, 4000) nload

      read(nload, 1000, end=999) iter, nfree, lfdset, lvldif, unitQ
      do 110, j = 1, n
         read(nload, 1010, end=999) jold, kx(j), istate(j), x(j)
  110 continue

      if (jold .ne. n) then
         if (iPrint .gt. 0) write(iPrint, 9000)
         inform = 1
         return
      end if

      do 120, j = n+1, n+nclin
         read(nload, 1020, end=999) jold, istate(j)
  120 continue

      if (jold .ne. n+nclin) then
         if (iPrint .gt. 0) write(iPrint, 9000)
         inform = 1
         return
      end if

      if (ncnln .gt. 0) then
         k = 1
         do 130, j = n+nclin+1, n+nclin+ncnln
            read(nload, 1030, end=999) jold, istate(j), cmul(k),
     $                                 rho(k)
            k = k + 1
  130    continue

         if (jold .ne. n+nclin+ncnln) then
            if (iPrint .gt. 0) write(iPrint, 9000)
            inform = 1
            return
         end if
     
         read(nload, 1040, end=999) rhomax, rhonrm, rhodmp, scale,
     $                              incrun
      end if

*     ------------------------------------------------------------------
*     Read   Q(free)  and the factor of  Q'HQ.
*     ------------------------------------------------------------------
      if (.not. unitQ) then
         do 160, j = 1, nfree
            do 150, i = 1, nfree
              read(nload, 1050, end=999) iold, jold, Q(i,j)
  150       continue     
  160    continue
      end if 

      do 180, j = 1, n
         do 170, i = 1, j
           read(nload, 1050, end=999) iold, jold, R(i,j)
  170    continue
  180 continue

      if (jold .ne. n  .or. iold .ne. n) then
         if (iPrint .gt. 0) write(iPrint, 9000)
         inform = 1
         return
      end if

*     ------------------------------------------------------------------
*     Read the finite-difference intervals.  
*     ------------------------------------------------------------------
      if (lvldif .gt. 0) then
         if (lfdset .eq. 0  .or.  lfdset .eq. 2) then
            do 190, j = 1, n
               read(nload, 1090, end=999) jold, hforwd(j), hcntrl(j)
  190       continue
         end if
      end if                                            

      if (iPrint .gt. 0) write(iPrint, 4010) n, iter

      call mcclos( nload )

*     ------------------------------------------------------------------
*     Now that all values have been set, we can print the parameters.
*     ------------------------------------------------------------------
      if (iPrint .gt. 0  .and.  msgNP .gt. 0) then
         epsmch = wmach( 3)
         write(iPrint, 2000)
         write(iPrint, 2100) nclin , tolfea, icrsh(lcrash) ,
     $                       n     , bigbnd, tolact,
     $                       dxlim , bigdx
         write(iPrint, 2200) ncnln , ftol  , epsrf ,
     $                       nlnj  , ctol  , epsmch,
     $                       nlnf  , eta   , iPrnt ,
     $                       lvlder, lverfy, iSumry
         write(iPrint, 2300) nmajor, msgNP ,
     $                       nminor, msgQP ,
     $                       nload , nsave , ksave

         if (lvlder .lt. 3) then
            if      (lfdset .eq. 0) then
               write(iPrint, 2400)
            else if (lfdset .eq. 1) then
               write(iPrint, 2401) fdint, cdint
            else if (lfdset .eq. 2) then
               write(iPrint, 2402)
            end if
         end if
      end if

      return
                       
 999  if (iPrint .gt. 0) write(iPrint, 9000)
      inform = 1 
      return

 1000 format(2i8, 1x, 2i2, 1x, l1 )
 1010 format(2i8, 1x,  i2, 1p, 2e24.14 )
 1020 format( i8, 1x,  i2 )
 1030 format( i8, 1x,  i2, 1p, 2e24.14 )
 1040 format( 1p, 4e24.14, 1x, l1 )
 1050 format(2i8, 1x, 1p,  e24.14 )
 1090 format( i8, 1x, 1p, 2e24.14 )

 2000 format(
     $//' Parameters'
     $/ ' ----------' )
 2100 format(
     $/ ' Linear constraints.....',     i10  , 6x,
     $  ' Linear feasibility.....', 1p, e10.2, 6x,
     $  1x, a4, ' start.............'
     $/ ' Variables..............',     i10,   6x,
     $  ' Infinite bound size....',     e10.2, 6x,
     $  ' Crash tolerance........',     e10.2
     $/ ' Step limit.............',     e10.2, 6x,
     $  ' Infinite step size.....',     e10.2  )
 2200 format(
     $/ ' Nonlinear constraints..',     i10,   6x,
     $  ' Optimality tolerance...', 1p, e10.2, 6x,
     $  ' Function precision.....',     e10.2
     $/ ' Nonlinear Jacobian vars',     i10,   6x,
     $  ' Nonlinear feasibility..', 1p, e10.2, 6x,
     $  ' eps (machine precision)', 1p, e10.2
     $/ ' Nonlinear objectiv vars',     i10,   6x,
     $  ' Linesearch tolerance...', 1p, e10.2, 6x,
     $  ' Print file.............',     i10
     $/ ' Derivative level.......',     i10,   6x,
     $  ' Verify level...........',     i10,   6x,
     $  ' Summary file...........',     i10)
 2300 format(
     $/ ' Major iterations limit.',     i10, 6x,
     $  ' Major print level......',     i10
     $/ ' Minor iterations limit.',     i10, 6x,
     $  ' Minor print level......',     i10
     $/ ' RUN loaded from file...',     i10, 6x,
     $  ' RUN to be saved on file',     i10, 6x,
     $  ' Save frequency.........',     i10)

 2400 format(/ ' Difference intervals to be computed.' )
 2401 format(/ ' Difference interval....', 1p, e10.2, 6x,
     $         ' Central diffce interval',     e10.2 )
 2402 format(/ ' User-supplied difference intervals.' )
 4000 format(/ ' OLD RUN to be loaded from file', i4)
 4010 format(/ ' Number of variables loaded = ', i5, '    Itn = ', i8)
 9000 format(/ ' Exit NPSOL - The OLD RUN file dimensions do not',
     $         ' match this problem.' )

*     end of npgetr
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npiqp ( feasqp, unitQ, nQPerr, majIts, minIts,
     $                   n, nclin, ncnln, ldcJ, ldAqp, ldR,
     $                   linact, nlnact, nactiv, nfree, nZ, numinf,
     $                   istate, kactiv, kx,
     $                   dxnorm, gdx, qpcurv,
     $                   Aqp, Adx, Anorm, Ax, bl, bu,
     $                   c, cJac, clamda, cmul, cs,
     $                   dlam, dslk, dx, qpbl, qpbu, qptol,
     $                   R, rho, slk, violn, x,
     $                   wtinf, iw, w )

      implicit           double precision (a-h,o-z)
      logical            feasqp, unitQ
      integer            istate(*), kactiv(n), kx(n)
      integer            iw(*)
      double precision   Aqp(ldAqp,*), Adx(*), Anorm(*), Ax(*),
     $                   bl(*), bu(*),
     $                   c(*), cJac(ldcJ,*), clamda(*), cmul(*), cs(*)
      double precision   dlam(*), dslk(*), dx(n)
      double precision   qpbl(*), qpbu(*),
     $                   qptol(*), R(ldR,*), rho(*), slk(*),
     $                   violn(*), x(n), wtinf(*)
      double precision   w(*)

*     ==================================================================
*     npiqp   does the following:
*
*     (1)  Generate the upper and lower bounds for the QP  subproblem.
*
*     (2)  Compute the  TQ  factors of the rows of  AQP  specified by
*          the array  ISTATE.  The part of the factorization defined by
*          the first contiguous group of linear constraints does not
*          need to be recomputed.  The remaining rows (which could be
*          comprised of both linear and nonlinear constraints) are
*          included as new rows of the  TQ  factorization stored in
*          T and Q.  Note that if there are no nonlinear constraints,
*          no factorization is required.
*
*     (3)  Solve the  QP  subproblem.
*                 minimize     1/2 (W p - d)'(Wp - d) + g'p
*
*                 subject to   qpbl .le. (  p ) .le. qpbu,
*                                        ( Ap )
*
*          where  W  is a matrix (not stored) such that  W'W = H  and
*          WQ = R,  d  is the zero vector,  and  g  is the gradient.
*          If the subproblem is infeasible, compute the point which
*          minimizes the sum of infeasibilities.
*
*    (4)   Find the value of each slack variable for which the merit
*          function is minimized.
*
*    (5)   Compute  dslk,  dlam  and  dx,  the search directions for
*          the slack variables, the multipliers and the variables.
*
*     Systems Optimization Laboratory, Stanford University.
*     Fortran 66 version written 10-January-1983.
*     This version of npiqp dated 23-Dec-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT   , ncolT , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize , dTmax , dTmin
      common    /sol6cm/ Rcndbd, Rfrobn, dRmax , dRmin
      
      integer            locls
      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            cmdbg, lsdbg, npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformH, lvlder, lverfy, msgNP , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgQP )

      character*8        names(1)
      logical            linobj, overfl, qpnamd, vertex
      parameter         (qpnamd = .false., vertex = .false. )
      parameter         (zero   = 0.0d+0, one = 1.0d+0, two = 2.0d+0)
      parameter         (hundrd = 1.0d+2                            )
      character*4        line
      data               line /'----'/

      idbgsv = idbg
      if (npdbg) then
         idbg   = 0
      else
         idbg = nminor + 1
      end if
      lsdbg  = npdbg
      cmdbg  = npdbg
      call icopy ( ldbg, ilsdbg, 1, icmdbg, 1 )

      lRpq   = locls( 5)
      lRpq0  = locls( 6)
      lHpq   = locls( 8)
      lgq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwrk1  = locls(14)

      nRpq   = 0
      ngq    = 1

      feasqp =  .true.
      linobj =  .true.

      biglow = - bigbnd
      bigupp =   bigbnd
      ssq1   =   zero

      nplin  = n     + nclin
      nctotl = nplin + ncnln
      ncqp   = nclin + ncnln
      nrank  = n
      nrejtd = 0                           

      if (msgQP .gt. 0) then
         if (iPrint .gt. 0) then
            write(iPrint, 1010) (line, j=1,32), majIts
         end if
         if (iSumm  .gt. 0) then
            write(iSumm , 1020) (line, j=1,16), majIts
         end if
      end if

*     ==================================================================
*     Generate the upper and lower bounds upon the search direction, the
*     weights on the sum of infeasibilities and the nonlinear constraint
*     violations.
*     ==================================================================
      wscale = - one
      do 170 j = 1, nctotl

         if (j .le. n) then                
            con = x(j)
         else if (j .le. nplin) then
            con = Ax(j-n)
         else
            con = c(j-nplin)
         end if

         blj = bl(j)
         buj = bu(j)
         if (blj .gt. biglow) blj = blj - con
         if (buj .lt. bigupp) buj = buj - con

         weight = one
         if (j .le. nplin) then
            if (abs(blj) .le. qptol(j)) blj = zero
            if (abs(buj) .le. qptol(j)) buj = zero
         else
            i    = j - nplin
            viol = zero
            if (bl(j) .gt. biglow) then
               if (blj .gt. zero) then
                  viol   = blj
                  if (rho(i) .gt. zero) then
                     weight =   viol*rho(i)
                  else
                     weight =   viol
                  end if
                  wscale = max( wscale,   weight )
                  go to 160
               end if
            end if

            if (bu(j) .lt. bigupp) then
               if (buj .lt. zero) then
                  viol   =   buj
                  if (rho(i) .gt. zero) then
                     weight = - viol*rho(i)
                  else
                     weight = - viol
                  end if
                  wscale = max( wscale, weight )
               end if
            end if

*           Set the vector of nonlinear constraint violations.

  160       violn(i) = viol
         end if

         wtinf(j) = weight
         qpbl(j)  = blj
         qpbu(j)  = buj
  170 continue

      if (wscale .gt. zero) then
         wscale = one/wscale
         call dscal ( nctotl, (wscale), wtinf, 1 )
      end if

      call dcond ( nctotl, wtinf, 1, wtmax, wtmin )
      wtmin  = epspt9*wtmax
      do 180, j = 1, nctotl
         wtinf(j) = max( wtinf(j), wtmin )
  180 continue

*     Set the maximum allowable condition estimator of the constraints
*     in the working set.  Note that a relatively well-conditioned
*     working set is used to start the QP iterations.

      condmx = max( one/epspt3, hundrd )

      if (ncnln .gt. 0) then
*        ===============================================================
*        Refactorize part of the  QP  constraint matrix.
*        ===============================================================
*        Load the new Jacobian into the  QP  matrix  A.  Compute the
*        2-norms of the rows of the Jacobian.

         call f06qff( 'General', ncnln, n, cJac, ldcJ,
     $                Aqp(nclin+1,1), ldAqp )

         do 190 j = nclin+1, ncqp
            Anorm(j) = dnrm2 ( n, Aqp(j,1), ldAqp )
  190    continue

*        Count the number of linear constraints in the working set and
*        move them to the front of kactiv.  Compute the norm of the
*        matrix of constraints in the working set.
*        Let k1  point to the first nonlinear constraint.  Constraints
*        with indices kactiv(k1),..., kactiv(nactiv)  must be
*        refactorized.

         Asize  = zero
         linact = 0
         k1     = nactiv + 1
         do 200 k = 1, nactiv
            i     = kactiv(k)
            Asize = max( Asize, Anorm(i) )

            if (i .le. nclin) then
               linact = linact + 1
               if (linact .ne. k) then
                  iswap  = kactiv(linact)
                  kactiv(linact) = i
                  kactiv(k)      = iswap
               end if
            else

*              Record the old position of the 1st. nonlinear constraint.

               if (k1 .gt. nactiv) k1 = k
            end if
  200    continue

         if (nactiv .le. 1 )
     $      call dcond ( ncqp, Anorm, 1, Asize, Amin )

*        Compute the absolute values of the nonlinear constraints in
*        the working set.  Use DX as workspace.

         do 210 k = linact+1, nactiv
            j = n + kactiv(k)
            if (istate(j) .eq. 1) dx(k) = abs( qpbl(j) )
            if (istate(j) .ge. 2) dx(k) = abs( qpbu(j) )
  210    continue

*        Sort the elements of kactiv corresponding to nonlinear
*        constraints in descending order of violation (i.e.,
*        the first element of kactiv for a nonlinear constraint
*        is associated with the most violated constraint.)
*        In this way, the rows of the Jacobian corresponding
*        to the more violated constraints tend to be included
*        in the  TQ  factorization.

*        The sorting procedure is taken from the simple insertion
*        sort in D. Knuth, ACP Volume 3, Sorting and Searching,
*        Page 81.  It should be replaced by a faster sort if the
*        number of active nonlinear constraints becomes large.

         do 230 k = linact+2, nactiv
            l     = k
            viol  = dx(l)
            kviol = kactiv(l)
*           while (l .gt. linact+1  .and.  dx(l-1) .lt. viol) do
  220       if    (l .gt. linact+1                          ) then
               if (                        dx(l-1) .lt. viol) then
                  dx(l)     = dx(l-1)
                  kactiv(l) = kactiv(l-1)
                  l         = l - 1
                  go to 220
               end if
*           end while
            end if
            dx(l)     = viol
            kactiv(l) = kviol
  230    continue

         k2     = nactiv
         nactiv = k1     - 1
         nZ     = nfree  - nactiv

*        Update the factors  R,  T  and  Q  to include constraints
*        k1  through  k2.

         if (k1 .le. k2)
     $      call lsadds( unitQ, vertex,
     $                   inform, k1, k2, nactiv, nartif, nZ, nfree,
     $                   nrank, nrejtd, nRpq, ngq,
     $                   n, ldQ, ldAqp, ldR, ldT,
     $                   istate, kactiv, kx, condmx,
     $                   Aqp, R, w(lT), w(lRpq), w(lgq), w(lQ),
     $                   w(lwrk1), dx, w(lrlam) )
      end if

*     ==================================================================
*     Solve for dx, the vector of minimum two-norm that satisfies the
*     constraints in the working set.
*     ==================================================================
      call npsetx( unitQ,
     $             ncqp, nactiv, nfree, nZ,
     $             n, nlnx, nctotl, ldQ, ldAqp, ldR, ldT,
     $             istate, kactiv, kx,
     $             dxnorm, gdx,
     $             Aqp, Adx, qpbl, qpbu, w(lRpq), w(lRpq0), dx, w(lgq),
     $             R, w(lT), w(lQ), w(lwrk1) )

*     ==================================================================
*     Solve a quadratic program for the search direction  DX  and
*     multiplier estimates  clamda.
*     ==================================================================
*     If there is no feasible point for the subproblem,  the sum of
*     infeasibilities is minimized subject to the linear constraints
*     (1  thru  jinf)  being satisfied.

      jinf  = n + nclin

      ntry  = 1
*+    repeat
  450    call lscore( 'QP subproblem', qpnamd, names, linobj, unitQ,
     $                nQPerr, minIts, jinf, ncqp, nctotl,
     $                nactiv, nfree, nrank, nZ, nZr,
     $                n, ldAqp, ldR,
     $                istate, kactiv, kx,
     $                gdx, ssq, ssq1, suminf, numinf, dxnorm,
     $                qpbl, qpbu, Aqp, clamda, Adx,
     $                qptol, R, dx, w )

         if (npdbg  .and.  inpdbg(1) .gt. 0)
     $      write(iPrint, 8000) nQPerr

         nviol = 0
         if (numinf .gt. 0) then

*           Count the violated linear constraints.

            do 460 j = 1, nplin
               if (istate(j) .lt. 0) nviol = nviol + 1
  460       continue

            if (nviol .gt. 0) then
               ntry   = ntry + 1
               unitQ  = .true.
               nactiv = 0
               nfree  = n
               nZ     = n
               call iload ( nctotl, (0), istate, 1 )

               call npsetx( unitQ,
     $                      ncqp, nactiv, nfree, nZ,
     $                      n, nlnx, nctotl, ldQ, ldAqp, ldR, ldT,
     $                      istate, kactiv, kx,
     $                      dxnorm, gdx,
     $                      Aqp, Adx, qpbl, qpbu, w(lRpq), w(lRpq0),
     $                      dx, w(lgq), R, w(lT), w(lQ), w(lwrk1) )
            end if
         end if
      if (.not. (nviol .eq. 0  .or.  ntry .gt. 2)) go to 450
*+    until (    nviol .eq. 0  .or.  ntry .gt. 2)

*     ==================================================================
*     Count the number of nonlinear constraint gradients in the  QP
*     working set.  Make sure that all small  QP  multipliers associated
*     with nonlinear inequality constraints have the correct sign.
*     ==================================================================
      nlnact  = 0
      if (nactiv .gt. 0  .and.  ncnln .gt. 0) then
         do 500 k = 1, nactiv
            l     = kactiv(k)
            if (l .gt. nclin) then
               nlnact = nlnact + 1
               j      = n      + l
               if (istate(j) .eq. 1) clamda(j) = max( zero, clamda(j) )
               if (istate(j) .eq. 2) clamda(j) = min( zero, clamda(j) )
            end if
  500    continue
      end if

      linact = nactiv - nlnact

*     ------------------------------------------------------------------
*     Extract various useful quantities from the QP solution.
*     ------------------------------------------------------------------
*     Compute  HPQ = R'R(pq)  from the transformed gradient of the QP
*     objective function and  R(pq)  from the transformed residual.

      call dscal ( n, (-one), w(lRpq), 1 )
      call daxpy ( n, (-one), w(lgq) , 1, w(lHpq), 1 )
      qpcurv = two*ssq

      if (ncnln .gt. 0) then
         if (numinf .gt. 0) then
            feasqp = .false.
            call dload ( nctotl, (zero), clamda, 1 )

            if (nZ .gt. 0) then
*              ---------------------------------------------------------
*              Compute a null space component for the search direction
*              as the solution of  Z'HZ(pz) = -Z'g - Z'HY(py).
*              ---------------------------------------------------------
*              Overwrite DX with the transformed search direction
*              Q'(dx).  The first NZ components of DX are zero.

               call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $                      kx, dx, w(lQ), w(lwrk1) )

*              Overwrite the first NZ components of DX with the solution
*              of  (Rz)u = -(v + w),  where  (Rz)'w = Z'g  and  v  is
*              vector of first NZ components of  R(pq).

               call dcopy ( nZ, w(lgq), 1, dx, 1 )
               call dtrsv ( 'U', 'T', 'N', nZ, R, ldR, dx, 1 )

               call daxpy ( nZ, (one), w(lRpq), 1, dx, 1 )

               call dtrsv ( 'U', 'N', 'N', nZ, R, ldR, dx, 1 )
               call dscal ( nZ, (-one), dx, 1 )

*              Recompute Rpq, Hpq, gdx and qpcurv.

               call dcopy ( nlnx, dx, 1, w(lRpq), 1 )
               call dtrmv ( 'U', 'N', 'N', nlnx, R, ldR, w(lRpq), 1 )
               if (nlnx .lt. n)
     $            call dgemv( 'N', nlnx, n-nlnx, one, R(1,nlnx+1), ldR,
     $                        dx(nlnx+1), 1, one, w(lRpq), 1 )

               gdx    = ddot  ( n, w(lgq) , 1, dx     , 1 )
               qpcurv = ddot  ( n, w(lRpq), 1, w(lRpq), 1 )

               call cmqmul( 3, n, nZ, nfree, ldQ, unitQ,
     $                      kx, dx, w(lQ), w(lwrk1) )

*              ---------------------------------------------------------
*              Recompute Adx and the 2-norm of dx.
*              ---------------------------------------------------------
               dxnorm  = dnrm2 ( n, dx, 1 )
               if (ncqp .gt. 0)
     $            call dgemv ( 'N', ncqp, n, one, Aqp, ldAqp,
     $                         dx, 1, zero, Adx, 1 )

               if (npdbg  .and.  inpdbg(2) .gt. 0)
     $            write(iPrint, 8100) (dx(j), j = 1, n)
            end if

            call dcopy ( nlnx, w(lRpq), 1, w(lHpq), 1 )
            call dtrmv ( 'U', 'T', 'N', nlnx, R, ldR, w(lHpq), 1 )
            if (nlnx .lt. n)
     $         call dgemv ( 'T', nlnx, n-nlnx, one, R(1,nlnx+1), ldR,
     $                      w(lRpq), 1, zero, w(lHpq+nlnx), 1 )
         end if

*        ===============================================================
*        For given values of the objective function and constraints,
*        attempt to minimize the merit function with respect to each
*        slack variable.
*        ===============================================================
         do 600 i = 1, ncnln
            j      = nplin + i
            con    = c(i)

            if (      .not. feasqp  .and.
     $          violn(i) .ne. zero  .and.  rho(i) .le. zero )
     $         rho(i) = one

            quotnt = ddiv  ( cmul(i), scale*rho(i), overfl )

*           Define the slack variable to be  CON - MULT / RHO.
*           Force each slack to lie within its upper and lower bounds.

            if (bl(j) .gt. biglow) then
               if (qpbl(j) .ge. - quotnt) then
                  slk(i) = bl(j)
                  go to 550
               end if
            end if

            if (bu(j) .lt. bigupp) then
               if (qpbu(j) .le. - quotnt) then
                  slk(i) = bu(j)
                  go to 550
               end if
            end if

            slk(i) = con - quotnt

*           The slack has been set within its bounds.

  550       cs(i)  = con - slk(i)

*           ------------------------------------------------------------
*           Compute the search direction for the slacks and multipliers.
*           ------------------------------------------------------------
            dslk(i) = Adx(nclin+i) + cs(i)

            if (feasqp) then
*
*              If any constraint is such that  (DLAM)*(C - S)  is
*              positive,  the merit function may be reduced immediately
*              by substituting the QP multiplier.
*
               dlam(i)  = clamda(j) - cmul(i)
               if (dlam(i) * cs(i) .ge. zero) then
                  cmul(i) = clamda(j)
                  dlam(i) = zero
               end if
            else

*              The  QP  subproblem was infeasible.

               dlam(i) = zero

               if (istate(j) .lt. 0  .or.  violn(i) .ne. zero)
     $            dslk(i)  = zero

            end if
  600    continue

         if (.not. feasqp)
     $      rhonrm = dnrm2 ( ncnln, rho, 1 )

         if (npdbg  .and.  inpdbg(2) .gt. 0) then
            write(iPrint, 8200) (violn(i), i=1,ncnln)
            write(iPrint, 8300) (slk(i)  , i=1,ncnln)
         end if
      end if

      call icopy ( ldbg, inpdbg, 1, icmdbg, 1 )
      idbg   = idbgsv

      return

 1010 format(  1x, 32a4 / ' Start of major itn', i4)
 1020 format(/ 1x, 16a4 / ' Start of major itn', i4)
 8000 format(/ ' //npiqp // nQPerr'
     $       / ' //npiqp // ',  i6 )
 8100 format(/ ' //npiqp // dx recomputed with null space portion...'
     $       / (5g12.3))
 8200 format(/ ' //npiqp // Violations = '/ (1p, 5e15.6))
 8300 format(/ ' //npiqp // Slacks     = '/ (1p, 5e15.6))

*     end of npiqp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npkey ( nout, buffer, key )

      implicit           double precision (a-h,o-z)
      character*(*)      buffer

*     ==================================================================
*     npkey   decodes the option contained in  buffer  in order to set
*     a parameter value in the relevant element of the parameter arrays.
*
*
*     Input:
*
*     nout   A unit number for printing error messages.
*            if nout = 0 no error messages are printed.
*
*     Output:
*
*     key    The first keyword contained in buffer.
*
*
*     npkey  calls opnumb and the subprograms
*                 lookup, scannr, tokens, upcase
*     (now called oplook, opscan, optokn, opuppr)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of npkey  dated 21-Mar-93.
*     ==================================================================
*     +Include lsparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter         (mxparm = 30)
      integer            iprmls(mxparm), ipsvls
      double precision   rprmls(mxparm), rpsvls

      common    /lspar1/ ipsvls(mxparm),
     $                   idbgls, iPrnt , iSumry, itmax1, itmax2, lcrash,
     $	                 ldbgls, lprob , msgls , nn    , nnclin, nprob , 
     $                   ipadls(18)

      common    /lspar2/ rpsvls(mxparm),
     $                   bigbnd, bigdx , bndlow, bndupp, tolact, tolfea,
     $                   tolrnk, rpadls(23)

      equivalence       (iprmls(1), idbgls), (rprmls(1), bigbnd)

      save      /lspar1/, /lspar2/
*     +Include npparm+++++++++++++++++++++++++++++++++++++++++++++++++++
      integer            iprmnp(mxparm), ipsvnp
      double precision   rprmnp(mxparm), rpsvnp

      common    /nppar1/ ipsvnp(mxparm),
     $                   idbgnp, itmxnp, jvrfy1, jvrfy2, jvrfy3, jvrfy4,
     $                   ldbgnp, lformH, lvlder, lverfy, msgNP , nlnf  ,
     $                   nlnj  , nlnx  , nncnln, nsave , nload , ksave ,
     $                   ipadnp(12)

      common    /nppar2/ rpsvnp(mxparm),
     $                   cdint , ctol  , dxlim , epsrf , eta   , fdint ,
     $                   ftol  , hcndbd, rpadnp(22)

      equivalence       (iprmnp(1), idbgnp), (rprmnp(1), cdint)

      save      /nppar1/, /nppar2/
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence  (idbgnp, idbg  ), (itmxnp, nmajor), (itmax2, nminor)
      equivalence  (ldbgls, mnrdbg), (ldbgnp, mjrdbg), (msgls , msgQP )

      external           opnumb
      logical            first , more  , number, opnumb, sorted
      save               first

      parameter         (     maxkey = 43,  maxtie = 22,   maxtok = 10)
      character*16       keys(maxkey), ties(maxtie), token(maxtok)
      character*16       key, key2, key3, value

      parameter         (idummy = -11111,  rdummy = -11111.0d+0,
     $                   sorted = .true.,  zero   =  0.0     )

      data                first
     $                  /.true./
      data   keys
     $ / 'BEGIN           ',
     $   'CENTRAL         ', 'COLD            ', 'CONDITION       ',
     $   'CONSTRAINTS     ',
     $   'CRASH           ', 'DEBUG           ', 'DEFAULTS        ',
     $   'DERIVATIVE      ', 'DIFFERENCE      ', 'END             ',
     $   'FEASIBILITY     ', 'FUNCTION        ', 'HESSIAN         ',
     $   'HOT             ', 'INFINITE        ', 'IPRMLS          ',
     $   'ITERATIONS      ', 'ITERS:ITERATIONS', 'ITNS :ITERATIONS',
     $   'LINEAR          ', 'LINESEARCH      ', 'LIST            ',
     $   'LOAD            ',
     $   'LOWER           ', 'MAJOR           ', 'MINOR           ',
     $   'NOLIST          ', 'NONLINEAR       ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'ROW             ',
     $   'RPRMLS          ', 'SAVE            ', 'START           ',
     $   'STEP            ', 'STOP            ', 'SUMMARY         ',
     $   'UPPER           ', 'VARIABLES       ', 'VERIFY          ',
     $   'WARM            '/

      data   ties
     $ / 'BOUND           ', 'CONSTRAINTS     ', 'DEBUG           ',
     $   'FEASIBILITY     ', 'FILE            ', 'FREQUENCY       ', 
     $   'GRADIENTS       ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'JACOBIAN        ', 'LEVEL           ',
     $   'NO              ', 'NO.      :NUMBER',
     $   'NUMBER          ', 'OBJECTIVE       ', 'PRINT           ',
     $   'RUN             ', 'STEP            ', 'TOLERANCE       ',
     $   'VARIABLES       ', 'YES             '/
*-----------------------------------------------------------------------

      if (first) then
         first  = .false.
         do 10,     i = 1, mxparm
            rprmls(i) = rdummy
            iprmls(i) = idummy
            rprmnp(i) = rdummy
            iprmnp(i) = idummy
   10    continue
      end if

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

*     Convert the keywords to their most fundamental form
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
         else if (key .eq. 'DEBUG       ') then
            idbg   = rvalue
         else if (key .eq. 'DEFAULTS    ') then
            do 20,     i = 1, mxparm
               iprmls(i) = idummy
               rprmls(i) = rdummy
               iprmnp(i) = idummy
               rprmnp(i) = rdummy
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
              if (key2.eq. 'BOUND       ') bigbnd = rvalue * 0.99999
              if (key2.eq. 'STEP        ') bigdx  = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'IPRMLS      ') then
*           Allow things like  iprmls 21 = 100  to set iprmls(21) = 100
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
               read (key3, '(bn, i16)') iprmls(ivalue)
            else
               if (nout .gt. 0) write(nout, 2400) ivalue
            end if
         else if (key .eq. 'ITERATIONS  ') then
            nmajor = rvalue
         else if (key .eq. 'LINEAR      ') then
            if (key2  .eq. 'CONSTRAINTS ') nnclin = rvalue
            if (key2  .eq. 'FEASIBILITY ') tolfea = rvalue
            if (loc2 .eq.  0             ) then
               if (nout .gt. 0)            write(nout, 2320) key2
            end if
         else if (key .eq. 'LINESEARCH  ') then
            eta    = rvalue
         else if (key .eq. 'LOWER       ') then
            bndlow = rvalue
         else if (key .eq. 'LOAD        ') then
            nload  = rvalue
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'MAJOR       ') then
              if (key2.eq. 'DEBUG       ') mjrdbg = rvalue
              if (key2.eq. 'ITERATIONS  ') nmajor = rvalue
              if (key2.eq. 'PRINT       ') msgNP  = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'MINOR       ') then
              if (key2.eq. 'DEBUG       ') mnrdbg = rvalue
              if (key2.eq. 'ITERATIONS  ') nminor = rvalue
              if (key2.eq. 'PRINT       ') msgQP  = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'NONLINEAR   ') then
              if (key2.eq. 'CONSTRAINTS ') nncnln = rvalue
              if (key2.eq. 'FEASIBILITY ') ctol   = rvalue
              if (key2.eq. 'JACOBIAN    ') nlnj   = rvalue
              if (key2.eq. 'OBJECTIVE   ') nlnf   = rvalue
              if (key2.eq. 'VARIABLES   ') nlnx   = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
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
              if (key2.eq. 'FILE        ') iPrnt  = rvalue
              if (key2.eq. 'LEVEL       ') msgNP  = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'PROBLEM     ') then
              if (key2.eq. 'NUMBER      ') nprob  = rvalue
         else if (key .eq. 'ROW         ') then
              if (key2.eq. 'TOLERANCE   ') ctol   = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'RPRMLS      ') then
*           Allow things like  rprmls 21 = 2  to set rprmls(21) = 2.0
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
               read (key3, '(bn, e16.0)') rprmls(ivalue)
            else
               if (nout .gt. 0) write(nout, 2400) ivalue
            end if
         else if (key .eq. 'SAVE        ') then
              if (key2.eq. 'RUN         ') nsave  = rvalue
              if (key2.eq. 'FREQUENCY   ') ksave  = rvalue
         else if (key .eq. 'START       ') then
              if (key2.eq. 'CONSTRAINTS ') jvrfy3 = rvalue
              if (key2.eq. 'OBJECTIVE   ') jvrfy1 = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'STEP        ') then
            dxlim  = rvalue
         else if (key .eq. 'STOP        ') then
              if (key2.eq. 'CONSTRAINTS ') jvrfy4 = rvalue
              if (key2.eq. 'OBJECTIVE   ') jvrfy2 = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'SUMMARY     ') then
            iSumry = rvalue
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
            if (nout .gt. 0) write(nout, 2300) key
         end if
      end if

  900 return

 2300 format(' XXX  Keyword not recognized:         ', a)
 2320 format(' XXX  Second keyword not recognized:  ', a)
 2400 format(' XXX  The parm subscript is out of range:', i10)

*     end of npkey
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nploc ( n, nclin, ncnln, nctotl, litotl, lwtotl)

      implicit           double precision (a-h,o-z)

*     ==================================================================
*     nploc   allocates the addresses of the work arrays for npcore and
*     lscore.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version   14-February-1985.
*     This version of  nploc  dated 12-July-1986.
*     ==================================================================
      common    /sol3cm/ lennam, ldT   , ncolT , ldQ

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      parameter         (lennp = 35)
      common    /sol1np/ locnp(lennp)

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      miniw     = litotl + 1
      minw      = lwtotl + 1

*     Assign array lengths that depend upon the problem dimensions.

      if (nclin + ncnln .eq. 0) then
         lenT     = 0
         lenQ     = 0
      else
         lenT = ldT *ncolT
         lenQ = ldQ*ldQ
      end if

      if (ncnln .eq. 0) then
         lenAqp = 0
      else
         lenAqp = (nclin + ncnln)*n
      end if

      lkactv    = miniw
      lkx       = lkactv + n
      lneedc    = lkx    + n
      liperm    = lneedc + ncnln
      miniw     = liperm + nctotl

      lhfrwd    = minw
      lhctrl    = lhfrwd + n
      lAnorm    = lhctrl + n
      lqpgq     = lAnorm + nclin + ncnln
      lgq       = lqpgq  + n
      lrlam     = lgq    + n
      lT        = lrlam  + n
      lQ        = lT     + lenT
      minw      = lQ     + lenQ

      locls( 1) = lkactv
      locls( 2) = lAnorm
      locls( 8) = lqpgq
      locls( 9) = lgq
      locls(10) = lrlam
      locls(11) = lT
      locls(12) = lQ

*     Assign the addresses for the workspace arrays used by  NPIQP.

      lqpAdx    = minw
      lqpdx     = lqpAdx + nclin + ncnln
      lRpq      = lqpdx  + n
      lRpq0     = lRpq   + n
      lqphz     = lRpq0  + n
      lwtinf    = lqphz  + n
      lwrk1     = lwtinf + nctotl
      lqptol    = lwrk1  + nctotl
      minw      = lqptol + nctotl

      locls( 3) = lqpAdx
      locls( 4) = lqpdx
      locls( 5) = lRpq
      locls( 6) = lRpq0
      locls( 7) = lqphz
      locls(13) = lwtinf
      locls(14) = lwrk1
      locls(15) = lqptol

*     Assign the addresses for arrays used in NPCORE.

      lAqp      = minw
      lAdx      = lAqp   + lenAqp
      lbl       = lAdx   + nclin  + ncnln
      lbu       = lbl    + nctotl
      ldx       = lbu    + nctotl
      lgq1      = ldx    + n
      lfeatl    = lgq1   + n
      lx1       = lfeatl + nctotl
      lwrk2     = lx1    + n
      minw      = lwrk2  + nctotl

      locnp( 1) = lkx
      locnp( 2) = liperm
      locnp( 3) = lAqp
      locnp( 4) = lAdx
      locnp( 5) = lbl
      locnp( 6) = lbu
      locnp( 7) = ldx
      locnp( 8) = lgq1
      locnp(10) = lfeatl
      locnp(11) = lx1
      locnp(12) = lwrk2

      lcs1      = minw
      lcs2      = lcs1   + ncnln
      lc1mul    = lcs2   + ncnln
      lcmul     = lc1mul + ncnln
      lcJdx     = lcmul  + ncnln
      ldlam     = lcJdx  + ncnln
      ldslk     = ldlam  + ncnln
      lrho      = ldslk  + ncnln
      lwrk3     = lrho   + ncnln
      lslk1     = lwrk3  + ncnln
      lslk      = lslk1  + ncnln
      minw      = lslk   + ncnln

      locnp(13) = lcs1
      locnp(14) = lcs2
      locnp(15) = lc1mul
      locnp(16) = lcmul
      locnp(17) = lcJdx
      locnp(18) = ldlam
      locnp(19) = ldslk
      locnp(20) = lrho
      locnp(21) = lwrk3
      locnp(22) = lslk1
      locnp(23) = lslk
      locnp(24) = lneedc

      lcJac     = minw
      lgrad     = lcJac  + ncnln*n
      minw      = lgrad  + n

      locnp(25) = lhfrwd
      locnp(26) = lhctrl
      locnp(27) = lcJac
      locnp(28) = lgrad

      litotl    = miniw - 1
      lwtotl    = minw  - 1

*     end of nploc
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npmrt ( feasqp, n, nclin, ncnln,
     $                   objalf, grdalf, qpcurv,
     $                   istate,
     $                   cJdx, cmul, cs,
     $                   dlam, rho, violn,
     $                   work1, work2 )

      implicit           double precision (a-h,o-z)

      logical            feasqp

      integer            istate(*)

      double precision   cJdx(*), cmul(*), cs(*),
     $                   dlam(*), rho(*), violn(*)
      double precision   work1(*), work2(*)

*     ==================================================================
*     npmrt   computes the value and directional derivative of the
*     augmented Lagrangian merit function.  The penalty parameters
*     rho(j) are boosted if the directional derivative of the resulting
*     augmented Lagrangian function is not sufficiently negative.  If
*     rho needs to be increased,  the perturbation with minimum two-norm
*     is found that gives a directional derivative equal to  - p'Hp.
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version written  27-May-1985.
*     This version of  NPMRT  dated 14-November-1985.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            boost , overfl
      external           ddiv  , ddot  , dnrm2
      intrinsic          abs   , max   , min   , sqrt
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

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1000) qpcurv, grdalf

      if (feasqp) then

*        Find the quantities that define  rhomin, the vector of minimum
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
         if (abs( tscl ) .le. rhomax  .and.  .not. overfl) then
*           ------------------------------------------------------------
*           Bounded  rhomin  found.  The final value of  rho(J)  will
*           never be less than  rhomin(j).  If the  QP  was feasible,  a
*           trial value  rhonew  is computed that is equal to the
*           geometric mean of the previous  rho  and a damped value of
*           rhomin.  The new  rho  is defined as  rhonew  if it is less
*           than half the previous  rho  and greater than  rhomin.
*           ------------------------------------------------------------
            scale  = one
            do 400, i = 1, ncnln
               rhomin = max(  (work2(i)/qnorm)*tscl, zero )
               rhoi   = rho(i)

               rhonew = sqrt( rhoi*(rhodmp + rhomin) )
               if (rhonew .lt. half*rhoi  ) rhoi = rhonew
               if (rhoi   .lt.      rhomin) rhoi = rhomin
               rho(i) = rhoi
  400       continue

            rho1   = rhonrm
            rhonrm = dnrm2 ( ncnln, rho, 1 )

*           ------------------------------------------------------------
*           If  incrun = true,  there has been a run of iterations in
*           which the norm of  rho  has not decreased.  Conversely,
*           incrun = false  implies that there has been a run of
*           iterations in which the norm of rho has not increased.  If
*           incrun changes during this iteration the damping parameter
*           rhodmp is increased by a factor of two.  This ensures that
*           rho(j) will oscillate only a finite number of times.
*           ------------------------------------------------------------
            boost  = .false.
            if (      incrun  .and.  rhonrm .lt. rho1) boost = .true.
            if (.not. incrun  .and.  rhonrm .gt. rho1) boost = .true.
            if (boost) then
               rhodmp = two*rhodmp
               incrun = .not. incrun
            end if
         end if

         if (npdbg  .and.  inpdbg(2) .gt. 0)
     $      write(iPrint, 1200) (rho(l), l=1,ncnln)

      else

*        The  QP  was infeasible.  Do not alter the penalty parameters,
*        but compute the scale factor so that the constraint violations
*        are reduced.

         call ddscl ( ncnln, rho, 1, work1, 1 )
         pterm2 = ddot  ( ncnln, work1, 1, cs, 1 )

         scale  = rhomax
         tscl   = ddiv  ( grdalf, pterm2, overfl )
         if (tscl .gt. scale  .and.  tscl .le. rhomax/(one+rhonrm)
     $                        .and.  .not. overfl)
     $      scale = tscl

         call dcopy ( ncnln, cs, 1, work1, 1 )
      end if

*     ------------------------------------------------------------------
*     Compute the new value and directional derivative of the
*     merit function.
*     ------------------------------------------------------------------
      call ddscl ( ncnln, rho, 1, work1, 1 )

      pterm  = ddot  ( ncnln, work1, 1, cs, 1 )
      objalf = objalf + half*scale*pterm

      if (feasqp)
     $  pterm2 = pterm

      grdalf = grdalf -      scale*pterm2

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1100) scale, rhonrm, grdalf

      return

 1000 format(/ ' //npmrt //        qpcurv        grdalf '
     $       / ' //npmrt //', 1p, 2e14.2 )
 1100 format(/ ' //npmrt //         scale        rhonrm        grdalf '
     $       / ' //npmrt //', 1p, 3e14.2 )
 1200 format(/ ' //npmrt //  Penalty parameters =       '/ (1p, 5e15.6))

*     end of npmrt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npoptn( string )
      character*(*)      string

*     ==================================================================
*     npoptn  loads the option supplied in string into the relevant
*     element of iprmls, rprmls, iprmnp or rprmnp.
*     ==================================================================
      logical             newopt
      common     /sol7np/ newopt
      save       /sol7np/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar
      character*16        key
      character*72        buffer
      logical             first , prnt
      save                first , nout  , prnt
      data                first /.true./

*     If first time in, set nout.
*     newopt is true first time into npfile or npoptn
*     and just after a call to an optimization routine.
*     prnt is set to true whenever newopt is true.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout   =  wmach(11)
      end if
      buffer = string

*     Call npkey to decode the option and set the parameter value.
*     If newopt is true, reset prnt and test specially for nolist.

      if (newopt) then
         newopt = .false.
         prnt   = .true.
         call npkey ( nout, buffer, key )

         if (key .eq. 'NOLIST') then
            prnt   = .false.
         else
            write(nout, '(// a / a /)')
     $         ' Calls to Option Routine',
     $         ' -----------------------'
            write(nout, '( 6x, a )') buffer
         end if
      else
         if (prnt) then
            iPrint = nout
            write(nout, '( 6x, a )') buffer
         else
            iPrint = 0
         end if

         call npkey ( iPrint, buffer, key )
         if (key .eq.   'LIST') prnt = .true.
         if (key .eq. 'NOLIST') prnt = .false.
      end if

*     end of npoptn
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npprt ( KTcond, convrg, MjrMsg, msgNP, msgQP,
     $                   ldR, ldT, n, nclin, ncnln,
     $                   nctotl, nactiv, linact, nlnact, nZ, nfree,
     $                   majit0, majIts, minIts, istate, alfa, nfun,
     $                   condHz, condH, condT, objalf, objf,
     $                   gfnorm, gznorm, cvnorm,
     $                   Ax, c, R, T, violn, x, work )

      implicit           double precision (a-h,o-z)
      character*5        MjrMsg
      logical            KTcond(2), convrg
      integer            istate(nctotl)
      double precision   Ax(*), c(*), R(ldR,*), T(ldT,*), violn(*)
      double precision   x(n)
      double precision   work(n)

*     ==================================================================
*     npprt  prints various levels of output for npcore.
*
*           Msg        Cumulative result
*           ---        -----------------
*
*        le   0        no output.
*
*        eq   1        nothing now (but full output later).
*
*        eq   5        one terse line of output.
*
*        ge  10        same as 5 (but full output later).
*
*        ge  20        objective function,  x,  Ax  and  c.
*
*        ge  30        diagonals of  T  and  R.
*
*     Debug print is performed depending on the logical variable npdbg.
*     npdbg is set true when idbg major iterations have been performed.
*     At this point,  printing is done according to a string of binary
*     digits of the form clsvt (stored in the integer array inpdbg).
*
*     C  set 'on' gives detailed information from the checking routines.
*     L  set 'on' gives information from the linesearch.
*     S  set 'on' gives information from the maximum step routine npalf.
*     V  set 'on' gives various vectors in  npcore  and its auxiliaries.
*     T  set 'on' gives a trace of which routine was called and an
*                 indication of the progress of the run.
*     For example, `Major debug level 11000' gives much output from 
*                  the checking routines and the linesearch.
*
*     Note: Much useless output can be generated by the uninitiated.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written November-1982.
*     This version of  npprt  dated  23-Dec-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm, lines1, lines2

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            first, nlnCon, newSet, prtHdr
      intrinsic          mod
      external           dnrm2

      if (msgNP  .ge. 20  .and.  iPrint .gt. 0) then
         write(iPrint, 1000) majIts
      end if

      if (msgNP  .ge. 5) then

         Mjr    = mod( majIts, 10000 )
         Mnr    = mod( minIts, 10000 )
         neval  = mod( nfun  , 10000 )
         ndf    = mod( nZ    , 10000 )
         nlnCon = ncnln  .gt. 0
         first  = majIts .eq. majIt0

*        ---------------------------------------------------------------
*        If necessary, print a header. 
*        Print a single line of information.
*        ---------------------------------------------------------------
         if (iPrint .gt. 0) then
*           ------------------------------
*           Terse line for the Print file.
*           ------------------------------
            newSet = lines1 .ge. 40
            prtHdr = msgQP  .gt. 0   .or.  first   .or.
     $               msgNP  .ge. 20  .or.  newSet

            if (prtHdr) then
               if (nlnCon) then 
                  write(iPrint, 2000)
               else
                  write(iPrint, 3000)
               end if
               lines1 = 0
            end if

            if (nlnCon) then
               write(iPrint, 2100) Mjr, Mnr, alfa, neval, objalf,
     $                              cvnorm, gznorm, ndf, 
     $                              n-nfree, linact,nlnact,scale*rhonrm,
     $                              gfnorm, condH, condHz, condT,
     $                              convrg, KTcond(1), KTcond(2), 
     $                              MjrMsg
            else
               write(iPrint, 3100) Mjr, Mnr, alfa, neval, objalf,
     $                              gznorm, ndf, n-nfree, linact, 
     $                              gfnorm, condH, condHz, condT,
     $                              convrg, KTcond(1), KTcond(2), 
     $                              MjrMsg
            end if
            lines1 = lines1 + 1
         end if

         if (iSumm .gt. 0) then
*           --------------------------------
*           Terse line for the Summary file.
*           --------------------------------
            newSet = lines2 .ge. 5
            prtHdr = msgQP  .gt. 0  .or.  first 
     $                              .or.  newSet

            if (prtHdr) then
               if (nlnCon) then
                  write(iSumm , 4000)
               else
                  write(iSumm , 5000)
               end if
               lines2 = 0
            end if

            if (nlnCon) then
               write(iSumm , 4100) Mjr, Mnr, alfa, neval, objalf,
     $                              cvnorm, gznorm, ndf, scale*rhonrm,
     $                              convrg, KTcond(1), KTcond(2), MjrMsg
            else
               write(iSumm , 5100) Mjr, Mnr, alfa, neval, objalf,
     $                              gznorm, ndf, 
     $                              convrg, KTcond(1), KTcond(2), MjrMsg
            end if
            lines2 = lines2 + 1
         end if

         if (msgNP .ge. 20  .and.  iPrint .gt. 0) then
            if (ncnln .eq. 0) then
               write(iPrint, 9100) objf
            else
               cviols = dnrm2 ( ncnln, violn, 1 )
               write(iPrint, 9200) objf, cviols
            end if

*           ------------------------------------------------------------
*           Print the constraint values.
*           ------------------------------------------------------------
            write(iPrint, 9300)
            write(iPrint, 9400) (x(j), istate(j), j=1,n)
            if (nclin .gt. 0)
     $         write(iPrint, 9500) (Ax(k), istate(n+k),      k=1,nclin)
            if (ncnln .gt. 0)
     $         write(iPrint, 9600) (c(k) , istate(n+nclin+k),k=1,ncnln)

            if (msgNP .ge. 30) then
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               incT   = ldT - 1
               if (nactiv .gt. 0) then
                  call dcopy( nactiv, T(nactiv,nZ+1), incT, work, 1 )
                  write(iPrint, 9700) (work(j), j=1,nactiv)
               end if
               write(iPrint, 9800) (R(j,j), j=1,n)
            end if
         end if
      end if

      if (msgNP .ge. 20  .and.  iPrint .gt. 0) then
         write(iPrint, 9900)
      end if

      MjrMsg(1:2) = '  '
      MjrMsg(4:5) = '  '

      return

 1000 format(/// ' Major iteration', i5
     $       /   ' ====================' )
 2000 format(//  '  Maj  Mnr    Step  Fun  Merit function',
     $           '  Violtn Norm gZ   nZ  Bnd  Lin  Nln', 
     $           ' Penalty Norm Gf  Cond H Cond Hz',
     $           '  Cond T Conv' )
 2100 format(    2i5, 1p, e8.1, i5, e16.8, 2e8.1, 4i5, 2e8.1, 3e8.0,
     $           1x, l1, 1x, 2l1, a5 )
 3000 format(//  '  Maj  Mnr    Step  Fun       Objective',
     $                       ' Norm gZ   nZ  Bnd  Lin',
     $                       ' Norm Gf  Cond H Cond Hz', 
     $           '  Cond T Conv' )
 3100 format(    2i5, 1p, e8.1, i5, e16.8,  e8.1, 3i5,  e8.1, 3e8.0,
     $           1x, l1, 1x, 2l1, a5 )
 4000 format(//  '  Maj  Mnr    Step  Fun  Merit function',
     $           '  Violtn Norm gZ   nZ Penalty Conv' )
 4100 format(    2i5, 1p, e8.1, i5, e16.8, 2e8.1, i5, e8.1, 
     $           1x, l1, 1x, 2l1, a5 )
 5000 format(//  '  Maj  Mnr    Step  Fun       Objective',
     $           ' Norm gZ   nZ Conv' )
 5100 format(    2i5, 1p, e8.1, i5, e16.8,  e8.1, i5, 
     $           1x, l1, 1x, 2l1, a5 )
 9100 format(/ ' Nonlinear objective value = ', 1p, e15.6 )
 9200 format(/ ' Nonlinear objective value = ', 1p, e15.6, '   Norm of',
     $         ' the nonlinear constraint violations = ', e15.6 )
 9300 format(/ ' Values of the constraints and their predicted status'
     $       / ' ----------------------------------------------------')
 9400 format(/ ' Variables                  '/ (1x, 5(1p, e15.6, i4)))
 9500 format(/ ' General linear constraints '/ (1x, 5(1p, e15.6, i4)))
 9600 format(/ ' Nonlinear constraints      '/ (1x, 5(1p, e15.6, i4)))
 9700 format(/ ' Diagonals of  T  =         '
     $       /   (1p, 5e15.6))
 9800 format(/ ' Diagonals of  R  =         '
     $       /   (1p, 5e15.6))
 9900 format(  ' ==================================================',
     $         '======================================='///)

*     end of npprt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nprset( unitQ,
     $                   n, nfree, nZ, ldQ, ldR,
     $                   iperm, kx,
     $                   gq, R, Q, work, qrwork )

      implicit           double precision (a-h,o-z)
      logical            unitQ
      integer            iperm(n), kx(n)
      double precision   gq(n), R(ldR,*), Q(ldQ,*)
      double precision   work(n), qrwork(2*n)

*     ==================================================================
*     nprset  bounds the condition estimator of the transformed Hessian.
*     On exit, R is of the form
*                  ( DRz   0     )
*                  (  0  sigma*I )
*     where D is a diagonal matrix such that DRz has a bounded condition
*     number,  I is the identity matrix and sigma  is the geometric mean
*     of the largest and smallest elements of DRz. The QR factorization
*     with interchanges is used to give diagonals of DRz that are
*     decreasing in modulus.
*
*     Systems Optimization Laboratory, Stanford University.
*
*     Original version of nprset dated  4-August-1986.
*     This version dated  14-Sep-92.
*     ==================================================================
      common    /sol6cm/ Rcndbd, Rfrobn, dRmax, dRmin

      logical            npdbg
      parameter         (ldbg   = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      intrinsic          max   , min   , log   , real  , sqrt
      external           dnorm , idrank
      parameter        ( zero   =0.0d+0, half =0.5d+0, one    =1.0d+0 )

*     ==================================================================
*     Bound the condition estimator of Q'HQ.
*     The scheme used here reduces the modulus of the larger
*     diagonals while increasing the modulus of the smaller ones.
*     ==================================================================
      if (nZ .gt. 1) then
*        ---------------------------------------------------------------
*        Refactorize Rz.  Interchanges are used to give diagonals
*        of decreasing magnitude.
*        ---------------------------------------------------------------
         do 100, j = 1, nZ-1
            call dload ( nZ-j, zero, R(j+1,j), 1 )
  100    continue

         call dgeqrp( 'Column iterchanges', nZ, nZ, R, ldR,
     $                work, iperm, qrwork, info )

         do 110, j = 1, nZ
            jmax = iperm(j)
            if (jmax .gt. j) then
               if (unitQ) then
                  jsave    = kx(jmax)
                  kx(jmax) = kx(j)
                  kx(j)    = jsave
               else
                  call dswap ( nfree, Q(1,jmax), 1, Q(1,j), 1 )
               end if

               gjmax    = gq(jmax)
               gq(jmax) = gq(j)
               gq(j)    = gjmax
            end if
  110    continue
      end if

      drgm  = one

      if (nZ .gt. 0) then
         nrank = idrank( nZ, R, ldR+1, one/Rcndbd )
         drgm  = half*sqrt(abs( R(1,1)*R(nrank,nrank) ))
         drgs  = abs( R(1,1) ) / Rcndbd

         if (nZ .gt. nrank) then
            do 120, j = nrank+1, nZ
               call dload ( j-1, zero, R(1,j), 1 )
  120       continue
            call dload ( nZ-nrank, drgs, R(nrank+1,nrank+1), ldR+1 )
         end if
      end if

*     ------------------------------------------------------------------
*     Reset the range-space partition of the Hessian.
*     ------------------------------------------------------------------
      if (nZ .lt. n) then
         do 130, j = nZ+1, n
            call dload ( j, zero, R(1,j), 1 )
  130    continue
         call dload ( n- nZ, drgm, R(nZ+1,nZ+1), ldR+1 )
      end if

*     Recompute the Frobenius norm of R.

      scle  = sqrt(real(n - nZ))*drgm
      sumsq = one
      do 140, j = 1, nZ
         call dssq  ( j, R(1,j), 1, scle, sumsq )
  140 continue
      Rfrobn = dnorm( scle, sumsq )

      return

*     end of nprset
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsavr( unitQ, n, nclin, ncnln, ldR, ldQ,
     $                   nfree, nsave, iter, istate, kx,
     $                   hforwd, hcntrl,
     $                   cmul, R, rho, x, Q )

      implicit           double precision (a-h,o-z)
      logical            unitQ
      integer            istate(n+nclin+ncnln), kx(n)
      double precision   R(ldR,*), x(n), Q(ldQ,*)
      double precision   hforwd(*), hcntrl(*)
      double precision   cmul(*), rho(*)

*     ==================================================================
*     npsavr  saves the details of this run on unit nsave.
*
*     Original version   24-Nov-89.
*     This version of  npsavr  dated  1-Dec-89.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      if (nsave .le. 0) return

*     We always save the computed difference intervals, if defined.

      if (lvldif .gt. 0  .and.  lfdset .eq. 0) then
         lfdnew = 2
      else 
         lfdnew = lfdset
      end if
 
      write(nsave, 1000) iter, nfree, lfdnew, lvldif, unitQ
      do 110, j = 1, n
         write(nsave, 1010) j, kx(j), istate(j), x(j)
  110 continue

      do 120, j = n+1, n+nclin
         write(nsave, 1020) j, istate(j)
  120 continue

      if (ncnln .gt. 0) then
         k = 1
         do 130, j = n+nclin+1, n+nclin+ncnln
            write(nsave, 1030) j, istate(j), cmul(k), rho(k)
            k = k + 1
  130    continue

         write(nsave, 1040) rhomax, rhonrm, rhodmp, scale, incrun
      end if

*     ------------------------------------------------------------------
*     Write  Q(free)  and the factor of  Q'HQ.
*     ------------------------------------------------------------------
      if (.not. unitQ) then
         do 210, j = 1, nfree
            do 200, i = 1, nfree
              write(nsave, 2000) i, j, Q(i,j)
  200       continue
  210    continue
      end if

      do 240, j = 1, n
         do 230, i = 1, j
           write(nsave, 2000) i, j, R(i,j)
  230    continue
  240 continue

*     ------------------------------------------------------------------
*     Read the finite-difference intervals.  
*     ------------------------------------------------------------------
      if (lvldif .gt. 0) then
         if (lfdnew .ne. 1) then
            do 250, j = 1, n
               write(nsave, 1050) j, hforwd(j), hcntrl(j)
  250       continue
         end if
      end if                                            
             
      if (iPrint .gt.      0) write(iPrint, 4000) nsave, iter
      if (nsave  .ne. iPrint) rewind nsave
      return

 1000 format(2i8, 1x, 2i2, 1x, l1 )
 1010 format(2i8, 1x,  i2, 1p, 2e24.14 )
 1020 format( i8, 1x,  i2 )
 1030 format( i8, 1x,  i2, 1p, 2e24.14 )
 1040 format( 1p, 4e24.14, 1x, l1 )
 1050 format( i8, 1x, 1p, 2e24.14 )
 2000 format(2i8, 1x, 1p,  e24.14 )
 4000 format(  ' Current run saved on file', I4, '    Itn = ', I7 )

*     end of npsavr
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsetx( unitQ,
     $                   ncqp, nactiv, nfree, nZ,
     $                   n, nlnx, nctotl, ldQ, ldAqp, ldR, ldT,
     $                   istate, kactiv, kx,
     $                   dxnorm, gdx,
     $                   Aqp, Adx, bl, bu, Rpq, Rpq0, dx, gq,
     $                   R, T, Q, work )

      implicit           double precision (a-h,o-z)
      logical            unitQ
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   Aqp(ldAqp,*), Adx(*), bl(nctotl), bu(nctotl),
     $                   Rpq(nlnx), Rpq0(nlnx), gq(n), R(ldR,*),
     $                   T(ldT,*), Q(ldQ,*), dx(n), work(n)

*     ==================================================================
*     npsetx   defines a point which lies on the initial working set for
*     the QP subproblem.  This routine is a similar to LSSETX except
*     that advantage is taken of the fact that the initial estimate of
*     the solution of the least-squares subproblem is zero.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 BLAS added 12-June-1986.
*     This version of npsetx dated 11-June-1986.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            npdbg
      parameter         (ldbg = 5)
      common    /npdebg/ inpdbg(ldbg), npdbg

      external           ddot, dnrm2
      intrinsic          abs , min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      nfixed = n - nfree

      gdx    = zero
      call dload ( n   , zero, dx  , 1 )
      call dload ( nlnx, zero, Rpq , 1 )
      call dload ( nlnx, zero, Rpq0, 1 )

      if (nactiv + nfixed .gt. 0) then

*        Set  work = residuals for constraints in the working set.
*        Solve for  dx,  the smallest correction to  x  that gives a
*        point on the constraints in the working set.
*        Set the fixed variables on their bounds,  solve the triangular
*        system  T*(dxy) = residuals,  and define  dx = Y*(dxy).
*        Use  (dxy)  to update  d(=Pr)  as  d = d - R'(  0  ).
*                                                     ( dxy )

         do 100, i = 1, nfixed
            j   = kx(nfree+i)
            if (istate(j) .le. 3) then
               bnd   = bl(j)
               if (istate(j) .eq. 2) bnd = bu(j)
               dx(j) = bnd
               work(nfree+i) = bnd
            else
               work(nfree+i) = zero
            end if
  100    continue

         do 110, i = 1, nactiv
            k   = kactiv(i)
            j   = n + k
            bnd = bl(j)
            if (istate(j) .eq. 2) bnd = bu(j)
            work(nZ+i) = bnd - ddot  ( n, Aqp(k,1), ldAqp, dx, 1 )
  110    continue

         if (nactiv .gt. 0)
     $      call cmtsol( 1, ldT, nactiv, T(1,nZ+1), work(nZ+1) )
         call dcopy ( nactiv+nfixed, work(nZ+1), 1, dx(nZ+1), 1 )
         if (nZ .gt. 0)
     $      call dload ( nZ, zero, dx, 1 )

         gdx  = ddot  ( nactiv+nfixed, gq(nZ+1), 1, dx(nZ+1), 1 )

         if (nZ .lt. n) then
            call dgemv ('N', nZ, n-nZ, -one, R(1,nZ+1), ldR,
     $                  dx(nZ+1), 1, one, Rpq, 1 )
            if (nZ .lt. nlnx) then
               nr  = ldR
               if (nZ+1 .eq. n) nr = 1
               call dcopy ( nlnx-nZ, dx(nZ+1), 1, Rpq(nZ+1), 1 )
               call dscal ( nlnx-nZ, (-one),      Rpq(nZ+1), 1 )
               call dtrmv ( 'U', 'N', 'N', nlnx-nZ, R(nZ+1,nZ+1), nr,
     $                      Rpq(nZ+1), 1 )
               if (nlnx .lt. n) then
                  nR = ldR
                  if (nlnx+1 .eq. n) nR = n - nZ
                  call dgemv( 'N', nlnx-nZ, n-nlnx, -one,R(nZ+1,nlnx+1),
     $                        nR, dx(nlnx+1), 1, one, Rpq(nZ+1), 1 )
               end if
            end if
         end if

         call cmqmul( 2, n, nZ, nfree, ldQ, unitQ, kx, dx, Q, work )
      end if

*     ------------------------------------------------------------------
*     Compute the 2-norm of  dx.
*     Initialize  A*dx.
*     ------------------------------------------------------------------
      dxnorm  = dnrm2 ( n, dx, 1 )
      if (ncqp .gt. 0)
     $   call dgemv ( 'N', ncqp, n, one, Aqp, ldAqp, dx, 1, zero,Adx,1)

      if (npdbg  .and.  inpdbg(2) .gt. 0)
     $   write(iPrint, 1200) (dx(j), j = 1, n)

      return

 1200 format(/ ' //npsetx// Variables after npsetx ... '/ (5g12.3))

*     end of npsetx
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npsrch( needfd, inform, n, ncnln,
     $                   ldcJ, ldcJu, nfun, ngrad,
     $                   needc, funcon, funobj,
     $                   alfa, alfbnd, alfmax, alfsml, dxnorm,
     $                   epsrf, eta, gdx, grdalf, glf1, glf,
     $                   objf, objalf, qpcurv, xnorm,
     $                   c, c2, cJac, cJacu, cJdx, cJdx2,
     $                   cmul1, cmul, cs1, cs, dx, dlam, dslk,
     $                   grad, gradu, qpmul, rho, slk1, slk, x1, x,
     $                   work, w, lenw )

      implicit           double precision (a-h,o-z)
      logical            needfd
      integer            needc(*)
      double precision   dx(n), grad(n), gradu(n), x1(n), x(n)
      double precision   c(*), c2(*), cJac(ldcJ,*), cJacu(ldcJu,*),
     $                   cJdx(*), cJdx2(*), cmul1(*), cmul(*)
      double precision   cs1(*), cs(*), dlam(*), dslk(*), qpmul(*),
     $                   rho(*), slk1(*), slk(*),  work(*)
      double precision   w(lenw)
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
*     Level 2 BLAS added 12-June-1986.
*     This version of npsrch dated  14-Sep-92.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol4np/ lvldif, ncdiff, nfdiff, lfdset
      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter        ( ldbg = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            debug , done  , first , imprvd
      external           ddot
      intrinsic          abs   , max   , min   , sqrt
      parameter        ( zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0 )
      parameter        ( two    =2.0d+0                                )
      parameter        ( tolg   =1.0d-1                                )
      parameter        ( rmu    =1.0d-4                                )

      epsmch = wmach(3)

      if (.not. needfd  .and.  ncnln .gt. 0)
     $   cs1jdx = ddot( ncnln, cs1, 1, cJdx, 1 )

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
*             if  M(tolabs) - M(0) .gt. epsaf,  the linesearch tries
*             steps in the range  toltny .le. alfa .le. tolabs.
*     ------------------------------------------------------------------
      nstate = 0
      debug  = npdbg  .and.  inpdbg(4) .gt. 0

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

      if (ncnln .gt. 0) call iload ( ncnln, (1), needc, 1 )

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
               gdx    = tgdx
               glf    = tglf

               if (ncnln .gt. 0) then
                  call dcopy ( ncnln, cJdx2, 1, cJdx, 1 )
                  call f06qff( 'General', ncnln, n, cJacu, ldcJu,
     $                         cJac, ldcJ )
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
               call funcon( mode, ncnln, n, ldcJu,
     $                      needc, x, c2, cJacu, nstate )
               if (mode .lt. 0) go to 999

               call dcopy ( ncnln,         c2 , 1, cs, 1 )
               call daxpy ( ncnln, (-one), slk, 1, cs, 1 )

               call dcopy ( ncnln, cs , 1,  work, 1 )
               call ddscl ( ncnln, rho, 1,  work, 1 )

               fterm  =            ddot( ncnln, cmul, 1, cs, 1 ) -
     $                  half*scale*ddot( ncnln, work, 1, cs, 1 )
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
               gtry   = ddot( n, gradu, 1, dx, 1 )
               tgdx   = gtry
               tglf   = gtry
               if (ncnln .gt. 0) then

*                 Compute the Jacobian times the search direction.

                  call dgemv ( 'n', ncnln, n, one, cJacu, ldcJu, dx, 1,
     $                         zero, cJdx2, 1 )

                  call dcopy ( ncnln,         cJdx2, 1, work, 1 )
                  call daxpy ( ncnln, (-one), dslk , 1, work, 1 )

                  gtry   = gtry - ddot( ncnln, cmul, 1, work, 1 )
                  if (alfa .le. one)
     $               gtry   = gtry - ddot( ncnln, dlam, 1, cs      , 1 )

                  call ddscl ( ncnln, rho , 1, work, 1 )
                  gtry = gtry  +
     $                     scale*ddot( ncnln, work, 1, cs   , 1 )

                  tglf = tgdx  - ddot( ncnln, cJdx2, 1, qpmul, 1 )

*                 ------------------------------------------------------
*                 If alfbnd .le. alfa .lt. alfmax and the norm of the
*                 quasi-Newton update is bounded, set alfmax to be alfa.
*                 This will cause the linesearch to stop if the merit
*                 function is decreasing at the boundary.
*                 ------------------------------------------------------
                  if (alfbnd .le. alfa  .and.  alfa .lt. alfmax) then
                     csjdx  = ddot   ( ncnln, cs, 1, cJdx2, 1 )
                     curvlf = tglf  - glf1
                     curvc  = abs( csjdx - cs1jdx )
                     rhobfs = max( qpcurv*tolg - curvlf, zero )
                     if (rhobfs .le. curvc*rhomax) then
                        alfmax = alfa
                     else
                        alfbnd = min( two*alfa, alfmax )
                     end if

                     if (npdbg  .and.  inpdbg(1) .gt. 0) then
                        write(iPrint, 1400) csjdx , cs1jdx, curvlf
                        write(iPrint, 1300) alfbnd, alfa  , alfmax
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

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1200) inform

      return

*     The user wants to stop.  Who am I to object?

  999 inform = mode
      return

 1200 format(/ ' //npsrch// inform  = ', i4 )
 1300 format(/ ' //npsrch//        alfbnd          alfa        alfmax'
     $       / ' //npsrch//', 1p, 3e14.2 )
 1400 format(/ ' //npsrch//         csjdx        cs1jdx        curvlf'
     $       / ' //npsrch//', 1p, 3e14.2 )

*     end of npsrch
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine npupdt( MjrMsg, unitQ,
     $                   n, ncnln, nfree, nZ,
     $                   ldcJ1, ldcJ2, ldQ, ldR, kx,
     $                   alfa, glf1, glf2, qpcurv,
     $                   cJac1, cJac2, cJdx1, cJdx2,
     $                   cs1, cs2, gq1, gq2, Hpq, Rpq,
     $                   qpmul, R, omega, Q, wrk1, wrk2 )

      implicit           double precision (a-h,o-z)
      character*5        MjrMsg
      logical            unitQ
      integer            kx(n)
      double precision   cJac1(ldcJ1,*), cJac2(ldcJ2,*),
     $                   cJdx1(*), cJdx2(*), cs1(*), cs2(*),
     $                   gq1(n), gq2(n), Hpq(n), Rpq(n), qpmul(*),
     $                   R(ldR,*), omega(*), Q(ldQ,*)
      double precision   wrk1(n+ncnln), wrk2(n)

*     ==================================================================
*     npupdt  computes the BFGS update for the approximate Hessian of
*     the Lagrangian.  If the approximate curvature of the Lagrangian
*     function is negative,  a nonnegative penalty vector OMEGA(i) of
*     minimum two norm is computed such that the approximate curvature
*     of the augmented Lagrangian will be positive. If no finite penalty
*     vector exists,  the BFGS update is performed with the approximate
*     curvature modified to be a small positive value.
*
*     On entry,  gq1 and gq2 contain the transformed objective gradients
*     at x1 and x2,  Hpq contains  R'R(pq), the transformed Hessian
*     times the transformed search direction.  The vectors gq1 and Hpq
*     are not saved.  If the regular BFGS quasi-Newton update could not
*     be performed, the first character of MjrMsg is loaded with 'M'.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 66 version written April 1984.
*     Level 2 BLAS added 12-June-1986.
*     Level 2 matrix routines added 22-Apr-1988.
*     This version of npuptd dated  22-Apr-1988.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol6cm/ Rcndbd, Rfrobn, dRmax, dRmin

      logical            incrun
      common    /sol6np/ rhomax, rhonrm, rhodmp, scale, incrun

      logical            npdbg
      parameter        ( ldbg   = 5 )
      common    /npdebg/ inpdbg(ldbg), npdbg

      logical            overfl, ssbfgs
      intrinsic          max   , min   , sqrt
      external           idamax, ddiv  , dnrm2
      parameter        ( zero   = 0.0d+0, one    = 1.0d+0 )
      parameter        ( tolg   = 1.0d-1                  )

      if (ncnln .gt. 0) call dload ( ncnln, zero, omega, 1 )

*     ------------------------------------------------------------------
*     Set curvl = (g2 - g1)'dx,  the approximate curvature along dx of
*     the (augmented) Lagrangian.  At first, the curvature is not scaled
*     by the steplength alfa.
*     ------------------------------------------------------------------
      curvl  = glf2 -   glf1
      tinycl =        qpcurv * tolg
      ssbfgs = curvl .le. alfa*tinycl
      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1000) ssbfgs, tinycl, curvl

*     ------------------------------------------------------------------
*     Test if curvl is sufficiently positive.  If there are no nonlinear
*     constraints,  no update can be performed.
*     ------------------------------------------------------------------
      if (curvl  .lt. tinycl) then
         MjrMsg(1:1) = 'modified BFGS'
         if (ncnln .gt. 0) then
            qmax = zero
            do 200 i = 1, ncnln
               qi    = cJdx2(i)*cs2(i) - cJdx1(i)*cs1(i)
               qmax  = max( qmax, qi )
               if (qi .le. zero) then
                  wrk1(i) = zero
               else
                  wrk1(i) = qi
               end if
  200       continue

            qnorm = dnrm2 ( ncnln, wrk1, 1 )

            test  = max( tinycl - curvl, zero )
            beta  = ddiv  ( qmax*test, qnorm*qnorm, overfl )
            if (beta .lt. rhomax  .and.  .not. overfl) then
               MjrMsg(1:1) = ' '
               beta  = test/(qnorm*qnorm)
               do 210 i = 1, ncnln
                  qi       = wrk1(i)
                  omega(i) =            beta*qi
                  curvl    = curvl    + beta*qi*qi
  210          continue

               if (npdbg) then
                  imax = idamax( ncnln, omega, 1 )
                  if (inpdbg(1) .gt. 0)
     $               write(iPrint, 1250) omega(imax)

                  if (inpdbg(2) .gt. 0)
     $               write(iPrint, 1300) (omega(i), i=1,ncnln)
               end if
            end if
         end if
      end if

*     ------------------------------------------------------------------
*     Compute the difference in the augmented Lagrangian gradient.
*     ------------------------------------------------------------------
*     Update gq1 to include the augmented Lagrangian terms.

      if (ncnln .gt. 0) then

         do 310 i = 1, ncnln
            wrk1(i) = - qpmul(i) + omega(i) * cs1(i)
  310    continue
         call dgemv ( 'T', ncnln, n, one, cJac1, ldcJ1, wrk1, 1,
     $                zero, wrk2, 1 )

         do 320 i = 1, ncnln
            wrk1(i) =   qpmul(i) - omega(i) * cs2(i)
  320    continue
         call dgemv ( 'T', ncnln, n, one, cJac2, ldcJ2, wrk1, 1,
     $                one, wrk2, 1 )

         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ, kx, wrk2, Q, wrk1 )
         call daxpy ( n, one, wrk2, 1, gq1, 1 )
      end if

      if (npdbg  .and.  inpdbg(1) .gt. 0)
     $   write(iPrint, 1100) alfa  , curvl

      if (curvl .lt. tinycl) curvl  = tinycl

      do 330 j = 1, n
         wrk2(j) = gq2(j) - gq1(j)
  330 continue

      rtgtp  = sqrt(qpcurv)
      rtyts  = sqrt(alfa*curvl)
      eta    = one
      if (ssbfgs)
     $   eta = rtyts / (rtgtp*alfa)

      trace1 = dnrm2 ( n,  Hpq, 1 ) /  rtgtp
      trace2 = dnrm2 ( n, wrk2, 1 ) / (rtyts*eta)
      Rfrobn = eta*sqrt( abs(  (Rfrobn - trace1)*(Rfrobn + trace1)
     $                                 + trace2**2) )

*     ==================================================================
*     Update the Cholesky factor of  Q'HQ.
*     ==================================================================
*     Normalize the vector  RPQ ( = R(pq) ).

      call dscal ( n, (one / rtgtp), Rpq, 1 )

*     Do the self-scaled or regular BFGS update.
*     Form the vector wrk1 = gamma * (gq2 - gq1) - beta * R'R*pq,
*     where  gamma = 1/sqrt( curv ) = 1/sqrt( (gq2 - gq1)'sq )

      call dscal ( n, (one / rtgtp), Hpq, 1 )

      if (ssbfgs) then
         do 410 j   = 1, n
            call dscal ( j, eta, R(1,j), 1 )
            wrk1(j) = wrk2(j)/rtyts  -  eta * Hpq(j)
  410    continue
      else
         do 420 j   = 1, n
            wrk1(j) = wrk2(j)/rtyts  -        Hpq(j)
  420    continue
      end if

*     Perform the update to  R = R + Rpq*wrk1'.
*     Rpq is overwritten. Arrays  gq1  and  Hpq  are used to store the
*     sines and cosines defined by the plane rotations.

      call cmr1md( n, 0, n, ldR, n, n, R, Hpq, Rpq, wrk1, gq1, Hpq )

      return

 1000 format(/ ' //npupdt// ssbfgs    min. curvl         curvl '
     $       / ' //npupdt//   ', l4, 1p, 2e14.2 )
 1100 format(/ ' //npupdt//          alfa         curvl '
     $       / ' //npupdt//', 1p, 2e14.2 )
 1250 format(/ ' //npupdt//   omega(imax)'
     $       / ' //npupdt//', 1pe14.2 )
 1300 format(/ ' //npupdt//  Penalty parameters = '  / (1p, 5e15.6))

*     end of npupdt
      end

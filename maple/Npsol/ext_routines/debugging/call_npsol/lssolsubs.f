*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*     file  lssubs.f
*
*     lssol    lsadd    lsadds   lsbnds   lschol   lscore   lscrsh
*     lsdel    lsdflt   lsfeas   lsfile   lsgetp   lsgset   lskey
*     lsloc    lsmove   lsmuls   lsoptn   lsprt    lssetx
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lssol ( mm, n,
     $                   nclin, ldA, ldR,
     $                   A, bl, bu, cvec,
     $                   istate, kx, x, R, b,
     $                   inform, iter, obj, clamda,
     $                   iw, leniw, w, lenw )

      implicit           double precision(a-h,o-z)
      integer            leniw, lenw
      integer            istate(n+nclin), kx(n)
      integer            iw(leniw)
      double precision   bl(n+nclin), bu(n+nclin), A(ldA,*)
      double precision   clamda(n+nclin), cvec(*)
      double precision   R(ldR,*), x(n), b(*)
      double precision   w(lenw)

*     ==================================================================
*     lssol  solves problems of the form
*
*           Minimize               f(x)
*              x
*                                 (  x )
*           subject to    bl  .le.(    ).ge.  bu,
*                                 ( Ax )
*
*     where  '  denotes the transpose of a column vector,  x  denotes
*     the n-vector of parameters and  f(x) is one of the following
*     functions...
*
*    FP = None                        (find a feasible point).
*    LP = c'x
*    QP1=       1/2 x'Rx               R  n times n, symmetric pos. def.
*    QP2= c'x + 1/2 x'Rx               .  .   ..        ..       ..  ..
*    QP3=       1/2 x'R'Rx             R  m times n, upper triangular.
*    QP4= c'x + 1/2 x'R'Rx             .  .   ..  .   ..      ...
*    LS1=       1/2 (b - Rx)'(b - Rx)  R  m times n, rectangular.
*    LS2= c'x + 1/2 (b - Rx)'(b - Rx)  .  .   ..  .     ...
*    LS3=       1/2 (b - Rx)'(b - Rx)  R  m times n, upper triangular.
*    LS4= c'x + 1/2 (b - Rx)'(b - Rx)  .  .   ..  .   ..      ...
*
*     The matrix  R  is entered as the two-dimensional array  R  (of
*     row dimension  ldR).  If  ldR = 0,  R  is not accessed.
*
*     The vector  c  is entered in the one-dimensional array  cvec.
*
*     nclin  is the number of general linear constraints (rows of  A).
*     (nclin may be zero.)
*
*     The first  n  components of  bl  and   bu  are lower and upper
*      bounds on the variables.  The next  nclin  components are
*     lower and upper bounds on the general linear constraints.
*
*     The matrix  A  of coefficients in the general linear constraints
*     is entered as the two-dimensional array  A  (of dimension
*     ldA by n).  If nclin = 0, A is not accessed.
*
*     The vector  x  must contain an initial estimate of the solution,
*     and will contain the computed solution on output.
*
*
*     Complete documentation for LSSOL is contained in Report SOL 86-1,
*     Users Guide for LSSOL (Version 1.0), by P.E. Gill,
*     S.J. Hammarling, W. Murray, M.A. Saunders and M.H. Wright,
*     Department of Operations Research, Stanford University, Stanford, 
*     California 94305.
*
*     Systems Optimization Laboratory, Stanford University.
*     Version 1.00 Dated  30-Jan-1986.
*     Version 1.01 Dated  30-Jun-1986.   Level-2 BLAS added
*     Version 1.02 Dated  13-May-1988.   Level-2 matrix routines added.
*     Version 1.03 Dated  19-Jun-1989.   Some obscure bugs fixed.      
*     Version 1.04 Dated  26-Aug-1991.   nrank bug fixed.      
*     Version 1.05 Dated  20-Sep-1992.   Output modified. 
*                         20-Oct-1992    Summary file included.
*                         
*     Copyright  1984/1993  Stanford University.
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
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT, ncolT, ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize, dTmax, dTmin

      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
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
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      intrinsic          max, min

*     Local variables.

      logical            cold  , factrz, linObj, named , rowerr,
     $                   unitQ , vertex
      character*2        prbtyp
      character*8        names(1)
      parameter         (zero   =0.0d+0, point3 =3.3d-1, point8 =0.8d+0)
      parameter         (point9 =0.9d+0, one    =1.0d+0, hundrd =1.0d+2)

      character*40       title
      data               title
     $                 / 'SOL/LSSOL  ---  Version 1.05-2  April 93' /
*                         123456789|123456789|123456789|123456789|
*
*     Set the machine-dependent constants.

      call mchpar()

      epsmch = wmach( 3)
      rteps  = wmach( 4)
      nout   = wmach(11)

      epspt3 = epsmch**point3
      epspt5 = rteps
      epspt8 = epsmch**point8
      epspt9 = epsmch**point9

      named  = .false.

      inform = 0
      iter   = 0

C-->  condmx = max( one/epspt5, hundrd )
C-->  condmx = max( one/epspt3, hundrd )
      condmx = max( one/epspt5, hundrd )

      nctotl = n + nclin

*     Set the default values of the parameters.

      call lsdflt( mm, n, nclin, title )

*     Set all parameters determined by the problem type.

      if      (lprob .eq. 1 ) then
         prbtyp    = 'FP'
         m      = 0
         linObj = .false.
         factrz = .true.
      else if (lprob .eq. 2 ) then
         prbtyp    = 'LP'
         m      = 0
         linObj = .true.
         factrz = .true.
      else if (lprob .eq. 3 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .false.
         factrz = .true.
      else if (lprob .eq. 4 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .true.
         factrz = .true.
      else if (lprob .eq. 5 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .false.
         factrz = .false.
      else if (lprob .eq. 6 ) then
         prbtyp    = 'QP'
         m      = mm
         linObj = .true.
         factrz = .false.
      else if (lprob .eq. 7 ) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .false.
         factrz = .true.
      else if (lprob .eq. 8 ) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .true.
         factrz = .true.
      else if (lprob .eq. 9 ) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .false.
         factrz = .false.
      else if (lprob .eq. 10) then
         prbtyp    = 'LS'
         m      = mm
         linObj = .true.
         factrz = .false.
      end if

*     Assign the dimensions of arrays in the parameter list of lscore.
*     Economies of storage are possible if the minimum number of active
*     constraints and the minimum number of fixed variables are known in
*     advance.  The expert user should alter minact and minfxd
*     accordingly.
*     If a linear program is being solved and the matrix of general
*     constraints is fat,  i.e.,  nclin .lt. n,  a non-zero value is
*     known for minfxd.  Note that in this case, vertex must be
*     set  .true..

      minact = 0
      minfxd = 0

      vertex = .false.
      if (      (prbtyp .eq. 'LP'  .or.  prbtyp .eq. 'FP')
     $    .and.  nclin  .lt. n   ) then
         minfxd = n - nclin - 1
         vertex = .true.
      end if

      mxfree = n - minfxd
      maxact = max( 1, min( n, nclin ) )
      maxnZ  = n - ( minfxd + minact )

      if (nclin .eq. 0) then
         ldQ    = 1
         ldT    = 1
         ncolt  = 1
         vertex = .false.
      else
         ldQ    = max( 1, mxfree )
         ldT    = max( maxnZ, maxact )
         ncolt  = mxfree
      end if

*     Allocate certain arrays that are not done in lsloc.

      litotl = 0

      lAx    = 1
      lwtotl = lAx + nclin  - 1

*     Allocate remaining work arrays.

      call lsloc ( lprob, n, nclin, litotl, lwtotl )

      cold  = lcrash .eq. 0

*     Check input parameters and storage limits.

      ncnln  = 0
      lennam = 1

      call cmchk ( nerror, msglvl, lcrash, (.not.factrz),
     $             leniw, lenw, litotl, lwtotl,
     $             n, nclin, ncnln,
     $             istate, kx, named, names,
     $             bigbnd, bl, bu, x )

      if (nerror .gt. 0) then
         inform = 6
         go to 800
      end if

      lkactv = locls( 1)

      lanorm = locls( 2)
      lpx    = locls( 4)
      lres   = locls( 5)
      lres0  = locls( 6)
      lgq    = locls( 8)
      lcq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk   = locls(14)
      lfeatl = locls(15)

      if (tolfea .gt. zero)
     $   call dload ( n+nclin, (tolfea), w(lfeatl), 1 )

      ianrmj = lanorm
      do 200, j = 1, nclin
         w(ianrmj) = dnrm2 ( n, A(j,1), ldA )
         ianrmj    = ianrmj + 1
  200 continue
      if (nclin .gt. 0)
     $   call dcond ( nclin, w(lanorm), 1, Asize, amin )

      call dcond ( nctotl, w(lfeatl), 1, feamax, feamin )
      call dcopy ( nctotl, w(lfeatl), 1, w(lwtinf), 1 )
      call dscal ( nctotl, (one/feamin), w(lwtinf), 1 )

      ssq1   = zero

      if (factrz) then
*        ===============================================================
*        Factorize R using QR or Cholesky.  kx must be initialized.
*        ===============================================================
         do 210, i = 1, n
            kx(i) = i
  210    continue

         if      (prbtyp .eq. 'LP'  .or.  prbtyp .eq. 'FP') then
            nrank = 0
         else if (prbtyp .eq. 'QP') then
*           ------------------------------------------------------------
*           Compute the Cholesky factorization of R.  The Hessian is
*           m by m and resides in the upper left-hand corner of R.
*           ------------------------------------------------------------
            do 220, j = m+1, n
               call dload ( m, (zero), R(1,j), 1 )
  220       continue

            call lschol( ldR, m, nrank, tolrnk, kx, R, info )

            if (nrank .gt. 0)
     $         call dload ( nrank, (zero), w(lres0), 1 )
         else if (prbtyp .eq. 'LS') then
*           ------------------------------------------------------------
*           Compute the orthogonal factorization PRQ = ( U ),  where P
*                                                      ( 0 )
*           is an orthogonal matrix and Q is a permutation matrix.
*           Overwrite R with the upper-triangle U.  The orthogonal
*           matrix P is applied to the residual and discarded.  The
*           permutation is stored in the array KX.  Once U has been
*           computed we need only work with vectors of length N within
*           lscore.  However, it is necessary to store the sum of
*           squares of the terms  b(nrank+1),...,b(m),  where B = Pr.
*           ------------------------------------------------------------
            call dgeqrp( 'Column iterchanges', m, n, R, ldR,
     $                   w(lwrk), iw(lkactv), w(lgq), info )

            lj  = lkactv
            do 230, j = 1, n
               jmax = iw(lj)
               if (jmax .gt. j) then
                  jsave    = kx(jmax)
                  kx(jmax) = kx(j)
                  kx(j)    = jsave
               end if
               lj = lj + 1
  230       continue

            call dgeapq( 'Transpose', 'Separate', m, min( n,m-1 ), 
     $                   R, ldR, w(lwrk), 1, b, m, w(lgq), info )

            rownrm = dnrm2 ( n, R(1,1), ldR )
            if (          rownrm  .le.        tolrnk
     $          .or.  abs(R(1,1)) .le. rownrm*tolrnk) then
               nrank = 0
            else
               nrank = idrank( min(n, m), R, ldR+1, tolrnk )
            end if

            if (m .gt. nrank) ssq1 = dnrm2 ( m-nrank, b(nrank+1), 1 )

            if (nrank .gt. 0)
     $         call dcopy ( nrank, b, 1, w(lres0), 1 )
         end if
      else
*        ===============================================================
*        R is input as an upper-triangular matrix with m rows.
*        ===============================================================
         nrank = m
         if (nrank .gt. 0) then
            if      (prbtyp .eq. 'QP') then
               call dload ( nrank, (zero), w(lres0), 1 )
            else if (prbtyp .eq. 'LS') then
               call dcopy ( nrank, b, 1, w(lres0), 1 )
            end if
         end if
      end if

      if (       msglvl .gt. 0     .and.  nrank  .lt. n
     $    .and.  prbtyp .ne. 'LP'  .and.  prbtyp .ne. 'FP') then
         if (iPrint .gt. 0) write(iPrint, 9000) nrank
      end if
*     ------------------------------------------------------------------
*     Find an initial working set.
*     ------------------------------------------------------------------
      call lscrsh( cold, vertex,
     $             nclin, nctotl, nactiv, nartif,
     $             nfree, n, ldA,
     $             istate, iw(lkactv),
     $             bigbnd, tolact,
     $             A, w(lAx), bl, bu, x, w(lgq), w(lwrk) )

*     ------------------------------------------------------------------
*     Compute the TQ factorization of the constraints while keeping R in
*     upper-triangular form.  Transformations associated with Q are
*     applied to cq.  Transformations associated with P are applied to
*     res0.  If some simple bounds are in the working set,  kx is
*     re-ordered so that the free variables come first.
*     ------------------------------------------------------------------
*     First, add the bounds. To save a bit of work, cq is not loaded
*     until after kx has been re-ordered.

      ngq   = 0
      nres  = 0
      if (nrank .gt. 0) nres = 1
      unitQ = .true.

      call lsbnds( unitQ,
     $             inform, nZ, nfree, nrank, nres, ngq,
     $             n, ldQ, ldA, ldR, ldT,
     $             istate, kx, condmx,
     $             A, R, w(lT), w(lres0), w(lcq), w(lQ),
     $             w(lwrk), w(lpx), w(lrlam) )

      if (linObj) then

*        Install the transformed linear term in cq.
*        cmqmul applies the permutations in kx to cvec.

         ngq = 1
         call dcopy ( n, cvec, 1, w(lcq), 1 )
         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ,
     $                kx, w(lcq), w(lQ), w(lwrk) )
      end if

      if (nactiv .gt. 0) then
         nact1  = nactiv
         nactiv = 0

         call lsadds( unitQ, vertex,
     $                inform, 1, nact1, nactiv, nartif, nZ, nfree,
     $                nrank, nrejtd, nres, ngq,
     $                n, ldQ, ldA, ldR, ldT,
     $                istate, iw(lkactv), kx, condmx,
     $                A, R, w(lT), w(lres0), w(lcq), w(lQ),
     $                w(lwrk), w(lpx), w(lrlam) )
      end if

*     ------------------------------------------------------------------
*     Move the initial  x  onto the constraints in the working set.
*     Compute the transformed residual vector  Pr = Pb - RQ'x.
*     ------------------------------------------------------------------
      call lssetx( linObj, rowerr, unitQ,
     $             nclin, nactiv, nfree, nrank, nZ,
     $             n, nctotl, ldQ, ldA, ldR, ldT,
     $             istate, iw(lkactv), kx,
     $             jmax, errmax, ctx, xnorm,
     $             A, w(lAx), bl, bu, w(lcq), w(lres), w(lres0),
     $             w(lfeatl), R, w(lT), x, w(lQ), w(lpx), w(lwrk) )

      if (rowerr) then
         if (iPrint .gt. 0) write(iPrint, 9010)
         inform = 3
         numinf = 1
         suminf = errmax
         go to 800
      end if

      jinf = 0

      call lscore( prbtyp, named, names, linObj, unitQ,
     $             inform, iter, jinf, nclin, nctotl,
     $             nactiv, nfree, nrank, nZ, nZr,
     $             n, ldA, ldR,
     $             istate, iw(lkactv), kx,
     $             ctx, obj, ssq1,
     $             suminf, numinf, xnorm,
     $             bl, bu, A, clamda, w(lAx),
     $             w(lfeatl), R, x, w )

      obj    = obj    + ctx
      if (prbtyp .eq. 'LS'  .and.  nrank .gt. 0)
     $   call dcopy ( nrank, w(lres), 1, b, 1 )

*     ==================================================================
*     Print messages if required.
*     ==================================================================
  800 if (msglvl .gt.   0) then
         if (inform .eq.   0) then
            if (prbtyp .eq. 'FP') then
               if (iPrint .gt. 0) write(iPrint, 2001)
               if (iSumm  .gt. 0) write(iSumm , 2001)
            else
               if (iPrint .gt. 0) write(iPrint, 2002) prbtyp
               if (iSumm  .gt. 0) write(iSumm , 2002) prbtyp
            end if
         end if
         if (iPrint .gt. 0) then
            if (inform .eq.   1) write(iPrint, 2010) prbtyp
            if (inform .eq.   2) write(iPrint, 2020) prbtyp
            if (inform .eq.   3) write(iPrint, 2030)
            if (inform .eq.   4) write(iPrint, 2040)
            if (inform .eq.   5) write(iPrint, 2050)
            if (inform .eq.   6) write(iPrint, 2060) nerror
         end if

         if (iSumm  .gt. 0) then
            if (inform .eq.   1) write(iSumm , 2010) prbtyp
            if (inform .eq.   2) write(iSumm , 2020) prbtyp
            if (inform .eq.   3) write(iSumm , 2030)
            if (inform .eq.   4) write(iSumm , 2040)
            if (inform .eq.   5) write(iSumm , 2050)
            if (inform .eq.   6) write(iSumm , 2060) nerror
         end if
            
         if (inform .lt.   6) then
            if      (numinf .eq. 0) then
                if (prbtyp .ne. 'FP') then
                   if (iPrint .gt. 0) write(iPrint, 3000) prbtyp, obj
                   if (iSumm  .gt. 0) write(iSumm , 3000) prbtyp, obj
                end if
            else if (inform .eq. 3) then
               if (iPrint .gt. 0) write(iPrint, 3010) suminf
               if (iSumm  .gt. 0) write(iSumm , 3010) suminf
            else
               if (iPrint .gt. 0) write(iPrint, 3020) suminf
               if (iSumm  .gt. 0) write(iSumm , 3020) suminf
            end if
            if (numinf .gt. 0) obj = suminf
         end if
      end if

*     Recover the optional parameters set by the user.

      call icopy ( mxparm, ipsvls, 1, iprmls, 1 )
      call dcopy ( mxparm, rpsvls, 1, rprmls, 1 )

      return

 2001 format(/ ' Exit LSSOL - Feasible point found.     ')
 2002 format(/ ' Exit LSSOL - Optimal ', A2, ' solution.')
 2010 format(/ ' Exit LSSOL - Weak ',    A2, ' solution.')
 2020 format(/ ' Exit LSSOL - ', A2,         ' solution is unbounded.' )
 2030 format(/ ' Exit LSSOL - Cannot satisfy the linear constraints. ' )
 2040 format(/ ' Exit LSSOL - Too many iterations.')
 2050 format(/ ' Exit LSSOL - Too many iterations without changing X.' )
 2060 format(/ ' Exit LSSOL - ', I10, ' errors found in the input',
     $         ' parameters.  Problem abandoned.'         )
 3000 format(/ ' Final ', A2, ' objective value =', G16.7 )
 3010 format(/ ' Minimum sum of infeasibilities =', G16.7 )
 3020 format(/ ' Final sum of infeasibilities =',   G16.7 )

 9000 format(/ ' Rank of the objective function data matrix = ', I5 )
 9010 format(  ' XXX  Cannot satisfy the working-set constraints to',
     $         ' the accuracy requested.')

*     end of lssol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsadd ( unitQ,
     $                   inform, ifix, iadd, jadd,
     $                   nactiv, nZ, nfree, nrank, nres, ngq,
     $                   n, ldA, ldQ, ldR, ldT,
     $                   kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kx(n)
      double precision   A(ldA,*), R(ldR,*), T(ldT,*),
     $                   res(n,*), gqm(n,*), Q(ldQ,*)
      double precision   w(n), c(n), s(n)

*     ==================================================================
*     lsadd   updates the factorization,  A(free) * (Z Y) = (0 T),  when
*     a constraint is added to the working set.  If  nrank .gt. 0, the
*     factorization  ( R ) = PCQ  is also updated,  where  C  is the
*                    ( 0 )
*     least squares matrix,  R  is upper-triangular,  and  P  is an
*     orthogonal matrix.  The matrices  C  and  P  are not stored.
*
*     There are three separate cases to consider (although each case
*     shares code with another)...
*
*     (1) A free variable becomes fixed on one of its bounds when there
*         are already some general constraints in the working set.
*
*     (2) A free variable becomes fixed on one of its bounds when there
*         are only bound constraints in the working set.
*
*     (3) A general constraint (corresponding to row  iadd  of  A) is
*         added to the working set.
*
*     In cases (1) and (2), we assume that  kx(ifix) = jadd.
*     In all cases,  jadd  is the index of the constraint being added.
*
*     If there are no general constraints in the working set,  the
*     matrix  Q = (Z Y)  is the identity and will not be touched.
*
*     If  nres .gt. 0,  the row transformations are applied to the rows
*     of the  (n by nres)  matrix  res.
*     If  ngq .gt. 0,  the column transformations are applied to the
*     columns of the  (ngq by n)  matrix  gqm'.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October--1984.
*     Level-2 matrix routines added 25-Apr-1988.
*     This version of  lsadd  dated 14-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize, dTmax, dTmin

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      logical            bound , overfl
      external           ddiv  , dnrm2
      intrinsic          max   , min
      parameter         (zero = 0.0d+0, one = 1.0d+0)

*     if the condition estimator of the updated factors is greater than
*     condbd,  a warning message is printed.

      condbd = one / epspt9

      overfl = .false.
      bound  = jadd .le. n

      if (bound) then
*        ===============================================================
*        A simple bound has entered the working set.  iadd  is not used.
*        ===============================================================
         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $      write(iPrint, 1010) nactiv, nZ, nfree, ifix, jadd, unitQ
         nanew = nactiv

         if (unitQ) then

*           Q  is not stored, but kx defines an ordering of the columns
*           of the identity matrix that implicitly define  Q.
*           Define the sequence of pairwise interchanges p that moves
*           the newly-fixed variable to position nfree.
*           Reorder kx accordingly.

            do 100, i = 1, nfree-1
               if (i .ge. ifix) then
                  w (i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if
  100       continue
         else
*           ------------------------------------------------------------
*           Q  is stored explicitly.
*           ------------------------------------------------------------
*           Set  w = the  (ifix)-th  row of  Q.
*           Move the  (nfree)-th  row of  Q  to position  ifix.

            call dcopy ( nfree, Q(ifix,1), ldQ, w, 1 )
            if (ifix .lt. nfree) then
               call dcopy ( nfree, Q(nfree,1), ldQ, Q(ifix,1), ldQ )
               kx(ifix) = kx(nfree)
            end if
         end if
         kx(nfree) = jadd
      else
*        ===============================================================
*        A general constraint has entered the working set.
*        ifix  is not used.
*        ===============================================================
         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $      write(iPrint, 1020) nactiv, nZ, nfree, iadd, jadd, unitQ

         nanew  = nactiv + 1

*        Transform the incoming row of  A  by  Q'.  Use c as workspace.

         call dcopy ( n, A(iadd,1), ldA, w, 1 )
         call cmqmul( 8, n, nZ, nfree, ldQ, unitQ, kx, w, Q, c )

*        Check that the incoming row is not dependent upon those
*        already in the working set.

         dTnew  = dnrm2 ( nZ, w, 1 )
         if (nactiv .eq. 0) then

*           This is the only general constraint in the working set.

            cond   = ddiv  ( Asize, dTnew, overfl )
            tdTmax = dTnew
            tdTmin = dTnew
         else

*           There are already some general constraints in the working
*           set. Update the estimate of the condition number.

            tdTmax = max( dTnew, dTmax )
            tdTmin = min( dTnew, dTmin )
            cond   = ddiv  ( tdTmax, tdTmin, overfl )
         end if

         if (cond .gt. condmx  .or.  overfl) go to 900

         if (unitQ) then

*           First general constraint added.  Set  Q = I.

            call f06qhf( 'General', nfree, nfree, zero, one, Q, ldQ )
            unitQ  = .false.
         end if
      end if

      if (bound) then
         npiv  = nfree
      else
         npiv  = nZ
      end if

      nT = min( nrank, npiv )

      if (unitQ) then
*        ---------------------------------------------------------------
*        Q (i.e., Q) is not stored explicitly.
*        Apply the sequence of pairwise interchanges P that moves the
*        newly-fixed variable to position nfree.
*        ---------------------------------------------------------------
         if (ngq .gt. 0)
     $      call f06qkf( 'Left', 'Transpose', nfree-1, w, ngq, gqm, n )
            
         if (nrank .gt. 0) then

*           Apply the pairwise interchanges to the triangular part of R.
*           The subdiagonal elements generated by this process are
*           stored in  s(1), s(2), ..., s(nt-1).

            call f06qnf( 'Right', n, ifix, nT, s, R, ldR )

            if (nt .lt. npiv) then

*              R is upper trapezoidal.  Apply the interchanges in
*              columns  nT  thru  npiv.

               do 200, i = ifix, nT-1
                  w(i) = i
  200          continue

               call f06qkf( 'Right', 'Normal', nfree-1, w, nT, R, ldR )
            end if
            
*           Eliminate the subdiagonal elements of R with a left-hand
*           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
*           Apply P2 to res.

            call f06qrf( 'Left ', n, ifix, nT, c, s, R, ldR )
            if (nres .gt. 0) 
     $         call f06qxf( 'Left', 'Variable', 'Forwards', nT, nres,
     $                      ifix, nT, c, s, res, n )
         end if
      else
*        ---------------------------------------------------------------
*        Full matrix Q.  Define a sweep of plane rotations P such that
*                           Pw = beta*e(npiv).
*        The rotations are applied in the planes (1,2), (2,3), ...,
*        (npiv-1,npiv).  The rotations must be applied to Q, R, T
*        and GQM'.
*        ---------------------------------------------------------------
         call f06fqf( 'Varble', 'Forwrds', npiv-1, w(npiv), w, 1, c, s )

         if (bound  .and.  nactiv .gt. 0) then

            call dcopy ( nactiv, s(nZ), 1, w(nZ), 1 )
         
            s(       nZ  ) = s(nZ)*T(nactiv,nZ+1)
            T(nactiv,nZ+1) = c(nZ)*T(nactiv,nZ+1)
         
            call f06qzf( 'Create', nactiv, 1, nactiv, c(nZ+1), s(nZ+1),
     $                   T(1,nZ+1), ldT )
            call dcopy ( nactiv, s(nZ), 1, T(nactiv,nZ), ldT-1 )
         
            call dcopy ( nactiv, w(nZ), 1, s(nZ), 1 )
         end if

         if (ngq .gt. 0)
     $      call f06qxf( 'Left ', 'Variable', 'Forwards', npiv , ngq,
     $                   1, npiv, c, s, gqm, n )
         call f06qxf( 'Right', 'Variable', 'Forwards', nfree, nfree,
     $                1, npiv, c, s, Q, ldQ )

         if (nrank .gt. 0) then

*           Apply the rotations to the triangular part of R.
*           The subdiagonal elements generated by this process are
*           stored in  s(1),  s(2), ..., s(nt-1).

            nT = min( nrank, npiv )
            call f06qvf( 'Right', n, 1, nT, c, s, R, ldR )

            if (nt .lt. npiv) then

*              R is upper trapezoidal.  Pretend R is (nt x n) and
*              apply the rotations in columns  nT  thru  npiv.

               call f06qxf( 'Right', 'Variable', 'Forwards', nT, n,
     $                      nT, npiv, c, s, R, ldR )
            end if

*           Eliminate the subdiagonal elements of R with a left-hand
*           sweep of rotations P2 in planes (1,2), (2,3), ...,(nt-1,nt).
*           Apply P2 to res.

            call f06qrf( 'Left ', n, 1, nT, c, s, R, ldR )
            if (nres .gt. 0)
     $         call f06qxf( 'Left', 'Variable', 'Forwards', nT, nres,
     $                      1, nT, c, s, res, n )
         end if

         if (bound) then

*           The last row and column of Q has been transformed to plus
*           or minus the unit vector e(nfree).  We can reconstitute the
*           columns of GQM and R corresponding to the new fixed variable.

            if (w(nfree) .lt. zero) then
               nf = min( nrank, nfree )
               if (nf  .gt. 0) call dscal ( nf , -one,   R(1,nfree), 1 )
               if (ngq .gt. 0) call dscal ( ngq, -one, gqm(nfree,1), n )
            end if

*           ------------------------------------------------------------
*           The diagonals of T have been altered.  Recompute the
*           largest and smallest values.
*           ------------------------------------------------------------
            if (nactiv .gt. 0) then
               call dcond( nactiv, T(nactiv,nZ), ldT-1, tdTmax, tdTmin )
               cond   = ddiv  ( tdTmax, tdTmin, overfl )
            end if
         else
*           ------------------------------------------------------------
*           General constraint.  Install the new row of T.
*           ------------------------------------------------------------
            call dcopy ( nanew, w(nZ), 1, T(nanew,nZ), ldT )
         end if
      end if

*     ==================================================================
*     Prepare to exit.  Check the magnitude of the condition estimator.
*     ==================================================================
  900 if (nanew .gt. 0) then
         if (cond .lt. condmx  .and.  .not. overfl) then

*           The factorization has been successfully updated.

            inform = 0
            dTmax  = tdTmax
            dTmin  = tdTmin
            if (cond .ge. condbd) then
               if (iPrint .gt. 0) write(iPrint, 2000) jadd
            end if
         else

*           The proposed working set appears to be linearly dependent.

            inform = 1
            if (lsdbg  .and.  ilsdbg(1) .gt. 0) then
               write( iPrint, 3000 )
               if (bound) then
                  write(iPrint, 3010) Asize, dTmax, dTmin
               else
                  if (nactiv .gt. 0) then
                     write(iPrint, 3020) Asize, dTmax, dTmin, dTnew
                  else
                     write(iPrint, 3030) Asize, dTnew
                  end if
               end if
            end if
         end if
      end if

      return

 1010 format(/ ' //lsadd //  Simple bound added.'
     $       / ' //lsadd //  nactiv    nZ nfree  ifix  jadd unitQ'
     $       / ' //lsadd //  ', 5i6, l6 )
 1020 format(/ ' //lsadd //  General constraint added.           '
     $       / ' //lsadd //  nactiv    nZ nfree  iadd  jadd unitQ'
     $       / ' //lsadd //  ', 5i6, l6 )
 2000 format(/ ' XXX  Serious ill-conditioning in the working set',
     $         ' after adding constraint ',  i5
     $       / ' XXX  Overflow may occur in subsequent iterations.'//)
 3000 format(/ ' //lsadd //  Dependent constraint rejected.' )
 3010 format(/ ' //lsadd //     Asize     dTmax     dTmin        '
     $       / ' //lsadd //', 1p, 3e10.2 )
 3020 format(/ ' //lsadd //     Asize     dTmax     dTmin     dTnew'
     $       / ' //lsadd //', 1p, 4e10.2 )
 3030 format(/ ' //lsadd //     Asize     dTnew'
     $       / ' //lsadd //', 1p, 2e10.2 )

*     end of lsadd
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsadds( unitQ, vertex,
     $                   inform, k1, k2, nactiv, nartif, nZ, nfree,
     $                   nrank, nrejtd, nres, ngq,
     $                   n, ldQ, ldA, ldR, ldT,
     $                   istate, kactiv, kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ, vertex
      integer            istate(*), kactiv(n), kx(n)
      double precision   condmx
      double precision   A(ldA,*), R(ldR,*),
     $                   T(ldT,*), res(n,*), gqm(n,*), Q(ldQ,*)
      double precision   w(n), c(n), s(n)

*     ==================================================================
*     lsadds  includes general constraints k1 thru k2 as new rows of
*     the TQ factorization stored in T, Q.  If nrank is nonZero, the
*     changes in Q are reflected in nrank by n triangular factor R such
*     that
*                         C  =  P ( R ) Q,
*                                 ( 0 )
*     where  P  is orthogonal.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  October-31-1984.
*     This version of lsadds dated  16-May-1988.
*     ==================================================================

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol5cm/ Asize, dTmax, dTmin

      external           dnrm2
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      rtmax  = wmach(8)

*     Estimate the condition number of the constraints that are not
*     to be refactorized.

      if (nactiv .eq. 0) then
         dTmax = zero
         dTmin = one
      else
         call dcond ( nactiv, T(nactiv,nZ+1), ldT-1, dTmax, dTmin )
      end if

      do 200, k = k1, k2
         iadd = kactiv(k)
         jadd = n + iadd
         if (nactiv .lt. nfree) then

            call lsadd ( unitQ,
     $                   inform, ifix, iadd, jadd,
     $                   nactiv, nZ, nfree, nrank, nres, ngq,
     $                   n, ldA, ldQ, ldR, ldT,
     $                   kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

            if (inform .eq. 0) then
               nactiv = nactiv + 1
               nZ     = nZ     - 1
            else
               istate(jadd) =   0
               kactiv(k)    = - kactiv(k)
            end if
         end if
  200 continue

      if (nactiv .lt. k2) then

*        Some of the constraints were classed as dependent and not
*        included in the factorization.  Re-order the part of kactiv
*        that holds the indices of the general constraints in the
*        working set.  Move accepted indices to the front and shift
*        rejected indices (with negative values) to the end.

         l      = k1 - 1
         do 300, k = k1, k2
            i         = kactiv(k)
            if (i .ge. 0) then
               l      = l + 1
               if (l .ne. k) then
                  iswap     = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
  300    continue

*        If a vertex is required, add some temporary bounds.
*        We must accept the resulting condition number of the working
*        set.

         if (vertex) then
            cndmax = rtmax
            nZadd  = nZ
            do 320, iartif = 1, nZadd
               if (unitQ) then
                  ifix = nfree
                  jadd = kx(ifix)
               else
                  rowmax = zero
                  do 310, i = 1, nfree
                     rnorm = dnrm2 ( nZ, Q(i,1), ldQ )
                     if (rowmax .lt. rnorm) then
                        rowmax = rnorm
                        ifix   = i
                     end if
  310             continue
                  jadd = kx(ifix)

                  call lsadd ( unitQ,
     $                         inform, ifix, iadd, jadd,
     $                         nactiv, nZ, nfree, nrank, nres, ngq,
     $                         n, ldA, ldQ, ldR, ldT,
     $                         kx, cndmax,
     $                         A, R, T, res, gqm, Q,
     $                         w, c, s  )
               end if
               nfree  = nfree  - 1
               nZ     = nZ     - 1
               nartif = nartif + 1
               istate(jadd) = 4
  320       continue
         end if
      end if

      nrejtd = k2 - nactiv

*     end of lsadds
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsbnds( unitQ,
     $                   inform, nZ, nfree, nrank, nres, ngq,
     $                   n, ldQ, ldA, ldR, ldT,
     $                   istate, kx, condmx,
     $                   A, R, T, res, gqm, Q,
     $                   w, c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            istate(*), kx(n)
      double precision   condmx
      double precision   A(ldA,*), R(ldR,*),
     $                   T(ldT,*), res(n,*), gqm(n,*), Q(ldQ,*)
      double precision   w(n), c(n), s(n)

*     ==================================================================
*     lsbnds updates the factor R as kx is reordered to reflect the
*     status of the bound constraints given by istate.  kx is reordered
*     so that the fixed variables come last.  One of two alternative
*     methods are used to reorder kx. One method needs fewer accesses 
*     to kx, the other gives a matrix  Rz  with more rows and columns.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  30-December-1985.
*     This version of lsbnds dated 13-May-88.
*     ==================================================================

      nfixed = n - nfree

      if (nrank .lt. n  .and.  nrank .gt. 0) then
*        ---------------------------------------------------------------
*        R is specified but singular.  Try and keep the dimension of Rz
*        as large as possible.
*        ---------------------------------------------------------------
         nactv = 0
         nfree = n
         nZ    = n

         j     = n
*+       while (j .gt. 0  .and.  n-nfree .lt. nfixed) do
  100    if    (j .gt. 0  .and.  n-nfree .lt. nfixed) then
            if (istate(j) .gt. 0) then
               jadd = j
               do 110, ifix = nfree, 1, -1
                  if (kx(ifix) .eq. jadd) go to 120
  110          continue

*              Add bound jadd.

  120          call lsadd ( unitQ,
     $                      inform, ifix, iadd, jadd,
     $                      nactv, nZ, nfree, nrank, nres, ngq,
     $                      n, ldA, ldQ, ldR, ldT,
     $                      kx, condmx,
     $                      A, R, T, res, gqm, Q,
     $                      w, c, s )

               nfree = nfree - 1
               nZ    = nZ    - 1
            end if
            j = j - 1
            go to 100
*+       end while
         end if
      else     
*        ---------------------------------------------------------------
*        R is of full rank,  or is not specified.
*        ---------------------------------------------------------------
         if (nfixed .gt. 0) then

*           Order kx so that the free variables come first.

            lstart = nfree + 1
            do 250, k = 1, nfree
               j = kx(k)
               if (istate(j) .gt. 0) then
                  do 220, l = lstart, n
                     j2 = kx(l)
                     if (istate(j2) .eq. 0) go to 230
  220             continue

  230             kx(k)  = j2
                  kx(l)  = j
                  lstart = l + 1

                  if (nrank .gt. 0)
     $               call cmrswp( n, nres, nrank, ldR, k, l,
     $                            R, res, c, s )
               end if
  250       continue

         end if
         nZ = nfree
      end if

*     end of lsbnds
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lschol( ldh, n, nrank, tolrnk, kx, h, inform )

      implicit           double precision (a-h,o-z)
      integer            kx(*)
      double precision   h(ldh,*)

*     ==================================================================
*     LSCHOL  forms the Cholesky factorization of the positive
*     semi-definite matrix H such that
*                   PHP'  =  R'R
*     where  P  is a permutation matrix and  R  is upper triangular.
*     The permutation P is chosen to maximize the diagonal of R at each
*     stage.  Only the diagonal and super-diagonal elements of H are
*     used.
*
*     Output:
*
*         inform = 0   the factorization was computed successfully,
*                      with the Cholesky factor written in the upper
*                      triangular part of H and P stored in kx.
*                  1   the matrix H was indefinite.
*
*     Original version of lschol dated  2-February-1981.
*     Level 2 Blas added 29-June-1986.
*     This version of lschol dated  26-Jun-1989. 
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      intrinsic          abs   , max   , sqrt
      external           idamax
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      inform = 0
      nrank  = 0

*     Main loop for computing rows of  R.

      do 200, j = 1, n

*        Find maximum available diagonal.    

         kmax = j - 1 + idamax( n-j+1, h(j,j), ldh+1 )
         dmax = h(kmax,kmax)

         if (dmax .le. tolrnk*abs(h(1,1))) go to 300

*        Perform a symmetric interchange if necessary.

         if (kmax .ne. j) then
            k        = kx(kmax)
            kx(kmax) = kx(j)
            kx(j)    = k

            call dswap ( kmax-j, h(j+1,kmax), 1, h(j,j+1 ), ldh )
            call dswap ( j     , h(1  ,j   ), 1, h(1,kmax), 1   )
            call dswap ( n-kmax+1, h(kmax,kmax), ldh,
     $                             h(j,kmax)   , ldh )

         end if

*        Set the diagonal of  R.

         d      = sqrt( dmax )
         h(j,j) = d
         nrank  = nrank + 1

         if (j .lt. n) then

*           Set the super-diagonal elements of this row of R and update
*           the elements of the block that is yet to be factorized.
                                                          
            call dscal ( n-j,   (one/d), h(j  ,j+1), ldh )
            call dsyr  ( 'u', n-j, -one, h(j  ,j+1), ldh,
     $                                   h(j+1,j+1), ldh )
         end if

  200 continue
*     ------------------------------------------------------------------
*     Check for the semi-definite case.
*     ------------------------------------------------------------------
  300 if (nrank .lt. n) then

*        Find the largest element in the unfactorized block.

         supmax = zero
         do 310, i = j, n-1
            k      = i + idamax( n-i, h(i,i+1), ldh )
            supmax = max( supmax, abs(h(i,k)) )
  310    continue

         if (supmax .gt. tolrnk*abs(h(1,1))) then
            if (iPrint .gt. 0) write(iPrint, 1000) dmax, supmax
            inform = 1
         end if
      end if

      return

 1000 format(' XXX  Hessian appears to be indefinite.'
     $      /' XXX  Maximum diagonal and off-diagonal ignored',
     $             ' in the Cholesky factorization:', 1p, 2e22.14 )

*     end of lschol
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lscore( prbtyp, named, names, linObj, unitQ,
     $                   inform, iter, jinf, nclin, nctotl,
     $                   nactiv, nfree, nrank, nZ, nZr,
     $                   n, ldA, ldR,
     $                   istate, kactiv, kx,
     $                   ctx, ssq, ssq1, suminf, numinf, xnorm,
     $                   bl, bu, A, clamda, Ax,
     $                   featol, R, x, w )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      character*8        names(*)
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   bl(nctotl), bu(nctotl), A(ldA,*),
     $                   clamda(nctotl), Ax(*),
     $                   featol(nctotl), R(ldR,*), x(n)
      double precision   w(*)
      logical            named, linObj, unitQ

*     ==================================================================
*     lscore  is a subroutine for linearly constrained linear-least
*     squares.  On entry, it is assumed that an initial working set of
*     linear constraints and bounds is available.
*     The arrays istate, kactiv and kx will have been set accordingly
*     and the arrays t and Q will contain the TQ factorization of
*     the matrix whose rows are the gradients of the active linear
*     constraints with the columns corresponding to the active bounds
*     removed.  the TQ factorization of the resulting (nactiv by nfree)
*     matrix is  A(free)*Q = (0 T),  where Q is (nfree by nfree) and t
*     is reverse-triangular.
*
*     Values of istate(j) for the linear constraints.......
*
*     istate(j)
*     ---------
*          0    constraint j is not in the working set.
*          1    constraint j is in the working set at its lower bound.
*          2    constraint j is in the working set at its upper bound.
*          3    constraint j is in the working set as an equality.
*
*     Constraint j may be violated by as much as featol(j).
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  lscore  dated  14-Sep-92.
*
*     Copyright  1984/1993  Stanford University.
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
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT   , ncolt , ldQ
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9
      common    /sol5cm/ Asize, dTmax  , dTmin

      integer            locls
      parameter         (lenls = 20)
      common    /sol1ls/ locls(lenls)

      logical            cmdbg, lsdbg
      parameter         (ldbg = 5)
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
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      external           ddiv  , dnrm2
      intrinsic          abs   , max   , sqrt
      logical            convrg, cyclin, error , firstv, hitcon,
     $                   hitlow, needfg, overfl, prnt  , rowerr
      logical            singlr, stall , statpt, unbndd, uncon , unitgZ,
     $                   weak
      parameter        ( zero   =0.0d+0, half   =0.5d+0, one   =1.0d+0 )
      parameter        ( mrefn  =1     , mstall =50                    )

*     Specify the machine-dependent parameters.

      flmax  = wmach(7)

      lanorm = locls( 2)
      lAp    = locls( 3)
      lpx    = locls( 4)
      lres   = locls( 5)
      lres0  = locls( 6)
      lhZ    = locls( 7)
      lgq    = locls( 8)
      lcq    = locls( 9)
      lrlam  = locls(10)
      lT     = locls(11)
      lQ     = locls(12)
      lwtinf = locls(13)
      lwrk   = locls(14)

*     Set up the adresses of the contiguous arrays  ( res0, res )
*     and  ( gq, cq ).

      nres   = 0
      if (nrank .gt. 0) nres = 2
      ngq    = 1
      if (linObj) ngq = 2

*     Initialize.

      irefn  =   0
      iter   =   0
      itmax  = itmax1

      jadd   =   0
      jdel   =   0
      nphase =   1
      nstall =   0
      numinf = - 1
      nZr    =   0

      alfa   = zero
      condmx = flmax
      dRzmax = one
      dRzmin = one
      ssq    = zero

      cyclin = .false.
      error  = .false.
      firstv = .false.
      prnt   = .true.
      needfg = .true.
      stall  = .true.
      uncon  = .false.
      unbndd = .false.

*     If debug output is required,  print nothing until iteration IDBG.

      msgsvd = msglvl
      if (idbg .gt. 0  .and.  idbg .le. itmax) then
         msglvl = 0
      end if

*======================== start of the main loop =======================
*
*      cyclin = false
*      unbndd = false
*      error  = false
*      k      = 0
*
*      repeat
*            repeat
*                  compute Z'g,  print details of this iteration
*                  stat pt = (Z'g .eq. 0)
*                  if (not stat pt) then
*                     error =  k .ge. itmax
*                     if (not error) then
*                        compute p, alfa
*                        error = unbndd  or  cyclin
*                        if (not error) then
*                           k = k + 1
*                           x = x + alfa p
*                           if (feasible) update Z'g
*                           if necessary, add a constraint
*                        end if
*                     end if
*                  end if
*            until  stat pt  or  error
*
*            compute lam1, lam2, smllst
*            optmul =  smllst .gt. 0
*            if ( not (optmul .or. error) ) then
*                  delete an artificial or regular constraint
*            end if
*      until optmul  or  error
*
*=======================================================================

*     repeat
*        repeat
  100       if (needfg) then
               if (nrank .gt. 0) then
                  resnrm = dnrm2 ( nrank, w(lres), 1 )
                  ssq    = half*(ssq1**2 + resnrm**2 )
               end if

               if (numinf .ne. 0) then

*                 Compute the transformed gradient of either the sum of
*                 of infeasibilities or the objective.  Initialize
*                 singlr and unitgZ.

                  call lsgset( prbtyp, linObj, singlr, unitgZ, unitQ,
     $                         n, nclin, nfree,
     $                         ldA, ldQ, ldR, nrank, nZ, nZr,
     $                         istate, kx,
     $                         bigbnd, tolrnk, numinf, suminf,
     $                         bl, bu, A, w(lres), featol,
     $                         w(lgq), w(lcq), R, x, w(lwtinf),
     $                         w(lQ), w(lwrk) )

                  if (prbtyp .ne. 'FP'  .and.  numinf .eq. 0
     $                                  .and.  nphase .eq. 1) then
                     itmax  = iter + itmax2
                     nphase = 2
                  end if
               end if
            end if

            gznorm = zero
            if (nZ  .gt. 0 ) gznorm = dnrm2 ( nZ, w(lgq), 1 )

            if (nZr .eq. nZ) then
               gZrnrm = gznorm
            else
               gZrnrm = zero
               if (nZr .gt. 0) gZrnrm = dnrm2 ( nZr, w(lgq), 1 )
            end if

            gfnorm = gznorm
            if (nfree .gt. 0  .and.  nactiv .gt. 0)
     $         gfnorm = dnrm2 ( nfree, w(lgq), 1 )

*           ------------------------------------------------------------
*           Print the details of this iteration.
*           ------------------------------------------------------------
*           Define small quantities that reflect the size of x, R and
*           the constraints in the working set.  If feasible,  estimate
*           the rank and condition number of Rz.
*           Note that nZr .le. nrank + 1.

            if (nZr .eq. 0) then
               singlr = .false.
            else
               if (numinf .gt. 0  .or.  nZr .gt. nrank) then
                  absrzz = zero
                  singlr = .true.
               else
                  call dcond ( nZr, R, ldR+1, dRzmax, dRzmin )
                  absrzz = abs( R(nZr,nZr) )
                  rownrm = dnrm2 ( n, R(1,1), ldR )
                  singlr =       absrzz      .le. dRzmax*tolrnk 
     $                     .or.  rownrm      .le.        tolrnk 
     $                     .or.  abs(R(1,1)) .le. rownrm*tolrnk
               end if

               if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $            write(iPrint, 9100) singlr, absrzz, dRzmax, dRzmin
            end if

            condRz = ddiv  ( dRzmax, dRzmin, overfl )
            condT  = one
            if (nactiv .gt. 0)
     $         condT  = ddiv  ( dTmax , dTmin , overfl )

            if (prnt) then
               call lsprt ( prbtyp, isdel, iter, jadd, jdel,
     $                      msglvl, nactiv, nfree, n, nclin,
     $                      nrank, ldR, ldT, nZ, nZr, istate,
     $                      alfa, condRz, condT, gfnorm, gZrnrm,
     $                      numinf, suminf, ctx, ssq,
     $                      Ax, R, w(lT), x, w(lwrk) )
               jdel  = 0
               jadd  = 0
               alfa  = zero
            end if

            if (numinf .gt. 0) then
               dinky  = zero
            else
               objsiz = one  + abs( ssq + ctx )
               wssize = zero
               if (nactiv .gt. 0) wssize = dTmax
               dinky  = epspt8 * max( wssize, objsiz, gfnorm )
               if (uncon) then
                  unitgZ = gZrnrm .le. dinky
               end if
            end if

            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $         write(iPrint, 9000) unitgZ, irefn, gZrnrm, dinky

*           If the projected gradient  Z'g  is small and Rz is of full
*           rank, X is a minimum on the working set.  An additional
*           refinement step is allowed to take care of an inaccurate
*           value of dinky.

            statpt = .not. singlr  .and.  gZrnrm .le. dinky
     $                             .or.   irefn  .gt. mrefn

            if (.not. statpt) then
*              ---------------------------------------------------------
*              Compute a search direction.
*              ---------------------------------------------------------
               prnt  = .true.

               error = iter .ge. itmax
               if (.not. error) then

                  irefn = irefn + 1
                  iter  = iter  + 1

                  if (iter .eq. idbg  .and.  iPrint .gt. 0) then
                     lsdbg  = .true.
                     cmdbg  =  lsdbg
                     msglvl =  msgsvd
                  end if

                  call lsgetp( linObj, singlr, unitgZ, unitQ,
     $                         n, nclin, nfree,
     $                         ldA, ldQ, ldR, nrank, numinf, nZr,
     $                         kx, ctp, pnorm,
     $                         A, w(lAp), w(lres), w(lhZ), w(lpx),
     $                         w(lgq), w(lcq), R, w(lQ), w(lwrk) )

*                 ------------------------------------------------------
*                 Find the constraint we bump into along p.
*                 Update x and Ax if the step alfa is nonZero.
*                 ------------------------------------------------------
*                 alfhit is initialized to bigalf.  If it remains
*                 that way after the call to cmalf, it will be
*                 regarded as infinite.

                  bigalf = ddiv  ( bigdx, pnorm, overfl )

                  call cmalf ( firstv, hitlow,
     $                         istate, inform, jadd, n, nctotl, numinf,
     $                         alfhit, palfa, atphit, 
     $                         bigalf, bigbnd, pnorm,
     $                         w(lAnorm), w(lAp), Ax, bl, bu, 
     $                         featol, w(lpx), x )

*                 If  Rz  is nonsingular,  alfa = 1.0  will be the
*                 step to the least-squares minimizer on the
*                 current subspace. If the unit step does not violate
*                 the nearest constraint by more than featol,  the
*                 constraint is not added to the working set.

                  hitcon = singlr  .or.  palfa  .le. one
                  uncon  = .not. hitcon

                  if (hitcon) then
                     alfa = alfhit
                  else
                     jadd   = 0
                     alfa   = one
                  end if

*                 Check for an unbounded solution or negligible step.

                  unbndd =  alfa .ge. bigalf
                  stall  = abs( alfa*pnorm ) .le. epspt9*xnorm
                  if (stall) then
                     nstall = nstall + 1
                     cyclin = nstall .gt. mstall
                  else
                     nstall = 0
                  end if

                  error = unbndd  .or.  cyclin
                  if (.not.  error) then
*                    ---------------------------------------------------
*                    Set x = x + alfa*p.  Update Ax, gq, res and ctx.
*                    ---------------------------------------------------
                     if (alfa .ne. zero)
     $                  call lsmove( hitcon, hitlow, linObj, unitgZ,
     $                               nclin, nrank, nZr,
     $                               n, ldR, jadd, numinf,
     $                               alfa, ctp, ctx, xnorm,
     $                               w(lAp), Ax, bl, bu, w(lgq),
     $                               w(lhZ), w(lpx), w(lres),
     $                               R, x, w(lwrk) )

                     if (hitcon) then
*                       ------------------------------------------------
*                       Add a constraint to the working set.
*                       Update the TQ factors of the working set.
*                       Use p as temporary work space.
*                       ------------------------------------------------
*                       Update  istate.

                        if (bl(jadd) .eq. bu(jadd)) then
                           istate(jadd) = 3
                        else if (hitlow) then
                           istate(jadd) = 1
                        else
                           istate(jadd) = 2
                        end if
                        iadd = jadd - n
                        if (jadd .le. n) then

                           do 510, ifix = 1, nfree
                              if (kx(ifix) .eq. jadd) go to 520
  510                      continue
  520                   end if

                        call lsadd ( unitQ,
     $                               inform, ifix, iadd, jadd,
     $                               nactiv, nZ, nfree, nrank, nres,ngq,
     $                               n, ldA, ldQ, ldR, ldT,
     $                               kx, condmx,
     $                               A, R, w(lT), w(lres),w(lgq),w(lQ),
     $                               w(lwrk), w(lrlam), w(lpx) )

                        nZr    = nZr - 1
                        nZ     = nZ  - 1

                        if (jadd .le. n) then

*                          A simple bound has been added.

                           nfree  = nfree  - 1
                        else

*                          A general constraint has been added.

                           nactiv = nactiv + 1
                           kactiv(nactiv) = iadd
                        end if
                        irefn  = 0
                     end if
*                    ---------------------------------------------------
*                    Check the feasibility of constraints with non-
*                    negative istate values.  If some violations have
*                    occurred.  Refine the current x and set inform so
*                    that feasibility is checked in lsgset.
*                    ---------------------------------------------------
                     call lsfeas( n, nclin, istate,
     $                            bigbnd, cnorm, err1, jmax1, nviol,
     $                            Ax, bl, bu, featol, x, w(lwrk) )

                     if (err1 .gt. featol(jmax1)) then
                        call lssetx( linObj, rowerr, unitQ,
     $                               nclin, nactiv, nfree, nrank, nZ,
     $                               n, nctotl, ldQ, ldA, ldR, ldT,
     $                               istate, kactiv, kx,
     $                               jmax1, err2, ctx, xnorm,
     $                               A, Ax, bl, bu, w(lcq),
     $                               w(lres), w(lres0), featol, R,
     $                               w(lT), x, w(lQ), w(lpx), w(lwrk) )

                        if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $                     write(iPrint, 2100) err1, err2
                        if (rowerr) then
                           if (iPrint .gt. 0) write(iPrint, 2200)
                        end if
                        uncon  =   .false.
                        irefn  =   0
                        numinf = - 1
                     end if
                     needfg = alfa .ne. zero
                  end if
               end if
            end if
*        until      statpt  .or.  error
         if (.not. (statpt  .or.  error) ) go to 100

*        ===============================================================
*        Try and find the index jdel of a constraint to drop from
*        the working set.
*        ===============================================================
         jdel   = 0

         if (numinf .eq. 0  .and.  prbtyp .eq. 'FP') then
            if (n .gt. nZ)
     $         call dload ( n-nZ, (zero), w(lrlam), 1 )
            jtiny  = 0
            jsmlst = 0
            jbigst = 0
         else

            call lsmuls( prbtyp,
     $                   msglvl, n, nactiv, nfree,
     $                   ldA, ldT, numinf, nZ, nZr,
     $                   istate, kactiv, kx, dinky,
     $                   jsmlst, ksmlst, jinf, jtiny,
     $                   jbigst, kbigst, trulam,
     $                   A, w(lanorm), w(lgq), w(lrlam),
     $                   w(lT), w(lwtinf) )
         end if

         if (.not. error) then
            if (     jsmlst .gt. 0) then

*              LSMULS found a regular constraint with multiplier less
*              than (-dinky).

               jdel   = jsmlst
               kdel   = ksmlst
               isdel  = istate(jdel)
               istate(jdel) = 0

            else if (jsmlst .lt. 0) then

               jdel   = jsmlst

            else if (numinf .gt. 0  .and.  jbigst .gt. 0) then

*              No feasible point exists for the constraints but the
*              sum of the constraint violations may be reduced by
*              moving off constraints with multipliers greater than 1.

               jdel   = jbigst
               kdel   = kbigst
               isdel  = istate(jdel)
               if (trulam .le. zero) is = - 1
               if (trulam .gt. zero) is = - 2
               istate(jdel) = is
               firstv = .true.
               numinf = numinf + 1
            end if

            if      (jdel .ne. 0  .and.  singlr) then

*              Cannot delete a constraint when Rz is singular.
*              Probably a weak minimum.

               jdel = 0
            else if (jdel .ne. 0               ) then

*              Constraint jdel has been deleted.
*              Update the matrix factorizations.

               call lsdel ( unitQ,
     $                      n, nactiv, nfree, nres, ngq, nZ, nZr,
     $                      ldA, ldQ, ldR, ldT, nrank,
     $                      jdel, kdel, kactiv, kx,
     $                      A, w(lres), R, w(lT), w(lgq), w(lQ),
     $                      w(lwrk), w(lpx) )
            end if
         end if

         irefn  =  0
         convrg =  jdel .eq. 0
         prnt   = .false.
         uncon  = .false.
         needfg = .false.

*     until       convrg  .or.  error
      if (.not.  (convrg  .or.  error)) go to 100

*  .........................End of main loop............................
      weak = jtiny .gt. 0  .or.  singlr

      if (error) then
         if (unbndd) then
            inform = 2
            if (numinf .gt. 0) inform = 3
         else if (iter .ge. itmax) then
            inform = 4
         else if (cyclin) then
            inform = 5
         end if
      else if (convrg) then
         inform = 0
         if (numinf .gt. 0) then
            inform = 3
         else if (prbtyp .ne. 'FP'  .and.  weak) then
            inform = 1
         end if
      end if

*     ------------------------------------------------------------------
*     Set   clamda.  Print the full solution.
*     ------------------------------------------------------------------
      msglvl = msgsvd
      if (msglvl .gt. 0  .and.  iPrint .gt. 0)
     $     write(iPrint, 2000) prbtyp, iter, inform

      call cmprt ( msglvl, nfree, ldA,
     $             n, nclin, nctotl, bigbnd,
     $             named, names,
     $             nactiv, istate, kactiv, kx,
     $             A, bl, bu, x, clamda, w(lrlam), x )

      return

 2000 format(/ ' Exit from ', a2, ' problem after ', i4, ' iterations.',
     $         '  inform =', i3 )
 2100 format(  ' XXX  Iterative refinement.  Maximum errors before and',
     $         ' after refinement are ',  1p, 2e14.2 )
 2200 format(  ' XXX  Warning.  Cannot satisfy the constraints to the',
     $         ' accuracy requested.')
 9000 format(/ ' //lscore//  unitgZ irefn     gZrnrm      dinky'
     $       / ' //lscore//  ', l6, i6, 1p, 2e11.2 )
 9100 format(/ ' //lscore//  singlr    abs(Rzz)      dRzmax      dRzmin'
     $       / ' //lscore//  ', l6,     1p, 3e12.4 )

*     end of lscore
      end                         

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lscrsh( cold, vertex,
     $                   nclin, nctotl, nactiv, nartif,
     $                   nfree, n, ldA,
     $                   istate, kactiv,
     $                   bigbnd, tolact,
     $                   A, Ax, bl, bu, x, wx, work )

      implicit           double precision(a-h,o-z)
      logical            cold, vertex
      integer            istate(nctotl), kactiv(n)
      double precision   A(ldA,*), Ax(*), bl(nctotl), bu(nctotl),
     $                   x(n), wx(n), work(n)

*     ==================================================================
*     lscrsh  computes the quantities  istate (optionally), kactiv,
*     nactiv, nZ and nfree  associated with the working set at x.
*     The computation depends upon the value of the input parameter
*     cold,  as follows...
*
*     cold = true.  An initial working set will be selected. First,
*                   nearly-satisfied or violated bounds are added.
*                   Next,  general linear constraints are added that
*                   have small residuals.
*
*     cold = false. The quantities kactiv, nactiv, nZ and nfree are
*                   computed from istate,  specified by the user.
*
*     Values of istate(j)....
*
*        - 2         - 1         0           1          2         3
*     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lscrsh dated 14-May-1992.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           ddot
      intrinsic          abs, min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      flmax  =   wmach(7)
      biglow = - bigbnd
      bigupp =   bigbnd

*     ------------------------------------------------------------------
*     Move the variables inside their bounds.
*     ------------------------------------------------------------------
      do 10, j = 1, n
         b1    = bl(j)
         b2    = bu(j)

         if (b1 .gt. biglow) then
            if (x(j) .lt. b1) x(j) = b1
         end if
      
         if (b2 .lt. bigupp) then
            if (x(j) .gt. b2) x(j) = b2
         end if
   10 continue

      call dcopy ( n, x, 1, wx, 1 )

      if (lsdbg) then
         if (ilsdbg(1) .gt. 0)
     $      write(iPrint, 1000) cold, nclin, nctotl
         if (ilsdbg(2) .gt. 0)
     $      write(iPrint, 1100) (wx(j), j = 1, n)
      end if

      nfixed = 0
      nactiv = 0
      nartif = 0

*     If a cold start is being made, initialize  istate.
*     If  bl(j) = bu(j),  set  istate(j)=3  for all variables and linear
*     constraints.

      if (cold) then
         do 100, j = 1, nctotl
            istate(j) = 0
            if (bl(j) .eq. bu(j)) istate(j) = 3
  100    continue
      else
         do 110, j = 1, nctotl
            if (istate(j) .gt. 3  .or.  istate(j) .lt. 0) istate(j) = 0
  110    continue
      end if

*     Initialize nfixed, nfree and kactiv.
*     Ensure that the number of bounds and general constraints in the
*     working set does not exceed n.

      do 200, j = 1, nctotl
         if (nfixed + nactiv .eq. n) istate(j) = 0
         if (istate(j) .gt. 0) then
            if (j .le. n) then
               nfixed = nfixed + 1
               if (istate(j) .eq. 1) wx(j) = bl(j)
               if (istate(j) .ge. 2) wx(j) = bu(j)
            else
               nactiv = nactiv + 1
               kactiv(nactiv) = j - n
            end if
         end if
  200 continue

*     ------------------------------------------------------------------
*     If a cold start is required,  attempt to add as many
*     constraints as possible to the working set.
*     ------------------------------------------------------------------
      if (cold) then

*        See if any bounds are violated or nearly satisfied.
*        If so,  add these bounds to the working set and set the
*        variables exactly on their bounds.

         j = n
*+       while (j .ge. 1  .and.  nfixed + nactiv .lt. n) do
  300    if    (j .ge. 1  .and.  nfixed + nactiv .lt. n) then
            if (istate(j) .eq. 0) then
               b1     = bl(j)
               b2     = bu(j)
               is     = 0
               if (b1 .gt. biglow) then
                  if (wx(j) - b1 .le. (one + abs( b1 ))*tolact) is = 1
               end if
               if (b2 .lt. bigupp) then
                  if (b2 - wx(j) .le. (one + abs( b2 ))*tolact) is = 2
               end if
               if (is .gt. 0) then
                  istate(j) = is
                  if (is .eq. 1) wx(j) = b1
                  if (is .eq. 2) wx(j) = b2
                  nfixed = nfixed + 1
               end if
            end if
            j = j - 1
            go to 300
*+       end while
         end if

*        ---------------------------------------------------------------
*        The following loop finds the linear constraint (if any) with
*        smallest residual less than or equal to tolact  and adds it
*        to the working set.  This is repeated until the working set
*        is complete or all the remaining residuals are too large.
*        ---------------------------------------------------------------
*        First, compute the residuals for all the constraints not in the
*        working set.

         if (nclin .gt. 0  .and.  nactiv+nfixed .lt. n) then
            do 410, i = 1, nclin
               if (istate(n+i) .le. 0)
     $         Ax(i) = ddot  (n, A(i,1), ldA, wx, 1 )
  410       continue

            is     = 1
            toobig = tolact + tolact

*+          while (is .gt. 0  .and.  nfixed + nactiv .lt. n) do
  500       if    (is .gt. 0  .and.  nfixed + nactiv .lt. n) then
               is     = 0
               resmin = tolact

               do 520, i = 1, nclin
                  j      = n + i
                  if (istate(j) .eq. 0) then
                     b1     = bl(j)
                     b2     = bu(j)
                     resl   = toobig
                     resu   = toobig
                     if (b1 .gt. biglow)
     $                  resl  = abs( Ax(i) - b1 ) / (one + abs( b1 ))
                     if (b2 .lt. bigupp)
     $                  resu  = abs( Ax(i) - b2 ) / (one + abs( b2 ))
                     residl   = min( resl, resu )
                     if(residl .lt. resmin) then
                        resmin = residl
                        imin   = i
                        is     = 1
                        if (resl .gt. resu) is = 2
                     end if
                  end if
  520          continue

               if (is .gt. 0) then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j         = n + imin
                  istate(j) = is
               end if
               go to 500
*+          end while
            end if
         end if
      end if
            
      if (vertex  .and.  nactiv+nfixed .lt. n) then
*        ---------------------------------------------------------------
*        Find an initial vertex by temporarily fixing some variables.
*        ---------------------------------------------------------------
*        Compute lengths of columns of selected linear constraints
*        (just the ones corresponding to variables eligible to be
*        temporarily fixed).        

         do 630, j = 1, n
            if (istate(j) .eq. 0) then
               colsiz = zero
               do 620, k = 1, nclin
                  if (istate(n+k) .gt. 0)
     $            colsiz = colsiz + abs( A(k,j) )
  620          continue
               work(j) = colsiz
            end if
  630    continue
         
*        Find the  nartif  smallest such columns.
*        This is an expensive loop.  Later we can replace it by a
*        4-pass process (say), accepting the first col that is within
*        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).
*        (This comment written in 1980).
         
*+       while (nfixed + nactiv .lt. n) do
  640    if    (nfixed + nactiv .lt. n) then
            colmin = flmax
            do 650, j = 1, n
               if (istate(j) .eq. 0) then
                  if (nclin .eq. 0) go to 660
                  colsiz = work(j)
                  if (colmin .gt. colsiz) then
                     colmin = colsiz
                     jmin   = j
                  end if
               end if
  650       continue
            j         = jmin
  660       istate(j) = 4
            nartif    = nartif + 1
            nfixed    = nfixed + 1
            go to 640
*+       end while
         end if
      end if
      
      nfree = n - nfixed

      if (lsdbg) then
         if (ilsdbg(1) .gt. 0)
     $       write(iPrint, 1300) nfixed, nactiv, nartif
         if (ilsdbg(2) .gt. 0)
     $       write(iPrint, 1200) (wx(j), j = 1, n)
      end if

      return

 1000 format(/ ' //lscrsh// cold nclin nctotl'
     $       / ' //lscrsh// ', l4, i6, i7 )
 1100 format(/ ' //lscrsh// Variables before crash... '/ (5g12.3))
 1200 format(/ ' //lscrsh// Variables after  crash... '/ (5g12.3))
 1300 format(/ ' //lscrsh// Working set selected ...             '
     $       / ' //lscrsh// nfixed nactiv nartif      '
     $       / ' //lscrsh// ', i6, 2i7 )

*     end of lscrsh
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsdel ( unitQ,
     $                   n, nactiv, nfree, nres, ngq, nZ, nZr,
     $                   ldA, ldQ, ldR, ldT, nrank,
     $                   jdel, kdel, kactiv, kx,
     $                   A, res, R, T, gq, Q,
     $                   c, s )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kactiv(n), kx(n)
      double precision   A(ldA,*), res(n,*), R(ldR,*), T(ldT,*),
     $                   gq(n,*), Q(ldQ,*)
      double precision   c(n), s(n)

*     ==================================================================
*     lsdel   updates the least-squares factor R and the factorization
*     A(free) (Z Y) = (0 T) when a regular, temporary or artificial
*     constraint is deleted from the working set.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level-2 matrix routines added 25-Apr-1988.
*     This version of lsdel dated  10-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol5cm/ Asize, dTmax, dTmin

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      
      intrinsic          max   , min
      external           idamax
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      if (jdel .gt. 0) then
*        ---------------------------------------------------------------
*        Regular constraint or temporary bound deleted.
*        ---------------------------------------------------------------

         if (jdel .le. n) then

*           Case 1.  A simple bound has been deleted.
*           =======  Columns nfree+1 and ir of r must be swapped.

            ir     = nZ    + kdel
            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $         write(iPrint, 1100) nactiv, nZ, nfree, ir, jdel, unitQ

            itdel  = 1
            nfree  = nfree + 1

            if (nfree .lt. ir) then
               kx(ir)    = kx(nfree)
               kx(nfree) = jdel
               if (nrank .gt. 0)
     $            call cmrswp( n, nres, nrank, ldR, nfree, ir,
     $                         R, res, c, s )
               call dswap ( ngq, gq(nfree,1), n, gq(ir,1), n )
            end if

            if (.not. unitQ) then

*              Copy the incoming column of  A(free)  into the end of T.

               do 130, ka = 1, nactiv
                  i = kactiv(ka)
                  T(ka,nfree) = A(i,jdel)
  130          continue

*              Expand Q by adding a unit row and column.

               if (nfree .gt. 1) then
                  call dload ( nfree-1, zero, Q(nfree,1), ldQ )
                  call dload ( nfree-1, zero, Q(1,nfree), 1  )
               end if
               Q(nfree,nfree) = one
            end if
         else        

*           Case 2.  A general constraint has been deleted.
*           =======

            if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $         write(iPrint, 1200) nactiv, nZ, nfree, kdel, jdel, unitQ

            itdel  = kdel
            nactiv = nactiv - 1

*           Delete row  kdel  of T and move up the ones below it.
*           T becomes reverse lower Hessenberg.

            do 220, i = kdel, nactiv
               kactiv(i) = kactiv(i+1)
               ld        = nfree - i
               call dcopy ( i+1, T(i+1,ld), ldT, T(i,ld), ldT )
  220       continue
         end if

         nZ    = nZ     + 1

         if (nactiv .eq. 0) then
            dTmax = one
            dTmin = one
         else
*           ------------------------------------------------------------
*           Restore the nactiv by (nactiv+1) reverse-Hessenberg matrix 
*           T  to reverse-triangular form.  The last  nactiv  super-
*           diagonal elements are removed using a backward sweep of
*           plane rotations.  The rotation for the singleton in the 
*           first column is generated separately.
*           ------------------------------------------------------------
            nsup   = nactiv - itdel + 1
                                              
            if (nsup .gt. 0) then
               npiv   = nfree  - itdel + 1

               if (nsup .gt. 1) then
                  call dcopy ( nsup-1, T(nactiv-1,nZ+1), ldT-1, 
     $                         s(nZ+1), 1)
                  call f06qzf( 'Remove', nactiv, 1, nsup, 
     $                         c(nZ+1), s(nZ+1), T(1,nZ+1), ldT )
               end if

               call f06baf( T(nactiv,nZ+1), T(nactiv,nZ), cs, sn )
               T(nactiv,nZ) = zero
               s(nZ)   = - sn
               c(nZ)   =   cs
         
               call f06qxf( 'Right', 'Variable', 'Backwards', 
     $                      nfree, nfree, nZ, npiv, c, s, Q, ldQ )
               call f06qxf( 'Left ', 'Variable', 'Backwards', 
     $                      npiv , ngq  , nZ, npiv, c, s, gq, n    )
            
               nT = min( nrank, npiv )
               
               if (nT .lt. npiv  .and.  nT .gt. 0) then
               
*                 R is upper trapezoidal, pretend R is (nt x n) and 
*                 apply the rotations in columns  max(nt,nZ)  thru npiv.
               
                  call f06qxf( 'Right', 'Variable', 'Backwards', nT, n,
     $                         max(nt,nZ), npiv, c, s, R, ldR )
               end if
               
*              Apply the column transformations to the triangular part 
*              of  R.  The arrays  c  and  s  containing the column
*              rotations are overwritten by the row rotations that
*              restore  R  to upper-triangular form.
               
               if (nZ .lt. nT) then
                  call f06qtf( 'Right', nT, nZ, nT, c, s, R, ldR )
               end if
               
*              Apply the row rotations to the remaining rows of R.
               
               if (n .gt. nT)
     $            call f06qxf( 'Left', 'Variable', 'Backwards', 
     $                         nT, n-nT, nZ, nT, c, s, R(1,nT+1), ldR )
               
               if (nres .gt. 0)
     $            call f06qxf( 'Left', 'Variable', 'Backwards', 
     $                         nT, nres, nZ, nT, c, s, res, n )
            end if
            
            call dcond ( nactiv, T(nactiv,nZ+1), ldT-1, dTmax, dTmin )
         end if
      end if

      nZr1 = nZr + 1

      if (nZ .gt. nZr) then
         if (jdel .gt. 0) then
            jart =   nZr1 - 1 + idamax( nZ-nZr1+1, gq(nZr1,1), 1 )
         else
            jart = - jdel
         end if

         if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $      write( iPrint, 1000 ) nZ, nZr1, jart

         if (jart .gt. nZr1) then

*           Swap columns NZR1 and JART of R.

            if (unitQ) then
               k        = kx(nZr1)
               kx(nZr1) = kx(jart)
               kx(jart) = k
            else
               call dswap ( nfree, Q(1,nZr1), 1, Q(1,jart), 1 )
            end if

            call dswap ( ngq, gq(nZr1,1), n, gq(jart,1), n )
            if (nrank .gt. 0)
     $         call cmrswp( n, nres, nrank, ldR, nZr1, jart,
     $                      R, res, c, s )
         end if
      end if

      nZr = nZr1

      return

 1000 format(/ ' //lsdel //  Artificial constraint deleted.      '
     $       / ' //lsdel //      nZ   nZr   jart                 '
     $       / ' //lsdel //  ', 3i6 )
 1100 format(/ ' //lsdel //  Simple bound deleted.               '
     $       / ' //lsdel //  nactiv    nZ nfree    ir  jdel unitQ'
     $       / ' //lsdel //  ', 5i6, l6 )
 1200 format(/ ' //lsdel //  General constraint deleted.         '
     $       / ' //lsdel //  nactiv    nZ nfree  kdel  jdel unitQ'
     $       / ' //lsdel //  ', 5i6, l6 )

*     end of lsdel
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsdflt( m, n, nclin, title )

      implicit           double precision(a-h,o-z)

      character*(*)      title

*     ==================================================================
*     lsdflt  loads the default values of the parameters not set by 
*     the user.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original Fortran 77 version written 17-September-1985.
*     This version of lsdflt dated  21-Mar-93.
*     ==================================================================
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol4cm/ epspt3, epspt5, epspt8, epspt9

      logical            cmdbg, lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg
      common    /cmdebg/ icmdbg(ldbg), cmdbg

      logical            newopt
      common    /sol3ls/ newopt
      save      /sol3ls/

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
      equivalence   (msgls , msglvl), (idbgls, idbg), (ldbgls, msgdbg)

      character*4        icrsh(0:2)
      character*3        lstype(1:10)
      character*16       key
      parameter         (zero   =  0.0d+0, ten    = 10.0d+0)
      parameter         (hundrd =100.0d+0)
      parameter         (rdummy = -11111., idummy = -11111)
      parameter         (gigant = 1.0d+20*.99999          )
      parameter         (wrktol = 1.0d-2                  )
      data               icrsh(0), icrsh(1), icrsh(2)
     $                 /'Cold'   ,'Warm'   ,'Hot '   /
      data               lstype(1), lstype(2)
     $                 /' FP'     ,' LP'     /
      data               lstype(3), lstype(4), lstype(5), lstype(6)
     $                 /'QP1'     ,'QP2'     ,'QP3'     ,'QP4'     /
      data               lstype(7), lstype(8), lstype(9), lstype(10)
     $                 /'LS1'     ,'LS2'     ,'LS3'     ,'LS4'     /

      epsmch = wmach( 3)

*     Make a dummy call to lskey to ensure that the defaults are set.

      call lskey ( nout, '*', key )
      newopt = .true.

*     Save the optional parameters set by the user.  The values in
*     rprmls and iprmls may be changed to their default values.

      call icopy ( mxparm, iprmls, 1, ipsvls, 1 )
      call dcopy ( mxparm, rprmls, 1, rpsvls, 1 )

      if (       iPrnt  .lt. 0      )  iPrnt   = nout
      if (       iSumry .lt. 0      )  iSumry  = 6
      if (       iSumry .eq. iPrnt  )  iSumry  = 0
                                       iPrint  = iPrnt
                                       iSumm   = iSumry
      if (       lprob  .lt. 0      )  lprob   = 7
      if (       lcrash .lt. 0
     $    .or.   lcrash .gt. 2      )  lcrash  = 0
      if (       itmax1 .lt. 0      )  itmax1  = max(50, 5*(n+nclin))
      if (       itmax2 .lt. 0      )  itmax2  = max(50, 5*(n+nclin))
      if (       msglvl .eq. idummy )  msglvl  = 10
      if (       idbg   .lt. 0
     $    .or.   idbg   .gt. itmax1 + itmax2
     $                              )  idbg    = 0
      if (       msgdbg .lt. 0      )  msgdbg  = 0
      if (       msgdbg .eq. 0      )  idbg    = itmax1 + itmax2 + 1
      if (       tolact .lt. zero   )  tolact  = wrktol
      if (       tolfea .eq. rdummy
     $    .or.  (tolfea .ge. zero
     $    .and.  tolfea .lt. epsmch))  tolfea  = epspt5
      if (       tolrnk .le. zero
     $    .and. (lprob  .eq. 5  .or.
     $           lprob  .eq. 7  .or.
     $           lprob  .eq. 9)     )  tolrnk  = hundrd*epsmch
      if (       tolrnk .le. zero   )  tolrnk  =    ten*epspt5
      if (       bigbnd .le. zero   )  bigbnd  = gigant
      if (       bigdx  .le. zero   )  bigdx   = max(gigant, bigbnd)

      lsdbg = idbg .eq. 0  .and.  iPrint .gt. 0
      cmdbg = lsdbg
      k     = 1
      msg   = msgdbg
      do 200, i = 1, ldbg
         ilsdbg(i) = mod( msg/k, 10 )
         icmdbg(i) = ilsdbg(i)
         k = k*10
  200 continue

      if (msglvl .gt. 0) then

*        Print the title.

         lenT = len( title )
         if (lenT .gt. 0) then
            nspace = (81 - lenT)/2 + 1
            if (iPrint .gt. 0) then
               write(iPrint, '(///// (80a1) )')
     $            (' ', j=1, nspace), (title(j:j), j=1,lenT)
               write(iPrint, '(80a1 //)')
     $            (' ', j=1, nspace), ('='       , j=1,lenT)
            end if

            if (iSumm .gt. 0) then
               write(iSumm , '(///// (80a1) )')
     $            (' ', j=1, nspace), (title(j:j), j=1,lenT)
               write(iSumm , '(80a1 //)')
     $            (' ', j=1, nspace), ('='       , j=1,lenT)
            end if
         end if

         if (iPrint .gt. 0) then
            write(iPrint, 2000)
            write(iPrint, 2100) lstype(lprob),
     $                          nclin , tolfea, icrsh(lcrash),
     $                          n     , bigbnd, tolact,
     $                          m     , bigdx , tolrnk
            write(iPrint, 2200) msglvl, iPrnt , itmax1, 
     $                          epsmch, iSumry, itmax2
         end if
      end if

      return

 2000 format(
     $//' Parameters'
     $/ ' ----------' )
 2100 format(
     $/ ' Problem type...........', 7x, a3
     $/ ' Linear constraints.....',     i10,   6x,
     $  ' Feasibility tolerance..', 1p, e10.2, 6x,
     $  1x, a4, ' start.............'
     $/ ' Variables..............',     i10,   6x,
     $  ' Infinite bound size....', 1p, e10.2, 6x,
     $  ' Crash tolerance........',     e10.2
     $/ ' Objective matrix rows..',     i10,   6x,
     $  ' Infinite step size.....', 1p, e10.2, 6x,
     $  ' Rank tolerance.........',     e10.2 )
 2200 format(
     $/ ' Print level............',     i10,   6x,
     $  ' Print file.............',     i10,   6x,
     $  ' Feasibility phase itns.',     i10
     $/ ' eps (machine precision)', 1p, e10.2, 6x,
     $  ' Summary file...........',     i10,   6x,
     $  ' Optimality  phase itns.',     i10 )

*     end of lsdflt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsfeas( n, nclin, istate,
     $                   bigbnd, cvnorm, errmax, jmax, nviol,
     $                   Ax, bl, bu, featol, x, work )

      implicit           double precision(a-h,o-z)
      integer            istate(n+nclin)
      double precision   Ax(*), bl(n+nclin), bu(n+nclin)
      double precision   featol(n+nclin), x(n)
      double precision   work(n+nclin)

*     ==================================================================
*     lsfeas  computes the following...
*     (1)  The number of constraints that are violated by more
*          than  featol  and the 2-norm of the constraint violations.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version      April    1984.
*     This version of  lsfeas  dated  17-October-1985.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           idamax, dnrm2
      intrinsic          abs
      parameter        ( zero = 0.0d+0 )

      biglow = - bigbnd
      bigupp =   bigbnd

*     ==================================================================
*     Compute nviol,  the number of constraints violated by more than
*     featol,  and cvnorm,  the 2-norm of the constraint violations and
*     residuals of the constraints in the working set.
*     ==================================================================
      nviol  = 0

      do 200, j = 1, n+nclin
         feasj  = featol(j)
         is     = istate(j)
         res    = zero

         if (is .ge. 0  .and.  is .lt. 4) then
            if (j .le. n) then
               con =  x(j)
            else
               i   = j - n
               con = Ax(i)
            end if

            tolj   = feasj

*           Check for constraint violations.

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

*           This constraint is satisfied,  but count the residual as a
*           violation if the constraint is in the working set.

            if (is .le. 0) res = zero
            if (is .eq. 1) res = bl(j) - con
            if (is .ge. 2) res = bu(j) - con
            if (abs( res ) .gt. feasj) nviol = nviol + 1
         end if
  190    work(j) = res
  200 continue

      jmax   = idamax( n+nclin, work, 1 )
      errmax = abs ( work(jmax) )

      if (lsdbg  .and.  ilsdbg(1) .gt. 0)
     $   write(iPrint, 1000) errmax, jmax

      cvnorm  = dnrm2 ( n+nclin, work, 1 )

      return

 1000 format(/ ' //lsfeas//  The maximum violation is ', 1pe14.2,
     $                     ' in constraint', i5 )

*     end of lsfeas
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsfile( ioptns, inform )
      integer            ioptns, inform

*     ==================================================================
*     lsfile  reads the options file from unit  ioptns  and loads the
*     options into the relevant elements of  iprmls  and  rprmls.
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
      common     /sol3ls/ newopt
      save       /sol3ls/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar, lskey
      logical             first
      save                first , nout
      data                first /.true./

*     If first time in, set nout.
*     newopt is true first time into lsfile or lsoptn
*     and just after a call to lssol.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout = wmach(11)
      end if

      call opfile( ioptns, nout, inform, lskey )


*     end of lsfile
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsgetp( linObj, singlr, unitgZ, unitQ,
     $                   n, nclin, nfree,
     $                   ldA, ldQ, ldR, nrank, numinf, nZr,
     $                   kx, ctp, pnorm,
     $                   A, Ap, res, hZ, p,
     $                   gq, cq, R, Q, work )

      implicit           double precision(a-h,o-z)
      logical            linObj, singlr, unitgZ, unitQ
      integer            kx(n)
      double precision   A(ldA,*), Ap(*), res(*), hZ(*), p(n),
     $                   gq(n), cq(*), R(ldR,*), Q(ldQ,*)
      double precision   work(n)

*     ==================================================================
*     lsgetp  computes the following quantities for  lscore.
*     (1) The vector  (hZ) = (Rz)(pz).
*         If X is not yet feasible,  the product is computed directly.
*         If  Rz is singular,  hZ  is zero.  Otherwise  hZ  satisfies
*         the equations
*                        Rz'hZ = -gz,
*         where  g  is the total gradient.  If there is no linear term
*         in the objective,  hZ  is set to  dz  directly.
*     (2) The search direction P (and its 2-norm).  The vector P is
*         defined as  Z*(pz), where  (pz)  depends upon whether or
*         not x is feasible and the nonsingularity of  (Rz).
*         If  numinf .gt. 0,  (pz)  is the steepest-descent direction.
*         Otherwise,  x  is the solution of the  nZr*nZr  triangular
*         system   (Rz)*(pz) = (hZ).
*     (3) The vector Ap,  where A is the matrix of linear constraints.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of lsgetp dated  23-Oct-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           ddot  , dnrm2
      intrinsic          min
      parameter        ( zero = 0.0d+0, one  = 1.0d+0 )

      if (singlr) then
*        ---------------------------------------------------------------
*        The triangular factor for the current objective function is
*        singular,  i.e., the objective is linear along the last column
*        of Zr.  This can only occur when unitgZ is true.
*        ---------------------------------------------------------------
         if (nZr .gt. 1) then
            call dcopy ( nZr-1, R(1,nZr), 1, p, 1 )
            call dtrsv ( 'U', 'N', 'N', nZr-1, R, ldR, p, 1 )
         end if
         p(nZr) = - one

         gtp = ddot  ( nZr, gq, 1, p, 1 )
         if (gtp .gt. zero) call dscal ( nZr, (-one), p, 1 )

         if (nZr .le. nrank) then
            if (numinf .eq. 0) then
               if (unitgZ) then
                  hZ(nZr) = R(nZr,nZr)*p(nZr)
               else
                  call dload ( nZr, (zero), hZ, 1 )
               end if
            else
               hZ(1)   = R(1,1)*p(1)
            end if
         end if
      else
*        ---------------------------------------------------------------
*        The objective is quadratic in the space spanned by Zr.
*        ---------------------------------------------------------------
         if (linObj) then
            if (unitgZ) then
               if (nZr .gt. 1)
     $            call dload ( nZr-1, (zero), hZ, 1 )
               hZ(nZr) = - gq(nZr)/R(nZr,nZr)
            else
               call dcopy ( nZr, gq  , 1, hZ, 1 )
               call dscal ( nZr, (-one), hZ, 1 )
               call dtrsv ( 'U', 'T', 'N', nZr, R, ldR, hZ, 1 )
            end if
         else
            call dcopy ( nZr, res, 1, hZ, 1 )
         end if

*        Solve  Rz*pz = hZ.

         call dcopy ( nZr, hZ, 1, p, 1 )
         call dtrsv ( 'U', 'N', 'N', nZr, R, ldR, p, 1 )
      end if

*     Compute  p = Zr*pz  and its norm.

      if (linObj)
     $   ctp = ddot  ( nZr, cq, 1, p, 1 )
      pnorm  = dnrm2 ( nZr, p, 1 )

      call cmqmul( 1, n, nZr, nfree, ldQ, unitQ, kx, p, Q, work )

      if (lsdbg  .and.  ilsdbg(2) .gt. 0)
     $   write(iPrint, 1000) (p(j), j = 1, n)

*     Compute  Ap.

      if (nclin .gt. 0) then
         call dgemv ( 'No transpose', nclin, n, one, A, ldA,
     $                p, 1, zero, Ap, 1 )

         if (lsdbg  .and.  ilsdbg(2) .gt. 0)
     $      write(iPrint, 1100) (ap(i), i = 1, nclin)
      end if

      return

 1000 format(/ ' //lsgetp//   p ... ' / (1p, 5e15.5))
 1100 format(/ ' //lsgetp//  Ap ... ' / (1p, 5e15.5))

*     end of lsgetp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsgset( prbtyp, linObj, singlr, unitgZ, unitQ,
     $                   n, nclin, nfree,
     $                   ldA, ldQ, ldR, nrank, nZ, nZr,
     $                   istate, kx,
     $                   bigbnd, tolrnk, numinf, suminf,
     $                   bl, bu, A, res, featol,
     $                   gq, cq, R, x, wtinf, Q, wrk )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      logical            linObj, singlr, unitgZ, unitQ
      integer            istate(*), kx(n)
      double precision   bl(*), bu(*), A(ldA,*),
     $                   res(*), featol(*)
      double precision   gq(n), cq(*), R(ldR,*), x(n), wtinf(*),
     $                   Q(ldQ,*)
      double precision   wrk(n)

*     ==================================================================
*     lsgset  finds the number and weighted sum of infeasibilities for
*     the bounds and linear constraints.   An appropriate transformed
*     gradient vector is returned in  GQ.
*
*     Positive values of  ISTATE(j)  will not be altered.  These mean
*     the following...
*
*               1             2           3
*           a'x = bl      a'x = bu     bl = bu
*
*     Other values of  ISTATE(j)  will be reset as follows...
*           a'x lt bl     a'x gt bu     a'x free
*              - 2           - 1           0
*
*     If  x  is feasible,  LSGSET computes the vector Q(free)'g(free),
*     where  g  is the gradient of the the sum of squares plus the
*     linear term.  The matrix Q is of the form
*                    ( Q(free)  0       ),
*                    (   0      I(fixed))
*     where  Q(free)  is the orthogonal factor of  A(free)  and  A  is
*     the matrix of constraints in the working set.  The transformed
*     gradients are stored in GQ.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     Level 2 Blas added 11-June-1986.
*     This version of lsgset dated 14-Sep-92.
*     ==================================================================
      external           ddot  , idrank
      intrinsic          abs   , max   , min
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )
                     
      bigupp =   bigbnd
      biglow = - bigbnd

      numinf =   0
      suminf =   zero
      call dload ( n, zero, gq, 1 )

      do 200, j = 1, n+nclin
         if (istate(j) .le. 0) then
            feasj  = featol(j)
            if (j .le. n) then
               ctx = x(j)
            else
               k   = j - n
               ctx = ddot  ( n, A(k,1), ldA, x, 1 )
            end if
            istate(j) = 0

*           See if the lower bound is violated.

            if (bl(j) .gt. biglow) then
               s = bl(j) - ctx
               if (s     .gt. feasj ) then
                  istate(j) = - 2
                  weight    = - wtinf(j)
                  go to 160
               end if
            end if

*           See if the upper bound is violated.

            if (bu(j) .ge. bigupp) go to 200
            s = ctx - bu(j)
            if (s     .le. feasj ) go to 200
            istate(j) = - 1
            weight    =   wtinf(j)

*           Add the infeasibility.

  160       numinf = numinf + 1
            suminf = suminf + abs( weight ) * s
            if (j .le. n) then
               gq(j) = weight
            else
               call daxpy ( n, weight, A(k,1), ldA, gq, 1 )
            end if
         end if
  200 continue

*     ------------------------------------------------------------------
*     Install  gq,  the transformed gradient.
*     ------------------------------------------------------------------
      singlr = .false.
      unitgZ = .true.

      if (numinf .gt. 0) then
         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ, kx, gq, Q, wrk )
      else if (numinf .eq. 0  .and.  prbtyp .eq. 'FP') then
         call dload ( n, zero, gq, 1 )
      else

*        Ready for the Optimality Phase.
*        Set nZr so that Rz is nonsingular.

         if (nrank .eq. 0) then
            if (linObj) then
               call dcopy ( n, cq, 1, gq, 1 )
            else
               call dload ( n, zero, gq, 1 )
            end if
            nZr    = 0
         else

*           Compute  gq = - R' * (transformed residual)

            call dcopy ( nrank, res, 1, gq, 1 )
            call dscal ( nrank, (-one), gq, 1 )
            call dtrmv ( 'U', 'T', 'N', nrank, R, ldR, gq, 1 )
            if (nrank .lt. n)
     $         call dgemv( 'T', nrank, n-nrank, -one,R(1,nrank+1),ldR,
     $                      res, 1, zero, gq(nrank+1), 1 )

            if (linObj) call daxpy ( n, one, cq, 1, gq, 1 )
            unitgZ = .false.

            rownrm = dnrm2 ( n, R(1,1), ldR )
            if (         rownrm  .le.        tolrnk
     $          .or. abs(R(1,1)) .le. rownrm*tolrnk) then
               nZr = 0
            else
               nZr = idrank( min(nrank, nZ), R, ldR+1, tolrnk )
            end if
         end if
      end if

*     end of lsgset
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lskey ( nout, buffer, key )

      implicit           double precision(a-h,o-z)
      character*(*)      buffer

*     ==================================================================
*     lskey   decodes the option contained in  buffer  in order to set
*     a parameter value in the relevant element of  iprmls  or  rprmls.
*
*
*     Input:
*
*     nout   a unit number for printing error messages.
*            if nout = 0 no error messages are printed.
*
*     Output:
*
*     key    The first keyword contained in BUFFER.
*
*
*     lskey  calls opnumb and the subprograms
*                 lookup, scannr, tokens, upcase
*     (now called oplook, opscan, optokn, opuppr)
*     supplied by Informatics General, Inc., Palo Alto, California.
*
*     Systems Optimization Laboratory, Stanford University.
*     This version of  lskey dated 19-Oct-92.
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
*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      external           opnumb
      logical            first , more  , number, opnumb, sorted
      save               first

      parameter         (     maxkey = 28,  maxtie = 12,   maxtok = 10,
     $                        maxtyp = 16)
      character*16       keys(maxkey), ties(maxtie), token(maxtok),
     $                   type(maxtyp)
      character*16       key, key2, key3, value

      parameter         (idummy = -11111,  rdummy = -11111.0,
     $                   sorted = .true.,  zero   =  0.0     )

      data                first
     $                  /.true./
      data   keys
     $ / 'BEGIN           ',
     $   'COLD            ', 'CONSTRAINTS     ', 'CRASH           ',
     $   'DEBUG           ', 'DEFAULTS        ', 'END             ',
     $   'FEASIBILITY     ', 'HOT             ', 'INFINITE        ',
     $   'IPRMLS          ', 'ITERATIONS      ', 'ITERS:ITERATIONS',
     $   'ITNS :ITERATIONS', 'LINEAR          ', 'LIST            ',
     $   'LOWER           ', 'NOLIST          ', 'OPTIMALITY      ',
     $   'PRINT           ', 'PROBLEM         ', 'RANK            ',
     $   'RPRMLS          ', 'START           ', 'SUMMARY         ',
     $   'UPPER           ', 'VARIABLES       ', 'WARM            '/

      data   ties
     $ / 'BOUND           ', 'CONSTRAINTS     ', 'FILE            ',
     $   'LEVEL           ',
     $   'NO              ', 'NO.      :NUMBER', 'NUMBER          ',
     $   'PHASE           ', 'STEP            ',
     $   'TOLERANCE       ', 'TYPE            ', 'YES             '/

      data   type
     $ / 'FP              ',
     $   'LEAST       :LS1', 'LINEAR       :LP', 'LP              ',
     $   'LS          :LS1', 'LS1             ', 'LS2             ',
     $   'LS3             ', 'LS4             ', 'LSQ         :LS1',
     $   'QP          :QP2', 'QP1             ', 'QP2             ',
     $   'QP3             ', 'QP4             ', 'QUADRATIC   :QP2'/
*-----------------------------------------------------------------------

      if (first) then
         first  = .false.
         do 10, i = 1, mxparm
            iprmls(i) = idummy
            rprmls(i) = rdummy
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
*     LOCi   says where the keywords are in the dictionaries.
*     LOCi = 0 signals that the keyword wasn't there.

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
         if (key .eq. 'COLD        ') then
            lcrash = 0
         else if (key .eq. 'CONSTRAINTS ') then
            nnclin = rvalue
         else if (key .eq. 'CRASH       ') then
            tolact = rvalue
         else if (key .eq. 'DEBUG       ') then
            ldbgls = rvalue
         else if (key .eq. 'DEFAULTS    ') then
            do 20, i = 1, mxparm
               iprmls(i) = idummy
               rprmls(i) = rdummy
   20       continue
         else if (key .eq. 'FEASIBILITY ') then
              if (key2.eq. 'PHASE       ') itmax1 = rvalue
              if (key2.eq. 'TOLERANCE   ') tolfea = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if (key .eq. 'HOT         ') then
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
            itmax2 = rvalue
         else if (key .eq. 'LINEAR      ') then
            nnclin = rvalue
         else if (key .eq. 'LOWER       ') then
            bndlow = rvalue
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'OPTIMALITY  ') then
            itmax2 = rvalue
         else if (key .eq. 'PROBLEM     ') then
            if      (key2 .eq. 'NUMBER') then
               nprob  = rvalue
            else if (key2 .eq. 'TYPE  ') then

*              Recognize     Problem type = LP     etc.

               call oplook( maxtyp, type, sorted, key3, loc3 )
               if (key3 .eq. 'FP' ) lprob = 1
               if (key3 .eq. 'LP' ) lprob = 2
               if (key3 .eq. 'QP1') lprob = 3
               if (key3 .eq. 'QP2') lprob = 4
               if (key3 .eq. 'QP3') lprob = 5
               if (key3 .eq. 'QP4') lprob = 6
               if (key3 .eq. 'LS1') lprob = 7
               if (key3 .eq. 'LS2') lprob = 8
               if (key3 .eq. 'LS3') lprob = 9
               if (key3 .eq. 'LS4') lprob = 10
               if (loc3 .eq.  0   ) then
                  if (nout .gt. 0)  write(nout, 2330) key3
               end if
            else
               if (nout .gt. 0) write(nout, 2320) key2
            end if
         else
            more   = .true.
         end if
      end if

      if (more) then
         more   = .false.
         if      (key .eq. 'PRINT       ') then
              if (key2.eq. 'FILE        ') iPrnt  = rvalue
              if (key2.eq. 'LEVEL       ') msgls  = rvalue
              if (loc2.eq.  0            ) then
                 if (nout .gt. 0)          write(nout, 2320) key2
              end if
         else if (key .eq. 'RANK        ') then
            tolrnk = rvalue
         else if (key .eq. 'RPRMLS      ') then
*           Allow things like  rprmls 21 = 2  to set rprmls(21) = 2.0
            ivalue = rvalue
            if (ivalue .ge. 1  .and. ivalue .le. mxparm) then
               read (key3, '(bn, e16.0)') rprmls(ivalue)
            else
               if (nout .gt. 0) write(nout, 2400) ivalue
            end if
         else if (key .eq. 'START       ') then
            idbgls = rvalue
         else if (key .eq. 'SUMMARY     ') then
            iSumry = rvalue
         else if (key .eq. 'UPPER       ') then
            bndupp = rvalue
         else if (key .eq. 'VARIABLES   ') then
            nn     = rvalue
         else if (key .eq. 'WARM        ') then
            lcrash = 1
         else
            if (nout .gt. 0) write(nout, 2300) key
         end if
      end if

  900 return

 2300 format(' XXX  Keyword not recognized:         ', a)
 2320 format(' XXX  Second keyword not recognized:  ', a)
 2330 format(' XXX  Third  keyword not recognized:  ', a)
 2400 format(' XXX  The parm subscript is out of range:', i10)

*     end of lskey
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsloc ( lprob, n, nclin, litotl, lwtotl )

      implicit           double precision(a-h,o-z)

*     ==================================================================
*     lsloc   allocates the addresses of the work arrays for  lscore.
*
*     Note that the arrays  ( gq, cq )  and  ( res, res0, hz )  lie in
*     contiguous areas of workspace.
*     res, res0 and hZ are not needed for LP.
*     CQ is defined when the objective has an explicit linear term.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written  29-October-1984.
*     This version of lsloc dated 16-February-1986.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2
      common    /sol3cm/ lennam, ldT   , ncolt, ldQ

      parameter        ( lenls = 20 )
      common    /sol1ls/ locls(lenls)

      logical            lsdbg
      parameter        ( ldbg = 5 )
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      miniw     = litotl + 1
      minw      = lwtotl + 1


*     Assign array lengths that depend upon the problem dimensions.

      if (nclin .eq. 0) then
         lenT  = 0
         lenQ = 0
      else
         lenT  = ldT *ncolt
         lenQ = ldQ*ldQ
      end if

      lencq  = 0
      if (lprob .eq. 2*(lprob/2)) lencq  = n
      lenres = 0
      if (lprob .gt. 2          ) lenres = n

      lkactv    = miniw
      miniw     = lkactv + n

      lanorm    = minw
      lAp       = lanorm + nclin
      lpx       = lAp    + nclin
      lgq       = lpx    + n
      lcq       = lgq    + n
      lres      = lcq    + lencq
      lres0     = lres   + lenres
      lhZ       = lres0  + lenres
      lrlam     = lhZ    + lenres
      lT        = lrlam  + n
      lQ        = lT     + lenT
      lwtinf    = lQ     + lenQ
      lwrk      = lwtinf + n  + nclin
      lfeatl    = lwrk   + n  + nclin
      minw      = lfeatl + n  + nclin

      locls( 1) = lkactv
      locls( 2) = lanorm
      locls( 3) = lAp
      locls( 4) = lpx
      locls( 5) = lres
      locls( 6) = lres0
      locls( 7) = lhZ
      locls( 8) = lgq
      locls( 9) = lcq
      locls(10) = lrlam
      locls(11) = lT
      locls(12) = lQ
      locls(13) = lwtinf
      locls(14) = lwrk
      locls(15) = lfeatl

      litotl    = miniw - 1
      lwtotl    = minw  - 1

*     end of lsloc
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsmove( hitcon, hitlow, linObj, unitgZ,
     $                   nclin, nrank, nZr,
     $                   n, ldR, jadd, numinf,
     $                   alfa, ctp, ctx, xnorm,
     $                   Ap, Ax, bl, bu, gq, hZ, p, res,
     $                   R, x, work )

      implicit           double precision (a-h,o-z)
      logical            hitcon, hitlow, linObj, unitgZ
      double precision   Ap(*), Ax(*), bl(*), bu(*), gq(*), hZ(*),
     $                   p(n), res(*), R(ldR,*), x(n)
      double precision   work(*)

*     ==================================================================
*     lsmove  changes x to x + alfa*p and updates ctx, Ax, res and gq
*     accordingly.
*
*     If a bound was added to the working set,  move x exactly on to it,
*     except when a negative step was taken (cmalf may have had to move
*     to some other closer constraint.)
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 27-December-1985.
*     Level 2 BLAS added 11-June-1986.
*     This version of lsmove dated 14-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           dnrm2
      intrinsic          abs   , min
      parameter        ( zero  = 0.0d+0, one = 1.0d+0 )

      call daxpy ( n, alfa, p, 1, x, 1 )
      if (linObj) ctx = ctx + alfa*ctp

      if (hitcon  .and.  jadd .le. n) then
         bnd = bu(jadd)
         if (hitlow) bnd = bl(jadd)
         if (alfa .ge. zero) x(jadd) = bnd
      end if
      xnorm  = dnrm2 ( n, x, 1 )

      if (nclin .gt. 0)
     $   call daxpy ( nclin, alfa, Ap, 1, Ax, 1 )

      if (nZr .le. nrank) then
         if (unitgZ) then
            res(nZr) = res(nZr) - alfa*hZ(nZr)
         else
            call daxpy ( nZr, (-alfa), hZ, 1, res, 1  )
         end if

         if (numinf .eq. 0) then

*           Update the transformed gradient GQ so that
*           gq = gq + alfa*R'( HZ ).
*                            ( 0  )

            if (unitgZ) then
               call daxpy ( n-nZr+1, alfa*hZ(nZr), R(nZr,nZr), ldR,
     $                                             gq(nZr)   , 1      )
            else
               call dcopy ( nZr, hZ, 1, work, 1 )
               call dtrmv ( 'U', 'T', 'N', nZr, R, ldR, work, 1 )
               if (nZr .lt. n)
     $            call dgemv ( 'T', nZr, n-nZr, one, R(1,nZr+1), ldR,
     $                         hZ, 1, zero, work(nZr+1), 1 )
               call daxpy ( n, alfa, work, 1, gq, 1 )
            end if
         end if
      end if

*     end of lsmove
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsmuls( prbtyp,
     $                   msglvl, n, nactiv, nfree,
     $                   ldA, ldT, numinf, nZ, nZr,
     $                   istate, kactiv, kx, dinky,
     $                   jsmlst, ksmlst, jinf, jtiny,
     $                   jbigst, kbigst, trulam,
     $                   A, anorms, gq, rlamda, T, wtinf )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      integer            istate(*), kactiv(n), kx(n)
      double precision   A(ldA,*), anorms(*),
     $                   gq(n), rlamda(n), T(ldT,*), wtinf(*)

*     ==================================================================
*     lsmuls  first computes the Lagrange multiplier estimates for the
*     given working set.  It then determines the values and indices of
*     certain significant multipliers.  In this process, the multipliers
*     for inequalities at their upper bounds are adjusted so that a
*     negative multiplier for an inequality constraint indicates non-
*     optimality.  All adjusted multipliers are scaled by the 2-norm
*     of the associated constraint row.  In the following, the term
*     minimum refers to the ordering of numbers on the real line,  and
*     not to their magnitude.
*
*     jsmlst  is the index of the minimum of the set of adjusted
*             multipliers with values less than  - dinky.  A negative
*             jsmlst defines the index in Q'g of the artificial
*             constraint to be deleted.
*     ksmlst  marks the position of general constraint jsmlst in kactiv.
*
*     jbigst  is the index of the largest of the set of adjusted
*             multipliers with values greater than (1 + dinky).
*     kbigst  marks its position in kactiv.
*
*     On exit,  elements 1 thru nactiv of rlamda contain the unadjusted
*     multipliers for the general constraints.  Elements nactiv onwards
*     of rlamda contain the unadjusted multipliers for the bounds.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lsmuls dated  14-Sep-92.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      intrinsic          abs, min
      parameter        ( one    =1.0d+0 )

      nfixed =   n - nfree

      jsmlst =   0
      ksmlst =   0
      smllst = - dinky

      tinylm =   dinky
      jtiny  =   0

      jbigst =   0
      kbigst =   0
      biggst =   one + dinky

      if (nZr .lt. nZ) then
*        ---------------------------------------------------------------
*        Compute jsmlst for the artificial constraints.
*        ---------------------------------------------------------------
         do 100, j = nZr+1, nZ
            rlam = - abs( gq(j) )
            if (rlam .lt. smllst) then
               smllst =   rlam
               jsmlst = - j
            else if (rlam .lt. tinylm) then
               tinylm =   rlam
               jtiny  =   j
            end if
  100    continue

         if (msglvl .ge. 20) then
            if (iPrint .gt. 0) write(iPrint, 1000) (gq(k), k=nZr+1,nZ)
         end if
      end if

*     ------------------------------------------------------------------
*     Compute jsmlst for regular constraints and temporary bounds.
*     ------------------------------------------------------------------
*     First, compute the Lagrange multipliers for the general
*     constraints in the working set, by solving  T'*lamda = Y'g.

      if (n .gt. nZ)
     $   call dcopy ( n-nZ, gq(nZ+1), 1, rlamda, 1 )
      if (nactiv .gt. 0)
     $   call cmtsol( 2, ldT, nactiv, T(1,nZ+1), rlamda )

*     -----------------------------------------------------------------
*     Now set elements nactiv, nactiv+1,... of  rlamda  equal to
*     the multipliers for the bound constraints.
*     -----------------------------------------------------------------
      do 190, l = 1, nfixed
         j     = kx(nfree+l)
         blam  = rlamda(nactiv+l)
         do 170, k = 1, nactiv
            i    = kactiv(k)
            blam = blam - A(i,j)*rlamda(k)
  170    continue
         rlamda(nactiv+l) = blam
  190 continue

*     -----------------------------------------------------------------
*     Find jsmlst and ksmlst.
*     -----------------------------------------------------------------
      do 330, k = 1, n - nZ
         if (k .gt. nactiv) then
            j = kx(nZ+k)
         else
            j = kactiv(k) + n
         end if

         is   = istate(j)

         i    = j - n
         if (j .le. n) anormj = one
         if (j .gt. n) anormj = anorms(i)

         rlam = rlamda(k)

*        Change the sign of the estimate if the constraint is in
*        the working set at its upper bound.

         if (is .eq. 2) rlam =      - rlam
         if (is .eq. 3) rlam =   abs( rlam )
         if (is .eq. 4) rlam = - abs( rlam )

         if (is .ne. 3) then
            scdlam = rlam * anormj
            if      (scdlam .lt. smllst) then
               smllst = scdlam
               jsmlst = j
               ksmlst = k
            else if (scdlam .lt. tinylm) then
               tinylm = scdlam
               jtiny  = j
            end if
         end if

         if (numinf .gt. 0  .and.  j .gt. jinf) then
            scdlam = rlam/wtinf(j)
            if (scdlam .gt. biggst) then
               biggst = scdlam
               trulam = rlamda(k)
               jbigst = j
               kbigst = k
            end if
         end if
  330 continue

*     -----------------------------------------------------------------
*     If required, print the multipliers.
*     -----------------------------------------------------------------
      if (msglvl .ge. 20  .and.  iPrint .gt. 0) then
         if (nfixed .gt. 0)
     $      write(iPrint, 1100) prbtyp, (kx(nfree+k),
     $                         rlamda(nactiv+k), k=1,nfixed)
         if (nactiv .gt. 0)
     $      write(iPrint, 1200) prbtyp, (kactiv(k),
     $                         rlamda(k), k=1,nactiv)
      end if

      if (lsdbg  .and.  ilsdbg(1) .gt. 0) then
         write(iPrint, 9000) jsmlst, smllst, ksmlst
         write(iPrint, 9100) jbigst, biggst, kbigst
         write(iPrint, 9200) jtiny , tinylm
      end if

      return

 1000 format(/ ' Multipliers for the artificial constraints        '
     $       / 4(5x, 1pe11.2))
 1100 format(/ ' Multipliers for the ', a2, ' bound  constraints   '
     $       / 4(i5, 1pe11.2))
 1200 format(/ ' Multipliers for the ', a2, ' linear constraints   '
     $       / 4(i5, 1pe11.2))
 9000 format(/ ' //lsmuls//  jsmlst     smllst     ksmlst (scaled) '
     $       / ' //lsmuls//  ', i6, 1pe11.2, 5x, i6 )
 9100 format(  ' //lsmuls//  jbigst     biggst     kbigst (scaled) '
     $       / ' //lsmuls//  ', i6, 1pe11.2, 5x, i6 )
 9200 format(  ' //lsmuls//   jtiny     tinylm                     '
     $       / ' //lsmuls//  ', i6, 1pe11.2)

*     end of lsmuls
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsoptn( string )
      character*(*)      string

*     ==================================================================
*     lsoptn  loads the option supplied in  string  into the relevant
*     element of  iprmls  or  rprmls.
*     ==================================================================

      logical             newopt
      common     /sol3ls/ newopt
      save       /sol3ls/

      double precision    wmach(15)
      common     /solmch/ wmach
      save       /solmch/

      external            mchpar
      character*16        key
      character*72        buffer
      logical             first , prnt
      save                first , nout  , prnt
      data                first /.true./

*     If first time in, set  nout.
*     newopt  is true first time into  lsfile  or  lsoptn
*     and just after a call to  lssol.
*     prnt    is set to true whenever  newopt  is true.

      if (first) then
         first  = .false.
         newopt = .true.
         call mchpar()
         nout   =  wmach(11)
      end if
      buffer = string

*     Call  lskey   to decode the option and set the parameter value.
*     If newopt is true, reset prnt and test specially for nolist.

      if (newopt) then
         newopt = .false.
         prnt   = .true.
         call lskey ( nout, buffer, key )

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

         call lskey ( iPrint, buffer, key )
         if (key .eq.   'list') prnt = .true.
         if (key .eq. 'nolist') prnt = .false.
      end if

*     end of lsoptn
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lsprt ( prbtyp, isdel, iter, jadd, jdel,
     $                   msglvl, nactiv, nfree, n, nclin,
     $                   nrank, ldR, ldT, nZ, nZr, istate,
     $                   alfa, condRz, condT, gfnorm, gZrnrm,
     $                   numinf, suminf, ctx, ssq,
     $                   Ax, R, T, x, work )

      implicit           double precision(a-h,o-z)
      character*2        prbtyp
      integer            istate(*)
      double precision   Ax(*), R(ldR,*), T(ldT,*), x(n)
      double precision   work(n)

*     ==================================================================
*     lsprt  prints various levels of output for  lscore.
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
*        ge  20        constraint status,  x  and  Ax.
*
*        ge  30        diagonals of  T  and  R.
*
*
*     Debug printing is performed depending on the logical variable 
*     lsdbg, which is set true when  idbg  major iterations have been
*     performed.  At this point,  printing is done according to a string
*     of binary digits of the form  SVT  (stored in the integer array
*     ilsdbg).
*
*     S  set 'on' gives information from the maximum step routine cmalf.
*     V  set 'on' gives various vectors in  lscore  and its auxiliaries.
*     T  set 'on' gives a trace of which routine was called and an
*                 indication of the progress of the run.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lsprt dated 07-Jan-93.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      logical            first , linObj, newSet, prtHdr
      character*2        ladd, ldel
      character*2        lstate(0:5)
      parameter         (mLine1 = 40, mLine2 = 5)
      data               lstate(0), lstate(1), lstate(2)
     $                  /'  '     , 'L '     , 'U '     /
      data               lstate(3), lstate(4), lstate(5)
     $                  /'E '     , 'F '     , 'A '     /

      if (msglvl .ge. 15) then
         if (iPrint .gt. 0) write(iPrint, 1000) prbtyp, iter
      end if

      if (msglvl .ge. 5) then

         first  = iter  .eq. 0
         linObj = nrank .eq. 0

         Itn    = mod( iter, 1000  )
         ndf    = mod( nZr , 10000 )

         nArt = nZ - nZr

         if      (jdel .gt. 0) then
            kdel =   isdel
         else if (jdel .lt. 0) then
            jdel =   nArt + 1
            kdel =   5
         else
            kdel =   0
         end if

         if (jadd .gt. 0) then
            kadd = istate(jadd)
         else
            kadd = 0
         end if

         ldel   = lstate(kdel)
         ladd   = lstate(kadd)

         if (numinf .gt. 0) then
            obj    = suminf
         else
            obj    = ssq + ctx
         end if

*        ---------------------------------------------------------------
*        If necessary, print a header. 
*        Print a single line of information.
*        ---------------------------------------------------------------
         if (iPrint .gt. 0) then
*           ------------------------------
*           Terse line for the Print file.
*           ------------------------------
            newSet = lines1 .ge. mLine1
            prtHdr = msglvl .ge. 15  .or.  first 
     $                               .or.  newSet

            if (prtHdr) then
               if (linObj) then 
                  write(iPrint, 1200)
               else
                  write(iPrint, 1300)
               end if
               lines1 = 0
            end if

            if (linObj) then
               write(iPrint, 1700) Itn, jdel, ldel, jadd, ladd,
     $                             alfa, numinf, obj, gZrnrm, ndf, nArt,
     $                             n-nfree, nactiv, gfnorm, condT
            else
               write(iPrint, 1700) Itn, jdel, ldel, jadd, ladd,
     $                             alfa, numinf, obj, gZrnrm, ndf, nArt,
     $                             n-nfree, nactiv, gfnorm, condT,
     $                             CondRz
            end if
            lines1 = lines1 + 1
         end if

         if (iSumm .gt. 0) then
*           --------------------------------
*           Terse line for the Summary file.
*           --------------------------------
            newSet = lines2 .ge. mLine2
            prtHdr =                      first 
     $                              .or.  newSet

            if (prtHdr) then
               write(iSumm , 1100)
               lines2 = 0
            end if
            write(iSumm , 1700) Itn, jdel, ldel, jadd, ladd,
     $                          alfa, numinf, obj, gZrnrm, ndf, nArt
            lines2 = lines2 + 1
         end if

         if (msglvl .ge. 20  .and.  iPrint .gt. 0) then
            write(iPrint, 2000) prbtyp
            write(iPrint, 2100) (x(j) , istate(j)  ,  j=1,n)
            if (nclin .gt. 0)
     $      write(iPrint, 2200) (Ax(k), istate(n+k), k=1,nclin )

            if (msglvl .ge. 30) then
*              ---------------------------------------------------------
*              Print the diagonals of  T  and  R.
*              ---------------------------------------------------------
               if (nactiv .gt. 0) then
                  call dcopy ( nactiv, T(nactiv,nZ+1), ldT-1, work,1 )
                  write(iPrint, 3000) prbtyp, (work(j), j=1,nactiv)
               end if
               if (nrank  .gt. 0)
     $            write(iPrint, 3100) prbtyp, (R(j,j) , j=1,nrank )
            end if
            write(iPrint, 5000)
         end if
      end if

      return

 1000 format(/// ' ', a2, ' iteration', i5
     $         / ' =================' )
 1100 format(// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     $          ' Norm gZ   Zr  Art' )
 1200 format(// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     $          ' Norm gZ   Zr  Art ',
     $          ' Bnd  Lin Norm gf  Cond T' )
 1300 format(// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     $          ' Norm gZ   Zr  Art ',
     $          ' Bnd  Lin Norm gf  Cond T Cond Rz' )
 1700 format(    i4, i5, a1, i5, a1, 1p, e8.1, i5, e16.8, 
     $           e8.1, 2i5,
     $           2i5, e8.1, 2e8.0 )
 2000 format(/ ' Values and status of the ', a2, ' constraints'
     $       / ' ---------------------------------------' )
 2100 format(/ ' Variables...'                 / (1x, 5(1p, e15.6, i5)))
 2200 format(/ ' General linear constraints...'/ (1x, 5(1p, e15.6, i5)))
 3000 format(/ ' Diagonals of ' , a2,' working set factor T'
     $       /   (1p, 5e15.6))
 3100 format(/ ' Diagonals of ' , a2, ' triangle R         '
     $       /   (1p, 5e15.6))
 5000 format(/// ' ---------------------------------------------------',
     $           '--------------------------------------------' )

*     end of lsprt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine lssetx( linObj, rowerr, unitQ,
     $                   nclin, nactiv, nfree, nrank, nZ,
     $                   n, nctotl, ldQ, ldA, ldR, ldT,
     $                   istate, kactiv, kx,
     $                   jmax, errmax, ctx, xnorm,
     $                   A, Ax, bl, bu, cq, res, res0, featol,
     $                   R, T, x, Q, p, work )

      implicit           double precision (a-h,o-z)
      logical            linObj, rowerr, unitQ
      integer            istate(nctotl), kactiv(n), kx(n)
      double precision   A(ldA,*), Ax(*), bl(nctotl), bu(nctotl),
     $                   cq(*), res(*), res0(*), featol(nctotl), p(n),
     $                   R(ldR,*), T(ldT,*), Q(ldQ,*), x(n)
      double precision   work(nctotl)

*     ==================================================================
*     lssetx  computes the point on a working set that is closest to the
*     input vector  x  (in the least-squares sense).  The norm of  x, 
*     the transformed residual vector  Pr - RQ'x,  and the constraint
*     values  Ax  are also initialized.
*
*     If the computed point gives a row error of more than the
*     feasibility tolerance, an extra step of iterative refinement is
*     used.  If  x  is still infeasible,  rowerr is set to true.
*
*     Systems Optimization Laboratory, Stanford University.
*     Original version written 31-October-1984.
*     This version of lssetx dated 29-December-1985.
*     ==================================================================
      common    /sol1cm/ nout  , iPrint, iSumm , lines1, lines2

      logical            lsdbg
      parameter         (ldbg = 5)
      common    /lsdebg/ ilsdbg(ldbg), lsdbg

      external           idamax, ddot
      intrinsic          abs, min
      parameter        ( ntry  = 2 )
      parameter        ( zero  = 0.0d+0, one = 1.0d+0 )

*     ------------------------------------------------------------------
*     Move  x  onto the simple bounds in the working set.
*     ------------------------------------------------------------------
      do 100, k = nfree+1, n
          j   = kx(k)
          is  = istate(j)
          bnd = bl(j)
          if (is .ge. 2) bnd  = bu(j)
          if (is .ne. 4) x(j) = bnd
  100 continue

*     ------------------------------------------------------------------
*     Move  x  onto the general constraints in the working set.
*     We shall make  ntry  tries at getting acceptable row errors.
*     ------------------------------------------------------------------
      ktry   = 1
      jmax   = 1
      errmax = zero

*     repeat
  200    if (nactiv .gt. 0) then

*           Set  work = residuals for constraints in the working set.
*           Solve for p, the smallest correction to x that gives a point
*           on the constraints in the working set.  Define  p = Y*(py),
*           where  py  solves the triangular system  T*(py) = residuals.

            do 220, i = 1, nactiv
               k   = kactiv(i)
               j   = n + k
               bnd = bl(j)
               if (istate(j) .eq. 2) bnd = bu(j)
               work(i) = bnd - ddot  ( n, A(k,1), ldA, x, 1 )
  220       continue

            call cmtsol( 1, ldT, nactiv, T(1,nZ+1), work )
            call dload ( n, zero, p, 1 )
            call dcopy ( nactiv, work, 1, p(nZ+1), 1 )

            call cmqmul( 2, n, nZ, nfree, ldQ, unitQ, kx, p, Q, work )
            call daxpy ( n, one, p, 1, x, 1 )
         end if

*        ---------------------------------------------------------------
*        Compute the 2-norm of  x.
*        Initialize  Ax  for all the general constraints.
*        ---------------------------------------------------------------
         xnorm  = dnrm2 ( n, x, 1 )
         if (nclin .gt. 0)
     $      call dgemv ( 'N', nclin, n, one, A, ldA,
     $                   x, 1, zero, Ax, 1 )

*        ---------------------------------------------------------------
*        Check the row residuals.
*        ---------------------------------------------------------------
         if (nactiv .gt. 0) then
            do 300, k = 1, nactiv
               i   = kactiv(k)
               j   = n + i
               is  = istate(j)
               if (is .eq. 1) work(k) = bl(j) - Ax(i)
               if (is .ge. 2) work(k) = bu(j) - Ax(i)
  300       continue

            jmax   = idamax( nactiv, work, 1 )
            errmax = abs( work(jmax) )
         end if

         ktry = ktry + 1
*     until    (errmax .le. featol(jmax) .or. ktry .gt. ntry
      if (.not.(errmax .le. featol(jmax) .or. ktry .gt. ntry)) go to 200

      rowerr = errmax .gt. featol(jmax)

*     ==================================================================
*     Compute the linear objective value  c'x  and the transformed
*     residual  Pr  -  RQ'x = res0  -  RQ'x.
*     ==================================================================
      if (nrank .gt. 0  .or.  linObj) then
         call dcopy ( n, x, 1, p, 1 )
         call cmqmul( 6, n, nZ, nfree, ldQ, unitQ, kx, p, Q, work )
      end if

      ctx = zero
      if (linObj)
     $   ctx = ddot  ( n, cq, 1, p, 1 )

      if (nrank .gt. 0) then

         call dtrmv ( 'U', 'N', 'N', nrank, R, ldR, p, 1 )
         if (nrank .lt. n)
     $      call dgemv ( 'N', nrank, n-nrank, one, R(1,nrank+1), ldR,
     $                   p(nrank+1), 1, one, p, 1 )

         call dcopy ( nrank,         res0, 1, res, 1 )
         call daxpy ( nrank, (-one), p   , 1, res, 1 )

      end if

      if (lsdbg  .and.  ilsdbg(2) .gt. 0)
     $   write(iPrint, 2200) (x(j), j = 1, n)

      return

 2200 format(/ ' //lssetx// Variables after refinement ... '/ (5g12.3))

*     end of lssetx
      end

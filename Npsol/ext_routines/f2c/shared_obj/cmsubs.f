*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  cmsubs.f
*
*     cmprnt   cmqmul   cmr1md   cmrswp
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmprnt( msglvl, n, nclin, nctotl, bigbnd,
     $                   named, names, istate,
     $                   bl, bu, clamda, featol, r )

      implicit           double precision(a-h,o-z)
      character*16       names(*)
      logical            named
      integer            istate(nctotl)
      double precision   bl(nctotl), bu(nctotl),
     $                   clamda(nctotl), featol(nctotl), r(nctotl)

*     ==================================================================
*     cmprnt  prints r(x) (x,  A*x and c(x)), the bounds, the
*     multipliers, and the slacks (distance to the nearer bound).
*
*     Original Fortran 77 version written  October 1984.
*     This version of  cmprnt dated  11-May-95.
*     ==================================================================
      common    /sol1cm/ iPrint, iSumm , lines1, lines2
      save      /sol1cm/ 

      character*1        key
      character*2        lstate(-2:4), state
      character*8        name
      character*102      line

      parameter         (zero  = 0.0d+0)
      data               lstate(-2) / '--' /, lstate(-1) / '++' /
      data               lstate( 0) / 'FR' /, lstate( 1) / 'LL' /
      data               lstate( 2) / 'UL' /, lstate( 3) / 'EQ' /
      data               lstate( 4) / 'TF' /

      if (iPrint .eq. 0
     $    .or.  (msglvl .lt. 10  .and.  msglvl .ne. 1)) return

      write(iPrint, 1000) 'Variable       '
      name   = 'variable'
      nplin  = n + nclin

      do 500, j = 1, nctotl
         b1     = bl(j)
         b2     = bu(j)
         wlam   = clamda(j)
         rj     = r(j)

         if (j .le. n) then
            number = j
         else if (j .le. nplin) then
            number = j - n
            if (number .eq. 1) then
               write(iPrint, 1000) 'Linear constrnt'
               name = 'lincon  '
            end if
         else
            number = j - nplin
            if (number .eq. 1) then
               write(iPrint, 1000) 'Nonlin constrnt'
               name = 'nlncon  '
            end if
         end if

*        Print a line for the jth variable or constraint.
*        ------------------------------------------------
         is     = istate(j)
         state  = lstate(is)
         tol    = featol(j)
         slk1   = rj - b1
         slk2   = b2 - rj
         if (abs(slk1) .lt. abs(slk2)) then
            slk = slk1
            if (b1 .le. - bigbnd) slk = slk2 
         else
            slk = slk2
            if (b2 .ge.   bigbnd) slk = slk1 
         end if

*        Flag infeasibilities, primal and dual degeneracies, 
*        and active QP constraints that are loose in NP.
*      
         key    = ' ' 
         if (slk1 .lt. -tol  .or.       slk2  .lt. -tol) key = 'I'
         if (is   .eq.  0    .and.  abs(slk ) .le.  tol) key = 'D'
         if (is   .ge.  1    .and.  abs(wlam) .le.  tol) key = 'A'

         write(line, 2000) name, number, key, state, 
     $                     rj, b1, b2, wlam, slk

*        Reset special cases:
*           Infinite bounds
*           Zero bounds
*           Lagrange multipliers for inactive constraints
*           Lagrange multipliers for infinite bounds
*           Infinite slacks
*           Zero slacks

         if (      named      ) line( 2: 17) = names(j)
         if (b1  .le. - bigbnd) line(39: 54) = '      None      '
         if (b2  .ge.   bigbnd) line(55: 70) = '      None      '
         if (b1  .eq.   zero  ) line(39: 54) = '        .       '
         if (b2  .eq.   zero  ) line(55: 70) = '        .       '
         if (is  .eq.   0       .or.    
     $       wlam.eq.   zero  ) line(71: 86) = '        .       '
         if (b1  .le. - bigbnd  .and. 
     $       b2  .ge.   bigbnd) then
                                line(71: 86) = '                '
                                line(87:102) = '                '
         end if
         if (slk .eq.   zero  ) line(87:102) = '        .       '

         write(iPrint, '(a)') line
  500 continue

      return

 1000 format(//  1x,  a15, 2x, 'State', 6x, 'Value',
     $           7x, 'Lower bound', 5x, 'Upper bound',
     $           3x, 'Lagr multiplier', 4x, '   Slack' / )
 2000 format( 1x, a8, i6, 3x, a1, 1x, a2, 4g16.7, g16.4 )

*     end of cmprnt
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cmqmul( mode, n, nZ, nfree, ldQ, unitQ,
     $                   kx, v, Q, w )

      implicit           double precision(a-h,o-z)
      logical            unitQ
      integer            kx(n)
      double precision   v(n), Q(ldQ,*), w(n)

*     ==================================================================
*     cmqmul  transforms the vector  v  in various ways using the
*     matrix  Q = ( Z  Y )  defined by the input parameters.
*
*        Mode               Result
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
      if (mode .eq. 1  .or.  mode .eq. 4) j2 = nZ
      if (mode .eq. 2  .or.  mode .eq. 5  .or.  mode .eq. 7) j1 = nZ + 1
      lenv   = j2 - j1 + 1
      if (mode .le. 3) then
*        ===============================================================
*        Mode = 1, 2  or  3.
*        ===============================================================
         if (nfree .gt. 0) call dload ( nfree, zero, w, 1 )

*        Copy  v(fixed)  into the end of wrk.

         if (mode .ge. 2  .and.  nfixed .gt. 0)
     $      call dcopy ( nfixed, v(nfree+1), 1, w(nfree+1), 1 )

*        Set  w  =  relevant part of  Qv.

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
               j      = kx(nfree+l)
               v(j)   = w (nfree+l)
  320       continue
         end if

      else
*        ===============================================================
*        Mode = 4, 5, 6, 7  or  8.
*        ===============================================================
*        Put the fixed components of  v  into the end of w.

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

*        Copy the fixed components of  w  into the end of v.

         if (nfixed .gt. 0  .and.  (mode .eq. 5  .or.  mode .eq. 6))
     $      call dcopy ( nfixed, w(nfree+1), 1, v(nfree+1), 1 )
      end if

*     end of cmqmul
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

*     end of cmr1md
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
*     This version of  cmrswp  dated  26-Aug-1991.
*     ==================================================================
      parameter        ( zero = 0.0d+0 )

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


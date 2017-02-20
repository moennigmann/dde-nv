*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  f06.f
*     A subset of the NAG F06 Chapter with some modifications.
*     The routines perform the same function as the NAG F06 routines.
*
*                         Level 0  F06  Scalar routines
*                         -------  ---- ---------------
*     f06aaz+         f06baf/drot3g+  f06bcf/dcsg+    ddiv+  
*     f06bmf/dnorm+
*
*                         Level 1  F06  Vector routines
*                         -------  ---  ---------------
*     iload/f06dbf    dload/f06fbf    dddiv           ddscl/f06fcf 
*     icopy/f06dff    dssq/f06fjf+    dcond/f06flf
*     idrank/f06klf+  f06fqf/dsrotg   dgrfg/f06frf+
*
*                         Level 2  F06  Matrix routines
*                         -------  ---  ---------------
*     f06qff          f06qgf          f06qhf
*     f06qkf          f06qnf          f06qrf          f06qsf   
*     f06qtf          f06qvf          f06qwf          f06qxf   
*     f06qzf
*
*    +Differs from the Nag F06 version.
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06AAZ( SRNAME, INFO )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.
C     .. Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*13       SRNAME
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      CHARACTER*80       REC (1)
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine f06baf( x, y, cs, sn )

      double precision   x, y, cs, sn

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
C
C  Note: f06baf/drot3g is different from the Nag routine f06baf.
C
C  f06baf  generates a plane rotation that reduces the vector (X, Y) to
C  the vector (A, 0),  where A is defined as follows...
C
C     If both X and Y are negligibly small, or
C     if Y is negligible relative to Y,
C     then  A = X,  and the identity rotation is returned.
C
C     If X is negligible relative to Y,
C     then  A = Y,  and the swap rotation is returned.
C
C     Otherwise,  A = sign(X) * sqrt( X**2 + Y**2 ).
C
C  In all cases,  X and Y are overwritten by A and 0,  and CS will lie
C  in the closed interval (0, 1).  Also,  the absolute value of CS and
C  SN (if nonzero) will be no less than the machine precision,  EPS.
C
C  DROT3G  guards against overflow and underflow.
C  It is assumed that  FLMIN .lt. EPS**2  (i.e.  RTMIN .lt. EPS).
C
C  Systems Optimization Laboratory, Stanford University.
C  Original version dated January 1982.
C  F77 version dated 28-June-1986.
C  This version of DROT3G dated 28-June-1986.
C
      double precision   a, b, eps, one, rtmin, zero
      logical            first
      intrinsic          abs, max, sqrt
      parameter        ( zero = 0.0d+0, one = 1.0d+0 )

      save               first , eps   , rtmin
      data               first / .true. /

      if( first )then
         first = .false.
         eps    = wmach(3)
         rtmin  = wmach(6)
      end if

      if (y .eq. zero) then
         cs = one
         sn = zero
      else if (x .eq. zero) then
         cs = zero
         sn = one
         x  = y
      else
         a      = abs(x)
         b      = abs(y)
         if (max(a,b) .le. rtmin) then
            cs = one
            sn = zero
         else
            if (a .ge. b) then
               if (b .le. eps*a) then
                  cs = one
                  sn = zero
                  go to 900
               else
                  a  = a * sqrt( one + (b/a)**2 )
               end if
            else
               if (a .le. eps*b) then
                  cs = zero
                  sn = one
                  x  = y
                  go to 900
               else
                  a  = b * sqrt( one + (a/b)**2 )
               end if
            end if
            if (x .lt. zero) a = - a
            cs = x/a
            sn = y/a
            x  = a
         end if
      end if

  900 y  = zero

*     end of  f06baf (drot3g).
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine f06bcf( t, c, s )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-602 (MAR 1988).

C     .. Scalar Arguments ..
      double precision   c, s, t
C     ..
C
C  F06BCF returns values c and s such that
C
C     c = cos( theta ),   s = sin( theta )
C
C  for a given value of
C
C     t = tan( theta ).
C
C  c is always non-negative and s has the same sign as t, so that
C
C     c = 1.0/sqrt( 1.0 + t**2 ),   s = t/sqrt( 1.0 + t**2 ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 28-February-1986.
C     Sven Hammarling, Nag Central Office.
C  -- Modified 19-August-1987.
C     Sven Hammarling and Jeremy Du Croz, Nag Central Office.
C        No longer sets s to zero when t is less than eps.
C  -- Modified 24-July-1991.
C     Philip E. Gill, UCSD.
C        Modified to call mchpar instead of x02ajf
C
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
C     .. Parameters ..
      double precision   one
      parameter        ( one = 1.0d+0 )
C     .. Local Scalars ..
      double precision   abst, eps, rrteps, rteps
      logical            first
C     .. External Functions ..
C+    DOUBLE PRECISION   X02AJF
C+    EXTERNAL           X02AJF
C     .. Intrinsic Functions ..
      intrinsic          abs, sign, sqrt
C     .. Save statement ..
      save               first, eps, rteps, rrteps
C     .. Data statements ..
      data               first/ .true. /
C     ..
C     .. Executable Statements ..
      if( first )then
         first  = .false.
         eps    = wmach(3)
C+       eps    = x02ajf( )
         rteps  =  sqrt( eps )
         rrteps =  1/rteps
      end if
C
      abst = abs( t )
      if( abst.lt.rteps )then
         c = one
         s = t
      else if( abst.gt.rrteps )then
         c = 1/abst
         s = sign( one, t )
      else
         c = 1/sqrt( 1 + abst**2 )
         s = c*t 
      end if
C
C     end of f06bcf. ( scsg )
      end    

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function ddiv  ( a, b, fail )

C     DOUBLE PRECISION          F06BLF
C     ENTRY                     F06BLF( A, B, FAIL )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.

C     .. Scalar Arguments ..
      double precision                  a, b
      logical                           fail
C     ..
C
C  F06BLF returns the value div given by
C
C     div = ( a/b                 if a/b does not overflow,
C           (
C           ( 0.0                 if a .eq. 0.0,
C           (
C           ( sign( a/b )*flmax   if a .ne. 0.0  and a/b would overflow,
C
C  where  flmax  is a large value, via the function name. In addition if
C  a/b would overflow then  fail is returned as true, otherwise  fail is
C  returned as false.
C
C  Note that when  a and b  are both zero, fail is returned as true, but
C  div  is returned as  0.0. In all other cases of overflow  div is such
C  that  abs( div ) = flmax.
C
C  When  b = 0  then  sign( a/b )  is taken as  sign( a ).
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 26-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision      absb, div, flmax, flmin
      logical               first
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
C     .. Intrinsic Functions ..
      intrinsic             abs, sign
C     .. Save statement ..
      save                  first, flmin, flmax
C     .. Data statements ..
      data                  first/ .true. /
C     ..
C     .. Executable Statements ..
      if( a.eq.zero )then
         div = zero
         if( b.eq.zero )then
            fail = .true.
         else
            fail = .false.
         end if
      else
C
         if( first )then
            first  = .false.
            flmin  = wmach( 5 )
            flmax  = wmach( 7 )
         end if
C
         if( b.eq.zero )then
            div  =  sign( flmax, a )
            fail = .true.
         else
            absb = abs( b )
            if( absb.ge.one )then
               fail = .false.
               if( abs( a ).ge.absb*flmin )then
                  div = a/b
               else
                  div = zero
               end if
            else
               if( abs( a ).le.absb*flmax )then
                  fail = .false.
                  div  =  a/b
               else
                  fail = .true.
                  div  = flmax
                  if( ( ( a.lt.zero ).and.( b.gt.zero ) ).or.
     $                ( ( a.gt.zero ).and.( b.lt.zero ) )     )
     $               div = -div
               end if
            end if
         end if
      end if
C
      ddiv   = div
C     F06BLF = DIV
C
C     end of f06blf. ( ddiv )
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function f06bmf( scale, ssq )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     Modified by PEG 9/25/88.

C     .. Scalar Arguments ..
      double precision                  scale, ssq
C     ..
C
C  F06BMF returns the value norm given by
C
C     norm = ( scale*sqrt( ssq ), scale*sqrt( ssq ) .lt. flmax
C            (
C            ( flmax,             scale*sqrt( ssq ) .ge. flmax
C
C  via the function name.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      double precision      flmax, norm, sqt
      logical               first
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
C     .. Intrinsic Functions ..
      INTRINSIC             SQRT
C     .. Save statement ..
      save                  first, flmax
C     .. Data statements ..
      data                  first/ .true. /
C     ..
C     .. Executable Statements ..
      if( first )then
         first = .false.
         flmax = wmach( 7 )
      end if
C
      sqt = sqrt( ssq )
      if( scale.lt.flmax/sqt )then
         norm = scale*sqt
      else
         norm = flmax
      end if
C
      f06bmf = norm
C
C     end of f06bmf. ( dnorm )
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine iload ( n, const, x, incx )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     ENTRY      F06DBF( N, CONST, X, INCX )

C     .. Scalar Arguments ..
      integer            const, incx, n
C     .. Array Arguments ..
      integer            x( * )
C     ..
C                      
C  iload/f06dbf performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 18-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( const.ne.0 )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = const
   10       continue
         else
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = 0
   20       continue
         end if
      end if
C
C     end of iload (f06dbf)
C
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dload ( n, const, x, incx )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     ENTRY      F06FBF( N, CONST, X, INCX )

C     .. Scalar Arguments ..
      double precision   const
      integer            incx, n
C     .. Array Arguments ..
      double precision   x( * )
C     ..
C
C  dload/f06fbf performs the operation
C
C     x = const*e,   e' = ( 1  1 ... 1 ).
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C                      
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( const.ne.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = const
   10       continue
         else
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
   20       continue
         end if
      end if
C
C     end of dload (f06fbf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ddscl ( n, d, incd, x, incx )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     ENTRY      F06FCF( N, D, INCD, X, INCX )

C     .. Scalar Arguments ..
      integer            incd, incx, n
C     .. Array Arguments ..
      double precision   d( * ), x( * )
C     ..
C
C  ddscl/f06fcf performs the operation
C
C     x := diag( d )*x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-September-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      integer            i, id, ix
C     .. External Subroutines ..
      external           dscal
C     .. Intrinsic Functions ..
      intrinsic          abs
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( ( incd.eq.0 ).and.( incx.ne.0 ) )then
            call dscal( n, d( 1 ), x, abs( incx ) )
         else if( ( incd.eq.incx ).and.( incd.gt.0 ) )then
            do 10, id = 1, 1 + ( n - 1 )*incd, incd
               x( id ) = d( id )*x( id )
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incd.gt.0 )then
               do 20, id = 1, 1 + ( n - 1 )*incd, incd
                  x( ix ) = d( id )*x( ix )
                  ix      = ix              + incx
   20          continue
            else
               id = 1 - ( n - 1 )*incd
               do 30, i = 1, n
                  x( ix ) = d( id )*x( ix )
                  id      = id              + incd
                  ix      = ix              + incx
   30          continue
            end if
         end if
      end if
C
C     end of ddscl (f06fcf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dddiv ( n, d, incd, x, incx )

      implicit           double precision (a-h,o-z)
      double precision   d(*), x(*)

*     dddiv  performs the diagonal scaling  x  =  x / d.

      integer            i, id, ix
      external           dscal
      intrinsic          abs
      parameter        ( one = 1.0d+0 )

      if (n .gt. 0) then
         if (incd .eq. 0  .and.  incx .ne. 0) then
            call dscal ( n, one/d(1), x, abs(incx) )
         else if (incd .eq. incx  .and.  incd .gt. 0) then
            do 10 id = 1, 1 + (n - 1)*incd, incd
               x(id) = x(id) / d(id)
   10       continue
         else
            if (incx .ge. 0) then
               ix = 1
            else
               ix = 1 - (n - 1)*incx
            end if
            if (incd .gt. 0) then
               do 20 id = 1, 1 + (n - 1)*incd, incd
                  x(ix) = x(ix) / d(id)
                  ix    = ix   + incx
   20          continue
            else
               id = 1 - (n - 1)*incd
               do 30  i = 1, n
                  x(ix) = x(ix) / d(id)
                  id    = id + incd
                  ix    = ix + incx
   30          continue
            end if
         end if
      end if

*     end of dddiv
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine icopy ( n, x, incx, y, incy )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     ENTRY      F06DFF( N, X, INCX, Y, INCY )

C     .. Scalar Arguments ..
      integer            incx, incy, n
C     .. Array Arguments ..
      integer            x( * ), y( * )
C     ..
C
C  F06DFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 10-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      integer            i, ix, iy
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
            do 10, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( iy )
   10       continue
         else
            if( incx.ge.0 )then
               ix = 1
            else
               ix = 1 - ( n - 1 )*incx
            end if
            if( incy.gt.0 )then
               do 20, iy = 1, 1 + ( n - 1 )*incy, incy
                  y( iy ) = x( ix )
                  ix      = ix      + incx
   20          continue
            else
               iy = 1 - ( n - 1 )*incy
               do 30, i = 1, n
                  y( iy ) = x( ix )
                  iy      = iy      + incy
                  ix      = ix      + incx
   30          continue
            end if
         end if
      end if
C
C     end of icopy (icopy)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine f06fjf( n, x, incx, scale, sumsq )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.

C     .. Scalar Arguments ..
      double precision   scale, sumsq
      integer            incx, n
C     .. Array Arguments ..                   
      double precision   x( * )
C     ..
C
C  f06fjf   returns the values scl and smsq such that
C
C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
C
C  where x( i ) = X( 1 + ( i - 1 )*INCX ). The value of sumsq is assumed
C  to be at least unity and the value of smsq will then satisfy
C
C     1.0 .le. smsq .le. ( sumsq + n ) .
C
C  scale is assumed to be non-negative and scl returns the value
C
C     scl = max( scale, abs( x( i ) ) ) .
C
C  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
C
C  The routine makes only one pass through the vector X.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 22-October-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   absxi
      integer            ix
C     .. Intrinsic Functions ..
      intrinsic          abs
C     ..
C     .. Executable Statements ..
      if( n.gt.0 )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  sumsq = 1     + sumsq*( scale/absxi )**2
                  scale = absxi
               else
                  sumsq = sumsq +       ( absxi/scale )**2
               end if
            end if
   10    continue
      end if
C
C     end of f06fjf
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dcond ( n, x, incx, xmax, xmin )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     ENTRY      F06FLF( N, X, INCX, XMAX, XMIN )
C     .. Scalar Arguments ..
      double precision   xmax, xmin
      integer            incx, n
C     .. Array Arguments ..
      double precision   x( * )
C     ..
C
C  dcond/f06flf returns the values xmax and xmin given by
C
C     xmax = max( abs( x( i ) ) ),   xmin = min( abs( x( i ) ) ).
C             i                              i
C
C  If n is less than unity then xmax and xmin are returned as zero.
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
C     .. Local Scalars ..
      integer            ix
C     .. Intrinsic Functions ..
      intrinsic          abs, max, min
C     ..
C     .. Executable Statements ..
      if( n.lt.1 )then
         xmax = zero
         xmin = zero
      else
         xmax = abs( x( 1 ) )
         xmin = xmax
         do 10 ix = 1 + incx, 1 + ( n - 1 )*incx, incx
            xmax = max( xmax, abs( x( ix ) ) )
            xmin = min( xmin, abs( x( ix ) ) )
   10    continue
      end if
C
C     end of dcond (f06flf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function idrank( n, x, incx, tol )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     INTEGER          F06KLF                    
C     ENTRY            F06KLF( N, X, INCX, TOL )
C     Modified by PEG 9/25/88.

C     .. Scalar Arguments ..
      double precision         tol
      integer                  incx, n
C     .. Array Arguments ..
      double precision         x( * )
C     ..
C
C  idrank/ F06KLF finds the first element of the n-vector x for which
C
C     abs( x( k ) ).le.tol*max( abs( x( 1 ) ), ..., abs( x( k - 1 ) ) )
C
C  and returns the value ( k - 1 ) in the function name F06KLF. If no
C  such k exists then F06KLF is returned as n.
C
C  If tol is supplied as less than zero then the value epsmch, where
C  epsmch is the relative machine precision, is used in place of tol.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 27-February-1986.
C     Sven Hammarling, Nag Central Office.
C                        
C     .. Parameters ..
      double precision         zero
      parameter              ( zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision         tl, xmax
      integer                  ix, k
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
C     .. Intrinsic Functions ..
      intrinsic                abs, max
C     ..
C     .. Executable Statements ..
      k = 0
      if( n.ge.1 )then
         ix = 1
         if( tol.lt.zero )then
            tl = wmach(3)
         else
            tl = tol
         end if
         xmax = abs( x( ix ) )
C
C+       while( k.lt.n )loop
   10    if   ( k.lt.n )then
            if( abs( x( ix ) ).le.tl*xmax )
     $         go to 20
            xmax = max( xmax, abs( x( ix ) ) )
            k    = k  + 1
            ix   = ix + incx
            go to 10
         end if
C+       end while
C
      end if
C
   20 idrank = k
C     f06klf = k
C
C     end of idrank (f06klf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06FQF( PIVOT, DIRECT, N, ALPHA, X, INCX, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, N
      CHARACTER*1        DIRECT, PIVOT
C     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * ), X( * )
C     ..
C
C  F06FQF generates the parameters of an orthogonal matrix P such that
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'F' or 'f'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'B' or 'b'
C
C        P*( alpha ) = ( beta ),
C          (   x   )   (   0  )
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'B' or 'b'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'F' or 'f'
C
C        P*(   x   ) = (   0  ),
C          ( alpha ) = ( beta )
C
C  where alpha is a scalar and x is an n element vector.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C        
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( 1, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( 1 )*P( 2 )*...*P( n )
C
C     where P( k ) is a plane rotation matrix for the ( k, n + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P( n )*P( n - 1 )*...*P( 1 )
C
C     where P( k ) is a plane rotation matrix for the ( k, k + 1 ) plane
C     designed to annihilate the kth element of x.
C
C  The routine returns the cosine, c( k ), and sine, s( k ) that define
C  the matrix P( k ), such that the two by two rotation part of P( k ),
C  R( k ), has the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  On entry, ALPHA must contain  the scalar alpha and on exit, ALPHA is
C  overwritten by beta. The cosines and sines are returned in the arrays
C  C and S and the vector x is overwritten by the tangents of the plane
C  rotations ( t( k ) = s( k )/c( k ) ).
C
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 19-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
            IX = 1 + ( N - 1 )*INCX
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
               DO 10, I = N, 2, -1
                  CALL F06BAF( X( IX - INCX ), X( IX ), C( I ), S( I ) )
                  IX = IX - INCX
   10          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( 1 ), S( 1 ) )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( alpha ) := (  c  s )*( alpha  )
C                 (   0   )    ( -s  c ) ( x( i ) )
C
C              which is equivalent to
C
C                 (   0   ) := ( c  -s )*( x( i ) )
C                 ( alpha )    ( s   c ) ( alpha  )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 20, I = N, 1, -1
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      - INCX
   20          CONTINUE
            END IF
         ELSE IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
            IX = 1
            IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
C
C              Here we choose c and s so that
C
C                 ( x( i + 1 ) ) := (  c  s )*( x( i + 1 ) )
C                 (    0       )    ( -s  c ) ( x( i )     )
C
C              which is equivalent to
C
C                 (    0       ) := ( c  -s )*( x( i )     )
C                 ( x( i + 1 ) )    ( s   c ) ( x( i + 1 ) )
C
C              and so we need to return  s( i ) = -s  in order to make
C              R( i ) look like
C
C                 R( i ) = (  c( i )  s( i ) ).
C                          ( -s( i )  c( i ) )
C
               DO 30, I = 1, N - 1
                  CALL F06BAF( X( IX + INCX ), X( IX ), C( I ), S( I ) )
                  S( I )  = -S( I )
                  X( IX ) = -X( IX )
                  IX      =  IX      + INCX
   30          CONTINUE
               CALL F06BAF( ALPHA, X( IX ), C( N ), S( N ) )
               S( N )  = -S( N )
               X( IX ) = -X( IX )
            ELSE IF( ( PIVOT.EQ.'F' ).OR.( PIVOT.EQ.'f' ) )THEN
               DO 40, I = 1, N
                  CALL F06BAF( ALPHA, X( IX ), C( I ), S( I ) )
                  IX = IX + INCX
   40          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06FQF. ( SSROTG )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dgrfg ( n, alpha, x, incx, tol, zeta )

C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     ENTRY      F06FRF( N, ALPHA, X, INCX, TOL, ZETA )
C     Modified by PEG 9/25/88.

C     .. Scalar Arguments ..
      double precision   alpha, tol, zeta
      integer            incx, n
C     .. Array Arguments ..
      double precision   x( * )
C     ..
C
C  dgrfg/f06frf generates a generalized Householder reflection such that
C
C     P*( alpha ) = ( beta ),   P'*P = I.
C       (   x   )   (   0  )
C
C  P is given in the form
C
C     P = I - ( zeta )*( zeta  z' ),
C             (   z  )
C
C  where z is an n element vector and zeta is a scalar that satisfies
C
C     1.0 .le. zeta .le. sqrt( 2.0 ).
C
C  zeta is returned in ZETA unless x is such that
C
C     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
C
C  where eps is the relative machine precision and tol is the user
C  supplied value TOL, in which case ZETA is returned as 0.0 and P can
C  be taken to be the unit matrix.
C
C  beta is overwritten on alpha and z is overwritten on x.
C  the routine may be called with  n = 0  and advantage is taken of the
C  case where  n = 1.
C
C
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 30-August-1984.
C     Sven Hammarling, Nag Central Office.
C     This version dated 28-September-1984.
C
C
C     .. Parameters ..
      double precision   one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
C     .. Local Scalars ..
      double precision   beta, eps, scale, ssq
      logical            first
      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/
C     .. External Subroutines ..
      external           f06fjf, dscal
C     .. Intrinsic Functions ..
      intrinsic          abs, max, sign, sqrt
C     .. Save statement ..
      save               eps, first
C     .. Data statements ..
      data               first/ .true. /
C     ..
C     .. Executable Statements ..
      if( n.lt.1 )then
         zeta = zero
      else if( ( n.eq.1 ).and.( x( 1 ).eq.zero ) )then
         zeta = zero
      else
C
         if( first )then
            first = .false.
            eps   =  wmach(3)
         end if
C
C        Treat case where P is a 2 by 2 matrix specially.
C
         if( n.eq.1 )then
C
C           Deal with cases where  ALPHA = zero  and
C           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
C
            if( alpha.eq.zero )then
               zeta   =  one
               alpha  =  abs ( x( 1 ) )
               x( 1 ) = -sign( one, x( 1 ) )
            else if( abs( x( 1 ) ).le.max( eps*abs( alpha ), tol ) )then
               zeta   =  zero
            else
               if( abs( alpha ).ge.abs( x( 1 ) ) )then
                  beta = abs( alpha ) *sqrt( 1 + ( x( 1 )/alpha )**2 )
               else
                  beta = abs( x( 1 ) )*sqrt( 1 + ( alpha/x( 1 ) )**2 )
               end if
               zeta = sqrt( ( abs( alpha ) + beta )/beta )
               if( alpha.ge.zero )
     $            beta = -beta
               x( 1 ) = -x( 1 )/( zeta*beta )
               alpha  = beta
            end if
         else
C
C           Now P is larger than 2 by 2.
C
            ssq   = one
            scale = zero
            call f06fjf( n, x, incx, scale, ssq )
C
C           Treat cases where  SCALE = zero,
C           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
C           ALPHA = zero  specially.
C           Note that  SCALE = max( abs( X( i ) ) ).
C
            if( ( scale.eq.zero ).or.
     $          ( scale.le.max( eps*abs( alpha ), tol ) ) )then
               zeta  = zero
            else if( alpha.eq.zero )then
               zeta  = one
               alpha = scale*sqrt( ssq )
               call dscal( n, -1/alpha, x, incx )
            else
               if( scale.lt.abs( alpha ) )then
                  beta = abs( alpha )*sqrt( 1 + ssq*( scale/alpha )**2 )
               else
                  beta = scale       *sqrt( ssq +   ( alpha/scale )**2 )
               end if
               zeta = sqrt( ( beta + abs( alpha ) )/beta )
               if( alpha.gt.zero )
     $            beta = -beta
               call dscal( n, -1/( zeta*beta ), x, incx )
               alpha = beta
            end if
         end if
      end if
C
C     end of dgrfg (f06frf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QFF( MATRIX, M, N, A, LDA, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      INTEGER            M, N, LDA, LDB
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
C     ..
C
C  F06QFF  copies  the  m by n  matrix  A  into  the  m by n  matrix  B.
C
C  If   MATRIX = 'G' or 'g'   then  A  and  B  are  regarded as  general
C                             matrices,
C  if   MATRIX = 'U' or 'u'   then  A  and  B  are  regarded  as   upper
C                             triangular,  and only  elements  for which
C                             i.le.j  are referenced,
C  if   MATRIX = 'L' or 'l'   then  A  and  B  are  regarded  as   lower
C                             triangular,  and only  elements  for which
C                             i.ge.j  are referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 40 J = 1, N
            DO 30 I = 1, MIN( M, J )
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 60 J = 1, MIN( M, N )
            DO 50 I = J, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QFF. ( SMCOPY )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      DOUBLE PRECISION FUNCTION F06QGF( NORM, MATRIX, M, N, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER                           LDA, M, N
      CHARACTER*1                       MATRIX, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION                  A( LDA, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QGF  returns the value of the one norm,  or the Frobenius norm, or
C  the element of  largest absolute value of a real matrix A.  A  may be
C  rectangular,  or square,  or triangular,  or symmetric.
C
C  Description
C  ===========
C
C  F06QGF returns the value
C
C     F06QGF = ( max( abs( a( i, j ) ) ) , NORM = 'M' or 'm'
C              (
C              ( norm1( A )  ,             NORM = '1', 'O' or 'o'
C              (
C              ( normF( A ) ,              NORM = 'F', 'f', 'E' or 'e'
C
C  where norm1 denotes the one norm of a matrix (maximum column sum) and
C  normF denotes the  Frobenius norm of a matrix  (square root of sum of
C  squares).  Note that  max( abs( a( i, j ) ) )  is not a  matrix norm.
C
C  The type of matrix for which  F06QGF is returned is determined by the
C  parameter MATRIX.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded as  a general matrix,
C  If   MATRIX = 'U' or 'u'   then  A  is regarded as  upper triangular,
C  If   MATRIX = 'L' or 'l'   then  A  is regarded as  lower triangular,
C  If   MATRIX = 'S' or 's'   then  A  is regarded as symmetric and only
C             or 'H' or 'h'   the  upper triangular part of the array  A
C                             is referenced,
C  If   MATRIX = 'Y' or 'y'   then  A  is regarded as symmetric and only
C             or 'E' or 'e'   the  lower triangular part of the array  A
C                             is referenced.
C
C  Parameters
C  ==========
C
C  NORM  -  CHARACTER*1.
C
C           On entry,  NORM specifies the value to be returned in F06QGF
C           as described above.
C
C           Unchanged on exit.
C
C  MATRIX - CHARACTER*1.
C
C           On entry,  MATRIX  specifies the type of matrix and,  in the
C           case of a  symmetric matrix,  the part of the array in which
C           the matrix is stored as described above.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry,  M  specifies the number of rows of the matrix  A.
C           M  must be at least  zero and when the  matrix is  symmetric
C           then  M must be equal to  N. When  M = 0  then F06QGF is set
C           to zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N  must be at least zero. When  N = 0  then F06QGF is set to
C           zero and an immediate return is effected.
C
C           Unchanged on exit.
C
C  A      - REAL array of DIMENSION ( LDA, n ).
C
C           Before entry,  A  must contain the  m by n  matrix for which
C           F06QGF is required.
C
C           If  MATRIX = 'U' or 'u' or 'S' or 's' or 'H' or 'h' then the
C           strictly lower triangular part of A is not referenced.
C
C           If  MATRIX = 'L' or 'l' or 'Y' or 'y' or 'E' or 'e' then the
C           strictly upper triangular part of A is not referenced.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in  the  calling  (sub)  program.  LDA  must be at least  M.
C
C           Unchanged on exit.
C
C  Further comments
C  ================
C
C  If A is part of a matrix B partitioned as
C
C     B = ( B1  B2 ) ,
C         ( B3  A  )
C
C  where  B1 is an l by k matrix  ( l.ge.0, k.ge.0 ),  then this routine
C  may be called with the parameter  A as  b( l + 1, k + 1 ) and  LDA as
C  the first dimension of  B  as declared in the calling  (sub) program.
C
C  This routine  can be inefficient on  paged machines when the one norm
C  is required, the matrix is symmetric and N is large.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION                  ONE, ZERO
      PARAMETER                         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION                  SCALE, SUM, VALUE
      INTEGER                           I, J
C     .. External Functions ..
      DOUBLE PRECISION                  F06BMF
      EXTERNAL                          F06BMF
C     .. External Subroutines ..
      EXTERNAL                          F06FJF
C     .. Intrinsic Functions ..
      INTRINSIC                         ABS, MAX, MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( M, N ).EQ.0 )THEN
         VALUE = ZERO
      ELSE IF( ( NORM.EQ.'M' ).OR.( NORM.EQ.'m' ) )THEN
C
C        Find  max( abs( a( i, j ) ) ).
C
         VALUE = ZERO
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 20 J = 1, N
               DO 10 I = 1, M
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ).OR.
     $            ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
            DO 40 J = 1, N
               DO 30 I = 1, MIN( M, J )
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ).OR.
     $            ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            DO 60 J = 1, MIN( M, N )
               DO 50 I = J, M
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE IF( ( NORM.EQ.'1' ).OR.( NORM.EQ.'O' ).OR.
     $         ( NORM.EQ.'o' ) )THEN
C
C        Find  norm1( A ).
C
         VALUE = ZERO
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 80 J = 1, N
               SUM = ZERO
               DO 70 I = 1, M
                  SUM = SUM + ABS( A( I, J ) )
   70          CONTINUE
               VALUE = MAX( VALUE, SUM )
   80       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
            DO 100 J = 1, N
               SUM = ZERO
               DO 90 I = 1, MIN( M, J )
                  SUM = SUM + ABS( A( I, J ) )
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
            DO 120 J = 1, MIN( M, N )
               SUM = ZERO
               DO 110 I = J, M
                  SUM = SUM + ABS( A( I, J ) )
  110          CONTINUE
               VALUE = MAX( VALUE, SUM )
  120       CONTINUE
         ELSE IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
            DO 150 J = 1, N
               SUM = ZERO
               DO 130 I = 1, J
                  SUM = SUM + ABS( A( I, J ) )
  130          CONTINUE
               DO 140 I = J + 1, N
                  SUM = SUM + ABS( A( J, I ) )
  140          CONTINUE
               VALUE = MAX( VALUE, SUM )
  150       CONTINUE
         ELSE IF( ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            DO 180 J = 1, N
               SUM = ZERO
               DO 160 I = 1, J - 1
                  SUM = SUM + ABS( A( J, I ) )
  160          CONTINUE
               DO 170 I = J, N
                  SUM = SUM + ABS( A( I, J ) )
  170          CONTINUE
               VALUE = MAX( VALUE, SUM )
  180       CONTINUE
         END IF
      ELSE IF( ( NORM.EQ.'F' ).OR.( NORM.EQ.'f' ).OR.( NORM.EQ.'E' ).OR.
     $         ( NORM.EQ.'e' ) )THEN
C
C        Find  normF( A ).
C
         SCALE = ZERO
         SUM = ONE
         IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
            DO 190 J = 1, N
               CALL F06FJF( M, A( 1, J ), 1, SCALE, SUM )
  190       CONTINUE
         ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
            DO 200 J = 1, N
               CALL F06FJF( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
  200       CONTINUE
         ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
            DO 210 J = 1, MIN( M, N )
               CALL F06FJF( M - J + 1, A( J, J ), 1, SCALE, SUM )
  210       CONTINUE
         ELSE IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $            ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ).OR.
     $            ( MATRIX.EQ.'Y' ).OR.( MATRIX.EQ.'y' ).OR.
     $            ( MATRIX.EQ.'E' ).OR.( MATRIX.EQ.'e' ) )THEN
            IF( ( MATRIX.EQ.'S' ).OR.( MATRIX.EQ.'s' ).OR.
     $          ( MATRIX.EQ.'H' ).OR.( MATRIX.EQ.'h' ) )THEN
               DO 220 J = 2, N
                  CALL F06FJF( J - 1, A( 1, J ), 1, SCALE, SUM )
  220          CONTINUE
            ELSE
               DO 230 J = 1, N - 1
                  CALL F06FJF( N - J, A( J + 1, J ), 1, SCALE, SUM )
  230          CONTINUE
            END IF
            SUM = 2*SUM
            CALL F06FJF( N, A( 1, 1 ), LDA + 1, SCALE, SUM )
         END IF
         VALUE = F06BMF( SCALE, SUM )
      END IF
C
      F06QGF = VALUE
      RETURN
C
C     End of F06QGF. ( SMNRM )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QHF( MATRIX, M, N, CONST, DIAG, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      CHARACTER*1        MATRIX
      DOUBLE PRECISION   CONST, DIAG
      INTEGER            LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
C     ..
C
C  F06QHF forms the m by n matrix A given by
C
C     a( i, j ) = (  diag  i.eq.j,
C                 (
C                 ( const  i.ne.j.
C
C  If   MATRIX = 'G' or 'g'   then  A  is regarded  as a general matrix,
C  if   MATRIX = 'U' or 'u'   then  A  is regarded  as upper triangular,
C                             and only  elements  for which  i.le.j  are
C                             referenced,
C  if   MATRIX = 'L' or 'l'   then  A  is regarded  as lower triangular,
C                             and only  elements  for which  i.ge.j  are
C                             referenced.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 21-November-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MATRIX.EQ.'G' ).OR.( MATRIX.EQ.'g' ) )THEN
         DO 20 J = 1, N
            DO 10 I = 1, M
               A( I, J ) = CONST
   10       CONTINUE
   20    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 30 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   30       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'U' ).OR.( MATRIX.EQ.'u' ) )THEN
         DO 50 J = 1, N
            DO 40 I = 1, MIN( M, J )
               A( I, J ) = CONST
   40       CONTINUE
   50    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 60 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   60       CONTINUE
         END IF
      ELSE IF( ( MATRIX.EQ.'L' ).OR.( MATRIX.EQ.'l' ) )THEN
         DO 80 J = 1, MIN( M, N )
            DO 70 I = J, M
               A( I, J ) = CONST
   70       CONTINUE
   80    CONTINUE
         IF( CONST.NE.DIAG )THEN
            DO 90 I = 1, MIN( M, N )
               A( I, I ) = DIAG
   90       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QHF. ( SMLOAD )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QKF( SIDE, TRANS, N, PERM, K, B, LDB )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K, LDB, N
      CHARACTER*1        SIDE, TRANS
C     .. Array Arguments ..
      DOUBLE PRECISION   PERM( * ), B( LDB, * )
C     ..
C
C  Purpose
C  =======
C
C  F06QKF performs one of the transformations
C
C     B := P'*B   or   B := P*B,   where B is an m by k matrix,
C
C  or
C
C     B := B*P'   or   B := B*P,   where B is a k by m matrix,
C
C  P being an m by m permutation matrix of the form
C
C     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
C
C  where  P( i, index( i ) ) is the permutation matrix that interchanges
C  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
C  with rows and columns  i and  index( i )  interchanged. Of course, if
C  index( i ) = i  then  P( i, index( i ) ) = I.
C
C  This  routine is intended  for use in conjunction with  Nag auxiliary
C  routines  that  perform  interchange  operations,  such  as  sorting.
C
C  Parameters
C  ==========
C
C  SIDE   - CHARACTER*1.
C  TRANS
C           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
C           TRANS  ( Transpose, or No transpose )  specify the operation
C           to be performed as follows.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := P'*B.
C
C           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := P*B.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
C
C              Perform the operation   B := B*P'.
C
C           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
C
C              Perform the operation   B := B*P.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the value of n.  N must be at least
C           zero.  When  N = 0  then an  immediate  return  is effected.
C
C           Unchanged on exit.
C
C  PERM   - REAL             array of DIMENSION at least ( n ).
C
C           Before  entry,  PERM  must  contain  the  n indices  for the
C           permutation matrices. index( i ) must satisfy
C
C              1 .le. index( i ) .le. m.
C
C           It is usual for index( i ) to be at least i, but this is not
C           necessary for this routine. It is assumed that the statement
C           INDEX = PERM( I )  returns the correct integer in  INDEX, so
C           that,  if necessary,  PERM( I )  should contain a real value
C           slightly larger than  INDEX.
C
C           Unchanged on exit.
C
C  K      - INTEGER.
C
C           On entry with  SIDE = 'L' or 'l',  K must specify the number
C           of columns of B and on entry with  SIDE = 'R' or 'r', K must
C           specify the number of rows of  B.  K must be at least  zero.
C           When  K = 0  then an immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - REAL  array  of  DIMENSION ( LDB, ncolb ),  where  ncolb = k
C           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
C           'r'.
C
C           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
C           part  of  the  array   B  must  contain  the  matrix  to  be
C           transformed  and before  entry with  SIDE = 'R' or 'r',  the
C           leading  K by m part of the array  B must contain the matrix
C           to  be  transformed.  On exit,   B  is  overwritten  by  the
C           transformed matrix.
C
C  LDB    - INTEGER.
C
C           On entry,  LDB  must specify  the  leading dimension  of the
C           array  B  as declared  in the  calling  (sub) program.  When
C           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
C           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
C           Unchanged on exit.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 11-August-1987.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      LOGICAL            LEFT, NULL, RIGHT, TRNSP
      INTEGER            I, J, L
      DOUBLE PRECISION   TEMP
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( MIN( N, K ).EQ.0 )
     $   RETURN
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      NULL = ( TRANS.EQ.'N' ).OR.( TRANS.EQ.'n' )
      TRNSP = ( TRANS.EQ.'T' ).OR.( TRANS.EQ.'t' )
      IF( LEFT )THEN
         IF( TRNSP )THEN
            DO 20 I = 1, N
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 10 J = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( L, J )
                     B( L, J ) = TEMP
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE IF( NULL )THEN
            DO 40 I = N, 1, -1
               L = PERM( I )
               IF( L.NE.I )THEN
                  DO 30 J = 1, K
                     TEMP = B( L, J )
                     B( L, J ) = B( I, J )
                     B( I, J ) = TEMP
   30             CONTINUE
               END IF
   40       CONTINUE
         END IF
      ELSE IF( RIGHT )THEN
         IF( TRNSP )THEN
            DO 60 J = N, 1, -1
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 50 I = 1, K
                     TEMP = B( I, J )
                     B( I, J ) = B( I, L )
                     B( I, L ) = TEMP
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE IF( NULL )THEN
            DO 80 J = 1, N
               L = PERM( J )
               IF( L.NE.J )THEN
                  DO 70 I = 1, K
                     TEMP = B( I, L )
                     B( I, L ) = B( I, J )
                     B( I, J ) = TEMP
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06QKF. ( SGEAPR )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QNF( SIDE, N, K1, K2, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * )
C     ..
C
C  F06QNF applies a  sequence  of  pairwise interchanges to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The interchanges are
C  applied in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
C  The  two by two
C  interchange part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( 0  1 ).
C              ( 1  0 )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the permutations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ    = ZERO
               S( J ) = A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP          = A( I, J )
               A( I + 1, J ) = TEMP
               AIJ           = AIJ
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane interchanges to  columns  k1  up to
C        ( k2 - 1 ) and  form   the   additional  sub-diagonal
C        elements,   storing  h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            DO 30 I = 1, J
               TEMP = A( I, J + 1 )
               A( I, J + 1 ) = A( I, J )
               A( I, J )     = TEMP
   30       CONTINUE
            S( J )            = A( J + 1, J + 1 )
            A( J + 1, J + 1 ) = ZERO
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QNF. ( SUTSRH )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QRF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QRF restores an upper Hessenberg matrix H to upper triangular form
C  by  applying a sequence of  plane rotations  from either the left, or
C  the right.  The matrix  H  is assumed to have  non-zero  sub-diagonal
C  elements  in  positions  h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only  and  h( k + 1, k )  must  be  supplied  in  s( k ).
C
C  H is restored to the upper triangular matrix R either as
C
C     R = P*H,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  or as
C
C     R = H*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  in both cases  P( k )  being a  plane rotation  for the  ( k, k + 1 )
C  plane.  The cosine and sine that define P( k ) are returned in c( k )
C  and  s( k )  respectively.  The two by two  rotation part of  P( k ),
C  Q( k ), is of the form
C
C     Q( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The upper triangular part of the matrix  H  must be supplied in the n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, SUBH, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )     ).
C           (     0     )    ( -s  c ) ( h( j + 1, j ) )
C
C        Apply the rotations in columns k1 up to n.
C
         DO 20 J = K1, N
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J, K2 ) - 1
               TEMP = A( I + 1, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               AIJ = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            IF( J.LT.K2 )THEN
C
C              Set up the rotation.
C
               SUBH = S( J )
               CALL F06BAF( AIJ, SUBH, C( J ), S( J ) )
               A( J, J ) = AIJ
            ELSE
               A( K2, J ) = AIJ
            END IF
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore   H  to  upper  triangular  form  by  annihilating  the
C        sub-diagonal elements of H.  The jth rotation is chosen so that
C
C           ( h( j + 1, j + 1 ) ) := (  c  s )*( h( j + 1, j + 1 ) ),
C           (         0         )    ( -s  c ) ( h( j + 1, j )     )
C
C        which can be expressed as
C
C           ( 0  h( j + 1, j + 1 ) ) :=
C
C               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c  s ).
C                                                    ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 40 J = K2 - 1, K1, -1
            SUBH = S( J )
            CALL F06BAF( A( J + 1, J + 1 ), SUBH, CTEMP, STEMP )
            STEMP = -STEMP
            S( J ) = STEMP
            C( J ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30 I = J, 1, -1
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QRF. ( SUHQR )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QSF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QSF restores an upper spiked matrix  H to upper triangular form by
C  applying a sequence of plane rotations, in planes  k1 up to k2,  from
C  either the left, or the right.
C
C  The matrix  H is assumed to have non-zero elements only in the spiked
C  positions, h( k2, k ) for a row spike and h( k + 1, k1 ) for a column
C  spike, k = k1, k1 + 1, ..., k2 - 1, and these must be supplied in the
C  elements s( k ).
C
C  When  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     H  is  assumed  to have a  row spike  and is restored to the upper
C     triangular matrix  R as
C
C        R = P*H,
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C     P( k )  being a  plane rotation  matrix for the  ( k, k2 )  plane.
C
C  When  SIDE = 'R' or 'r'  ( Right-hand side )
C
C     H  is assumed to have a  column spike and is restored to the upper
C     triangular matrix R as
C
C        R = H*P',
C
C     where P is an orthogonal matrix of the form
C
C        P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C     P( k ) being a plane rotation matrix for the  ( k1, k + 1 ) plane.
C
C  The  two by two  rotation  part of  P( k ),  Q( k ),  is of  the form
C
C     Q( k ) = (  c( k )  s( k ) )
C              ( -s( k )  c( k ) )
C
C  and  c( k ) and s( k ) are returned in the kth elements of the arrays
C  C and S respectively.
C
C  The upper triangular part of the matrix  H must be supplied in the  n
C  by n  leading upper triangular part of  A, and this is overwritten by
C  the upper triangular matrix R.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Restore H to upper triangular form by annihilating the elements
C        in  the  spike  of  H.  The  jth rotation  is  chosen  so  that
C
C        ( h( j, j ) ) := (  c  s )*( h( j , j ) ).
C        (     0     )    ( -s  c ) ( h( k2, j ) )
C
C        Apply the rotations in columns k1 up to ( k2 - 1 ).
C
         DO 20 J = K1, K2 - 1
            SPIKE = S( J )
            DO 10 I = K1, J - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   10       CONTINUE
C
C           Set up the rotation.
C
            CALL F06BAF( A( J, J ), SPIKE, C( J ), S( J ) )
   20    CONTINUE
C
C        Apply the rotations to columns k2 up to n.
C
         DO 40 J = K2, N
            TEMP = A( K2, J )
            DO 30 I = K1, K2 - 1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   30       CONTINUE
            A( K2, J ) = TEMP
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Restore H to upper triangular form by annihilating the spike of
C        H. The jth rotation is chosen so that
C
C           ( h( j, j ) ) := (  c  s )*( h( j, j )  ),
C           (     0     )    ( -s  c ) ( h( j, k1 ) )
C
C        which can be expressed as
C
C           ( 0  h( j, j ) ) := ( h( j, k1 )  h( j, j ) )*(  c  s ).
C                                                         ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           Q( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 70 J = K2, K1 + 1, -1
            CALL F06BAF( A( J, J ), S( J - 1 ), CTEMP, STEMP )
            STEMP = -STEMP
            S( J - 1 ) = STEMP
            C( J - 1 ) = CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = J - 1, K1 + 1, -1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   50          CONTINUE
               DO 60 I = K1, 1, -1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   60          CONTINUE
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QSF. ( SUSQR )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QTF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QTF performs the transformation
C
C     R := P*U*Q'  when  SIDE = 'L' or 'l'  (  Left-hand side )
C
C     R := Q*U*P'  when  SIDE = 'R' or 'r'  ( Right-hand side ),
C
C  where  U and R  are  n by n  upper  triangular  matrices,   P  is  an
C  orthogonal matrix,  consisting of a given sequence of plane rotations
C  to be  applied  in  planes  k1 to k2,  and  Q  is  a  unitary  matrix
C  consisting of a sequence of plane rotations, applied in planes  k1 to
C  k2,  chosen to make  R  upper triangular.
C
C  When  SIDE = 'L' or 'l'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k2 - 1 )*...*Q( k1 + 1 )*Q( k1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  When  SIDE = 'R' or 'r'  then  P  is  given  as a  sequence of  plane
C  rotation matrices
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  where  P( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C  In this case the matrix Q is given as
C
C     Q = Q( k1 )*Q( k1 + 1 )*...*Q( k2 - 1 ),
C
C  where  Q( k ) is a plane rotation matrix for the  ( k, k + 1 ) plane.
C
C  The  upper  triangular  matrix  U  must  be  supplied  in the  n by n
C  leading upper triangular part of  A,  and this  is overwritten by the
C  upper triangular matrix  R.  The cosine  and  sine  that  define  the
C  plane rotation matrix  P( k )  must be supplied in  c( k ) and s( k )
C  respectively,  and  the two by two rotation part of  P( k ),  T( k ),
C  is assumed to be of the form
C
C     T( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The cosine  and  sine that define  Q( k )  are overwritten on  c( k )
C  and  s( k )  respectively and the two by two rotation part of  Q( k )
C  will have the form of  T( k )  above.
C
C  If  n or k1  are less  than  unity, or  k1  is not  less than  k2, or
C  k2  is greater than  n  then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 26-November-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, FILL, STEMP, TEMP
      INTEGER            I, I1, J
C     .. External Subroutines ..
      EXTERNAL           F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the left-hand transformations,  column by column,  to the
C        triangular part of  U,  but not to  anywhere  that would  cause
C        fill.
C
         DO 20 J = K1 + 1, N
C
C           Apply  P( k1 ) ... P( j - 1 )  to column j.
C
            AIJ = A( K1, J )
            DO 10 I = K1, MIN( J - 1, K2 - 1 )
               A( I, J ) = S( I )*A( I + 1, J ) + C( I )*AIJ
               AIJ = C( I )*A( I + 1, J ) - S( I )*AIJ
   10       CONTINUE
            A( I, J ) = AIJ
   20    CONTINUE
C
C           Now apply each  left-hand tranformation  to form the fill-in
C           elements and apply a  right-hand transformation to eliminate
C           the fill-in element.
C
         DO 40 J = K1, K2 - 1
C
C           Apply  P( j )  to the jth diagonal element  and the  fill-in
C           position.
C
            FILL = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
C
C           Now  set up  the rotation  Q( j )  to eliminate the  fill-in
C           element,  and  apply  Q( j )  to  the  jth  and  ( j + 1 )th
C           columns.
C
            CALL F06BAF( A( J + 1, J + 1 ), FILL, CTEMP, STEMP )
            C( J ) = CTEMP
            S( J ) = -STEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               STEMP = -STEMP
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        We intermingle the  left and right hand transformations so that
C        at the kth step we form
C
C           A := Q( k )*A*P( k )'.
C
C        First  apply  the  transformations  in  columns  k2 back to k1.
C
         DO 60 J = K2 - 1, K1, -1
C
C           First apply  P( j ).
C
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               CTEMP = C( J )
               STEMP = S( J )
               DO 50 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   50          CONTINUE
C
C              Next form the fill-in element  a( j + 1, j )  by applying
C              P( j ).
C
               FILL = S( J )*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = C( J )*A( J + 1, J + 1 )
C
C              Now set up the rotation  Q( j )  to eliminate the fill-in
C              element.
C
               CALL F06BAF( A( J, J ), FILL, C( J ), S( J ) )
            END IF
   60    CONTINUE
C
C        Finally  apply  Q( k2 - 1 ) ... Q( k1 )  to columns  n  back to
C        ( k1 + 1 ).
C
         DO 80 J = N, K1 + 1, -1
            I1 = MIN( K2, J )
            AIJ = A( I1, J )
            DO 70 I = I1 - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   70       CONTINUE
            A( K1, J ) = AIJ
   80    CONTINUE
      END IF
      RETURN
C
C     End of F06QTF. ( SUTSQR )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QVF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QVF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The rotations are applied
C  in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k1.
C
         DO 20 J = N, K1, -1
            IF( J.GE.K2 )THEN
               AIJ = A( K2, J )
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = C( J )*A( J, J )
               S( J ) = -S( J )*A( J, J )
            END IF
            DO 10 I = MIN( K2, J ) - 1, K1, -1
               TEMP = A( I, J )
               A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
               AIJ = S( I )*AIJ + C( I )*TEMP
   10       CONTINUE
            A( K1, J ) = AIJ
   20    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
C        and  form   the   additional  sub-diagonal  elements,   storing
C        h( j + 1, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP = S( J )
               CTEMP = C( J )
               DO 30 I = 1, J
                  TEMP = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
               S( J ) = STEMP*A( J + 1, J + 1 )
               A( J + 1, J + 1 ) = CTEMP*A( J + 1, J + 1 )
            END IF
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QVF. ( SUTSRH )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QWF( SIDE, N, K1, K2, C, S, A, LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, N
      CHARACTER*1        SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QWF applies a  given sequence  of  plane rotations  to either  the
C  left,  or the right,  of the  n by n  upper triangular  matrix  U  to
C  transform  U  to an upper spiked matrix. The rotations are applied in
C  planes k1 up to k2.
C
C  The upper spiked matrix, H, is formed as
C
C     H = P*U,   when   SIDE = 'L' or 'l',  ( Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation matrix for the ( k, k2 ) plane, and is
C  formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k )  being a  plane rotation matrix for the  ( k1, k + 1 )  plane.
C
C  The cosine and sine that define  P( k ), k = k1, k1 + 1, ..., k2 - 1,
C  must be  supplied  in  c( k ) and s( k ) respectively. The two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of H.
C
C  When  SIDE = 'L' or 'l'  then a  row spike  is  generated  in  H  and
C  when  SIDE = 'R' or 'r'  then a  column spike is generated. For a row
C  spike the elements  h( k2, k )  and for a  column spike  the elements
C  h( k + 1, k1 ), k = k1, k1 + 1, ..., k2 - 1, are returned in  s( k ).
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, SPIKE, STEMP, TEMP
      INTEGER            I, J
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $   ( K2.GT.N ) )RETURN
      IF( ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' ) )THEN
C
C        Apply the plane rotations to columns n back to k2.
C
         DO 20 J = N, K2, -1
            TEMP = A( K2, J )
            DO 10 I = K2 - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*TEMP + C( I )*AIJ
               TEMP = C( I )*TEMP - S( I )*AIJ
   10       CONTINUE
            A( K2, J ) = TEMP
   20    CONTINUE
C
C        Form  the spike  and apply the rotations in columns  ( k2 - 1 )
C        back to k1.
C
         DO 40 J = K2 - 1, K1, -1
            SPIKE = -S( J )*A( J, J )
            A( J, J ) = C( J )*A( J, J )
            DO 30 I = J - 1, K1, -1
               AIJ = A( I, J )
               A( I, J ) = S( I )*SPIKE + C( I )*AIJ
               SPIKE = C( I )*SPIKE - S( I )*AIJ
   30       CONTINUE
            S( J ) = SPIKE
   40    CONTINUE
      ELSE IF( ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' ) )THEN
C
C        Apply the  plane rotations to columns  ( k1 + 1 ) up to k2  and
C        form the spike.
C
         DO 70 J = K1 + 1, K2
            CTEMP = C( J - 1 )
            STEMP = S( J - 1 )
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 50 I = 1, K1
                  TEMP = A( I, K1 )
                  A( I, K1 ) = STEMP*A( I, J ) + CTEMP*TEMP
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*TEMP
   50          CONTINUE
               DO 60 I = K1 + 1, J - 1
                  SPIKE = S( I - 1 )
                  S( I - 1 ) = STEMP*A( I, J ) + CTEMP*SPIKE
                  A( I, J ) = CTEMP*A( I, J ) - STEMP*SPIKE
   60          CONTINUE
               S( J - 1 ) = STEMP*A( J, J )
               A( J, J ) = CTEMP*A( J, J )
            END IF
   70    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QWF. ( SUTSRS )
C
      END

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE F06QXF( SIDE, PIVOT, DIRECT, M, N, K1, K2, C, S, A,
     $                   LDA )
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER            K1, K2, LDA, M, N
      CHARACTER*1        DIRECT, PIVOT, SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
C     ..
C
C  F06QXF  performs the transformation
C
C     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
C
C     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
C
C  where A is an m by n matrix and P is an orthogonal matrix, consisting
C  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
C  determined by the parameters PIVOT and DIRECT as follows:
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'T' or 't'  ( Top pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
C     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
C
C        P is given as a sequence of plane rotation matrices
C
C           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
C        plane.
C
C  c( k ) and s( k )  must contain the  cosine and sine  that define the
C  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
C  P( k ), R( k ), is assumed to be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
C  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
C  and  k2  is greater than  n,  then an  immediate return  is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 20-November-1986.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   AIJ, CTEMP, STEMP, TEMP
      INTEGER            I, J
      LOGICAL            LEFT, RIGHT
C     .. Intrinsic Functions ..
      INTRINSIC          MIN
C     ..
C     .. Executable Statements ..
      LEFT = ( SIDE.EQ.'L' ).OR.( SIDE.EQ.'l' )
      RIGHT = ( SIDE.EQ.'R' ).OR.( SIDE.EQ.'r' )
      IF( ( MIN( M, N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.
     $    (  LEFT .AND. K2.GT.M  ).OR.
     $    (  RIGHT .AND. K2.GT.N  ) )RETURN
      IF( LEFT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 20 J = 1, N
                  AIJ = A( K1, J )
                  DO 10 I = K1, K2 - 1
                     TEMP = A( I + 1, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     AIJ = C( I )*TEMP - S( I )*AIJ
   10             CONTINUE
                  A( K2, J ) = AIJ
   20          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 40 J = 1, N
                  AIJ = A( K2, J )
                  DO 30 I = K2 - 1, K1, -1
                     TEMP = A( I, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     AIJ = S( I )*AIJ + C( I )*TEMP
   30             CONTINUE
                  A( K1, J ) = AIJ
   40          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 60 J = 1, N
                  TEMP = A( K1, J )
                  DO 50 I = K1, K2 - 1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   50             CONTINUE
                  A( K1, J ) = TEMP
   60          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 80 J = 1, N
                  TEMP = A( K1, J )
                  DO 70 I = K2 - 1, K1, -1
                     AIJ = A( I + 1, J )
                     A( I + 1, J ) = C( I )*AIJ - S( I )*TEMP
                     TEMP = S( I )*AIJ + C( I )*TEMP
   70             CONTINUE
                  A( K1, J ) = TEMP
   80          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 100 J = 1, N
                  TEMP = A( K2, J )
                  DO 90 I = K1, K2 - 1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
   90             CONTINUE
                  A( K2, J ) = TEMP
  100          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 120 J = 1, N
                  TEMP = A( K2, J )
                  DO 110 I = K2 - 1, K1, -1
                     AIJ = A( I, J )
                     A( I, J ) = S( I )*TEMP + C( I )*AIJ
                     TEMP = C( I )*TEMP - S( I )*AIJ
  110             CONTINUE
                  A( K2, J ) = TEMP
  120          CONTINUE
            END IF
         END IF
      ELSE IF( RIGHT )THEN
         IF( ( PIVOT.EQ.'V' ).OR.( PIVOT.EQ.'v' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 140 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 130 I = 1, M
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 160 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 150 I = M, 1, -1
                        TEMP = A( I, J + 1 )
                        A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'T' ).OR.( PIVOT.EQ.'t' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 180 J = K1 + 1, K2
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 200 J = K2, K1 + 1, -1
                  CTEMP = C( J - 1 )
                  STEMP = S( J - 1 )
                  IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
                     DO 190 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, K1 )
                        A( I, K1 ) = STEMP*TEMP + CTEMP*A( I, K1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( ( PIVOT.EQ.'B' ).OR.( PIVOT.EQ.'b' ) )THEN
            IF( ( DIRECT.EQ.'F' ).OR.( DIRECT.EQ.'f' ) )THEN
               DO 220 J = K1, K2 - 1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( ( DIRECT.EQ.'B' ).OR.( DIRECT.EQ.'b' ) )THEN
               DO 240 J = K2 - 1, K1, -1
                  IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
                     CTEMP = C( J )
                     STEMP = S( J )
                     DO 230 I = M, 1, -1
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, K2 ) + CTEMP*TEMP
                        A( I, K2 ) = CTEMP*A( I, K2 ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06QXF. ( SGESRC )
C
      END
** END OF F06QXFTEXT
      SUBROUTINE F06QZF( HESS, N, K1, K2, C, S, A, LDA )
*     .. Scalar Arguments ..
      CHARACTER*1        HESS
      INTEGER            K1, K2, LDA, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  F06QZF  either applies a  given sequence  of  plane rotations  to the
*  right of the n by n reverse lower triangular matrix T, to transform T
*  to a  reverse lower Hessenberg matrix  H, or restores a reverse lower
*  Hessenberg matrix H to reverse lower triangular form T, by applying a
*  sequence of plane rotations from the right.
*
*  The rotations are applied  in planes k1 up to k2.
*
*  When   HESS = 'C' or 'c',   ( Create ),  then   the   reverse   lower
*  Hessenberg matrix, H, is formed as
*
*     H = T*P',
*
*  where P is an orthogonal matrix of the form
*
*     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
*
*  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
*  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
*  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
*  rotation part of P( k ), R( k ), is assumed to have the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  The matrix  T must be supplied in the n by n reverse lower triangular
*  part  of the array  A,  and this is overwritten by the  reverse lower
*  triangular part of  H.
*
*  The super-diagonal elements of  H, h( n - k, k ), are returned in the
*  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
*
*  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
*  greater than n then an immediate return is effected.
*
*  When   HESS = 'R' or 'r',   ( Remove ),  then   the   reverse   lower
*  Hessenberg matrix  H  is  assumed  to  have  non-zero  super-diagonal
*  elements  in  positions  h( n - k, k ),  k = k1, k1 + 1, ..., k2 - 1,
*  only and  h( n - k, k ) must be supplied in  s( k ). H is restored to
*  the reverse lower triangular matrix T as
*
*     T = H*P',
*
*  where P is an orthogonal matrix of the form
*
*     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
*
*  P( k ) being a plane rotation for the  ( k, k + 1 ) plane. The cosine
*  and  sine  that  define  P( k )  are  returned  in  c( k ) and s( k )
*  respectively.  The  two by two  rotation part of  P( k ),  R( k ), is
*  of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  The reverse lower triangular part of the matrix H must be supplied in
*  the  n by n  reverse  lower  triangular  part  of  A,   and  this  is
*  overwritten by the reverse triangular matrix T.
*
*  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
*  greater than n then an immediate return is effected.
*
*  When   n = 7, k1 = 2 and k2 = 5   then  T  and  H  are  of  the  form
*
*     T = ( 0  0  0  0  0  0  X ),   H = ( 0  0  0  0  0  0  X ).
*         ( 0  0  0  0  0  X  X )        ( 0  0  0  0  X  X  X )
*         ( 0  0  0  0  X  X  X )        ( 0  0  0  X  X  X  X )
*         ( 0  0  0  X  X  X  X )        ( 0  0  X  X  X  X  X )
*         ( 0  0  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
*         ( 0  X  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
*         ( X  X  X  X  X  X  X )        ( X  X  X  X  X  X  X )
*
*
*  This routine  is  principally intended  for use  with the  non-linear
*  optimization routines such as E04UCF, in order to help vectorization.
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 10-May-1988.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     .. External Subroutines ..
      EXTERNAL           F06BAF
*     .. Local Scalars ..
      DOUBLE PRECISION   CTEMP, STEMP, SUPH, TEMP
      INTEGER            I, J
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
      IF( ( MIN( N, K1 ).LT.1 ).OR.( K2.LE.K1 ).OR.( K2.GT.N ) )
     $   RETURN
      IF( ( HESS.EQ.'C' ).OR.( HESS.EQ.'c' ) )THEN
*
*        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
*        and  form   the  additional  super-diagonal  elements,  storing
*        h( n - j, j ) in s( j ).
*
         DO 20, J = K1, K2 - 1
            IF( ( C( J ).NE.ONE ).OR.( S( J ).NE.ZERO ) )THEN
               STEMP             = S( J )
               CTEMP             = C( J )
               S( J )            = STEMP*A( N - J, J + 1 )
               A( N - J, J + 1 ) = CTEMP*A( N - J, J + 1 )
               DO 10, I = N - J + 1, N
                  TEMP          = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J )     = STEMP*TEMP + CTEMP*A( I, J )
   10          CONTINUE
            END IF
   20    CONTINUE
      ELSE IF( ( HESS.EQ.'R' ).OR.( HESS.EQ.'r' ) )THEN
*
*        Restore  H to reverse lower triangular form by annihilating the
*        super-diagonal elements of  H.  The  jth rotation  is chosen so
*        that
*
*          ( h( n - j, n - j ) ) := (  c  s )*( h( n - j, n - j     ) ),
*          (         0         )    ( -s  c ) ( h( n - j, n - j - 1 ) )
*
*        which can be expressed as
*
*           ( 0  h( n - j, n - j ) ) :=
*
*               ( h( n - j, n - j - 1 )  h( n - j, n - j ) )*(  c  s ).
*                                                            ( -s  c )
*
*        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
*        rotation matrix look like
*
*           R( j ) = (  c( j )  s( j ) ).
*                    ( -s( j )  c( j ) )
*
         DO 40, J = K2 - 1, K1, -1
            SUPH   =  S( J )
            CALL F06BAF( A( N - J, J + 1 ), SUPH, CTEMP, STEMP )
            STEMP  = -STEMP
            S( J ) =  STEMP
            C( J ) =  CTEMP
            IF( ( CTEMP.NE.ONE ).OR.( STEMP.NE.ZERO ) )THEN
               DO 30, I = N - J + 1, N
                  TEMP          = A( I, J + 1 )
                  A( I, J + 1 ) = CTEMP*TEMP - STEMP*A( I, J )
                  A( I, J )     = STEMP*TEMP + CTEMP*A( I, J )
   30          CONTINUE
            END IF
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of F06QZF.
*
      END
** END OF F06QZFTEXT

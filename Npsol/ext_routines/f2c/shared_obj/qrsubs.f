*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  qrsubs.f
*
*     dgeqr    dgeqrp   dgeap    dgeapq
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dgeqr ( m, n, A, ldA, zeta, inform )
      integer            m, n, ldA, inform
      double precision   A( ldA, * ), zeta( * )
C
C  1. Purpose
C     =======
C
C  DGEQR  reduces the  m by n, m.ge.n, matrix A to upper triangular form
C  by means of orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )   when   m.gt.n,
C           ( 0 )
C
C     A = Q*R       when   m = n,
C
C  where  Q  is an  m by m  orthogonal matrix and  R  is an n by n upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used to introduce zeros  into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of zeta and in the
C  kth column of A, such that zeta( k ) is in zeta( k ) and the elements
C  of z( k ) are in A( k + 1, k ), ..., A( m, k ). The elements of R are
C  returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( n, m - 1 ).
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at least  n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be  at  least zero. When  N = 0  then an immediate return is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  N by N upper triangular part of A will contain
C           the  upper  triangular  matrix  R  and the  M by N  strictly
C           lower triangular part of  A  will  contain  details  of  the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  zeta   - 'real' array of DIMENSION at least ( n ).
C
C           On  exit, zeta( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then   zeta( k ) = 0.0
C           otherwise  zeta( k )  contains  zeta( k ) as described above
C           and is always in the range ( 1.0, sqrt( 2.0 ) ).
C
C  INFORM - INTEGER.
C
C           On successful  exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        M   .lt. N
C        N   .lt. 0
C        LDA .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear  algebra routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, zeta,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, zeta,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be a  k  element array  that  is used  as
C  workspace. If  B  is a one-dimensional array (single column) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Operations involving the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRMV  and DTRSV . Note that no test for near
C  singularity of R is incorporated in this routine or in routine  DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular  then  the  auxiliary linear algebra  routine  DUTSV
C  can  be used to  determine  the  singular value decomposition  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      external           dgemv , dger  , dgrfg
      intrinsic          min
      integer            k     , la
      double precision   temp
      double precision   one   ,         zero
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

*     Check the input parameters.

      if( n.eq.0 )then
         inform = 0
         return
      end if
      if( ( m.lt.n ).or.( n.lt.0 ).or.( ldA.lt.m ) )then
         inform = 1
         return
      end if

*     Perform the factorization.

      la = ldA
      do 20, k = 1, min( m - 1, n )

*        Use a Householder reflection to zero the kth column of A.
*        First set up the reflection.

         call dgrfg ( m - k, A( k, k ), A( k + 1, k ), 1, zero,
     $                zeta( k ) )
         if( ( zeta( k ).gt.zero ).and.( k.lt.n ) )then
            if( ( k + 1 ).eq.n )
     $         la = m - k + 1
            temp      = A( k, k )
            A( k, k ) = zeta( k )

*           We now perform the operation  A := Q( k )*A.

*           Let B denote the bottom ( m - k + 1 ) by ( n - k ) part
*           of A.

*           First form  work = B'*u. ( work is stored in the elements
*           zeta( k + 1 ), ..., zeta( n ). )

            call dgemv ( 'Transpose', m - k + 1, n - k,
     $                   one, A( k, k + 1 ), la, A( k, k ), 1,
     $                   zero, zeta( k + 1 ), 1 )

*           Now form  B := B - u*work'.

            call dger  ( m - k + 1, n - k, -one, A( k, k ), 1,
     $                   zeta( k + 1 ), 1, A( k, k + 1 ), la )

*           Restore beta.

            A( k, k ) = temp
         end if
   20 continue

*     Store the final zeta when m.eq.n.

      if( m.eq.n )
     $   zeta( n ) = zero

      inform = 0

*     end of dgeqr
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dgeqrp( pivot, m, n, A, ldA, zeta, perm, work, inform )
      character*1        pivot
      integer            m, n, ldA, inform
      integer            perm( * )
      double precision   A( ldA, * ), zeta( * ), work( * )

      double precision   wmach
      common    /solmch/ wmach(15)
      save      /solmch/

C  1. Purpose
C     =======
C
C  DGEQRP reduces the  m by n matrix A to upper triangular form by means
C  of orthogonal transformations and column permutations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )*P'      when   m.gt.n,
C           ( 0 )
C
C     A = Q*R*P'          when   m = n,
C
C     A = Q*( R  X )*P'   when   m.lt.n,
C
C  where  Q  is  an  m by m  orthogonal matrix, R  is a  min( m, n )  by
C  min( m, n )  upper triangular matrix and  P is an  n by n permutation
C  matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ),  which is used to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k )  and  z( k ) are chosen to annhilate the elements below the
C  triangular part of  A.
C
C  The vector  u( k )  is returned in the kth element of zeta and in the
C  kth column of A, such that zeta( k ) is in zeta( k ) and the elements
C  of z( k ) are in A( k + 1, k ), ..., A( m, k ). The elements of R are
C  returned in the upper triangular part of A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  where p = min( m - 1, n ).
C
C  Two options are available for the column permutations. In either case
C  the column for which the  sub-diagonal elements are to be annihilated
C  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
C  The  particular column chosen as the pivot column is either that  for
C  which  the  unreduced  part  ( elements k onwards )  has the  largest
C  Euclidean  length, or  is that for  which the ratio of the  Euclidean
C  length  of the  unreduced part  to the  Euclidean length of the whole
C  column is a maximum.
C
C  3. Parameters
C     ==========
C
C  PIVOT  - CHARACTER*1.
C
C           On  entry, PIVOT  specifies  the  pivoting  strategy  to  be
C           performed as follows.
C
C           PIVOT = 'C' or 'c'
C
C              Column  interchanges  are  to be  incorporated  into  the
C              factorization, such that the  column whose unreduced part
C              has  maximum  Euclidean  length  is chosen  as the  pivot
C              column at each step.
C
C           PIVOT = 'S' or 's'
C
C              Scaled  column interchanges  are to be  incorporated into
C              the  factorization, such  that the  column for which  the
C              ratio  of the  Euclidean  length of the unreduced part of
C              the column to the original Euclidean length of the column
C              is a maximum is chosen as the  pivot column at each step.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of A will contain the upper triangular matrix R and the
C           M by min( M, N )  strictly lower triangular part of  A  will
C           contain details  of the  factorization  as  described above.
C           When m.lt.n then the remaining M by ( N - M ) part of A will
C           contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must  specify  the leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  zeta   - 'real' array of DIMENSION at least ( n ).
C
C           On exit, zeta( k )  contains the scalar  zeta  for  the  kth
C           transformation. If T( k ) = I then zeta( k) = 0.0, otherwise
C           zeta( k )  contains the scalar  zeta( k ) as described above
C           and  is  always  in  the  range  ( 1.0, sqrt( 2.0 ) ).  When
C           n .gt. m  the  elements  zeta( m + 1 ),  zeta( m + 2 ), ...,
C           zeta( n )  are used as internal workspace.
C
C  PERM   - INTEGER array of DIMENSION at least min( m, n ).
C
C           On exit, PERM  contains details of the permutation matrix P,
C           such  that  PERM( k ) = k  if no  column interchange occured
C           at  the  kth  step  and  PERM( k ) = j, ( k .lt. j .le. n ),
C           if  columns  k and j  were  interchanged at  the  kth  step.
C           Note  that, although  there are  min( m - 1, n )  orthogonal
C           transformations, there are min( m, n ) permutations.
C
C  WORK   - 'real' array of DIMENSION at least ( 2*n ).
C
C           Used as internal workspace.
C
C           On exit, WORK( j ), j = 1, 2, ..., n, contains the Euclidean
C           length  of the  jth  column  of the  permuted  matrix  A*P'.
C
C  INFORM - INTEGER.
C
C           On  successful exit, INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly supplied. See the next section for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        PIVOT .ne. 'C' or 'c' or 'S' or 's'
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C     B := Q'*B   and   B := Q*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary  linear algebra  routine  DGEAPQ. The  operation  B := Q'*B
C  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'Transpose', 'Separate', M, N, A, LDA, zeta,
C    $             K, B, LDB, WORK, INFORM )
C
C  and  B := Q*B  can be obtained by the call:
C
C     INFORM = 0
C     CALL DGEAPQ( 'No transpose', 'Separate', M, N, A, LDA, zeta,
C    $             K, B, LDB, WORK, INFORM )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine DGEAPQ for further
C  details.
C
C  Also following the use of this routine the operations
C
C     B := P'*B   and   B := P*B,
C
C  where B is an n by k matrix, and the operations
C
C     B := B*P    and   B := B*P',
C
C  where  B is a k by n  matrix, can  be performed by calls to the basic
C  linear  algebra  routine  DGEAP .  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL DGEAP ( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL DGEAP ( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  N  in the above two calls.
C  The operation  B := B*P  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL DGEAP ( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine DGEAP for further details.
C
C  Operations involving  the matrix  R  are performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV.  Note that no test for near
C  singularity of  R is incorporated in this routine or in routine DTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  DUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular then  the  auxiliary  linear algebra  routine  DUTSV
C  can  be  used  to  determine  the  singular value decomposition of R.
C  Operations  involving  the  matrix  X  can also be  performed  by the
C  Level 2  BLAS  routines.  Matrices  of  the  form   ( R  X )  can  be
C  factorized as
C
C     ( R  X ) = ( T  0 )*S',
C
C  where  T is upper triangular and S is orthogonal, using the auxiliary
C  linear algebra routine  DUTRQ .
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-December-1984.
C     Sven Hammarling, Nag Central Office.
C
      external           mchpar, dgemv , dger  , dgrfg , dnrm2 , dswap
      intrinsic          abs   , max   , min   , sqrt
      integer            j     , jmax  , k     , la
      double precision   eps   , maxnrm, norm  , dnrm2 , temp  , tol
      double precision   lamda
      parameter        ( lamda = 1.0d-2 )
      double precision   one   ,         zero
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

*     Check the input parameters.

      if( min( m, n ).eq.0 )then
         inform = 0
         return
      end if
      if( ( ( pivot.ne.'C' ).and.( pivot.ne.'c' ).and.
     $      ( pivot.ne.'S' ).and.( pivot.ne.'s' )      ).or.
     $    ( m.lt.0 ).or.( n.lt.0 ).or.( ldA.lt.m )           )then
         inform = 1
         return
      end if

*     Compute eps and the initial column norms.

      call mchpar()
      eps = wmach( 3 )
      do 10, j = 1, n
         work( j )     = dnrm2 ( m, A( 1, j ), 1 )
         work( j + n ) = work( j )
   10 continue

*     Perform the factorization. TOL is the tolerance for DGRFG .

      la = ldA
      do 50, k = 1, min( m, n )

*        Find the pivot column.

         maxnrm = zero
         jmax   = k
         do 20, j = k, n
            if( ( pivot.eq.'C' ).or.( pivot.eq.'c' ) )then
               if( work( j + n  ).gt.maxnrm )then
                  maxnrm = work( j + n )
                  jmax   = j
               end if
            else if( work( j ).gt.zero )then
               if( ( work( j + n )/work( j ) ).gt.maxnrm )then
                  maxnrm = work( j + n )/work( j )
                  jmax   = j
               end if
            end if
   20    continue
         perm( k ) = jmax
         if( jmax.gt.k )then
            call dswap ( m, A( 1, k ), 1, A( 1, jmax ), 1 )
            temp             = work( k )
            work( k )        = work( jmax )
            work( jmax )     = temp
            work( jmax + n ) = work( k + n )
            perm( k )        = jmax
         end if
         tol = eps*work( k )
         if( k.lt.m )then

*           Use a Householder reflection to zero the kth column of A.
*           First set up the reflection.

            call dgrfg ( m - k, A( k, k ), A( k + 1, k ), 1, tol,
     $                   zeta( k ) )
            if( k.lt.n )then
               if( zeta( k ).gt.zero )then
                  if( ( k + 1 ).eq.n )
     $               la = m - k + 1
                  temp      = A( k, k )
                  A( k, k ) = zeta( k )

*                 We now perform the operation  A := Q( k )*A.

*                 Let B denote the bottom ( m - k + 1 ) by ( n - k )
*                 part of A.

*                 First form  work = B'*u. ( work is stored in the
*                 elements zeta( k + 1 ), ..., zeta( n ). )

                  call dgemv ( 'Transpose', m - k + 1, n - k,
     $                         one, A( k, k + 1 ), la, A( k, k ), 1,
     $                         zero, zeta( k + 1 ), 1 )

*                 Now form  B := B - u*work'.

                  call dger  ( m - k + 1, n - k, -one, A( k, k ), 1,
     $                         zeta( k + 1 ), 1, A( k, k + 1 ), la )

*                 Restore beta.

                  A( k, k ) = temp
               end if

*              Update the unreduced column norms. Use the Linpack
*              criterion for when to recompute the norms, except that
*              we retain the original column lengths throughout and use
*              a smaller lamda.

               do 40, j = k + 1, n
                  if( work( j + n ).gt.zero )then
                     temp = abs( A( k, j ) )/work( j + n )
                     temp = max( ( one + temp )*( one - temp ), zero )
                     norm = temp
                     temp = one +
     $                      lamda*temp*( work( j + n )/work( j ) )**2
                     if( temp.gt.one )then
                        work( j + n ) = work( j + n )*sqrt( norm )
                     else
                        work( j + n ) = dnrm2 ( m - k,
     $                                          A( k + 1, j ), 1 )
                     end if
                  end if
   40          continue
            end if
         end if
   50 continue

*     Store the final zeta when m.le.n.

      if( m.le.n )
     $   zeta( m ) = zero

      inform = 0

*     end of dgeqrp
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dgeap ( side, trans, m, n, perm, k, b, ldB )
*     .. Scalar Arguments ..
      integer            k, ldB, m, n
      character*1        side, trans
*     .. Array Arguments ..
      double precision   B( ldB, * )
      integer            perm( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEAP  performs one of the transformations
*
*     B := P'*B   or   B := P*B,   where B is an m by k matrix,
*
*  or
*
*     B := B*P'   or   B := B*P,   where B is a k by m matrix,
*
*  P being an m by m permutation matrix of the form
*
*     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
*
*  where  P( i, index( i ) ) is the permutation matrix that interchanges
*  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
*  with rows and columns  i and index( i )  interchanged.  Of course, if
*  index( i ) = i  then  P( i, index( i ) ) = I.
*
*  This routine  is intended for use in  conjunction with  Nag auxiliary
*  routines that  perform  interchange  operations,  such  as  pivoting.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*  TRANS
*           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
*           TRANS  ( Transpose, or No transpose )  specify the operation
*           to be performed as follows.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := P'*B.
*
*           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := P*B.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
*
*              Perform the operation   B := B*P'.
*
*           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
*
*              Perform the operation   B := B*P.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*
*           On entry, M must specify the order of the permutation matrix
*           P.  M must be at least zero.  When  M = 0  then an immediate
*           return is effected.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*
*           On entry,  N must specify the value of n. N must be at least
*           zero.  When  N = 0  then an  immediate  return is  effected.
*
*           Unchanged on exit.
*
*  PERM   - INTEGER array of DIMENSION at least ( n ).
*
*           Before  entry,  PERM  must  contain the  n  indices  for the
*           permutation matrices. index( i ) must satisfy
*
*              1 .le. index( i ) .le. m.
*
*           It is usual for index( i ) to be at least i, but this is not
*           necessary for this routine.
*
*           Unchanged on exit.
*
*  K      - INTEGER.
*
*           On entry with  SIDE = 'L' or 'l',  K must specify the number
*           of columns of B and on entry with  SIDE = 'R' or 'r', K must
*           specify the  number of rows of B.  K must be at least  zero.
*           When  K = 0  then an immediate return is effected.
*
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of  DIMENSION  ( LDB, ncolB ),  where
*           ncolB = k   when   SIDE = 'L' or 'l'  and   ncolB = m   when
*           SIDE = 'R' or 'r'.
*
*           Before entry  with  SIDE = 'L' or 'l',  the  leading  M by K
*           part  of  the  array   B  must  contain  the  matrix  to  be
*           transformed  and  before entry with  SIDE = 'R' or 'r',  the
*           leading  K by M part of the array  B must contain the matrix
*           to  be  transformed.  On  exit,  B  is  overwritten  by  the
*           transformed matrix.
*
*  LDB    - INTEGER.
*
*           On entry,  LDB  must specify  the  leading dimension  of the
*           array  B  as declared  in the  calling  (sub) program.  When
*           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
*           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
*           Unchanged on exit.
*
*
*  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
*
*  -- Written on 13-January-1986.
*     Sven Hammarling, Nag Central Office.
*
*
*     .. Local Scalars ..
      double precision   temp
      integer            i, j, l
      logical            left, null, right, trnsp
*     .. Intrinsic Functions ..
      intrinsic          min
*     ..
*     .. Executable Statements ..
      if( min( m, n, k ).eq.0 )
     $   return
      left  = ( side .eq.'L' ).or.( side .eq.'l' )
      right = ( side .eq.'R' ).or.( side .eq.'r' )
      null  = ( trans.eq.'N' ).or.( trans.eq.'n' )
      trnsp = ( trans.eq.'T' ).or.( trans.eq.'t' )
      if( left )then
         if( trnsp )then
            do 20, i = 1, n
               if( perm( i ).ne.i )then
                  l = perm( i )
                  do 10, j = 1, k
                     temp      = B( i, j )
                     B( i, j ) = B( l, j )
                     B( l, j ) = temp
   10             continue
               end if
   20       continue
         else if( null )then
            do 40, i = n, 1, -1
               if( perm( i ).ne.i )then
                  l = perm( i )
                  do 30, j = 1, k
                     temp      = B( l, j )
                     B( l, j ) = B( i, j )
                     B( i, j ) = temp
   30             continue
               end if
   40       continue
         end if
      else if( right )then
         if( trnsp )then
            do 60, j = 1, n
               if( perm( j ).ne.j )then
                  l = perm( j )
                  do 50, i = 1, k
                     temp      = B( i, j )
                     B( i, j ) = B( l, j )
                     B( l, j ) = temp
   50             continue
               end if
   60       continue
         else if( null )then
            do 80, j = n, 1, -1
               if( perm( j ).ne.j )then
                  l = perm( j )
                  do 70, i = 1, k
                     temp      = B( l, j )
                     B( l, j ) = B( i, j )
                     B( i, j ) = temp
   70             continue
               end if
   80       continue
         end if
      end if
*
*     end of dgeap (f06qjf)
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine dgeapq( trans, wherez, m, n, A, ldA, zeta,
     $                   ncolB, b, ldB, work, inform )
      character*1        trans, wherez
      integer            m, n, ldA, ncolB, ldB, inform
      double precision   A( ldA, * ), zeta( * ), B( ldB, * ), work( * )
C
C  1. Purpose
C     =======
C
C  DGEAPQ performs one of the transformations
C
C     B := Q'*B   or   B := Q*B,
C
C  where B is an m by ncolB matrix and Q is an m by m orthogonal matrix,
C  given as the product of  Householder transformation matrices, details
C  of  which are stored in the  m by n ( m.ge.n )  array  A  and, if the
C  parameter  WHEREZ = 'S' or 's', in the array zeta.
C
C  This  routine is  intended for use following auxiliary linear algebra
C  routines such as  DGEQR , DGEHES and DSLTRI. ( See those routines for
C  example calls. )
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta( k ) ),
C                                             (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  A( k + 1, k ), ..., A( m, k )  and  zeta( k ) must be supplied either
C  in  A( k, k )  or in  zeta( k ), depending upon the parameter WHEREZ.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining Q'.
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = ' ' or 'N' or 'n'
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'T' or 't' or 'C' or 'c'
C
C              Perform the operation  B := Q'*B.
C
C           Unchanged on exit.
C
C  WHEREZ - CHARACTER*1.
C
C           On entry, WHEREZ specifies where the elements of zeta are to
C           be found as follows.
C
C           WHEREZ = 'I' or 'i'
C
C              The elements of zeta are in A.
C
C           WHEREZ = 'S' or 's'
C
C              The elements of zeta are separate from A, in zeta.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - 'real' array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHEREZ = 'I' or 'i'  then  the  diagonal
C           elements of A must contain the elements of zeta.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  zeta   - 'real' array of DIMENSION at least min( m - 1, n ).
C
C           Before entry with  WHEREZ = 'S' or 's', the array  zeta must
C           contain the elements of the vector  zeta.
C
C           When  WHEREZ = 'I' or 'i', the array zeta is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - 'real' array of DIMENSION ( LDB, ncolB ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - 'real' array of DIMENSION at least ( ncolB ).
C
C           Used as internal workspace.
C
C  INFORM - INTEGER.
C
C           On  successful exit  INFORM  will be zero, otherwise  INFORM
C           will  be set to unity indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  INFORM = 1
C
C     One or more of the following conditions holds:
C
C        TRANS  .ne. ' ' or 'N' or 'n' or 'T' or 't' or 'C' or 'c'
C        WHEREZ .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 15-November-1984.
C     Sven Hammarling, Nag Central Office.
C
      external           dgemv , dger
      intrinsic          min
      integer            j     , k     , kk    , lb
      double precision   temp
      double precision   one   ,         zero
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

*     Check the input parameters.

      if( min( n, ncolB ).eq.0 )then
         inform = 0
         return
      end if
      if( ( ( trans .ne.' ' ).and.
     $      ( trans .ne.'N' ).and.( trans .ne.'n' ).and.
     $      ( trans .ne.'T' ).and.( trans .ne.'t' ).and.
     $      ( trans .ne.'C' ).and.( trans .ne.'c' )      ).or.
     $    ( ( wherez.ne.'I' ).and.( wherez.ne.'i' ).and.
     $      ( wherez.ne.'S' ).and.( wherez.ne.'s' )      ).or.
     $    ( m.lt.n ).or.( n.lt.0 ).or.( ldA.lt.m ).or.
     $    ( ncolB.lt.0 ).or.( ldB.lt.m )                      )then
         inform = 1
         return
      end if

*     Perform the transformation.

      lb = ldB
      do 20, kk = 1, min( m - 1, n )
         if( ( trans.eq.'T' ).or.( trans.eq.'t' ).or.
     $       ( trans.eq.'C' ).or.( trans.eq.'c' )     )then

*           Q'*B = Q( p )*...*Q( 2 )*Q( 1 )*B,     p = min( m - 1, n ).

            k = kk
         else

*           Q*B  = Q( 1 )'*Q( 2 )'*...*Q( p )'*B,  p = min( m - 1, n ).
*           Note that  Q( k )' = Q( k ).

            k = min( n, m - 1 ) + 1 - kk
         end if
         if( ( wherez.eq.'S' ).or.( wherez.eq.'s' ) )then
            temp      = A( k, k )
            A( k, k ) = zeta( k )
         end if

*        If zeta( k ) is zero then Q( k ) = I and we can skip the kth
*        transformation.

         if( A( k, k ).gt.zero )then
            if( ncolB.eq.1 )
     $         lb = m - k + 1

*           Let C denote the bottom ( m - k + 1 ) by ncolB part of B.

*           First form  work = C'*u.

            do 10, j = 1, ncolB
               work( j ) = zero
   10       continue
            call dgemv ( 'Transpose', m - k + 1, ncolB,
     $                   one, B( k, 1 ), lb, A( k, k ), 1,
     $                   zero, work, 1 )

*           Now form  C := C - u*work'.

            call dger  ( m - k + 1, ncolB, -one, A( k, k ), 1,
     $                   work, 1, B( k, 1 ), lb )
         end if

*        Restore the diagonal element of A.

         if( ( wherez.eq.'S' ).or.( wherez.eq.'s' ) )
     $      A( k, k ) = temp
   20 continue

      inform = 0

*     end of dgeapq
      end


      SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
     $                   INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEBAK forms the right or left eigenvectors of a real general matrix
*  by backward transformation on the computed eigenvectors of the
*  balanced matrix output by DGEBAL.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the type of backward transformation required:
*          = 'N', do nothing, return immediately;
*          = 'P', do backward transformation for permutation only;
*          = 'S', do backward transformation for scaling only;
*          = 'B', do backward transformations for both permutation and
*                 scaling.
*          JOB must be the same as the argument JOB supplied to DGEBAL.
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  V contains right eigenvectors;
*          = 'L':  V contains left eigenvectors.
*
*  N       (input) INTEGER
*          The number of rows of the matrix V.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          The integers ILO and IHI determined by DGEBAL.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  SCALE   (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutation and scaling factors, as returned
*          by DGEBAL.
*
*  M       (input) INTEGER
*          The number of columns of the matrix V.  M >= 0.
*
*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)
*          On entry, the matrix of right or left eigenvectors to be
*          transformed, as returned by DHSEIN or DTREVC.
*          On exit, V is overwritten by the transformed eigenvectors.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V. LDV >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEBAK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( M.EQ.0 )
     $   RETURN
      IF( LSAME( JOB, 'N' ) )
     $   RETURN
*
      IF( ILO.EQ.IHI )
     $   GO TO 30
*
*     Backward balance
*
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               S = SCALE( I )
               CALL DSCAL( M, S, V( I, 1 ), LDV )
   10       CONTINUE
         END IF
*
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               S = ONE / SCALE( I )
               CALL DSCAL( M, S, V( I, 1 ), LDV )
   20       CONTINUE
         END IF
*
      END IF
*
*     Backward permutation
*
*     For  I = ILO-1 step -1 until 1,
*              IHI+1 step 1 until N do --
*
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
         IF( RIGHTV ) THEN
            DO 40 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI )
     $            GO TO 40
               IF( I.LT.ILO )
     $            I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I )
     $            GO TO 40
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
         END IF
*
         IF( LEFTV ) THEN
            DO 50 II = 1, N
               I = II
               IF( I.GE.ILO .AND. I.LE.IHI )
     $            GO TO 50
               IF( I.LT.ILO )
     $            I = ILO - II
               K = SCALE( I )
               IF( K.EQ.I )
     $            GO TO 50
               CALL DSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   50       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEBAK
*
      END
      SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEBAL balances a general real matrix A.  This involves, first,
*  permuting A by a similarity transformation to isolate eigenvalues
*  in the first 1 to ILO-1 and last IHI+1 to N elements on the
*  diagonal; and second, applying a diagonal similarity transformation
*  to rows and columns ILO to IHI to make the rows and columns as
*  close in norm as possible.  Both steps are optional.
*
*  Balancing may reduce the 1-norm of the matrix, and improve the
*  accuracy of the computed eigenvalues and/or eigenvectors.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the operations to be performed on A:
*          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
*                  for i = 1,...,N;
*          = 'P':  permute only;
*          = 'S':  scale only;
*          = 'B':  both permute and scale.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the input matrix A.
*          On exit,  A is overwritten by the balanced matrix.
*          If JOB = 'N', A is not referenced.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  ILO     (output) INTEGER
*  IHI     (output) INTEGER
*          ILO and IHI are set to integers such that on exit
*          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
*
*  SCALE   (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied to
*          A.  If P(j) is the index of the row and column interchanged
*          with row and column j and D(j) is the scaling factor
*          applied to row and column j, then
*          SCALE(j) = P(j)    for j = 1,...,ILO-1
*                   = D(j)    for j = ILO,...,IHI
*                   = P(j)    for j = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The permutations consist of row and column interchanges which put
*  the matrix in the form
*
*             ( T1   X   Y  )
*     P A P = (  0   B   Z  )
*             (  0   0   T2 )
*
*  where T1 and T2 are upper triangular matrices whose eigenvalues lie
*  along the diagonal.  The column indices ILO and IHI mark the starting
*  and ending columns of the submatrix B. Balancing consists of applying
*  a diagonal similarity transformation inv(D) * B * D to make the
*  1-norms of each row of B and its corresponding column nearly equal.
*  The output matrix is
*
*     ( T1     X*D          Y    )
*     (  0  inv(D)*B*D  inv(D)*Z ).
*     (  0      0           T2   )
*
*  Information about the permutations P and the diagonal matrix D is
*  returned in the vector SCALE.
*
*  This subroutine is based on the EISPACK routine BALANC.
*
*  Modified by Tzu-Yi Chen, Computer Science Division, University of
*    California at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 2.0D+0 )
      DOUBLE PRECISION   FACTOR
      PARAMETER          ( FACTOR = 0.95D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1,
     $                   SFMIN2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEBAL', -INFO )
         RETURN
      END IF
*
      K = 1
      L = N
*
      IF( N.EQ.0 )
     $   GO TO 210
*
      IF( LSAME( JOB, 'N' ) ) THEN
         DO 10 I = 1, N
            SCALE( I ) = ONE
   10    CONTINUE
         GO TO 210
      END IF
*
      IF( LSAME( JOB, 'S' ) )
     $   GO TO 120
*
*     Permutation to isolate eigenvalues if possible
*
      GO TO 50
*
*     Row and column exchange.
*
   20 CONTINUE
      SCALE( M ) = J
      IF( J.EQ.M )
     $   GO TO 30
*
      CALL DSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL DSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
*
   30 CONTINUE
      GO TO ( 40, 80 )IEXC
*
*     Search for rows isolating an eigenvalue and push them down.
*
   40 CONTINUE
      IF( L.EQ.1 )
     $   GO TO 210
      L = L - 1
*
   50 CONTINUE
      DO 70 J = L, 1, -1
*
         DO 60 I = 1, L
            IF( I.EQ.J )
     $         GO TO 60
            IF( A( J, I ).NE.ZERO )
     $         GO TO 70
   60    CONTINUE
*
         M = L
         IEXC = 1
         GO TO 20
   70 CONTINUE
*
      GO TO 90
*
*     Search for columns isolating an eigenvalue and push them left.
*
   80 CONTINUE
      K = K + 1
*
   90 CONTINUE
      DO 110 J = K, L
*
         DO 100 I = K, L
            IF( I.EQ.J )
     $         GO TO 100
            IF( A( I, J ).NE.ZERO )
     $         GO TO 110
  100    CONTINUE
*
         M = K
         IEXC = 2
         GO TO 20
  110 CONTINUE
*
  120 CONTINUE
      DO 130 I = K, L
         SCALE( I ) = ONE
  130 CONTINUE
*
      IF( LSAME( JOB, 'P' ) )
     $   GO TO 210
*
*     Balance the submatrix in rows K to L.
*
*     Iterative loop for norm reduction
*
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
  140 CONTINUE
      NOCONV = .FALSE.
*
      DO 200 I = K, L
         C = ZERO
         R = ZERO
*
         DO 150 J = K, L
            IF( J.EQ.I )
     $         GO TO 150
            C = C + ABS( A( J, I ) )
            R = R + ABS( A( I, J ) )
  150    CONTINUE
         ICA = IDAMAX( L, A( 1, I ), 1 )
         CA = ABS( A( ICA, I ) )
         IRA = IDAMAX( N-K+1, A( I, K ), LDA )
         RA = ABS( A( I, IRA+K-1 ) )
*
*        Guard against zero C or R due to underflow.
*
         IF( C.EQ.ZERO .OR. R.EQ.ZERO )
     $      GO TO 200
         G = R / SCLFAC
         F = ONE
         S = C + R
  160    CONTINUE
         IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR.
     $       MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
         F = F*SCLFAC
         C = C*SCLFAC
         CA = CA*SCLFAC
         R = R / SCLFAC
         G = G / SCLFAC
         RA = RA / SCLFAC
         GO TO 160
*
  170    CONTINUE
         G = C / SCLFAC
  180    CONTINUE
         IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR.
     $       MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 190
         F = F / SCLFAC
         C = C / SCLFAC
         G = G / SCLFAC
         CA = CA / SCLFAC
         R = R*SCLFAC
         RA = RA*SCLFAC
         GO TO 180
*
*        Now balance.
*
  190    CONTINUE
         IF( ( C+R ).GE.FACTOR*S )
     $      GO TO 200
         IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
            IF( F*SCALE( I ).LE.SFMIN1 )
     $         GO TO 200
         END IF
         IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
            IF( SCALE( I ).GE.SFMAX1 / F )
     $         GO TO 200
         END IF
         G = ONE / F
         SCALE( I ) = SCALE( I )*F
         NOCONV = .TRUE.
*
         CALL DSCAL( N-K+1, G, A( I, K ), LDA )
         CALL DSCAL( L, F, A( 1, I ), 1 )
*
  200 CONTINUE
*
      IF( NOCONV )
     $   GO TO 140
*
  210 CONTINUE
      ILO = K
      IHI = L
*
      RETURN
*
*     End of DGEBAL
*
      END
      SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGECON estimates the reciprocal of the condition number of a general
*  real matrix A, in either the 1-norm or the infinity-norm, using
*  the LU factorization computed by DGETRF.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as
*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'O':  1-norm;
*          = 'I':         Infinity-norm.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  ANORM   (input) DOUBLE PRECISION
*          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*          If NORM = 'I', the infinity-norm of the original matrix A.
*
*  RCOND   (output) DOUBLE PRECISION
*          The reciprocal of the condition number of the matrix A,
*          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      DOUBLE PRECISION   AINVNM, SCALE, SL, SMLNUM, SU
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DLATRS, DRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGECON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Estimate the norm of inv(A).
*
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL DLACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(L).
*
            CALL DLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A,
     $                   LDA, WORK, SL, WORK( 2*N+1 ), INFO )
*
*           Multiply by inv(U).
*
            CALL DLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N,
     $                   A, LDA, WORK, SU, WORK( 3*N+1 ), INFO )
         ELSE
*
*           Multiply by inv(U').
*
            CALL DLATRS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N, A,
     $                   LDA, WORK, SU, WORK( 3*N+1 ), INFO )
*
*           Multiply by inv(L').
*
            CALL DLATRS( 'Lower', 'Transpose', 'Unit', NORMIN, N, A,
     $                   LDA, WORK, SL, WORK( 2*N+1 ), INFO )
         END IF
*
*        Divide X by 1/(SL*SU) if doing so will not cause overflow.
*
         SCALE = SL*SU
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
      RETURN
*
*     End of DGECON
*
      END
      SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEEQU computes row and column scalings intended to equilibrate an
*  M-by-N matrix A and reduce its condition number.  R returns the row
*  scale factors and C the column scale factors, chosen to try to make
*  the largest element in each row and column of the matrix B with
*  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
*
*  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
*  number and BIGNUM = largest safe number.  Use of these scaling
*  factors is not guaranteed to reduce the condition number of A but
*  works well in practice.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The M-by-N matrix whose equilibration factors are
*          to be computed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  R       (output) DOUBLE PRECISION array, dimension (M)
*          If INFO = 0 or INFO > M, R contains the row scale factors
*          for A.
*
*  C       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0,  C contains the column scale factors for A.
*
*  ROWCND  (output) DOUBLE PRECISION
*          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
*          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
*          AMAX is neither too large nor too small, it is not worth
*          scaling by R.
*
*  COLCND  (output) DOUBLE PRECISION
*          If INFO = 0, COLCND contains the ratio of the smallest
*          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
*          worth scaling by C.
*
*  AMAX    (output) DOUBLE PRECISION
*          Absolute value of largest matrix element.  If AMAX is very
*          close to overflow or very close to underflow, the matrix
*          should be scaled.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i,  and i is
*                <= M:  the i-th row of A is exactly zero
*                >  M:  the (i-M)-th column of A is exactly zero
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   BIGNUM, RCMAX, RCMIN, SMLNUM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
*     Get machine constants.
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
*     Compute row scale factors.
*
      DO 10 I = 1, M
         R( I ) = ZERO
   10 CONTINUE
*
*     Find the maximum element in each row.
*
      DO 30 J = 1, N
         DO 20 I = 1, M
            R( I ) = MAX( R( I ), ABS( A( I, J ) ) )
   20    CONTINUE
   30 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 50 I = 1, M
            IF( R( I ).EQ.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   50    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
*
*        Compute ROWCND = min(R(I)) / max(R(I))
*
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
*     Compute column scale factors
*
      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE
*
*     Find the maximum element in each column,
*     assuming the row scaling computed above.
*
      DO 90 J = 1, N
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), ABS( A( I, J ) )*R( I ) )
   80    CONTINUE
   90 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 110 J = 1, N
            IF( C( J ).EQ.ZERO ) THEN
               INFO = M + J
               RETURN
            END IF
  110    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
*
*        Compute COLCND = min(C(J)) / max(C(J))
*
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
      RETURN
*
*     End of DGEEQU
*
      END
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
*  eigenvalues and, optionally, the left and/or right eigenvectors.
*
*  The right eigenvector v(j) of A satisfies
*                   A * v(j) = lambda(j) * v(j)
*  where lambda(j) is its eigenvalue.
*  The left eigenvector u(j) of A satisfies
*                u(j)**H * A = lambda(j) * u(j)**H
*  where u(j)**H denotes the conjugate transpose of u(j).
*
*  The computed eigenvectors are normalized to have Euclidean norm
*  equal to 1 and largest component real.
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N': left eigenvectors of A are not computed;
*          = 'V': left eigenvectors of A are computed.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N': right eigenvectors of A are not computed;
*          = 'V': right eigenvectors of A are computed.
*
*  N       (input) INTEGER
*          The order of the matrix A. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.
*          On exit, A has been overwritten.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues.  Complex
*          conjugate pairs of eigenvalues appear consecutively
*          with the eigenvalue having the positive imaginary part
*          first.
*
*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
*          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*          after another in the columns of VL, in the same order
*          as their eigenvalues.
*          If JOBVL = 'N', VL is not referenced.
*          If the j-th eigenvalue is real, then u(j) = VL(:,j),
*          the j-th column of VL.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
*          u(j+1) = VL(:,j) - i*VL(:,j+1).
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.  LDVL >= 1; if
*          JOBVL = 'V', LDVL >= N.
*
*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
*          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*          after another in the columns of VR, in the same order
*          as their eigenvalues.
*          If JOBVR = 'N', VR is not referenced.
*          If the j-th eigenvalue is real, then v(j) = VR(:,j),
*          the j-th column of VR.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
*          v(j+1) = VR(:,j) - i*VR(:,j+1).
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= 1; if
*          JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,3*N), and
*          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
*          performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the QR algorithm failed to compute all the
*                eigenvalues, and no eigenvectors have been computed;
*                elements i+1:N of WR and WI contain eigenvalues which
*                have converged.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K,
     $                   MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM,
     $                   SN
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLABAD, DLACPY,
     $                   DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC,
     $                   XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2,
     $                   DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -9
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -11
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by DHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
            IF( WANTVL ) THEN
               MINWRK = 4*N
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1,
     $                       'DORGHR', ' ', N, 1, N, -1 ) )
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL,
     $                WORK, -1, INFO )
               HSWORK = WORK( 1 )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
               MAXWRK = MAX( MAXWRK, 4*N )
            ELSE IF( WANTVR ) THEN
               MINWRK = 4*N
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1,
     $                       'DORGHR', ' ', N, 1, N, -1 ) )
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR,
     $                WORK, -1, INFO )
               HSWORK = WORK( 1 )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
               MAXWRK = MAX( MAXWRK, 4*N )
            ELSE 
               MINWRK = 3*N
               CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR,
     $                WORK, -1, INFO )
               HSWORK = WORK( 1 )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
            END IF
            MAXWRK = MAX( MAXWRK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $   CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*
*     Balance the matrix
*     (Workspace: need N)
*
      IBAL = 1
      CALL DGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
*
*     Reduce to upper Hessenberg form
*     (Workspace: need 3*N, prefer 2*N+N*NB)
*
      ITAU = IBAL + N
      IWRK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
     $             LWORK-IWRK+1, IERR )
*
      IF( WANTVL ) THEN
*
*        Want left eigenvectors
*        Copy Householder vectors to VL
*
         SIDE = 'L'
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
*
*        Generate orthogonal matrix in VL
*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
*
         IF( WANTVR ) THEN
*
*           Want left and right eigenvectors
*           Copy Schur vectors to VR
*
            SIDE = 'B'
            CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
*
      ELSE IF( WANTVR ) THEN
*
*        Want right eigenvectors
*        Copy Householder vectors to VR
*
         SIDE = 'R'
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
*
*        Generate orthogonal matrix in VR
*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
*
      ELSE
*
*        Compute eigenvalues only
*        (Workspace: need N+1, prefer N+HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
*
*     If INFO > 0 from DHSEQR, then quit
*
      IF( INFO.GT.0 )
     $   GO TO 50
*
      IF( WANTVL .OR. WANTVR ) THEN
*
*        Compute left and/or right eigenvectors
*        (Workspace: need 4*N)
*
         CALL DTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                N, NOUT, WORK( IWRK ), IERR )
      END IF
*
      IF( WANTVL ) THEN
*
*        Undo balancing of left eigenvectors
*        (Workspace: need N)
*
         CALL DGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL,
     $                IERR )
*
*        Normalize left eigenvectors and make largest component real
*
         DO 20 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ),
     $               DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = IDAMAX( N, WORK( IWRK ), 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   20    CONTINUE
      END IF
*
      IF( WANTVR ) THEN
*
*        Undo balancing of right eigenvectors
*        (Workspace: need N)
*
         CALL DGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR,
     $                IERR )
*
*        Normalize right eigenvectors and make largest component real
*
         DO 40 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ),
     $               DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK( IWRK ), 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
   40    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N,
     $                   IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N,
     $                   IERR )
         END IF
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
*     End of DGEEV
*
      END
      SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,
     $                   VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          BALANC, JOBVL, JOBVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ),
     $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEEVX computes for an N-by-N real nonsymmetric matrix A, the
*  eigenvalues and, optionally, the left and/or right eigenvectors.
*
*  Optionally also, it computes a balancing transformation to improve
*  the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
*  SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
*  (RCONDE), and reciprocal condition numbers for the right
*  eigenvectors (RCONDV).
*
*  The right eigenvector v(j) of A satisfies
*                   A * v(j) = lambda(j) * v(j)
*  where lambda(j) is its eigenvalue.
*  The left eigenvector u(j) of A satisfies
*                u(j)**H * A = lambda(j) * u(j)**H
*  where u(j)**H denotes the conjugate transpose of u(j).
*
*  The computed eigenvectors are normalized to have Euclidean norm
*  equal to 1 and largest component real.
*
*  Balancing a matrix means permuting the rows and columns to make it
*  more nearly upper triangular, and applying a diagonal similarity
*  transformation D * A * D**(-1), where D is a diagonal matrix, to
*  make its rows and columns closer in norm and the condition numbers
*  of its eigenvalues and eigenvectors smaller.  The computed
*  reciprocal condition numbers correspond to the balanced matrix.
*  Permuting rows and columns will not change the condition numbers
*  (in exact arithmetic) but diagonal scaling will.  For further
*  explanation of balancing, see section 4.10.2 of the LAPACK
*  Users' Guide.
*
*  Arguments
*  =========
*
*  BALANC  (input) CHARACTER*1
*          Indicates how the input matrix should be diagonally scaled
*          and/or permuted to improve the conditioning of its
*          eigenvalues.
*          = 'N': Do not diagonally scale or permute;
*          = 'P': Perform permutations to make the matrix more nearly
*                 upper triangular. Do not diagonally scale;
*          = 'S': Diagonally scale the matrix, i.e. replace A by
*                 D*A*D**(-1), where D is a diagonal matrix chosen
*                 to make the rows and columns of A more equal in
*                 norm. Do not permute;
*          = 'B': Both diagonally scale and permute A.
*
*          Computed reciprocal condition numbers will be for the matrix
*          after balancing and/or permuting. Permuting does not change
*          condition numbers (in exact arithmetic), but balancing does.
*
*  JOBVL   (input) CHARACTER*1
*          = 'N': left eigenvectors of A are not computed;
*          = 'V': left eigenvectors of A are computed.
*          If SENSE = 'E' or 'B', JOBVL must = 'V'.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N': right eigenvectors of A are not computed;
*          = 'V': right eigenvectors of A are computed.
*          If SENSE = 'E' or 'B', JOBVR must = 'V'.
*
*  SENSE   (input) CHARACTER*1
*          Determines which reciprocal condition numbers are computed.
*          = 'N': None are computed;
*          = 'E': Computed for eigenvalues only;
*          = 'V': Computed for right eigenvectors only;
*          = 'B': Computed for eigenvalues and right eigenvectors.
*
*          If SENSE = 'E' or 'B', both left and right eigenvectors
*          must also be computed (JOBVL = 'V' and JOBVR = 'V').
*
*  N       (input) INTEGER
*          The order of the matrix A. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.
*          On exit, A has been overwritten.  If JOBVL = 'V' or
*          JOBVR = 'V', A contains the real Schur form of the balanced
*          version of the input matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues.  Complex
*          conjugate pairs of eigenvalues will appear consecutively
*          with the eigenvalue having the positive imaginary part
*          first.
*
*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
*          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*          after another in the columns of VL, in the same order
*          as their eigenvalues.
*          If JOBVL = 'N', VL is not referenced.
*          If the j-th eigenvalue is real, then u(j) = VL(:,j),
*          the j-th column of VL.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
*          u(j+1) = VL(:,j) - i*VL(:,j+1).
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.  LDVL >= 1; if
*          JOBVL = 'V', LDVL >= N.
*
*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
*          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*          after another in the columns of VR, in the same order
*          as their eigenvalues.
*          If JOBVR = 'N', VR is not referenced.
*          If the j-th eigenvalue is real, then v(j) = VR(:,j),
*          the j-th column of VR.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
*          v(j+1) = VR(:,j) - i*VR(:,j+1).
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= 1, and if
*          JOBVR = 'V', LDVR >= N.
*
*  ILO     (output) INTEGER
*  IHI     (output) INTEGER
*          ILO and IHI are integer values determined when A was
*          balanced.  The balanced A(i,j) = 0 if I > J and
*          J = 1,...,ILO-1 or I = IHI+1,...,N.
*
*  SCALE   (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          when balancing A.  If P(j) is the index of the row and column
*          interchanged with row and column j, and D(j) is the scaling
*          factor applied to row and column j, then
*          SCALE(J) = P(J),    for J = 1,...,ILO-1
*                   = D(J),    for J = ILO,...,IHI
*                   = P(J)     for J = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  ABNRM   (output) DOUBLE PRECISION
*          The one-norm of the balanced matrix (the maximum
*          of the sum of absolute values of elements of any column).
*
*  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
*          RCONDE(j) is the reciprocal condition number of the j-th
*          eigenvalue.
*
*  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
*          RCONDV(j) is the reciprocal condition number of the j-th
*          right eigenvector.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.   If SENSE = 'N' or 'E',
*          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',
*          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (2*N-2)
*          If SENSE = 'N' or 'E', not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the QR algorithm failed to compute all the
*                eigenvalues, and no eigenvectors or condition numbers
*                have been computed; elements 1:ILO-1 and i+1:N of WR
*                and WI contain eigenvalues which have converged.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE,
     $                   WNTSNN, WNTSNV
      CHARACTER          JOB, SIDE
      INTEGER            HSWORK, I, ICOND, IERR, ITAU, IWRK, K, MAXWRK,
     $                   MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM,
     $                   SN
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLABAD, DLACPY,
     $                   DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC,
     $                   DTRSNA, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2,
     $                   DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      WNTSNN = LSAME( SENSE, 'N' )
      WNTSNE = LSAME( SENSE, 'E' )
      WNTSNV = LSAME( SENSE, 'V' )
      WNTSNB = LSAME( SENSE, 'B' )
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. LSAME( BALANC,
     $    'S' ) .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' ) ) )
     $     THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( WNTSNN .OR. WNTSNE .OR. WNTSNB .OR. WNTSNV ) .OR.
     $         ( ( WNTSNE .OR. WNTSNB ) .AND. .NOT.( WANTVL .AND.
     $         WANTVR ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by DHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
*
            IF( WANTVL ) THEN
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL,
     $                WORK, -1, INFO )
            ELSE IF( WANTVR ) THEN
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR,
     $                WORK, -1, INFO )
            ELSE
               IF( WNTSNN ) THEN
                  CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR,
     $                LDVR, WORK, -1, INFO )
               ELSE
                  CALL DHSEQR( 'S', 'N', N, 1, N, A, LDA, WR, WI, VR,
     $                LDVR, WORK, -1, INFO )
               END IF
            END IF
            HSWORK = WORK( 1 )
*
            IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
               MINWRK = 2*N
               IF( .NOT.WNTSNN )
     $            MINWRK = MAX( MINWRK, N*N+6*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               IF( .NOT.WNTSNN )
     $            MAXWRK = MAX( MAXWRK, N*N + 6*N )
            ELSE
               MINWRK = 3*N
               IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) )
     $            MINWRK = MAX( MINWRK, N*N + 6*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'DORGHR',
     $                       ' ', N, 1, N, -1 ) )
               IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) )
     $            MAXWRK = MAX( MAXWRK, N*N + 6*N )
               MAXWRK = MAX( MAXWRK, 3*N )
            END IF
            MAXWRK = MAX( MAXWRK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -21
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ICOND = 0
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA )
     $   CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*
*     Balance the matrix and compute ABNRM
*
      CALL DGEBAL( BALANC, N, A, LDA, ILO, IHI, SCALE, IERR )
      ABNRM = DLANGE( '1', N, N, A, LDA, DUM )
      IF( SCALEA ) THEN
         DUM( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
         ABNRM = DUM( 1 )
      END IF
*
*     Reduce to upper Hessenberg form
*     (Workspace: need 2*N, prefer N+N*NB)
*
      ITAU = 1
      IWRK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
     $             LWORK-IWRK+1, IERR )
*
      IF( WANTVL ) THEN
*
*        Want left eigenvectors
*        Copy Householder vectors to VL
*
         SIDE = 'L'
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
*
*        Generate orthogonal matrix in VL
*        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (Workspace: need 1, prefer HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
*
         IF( WANTVR ) THEN
*
*           Want left and right eigenvectors
*           Copy Schur vectors to VR
*
            SIDE = 'B'
            CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
*
      ELSE IF( WANTVR ) THEN
*
*        Want right eigenvectors
*        Copy Householder vectors to VR
*
         SIDE = 'R'
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
*
*        Generate orthogonal matrix in VR
*        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (Workspace: need 1, prefer HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
*
      ELSE
*
*        Compute eigenvalues only
*        If condition numbers desired, compute Schur form
*
         IF( WNTSNN ) THEN
            JOB = 'E'
         ELSE
            JOB = 'S'
         END IF
*
*        (Workspace: need 1, prefer HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( JOB, 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
*
*     If INFO > 0 from DHSEQR, then quit
*
      IF( INFO.GT.0 )
     $   GO TO 50
*
      IF( WANTVL .OR. WANTVR ) THEN
*
*        Compute left and/or right eigenvectors
*        (Workspace: need 3*N)
*
         CALL DTREVC( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                N, NOUT, WORK( IWRK ), IERR )
      END IF
*
*     Compute condition numbers if desired
*     (Workspace: need N*N+6*N unless SENSE = 'E')
*
      IF( .NOT.WNTSNN ) THEN
         CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, IWORK,
     $                ICOND )
      END IF
*
      IF( WANTVL ) THEN
*
*        Undo balancing of left eigenvectors
*
         CALL DGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL,
     $                IERR )
*
*        Normalize left eigenvectors and make largest component real
*
         DO 20 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ),
     $               DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( K ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   20    CONTINUE
      END IF
*
      IF( WANTVR ) THEN
*
*        Undo balancing of right eigenvectors
*
         CALL DGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR,
     $                IERR )
*
*        Normalize right eigenvectors and make largest component real
*
         DO 40 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ),
     $               DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( K ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
   40    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ),
     $                MAX( N-INFO, 1 ), IERR )
         IF( INFO.EQ.0 ) THEN
            IF( ( WNTSNV .OR. WNTSNB ) .AND. ICOND.EQ.0 )
     $         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N,
     $                      IERR )
         ELSE
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N,
     $                   IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N,
     $                   IERR )
         END IF
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
*     End of DGEEVX
*
      END
      SUBROUTINE DGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEHD2 reduces a real general matrix A to upper Hessenberg form H by
*  an orthogonal similarity transformation:  Q' * A * Q = H .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows
*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*          set by a previous call to DGEBAL; otherwise they should be
*          set to 1 and N respectively. See Further Details.
*          1 <= ILO <= IHI <= max(1,N).
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the n by n general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          elements below the first subdiagonal, with the array TAU,
*          represent the orthogonal matrix Q as a product of elementary
*          reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of (ihi-ilo) elementary
*  reflectors
*
*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*  exit in A(i+2:ihi,i), and tau in TAU(i).
*
*  The contents of A are illustrated by the following example, with
*  n = 7, ilo = 2 and ihi = 6:
*
*  on entry,                        on exit,
*
*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*  (                         a )    (                          a )
*
*  where a denotes an element of the original matrix A, h denotes a
*  modified element of the upper Hessenberg matrix H, and vi denotes an
*  element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEHD2', -INFO )
         RETURN
      END IF
*
      DO 10 I = ILO, IHI - 1
*
*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         CALL DLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                TAU( I ) )
         AII = A( I+1, I )
         A( I+1, I ) = ONE
*
*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL DLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ),
     $               A( 1, I+1 ), LDA, WORK )
*
*        Apply H(i) to A(i+1:ihi,i+1:n) from the left
*
         CALL DLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ),
     $               A( I+1, I+1 ), LDA, WORK )
*
         A( I+1, I ) = AII
   10 CONTINUE
*
      RETURN
*
*     End of DGEHD2
*
      END
      SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2.1)                                  --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2009                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEHRD reduces a real general matrix A to upper Hessenberg form H by
*  an orthogonal similarity transformation:  Q' * A * Q = H .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows
*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*          set by a previous call to DGEBAL; otherwise they should be
*          set to 1 and N respectively. See Further Details.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          elements below the first subdiagonal, with the array TAU,
*          represent the orthogonal matrix Q as a product of elementary
*          reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
*          zero.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of (ihi-ilo) elementary
*  reflectors
*
*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*  exit in A(i+2:ihi,i), and tau in TAU(i).
*
*  The contents of A are illustrated by the following example, with
*  n = 7, ilo = 2 and ihi = 6:
*
*  on entry,                        on exit,
*
*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*  (                         a )    (                          a )
*
*  where a denotes an element of the original matrix A, h denotes a
*  modified element of the upper Hessenberg matrix H, and vi denotes an
*  element of the vector defining H(i).
*
*  This file is a slight modification of LAPACK-3.0's DGEHRD
*  subroutine incorporating improvements proposed by Quintana-Orti and
*  Van de Geijn (2006). (See DLAHR2.)
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, 
     $                     ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, LDWORK, LWKOPT, NB,
     $                   NBMIN, NH, NX
      DOUBLE PRECISION  EI
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION  T( LDT, NBMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DGEHD2, DGEMM, DLAHR2, DLARFB, DTRMM,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      NB = MIN( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEHRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
*
      DO 10 I = 1, ILO - 1
         TAU( I ) = ZERO
   10 CONTINUE
      DO 20 I = MAX( 1, IHI ), N - 1
         TAU( I ) = ZERO
   20 CONTINUE
*
*     Quick return if possible
*
      NH = IHI - ILO + 1
      IF( NH.LE.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size
*
      NB = MIN( NBMAX, ILAENV( 1, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
      NBMIN = 2
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.NH ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code)
*
         NX = MAX( NB, ILAENV( 3, 'DGEHRD', ' ', N, ILO, IHI, -1 ) )
         IF( NX.LT.NH ) THEN
*
*           Determine if workspace is large enough for blocked code
*
            IWS = N*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code
*
               NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N, ILO, IHI,
     $                 -1 ) )
               IF( LWORK.GE.N*NBMIN ) THEN
                  NB = LWORK / N
               ELSE
                  NB = 1
               END IF
            END IF
         END IF
      END IF
      LDWORK = N
*
      IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
*
*        Use unblocked code below
*
         I = ILO
*
      ELSE
*
*        Use blocked code
*
         DO 40 I = ILO, IHI - 1 - NX, NB
            IB = MIN( NB, IHI-I )
*
*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
*           matrices V and T of the block reflector H = I - V*T*V'
*           which performs the reduction, and also the matrix Y = A*V*T
*
            CALL DLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), T, LDT,
     $                   WORK, LDWORK )
*
*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
*           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set
*           to 1
*
            EI = A( I+IB, I+IB-1 )
            A( I+IB, I+IB-1 ) = ONE
            CALL DGEMM( 'No transpose', 'Transpose', 
     $                  IHI, IHI-I-IB+1,
     $                  IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE,
     $                  A( 1, I+IB ), LDA )
            A( I+IB, I+IB-1 ) = EI
*
*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
*           right
*
            CALL DTRMM( 'Right', 'Lower', 'Transpose',
     $                  'Unit', I, IB-1,
     $                  ONE, A( I+1, I ), LDA, WORK, LDWORK )
            DO 30 J = 0, IB-2
               CALL DAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1,
     $                     A( 1, I+J+1 ), 1 )
   30       CONTINUE
*
*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
*           left
*
            CALL DLARFB( 'Left', 'Transpose', 'Forward',
     $                   'Columnwise',
     $                   IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, T, LDT,
     $                   A( I+1, I+IB ), LDA, WORK, LDWORK )
   40    CONTINUE
      END IF
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL DGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
      WORK( 1 ) = IWS
*
      RETURN
*
*     End of DGEHRD
*
      END
      SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEQR2 computes a QR factorization of a real m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL DLARFP( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i) to A(i:m,i+1:n) from the left
*
            AII = A( I, I )
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END
      SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
     $                   X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DGERFS improves the computed solution to a system of linear
*  equations and provides error bounds and backward error estimates for
*  the solution.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The original N-by-N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDAF    (input) INTEGER
*          The leading dimension of the array AF.  LDAF >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          The right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          On entry, the solution matrix X, as computed by DGETRS.
*          On exit, the improved solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Internal Parameters
*  ===================
*
*  ITMAX is the maximum number of steps of iterative refinement.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
      DOUBLE PRECISION   THREE
      PARAMETER          ( THREE = 3.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      CHARACTER          TRANST
      INTEGER            COUNT, I, J, K, KASE, NZ
      DOUBLE PRECISION   EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DGETRS, DLACN2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERFS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = N + 1
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 140 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A, A**T, or A**H, depending on TRANS.
*
         CALL DCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL DGEMV( TRANS, N, N, -ONE, A, LDA, X( 1, J ), 1, ONE,
     $               WORK( N+1 ), 1 )
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         DO 30 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   30    CONTINUE
*
*        Compute abs(op(A))*abs(X) + abs(B).
*
         IF( NOTRAN ) THEN
            DO 50 K = 1, N
               XK = ABS( X( K, J ) )
               DO 40 I = 1, N
                  WORK( I ) = WORK( I ) + ABS( A( I, K ) )*XK
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 70 K = 1, N
               S = ZERO
               DO 60 I = 1, N
                  S = S + ABS( A( I, K ) )*ABS( X( I, J ) )
   60          CONTINUE
               WORK( K ) = WORK( K ) + S
   70       CONTINUE
         END IF
         S = ZERO
         DO 80 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) /
     $             ( WORK( I )+SAFE1 ) )
            END IF
   80    CONTINUE
         BERR( J ) = S
*
*        Test stopping criterion. Continue iterating if
*           1) The residual BERR(J) is larger than machine epsilon, and
*           2) BERR(J) decreased by at least a factor of 2 during the
*              last iteration, and
*           3) At most ITMAX iterations tried.
*
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND.
     $       COUNT.LE.ITMAX ) THEN
*
*           Update solution and try again.
*
            CALL DGETRS( TRANS, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,
     $                   INFO )
            CALL DAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(op(A)))*
*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(op(A)) is the inverse of op(A)
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*        Use DLACN2 to estimate the infinity-norm of the matrix
*           inv(op(A)) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         DO 90 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
   90    CONTINUE
*
         KASE = 0
  100    CONTINUE
         CALL DLACN2( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ),
     $                KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(op(A)**T).
*
               CALL DGETRS( TRANST, N, 1, AF, LDAF, IPIV, WORK( N+1 ),
     $                      N, INFO )
               DO 110 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  110          CONTINUE
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               DO 120 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  120          CONTINUE
               CALL DGETRS( TRANS, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,
     $                      INFO )
            END IF
            GO TO 100
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO )
     $      FERR( J ) = FERR( J ) / LSTRES
*
  140 CONTINUE
*
      RETURN
*
*     End of DGERFS
*
      END
      SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. External Subroutines ..
      EXTERNAL           DGETRF, DGETRS, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESV ', -INFO )
         RETURN
      END IF
*
*     Compute the LU factorization of A.
*
      CALL DGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        Solve the system A*X = B, overwriting B with X.
*
         CALL DGETRS( 'No transpose', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
*
*     End of DGESV
*
      END
      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), C( * ), FERR( * ), R( * ),
     $                   WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DGESVX uses the LU factorization to compute the solution to a real
*  system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*  Description
*  ===========
*
*  The following steps are performed:
*
*  1. If FACT = 'E', real scaling factors are computed to equilibrate
*     the system:
*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
*     Whether or not the system will be equilibrated depends on the
*     scaling of the matrix A, but if equilibration is used, A is
*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
*     or diag(C)*B (if TRANS = 'T' or 'C').
*
*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
*     matrix A (after equilibration if FACT = 'E') as
*        A = P * L * U,
*     where P is a permutation matrix, L is a unit lower triangular
*     matrix, and U is upper triangular.
*
*  3. If some U(i,i)=0, so that U is exactly singular, then the routine
*     returns with INFO = i. Otherwise, the factored form of A is used
*     to estimate the condition number of the matrix A.  If the
*     reciprocal of the condition number is less than machine precision,
*     INFO = N+1 is returned as a warning, but the routine still goes on
*     to solve for X and compute error bounds as described below.
*
*  4. The system of equations is solved for X using the factored form
*     of A.
*
*  5. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  6. If equilibration was used, the matrix X is premultiplied by
*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
*     that it solves the original system before equilibration.
*
*  Arguments
*  =========
*
*  FACT    (input) CHARACTER*1
*          Specifies whether or not the factored form of the matrix A is
*          supplied on entry, and if not, whether the matrix A should be
*          equilibrated before it is factored.
*          = 'F':  On entry, AF and IPIV contain the factored form of A.
*                  If EQUED is not 'N', the matrix A has been
*                  equilibrated with scaling factors given by R and C.
*                  A, AF, and IPIV are not modified.
*          = 'N':  The matrix A will be copied to AF and factored.
*          = 'E':  The matrix A will be equilibrated if necessary, then
*                  copied to AF and factored.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Transpose)
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
*          not 'N', then A must have been equilibrated by the scaling
*          factors in R and/or C.  A is not modified if FACT = 'F' or
*          'N', or if FACT = 'E' and EQUED = 'N' on exit.
*
*          On exit, if EQUED .ne. 'N', A is scaled as follows:
*          EQUED = 'R':  A := diag(R) * A
*          EQUED = 'C':  A := A * diag(C)
*          EQUED = 'B':  A := diag(R) * A * diag(C).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)
*          If FACT = 'F', then AF is an input argument and on entry
*          contains the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then
*          AF is the factored form of the equilibrated matrix A.
*
*          If FACT = 'N', then AF is an output argument and on exit
*          returns the factors L and U from the factorization A = P*L*U
*          of the original matrix A.
*
*          If FACT = 'E', then AF is an output argument and on exit
*          returns the factors L and U from the factorization A = P*L*U
*          of the equilibrated matrix A (see the description of A for
*          the form of the equilibrated matrix).
*
*  LDAF    (input) INTEGER
*          The leading dimension of the array AF.  LDAF >= max(1,N).
*
*  IPIV    (input or output) INTEGER array, dimension (N)
*          If FACT = 'F', then IPIV is an input argument and on entry
*          contains the pivot indices from the factorization A = P*L*U
*          as computed by DGETRF; row i of the matrix was interchanged
*          with row IPIV(i).
*
*          If FACT = 'N', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization A = P*L*U
*          of the original matrix A.
*
*          If FACT = 'E', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization A = P*L*U
*          of the equilibrated matrix A.
*
*  EQUED   (input or output) CHARACTER*1
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration (always true if FACT = 'N').
*          = 'R':  Row equilibration, i.e., A has been premultiplied by
*                  diag(R).
*          = 'C':  Column equilibration, i.e., A has been postmultiplied
*                  by diag(C).
*          = 'B':  Both row and column equilibration, i.e., A has been
*                  replaced by diag(R) * A * diag(C).
*          EQUED is an input argument if FACT = 'F'; otherwise, it is an
*          output argument.
*
*  R       (input or output) DOUBLE PRECISION array, dimension (N)
*          The row scale factors for A.  If EQUED = 'R' or 'B', A is
*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
*          is not accessed.  R is an input argument if FACT = 'F';
*          otherwise, R is an output argument.  If FACT = 'F' and
*          EQUED = 'R' or 'B', each element of R must be positive.
*
*  C       (input or output) DOUBLE PRECISION array, dimension (N)
*          The column scale factors for A.  If EQUED = 'C' or 'B', A is
*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
*          is not accessed.  C is an input argument if FACT = 'F';
*          otherwise, C is an output argument.  If FACT = 'F' and
*          EQUED = 'C' or 'B', each element of C must be positive.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit,
*          if EQUED = 'N', B is not modified;
*          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
*          diag(R)*B;
*          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
*          overwritten by diag(C)*B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
*          to the original system of equations.  Note that A and B are
*          modified on exit if EQUED .ne. 'N', and the solution to the
*          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
*          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
*          and EQUED = 'R' or 'B'.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  RCOND   (output) DOUBLE PRECISION
*          The estimate of the reciprocal condition number of the matrix
*          A after equilibration (if done).  If RCOND is less than the
*          machine precision (in particular, if RCOND = 0), the matrix
*          is singular to working precision.  This condition is
*          indicated by a return code of INFO > 0.
*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (4*N)
*          On exit, WORK(1) contains the reciprocal pivot growth
*          factor norm(A)/norm(U). The "max absolute element" norm is
*          used. If WORK(1) is much less than 1, then the stability
*          of the LU factorization of the (equilibrated) matrix A
*          could be poor. This also means that the solution X, condition
*          estimator RCOND, and forward error bound FERR could be
*          unreliable. If factorization fails with 0<INFO<=N, then
*          WORK(1) contains the reciprocal pivot growth factor for the
*          leading INFO columns of A.
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, and i is
*                <= N:  U(i,i) is exactly zero.  The factorization has
*                       been completed, but the factor U is exactly
*                       singular, so the solution and error bounds
*                       could not be computed. RCOND = 0 is returned.
*                = N+1: U is nonsingular, but RCOND is less than machine
*                       precision, meaning that the matrix is singular
*                       to working precision.  Nevertheless, the
*                       solution and error bounds are computed because
*                       there are a number of situations where the
*                       computed solution can be more accurate than the
*                       value of RCOND would suggest.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU
      CHARACTER          NORM
      INTEGER            I, INFEQU, J
      DOUBLE PRECISION   AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN,
     $                   ROWCND, RPVGRW, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANTR
      EXTERNAL           LSAME, DLAMCH, DLANGE, DLANTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGECON, DGEEQU, DGERFS, DGETRF, DGETRS, DLACPY,
     $                   DLAQGE, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      ELSE
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF
*
*     Test the input parameters.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) )
     $     THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT.
     $         ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -10
      ELSE
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
   10       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
   20       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -12
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               COLCND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -14
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -16
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESVX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL DGEEQU( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL DLAQGE( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
      END IF
*
*     Scale the right hand side.
*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) THEN
            DO 40 J = 1, NRHS
               DO 30 I = 1, N
                  B( I, J ) = R( I )*B( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( COLEQU ) THEN
         DO 60 J = 1, NRHS
            DO 50 I = 1, N
               B( I, J ) = C( I )*B( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the LU factorization of A.
*
         CALL DLACPY( 'Full', N, N, A, LDA, AF, LDAF )
         CALL DGETRF( N, N, AF, LDAF, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.GT.0 ) THEN
*
*           Compute the reciprocal pivot growth factor of the
*           leading rank-deficient INFO columns of A.
*
            RPVGRW = DLANTR( 'M', 'U', 'N', INFO, INFO, AF, LDAF,
     $               WORK )
            IF( RPVGRW.EQ.ZERO ) THEN
               RPVGRW = ONE
            ELSE
               RPVGRW = DLANGE( 'M', N, INFO, A, LDA, WORK ) / RPVGRW
            END IF
            WORK( 1 ) = RPVGRW
            RCOND = ZERO
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A and the
*     reciprocal pivot growth factor RPVGRW.
*
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = DLANGE( NORM, N, N, A, LDA, WORK )
      RPVGRW = DLANTR( 'M', 'U', 'N', N, N, AF, LDAF, WORK )
      IF( RPVGRW.EQ.ZERO ) THEN
         RPVGRW = ONE
      ELSE
         RPVGRW = DLANGE( 'M', N, N, A, LDA, WORK ) / RPVGRW
      END IF
*
*     Compute the reciprocal of the condition number of A.
*
      CALL DGECON( NORM, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO )
*
*     Compute the solution matrix X.
*
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DGETRS( TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X,
     $             LDX, FERR, BERR, WORK, IWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      IF( NOTRAN ) THEN
         IF( COLEQU ) THEN
            DO 80 J = 1, NRHS
               DO 70 I = 1, N
                  X( I, J ) = C( I )*X( I, J )
   70          CONTINUE
   80       CONTINUE
            DO 90 J = 1, NRHS
               FERR( J ) = FERR( J ) / COLCND
   90       CONTINUE
         END IF
      ELSE IF( ROWEQU ) THEN
         DO 110 J = 1, NRHS
            DO 100 I = 1, N
               X( I, J ) = R( I )*X( I, J )
  100       CONTINUE
  110    CONTINUE
         DO 120 J = 1, NRHS
            FERR( J ) = FERR( J ) / ROWCND
  120    CONTINUE
      END IF
*
      WORK( 1 ) = RPVGRW
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) )
     $   INFO = N + 1
      RETURN
*
*     End of DGESVX
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN 
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH      
      INTEGER            IDAMAX
      EXTERNAL           DLAMCH, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum 
* 
      SFMIN = DLAMCH('S')  
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN 
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN 
                  CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 ) 
               ELSE 
                 DO 20 I = 1, M-J 
                    A( J+I, J ) = A( J+I, J ) / A( J, J ) 
   20            CONTINUE 
               END IF 
            END IF 
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dgetv0
c
c\Description: 
c  Generate a random initial residual vector for the Arnoldi process.
c  Force the residual vector to be in the range of the operator OP.  
c
c\Usage:
c  call dgetv0
c     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
c       IPNTR, WORKD, IERR )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to dgetv0.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B in the (generalized)
c          eigenvalue problem A*x = lambda*B*x.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  ITRY    Integer.  (INPUT)
c          ITRY counts the number of times that dgetv0 is called.  
c          It should be set to 1 on the initial call to dgetv0.
c
c  INITV   Logical variable.  (INPUT)
c          .TRUE.  => the initial residual vector is given in RESID.
c          .FALSE. => generate a random initial residual vector.
c
c  N       Integer.  (INPUT)
c          Dimension of the problem.
c
c  J       Integer.  (INPUT)
c          Index of the residual vector to be generated, with respect to
c          the Arnoldi process.  J > 1 in case of a "restart".
c
c  V       Double precision N by J array.  (INPUT)
c          The first J-1 columns of V contain the current Arnoldi basis
c          if this is a "restart".
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          Initial residual vector to be generated.  If RESID is 
c          provided, force RESID into the range of the operator OP.
c
c  RNORM   Double precision scalar.  (OUTPUT)
c          B-norm of the generated residual.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c
c  WORKD   Double precision work array of length 2*N.  (REVERSE COMMUNICATION).
c          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
c
c  IERR    Integer.  (OUTPUT)
c          =  0: Normal exit.
c          = -1: Cannot generate a nontrivial restarted residual vector
c                in the range of the operator OP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine for vector output.
c     dlarnv  LAPACK routine for generating a random vector.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\SCCS Information: @(#) 
c FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dgetv0 
     &   ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, 
     &     ipntr, workd, ierr )
c 
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           resid(n), v(ldv,j), workd(2*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    first, inits, orth
      integer    idist, iseed(4), iter, msglvl, jj
      Double precision
     &           rnorm0
      save       first, iseed, inits, iter, msglvl, orth, rnorm0
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dlarnv, dvout, dcopy, dgemv, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2
      external   ddot, dnrm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %-----------------%
c     | Data Statements |
c     %-----------------%
c
      data       inits /.true./
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-----------------------------------%
c     | Initialize the seed of the LAPACK |
c     | random number generator           |
c     %-----------------------------------%
c
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
c
      if (ido .eq.  0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call second (t0)
         msglvl = mgetv0
c 
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
c
c        %-----------------------------------------------------%
c        | Possibly generate a random starting vector in RESID |
c        | Use a LAPACK random number generator used by the    |
c        | matrix generation routines.                         |
c        |    idist = 1: uniform (0,1)  distribution;          |
c        |    idist = 2: uniform (-1,1) distribution;          |
c        |    idist = 3: normal  (0,1)  distribution;          |
c        %-----------------------------------------------------%
c
         if (.not.initv) then
            idist = 2
            call dlarnv (idist, iseed, n, resid)
         end if
c 
c        %----------------------------------------------------------%
c        | Force the starting vector into the range of OP to handle |
c        | the generalized problem when B is possibly (singular).   |
c        %----------------------------------------------------------%
c
         call second (t2)
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call dcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
c 
c     %-----------------------------------------%
c     | Back from computing OP*(initial-vector) |
c     %-----------------------------------------%
c
      if (first) go to 20
c
c     %-----------------------------------------------%
c     | Back from computing B*(orthogonalized-vector) |
c     %-----------------------------------------------%
c
      if (orth)  go to 40
c 
      if (bmat .eq. 'G') then
         call second (t3)
         tmvopx = tmvopx + (t3 - t2)
      end if
c 
c     %------------------------------------------------------%
c     | Starting vector is now in the range of OP; r = OP*r; |
c     | Compute B-norm of starting vector.                   |
c     %------------------------------------------------------%
c
      call second (t2)
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
c 
   20 continue
c
      if (bmat .eq. 'G') then
         call second (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
c 
      first = .FALSE.
      if (bmat .eq. 'G') then
          rnorm0 = ddot (n, resid, 1, workd, 1)
          rnorm0 = sqrt(abs(rnorm0))
      else if (bmat .eq. 'I') then
           rnorm0 = dnrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
c
c     %---------------------------------------------%
c     | Exit if this is the very first Arnoldi step |
c     %---------------------------------------------%
c
      if (j .eq. 1) go to 50
c 
c     %----------------------------------------------------------------
c     | Otherwise need to B-orthogonalize the starting vector against |
c     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
c     | This is the case where an invariant subspace is encountered   |
c     | in the middle of the Arnoldi factorization.                   |
c     |                                                               |
c     |       s = V^{T}*B*r;   r = r - V*s;                           |
c     |                                                               |
c     | Stopping criteria used for iter. ref. is discussed in         |
c     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
c     %---------------------------------------------------------------%
c
      orth = .TRUE.
   30 continue
c
      call dgemv ('T', n, j-1, one, v, ldv, workd, 1, 
     &            zero, workd(n+1), 1)
      call dgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, 
     &            one, resid, 1)
c 
c     %----------------------------------------------------------%
c     | Compute the B-norm of the orthogonalized starting vector |
c     %----------------------------------------------------------%
c
      call second (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call dcopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat .eq. 'I') then
         call dcopy (n, resid, 1, workd, 1)
      end if
c 
   40 continue
c
      if (bmat .eq. 'G') then
         call second (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
c 
      if (bmat .eq. 'G') then
         rnorm = ddot (n, resid, 1, workd, 1)
         rnorm = sqrt(abs(rnorm))
      else if (bmat .eq. 'I') then
         rnorm = dnrm2(n, resid, 1)
      end if
c
c     %--------------------------------------%
c     | Check for further orthogonalization. |
c     %--------------------------------------%
c
      if (msglvl .gt. 2) then
          call dvout (logfil, 1, rnorm0, ndigit, 
     &                '_getv0: re-orthonalization ; rnorm0 is')
          call dvout (logfil, 1, rnorm, ndigit, 
     &                '_getv0: re-orthonalization ; rnorm is')
      end if
c
      if (rnorm .gt. 0.717*rnorm0) go to 50
c 
      iter = iter + 1
      if (iter .le. 5) then
c
c        %-----------------------------------%
c        | Perform iterative refinement step |
c        %-----------------------------------%
c
         rnorm0 = rnorm
         go to 30
      else
c
c        %------------------------------------%
c        | Iterative refinement step "failed" |
c        %------------------------------------%
c
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = zero
         ierr = -1
      end if
c 
   50 continue
c
      if (msglvl .gt. 0) then
         call dvout (logfil, 1, rnorm, ndigit,
     &        '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl .gt. 3) then
         call dvout (logfil, n, resid, ndigit,
     &        '_getv0: initial / restarted starting vector')
      end if
      ido = 99
c 
      call second (t1)
      tgetv0 = tgetv0 + (t1 - t0)
c 
 9000 continue
      return
c
c     %---------------%
c     | End of dgetv0 |
c     %---------------%
c
      end
      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,
     $                   LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   Z( LDZ, * )
*     ..
*     Purpose
*     =======
*
*     DHSEQR computes the eigenvalues of a Hessenberg matrix H
*     and, optionally, the matrices T and Z from the Schur decomposition
*     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
*     Schur form), and Z is the orthogonal matrix of Schur vectors.
*
*     Optionally Z may be postmultiplied into an input orthogonal
*     matrix Q so that this routine can give the Schur factorization
*     of a matrix A which has been reduced to the Hessenberg form H
*     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
*
*     Arguments
*     =========
*
*     JOB   (input) CHARACTER*1
*           = 'E':  compute eigenvalues only;
*           = 'S':  compute eigenvalues and the Schur form T.
*
*     COMPZ (input) CHARACTER*1
*           = 'N':  no Schur vectors are computed;
*           = 'I':  Z is initialized to the unit matrix and the matrix Z
*                   of Schur vectors of H is returned;
*           = 'V':  Z must contain an orthogonal matrix Q on entry, and
*                   the product Q*Z is returned.
*
*     N     (input) INTEGER
*           The order of the matrix H.  N .GE. 0.
*
*     ILO   (input) INTEGER
*     IHI   (input) INTEGER
*           It is assumed that H is already upper triangular in rows
*           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*           set by a previous call to DGEBAL, and then passed to DGEHRD
*           when the matrix output by DGEBAL is reduced to Hessenberg
*           form. Otherwise ILO and IHI should be set to 1 and N
*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*           If N = 0, then ILO = 1 and IHI = 0.
*
*     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*           On entry, the upper Hessenberg matrix H.
*           On exit, if INFO = 0 and JOB = 'S', then H contains the
*           upper quasi-triangular matrix T from the Schur decomposition
*           (the Schur form); 2-by-2 diagonal blocks (corresponding to
*           complex conjugate pairs of eigenvalues) are returned in
*           standard form, with H(i,i) = H(i+1,i+1) and
*           H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and JOB = 'E', the
*           contents of H are unspecified on exit.  (The output value of
*           H when INFO.GT.0 is given under the description of INFO
*           below.)
*
*           Unlike earlier versions of DHSEQR, this subroutine may
*           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
*           or j = IHI+1, IHI+2, ... N.
*
*     LDH   (input) INTEGER
*           The leading dimension of the array H. LDH .GE. max(1,N).
*
*     WR    (output) DOUBLE PRECISION array, dimension (N)
*     WI    (output) DOUBLE PRECISION array, dimension (N)
*           The real and imaginary parts, respectively, of the computed
*           eigenvalues. If two eigenvalues are computed as a complex
*           conjugate pair, they are stored in consecutive elements of
*           WR and WI, say the i-th and (i+1)th, with WI(i) .GT. 0 and
*           WI(i+1) .LT. 0. If JOB = 'S', the eigenvalues are stored in
*           the same order as on the diagonal of the Schur form returned
*           in H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2
*           diagonal block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
*           WI(i+1) = -WI(i).
*
*     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*           If COMPZ = 'N', Z is not referenced.
*           If COMPZ = 'I', on entry Z need not be set and on exit,
*           if INFO = 0, Z contains the orthogonal matrix Z of the Schur
*           vectors of H.  If COMPZ = 'V', on entry Z must contain an
*           N-by-N matrix Q, which is assumed to be equal to the unit
*           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
*           if INFO = 0, Z contains Q*Z.
*           Normally Q is the orthogonal matrix generated by DORGHR
*           after the call to DGEHRD which formed the Hessenberg matrix
*           H. (The output value of Z when INFO.GT.0 is given under
*           the description of INFO below.)
*
*     LDZ   (input) INTEGER
*           The leading dimension of the array Z.  if COMPZ = 'I' or
*           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
*
*     WORK  (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*           On exit, if INFO = 0, WORK(1) returns an estimate of
*           the optimal value for LWORK.
*
*     LWORK (input) INTEGER
*           The dimension of the array WORK.  LWORK .GE. max(1,N)
*           is sufficient and delivers very good and sometimes
*           optimal performance.  However, LWORK as large as 11*N
*           may be required for optimal performance.  A workspace
*           query is recommended to determine the optimal workspace
*           size.
*
*           If LWORK = -1, then DHSEQR does a workspace query.
*           In this case, DHSEQR checks the input parameters and
*           estimates the optimal workspace size for the given
*           values of N, ILO and IHI.  The estimate is returned
*           in WORK(1).  No error message related to LWORK is
*           issued by XERBLA.  Neither H nor Z are accessed.
*
*
*     INFO  (output) INTEGER
*             =  0:  successful exit
*           .LT. 0:  if INFO = -i, the i-th argument had an illegal
*                    value
*           .GT. 0:  if INFO = i, DHSEQR failed to compute all of
*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*                and WI contain those eigenvalues which have been
*                successfully computed.  (Failures are rare.)
*
*                If INFO .GT. 0 and JOB = 'E', then on exit, the
*                remaining unconverged eigenvalues are the eigen-
*                values of the upper Hessenberg matrix rows and
*                columns ILO through INFO of the final, output
*                value of H.
*
*                If INFO .GT. 0 and JOB   = 'S', then on exit
*
*           (*)  (initial value of H)*U  = U*(final value of H)
*
*                where U is an orthogonal matrix.  The final
*                value of H is upper Hessenberg and quasi-triangular
*                in rows and columns INFO+1 through IHI.
*
*                If INFO .GT. 0 and COMPZ = 'V', then on exit
*
*                  (final value of Z)  =  (initial value of Z)*U
*
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'I', then on exit
*                      (final value of Z)  = U
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of JOB.)
*
*                If INFO .GT. 0 and COMPZ = 'N', then Z is not
*                accessed.
*
*     ================================================================
*             Default values supplied by
*             ILAENV(ISPEC,'DHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
*             It is suggested that these defaults be adjusted in order
*             to attain best performance in each particular
*             computational environment.
*
*            ISPEC=12: The DLAHQR vs DLAQR0 crossover point.
*                      Default: 75. (Must be at least 11.)
*
*            ISPEC=13: Recommended deflation window size.
*                      This depends on ILO, IHI and NS.  NS is the
*                      number of simultaneous shifts returned
*                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
*                      The default for (IHI-ILO+1).LE.500 is NS.
*                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*            ISPEC=14: Nibble crossover point. (See IPARMQ for
*                      details.)  Default: 14% of deflation window
*                      size.
*
*            ISPEC=15: Number of simultaneous shifts in a multishift
*                      QR iteration.
*
*                      If IHI-ILO+1 is ...
*
*                      greater than      ...but less    ... the
*                      or equal to ...      than        default is
*
*                           1               30          NS =   2(+)
*                          30               60          NS =   4(+)
*                          60              150          NS =  10(+)
*                         150              590          NS =  **
*                         590             3000          NS =  64
*                        3000             6000          NS = 128
*                        6000             infinity      NS = 256
*
*                  (+)  By default some or all matrices of this order
*                       are passed to the implicit double shift routine
*                       DLAHQR and this parameter is ignored.  See
*                       ISPEC=12 above and comments in IPARMQ for
*                       details.
*
*                 (**)  The asterisks (**) indicate an ad-hoc
*                       function of N increasing from 10 to 64.
*
*            ISPEC=16: Select structured matrix multiply.
*                      If the number of simultaneous shifts (specified
*                      by ISPEC=15) is less than 14, then the default
*                      for ISPEC=16 is 0.  Otherwise the default for
*                      ISPEC=16 is 2.
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*     References:
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*       929--947, 2002.
*
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*       of Matrix Analysis, volume 23, pages 948--973, 2002.
*
*     ================================================================
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
*
*     ==== NL allocates some local workspace to help small matrices
*     .    through a rare DLAHQR failure.  NL .GT. NTINY = 11 is
*     .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
*     .    allows up to six simultaneous shifts and a 16-by-16
*     .    deflation window.  ====
      INTEGER            NL
      PARAMETER          ( NL = 49 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   HL( NL, NL ), WORKL( NL )
*     ..
*     .. Local Scalars ..
      INTEGER            I, KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      LOGICAL            LSAME
      EXTERNAL           ILAENV, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLAHQR, DLAQR0, DLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     ==== Decode and check the input parameters. ====
*
      WANTT = LSAME( JOB, 'S' )
      INITZ = LSAME( COMPZ, 'I' )
      WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
      WORK( 1 ) = DBLE( MAX( 1, N ) )
      LQUERY = LWORK.EQ.-1
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.NE.0 ) THEN
*
*        ==== Quick return in case of invalid argument. ====
*
         CALL XERBLA( 'DHSEQR', -INFO )
         RETURN
*
      ELSE IF( N.EQ.0 ) THEN
*
*        ==== Quick return in case N = 0; nothing to do. ====
*
         RETURN
*
      ELSE IF( LQUERY ) THEN
*
*        ==== Quick return in case of a workspace query ====
*
         CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO,
     $                IHI, Z, LDZ, WORK, LWORK, INFO )
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
         RETURN
*
      ELSE
*
*        ==== copy eigenvalues isolated by DGEBAL ====
*
         DO 10 I = 1, ILO - 1
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   10    CONTINUE
         DO 20 I = IHI + 1, N
            WR( I ) = H( I, I )
            WI( I ) = ZERO
   20    CONTINUE
*
*        ==== Initialize Z, if requested ====
*
         IF( INITZ )
     $      CALL DLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
*
*        ==== Quick return if possible ====
*
         IF( ILO.EQ.IHI ) THEN
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         NMIN = ILAENV( 12, 'DHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N,
     $          ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== DLAQR0 for big matrices; DLAHQR for small ones ====
*
         IF( N.GT.NMIN ) THEN
            CALL DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO,
     $                   IHI, Z, LDZ, WORK, LWORK, INFO )
         ELSE
*
*           ==== Small matrix ====
*
            CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO,
     $                   IHI, Z, LDZ, INFO )
*
            IF( INFO.GT.0 ) THEN
*
*              ==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
*              .    when DLAHQR fails. ====
*
               KBOT = INFO
*
               IF( N.GE.NL ) THEN
*
*                 ==== Larger matrices have enough subdiagonal scratch
*                 .    space to call DLAQR0 directly. ====
*
                  CALL DLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR,
     $                         WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
*
               ELSE
*
*                 ==== Tiny matrices don't have enough subdiagonal
*                 .    scratch space to benefit from DLAQR0.  Hence,
*                 .    tiny matrices must be copied into a larger
*                 .    array before calling DLAQR0. ====
*
                  CALL DLACPY( 'A', N, N, H, LDH, HL, NL )
                  HL( N+1, N ) = ZERO
                  CALL DLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ),
     $                         NL )
                  CALL DLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR,
     $                         WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )
                  IF( WANTT .OR. INFO.NE.0 )
     $               CALL DLACPY( 'A', N, N, HL, NL, H, LDH )
               END IF
            END IF
         END IF
*
*        ==== Clear out the trash, if necessary. ====
*
         IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 )
     $      CALL DLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
*
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
*
         WORK( 1 ) = MAX( DBLE( MAX( 1, N ) ), WORK( 1 ) )
      END IF
*
*     ==== End of DHSEQR ====
*
      END
      LOGICAL FUNCTION DISNAN(DIN)
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN
*     ..
*
*  Purpose
*  =======
*
*  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  DIN      (input) DOUBLE PRECISION
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
*  ..
*  .. Executable Statements ..
      DISNAN = DLAISNAN(DIN,DIN)
      RETURN
      END
      SUBROUTINE DLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  Purpose
*  =======
*
*  DLABAD takes as input the values computed by DLAMCH for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of LARGE is sufficiently large.  This subroutine is intended to
*  identify machines with a large exponent range, such as the Crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by DLAMCH.  This subroutine is needed because
*  DLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  Arguments
*  =========
*
*  SMALL   (input/output) DOUBLE PRECISION
*          On entry, the underflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (input/output) DOUBLE PRECISION
*          On entry, the overflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
*
      RETURN
*
*     End of DLABAD
*
      END
      SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..
      INTEGER            ISGN( * ), ISAVE( 3 )
      DOUBLE PRECISION   V( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLACN2 estimates the 1-norm of a square, real matrix A.
*  Reverse communication is used for evaluating matrix-vector products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  V      (workspace) DOUBLE PRECISION array, dimension (N)
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and DLACN2 must be re-called with all the other parameters
*         unchanged.
*
*  ISGN   (workspace) INTEGER array, dimension (N)
*
*  EST    (input/output) DOUBLE PRECISION
*         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
*         unchanged from the previous call to DLACN2.
*         On exit, EST is an estimate (a lower bound) for norm(A). 
*
*  KASE   (input/output) INTEGER
*         On the initial call to DLACN2, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from DLACN2, KASE will again be 0.
*
*  ISAVE  (input/output) INTEGER array, dimension (3)
*         ISAVE is used to save variables between calls to DLACN2
*
*  Further Details
*  ======= =======
*
*  Contributed by Nick Higham, University of Manchester.
*  Originally named SONEST, dated March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*  This is a thread safe version of DLACON, which uses the array ISAVE
*  in place of a SAVE statement, as follows:
*
*     DLACON     DLACN2
*      JUMP     ISAVE(1)
*      J        ISAVE(2)
*      ITER     ISAVE(3)
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, JLAST
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM
      EXTERNAL           IDAMAX, DASUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
*     ..
*     .. Executable Statements ..
*
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      END IF
*
      GO TO ( 20, 40, 70, 110, 140 )ISAVE( 1 )
*
*     ................ ENTRY   (ISAVE( 1 ) = 1)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
*
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
*        ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
*
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN
*
*     ................ ENTRY   (ISAVE( 1 ) = 2)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*
   40 CONTINUE
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      ISAVE( 3 ) = 2
*
*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
*
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( ISAVE( 2 ) ) = ONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN
*
*     ................ ENTRY   (ISAVE( 1 ) = 3)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) )
     $      GO TO 90
   80 CONTINUE
*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
*
   90 CONTINUE
*     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )
     $   GO TO 120
*
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN
*
*     ................ ENTRY   (ISAVE( 1 ) = 4)
*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*
  110 CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND.
     $    ( ISAVE( 3 ).LT.ITMAX ) ) THEN
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      END IF
*
*     ITERATION COMPLETE.  FINAL STAGE.
*
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN
*
*     ................ ENTRY   (ISAVE( 1 ) = 5)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
*
  150 CONTINUE
      KASE = 0
      RETURN
*
*     End of DLACN2
*
      END
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of DLACPY
*
      END
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   E, F
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
*
      RETURN
*
*     End of DLADIV
*
      END
      SUBROUTINE DLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK,
     $                   INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
*  an upper quasi-triangular matrix T by an orthogonal similarity
*  transformation.
*
*  T must be in Schur canonical form, that is, block upper triangular
*  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
*  has its diagonal elemnts equal and its off-diagonal elements of
*  opposite sign.
*
*  Arguments
*  =========
*
*  WANTQ   (input) LOGICAL
*          = .TRUE. : accumulate the transformation in the matrix Q;
*          = .FALSE.: do not accumulate the transformation.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*          On exit, the updated matrix T, again in Schur canonical form.
*
*  LDT     (input)  INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
*          On exit, if WANTQ is .TRUE., the updated matrix Q.
*          If WANTQ is .FALSE., Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
*
*  J1      (input) INTEGER
*          The index of the first row of the first block T11.
*
*  N1      (input) INTEGER
*          The order of the first block T11. N1 = 0, 1 or 2.
*
*  N2      (input) INTEGER
*          The order of the second block T22. N2 = 0, 1 or 2.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          = 1: the transformed matrix T would be too far from Schur
*               form; the blocks are not swapped and T and Q are
*               unchanged.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 1.0D+1 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      DOUBLE PRECISION   CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22,
     $                   T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2,
     $                   WR1, WR2, XNORM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ),
     $                   X( LDX, 2 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLANV2, DLARFG, DLARFX, DLARTG, DLASY2,
     $                   DROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. N1.EQ.0 .OR. N2.EQ.0 )
     $   RETURN
      IF( J1+N1.GT.N )
     $   RETURN
*
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
*
      IF( N1.EQ.1 .AND. N2.EQ.1 ) THEN
*
*        Swap two 1-by-1 blocks.
*
         T11 = T( J1, J1 )
         T22 = T( J2, J2 )
*
*        Determine the transformation to perform the interchange.
*
         CALL DLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
*
*        Apply transformation to the matrix T.
*
         IF( J3.LE.N )
     $      CALL DROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS,
     $                 SN )
         CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
*
         T( J1, J1 ) = T22
         T( J2, J2 ) = T11
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
*
      ELSE
*
*        Swapping involves at least one 2-by-2 block.
*
*        Copy the diagonal block of order N1+N2 to the local array D
*        and compute its norm.
*
         ND = N1 + N2
         CALL DLACPY( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
         DNORM = DLANGE( 'Max', ND, ND, D, LDD, WORK )
*
*        Compute machine-dependent threshold for test for accepting
*        swap.
*
         EPS = DLAMCH( 'P' )
         SMLNUM = DLAMCH( 'S' ) / EPS
         THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
*
*        Solve T11*X - X*T22 = scale*T12 for X.
*
         CALL DLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD,
     $                D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X,
     $                LDX, XNORM, IERR )
*
*        Swap the adjacent diagonal blocks.
*
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
*
   10    CONTINUE
*
*        N1 = 1, N2 = 2: generate elementary reflector H so that:
*
*        ( scale, X11, X12 ) H = ( 0, 0, * )
*
         U( 1 ) = SCALE
         U( 2 ) = X( 1, 1 )
         U( 3 ) = X( 1, 2 )
         CALL DLARFG( 3, U( 3 ), U, 1, TAU )
         U( 3 ) = ONE
         T11 = T( J1, J1 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3,
     $       3 )-T11 ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL DLARFX( 'L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK )
*
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J3, J3 ) = T11
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
*
   20    CONTINUE
*
*        N1 = 2, N2 = 1: generate elementary reflector H so that:
*
*        H (  -X11 ) = ( * )
*          (  -X21 ) = ( 0 )
*          ( scale ) = ( 0 )
*
         U( 1 ) = -X( 1, 1 )
         U( 2 ) = -X( 2, 1 )
         U( 3 ) = SCALE
         CALL DLARFG( 3, U( 1 ), U( 2 ), 1, TAU )
         U( 1 ) = ONE
         T33 = T( J3, J3 )
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL DLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL DLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1,
     $       1 )-T33 ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL DLARFX( 'R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK )
*
         T( J1, J1 ) = T33
         T( J2, J1 ) = ZERO
         T( J3, J1 ) = ZERO
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
*
   30    CONTINUE
*
*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
*        that:
*
*        H(2) H(1) (  -X11  -X12 ) = (  *  * )
*                  (  -X21  -X22 )   (  0  * )
*                  ( scale    0  )   (  0  0 )
*                  (    0  scale )   (  0  0 )
*
         U1( 1 ) = -X( 1, 1 )
         U1( 2 ) = -X( 2, 1 )
         U1( 3 ) = SCALE
         CALL DLARFG( 3, U1( 1 ), U1( 2 ), 1, TAU1 )
         U1( 1 ) = ONE
*
         TEMP = -TAU1*( X( 1, 2 )+U1( 2 )*X( 2, 2 ) )
         U2( 1 ) = -TEMP*U1( 2 ) - X( 2, 2 )
         U2( 2 ) = -TEMP*U1( 3 )
         U2( 3 ) = SCALE
         CALL DLARFG( 3, U2( 1 ), U2( 2 ), 1, TAU2 )
         U2( 1 ) = ONE
*
*        Perform swap provisionally on diagonal block in D.
*
         CALL DLARFX( 'L', 3, 4, U1, TAU1, D, LDD, WORK )
         CALL DLARFX( 'R', 4, 3, U1, TAU1, D, LDD, WORK )
         CALL DLARFX( 'L', 3, 4, U2, TAU2, D( 2, 1 ), LDD, WORK )
         CALL DLARFX( 'R', 4, 3, U2, TAU2, D( 1, 2 ), LDD, WORK )
*
*        Test whether to reject swap.
*
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ),
     $       ABS( D( 4, 2 ) ) ).GT.THRESH )GO TO 50
*
*        Accept swap: apply transformation to the entire matrix T.
*
         CALL DLARFX( 'L', 3, N-J1+1, U1, TAU1, T( J1, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J4, 3, U1, TAU1, T( 1, J1 ), LDT, WORK )
         CALL DLARFX( 'L', 3, N-J1+1, U2, TAU2, T( J2, J1 ), LDT, WORK )
         CALL DLARFX( 'R', J4, 3, U2, TAU2, T( 1, J2 ), LDT, WORK )
*
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J4, J1 ) = ZERO
         T( J4, J2 ) = ZERO
*
         IF( WANTQ ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL DLARFX( 'R', N, 3, U1, TAU1, Q( 1, J1 ), LDQ, WORK )
            CALL DLARFX( 'R', N, 3, U2, TAU2, Q( 1, J2 ), LDQ, WORK )
         END IF
*
   40    CONTINUE
*
         IF( N2.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T11
*
            CALL DLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ),
     $                   T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
            CALL DROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT,
     $                 CS, SN )
            CALL DROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
            IF( WANTQ )
     $         CALL DROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
*
         IF( N1.EQ.2 ) THEN
*
*           Standardize new 2-by-2 block T22
*
            J3 = J1 + N2
            J4 = J3 + 1
            CALL DLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ),
     $                   T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
            IF( J3+2.LE.N )
     $         CALL DROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ),
     $                    LDT, CS, SN )
            CALL DROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
            IF( WANTQ )
     $         CALL DROT( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN )
         END IF
*
      END IF
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   50 CONTINUE
      INFO = 1
      RETURN
*
*     End of DLAEXC
*
      END
      SUBROUTINE DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
*     ..
*
*     Purpose
*     =======
*
*     DLAHQR is an auxiliary routine called by DHSEQR to update the
*     eigenvalues and Schur decomposition already computed by DHSEQR, by
*     dealing with the Hessenberg submatrix in rows and columns ILO to
*     IHI.
*
*     Arguments
*     =========
*
*     WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*     WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*     N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*     ILO     (input) INTEGER
*     IHI     (input) INTEGER
*          It is assumed that H is already upper quasi-triangular in
*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
*          ILO = 1). DLAHQR works primarily with the Hessenberg
*          submatrix in rows and columns ILO to IHI, but applies
*          transformations to all of H if WANTT is .TRUE..
*          1 <= ILO <= max(1,IHI); IHI <= N.
*
*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*          On entry, the upper Hessenberg matrix H.
*          On exit, if INFO is zero and if WANTT is .TRUE., H is upper
*          quasi-triangular in rows and columns ILO:IHI, with any
*          2-by-2 diagonal blocks in standard form. If INFO is zero
*          and WANTT is .FALSE., the contents of H are unspecified on
*          exit.  The output state of H if INFO is nonzero is given
*          below under the description of INFO.
*
*     LDH     (input) INTEGER
*          The leading dimension of the array H. LDH >= max(1,N).
*
*     WR      (output) DOUBLE PRECISION array, dimension (N)
*     WI      (output) DOUBLE PRECISION array, dimension (N)
*          The real and imaginary parts, respectively, of the computed
*          eigenvalues ILO to IHI are stored in the corresponding
*          elements of WR and WI. If two eigenvalues are computed as a
*          complex conjugate pair, they are stored in consecutive
*          elements of WR and WI, say the i-th and (i+1)th, with
*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
*          eigenvalues are stored in the same order as on the diagonal
*          of the Schur form returned in H, with WR(i) = H(i,i), and, if
*          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
*
*     ILOZ    (input) INTEGER
*     IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE..
*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*
*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          If WANTZ is .TRUE., on entry Z must contain the current
*          matrix Z of transformations accumulated by DHSEQR, and on
*          exit Z has been updated; transformations are applied only to
*          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*          If WANTZ is .FALSE., Z is not referenced.
*
*     LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= max(1,N).
*
*     INFO    (output) INTEGER
*           =   0: successful exit
*          .GT. 0: If INFO = i, DLAHQR failed to compute all the
*                  eigenvalues ILO to IHI in a total of 30 iterations
*                  per eigenvalue; elements i+1:ihi of WR and WI
*                  contain those eigenvalues which have been
*                  successfully computed.
*
*                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,
*                  the remaining unconverged eigenvalues are the
*                  eigenvalues of the upper Hessenberg matrix rows
*                  and columns ILO thorugh INFO of the final, output
*                  value of H.
*
*                  If INFO .GT. 0 and WANTT is .TRUE., then on exit
*          (*)       (initial value of H)*U  = U*(final value of H)
*                  where U is an orthognal matrix.    The final
*                  value of H is upper Hessenberg and triangular in
*                  rows and columns INFO+1 through IHI.
*
*                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit
*                      (final value of Z)  = (initial value of Z)*U
*                  where U is the orthogonal matrix in (*)
*                  (regardless of the value of WANTT.)
*
*     Further Details
*     ===============
*
*     02-96 Based on modifications by
*     David Day, Sandia National Laboratory, USA
*
*     12-04 Further modifications by
*     Ralph Byers, University of Kansas, USA
*     This is a modified version of DLAHQR from LAPACK version 3.0.
*     It is (1) more robust against overflow and underflow and
*     (2) adopts the more conservative Ahues & Tisseur stopping
*     criterion (LAWN 122, 1997).
*
*     =========================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0 )
      DOUBLE PRECISION   DAT1, DAT2
      PARAMETER          ( DAT1 = 3.0d0 / 4.0d0, DAT2 = -0.4375d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S,
     $                   H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX,
     $                   SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST,
     $                   ULP, V2, V3
      INTEGER            I, I1, I2, ITS, J, K, L, M, NH, NR, NZ
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   V( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLABAD, DLANV2, DLARFG, DROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( ILO.EQ.IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 J = ILO, IHI - 3
         H( J+2, J ) = ZERO
         H( J+3, J ) = ZERO
   10 CONTINUE
      IF( ILO.LE.IHI-2 )
     $   H( IHI, IHI-2 ) = ZERO
*
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( NH ) / ULP )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
*     H(L,L-1) is negligible so that the matrix splits.
*
      I = IHI
   20 CONTINUE
      L = ILO
      IF( I.LT.ILO )
     $   GO TO 160
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 or 2 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      DO 140 ITS = 0, ITMAX
*
*        Look for a single small subdiagonal element.
*
         DO 30 K = I, L + 1, -1
            IF( ABS( H( K, K-1 ) ).LE.SMLNUM )
     $         GO TO 40
            TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST.EQ.ZERO ) THEN
               IF( K-2.GE.ILO )
     $            TST = TST + ABS( H( K-1, K-2 ) )
               IF( K+1.LE.IHI )
     $            TST = TST + ABS( H( K+1, K ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some cases.  ====
            IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
               AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
               AA = MAX( ABS( H( K, K ) ),
     $              ABS( H( K-1, K-1 )-H( K, K ) ) )
               BB = MIN( ABS( H( K, K ) ),
     $              ABS( H( K-1, K-1 )-H( K, K ) ) )
               S = AA + AB
               IF( BA*( AB / S ).LE.MAX( SMLNUM,
     $             ULP*( BB*( AA / S ) ) ) )GO TO 40
            END IF
   30    CONTINUE
   40    CONTINUE
         L = K
         IF( L.GT.ILO ) THEN
*
*           H(L,L-1) is negligible
*
            H( L, L-1 ) = ZERO
         END IF
*
*        Exit from loop if a submatrix of order 1 or 2 has split off.
*
         IF( L.GE.I-1 )
     $      GO TO 150
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
*
         IF( ITS.EQ.10 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( L+1, L ) ) + ABS( H( L+2, L+1 ) )
            H11 = DAT1*S + H( L, L )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE IF( ITS.EQ.20 ) THEN
*
*           Exceptional shift.
*
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H11 = DAT1*S + H( I, I )
            H12 = DAT2*S
            H21 = S
            H22 = H11
         ELSE
*
*           Prepare to use Francis' double shift
*           (i.e. 2nd degree generalized Rayleigh quotient)
*
            H11 = H( I-1, I-1 )
            H21 = H( I, I-1 )
            H12 = H( I-1, I )
            H22 = H( I, I )
         END IF
         S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
         IF( S.EQ.ZERO ) THEN
            RT1R = ZERO
            RT1I = ZERO
            RT2R = ZERO
            RT2I = ZERO
         ELSE
            H11 = H11 / S
            H21 = H21 / S
            H12 = H12 / S
            H22 = H22 / S
            TR = ( H11+H22 ) / TWO
            DET = ( H11-TR )*( H22-TR ) - H12*H21
            RTDISC = SQRT( ABS( DET ) )
            IF( DET.GE.ZERO ) THEN
*
*              ==== complex conjugate shifts ====
*
               RT1R = TR*S
               RT2R = RT1R
               RT1I = RTDISC*S
               RT2I = -RT1I
            ELSE
*
*              ==== real shifts (use only one of them)  ====
*
               RT1R = TR + RTDISC
               RT2R = TR - RTDISC
               IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
                  RT1R = RT1R*S
                  RT2R = RT1R
               ELSE
                  RT2R = RT2R*S
                  RT1R = RT2R
               END IF
               RT1I = ZERO
               RT2I = ZERO
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 50 M = I - 2, L, -1
*           Determine the effect of starting the double-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.  (The following uses scaling to avoid
*           overflows and most underflows.)
*
            H21S = H( M+1, M )
            S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
            H21S = H( M+1, M ) / S
            V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )*
     $               ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
            V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
            V( 3 ) = H21S*H( M+2, M+1 )
            S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
            V( 1 ) = V( 1 ) / S
            V( 2 ) = V( 2 ) / S
            V( 3 ) = V( 3 ) / S
            IF( M.EQ.L )
     $         GO TO 60
            IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE.
     $          ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M,
     $          M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
   50    CONTINUE
   60    CONTINUE
*
*        Double-shift QR step
*
         DO 130 K = M, I - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix. NR is the order of G.
*
            NR = MIN( 3, I-K+1 )
            IF( K.GT.M )
     $         CALL DCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL DLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K.GT.M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K.LT.I-1 )
     $            H( K+2, K-1 ) = ZERO
            ELSE IF( M.GT.L ) THEN
*               ==== Use the following instead of
*               .    H( K, K-1 ) = -H( K, K-1 ) to
*               .    avoid a bug when v(2) and v(3)
*               .    underflow. ====
               H( K, K-1 ) = H( K, K-1 )*( ONE-T1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR.EQ.3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 70 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   70          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 80 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   80          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 90 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   90             CONTINUE
               END IF
            ELSE IF( NR.EQ.2 ) THEN
*
*              Apply G from the left to transform the rows of the matrix
*              in columns K to I2.
*
               DO 100 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
  100          CONTINUE
*
*              Apply G from the right to transform the columns of the
*              matrix in rows I1 to min(K+3,I).
*
               DO 110 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  110          CONTINUE
*
               IF( WANTZ ) THEN
*
*                 Accumulate transformations in the matrix Z
*
                  DO 120 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  120             CONTINUE
               END IF
            END IF
  130    CONTINUE
*
  140 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      INFO = I
      RETURN
*
  150 CONTINUE
*
      IF( L.EQ.I ) THEN
*
*        H(I,I-1) is negligible: one eigenvalue has converged.
*
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L.EQ.I-1 ) THEN
*
*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
*
*        Transform the 2-by-2 submatrix to standard Schur form,
*        and compute and store the eigenvalues.
*
         CALL DLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ),
     $                H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ),
     $                CS, SN )
*
         IF( WANTT ) THEN
*
*           Apply the transformation to the rest of H.
*
            IF( I2.GT.I )
     $         CALL DROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH,
     $                    CS, SN )
            CALL DROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
*
*           Apply the transformation to Z.
*
            CALL DROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
      END IF
*
*     return to start of the main loop with new value of I.
*
      I = L - 1
      GO TO 20
*
  160 CONTINUE
      RETURN
*
*     End of DLAHQR
*
      END
      SUBROUTINE DLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2009                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  A( LDA, * ), T( LDT, NB ), TAU( NB ),
     $                   Y( LDY, NB )
*     ..
*
*  Purpose
*  =======
*
*  DLAHR2 reduces the first NB columns of A real general n-BY-(n-k+1)
*  matrix A so that elements below the k-th subdiagonal are zero. The
*  reduction is performed by an orthogonal similarity transformation
*  Q' * A * Q. The routine returns the matrices V and T which determine
*  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T.
*
*  This is an auxiliary routine called by DGEHRD.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  K       (input) INTEGER
*          The offset for the reduction. Elements below the k-th
*          subdiagonal in the first NB columns are reduced to zero.
*          K < N.
*
*  NB      (input) INTEGER
*          The number of columns to be reduced.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N-K+1)
*          On entry, the n-by-(n-k+1) general matrix A.
*          On exit, the elements on and above the k-th subdiagonal in
*          the first NB columns are overwritten with the corresponding
*          elements of the reduced matrix; the elements below the k-th
*          subdiagonal, with the array TAU, represent the matrix Q as a
*          product of elementary reflectors. The other columns of A are
*          unchanged. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (NB)
*          The scalar factors of the elementary reflectors. See Further
*          Details.
*
*  T       (output) DOUBLE PRECISION array, dimension (LDT,NB)
*          The upper triangular matrix T.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T.  LDT >= NB.
*
*  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
*          The n-by-nb matrix Y.
*
*  LDY     (input) INTEGER
*          The leading dimension of the array Y. LDY >= N.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of nb elementary reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
*  A(i+k+1:n,i), and tau in TAU(i).
*
*  The elements of the vectors v together form the (n-k+1)-by-nb matrix
*  V which is needed, with T and Y, to apply the transformation to the
*  unreduced part of the matrix, using an update of the form:
*  A := (I - V*T*V') * (A - Y*V').
*
*  The contents of A on exit are illustrated by the following example
*  with n = 7, k = 3 and nb = 2:
*
*     ( a   a   a   a   a )
*     ( a   a   a   a   a )
*     ( a   a   a   a   a )
*     ( h   h   a   a   a )
*     ( v1  h   a   a   a )
*     ( v1  v2  a   a   a )
*     ( v1  v2  a   a   a )
*
*  where a denotes an element of the original matrix A, h denotes a
*  modified element of the upper Hessenberg matrix H, and vi denotes an
*  element of the vector defining H(i).
*
*  This subroutine is a slight modification of LAPACK-3.0's DLAHRD
*  incorporating improvements proposed by Quintana-Orti and Van de
*  Gejin. Note that the entries of A(1:K,2:NB) differ from those
*  returned by the original LAPACK-3.0's DLAHRD routine. (This
*  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)
*
*  References
*  ==========
*
*  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
*  performance of reduction to Hessenberg form," ACM Transactions on
*  Mathematical Software, 32(2):180-194, June 2006.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, 
     $                     ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION  EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMM, DGEMV, DLACPY,
     $                   DLARFG, DSCAL, DTRMM, DTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
*
*           Update A(K+1:N,I)
*
*           Update I-th column of A - Y * V'
*
            CALL DGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY,
     $                  A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
*
*           Apply I - V * T' * V' to this column (call it b) from the
*           left, using the last column of T as workspace
*
*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
*                    ( V2 )             ( b2 )
*
*           where V1 is unit lower triangular
*
*           w := V1' * b1
*
            CALL DCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL DTRMV( 'Lower', 'Transpose', 'UNIT', 
     $                  I-1, A( K+1, 1 ),
     $                  LDA, T( 1, NB ), 1 )
*
*           w := w + V2'*b2
*
            CALL DGEMV( 'Transpose', N-K-I+1, I-1, 
     $                  ONE, A( K+I, 1 ),
     $                  LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
*
*           w := T'*w
*
            CALL DTRMV( 'Upper', 'Transpose', 'NON-UNIT', 
     $                  I-1, T, LDT,
     $                  T( 1, NB ), 1 )
*
*           b2 := b2 - V2*w
*
            CALL DGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE, 
     $                  A( K+I, 1 ),
     $                  LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
*
*           b1 := b1 - V1*w
*
            CALL DTRMV( 'Lower', 'NO TRANSPOSE', 
     $                  'UNIT', I-1,
     $                  A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL DAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
*
            A( K+I-1, I-1 ) = EI
         END IF
*
*        Generate the elementary reflector H(I) to annihilate
*        A(K+I+1:N,I)
*
         CALL DLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1,
     $                TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
*
*        Compute  Y(K+1:N,I)
*
         CALL DGEMV( 'NO TRANSPOSE', N-K, N-K-I+1, 
     $               ONE, A( K+1, I+1 ),
     $               LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
         CALL DGEMV( 'Transpose', N-K-I+1, I-1, 
     $               ONE, A( K+I, 1 ), LDA,
     $               A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL DGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, 
     $               Y( K+1, 1 ), LDY,
     $               T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
         CALL DSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
*
*        Compute T(1:I,I)
*
         CALL DSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL DTRMV( 'Upper', 'No Transpose', 'NON-UNIT', 
     $               I-1, T, LDT,
     $               T( 1, I ), 1 )
         T( I, I ) = TAU( I )
*
   10 CONTINUE
      A( K+NB, NB ) = EI
*
*     Compute Y(1:K,1:NB)
*
      CALL DLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
      CALL DTRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE', 
     $            'UNIT', K, NB,
     $            ONE, A( K+1, 1 ), LDA, Y, LDY )
      IF( N.GT.K+NB )
     $   CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, 
     $               NB, N-K-NB, ONE,
     $               A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y,
     $               LDY )
      CALL DTRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE', 
     $            'NON-UNIT', K, NB,
     $            ONE, T, LDT, Y, LDY )
*
      RETURN
*
*     End of DLAHR2
*
      END
      LOGICAL FUNCTION DLAISNAN(DIN1,DIN2)
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DIN1,DIN2
*     ..
*
*  Purpose
*  =======
*
*  This routine is not for general use.  It exists solely to avoid
*  over-optimization in DISNAN.
*
*  DLAISNAN checks for NaNs by comparing its two arguments for
*  inequality.  NaN is the only floating-point value where NaN != NaN
*  returns .TRUE.  To check for NaNs, pass the same variable as both
*  arguments.
*
*  A compiler must assume that the two arguments are
*  not the same variable, and the test will not be optimized away.
*  Interprocedural or whole-program optimization may delete this
*  test.  The ISNAN functions will be replaced by the correct
*  Fortran 03 intrinsic once the intrinsic is widely available.
*
*  Arguments
*  =========
*
*  DIN1     (input) DOUBLE PRECISION
*  DIN2     (input) DOUBLE PRECISION
*          Two numbers to compare for inequality.
*
*  =====================================================================
*
*  .. Executable Statements ..
      DLAISNAN = (DIN1.NE.DIN2)
      RETURN
      END
      SUBROUTINE DLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B,
     $                   LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      DOUBLE PRECISION   CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DLALN2 solves a system of the form  (ca A - w D ) X = s B
*  or (ca A' - w D) X = s B   with possible scaling ("s") and
*  perturbation of A.  (A' means A-transpose.)
*
*  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
*  real diagonal matrix, w is a real or complex value, and X and B are
*  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
*  may be 1 or 2.
*
*  If w is complex, X and B are represented as NA x 2 matrices,
*  the first column of each being the real part and the second
*  being the imaginary part.
*
*  "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
*  so chosen that X can be computed without overflow.  X is further
*  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
*  than overflow.
*
*  If both singular values of (ca A - w D) are less than SMIN,
*  SMIN*identity will be used instead of (ca A - w D).  If only one
*  singular value is less than SMIN, one element of (ca A - w D) will be
*  perturbed enough to make the smallest singular value roughly SMIN.
*  If both singular values are at least SMIN, (ca A - w D) will not be
*  perturbed.  In any case, the perturbation will be at most some small
*  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
*  are computed by infinity-norm approximations, and thus will only be
*  correct to a factor of 2 or so.
*
*  Note: all input quantities are assumed to be smaller than overflow
*  by a reasonable factor.  (See BIGNUM.)
*
*  Arguments
*  ==========
*
*  LTRANS  (input) LOGICAL
*          =.TRUE.:  A-transpose will be used.
*          =.FALSE.: A will be used (not transposed.)
*
*  NA      (input) INTEGER
*          The size of the matrix A.  It may (only) be 1 or 2.
*
*  NW      (input) INTEGER
*          1 if "w" is real, 2 if "w" is complex.  It may only be 1
*          or 2.
*
*  SMIN    (input) DOUBLE PRECISION
*          The desired lower bound on the singular values of A.  This
*          should be a safe distance away from underflow or overflow,
*          say, between (underflow/machine precision) and  (machine
*          precision * overflow ).  (See BIGNUM and ULP.)
*
*  CA      (input) DOUBLE PRECISION
*          The coefficient c, which A is multiplied by.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,NA)
*          The NA x NA matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least NA.
*
*  D1      (input) DOUBLE PRECISION
*          The 1,1 element in the diagonal matrix D.
*
*  D2      (input) DOUBLE PRECISION
*          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NW)
*          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
*          complex), column 1 contains the real part of B and column 2
*          contains the imaginary part.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  It must be at least NA.
*
*  WR      (input) DOUBLE PRECISION
*          The real part of the scalar "w".
*
*  WI      (input) DOUBLE PRECISION
*          The imaginary part of the scalar "w".  Not used if NW=1.
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NW)
*          The NA x NW matrix X (unknowns), as computed by DLALN2.
*          If NW=2 ("w" is complex), on exit, column 1 will contain
*          the real part of X and column 2 will contain the imaginary
*          part.
*
*  LDX     (input) INTEGER
*          The leading dimension of X.  It must be at least NA.
*
*  SCALE   (output) DOUBLE PRECISION
*          The scale factor that B must be multiplied by to insure
*          that overflow does not occur when computing X.  Thus,
*          (ca A - w D) X  will be SCALE*B, not B (ignoring
*          perturbations of A.)  It will be at most 1.
*
*  XNORM   (output) DOUBLE PRECISION
*          The infinity-norm of X, when X is regarded as an NA x NW
*          real matrix.
*
*  INFO    (output) INTEGER
*          An error flag.  It will be set to zero if no error occurs,
*          a negative number if an argument is in error, or a positive
*          number if  ca A - w D  had to be perturbed.
*          The possible values are:
*          = 0: No error occurred, and (ca A - w D) did not have to be
*                 perturbed.
*          = 1: (ca A - w D) had to be perturbed to make its smallest
*               (or only) singular value greater than SMIN.
*          NOTE: In the interests of speed, this routine does not
*                check the inputs for errors.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICMAX, J
      DOUBLE PRECISION   BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21,
     $                   CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21,
     $                   LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R,
     $                   UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S,
     $                   UR22, XI1, XI2, XR1, XR2
*     ..
*     .. Local Arrays ..
      LOGICAL            RSWAP( 4 ), ZSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      DOUBLE PRECISION   CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Equivalences ..
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ),
     $                   ( CR( 1, 1 ), CRV( 1 ) )
*     ..
*     .. Data statements ..
      DATA               ZSWAP / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               RSWAP / .FALSE., .TRUE., .FALSE., .TRUE. /
      DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4,
     $                   3, 2, 1 /
*     ..
*     .. Executable Statements ..
*
*     Compute BIGNUM
*
      SMLNUM = TWO*DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMINI = MAX( SMIN, SMLNUM )
*
*     Don't check for input errors
*
      INFO = 0
*
*     Standard Initializations
*
      SCALE = ONE
*
      IF( NA.EQ.1 ) THEN
*
*        1 x 1  (i.e., scalar) system   C X = B
*
         IF( NW.EQ.1 ) THEN
*
*           Real 1x1 system.
*
*           C = ca A - w D
*
            CSR = CA*A( 1, 1 ) - WR*D1
            CNORM = ABS( CSR )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CNORM = SMINI
               INFO = 1
            END IF
*
*           Check scaling for  X = B / C
*
            BNORM = ABS( B( 1, 1 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM )
     $            SCALE = ONE / BNORM
            END IF
*
*           Compute X
*
            X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR
            XNORM = ABS( X( 1, 1 ) )
         ELSE
*
*           Complex 1x1 system (w is complex)
*
*           C = ca A - w D
*
            CSR = CA*A( 1, 1 ) - WR*D1
            CSI = -WI*D1
            CNORM = ABS( CSR ) + ABS( CSI )
*
*           If | C | < SMINI, use C = SMINI
*
            IF( CNORM.LT.SMINI ) THEN
               CSR = SMINI
               CSI = ZERO
               CNORM = SMINI
               INFO = 1
            END IF
*
*           Check scaling for  X = B / C
*
            BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) )
            IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
               IF( BNORM.GT.BIGNUM*CNORM )
     $            SCALE = ONE / BNORM
            END IF
*
*           Compute X
*
            CALL DLADIV( SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI,
     $                   X( 1, 1 ), X( 1, 2 ) )
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
         END IF
*
      ELSE
*
*        2x2 System
*
*        Compute the real part of  C = ca A - w D  (or  ca A' - w D )
*
         CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1
         CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2
         IF( LTRANS ) THEN
            CR( 1, 2 ) = CA*A( 2, 1 )
            CR( 2, 1 ) = CA*A( 1, 2 )
         ELSE
            CR( 2, 1 ) = CA*A( 2, 1 )
            CR( 1, 2 ) = CA*A( 1, 2 )
         END IF
*
         IF( NW.EQ.1 ) THEN
*
*           Real 2x2 system  (w is real)
*
*           Find the largest element in C
*
            CMAX = ZERO
            ICMAX = 0
*
            DO 10 J = 1, 4
               IF( ABS( CRV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) )
                  ICMAX = J
               END IF
   10       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI )
     $               SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            UR11 = CRV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            UR11R = ONE / UR11
            LR21 = UR11R*CR21
            UR22 = CR22 - UR12*LR21
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( ABS( UR22 ).LT.SMINI ) THEN
               UR22 = SMINI
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR1 = B( 2, 1 )
               BR2 = B( 1, 1 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
            END IF
            BR2 = BR2 - LR21*BR1
            BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) )
            IF( BBND.GT.ONE .AND. ABS( UR22 ).LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*ABS( UR22 ) )
     $            SCALE = ONE / BBND
            END IF
*
            XR2 = ( BR2*SCALE ) / UR22
            XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 )
            IF( ZSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
            END IF
            XNORM = MAX( ABS( XR1 ), ABS( XR2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         ELSE
*
*           Complex 2x2 system  (w is complex)
*
*           Find the largest element in C
*
            CI( 1, 1 ) = -WI*D1
            CI( 2, 1 ) = ZERO
            CI( 1, 2 ) = ZERO
            CI( 2, 2 ) = -WI*D2
            CMAX = ZERO
            ICMAX = 0
*
            DO 20 J = 1, 4
               IF( ABS( CRV( J ) )+ABS( CIV( J ) ).GT.CMAX ) THEN
                  CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) )
                  ICMAX = J
               END IF
   20       CONTINUE
*
*           If norm(C) < SMINI, use SMINI*identity.
*
            IF( CMAX.LT.SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ),
     $                 ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
               IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                  IF( BNORM.GT.BIGNUM*SMINI )
     $               SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               X( 1, 2 ) = TEMP*B( 1, 2 )
               X( 2, 2 ) = TEMP*B( 2, 2 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
*
*           Gaussian elimination with complete pivoting.
*
            UR11 = CRV( ICMAX )
            UI11 = CIV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            CI21 = CIV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            UI12 = CIV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            CI22 = CIV( IPIVOT( 4, ICMAX ) )
            IF( ICMAX.EQ.1 .OR. ICMAX.EQ.4 ) THEN
*
*              Code when off-diagonals of pivoted C are real
*
               IF( ABS( UR11 ).GT.ABS( UI11 ) ) THEN
                  TEMP = UI11 / UR11
                  UR11R = ONE / ( UR11*( ONE+TEMP**2 ) )
                  UI11R = -TEMP*UR11R
               ELSE
                  TEMP = UR11 / UI11
                  UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) )
                  UR11R = -TEMP*UI11R
               END IF
               LR21 = CR21*UR11R
               LI21 = CR21*UI11R
               UR12S = UR12*UR11R
               UI12S = UR12*UI11R
               UR22 = CR22 - UR12*LR21
               UI22 = CI22 - UR12*LI21
            ELSE
*
*              Code when diagonals of pivoted C are real
*
               UR11R = ONE / UR11
               UI11R = ZERO
               LR21 = CR21*UR11R
               LI21 = CI21*UR11R
               UR12S = UR12*UR11R
               UI12S = UI12*UR11R
               UR22 = CR22 - UR12*LR21 + UI12*LI21
               UI22 = -UR12*LI21 - UI12*LR21
            END IF
            U22ABS = ABS( UR22 ) + ABS( UI22 )
*
*           If smaller pivot < SMINI, use SMINI
*
            IF( U22ABS.LT.SMINI ) THEN
               UR22 = SMINI
               UI22 = ZERO
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR2 = B( 1, 1 )
               BR1 = B( 2, 1 )
               BI2 = B( 1, 2 )
               BI1 = B( 2, 2 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
               BI1 = B( 1, 2 )
               BI2 = B( 2, 2 )
            END IF
            BR2 = BR2 - LR21*BR1 + LI21*BI1
            BI2 = BI2 - LI21*BR1 - LR21*BI1
            BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )*
     $             ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ),
     $             ABS( BR2 )+ABS( BI2 ) )
            IF( BBND.GT.ONE .AND. U22ABS.LT.ONE ) THEN
               IF( BBND.GE.BIGNUM*U22ABS ) THEN
                  SCALE = ONE / BBND
                  BR1 = SCALE*BR1
                  BI1 = SCALE*BI1
                  BR2 = SCALE*BR2
                  BI2 = SCALE*BI2
               END IF
            END IF
*
            CALL DLADIV( BR2, BI2, UR22, UI22, XR2, XI2 )
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
            IF( ZSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
               X( 1, 2 ) = XI2
               X( 2, 2 ) = XI1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
               X( 1, 2 ) = XI1
               X( 2, 2 ) = XI2
            END IF
            XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) )
*
*           Further scaling if  norm(A) norm(X) > overflow
*
            IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
               IF( XNORM.GT.BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  X( 1, 2 ) = TEMP*X( 1, 2 )
                  X( 2, 2 ) = TEMP*X( 2, 2 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DLALN2
*
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
         FIRST = .FALSE.
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) INTEGER 
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix A.
*
*  Description
*  ===========
*
*  DLANGE returns the value
*
*     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANGE as described
*          above.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.  When M = 0,
*          DLANGE is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.  When N = 0,
*          DLANGE is set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANGE = VALUE
      RETURN
*
*     End of DLANGE
*
      END
      DOUBLE PRECISION FUNCTION DLANHS( NORM, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANHS  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  Hessenberg matrix A.
*
*  Description
*  ===========
*
*  DLANHS returns the value
*
*     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANHS as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
*          set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The n by n upper Hessenberg matrix A; the part of A below the
*          first sub-diagonal is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANHS = VALUE
      RETURN
*
*     End of DLANHS
*
      END
      DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
     $                 WORK )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANTR  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  trapezoidal or triangular matrix A.
*
*  Description
*  ===========
*
*  DLANTR returns the value
*
*     DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANTR as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower trapezoidal.
*          = 'U':  Upper trapezoidal
*          = 'L':  Lower trapezoidal
*          Note that A is triangular instead of trapezoidal if M = N.
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A has unit diagonal.
*          = 'N':  Non-unit diagonal
*          = 'U':  Unit diagonal
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0, and if
*          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0, and if
*          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The trapezoidal matrix A (A is triangular if M = N).
*          If UPLO = 'U', the leading m by n upper trapezoidal part of
*          the array A contains the upper trapezoidal matrix, and the
*          strictly lower triangular part of A is not referenced.
*          If UPLO = 'L', the leading m by n lower trapezoidal part of
*          the array A contains the lower trapezoidal matrix, and the
*          strictly upper triangular part of A is not referenced.  Note
*          that when DIAG = 'U', the diagonal elements of A are not
*          referenced and are assumed to be one.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UDIAG
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         IF( LSAME( DIAG, 'U' ) ) THEN
            VALUE = ONE
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 20 J = 1, N
                  DO 10 I = 1, MIN( M, J-1 )
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = J + 1, M
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 60 J = 1, N
                  DO 50 I = 1, MIN( M, J )
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 70 I = J, M
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         UDIAG = LSAME( DIAG, 'U' )
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 1, N
               IF( ( UDIAG ) .AND. ( J.LE.M ) ) THEN
                  SUM = ONE
                  DO 90 I = 1, J - 1
                     SUM = SUM + ABS( A( I, J ) )
   90             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 100 I = 1, MIN( M, J )
                     SUM = SUM + ABS( A( I, J ) )
  100             CONTINUE
               END IF
               VALUE = MAX( VALUE, SUM )
  110       CONTINUE
         ELSE
            DO 140 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 120 I = J + 1, M
                     SUM = SUM + ABS( A( I, J ) )
  120             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 130 I = J, M
                     SUM = SUM + ABS( A( I, J ) )
  130             CONTINUE
               END IF
               VALUE = MAX( VALUE, SUM )
  140       CONTINUE
         END IF
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 150 I = 1, M
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  DO 160 I = 1, MIN( M, J-1 )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  160             CONTINUE
  170          CONTINUE
            ELSE
               DO 180 I = 1, M
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  DO 190 I = 1, MIN( M, J )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  190             CONTINUE
  200          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 210 I = 1, N
                  WORK( I ) = ONE
  210          CONTINUE
               DO 220 I = N + 1, M
                  WORK( I ) = ZERO
  220          CONTINUE
               DO 240 J = 1, N
                  DO 230 I = J + 1, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  230             CONTINUE
  240          CONTINUE
            ELSE
               DO 250 I = 1, M
                  WORK( I ) = ZERO
  250          CONTINUE
               DO 270 J = 1, N
                  DO 260 I = J, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  260             CONTINUE
  270          CONTINUE
            END IF
         END IF
         VALUE = ZERO
         DO 280 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
  280    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = MIN( M, N )
               DO 290 J = 2, N
                  CALL DLASSQ( MIN( M, J-1 ), A( 1, J ), 1, SCALE, SUM )
  290          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 300 J = 1, N
                  CALL DLASSQ( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
  300          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = MIN( M, N )
               DO 310 J = 1, N
                  CALL DLASSQ( M-J, A( MIN( M, J+1 ), J ), 1, SCALE,
     $                         SUM )
  310          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 320 J = 1, N
                  CALL DLASSQ( M-J+1, A( J, J ), 1, SCALE, SUM )
  320          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANTR = VALUE
      RETURN
*
*     End of DLANTR
*
      END
      SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
*  matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
*
*  where either
*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
*  conjugate eigenvalues.
*
*  Arguments
*  =========
*
*  A       (input/output) DOUBLE PRECISION
*  B       (input/output) DOUBLE PRECISION
*  C       (input/output) DOUBLE PRECISION
*  D       (input/output) DOUBLE PRECISION
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1R    (output) DOUBLE PRECISION
*  RT1I    (output) DOUBLE PRECISION
*  RT2R    (output) DOUBLE PRECISION
*  RT2I    (output) DOUBLE PRECISION
*          The real and imaginary parts of the eigenvalues. If the
*          eigenvalues are a complex conjugate pair, RT1I > 0.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) DOUBLE PRECISION
*          Parameters of the rotation matrix.
*
*  Further Details
*  ===============
*
*  Modified by V. Sima, Research Institute for Informatics, Bucharest,
*  Romania, to reduce the risk of cancellation errors,
*  when computing real eigenvalues, and to ensure, if possible, that
*  abs(RT1R) >= abs(RT2R).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   MULTPL
      PARAMETER          ( MULTPL = 4.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB,
     $                   SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'P' )
      IF( C.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         GO TO 10
*
      ELSE IF( B.EQ.ZERO ) THEN
*
*        Swap rows and columns
*
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( ( A-D ).EQ.ZERO .AND. SIGN( ONE, B ).NE.SIGN( ONE, C ) )
     $          THEN
         CS = ONE
         SN = ZERO
         GO TO 10
      ELSE
*
         TEMP = A - D
         P = HALF*TEMP
         BCMAX = MAX( ABS( B ), ABS( C ) )
         BCMIS = MIN( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C )
         SCALE = MAX( ABS( P ), BCMAX )
         Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS
*
*        If Z is of the order of the machine accuracy, postpone the
*        decision on the nature of eigenvalues
*
         IF( Z.GE.MULTPL*EPS ) THEN
*
*           Real eigenvalues. Compute A and D.
*
            Z = P + SIGN( SQRT( SCALE )*SQRT( Z ), P )
            A = D + Z
            D = D - ( BCMAX / Z )*BCMIS
*
*           Compute B and the rotation matrix
*
            TAU = DLAPY2( C, Z )
            CS = Z / TAU
            SN = C / TAU
            B = B - C
            C = ZERO
         ELSE
*
*           Complex eigenvalues, or real (almost) equal eigenvalues.
*           Make diagonal elements equal.
*
            SIGMA = B + C
            TAU = DLAPY2( SIGMA, TEMP )
            CS = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
            SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA )
*
*           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
*                   [ CC  DD ]   [ C  D ] [ SN  CS ]
*
            AA = A*CS + B*SN
            BB = -A*SN + B*CS
            CC = C*CS + D*SN
            DD = -C*SN + D*CS
*
*           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
*                   [ C  D ]   [-SN  CS ] [ CC  DD ]
*
            A = AA*CS + CC*SN
            B = BB*CS + DD*SN
            C = -AA*SN + CC*CS
            D = -BB*SN + DD*CS
*
            TEMP = HALF*( A+D )
            A = TEMP
            D = TEMP
*
            IF( C.NE.ZERO ) THEN
               IF( B.NE.ZERO ) THEN
                  IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
*
*                    Real eigenvalues: reduce to upper triangular form
*
                     SAB = SQRT( ABS( B ) )
                     SAC = SQRT( ABS( C ) )
                     P = SIGN( SAB*SAC, C )
                     TAU = ONE / SQRT( ABS( B+C ) )
                     A = TEMP + P
                     D = TEMP - P
                     B = B - C
                     C = ZERO
                     CS1 = SAB*TAU
                     SN1 = SAC*TAU
                     TEMP = CS*CS1 - SN*SN1
                     SN = CS*SN1 + SN*CS1
                     CS = TEMP
                  END IF
               ELSE
                  B = -C
                  C = ZERO
                  TEMP = CS
                  CS = -SN
                  SN = TEMP
               END IF
            END IF
         END IF
*
      END IF
*
   10 CONTINUE
*
*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
*
      RT1R = A
      RT2R = D
      IF( C.EQ.ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN
*
*     End of DLANV2
*
      END
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
      SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   EQUED )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED
      INTEGER            LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAQGE equilibrates a general M by N matrix A using the row and
*  column scaling factors in the vectors R and C.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M by N matrix A.
*          On exit, the equilibrated matrix.  See EQUED for the form of
*          the equilibrated matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  R       (input) DOUBLE PRECISION array, dimension (M)
*          The row scale factors for A.
*
*  C       (input) DOUBLE PRECISION array, dimension (N)
*          The column scale factors for A.
*
*  ROWCND  (input) DOUBLE PRECISION
*          Ratio of the smallest R(i) to the largest R(i).
*
*  COLCND  (input) DOUBLE PRECISION
*          Ratio of the smallest C(i) to the largest C(i).
*
*  AMAX    (input) DOUBLE PRECISION
*          Absolute value of largest matrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration
*          = 'R':  Row equilibration, i.e., A has been premultiplied by
*                  diag(R).
*          = 'C':  Column equilibration, i.e., A has been postmultiplied
*                  by diag(C).
*          = 'B':  Both row and column equilibration, i.e., A has been
*                  replaced by diag(R) * A * diag(C).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if row or column scaling
*  should be done based on the ratio of the row or column scaling
*  factors.  If ROWCND < THRESH, row scaling is done, and if
*  COLCND < THRESH, column scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if row scaling
*  should be done based on the absolute size of the largest matrix
*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( ROWCND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE )
     $     THEN
*
*        No row scaling
*
         IF( COLCND.GE.THRESH ) THEN
*
*           No column scaling
*
            EQUED = 'N'
         ELSE
*
*           Column scaling
*
            DO 20 J = 1, N
               CJ = C( J )
               DO 10 I = 1, M
                  A( I, J ) = CJ*A( I, J )
   10          CONTINUE
   20       CONTINUE
            EQUED = 'C'
         END IF
      ELSE IF( COLCND.GE.THRESH ) THEN
*
*        Row scaling, no column scaling
*
         DO 40 J = 1, N
            DO 30 I = 1, M
               A( I, J ) = R( I )*A( I, J )
   30       CONTINUE
   40    CONTINUE
         EQUED = 'R'
      ELSE
*
*        Row and column scaling
*
         DO 60 J = 1, N
            CJ = C( J )
            DO 50 I = 1, M
               A( I, J ) = CJ*R( I )*A( I, J )
   50       CONTINUE
   60    CONTINUE
         EQUED = 'B'
      END IF
*
      RETURN
*
*     End of DLAQGE
*
      END
      SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   Z( LDZ, * )
*     ..
*
*     Purpose
*     =======
*
*     DLAQR0 computes the eigenvalues of a Hessenberg matrix H
*     and, optionally, the matrices T and Z from the Schur decomposition
*     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
*     Schur form), and Z is the orthogonal matrix of Schur vectors.
*
*     Optionally Z may be postmultiplied into an input orthogonal
*     matrix Q so that this routine can give the Schur factorization
*     of a matrix A which has been reduced to the Hessenberg form H
*     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
*
*     Arguments
*     =========
*
*     WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*     WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*     N     (input) INTEGER
*           The order of the matrix H.  N .GE. 0.
*
*     ILO   (input) INTEGER
*     IHI   (input) INTEGER
*           It is assumed that H is already upper triangular in rows
*           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
*           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
*           previous call to DGEBAL, and then passed to DGEHRD when the
*           matrix output by DGEBAL is reduced to Hessenberg form.
*           Otherwise, ILO and IHI should be set to 1 and N,
*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*           If N = 0, then ILO = 1 and IHI = 0.
*
*     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*           On entry, the upper Hessenberg matrix H.
*           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
*           the upper quasi-triangular matrix T from the Schur
*           decomposition (the Schur form); 2-by-2 diagonal blocks
*           (corresponding to complex conjugate pairs of eigenvalues)
*           are returned in standard form, with H(i,i) = H(i+1,i+1)
*           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is
*           .FALSE., then the contents of H are unspecified on exit.
*           (The output value of H when INFO.GT.0 is given under the
*           description of INFO below.)
*
*           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
*           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
*
*     LDH   (input) INTEGER
*           The leading dimension of the array H. LDH .GE. max(1,N).
*
*     WR    (output) DOUBLE PRECISION array, dimension (IHI)
*     WI    (output) DOUBLE PRECISION array, dimension (IHI)
*           The real and imaginary parts, respectively, of the computed
*           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
*           and WI(ILO:IHI). If two eigenvalues are computed as a
*           complex conjugate pair, they are stored in consecutive
*           elements of WR and WI, say the i-th and (i+1)th, with
*           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then
*           the eigenvalues are stored in the same order as on the
*           diagonal of the Schur form returned in H, with
*           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
*           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
*           WI(i+1) = -WI(i).
*
*     ILOZ     (input) INTEGER
*     IHIZ     (input) INTEGER
*           Specify the rows of Z to which transformations must be
*           applied if WANTZ is .TRUE..
*           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
*
*     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
*           If WANTZ is .FALSE., then Z is not referenced.
*           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
*           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
*           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
*           (The output value of Z when INFO.GT.0 is given under
*           the description of INFO below.)
*
*     LDZ   (input) INTEGER
*           The leading dimension of the array Z.  if WANTZ is .TRUE.
*           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
*
*     WORK  (workspace/output) DOUBLE PRECISION array, dimension LWORK
*           On exit, if LWORK = -1, WORK(1) returns an estimate of
*           the optimal value for LWORK.
*
*     LWORK (input) INTEGER
*           The dimension of the array WORK.  LWORK .GE. max(1,N)
*           is sufficient, but LWORK typically as large as 6*N may
*           be required for optimal performance.  A workspace query
*           to determine the optimal workspace size is recommended.
*
*           If LWORK = -1, then DLAQR0 does a workspace query.
*           In this case, DLAQR0 checks the input parameters and
*           estimates the optimal workspace size for the given
*           values of N, ILO and IHI.  The estimate is returned
*           in WORK(1).  No error message related to LWORK is
*           issued by XERBLA.  Neither H nor Z are accessed.
*
*
*     INFO  (output) INTEGER
*             =  0:  successful exit
*           .GT. 0:  if INFO = i, DLAQR0 failed to compute all of
*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*                and WI contain those eigenvalues which have been
*                successfully computed.  (Failures are rare.)
*
*                If INFO .GT. 0 and WANT is .FALSE., then on exit,
*                the remaining unconverged eigenvalues are the eigen-
*                values of the upper Hessenberg matrix rows and
*                columns ILO through INFO of the final, output
*                value of H.
*
*                If INFO .GT. 0 and WANTT is .TRUE., then on exit
*
*           (*)  (initial value of H)*U  = U*(final value of H)
*
*                where U is an orthogonal matrix.  The final
*                value of H is upper Hessenberg and quasi-triangular
*                in rows and columns INFO+1 through IHI.
*
*                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
*
*                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
*                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
*
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of WANTT.)
*
*                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
*                accessed.
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*     References:
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*       929--947, 2002.
*
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*       of Matrix Analysis, volume 23, pages 948--973, 2002.
*
*     ================================================================
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      PARAMETER          ( KEXNW = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 6 )
*
*     ==== The constants WILK1 and WILK2 are used to form the
*     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
     $                   LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS,
     $                   NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR3, DLAQR4, DLAQR5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      INFO = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
*
      IF( N.LE.NTINY ) THEN
*
*        ==== Tiny matrices must use DLAHQR. ====
*
         LWKOPT = 1
         IF( LWORK.NE.-1 )
     $      CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         INFO = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( WANTT ) THEN
            JBCMPZ( 1: 1 ) = 'S'
         ELSE
            JBCMPZ( 1: 1 ) = 'E'
         END IF
         IF( WANTZ ) THEN
            JBCMPZ( 2: 2 ) = 'V'
         ELSE
            JBCMPZ( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 11, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.3.) ====
*
         NWR = ILAENV( 13, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 11, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         NSR = ILAENV( 15, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to DLAQR3 ====
*
         CALL DLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,
     $                IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH,
     $                N, H, LDH, WORK, -1 )
*
*        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
*
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         NMIN = ILAENV( 12, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== Nibble crossover point ====
*
         NIBBLE = ILAENV( 14, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         KACC22 = ILAENV( 16, 'DLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
         NW = NWMAX
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         NDFL = 1
*
*        ==== ITMAX = iteration limit ====
*
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
*
*        ==== Last row and column in the active block ====
*
         KBOT = IHI
*
*        ==== Main Loop ====
*
         DO 80 IT = 1, ITMAX
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( KBOT.LT.ILO )
     $         GO TO 90
*
*           ==== Locate active block ====
*
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO )
     $            GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            NH = KBOT - KTOP + 1
            NWUPBD = MIN( NH, NWMAX )
            IF( NDFL.LT.KEXNW ) THEN
               NW = MIN( NWUPBD, NWR )
            ELSE
               NW = MIN( NWUPBD, 2*NW )
            END IF
            IF( NW.LT.NWMAX ) THEN
               IF( NW.GE.NH-1 ) THEN
                  NW = NH
               ELSE
                  KWTOP = KBOT - NW + 1
                  IF( ABS( H( KWTOP, KWTOP-1 ) ).GT.
     $                ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 )
     $            NDEC = 0
               NW = NW - NDEC
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
*
*           ==== Aggressive early deflation ====
*
            CALL DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH,
     $                   NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH,
     $                   WORK, LWORK )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            KBOT = KBOT - LD
*
*           ==== KS points to the shifts. ====
*
            KS = KBOT - LS + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
     $          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if DLAQR3
*              .    did not provide that many shifts. ====
*
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    DLAQR3 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
                     SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                     AA = WILK1*SS + H( I, I )
                     BB = SS
                     CC = WILK2*SS
                     DD = AA
                     CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ),
     $                            WR( I ), WI( I ), CS, SN )
   30             CONTINUE
                  IF( KS.EQ.KTOP ) THEN
                     WR( KS+1 ) = H( KS+1, KS+1 )
                     WI( KS+1 ) = ZERO
                     WR( KS ) = WR( KS+1 )
                     WI( KS ) = WI( KS+1 )
                  END IF
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
*                 .    DLAHQR on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH,
     $                            H( KT, 1 ), LDH )
                     IF( NS.GT.NMIN ) THEN
                        CALL DLAQR4( .false., .false., NS, 1, NS,
     $                               H( KT, 1 ), LDH, WR( KS ),
     $                               WI( KS ), 1, 1, ZDUM, 1, WORK,
     $                               LWORK, INF )
                     ELSE
                        CALL DLAHQR( .false., .false., NS, 1, NS,
     $                               H( KT, 1 ), LDH, WR( KS ),
     $                               WI( KS ), 1, 1, ZDUM, 1, INF )
                     END IF
                     KS = KS + INF
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  ====
*
                     IF( KS.GE.KBOT ) THEN
                        AA = H( KBOT-1, KBOT-1 )
                        CC = H( KBOT, KBOT-1 )
                        BB = H( KBOT-1, KBOT )
                        DD = H( KBOT, KBOT )
                        CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ),
     $                               WI( KBOT-1 ), WR( KBOT ),
     $                               WI( KBOT ), CS, SN )
                        KS = KBOT - 1
                     END IF
                  END IF
*
                  IF( KBOT-KS+1.GT.NS ) THEN
*
*                    ==== Sort the shifts (Helps a little)
*                    .    Bubble sort keeps complex conjugate
*                    .    pairs together. ====
*
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED )
     $                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( ABS( WR( I ) )+ABS( WI( I ) ).LT.
     $                         ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                              SORTED = .false.
*
                              SWAP = WR( I )
                              WR( I ) = WR( I+1 )
                              WR( I+1 ) = SWAP
*
                              SWAP = WI( I )
                              WI( I ) = WI( I+1 )
                              WI( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
*
*                 ==== Shuffle shifts into pairs of real shifts
*                 .    and pairs of complex conjugate shifts
*                 .    assuming complex conjugate shifts are
*                 .    already adjacent to one another. (Yes,
*                 .    they are.)  ====
*
                  DO 70 I = KBOT, KS + 2, -2
                     IF( WI( I ).NE.-WI( I-1 ) ) THEN
*
                        SWAP = WR( I )
                        WR( I ) = WR( I-1 )
                        WR( I-1 ) = WR( I-2 )
                        WR( I-2 ) = SWAP
*
                        SWAP = WI( I )
                        WI( I ) = WI( I-1 )
                        WI( I-1 ) = WI( I-2 )
                        WI( I-2 ) = SWAP
                     END IF
   70             CONTINUE
               END IF
*
*              ==== If there are only two shifts and both are
*              .    real, then use only one.  ====
*
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( WI( KBOT ).EQ.ZERO ) THEN
                     IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT.
     $                   ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                        WR( KBOT-1 ) = WR( KBOT )
                     ELSE
                        WR( KBOT ) = WR( KBOT-1 )
                     END IF
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnatiude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               KDU = 3*NS - 3
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,
     $                      WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z,
     $                      LDZ, WORK, 3, H( KU, 1 ), LDH, NVE,
     $                      H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
*
*           ==== End of main loop ====
   80    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         INFO = KBOT
   90    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     ==== End of DLAQR0 ====
*
      END
      SUBROUTINE DLAQR1( N, H, LDH, SR1, SI1, SR2, SI2, V )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   SI1, SI2, SR1, SR2
      INTEGER            LDH, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), V( * )
*     ..
*
*       Given a 2-by-2 or 3-by-3 matrix H, DLAQR1 sets v to a
*       scalar multiple of the first column of the product
*
*       (*)  K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)
*
*       scaling to avoid overflows and most underflows. It
*       is assumed that either
*
*               1) sr1 = sr2 and si1 = -si2
*           or
*               2) si1 = si2 = 0.
*
*       This is useful for starting double implicit shift bulges
*       in the QR algorithm.
*
*
*       N      (input) integer
*              Order of the matrix H. N must be either 2 or 3.
*
*       H      (input) DOUBLE PRECISION array of dimension (LDH,N)
*              The 2-by-2 or 3-by-3 matrix H in (*).
*
*       LDH    (input) integer
*              The leading dimension of H as declared in
*              the calling procedure.  LDH.GE.N
*
*       SR1    (input) DOUBLE PRECISION
*       SI1    The shifts in (*).
*       SR2
*       SI2
*
*       V      (output) DOUBLE PRECISION array of dimension N
*              A scalar multiple of the first column of the
*              matrix K in (*).
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   H21S, H31S, S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
      IF( N.EQ.2 ) THEN
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            V( 1 ) = H21S*H( 1, 2 ) + ( H( 1, 1 )-SR1 )*
     $               ( ( H( 1, 1 )-SR2 ) / S ) - SI1*( SI2 / S )
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 )
         END IF
      ELSE
         S = ABS( H( 1, 1 )-SR2 ) + ABS( SI2 ) + ABS( H( 2, 1 ) ) +
     $       ABS( H( 3, 1 ) )
         IF( S.EQ.ZERO ) THEN
            V( 1 ) = ZERO
            V( 2 ) = ZERO
            V( 3 ) = ZERO
         ELSE
            H21S = H( 2, 1 ) / S
            H31S = H( 3, 1 ) / S
            V( 1 ) = ( H( 1, 1 )-SR1 )*( ( H( 1, 1 )-SR2 ) / S ) -
     $               SI1*( SI2 / S ) + H( 1, 2 )*H21S + H( 1, 3 )*H31S
            V( 2 ) = H21S*( H( 1, 1 )+H( 2, 2 )-SR1-SR2 ) +
     $               H( 2, 3 )*H31S
            V( 3 ) = H31S*( H( 1, 1 )+H( 3, 3 )-SR1-SR2 ) +
     $               H21S*H( 3, 2 )
         END IF
      END IF
      END
      SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
     $                   LDT, NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*  -- April 2009                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
     $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
     $                   Z( LDZ, * )
*     ..
*
*     This subroutine is identical to DLAQR3 except that it avoids
*     recursion by calling DLAHQR instead of DLAQR4.
*
*
*     ******************************************************************
*     Aggressive early deflation:
*
*     This subroutine accepts as input an upper Hessenberg matrix
*     H and performs an orthogonal similarity transformation
*     designed to detect and deflate fully converged eigenvalues from
*     a trailing principal submatrix.  On output H has been over-
*     written by a new Hessenberg matrix that is a perturbation of
*     an orthogonal similarity transformation of H.  It is to be
*     hoped that the final version of H has many zero subdiagonal
*     entries.
*
*     ******************************************************************
*     WANTT   (input) LOGICAL
*          If .TRUE., then the Hessenberg matrix H is fully updated
*          so that the quasi-triangular Schur factor may be
*          computed (in cooperation with the calling subroutine).
*          If .FALSE., then only enough of H is updated to preserve
*          the eigenvalues.
*
*     WANTZ   (input) LOGICAL
*          If .TRUE., then the orthogonal matrix Z is updated so
*          so that the orthogonal Schur factor may be computed
*          (in cooperation with the calling subroutine).
*          If .FALSE., then Z is not referenced.
*
*     N       (input) INTEGER
*          The order of the matrix H and (if WANTZ is .TRUE.) the
*          order of the orthogonal matrix Z.
*
*     KTOP    (input) INTEGER
*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*          KBOT and KTOP together determine an isolated block
*          along the diagonal of the Hessenberg matrix.
*
*     KBOT    (input) INTEGER
*          It is assumed without a check that either
*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*          determine an isolated block along the diagonal of the
*          Hessenberg matrix.
*
*     NW      (input) INTEGER
*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*
*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*          On input the initial N-by-N section of H stores the
*          Hessenberg matrix undergoing aggressive early deflation.
*          On output H has been transformed by an orthogonal
*          similarity transformation, perturbed, and the returned
*          to Hessenberg form that (it is to be hoped) has some
*          zero subdiagonal entries.
*
*     LDH     (input) integer
*          Leading dimension of H just as declared in the calling
*          subroutine.  N .LE. LDH
*
*     ILOZ    (input) INTEGER
*     IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*
*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          IF WANTZ is .TRUE., then on output, the orthogonal
*          similarity transformation mentioned above has been
*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
*          If WANTZ is .FALSE., then Z is unreferenced.
*
*     LDZ     (input) integer
*          The leading dimension of Z just as declared in the
*          calling subroutine.  1 .LE. LDZ.
*
*     NS      (output) integer
*          The number of unconverged (ie approximate) eigenvalues
*          returned in SR and SI that may be used as shifts by the
*          calling subroutine.
*
*     ND      (output) integer
*          The number of converged eigenvalues uncovered by this
*          subroutine.
*
*     SR      (output) DOUBLE PRECISION array, dimension KBOT
*     SI      (output) DOUBLE PRECISION array, dimension KBOT
*          On output, the real and imaginary parts of approximate
*          eigenvalues that may be used for shifts are stored in
*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
*          The real and imaginary parts of converged eigenvalues
*          are stored in SR(KBOT-ND+1) through SR(KBOT) and
*          SI(KBOT-ND+1) through SI(KBOT), respectively.
*
*     V       (workspace) DOUBLE PRECISION array, dimension (LDV,NW)
*          An NW-by-NW work array.
*
*     LDV     (input) integer scalar
*          The leading dimension of V just as declared in the
*          calling subroutine.  NW .LE. LDV
*
*     NH      (input) integer scalar
*          The number of columns of T.  NH.GE.NW.
*
*     T       (workspace) DOUBLE PRECISION array, dimension (LDT,NW)
*
*     LDT     (input) integer
*          The leading dimension of T just as declared in the
*          calling subroutine.  NW .LE. LDT
*
*     NV      (input) integer
*          The number of rows of work array WV available for
*          workspace.  NV.GE.NW.
*
*     WV      (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)
*
*     LDWV    (input) integer
*          The leading dimension of W just as declared in the
*          calling subroutine.  NW .LE. LDV
*
*     WORK    (workspace) DOUBLE PRECISION array, dimension LWORK.
*          On exit, WORK(1) is set to an estimate of the optimal value
*          of LWORK for the given values of N, NW, KTOP and KBOT.
*
*     LWORK   (input) integer
*          The dimension of the work array WORK.  LWORK = 2*NW
*          suffices, but greater efficiency may result from larger
*          values of LWORK.
*
*          If LWORK = -1, then a workspace query is assumed; DLAQR2
*          only estimates the optimal workspace size for the given
*          values of N, NW, KTOP and KBOT.  The estimate is returned
*          in WORK(1).  No error message related to LWORK is issued
*          by XERBLA.  Neither H nor Z are accessed.
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
     $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
     $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2,
     $                   LWKOPT
      LOGICAL            BULGE, SORTED
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEHRD, DGEMM, DLABAD, DLACPY, DLAHQR,
     $                   DLANV2, DLARF, DLARFG, DLASET, DORMHR, DTREXC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
*
*        ==== Workspace query call to DGEHRD ====
*
         CALL DGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
*
*        ==== Workspace query call to DORMHR ====
*
         CALL DORMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV,
     $                WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
*
*        ==== Optimal workspace ====
*
         LWKOPT = JW + MAX( LWK1, LWK2 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( NW.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     ==== Setup deflation window ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
*
      IF( KBOT.EQ.KWTOP ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         SR( KWTOP ) = H( KWTOP, KWTOP )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) )
     $        THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP )
     $         H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         WORK( 1 ) = ONE
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL DLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL DCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
*
      CALL DLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      CALL DLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ),
     $             SI( KWTOP ), 1, JW, V, LDV, INFQR )
*
*     ==== DTREXC needs a clean margin near the diagonal ====
*
      DO 10 J = 1, JW - 3
         T( J+2, J ) = ZERO
         T( J+3, J ) = ZERO
   10 CONTINUE
      IF( JW.GT.2 )
     $   T( JW, JW-2 ) = ZERO
*
*     ==== Deflation detection loop ====
*
      NS = JW
      ILST = INFQR + 1
   20 CONTINUE
      IF( ILST.LE.NS ) THEN
         IF( NS.EQ.1 ) THEN
            BULGE = .FALSE.
         ELSE
            BULGE = T( NS, NS-1 ).NE.ZERO
         END IF
*
*        ==== Small spike tip test for deflation ====
*
         IF( .NOT.BULGE ) THEN
*
*           ==== Real eigenvalue ====
*
            FOO = ABS( T( NS, NS ) )
            IF( FOO.EQ.ZERO )
     $         FOO = ABS( S )
            IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
*
*              ==== Deflatable ====
*
               NS = NS - 1
            ELSE
*
*              ==== Undeflatable.   Move it up out of the way.
*              .    (DTREXC can not fail in this case.) ====
*
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK,
     $                      INFO )
               ILST = ILST + 1
            END IF
         ELSE
*
*           ==== Complex conjugate pair ====
*
            FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )*
     $            SQRT( ABS( T( NS-1, NS ) ) )
            IF( FOO.EQ.ZERO )
     $         FOO = ABS( S )
            IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE.
     $          MAX( SMLNUM, ULP*FOO ) ) THEN
*
*              ==== Deflatable ====
*
               NS = NS - 2
            ELSE
*
*              ==== Undeflatable. Move them up out of the way.
*              .    Fortunately, DTREXC does the right thing with
*              .    ILST in case of a rare exchange failure. ====
*
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK,
     $                      INFO )
               ILST = ILST + 2
            END IF
         END IF
*
*        ==== End deflation detection loop ====
*
         GO TO 20
      END IF
*
*        ==== Return to Hessenberg form ====
*
      IF( NS.EQ.0 )
     $   S = ZERO
*
      IF( NS.LT.JW ) THEN
*
*        ==== sorting diagonal blocks of T improves accuracy for
*        .    graded matrices.  Bubble sort deals well with
*        .    exchange failures. ====
*
         SORTED = .false.
         I = NS + 1
   30    CONTINUE
         IF( SORTED )
     $      GO TO 50
         SORTED = .true.
*
         KEND = I - 1
         I = INFQR + 1
         IF( I.EQ.NS ) THEN
            K = I + 1
         ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
            K = I + 1
         ELSE
            K = I + 2
         END IF
   40    CONTINUE
         IF( K.LE.KEND ) THEN
            IF( K.EQ.I+1 ) THEN
               EVI = ABS( T( I, I ) )
            ELSE
               EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )*
     $               SQRT( ABS( T( I, I+1 ) ) )
            END IF
*
            IF( K.EQ.KEND ) THEN
               EVK = ABS( T( K, K ) )
            ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
               EVK = ABS( T( K, K ) )
            ELSE
               EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )*
     $               SQRT( ABS( T( K, K+1 ) ) )
            END IF
*
            IF( EVI.GE.EVK ) THEN
               I = K
            ELSE
               SORTED = .false.
               IFST = I
               ILST = K
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK,
     $                      INFO )
               IF( INFO.EQ.0 ) THEN
                  I = ILST
               ELSE
                  I = K
               END IF
            END IF
            IF( I.EQ.KEND ) THEN
               K = I + 1
            ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
               K = I + 1
            ELSE
               K = I + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      I = JW
   60 CONTINUE
      IF( I.GE.INFQR+1 ) THEN
         IF( I.EQ.INFQR+1 ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE
            AA = T( I-1, I-1 )
            CC = T( I, I-1 )
            BB = T( I-1, I )
            DD = T( I, I )
            CALL DLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ),
     $                   SI( KWTOP+I-2 ), SR( KWTOP+I-1 ),
     $                   SI( KWTOP+I-1 ), CS, SN )
            I = I - 2
         END IF
         GO TO 60
      END IF
*
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL DCOPY( NS, V, LDV, WORK, 1 )
            BETA = WORK( 1 )
            CALL DLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
*
            CALL DLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
*
            CALL DLARF( 'L', NS, JW, WORK, 1, TAU, T, LDT,
     $                  WORK( JW+1 ) )
            CALL DLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT,
     $                  WORK( JW+1 ) )
            CALL DLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV,
     $                  WORK( JW+1 ) )
*
            CALL DGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
     $                   LWORK-JW, INFO )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( KWTOP.GT.1 )
     $      H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
         CALL DLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL DCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),
     $               LDH+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( NS.GT.1 .AND. S.NE.ZERO )
     $      CALL DORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV,
     $                   WORK( JW+1 ), LWORK-JW, INFO )
*
*        ==== Update vertical slab in H ====
*
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 70 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),
     $                  LDH, V, LDV, ZERO, WV, LDWV )
            CALL DLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   70    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( WANTT ) THEN
            DO 80 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL DGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,
     $                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL DLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),
     $                      LDH )
   80       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( WANTZ ) THEN
            DO 90 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ),
     $                     LDZ, V, LDV, ZERO, WV, LDWV )
               CALL DLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),
     $                      LDZ )
   90       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      ND = JW - NS
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      NS = NS - INFQR
*
*      ==== Return optimal workspace. ====
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     ==== End of DLAQR2 ====
*
      END
      SUBROUTINE DLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
     $                   LDT, NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*  -- April 2009                                                      --
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
     $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
     $                   Z( LDZ, * )
*     ..
*
*     ******************************************************************
*     Aggressive early deflation:
*
*     This subroutine accepts as input an upper Hessenberg matrix
*     H and performs an orthogonal similarity transformation
*     designed to detect and deflate fully converged eigenvalues from
*     a trailing principal submatrix.  On output H has been over-
*     written by a new Hessenberg matrix that is a perturbation of
*     an orthogonal similarity transformation of H.  It is to be
*     hoped that the final version of H has many zero subdiagonal
*     entries.
*
*     ******************************************************************
*     WANTT   (input) LOGICAL
*          If .TRUE., then the Hessenberg matrix H is fully updated
*          so that the quasi-triangular Schur factor may be
*          computed (in cooperation with the calling subroutine).
*          If .FALSE., then only enough of H is updated to preserve
*          the eigenvalues.
*
*     WANTZ   (input) LOGICAL
*          If .TRUE., then the orthogonal matrix Z is updated so
*          so that the orthogonal Schur factor may be computed
*          (in cooperation with the calling subroutine).
*          If .FALSE., then Z is not referenced.
*
*     N       (input) INTEGER
*          The order of the matrix H and (if WANTZ is .TRUE.) the
*          order of the orthogonal matrix Z.
*
*     KTOP    (input) INTEGER
*          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*          KBOT and KTOP together determine an isolated block
*          along the diagonal of the Hessenberg matrix.
*
*     KBOT    (input) INTEGER
*          It is assumed without a check that either
*          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*          determine an isolated block along the diagonal of the
*          Hessenberg matrix.
*
*     NW      (input) INTEGER
*          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*
*     H       (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*          On input the initial N-by-N section of H stores the
*          Hessenberg matrix undergoing aggressive early deflation.
*          On output H has been transformed by an orthogonal
*          similarity transformation, perturbed, and the returned
*          to Hessenberg form that (it is to be hoped) has some
*          zero subdiagonal entries.
*
*     LDH     (input) integer
*          Leading dimension of H just as declared in the calling
*          subroutine.  N .LE. LDH
*
*     ILOZ    (input) INTEGER
*     IHIZ    (input) INTEGER
*          Specify the rows of Z to which transformations must be
*          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*
*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          IF WANTZ is .TRUE., then on output, the orthogonal
*          similarity transformation mentioned above has been
*          accumulated into Z(ILOZ:IHIZ,ILO:IHI) from the right.
*          If WANTZ is .FALSE., then Z is unreferenced.
*
*     LDZ     (input) integer
*          The leading dimension of Z just as declared in the
*          calling subroutine.  1 .LE. LDZ.
*
*     NS      (output) integer
*          The number of unconverged (ie approximate) eigenvalues
*          returned in SR and SI that may be used as shifts by the
*          calling subroutine.
*
*     ND      (output) integer
*          The number of converged eigenvalues uncovered by this
*          subroutine.
*
*     SR      (output) DOUBLE PRECISION array, dimension KBOT
*     SI      (output) DOUBLE PRECISION array, dimension KBOT
*          On output, the real and imaginary parts of approximate
*          eigenvalues that may be used for shifts are stored in
*          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
*          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
*          The real and imaginary parts of converged eigenvalues
*          are stored in SR(KBOT-ND+1) through SR(KBOT) and
*          SI(KBOT-ND+1) through SI(KBOT), respectively.
*
*     V       (workspace) DOUBLE PRECISION array, dimension (LDV,NW)
*          An NW-by-NW work array.
*
*     LDV     (input) integer scalar
*          The leading dimension of V just as declared in the
*          calling subroutine.  NW .LE. LDV
*
*     NH      (input) integer scalar
*          The number of columns of T.  NH.GE.NW.
*
*     T       (workspace) DOUBLE PRECISION array, dimension (LDT,NW)
*
*     LDT     (input) integer
*          The leading dimension of T just as declared in the
*          calling subroutine.  NW .LE. LDT
*
*     NV      (input) integer
*          The number of rows of work array WV available for
*          workspace.  NV.GE.NW.
*
*     WV      (workspace) DOUBLE PRECISION array, dimension (LDWV,NW)
*
*     LDWV    (input) integer
*          The leading dimension of W just as declared in the
*          calling subroutine.  NW .LE. LDV
*
*     WORK    (workspace) DOUBLE PRECISION array, dimension LWORK.
*          On exit, WORK(1) is set to an estimate of the optimal value
*          of LWORK for the given values of N, NW, KTOP and KBOT.
*
*     LWORK   (input) integer
*          The dimension of the work array WORK.  LWORK = 2*NW
*          suffices, but greater efficiency may result from larger
*          values of LWORK.
*
*          If LWORK = -1, then a workspace query is assumed; DLAQR3
*          only estimates the optimal workspace size for the given
*          values of N, NW, KTOP and KBOT.  The estimate is returned
*          in WORK(1).  No error message related to LWORK is issued
*          by XERBLA.  Neither H nor Z are accessed.
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S,
     $                   SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL,
     $                   KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3,
     $                   LWKOPT, NMIN
      LOGICAL            BULGE, SORTED
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            ILAENV
      EXTERNAL           DLAMCH, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEHRD, DGEMM, DLABAD, DLACPY, DLAHQR,
     $                   DLANV2, DLAQR4, DLARF, DLARFG, DLASET, DORMHR,
     $                   DTREXC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      IF( JW.LE.2 ) THEN
         LWKOPT = 1
      ELSE
*
*        ==== Workspace query call to DGEHRD ====
*
         CALL DGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
         LWK1 = INT( WORK( 1 ) )
*
*        ==== Workspace query call to DORMHR ====
*
         CALL DORMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV,
     $                WORK, -1, INFO )
         LWK2 = INT( WORK( 1 ) )
*
*        ==== Workspace query call to DLAQR4 ====
*
         CALL DLAQR4( .true., .true., JW, 1, JW, T, LDT, SR, SI, 1, JW,
     $                V, LDV, WORK, -1, INFQR )
         LWK3 = INT( WORK( 1 ) )
*
*        ==== Optimal workspace ====
*
         LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( LWORK.EQ.-1 ) THEN
         WORK( 1 ) = DBLE( LWKOPT )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      NS = 0
      ND = 0
      WORK( 1 ) = ONE
      IF( KTOP.GT.KBOT )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( NW.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     ==== Setup deflation window ====
*
      JW = MIN( NW, KBOT-KTOP+1 )
      KWTOP = KBOT - JW + 1
      IF( KWTOP.EQ.KTOP ) THEN
         S = ZERO
      ELSE
         S = H( KWTOP, KWTOP-1 )
      END IF
*
      IF( KBOT.EQ.KWTOP ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         SR( KWTOP ) = H( KWTOP, KWTOP )
         SI( KWTOP ) = ZERO
         NS = 1
         ND = 0
         IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) )
     $        THEN
            NS = 0
            ND = 1
            IF( KWTOP.GT.KTOP )
     $         H( KWTOP, KWTOP-1 ) = ZERO
         END IF
         WORK( 1 ) = ONE
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL DLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
      CALL DCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
*
      CALL DLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
      NMIN = ILAENV( 12, 'DLAQR3', 'SV', JW, 1, JW, LWORK )
      IF( JW.GT.NMIN ) THEN
         CALL DLAQR4( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ),
     $                SI( KWTOP ), 1, JW, V, LDV, WORK, LWORK, INFQR )
      ELSE
         CALL DLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ),
     $                SI( KWTOP ), 1, JW, V, LDV, INFQR )
      END IF
*
*     ==== DTREXC needs a clean margin near the diagonal ====
*
      DO 10 J = 1, JW - 3
         T( J+2, J ) = ZERO
         T( J+3, J ) = ZERO
   10 CONTINUE
      IF( JW.GT.2 )
     $   T( JW, JW-2 ) = ZERO
*
*     ==== Deflation detection loop ====
*
      NS = JW
      ILST = INFQR + 1
   20 CONTINUE
      IF( ILST.LE.NS ) THEN
         IF( NS.EQ.1 ) THEN
            BULGE = .FALSE.
         ELSE
            BULGE = T( NS, NS-1 ).NE.ZERO
         END IF
*
*        ==== Small spike tip test for deflation ====
*
         IF( .NOT.BULGE ) THEN
*
*           ==== Real eigenvalue ====
*
            FOO = ABS( T( NS, NS ) )
            IF( FOO.EQ.ZERO )
     $         FOO = ABS( S )
            IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
*
*              ==== Deflatable ====
*
               NS = NS - 1
            ELSE
*
*              ==== Undeflatable.   Move it up out of the way.
*              .    (DTREXC can not fail in this case.) ====
*
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK,
     $                      INFO )
               ILST = ILST + 1
            END IF
         ELSE
*
*           ==== Complex conjugate pair ====
*
            FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )*
     $            SQRT( ABS( T( NS-1, NS ) ) )
            IF( FOO.EQ.ZERO )
     $         FOO = ABS( S )
            IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE.
     $          MAX( SMLNUM, ULP*FOO ) ) THEN
*
*              ==== Deflatable ====
*
               NS = NS - 2
            ELSE
*
*              ==== Undeflatable. Move them up out of the way.
*              .    Fortunately, DTREXC does the right thing with
*              .    ILST in case of a rare exchange failure. ====
*
               IFST = NS
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK,
     $                      INFO )
               ILST = ILST + 2
            END IF
         END IF
*
*        ==== End deflation detection loop ====
*
         GO TO 20
      END IF
*
*        ==== Return to Hessenberg form ====
*
      IF( NS.EQ.0 )
     $   S = ZERO
*
      IF( NS.LT.JW ) THEN
*
*        ==== sorting diagonal blocks of T improves accuracy for
*        .    graded matrices.  Bubble sort deals well with
*        .    exchange failures. ====
*
         SORTED = .false.
         I = NS + 1
   30    CONTINUE
         IF( SORTED )
     $      GO TO 50
         SORTED = .true.
*
         KEND = I - 1
         I = INFQR + 1
         IF( I.EQ.NS ) THEN
            K = I + 1
         ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
            K = I + 1
         ELSE
            K = I + 2
         END IF
   40    CONTINUE
         IF( K.LE.KEND ) THEN
            IF( K.EQ.I+1 ) THEN
               EVI = ABS( T( I, I ) )
            ELSE
               EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )*
     $               SQRT( ABS( T( I, I+1 ) ) )
            END IF
*
            IF( K.EQ.KEND ) THEN
               EVK = ABS( T( K, K ) )
            ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
               EVK = ABS( T( K, K ) )
            ELSE
               EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )*
     $               SQRT( ABS( T( K, K+1 ) ) )
            END IF
*
            IF( EVI.GE.EVK ) THEN
               I = K
            ELSE
               SORTED = .false.
               IFST = I
               ILST = K
               CALL DTREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK,
     $                      INFO )
               IF( INFO.EQ.0 ) THEN
                  I = ILST
               ELSE
                  I = K
               END IF
            END IF
            IF( I.EQ.KEND ) THEN
               K = I + 1
            ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
               K = I + 1
            ELSE
               K = I + 2
            END IF
            GO TO 40
         END IF
         GO TO 30
   50    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      I = JW
   60 CONTINUE
      IF( I.GE.INFQR+1 ) THEN
         IF( I.EQ.INFQR+1 ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
            SR( KWTOP+I-1 ) = T( I, I )
            SI( KWTOP+I-1 ) = ZERO
            I = I - 1
         ELSE
            AA = T( I-1, I-1 )
            CC = T( I, I-1 )
            BB = T( I-1, I )
            DD = T( I, I )
            CALL DLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ),
     $                   SI( KWTOP+I-2 ), SR( KWTOP+I-1 ),
     $                   SI( KWTOP+I-1 ), CS, SN )
            I = I - 2
         END IF
         GO TO 60
      END IF
*
      IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
         IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL DCOPY( NS, V, LDV, WORK, 1 )
            BETA = WORK( 1 )
            CALL DLARFG( NS, BETA, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
*
            CALL DLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
*
            CALL DLARF( 'L', NS, JW, WORK, 1, TAU, T, LDT,
     $                  WORK( JW+1 ) )
            CALL DLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT,
     $                  WORK( JW+1 ) )
            CALL DLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV,
     $                  WORK( JW+1 ) )
*
            CALL DGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ),
     $                   LWORK-JW, INFO )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( KWTOP.GT.1 )
     $      H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
         CALL DLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
         CALL DCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ),
     $               LDH+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( NS.GT.1 .AND. S.NE.ZERO )
     $      CALL DORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV,
     $                   WORK( JW+1 ), LWORK-JW, INFO )
*
*        ==== Update vertical slab in H ====
*
         IF( WANTT ) THEN
            LTOP = 1
         ELSE
            LTOP = KTOP
         END IF
         DO 70 KROW = LTOP, KWTOP - 1, NV
            KLN = MIN( NV, KWTOP-KROW )
            CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ),
     $                  LDH, V, LDV, ZERO, WV, LDWV )
            CALL DLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
   70    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( WANTT ) THEN
            DO 80 KCOL = KBOT + 1, N, NH
               KLN = MIN( NH, N-KCOL+1 )
               CALL DGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV,
     $                     H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
               CALL DLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ),
     $                      LDH )
   80       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( WANTZ ) THEN
            DO 90 KROW = ILOZ, IHIZ, NV
               KLN = MIN( NV, IHIZ-KROW+1 )
               CALL DGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ),
     $                     LDZ, V, LDV, ZERO, WV, LDWV )
               CALL DLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ),
     $                      LDZ )
   90       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      ND = JW - NS
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      NS = NS - INFQR
*
*      ==== Return optimal workspace. ====
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     ==== End of DLAQR3 ====
*
      END
      SUBROUTINE DLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
     $                   Z( LDZ, * )
*     ..
*
*     This subroutine implements one level of recursion for DLAQR0.
*     It is a complete implementation of the small bulge multi-shift
*     QR algorithm.  It may be called by DLAQR0 and, for large enough
*     deflation window size, it may be called by DLAQR3.  This
*     subroutine is identical to DLAQR0 except that it calls DLAQR2
*     instead of DLAQR3.
*
*     Purpose
*     =======
*
*     DLAQR4 computes the eigenvalues of a Hessenberg matrix H
*     and, optionally, the matrices T and Z from the Schur decomposition
*     H = Z T Z**T, where T is an upper quasi-triangular matrix (the
*     Schur form), and Z is the orthogonal matrix of Schur vectors.
*
*     Optionally Z may be postmultiplied into an input orthogonal
*     matrix Q so that this routine can give the Schur factorization
*     of a matrix A which has been reduced to the Hessenberg form H
*     by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
*
*     Arguments
*     =========
*
*     WANTT   (input) LOGICAL
*          = .TRUE. : the full Schur form T is required;
*          = .FALSE.: only eigenvalues are required.
*
*     WANTZ   (input) LOGICAL
*          = .TRUE. : the matrix of Schur vectors Z is required;
*          = .FALSE.: Schur vectors are not required.
*
*     N     (input) INTEGER
*           The order of the matrix H.  N .GE. 0.
*
*     ILO   (input) INTEGER
*     IHI   (input) INTEGER
*           It is assumed that H is already upper triangular in rows
*           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
*           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
*           previous call to DGEBAL, and then passed to DGEHRD when the
*           matrix output by DGEBAL is reduced to Hessenberg form.
*           Otherwise, ILO and IHI should be set to 1 and N,
*           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*           If N = 0, then ILO = 1 and IHI = 0.
*
*     H     (input/output) DOUBLE PRECISION array, dimension (LDH,N)
*           On entry, the upper Hessenberg matrix H.
*           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
*           the upper quasi-triangular matrix T from the Schur
*           decomposition (the Schur form); 2-by-2 diagonal blocks
*           (corresponding to complex conjugate pairs of eigenvalues)
*           are returned in standard form, with H(i,i) = H(i+1,i+1)
*           and H(i+1,i)*H(i,i+1).LT.0. If INFO = 0 and WANTT is
*           .FALSE., then the contents of H are unspecified on exit.
*           (The output value of H when INFO.GT.0 is given under the
*           description of INFO below.)
*
*           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
*           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
*
*     LDH   (input) INTEGER
*           The leading dimension of the array H. LDH .GE. max(1,N).
*
*     WR    (output) DOUBLE PRECISION array, dimension (IHI)
*     WI    (output) DOUBLE PRECISION array, dimension (IHI)
*           The real and imaginary parts, respectively, of the computed
*           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
*           and WI(ILO:IHI). If two eigenvalues are computed as a
*           complex conjugate pair, they are stored in consecutive
*           elements of WR and WI, say the i-th and (i+1)th, with
*           WI(i) .GT. 0 and WI(i+1) .LT. 0. If WANTT is .TRUE., then
*           the eigenvalues are stored in the same order as on the
*           diagonal of the Schur form returned in H, with
*           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
*           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
*           WI(i+1) = -WI(i).
*
*     ILOZ     (input) INTEGER
*     IHIZ     (input) INTEGER
*           Specify the rows of Z to which transformations must be
*           applied if WANTZ is .TRUE..
*           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
*
*     Z     (input/output) DOUBLE PRECISION array, dimension (LDZ,IHI)
*           If WANTZ is .FALSE., then Z is not referenced.
*           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
*           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
*           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
*           (The output value of Z when INFO.GT.0 is given under
*           the description of INFO below.)
*
*     LDZ   (input) INTEGER
*           The leading dimension of the array Z.  if WANTZ is .TRUE.
*           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
*
*     WORK  (workspace/output) DOUBLE PRECISION array, dimension LWORK
*           On exit, if LWORK = -1, WORK(1) returns an estimate of
*           the optimal value for LWORK.
*
*     LWORK (input) INTEGER
*           The dimension of the array WORK.  LWORK .GE. max(1,N)
*           is sufficient, but LWORK typically as large as 6*N may
*           be required for optimal performance.  A workspace query
*           to determine the optimal workspace size is recommended.
*
*           If LWORK = -1, then DLAQR4 does a workspace query.
*           In this case, DLAQR4 checks the input parameters and
*           estimates the optimal workspace size for the given
*           values of N, ILO and IHI.  The estimate is returned
*           in WORK(1).  No error message related to LWORK is
*           issued by XERBLA.  Neither H nor Z are accessed.
*
*
*     INFO  (output) INTEGER
*             =  0:  successful exit
*           .GT. 0:  if INFO = i, DLAQR4 failed to compute all of
*                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*                and WI contain those eigenvalues which have been
*                successfully computed.  (Failures are rare.)
*
*                If INFO .GT. 0 and WANT is .FALSE., then on exit,
*                the remaining unconverged eigenvalues are the eigen-
*                values of the upper Hessenberg matrix rows and
*                columns ILO through INFO of the final, output
*                value of H.
*
*                If INFO .GT. 0 and WANTT is .TRUE., then on exit
*
*           (*)  (initial value of H)*U  = U*(final value of H)
*
*                where U is an orthogonal matrix.  The final
*                value of H is upper Hessenberg and quasi-triangular
*                in rows and columns INFO+1 through IHI.
*
*                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
*
*                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
*                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
*
*                where U is the orthogonal matrix in (*) (regard-
*                less of the value of WANTT.)
*
*                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
*                accessed.
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*     References:
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*       929--947, 2002.
*
*       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*       of Matrix Analysis, volume 23, pages 948--973, 2002.
*
*     ================================================================
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    DLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      PARAMETER          ( NTINY = 11 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      PARAMETER          ( KEXNW = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      PARAMETER          ( KEXSH = 6 )
*
*     ==== The constants WILK1 and WILK2 are used to form the
*     .    exceptional shifts. ====
      DOUBLE PRECISION   WILK1, WILK2
      PARAMETER          ( WILK1 = 0.75d0, WILK2 = -0.4375d0 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS, DD, SN, SS, SWAP
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS,
     $                   LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS,
     $                   NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLAHQR, DLANV2, DLAQR2, DLAQR5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
      INFO = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = ONE
         RETURN
      END IF
*
      IF( N.LE.NTINY ) THEN
*
*        ==== Tiny matrices must use DLAHQR. ====
*
         LWKOPT = 1
         IF( LWORK.NE.-1 )
     $      CALL DLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
     $                   ILOZ, IHIZ, Z, LDZ, INFO )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         INFO = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( WANTT ) THEN
            JBCMPZ( 1: 1 ) = 'S'
         ELSE
            JBCMPZ( 1: 1 ) = 'E'
         END IF
         IF( WANTZ ) THEN
            JBCMPZ( 2: 2 ) = 'V'
         ELSE
            JBCMPZ( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 11, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.3.) ====
*
         NWR = ILAENV( 13, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NWR = MAX( 2, NWR )
         NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 11, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         NSR = ILAENV( 15, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
         NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to DLAQR2 ====
*
         CALL DLAQR2( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ,
     $                IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH,
     $                N, H, LDH, WORK, -1 )
*
*        ==== Optimal workspace = MAX(DLAQR5, DLAQR2) ====
*
         LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = DBLE( LWKOPT )
            RETURN
         END IF
*
*        ==== DLAHQR/DLAQR0 crossover point ====
*
         NMIN = ILAENV( 12, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NMIN = MAX( NTINY, NMIN )
*
*        ==== Nibble crossover point ====
*
         NIBBLE = ILAENV( 14, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         NIBBLE = MAX( 0, NIBBLE )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         KACC22 = ILAENV( 16, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK )
         KACC22 = MAX( 0, KACC22 )
         KACC22 = MIN( 2, KACC22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
         NW = NWMAX
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
         NSMAX = NSMAX - MOD( NSMAX, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         NDFL = 1
*
*        ==== ITMAX = iteration limit ====
*
         ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
*
*        ==== Last row and column in the active block ====
*
         KBOT = IHI
*
*        ==== Main Loop ====
*
         DO 80 IT = 1, ITMAX
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( KBOT.LT.ILO )
     $         GO TO 90
*
*           ==== Locate active block ====
*
            DO 10 K = KBOT, ILO + 1, -1
               IF( H( K, K-1 ).EQ.ZERO )
     $            GO TO 20
   10       CONTINUE
            K = ILO
   20       CONTINUE
            KTOP = K
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            NH = KBOT - KTOP + 1
            NWUPBD = MIN( NH, NWMAX )
            IF( NDFL.LT.KEXNW ) THEN
               NW = MIN( NWUPBD, NWR )
            ELSE
               NW = MIN( NWUPBD, 2*NW )
            END IF
            IF( NW.LT.NWMAX ) THEN
               IF( NW.GE.NH-1 ) THEN
                  NW = NH
               ELSE
                  KWTOP = KBOT - NW + 1
                  IF( ABS( H( KWTOP, KWTOP-1 ) ).GT.
     $                ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
               END IF
            END IF
            IF( NDFL.LT.KEXNW ) THEN
               NDEC = -1
            ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
               NDEC = NDEC + 1
               IF( NW-NDEC.LT.2 )
     $            NDEC = 0
               NW = NW - NDEC
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            KV = N - NW + 1
            KT = NW + 1
            NHO = ( N-NW-1 ) - KT + 1
            KWV = NW + 2
            NVE = ( N-NW ) - KWV + 1
*
*           ==== Aggressive early deflation ====
*
            CALL DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH,
     $                   NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH,
     $                   WORK, LWORK )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            KBOT = KBOT - LD
*
*           ==== KS points to the shifts. ====
*
            KS = KBOT - LS + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT-
     $          KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if DLAQR2
*              .    did not provide that many shifts. ====
*
               NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
               NS = NS - MOD( NS, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    DLAQR2 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                  KS = KBOT - NS + 1
                  DO 30 I = KBOT, MAX( KS+1, KTOP+2 ), -2
                     SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                     AA = WILK1*SS + H( I, I )
                     BB = SS
                     CC = WILK2*SS
                     DD = AA
                     CALL DLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ),
     $                            WR( I ), WI( I ), CS, SN )
   30             CONTINUE
                  IF( KS.EQ.KTOP ) THEN
                     WR( KS+1 ) = H( KS+1, KS+1 )
                     WI( KS+1 ) = ZERO
                     WR( KS ) = WR( KS+1 )
                     WI( KS ) = WI( KS+1 )
                  END IF
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use DLAHQR
*                 .    on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( KBOT-KS+1.LE.NS / 2 ) THEN
                     KS = KBOT - NS + 1
                     KT = N - NS + 1
                     CALL DLACPY( 'A', NS, NS, H( KS, KS ), LDH,
     $                            H( KT, 1 ), LDH )
                     CALL DLAHQR( .false., .false., NS, 1, NS,
     $                            H( KT, 1 ), LDH, WR( KS ), WI( KS ),
     $                            1, 1, ZDUM, 1, INF )
                     KS = KS + INF
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  ====
*
                     IF( KS.GE.KBOT ) THEN
                        AA = H( KBOT-1, KBOT-1 )
                        CC = H( KBOT, KBOT-1 )
                        BB = H( KBOT-1, KBOT )
                        DD = H( KBOT, KBOT )
                        CALL DLANV2( AA, BB, CC, DD, WR( KBOT-1 ),
     $                               WI( KBOT-1 ), WR( KBOT ),
     $                               WI( KBOT ), CS, SN )
                        KS = KBOT - 1
                     END IF
                  END IF
*
                  IF( KBOT-KS+1.GT.NS ) THEN
*
*                    ==== Sort the shifts (Helps a little)
*                    .    Bubble sort keeps complex conjugate
*                    .    pairs together. ====
*
                     SORTED = .false.
                     DO 50 K = KBOT, KS + 1, -1
                        IF( SORTED )
     $                     GO TO 60
                        SORTED = .true.
                        DO 40 I = KS, K - 1
                           IF( ABS( WR( I ) )+ABS( WI( I ) ).LT.
     $                         ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                              SORTED = .false.
*
                              SWAP = WR( I )
                              WR( I ) = WR( I+1 )
                              WR( I+1 ) = SWAP
*
                              SWAP = WI( I )
                              WI( I ) = WI( I+1 )
                              WI( I+1 ) = SWAP
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
*
*                 ==== Shuffle shifts into pairs of real shifts
*                 .    and pairs of complex conjugate shifts
*                 .    assuming complex conjugate shifts are
*                 .    already adjacent to one another. (Yes,
*                 .    they are.)  ====
*
                  DO 70 I = KBOT, KS + 2, -2
                     IF( WI( I ).NE.-WI( I-1 ) ) THEN
*
                        SWAP = WR( I )
                        WR( I ) = WR( I-1 )
                        WR( I-1 ) = WR( I-2 )
                        WR( I-2 ) = SWAP
*
                        SWAP = WI( I )
                        WI( I ) = WI( I-1 )
                        WI( I-1 ) = WI( I-2 )
                        WI( I-2 ) = SWAP
                     END IF
   70             CONTINUE
               END IF
*
*              ==== If there are only two shifts and both are
*              .    real, then use only one.  ====
*
               IF( KBOT-KS+1.EQ.2 ) THEN
                  IF( WI( KBOT ).EQ.ZERO ) THEN
                     IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT.
     $                   ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                        WR( KBOT-1 ) = WR( KBOT )
                     ELSE
                        WR( KBOT ) = WR( KBOT-1 )
                     END IF
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnatiude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               NS = MIN( NS, KBOT-KS+1 )
               NS = NS - MOD( NS, 2 )
               KS = KBOT - NS + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               KDU = 3*NS - 3
               KU = N - KDU + 1
               KWH = KDU + 1
               NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
               KWV = KDU + 4
               NVE = N - KDU - KWV + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS,
     $                      WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z,
     $                      LDZ, WORK, 3, H( KU, 1 ), LDH, NVE,
     $                      H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( LD.GT.0 ) THEN
               NDFL = 1
            ELSE
               NDFL = NDFL + 1
            END IF
*
*           ==== End of main loop ====
   80    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         INFO = KBOT
   90    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      WORK( 1 ) = DBLE( LWKOPT )
*
*     ==== End of DLAQR4 ====
*
      END
      SUBROUTINE DLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS,
     $                   SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U,
     $                   LDU, NV, WV, LDWV, NH, WH, LDWH )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
     $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), U( LDU, * ),
     $                   V( LDV, * ), WH( LDWH, * ), WV( LDWV, * ),
     $                   Z( LDZ, * )
*     ..
*
*     This auxiliary subroutine called by DLAQR0 performs a
*     single small-bulge multi-shift QR sweep.
*
*      WANTT  (input) logical scalar
*             WANTT = .true. if the quasi-triangular Schur factor
*             is being computed.  WANTT is set to .false. otherwise.
*
*      WANTZ  (input) logical scalar
*             WANTZ = .true. if the orthogonal Schur factor is being
*             computed.  WANTZ is set to .false. otherwise.
*
*      KACC22 (input) integer with value 0, 1, or 2.
*             Specifies the computation mode of far-from-diagonal
*             orthogonal updates.
*        = 0: DLAQR5 does not accumulate reflections and does not
*             use matrix-matrix multiply to update far-from-diagonal
*             matrix entries.
*        = 1: DLAQR5 accumulates reflections and uses matrix-matrix
*             multiply to update the far-from-diagonal matrix entries.
*        = 2: DLAQR5 accumulates reflections, uses matrix-matrix
*             multiply to update the far-from-diagonal matrix entries,
*             and takes advantage of 2-by-2 block structure during
*             matrix multiplies.
*
*      N      (input) integer scalar
*             N is the order of the Hessenberg matrix H upon which this
*             subroutine operates.
*
*      KTOP   (input) integer scalar
*      KBOT   (input) integer scalar
*             These are the first and last rows and columns of an
*             isolated diagonal block upon which the QR sweep is to be
*             applied. It is assumed without a check that
*                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
*             and
*                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
*
*      NSHFTS (input) integer scalar
*             NSHFTS gives the number of simultaneous shifts.  NSHFTS
*             must be positive and even.
*
*      SR     (input/output) DOUBLE PRECISION array of size (NSHFTS)
*      SI     (input/output) DOUBLE PRECISION array of size (NSHFTS)
*             SR contains the real parts and SI contains the imaginary
*             parts of the NSHFTS shifts of origin that define the
*             multi-shift QR sweep.  On output SR and SI may be
*             reordered.
*
*      H      (input/output) DOUBLE PRECISION array of size (LDH,N)
*             On input H contains a Hessenberg matrix.  On output a
*             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
*             to the isolated diagonal block in rows and columns KTOP
*             through KBOT.
*
*      LDH    (input) integer scalar
*             LDH is the leading dimension of H just as declared in the
*             calling procedure.  LDH.GE.MAX(1,N).
*
*      ILOZ   (input) INTEGER
*      IHIZ   (input) INTEGER
*             Specify the rows of Z to which transformations must be
*             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
*
*      Z      (input/output) DOUBLE PRECISION array of size (LDZ,IHI)
*             If WANTZ = .TRUE., then the QR Sweep orthogonal
*             similarity transformation is accumulated into
*             Z(ILOZ:IHIZ,ILO:IHI) from the right.
*             If WANTZ = .FALSE., then Z is unreferenced.
*
*      LDZ    (input) integer scalar
*             LDA is the leading dimension of Z just as declared in
*             the calling procedure. LDZ.GE.N.
*
*      V      (workspace) DOUBLE PRECISION array of size (LDV,NSHFTS/2)
*
*      LDV    (input) integer scalar
*             LDV is the leading dimension of V as declared in the
*             calling procedure.  LDV.GE.3.
*
*      U      (workspace) DOUBLE PRECISION array of size
*             (LDU,3*NSHFTS-3)
*
*      LDU    (input) integer scalar
*             LDU is the leading dimension of U just as declared in the
*             in the calling subroutine.  LDU.GE.3*NSHFTS-3.
*
*      NH     (input) integer scalar
*             NH is the number of columns in array WH available for
*             workspace. NH.GE.1.
*
*      WH     (workspace) DOUBLE PRECISION array of size (LDWH,NH)
*
*      LDWH   (input) integer scalar
*             Leading dimension of WH just as declared in the
*             calling procedure.  LDWH.GE.3*NSHFTS-3.
*
*      NV     (input) integer scalar
*             NV is the number of rows in WV agailable for workspace.
*             NV.GE.1.
*
*      WV     (workspace) DOUBLE PRECISION array of size
*             (LDWV,3*NSHFTS-3)
*
*      LDWV   (input) integer scalar
*             LDWV is the leading dimension of WV as declared in the
*             in the calling subroutine.  LDWV.GE.NV.
*
*     ================================================================
*     Based on contributions by
*        Karen Braman and Ralph Byers, Department of Mathematics,
*        University of Kansas, USA
*
*     ================================================================
*     Reference:
*
*     K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*     Algorithm Part I: Maintaining Well Focused Shifts, and
*     Level 3 Performance, SIAM Journal of Matrix Analysis,
*     volume 23, pages 929--947, 2002.
*
*     ================================================================
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0d0, ONE = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, BETA, H11, H12, H21, H22, REFSUM,
     $                   SAFMAX, SAFMIN, SCL, SMLNUM, SWAP, TST1, TST2,
     $                   ULP
      INTEGER            I, I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN,
     $                   JROW, JTOP, K, K1, KDU, KMS, KNZ, KRCOL, KZS,
     $                   M, M22, MBOT, MEND, MSTART, MTOP, NBMPS, NDCOL,
     $                   NS, NU
      LOGICAL            ACCUM, BLK22, BMP22
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          ABS, DBLE, MAX, MIN, MOD
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   VT( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DLABAD, DLACPY, DLAQR1, DLARFG, DLASET,
     $                   DTRMM
*     ..
*     .. Executable Statements ..
*
*     ==== If there are no shifts, then there is nothing to do. ====
*
      IF( NSHFTS.LT.2 )
     $   RETURN
*
*     ==== If the active block is empty or 1-by-1, then there
*     .    is nothing to do. ====
*
      IF( KTOP.GE.KBOT )
     $   RETURN
*
*     ==== Shuffle shifts into pairs of real shifts and pairs
*     .    of complex conjugate shifts assuming complex
*     .    conjugate shifts are already adjacent to one
*     .    another. ====
*
      DO 10 I = 1, NSHFTS - 2, 2
         IF( SI( I ).NE.-SI( I+1 ) ) THEN
*
            SWAP = SR( I )
            SR( I ) = SR( I+1 )
            SR( I+1 ) = SR( I+2 )
            SR( I+2 ) = SWAP
*
            SWAP = SI( I )
            SI( I ) = SI( I+1 )
            SI( I+1 ) = SI( I+2 )
            SI( I+2 ) = SWAP
         END IF
   10 CONTINUE
*
*     ==== NSHFTS is supposed to be even, but if it is odd,
*     .    then simply reduce it by one.  The shuffle above
*     .    ensures that the dropped shift is real and that
*     .    the remaining shifts are paired. ====
*
      NS = NSHFTS - MOD( NSHFTS, 2 )
*
*     ==== Machine constants for deflation ====
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N ) / ULP )
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries ? ====
*
      ACCUM = ( KACC22.EQ.1 ) .OR. ( KACC22.EQ.2 )
*
*     ==== If so, exploit the 2-by-2 block structure? ====
*
      BLK22 = ( NS.GT.2 ) .AND. ( KACC22.EQ.2 )
*
*     ==== clear trash ====
*
      IF( KTOP+2.LE.KBOT )
     $   H( KTOP+2, KTOP ) = ZERO
*
*     ==== NBMPS = number of 2-shift bulges in the chain ====
*
      NBMPS = NS / 2
*
*     ==== KDU = width of slab ====
*
      KDU = 6*NBMPS - 3
*
*     ==== Create and chase chains of NBMPS bulges ====
*
      DO 220 INCOL = 3*( 1-NBMPS ) + KTOP - 1, KBOT - 2, 3*NBMPS - 2
         NDCOL = INCOL + KDU
         IF( ACCUM )
     $      CALL DLASET( 'ALL', KDU, KDU, ZERO, ONE, U, LDU )
*
*        ==== Near-the-diagonal bulge chase.  The following loop
*        .    performs the near-the-diagonal part of a small bulge
*        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
*        .    chunk extends from column INCOL to column NDCOL
*        .    (including both column INCOL and column NDCOL). The
*        .    following loop chases a 3*NBMPS column long chain of
*        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
*        .    may be less than KTOP and and NDCOL may be greater than
*        .    KBOT indicating phantom columns from which to chase
*        .    bulges before they are actually introduced or to which
*        .    to chase bulges beyond column KBOT.)  ====
*
         DO 150 KRCOL = INCOL, MIN( INCOL+3*NBMPS-3, KBOT-2 )
*
*           ==== Bulges number MTOP to MBOT are active double implicit
*           .    shift bulges.  There may or may not also be small
*           .    2-by-2 bulge, if there is room.  The inactive bulges
*           .    (if any) must wait until the active bulges have moved
*           .    down the diagonal to make room.  The phantom matrix
*           .    paradigm described above helps keep track.  ====
*
            MTOP = MAX( 1, ( ( KTOP-1 )-KRCOL+2 ) / 3+1 )
            MBOT = MIN( NBMPS, ( KBOT-KRCOL ) / 3 )
            M22 = MBOT + 1
            BMP22 = ( MBOT.LT.NBMPS ) .AND. ( KRCOL+3*( M22-1 ) ).EQ.
     $              ( KBOT-2 )
*
*           ==== Generate reflections to chase the chain right
*           .    one column.  (The minimum value of K is KTOP-1.) ====
*
            DO 20 M = MTOP, MBOT
               K = KRCOL + 3*( M-1 )
               IF( K.EQ.KTOP-1 ) THEN
                  CALL DLAQR1( 3, H( KTOP, KTOP ), LDH, SR( 2*M-1 ),
     $                         SI( 2*M-1 ), SR( 2*M ), SI( 2*M ),
     $                         V( 1, M ) )
                  ALPHA = V( 1, M )
                  CALL DLARFG( 3, ALPHA, V( 2, M ), 1, V( 1, M ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M ) = H( K+2, K )
                  V( 3, M ) = H( K+3, K )
                  CALL DLARFG( 3, BETA, V( 2, M ), 1, V( 1, M ) )
*
*                 ==== A Bulge may collapse because of vigilant
*                 .    deflation or destructive underflow.  In the
*                 .    underflow case, try the two-small-subdiagonals
*                 .    trick to try to reinflate the bulge.  ====
*
                  IF( H( K+3, K ).NE.ZERO .OR. H( K+3, K+1 ).NE.
     $                ZERO .OR. H( K+3, K+2 ).EQ.ZERO ) THEN
*
*                    ==== Typical case: not collapsed (yet). ====
*
                     H( K+1, K ) = BETA
                     H( K+2, K ) = ZERO
                     H( K+3, K ) = ZERO
                  ELSE
*
*                    ==== Atypical case: collapsed.  Attempt to
*                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
*                    .    If the fill resulting from the new
*                    .    reflector is too large, then abandon it.
*                    .    Otherwise, use the new one. ====
*
                     CALL DLAQR1( 3, H( K+1, K+1 ), LDH, SR( 2*M-1 ),
     $                            SI( 2*M-1 ), SR( 2*M ), SI( 2*M ),
     $                            VT )
                     ALPHA = VT( 1 )
                     CALL DLARFG( 3, ALPHA, VT( 2 ), 1, VT( 1 ) )
                     REFSUM = VT( 1 )*( H( K+1, K )+VT( 2 )*
     $                        H( K+2, K ) )
*
                     IF( ABS( H( K+2, K )-REFSUM*VT( 2 ) )+
     $                   ABS( REFSUM*VT( 3 ) ).GT.ULP*
     $                   ( ABS( H( K, K ) )+ABS( H( K+1,
     $                   K+1 ) )+ABS( H( K+2, K+2 ) ) ) ) THEN
*
*                       ==== Starting a new bulge here would
*                       .    create non-negligible fill.  Use
*                       .    the old one with trepidation. ====
*
                        H( K+1, K ) = BETA
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                     ELSE
*
*                       ==== Stating a new bulge here would
*                       .    create only negligible fill.
*                       .    Replace the old reflector with
*                       .    the new one. ====
*
                        H( K+1, K ) = H( K+1, K ) - REFSUM
                        H( K+2, K ) = ZERO
                        H( K+3, K ) = ZERO
                        V( 1, M ) = VT( 1 )
                        V( 2, M ) = VT( 2 )
                        V( 3, M ) = VT( 3 )
                     END IF
                  END IF
               END IF
   20       CONTINUE
*
*           ==== Generate a 2-by-2 reflection, if needed. ====
*
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 ) THEN
               IF( K.EQ.KTOP-1 ) THEN
                  CALL DLAQR1( 2, H( K+1, K+1 ), LDH, SR( 2*M22-1 ),
     $                         SI( 2*M22-1 ), SR( 2*M22 ), SI( 2*M22 ),
     $                         V( 1, M22 ) )
                  BETA = V( 1, M22 )
                  CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
               ELSE
                  BETA = H( K+1, K )
                  V( 2, M22 ) = H( K+2, K )
                  CALL DLARFG( 2, BETA, V( 2, M22 ), 1, V( 1, M22 ) )
                  H( K+1, K ) = BETA
                  H( K+2, K ) = ZERO
               END IF
            END IF
*
*           ==== Multiply H by reflections from the left ====
*
            IF( ACCUM ) THEN
               JBOT = MIN( NDCOL, KBOT )
            ELSE IF( WANTT ) THEN
               JBOT = N
            ELSE
               JBOT = KBOT
            END IF
            DO 40 J = MAX( KTOP, KRCOL ), JBOT
               MEND = MIN( MBOT, ( J-KRCOL+2 ) / 3 )
               DO 30 M = MTOP, MEND
                  K = KRCOL + 3*( M-1 )
                  REFSUM = V( 1, M )*( H( K+1, J )+V( 2, M )*
     $                     H( K+2, J )+V( 3, M )*H( K+3, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M )
                  H( K+3, J ) = H( K+3, J ) - REFSUM*V( 3, M )
   30          CONTINUE
   40       CONTINUE
            IF( BMP22 ) THEN
               K = KRCOL + 3*( M22-1 )
               DO 50 J = MAX( K+1, KTOP ), JBOT
                  REFSUM = V( 1, M22 )*( H( K+1, J )+V( 2, M22 )*
     $                     H( K+2, J ) )
                  H( K+1, J ) = H( K+1, J ) - REFSUM
                  H( K+2, J ) = H( K+2, J ) - REFSUM*V( 2, M22 )
   50          CONTINUE
            END IF
*
*           ==== Multiply H by reflections from the right.
*           .    Delay filling in the last row until the
*           .    vigilant deflation check is complete. ====
*
            IF( ACCUM ) THEN
               JTOP = MAX( KTOP, INCOL )
            ELSE IF( WANTT ) THEN
               JTOP = 1
            ELSE
               JTOP = KTOP
            END IF
            DO 90 M = MTOP, MBOT
               IF( V( 1, M ).NE.ZERO ) THEN
                  K = KRCOL + 3*( M-1 )
                  DO 60 J = JTOP, MIN( KBOT, K+3 )
                     REFSUM = V( 1, M )*( H( J, K+1 )+V( 2, M )*
     $                        H( J, K+2 )+V( 3, M )*H( J, K+3 ) )
                     H( J, K+1 ) = H( J, K+1 ) - REFSUM
                     H( J, K+2 ) = H( J, K+2 ) - REFSUM*V( 2, M )
                     H( J, K+3 ) = H( J, K+3 ) - REFSUM*V( 3, M )
   60             CONTINUE
*
                  IF( ACCUM ) THEN
*
*                    ==== Accumulate U. (If necessary, update Z later
*                    .    with with an efficient matrix-matrix
*                    .    multiply.) ====
*
                     KMS = K - INCOL
                     DO 70 J = MAX( 1, KTOP-INCOL ), KDU
                        REFSUM = V( 1, M )*( U( J, KMS+1 )+V( 2, M )*
     $                           U( J, KMS+2 )+V( 3, M )*U( J, KMS+3 ) )
                        U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                        U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*V( 2, M )
                        U( J, KMS+3 ) = U( J, KMS+3 ) - REFSUM*V( 3, M )
   70                CONTINUE
                  ELSE IF( WANTZ ) THEN
*
*                    ==== U is not accumulated, so update Z
*                    .    now by multiplying by reflections
*                    .    from the right. ====
*
                     DO 80 J = ILOZ, IHIZ
                        REFSUM = V( 1, M )*( Z( J, K+1 )+V( 2, M )*
     $                           Z( J, K+2 )+V( 3, M )*Z( J, K+3 ) )
                        Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                        Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*V( 2, M )
                        Z( J, K+3 ) = Z( J, K+3 ) - REFSUM*V( 3, M )
   80                CONTINUE
                  END IF
               END IF
   90       CONTINUE
*
*           ==== Special case: 2-by-2 reflection (if needed) ====
*
            K = KRCOL + 3*( M22-1 )
            IF( BMP22 .AND. ( V( 1, M22 ).NE.ZERO ) ) THEN
               DO 100 J = JTOP, MIN( KBOT, K+3 )
                  REFSUM = V( 1, M22 )*( H( J, K+1 )+V( 2, M22 )*
     $                     H( J, K+2 ) )
                  H( J, K+1 ) = H( J, K+1 ) - REFSUM
                  H( J, K+2 ) = H( J, K+2 ) - REFSUM*V( 2, M22 )
  100          CONTINUE
*
               IF( ACCUM ) THEN
                  KMS = K - INCOL
                  DO 110 J = MAX( 1, KTOP-INCOL ), KDU
                     REFSUM = V( 1, M22 )*( U( J, KMS+1 )+V( 2, M22 )*
     $                        U( J, KMS+2 ) )
                     U( J, KMS+1 ) = U( J, KMS+1 ) - REFSUM
                     U( J, KMS+2 ) = U( J, KMS+2 ) - REFSUM*V( 2, M22 )
  110             CONTINUE
               ELSE IF( WANTZ ) THEN
                  DO 120 J = ILOZ, IHIZ
                     REFSUM = V( 1, M22 )*( Z( J, K+1 )+V( 2, M22 )*
     $                        Z( J, K+2 ) )
                     Z( J, K+1 ) = Z( J, K+1 ) - REFSUM
                     Z( J, K+2 ) = Z( J, K+2 ) - REFSUM*V( 2, M22 )
  120             CONTINUE
               END IF
            END IF
*
*           ==== Vigilant deflation check ====
*
            MSTART = MTOP
            IF( KRCOL+3*( MSTART-1 ).LT.KTOP )
     $         MSTART = MSTART + 1
            MEND = MBOT
            IF( BMP22 )
     $         MEND = MEND + 1
            IF( KRCOL.EQ.KBOT-2 )
     $         MEND = MEND + 1
            DO 130 M = MSTART, MEND
               K = MIN( KBOT-1, KRCOL+3*( M-1 ) )
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( H( K+1, K ).NE.ZERO ) THEN
                  TST1 = ABS( H( K, K ) ) + ABS( H( K+1, K+1 ) )
                  IF( TST1.EQ.ZERO ) THEN
                     IF( K.GE.KTOP+1 )
     $                  TST1 = TST1 + ABS( H( K, K-1 ) )
                     IF( K.GE.KTOP+2 )
     $                  TST1 = TST1 + ABS( H( K, K-2 ) )
                     IF( K.GE.KTOP+3 )
     $                  TST1 = TST1 + ABS( H( K, K-3 ) )
                     IF( K.LE.KBOT-2 )
     $                  TST1 = TST1 + ABS( H( K+2, K+1 ) )
                     IF( K.LE.KBOT-3 )
     $                  TST1 = TST1 + ABS( H( K+3, K+1 ) )
                     IF( K.LE.KBOT-4 )
     $                  TST1 = TST1 + ABS( H( K+4, K+1 ) )
                  END IF
                  IF( ABS( H( K+1, K ) ).LE.MAX( SMLNUM, ULP*TST1 ) )
     $                 THEN
                     H12 = MAX( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
                     H21 = MIN( ABS( H( K+1, K ) ), ABS( H( K, K+1 ) ) )
                     H11 = MAX( ABS( H( K+1, K+1 ) ),
     $                     ABS( H( K, K )-H( K+1, K+1 ) ) )
                     H22 = MIN( ABS( H( K+1, K+1 ) ),
     $                     ABS( H( K, K )-H( K+1, K+1 ) ) )
                     SCL = H11 + H12
                     TST2 = H22*( H11 / SCL )
*
                     IF( TST2.EQ.ZERO .OR. H21*( H12 / SCL ).LE.
     $                   MAX( SMLNUM, ULP*TST2 ) )H( K+1, K ) = ZERO
                  END IF
               END IF
  130       CONTINUE
*
*           ==== Fill in the last row of each bulge. ====
*
            MEND = MIN( NBMPS, ( KBOT-KRCOL-1 ) / 3 )
            DO 140 M = MTOP, MEND
               K = KRCOL + 3*( M-1 )
               REFSUM = V( 1, M )*V( 3, M )*H( K+4, K+3 )
               H( K+4, K+1 ) = -REFSUM
               H( K+4, K+2 ) = -REFSUM*V( 2, M )
               H( K+4, K+3 ) = H( K+4, K+3 ) - REFSUM*V( 3, M )
  140       CONTINUE
*
*           ==== End of near-the-diagonal bulge chase. ====
*
  150    CONTINUE
*
*        ==== Use U (if accumulated) to update far-from-diagonal
*        .    entries in H.  If required, use U to update Z as
*        .    well. ====
*
         IF( ACCUM ) THEN
            IF( WANTT ) THEN
               JTOP = 1
               JBOT = N
            ELSE
               JTOP = KTOP
               JBOT = KBOT
            END IF
            IF( ( .NOT.BLK22 ) .OR. ( INCOL.LT.KTOP ) .OR.
     $          ( NDCOL.GT.KBOT ) .OR. ( NS.LE.2 ) ) THEN
*
*              ==== Updates not exploiting the 2-by-2 block
*              .    structure of U.  K1 and NU keep track of
*              .    the location and size of U in the special
*              .    cases of introducing bulges and chasing
*              .    bulges off the bottom.  In these special
*              .    cases and in case the number of shifts
*              .    is NS = 2, there is no 2-by-2 block
*              .    structure to exploit.  ====
*
               K1 = MAX( 1, KTOP-INCOL )
               NU = ( KDU-MAX( 0, NDCOL-KBOT ) ) - K1 + 1
*
*              ==== Horizontal Multiply ====
*
               DO 160 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
                  CALL DGEMM( 'C', 'N', NU, JLEN, NU, ONE, U( K1, K1 ),
     $                        LDU, H( INCOL+K1, JCOL ), LDH, ZERO, WH,
     $                        LDWH )
                  CALL DLACPY( 'ALL', NU, JLEN, WH, LDWH,
     $                         H( INCOL+K1, JCOL ), LDH )
  160          CONTINUE
*
*              ==== Vertical multiply ====
*
               DO 170 JROW = JTOP, MAX( KTOP, INCOL ) - 1, NV
                  JLEN = MIN( NV, MAX( KTOP, INCOL )-JROW )
                  CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE,
     $                        H( JROW, INCOL+K1 ), LDH, U( K1, K1 ),
     $                        LDU, ZERO, WV, LDWV )
                  CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV,
     $                         H( JROW, INCOL+K1 ), LDH )
  170          CONTINUE
*
*              ==== Z multiply (also vertical) ====
*
               IF( WANTZ ) THEN
                  DO 180 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
                     CALL DGEMM( 'N', 'N', JLEN, NU, NU, ONE,
     $                           Z( JROW, INCOL+K1 ), LDZ, U( K1, K1 ),
     $                           LDU, ZERO, WV, LDWV )
                     CALL DLACPY( 'ALL', JLEN, NU, WV, LDWV,
     $                            Z( JROW, INCOL+K1 ), LDZ )
  180             CONTINUE
               END IF
            ELSE
*
*              ==== Updates exploiting U's 2-by-2 block structure.
*              .    (I2, I4, J2, J4 are the last rows and columns
*              .    of the blocks.) ====
*
               I2 = ( KDU+1 ) / 2
               I4 = KDU
               J2 = I4 - I2
               J4 = KDU
*
*              ==== KZS and KNZ deal with the band of zeros
*              .    along the diagonal of one of the triangular
*              .    blocks. ====
*
               KZS = ( J4-J2 ) - ( NS+1 )
               KNZ = NS + 1
*
*              ==== Horizontal multiply ====
*
               DO 190 JCOL = MIN( NDCOL, KBOT ) + 1, JBOT, NH
                  JLEN = MIN( NH, JBOT-JCOL+1 )
*
*                 ==== Copy bottom of H to top+KZS of scratch ====
*                  (The first KZS rows get multiplied by zero.) ====
*
                  CALL DLACPY( 'ALL', KNZ, JLEN, H( INCOL+1+J2, JCOL ),
     $                         LDH, WH( KZS+1, 1 ), LDWH )
*
*                 ==== Multiply by U21' ====
*
                  CALL DLASET( 'ALL', KZS, JLEN, ZERO, ZERO, WH, LDWH )
                  CALL DTRMM( 'L', 'U', 'C', 'N', KNZ, JLEN, ONE,
     $                        U( J2+1, 1+KZS ), LDU, WH( KZS+1, 1 ),
     $                        LDWH )
*
*                 ==== Multiply top of H by U11' ====
*
                  CALL DGEMM( 'C', 'N', I2, JLEN, J2, ONE, U, LDU,
     $                        H( INCOL+1, JCOL ), LDH, ONE, WH, LDWH )
*
*                 ==== Copy top of H to bottom of WH ====
*
                  CALL DLACPY( 'ALL', J2, JLEN, H( INCOL+1, JCOL ), LDH,
     $                         WH( I2+1, 1 ), LDWH )
*
*                 ==== Multiply by U21' ====
*
                  CALL DTRMM( 'L', 'L', 'C', 'N', J2, JLEN, ONE,
     $                        U( 1, I2+1 ), LDU, WH( I2+1, 1 ), LDWH )
*
*                 ==== Multiply by U22 ====
*
                  CALL DGEMM( 'C', 'N', I4-I2, JLEN, J4-J2, ONE,
     $                        U( J2+1, I2+1 ), LDU,
     $                        H( INCOL+1+J2, JCOL ), LDH, ONE,
     $                        WH( I2+1, 1 ), LDWH )
*
*                 ==== Copy it back ====
*
                  CALL DLACPY( 'ALL', KDU, JLEN, WH, LDWH,
     $                         H( INCOL+1, JCOL ), LDH )
  190          CONTINUE
*
*              ==== Vertical multiply ====
*
               DO 200 JROW = JTOP, MAX( INCOL, KTOP ) - 1, NV
                  JLEN = MIN( NV, MAX( INCOL, KTOP )-JROW )
*
*                 ==== Copy right of H to scratch (the first KZS
*                 .    columns get multiplied by zero) ====
*
                  CALL DLACPY( 'ALL', JLEN, KNZ, H( JROW, INCOL+1+J2 ),
     $                         LDH, WV( 1, 1+KZS ), LDWV )
*
*                 ==== Multiply by U21 ====
*
                  CALL DLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV, LDWV )
                  CALL DTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,
     $                        U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),
     $                        LDWV )
*
*                 ==== Multiply by U11 ====
*
                  CALL DGEMM( 'N', 'N', JLEN, I2, J2, ONE,
     $                        H( JROW, INCOL+1 ), LDH, U, LDU, ONE, WV,
     $                        LDWV )
*
*                 ==== Copy left of H to right of scratch ====
*
                  CALL DLACPY( 'ALL', JLEN, J2, H( JROW, INCOL+1 ), LDH,
     $                         WV( 1, 1+I2 ), LDWV )
*
*                 ==== Multiply by U21 ====
*
                  CALL DTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,
     $                        U( 1, I2+1 ), LDU, WV( 1, 1+I2 ), LDWV )
*
*                 ==== Multiply by U22 ====
*
                  CALL DGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,
     $                        H( JROW, INCOL+1+J2 ), LDH,
     $                        U( J2+1, I2+1 ), LDU, ONE, WV( 1, 1+I2 ),
     $                        LDWV )
*
*                 ==== Copy it back ====
*
                  CALL DLACPY( 'ALL', JLEN, KDU, WV, LDWV,
     $                         H( JROW, INCOL+1 ), LDH )
  200          CONTINUE
*
*              ==== Multiply Z (also vertical) ====
*
               IF( WANTZ ) THEN
                  DO 210 JROW = ILOZ, IHIZ, NV
                     JLEN = MIN( NV, IHIZ-JROW+1 )
*
*                    ==== Copy right of Z to left of scratch (first
*                    .     KZS columns get multiplied by zero) ====
*
                     CALL DLACPY( 'ALL', JLEN, KNZ,
     $                            Z( JROW, INCOL+1+J2 ), LDZ,
     $                            WV( 1, 1+KZS ), LDWV )
*
*                    ==== Multiply by U12 ====
*
                     CALL DLASET( 'ALL', JLEN, KZS, ZERO, ZERO, WV,
     $                            LDWV )
                     CALL DTRMM( 'R', 'U', 'N', 'N', JLEN, KNZ, ONE,
     $                           U( J2+1, 1+KZS ), LDU, WV( 1, 1+KZS ),
     $                           LDWV )
*
*                    ==== Multiply by U11 ====
*
                     CALL DGEMM( 'N', 'N', JLEN, I2, J2, ONE,
     $                           Z( JROW, INCOL+1 ), LDZ, U, LDU, ONE,
     $                           WV, LDWV )
*
*                    ==== Copy left of Z to right of scratch ====
*
                     CALL DLACPY( 'ALL', JLEN, J2, Z( JROW, INCOL+1 ),
     $                            LDZ, WV( 1, 1+I2 ), LDWV )
*
*                    ==== Multiply by U21 ====
*
                     CALL DTRMM( 'R', 'L', 'N', 'N', JLEN, I4-I2, ONE,
     $                           U( 1, I2+1 ), LDU, WV( 1, 1+I2 ),
     $                           LDWV )
*
*                    ==== Multiply by U22 ====
*
                     CALL DGEMM( 'N', 'N', JLEN, I4-I2, J4-J2, ONE,
     $                           Z( JROW, INCOL+1+J2 ), LDZ,
     $                           U( J2+1, I2+1 ), LDU, ONE,
     $                           WV( 1, 1+I2 ), LDWV )
*
*                    ==== Copy the result back to Z ====
*
                     CALL DLACPY( 'ALL', JLEN, KDU, WV, LDWV,
     $                            Z( JROW, INCOL+1 ), LDZ )
  210             CONTINUE
               END IF
            END IF
         END IF
  220 CONTINUE
*
*     ==== End of DLAQR5 ====
*
      END
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dlaqrb
c
c\Description:
c  Compute the eigenvalues and the Schur decomposition of an upper 
c  Hessenberg submatrix in rows and columns ILO to IHI.  Only the
c  last component of the Schur vectors are computed.
c
c  This is mostly a modification of the LAPACK routine dlahqr.
c  
c\Usage:
c  call dlaqrb
c     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO )
c
c\Arguments
c  WANTT   Logical variable.  (INPUT)
c          = .TRUE. : the full Schur form T is required;
c          = .FALSE.: only eigenvalues are required.
c
c  N       Integer.  (INPUT)
c          The order of the matrix H.  N >= 0.
c
c  ILO     Integer.  (INPUT)
c  IHI     Integer.  (INPUT)
c          It is assumed that H is already upper quasi-triangular in
c          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
c          ILO = 1). SLAQRB works primarily with the Hessenberg
c          submatrix in rows and columns ILO to IHI, but applies
c          transformations to all of H if WANTT is .TRUE..
c          1 <= ILO <= max(1,IHI); IHI <= N.
c
c  H       Double precision array, dimension (LDH,N).  (INPUT/OUTPUT)
c          On entry, the upper Hessenberg matrix H.
c          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
c          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
c          standard form. If WANTT is .FALSE., the contents of H are
c          unspecified on exit.
c
c  LDH     Integer.  (INPUT)
c          The leading dimension of the array H. LDH >= max(1,N).
c
c  WR      Double precision array, dimension (N).  (OUTPUT)
c  WI      Double precision array, dimension (N).  (OUTPUT)
c          The real and imaginary parts, respectively, of the computed
c          eigenvalues ILO to IHI are stored in the corresponding
c          elements of WR and WI. If two eigenvalues are computed as a
c          complex conjugate pair, they are stored in consecutive
c          elements of WR and WI, say the i-th and (i+1)th, with
c          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
c          eigenvalues are stored in the same order as on the diagonal
c          of the Schur form returned in H, with WR(i) = H(i,i), and, if
c          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
c          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
c
c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit Z contains the last components of the Schur vectors.
c
c  INFO    Integer.  (OUPUT)
c          = 0: successful exit
c          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI
c               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
c               elements i+1:ihi of WR and WI contain those eigenvalues
c               which have been successfully computed.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlanv2  LAPACK routine that computes the Schur factorization of
c             2 by 2 nonsymmetric matrix in standard form.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     drot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix.

c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c               Modified from the LAPACK routine dlahqr so that only the
c               last component of the Schur vectors are computed.
c
c\SCCS Information: @(#) 
c FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dlaqrb ( wantt, n, ilo, ihi, h, ldh, wr, wi,
     &                    z, info )
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      logical    wantt
      integer    ihi, ilo, info, ldh, n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h( ldh, * ), wi( * ), wr( * ), z( * )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero, one, dat1, dat2
      parameter (zero = 0.0D+0, one = 1.0D+0, dat1 = 7.5D-1, 
     &           dat2 = -4.375D-1)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, i1, i2, itn, its, j, k, l, m, nh, nr
      Double precision
     &           cs, h00, h10, h11, h12, h21, h22, h33, h33s,
     &           h43h34, h44, h44s, ovfl, s, smlnum, sn, sum,
     &           t1, t2, t3, tst1, ulp, unfl, v1, v2, v3
      Double precision
     &           v( 3 ), work( 1 )
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch, dlanhs
      external   dlamch, dlanhs
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dlabad, dlanv2, dlarfg, drot
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      info = 0
c
c     %--------------------------%
c     | Quick return if possible |
c     %--------------------------%
c
      if( n.eq.0 )
     &   return
      if( ilo.eq.ihi ) then
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = zero
         return
      end if
c 
c     %---------------------------------------------%
c     | Initialize the vector of last components of |
c     | the Schur vectors for accumulation.         |
c     %---------------------------------------------%
c
      do 5 j = 1, n-1
         z(j) = zero
  5   continue 
      z(n) = one
c 
      nh = ihi - ilo + 1
c
c     %-------------------------------------------------------------%
c     | Set machine-dependent constants for the stopping criterion. |
c     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
c     %-------------------------------------------------------------%
c
      unfl = dlamch( 'safe minimum' )
      ovfl = one / unfl
      call dlabad( unfl, ovfl )
      ulp = dlamch( 'precision' )
      smlnum = unfl*( nh / ulp )
c
c     %---------------------------------------------------------------%
c     | I1 and I2 are the indices of the first row and last column    |
c     | of H to which transformations must be applied. If eigenvalues |
c     | only are computed, I1 and I2 are set inside the main loop.    |
c     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
c     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
c     %---------------------------------------------------------------%
c
      if( wantt ) then
         i1 = 1
         i2 = n
         do 8 i=1,i2-2
            h(i1+i+1,i) = zero
 8       continue
      else
         do 9 i=1, ihi-ilo-1
            h(ilo+i+1,ilo+i-1) = zero
 9       continue
      end if
c 
c     %---------------------------------------------------%
c     | ITN is the total number of QR iterations allowed. |
c     %---------------------------------------------------%
c
      itn = 30*nh
c 
c     ------------------------------------------------------------------
c     The main loop begins here. I is the loop index and decreases from
c     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
c     with the active submatrix in rows and columns L to I.
c     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
c     H(L,L-1) is negligible so that the matrix splits.
c     ------------------------------------------------------------------
c 
      i = ihi
   10 continue
      l = ilo
      if( i.lt.ilo )
     &   go to 150
 
c     %--------------------------------------------------------------%
c     | Perform QR iterations on rows and columns ILO to I until a   |
c     | submatrix of order 1 or 2 splits off at the bottom because a |
c     | subdiagonal element has become negligible.                   |
c     %--------------------------------------------------------------%
 
      do 130 its = 0, itn
c
c        %----------------------------------------------%
c        | Look for a single small subdiagonal element. |
c        %----------------------------------------------%
c
         do 20 k = i, l + 1, -1
            tst1 = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            if( tst1.eq.zero )
     &         tst1 = dlanhs( '1', i-l+1, h( l, l ), ldh, work )
            if( abs( h( k, k-1 ) ).le.max( ulp*tst1, smlnum ) )
     &         go to 30
   20    continue
   30    continue
         l = k
         if( l.gt.ilo ) then
c
c           %------------------------%
c           | H(L,L-1) is negligible |
c           %------------------------%
c
            h( l, l-1 ) = zero
         end if
c
c        %-------------------------------------------------------------%
c        | Exit from loop if a submatrix of order 1 or 2 has split off |
c        %-------------------------------------------------------------%
c
         if( l.ge.i-1 )
     &      go to 140
c
c        %---------------------------------------------------------%
c        | Now the active submatrix is in rows and columns L to I. |
c        | If eigenvalues only are being computed, only the active |
c        | submatrix need be transformed.                          |
c        %---------------------------------------------------------%
c
         if( .not.wantt ) then
            i1 = l
            i2 = i
         end if
c 
         if( its.eq.10 .or. its.eq.20 ) then
c
c           %-------------------%
c           | Exceptional shift |
c           %-------------------%
c
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h44 = dat1*s
            h33 = h44
            h43h34 = dat2*s*s
c
         else
c
c           %-----------------------------------------%
c           | Prepare to use Wilkinson's double shift |
c           %-----------------------------------------%
c
            h44 = h( i, i )
            h33 = h( i-1, i-1 )
            h43h34 = h( i, i-1 )*h( i-1, i )
         end if
c
c        %-----------------------------------------------------%
c        | Look for two consecutive small subdiagonal elements |
c        %-----------------------------------------------------%
c
         do 40 m = i - 2, l, -1
c
c           %---------------------------------------------------------%
c           | Determine the effect of starting the double-shift QR    |
c           | iteration at row M, and see if this would make H(M,M-1) |
c           | negligible.                                             |
c           %---------------------------------------------------------%
c
            h11 = h( m, m )
            h22 = h( m+1, m+1 )
            h21 = h( m+1, m )
            h12 = h( m, m+1 )
            h44s = h44 - h11
            h33s = h33 - h11
            v1 = ( h33s*h44s-h43h34 ) / h21 + h12
            v2 = h22 - h11 - h33s - h44s
            v3 = h( m+2, m+1 )
            s = abs( v1 ) + abs( v2 ) + abs( v3 )
            v1 = v1 / s
            v2 = v2 / s
            v3 = v3 / s
            v( 1 ) = v1
            v( 2 ) = v2
            v( 3 ) = v3
            if( m.eq.l )
     &         go to 50
            h00 = h( m-1, m-1 )
            h10 = h( m, m-1 )
            tst1 = abs( v1 )*( abs( h00 )+abs( h11 )+abs( h22 ) )
            if( abs( h10 )*( abs( v2 )+abs( v3 ) ).le.ulp*tst1 )
     &         go to 50
   40    continue
   50    continue
c
c        %----------------------%
c        | Double-shift QR step |
c        %----------------------%
c
         do 120 k = m, i - 1
c 
c           ------------------------------------------------------------
c           The first iteration of this loop determines a reflection G
c           from the vector V and applies it from left and right to H,
c           thus creating a nonzero bulge below the subdiagonal.
c
c           Each subsequent iteration determines a reflection G to
c           restore the Hessenberg form in the (K-1)th column, and thus
c           chases the bulge one step toward the bottom of the active
c           submatrix. NR is the order of G.
c           ------------------------------------------------------------
c 
            nr = min( 3, i-k+1 )
            if( k.gt.m )
     &         call dcopy( nr, h( k, k-1 ), 1, v, 1 )
            call dlarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            if( k.gt.m ) then
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
               if( k.lt.i-1 )
     &            h( k+2, k-1 ) = zero
            else if( m.gt.l ) then
               h( k, k-1 ) = -h( k, k-1 )
            end if
            v2 = v( 2 )
            t2 = t1*v2
            if( nr.eq.3 ) then
               v3 = v( 3 )
               t3 = t1*v3
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
               do 60 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
                  h( k+2, j ) = h( k+2, j ) - sum*t3
   60          continue
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
               do 70 j = i1, min( k+3, i )
                  sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
                  h( j, k+2 ) = h( j, k+2 ) - sum*t3
   70          continue
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
               sum      = z( k ) + v2*z( k+1 ) + v3*z( k+2 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
               z( k+2 ) = z( k+2 ) - sum*t3
 
            else if( nr.eq.2 ) then
c
c              %------------------------------------------------%
c              | Apply G from the left to transform the rows of |
c              | the matrix in columns K to I2.                 |
c              %------------------------------------------------%
c
               do 90 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
   90          continue
c
c              %----------------------------------------------------%
c              | Apply G from the right to transform the columns of |
c              | the matrix in rows I1 to min(K+3,I).               |
c              %----------------------------------------------------%
c
               do 100 j = i1, i
                  sum = h( j, k ) + v2*h( j, k+1 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
  100          continue
c
c              %----------------------------------%
c              | Accumulate transformations for Z |
c              %----------------------------------%
c
               sum      = z( k ) + v2*z( k+1 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
            end if
  120    continue
 
  130 continue
c
c     %-------------------------------------------------------%
c     | Failure to converge in remaining number of iterations |
c     %-------------------------------------------------------%
c
      info = i
      return
 
  140 continue
 
      if( l.eq.i ) then
c
c        %------------------------------------------------------%
c        | H(I,I-1) is negligible: one eigenvalue has converged |
c        %------------------------------------------------------%
c
         wr( i ) = h( i, i )
         wi( i ) = zero

      else if( l.eq.i-1 ) then
c
c        %--------------------------------------------------------%
c        | H(I-1,I-2) is negligible;                              |
c        | a pair of eigenvalues have converged.                  |
c        |                                                        |
c        | Transform the 2-by-2 submatrix to standard Schur form, |
c        | and compute and store the eigenvalues.                 |
c        %--------------------------------------------------------%
c
         call dlanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ),
     &                h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ),
     &                cs, sn )
 
         if( wantt ) then
c
c           %-----------------------------------------------------%
c           | Apply the transformation to the rest of H and to Z, |
c           | as required.                                        |
c           %-----------------------------------------------------%
c
            if( i2.gt.i )
     &         call drot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh,
     &                    cs, sn )
            call drot( i-i1-1, h( i1, i-1 ), 1, h( i1, i ), 1, cs, sn )
            sum      = cs*z( i-1 ) + sn*z( i )
            z( i )   = cs*z( i )   - sn*z( i-1 )
            z( i-1 ) = sum
         end if
      end if
c
c     %---------------------------------------------------------%
c     | Decrement number of remaining iterations, and return to |
c     | start of the main loop with new value of I.             |
c     %---------------------------------------------------------%
c
      itn = itn - its
      i = l - 1
      go to 10
 
  150 continue
      return
c
c     %---------------%
c     | End of dlaqrb |
c     %---------------%
c
      end
      SUBROUTINE DLAQTR( LTRAN, LREAL, N, T, LDT, B, W, SCALE, X, WORK,
     $                   INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            LREAL, LTRAN
      INTEGER            INFO, LDT, N
      DOUBLE PRECISION   SCALE, W
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( * ), T( LDT, * ), WORK( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAQTR solves the real quasi-triangular system
*
*               op(T)*p = scale*c,               if LREAL = .TRUE.
*
*  or the complex quasi-triangular systems
*
*             op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE.
*
*  in real arithmetic, where T is upper quasi-triangular.
*  If LREAL = .FALSE., then the first diagonal block of T must be
*  1 by 1, B is the specially structured matrix
*
*                 B = [ b(1) b(2) ... b(n) ]
*                     [       w            ]
*                     [           w        ]
*                     [              .     ]
*                     [                 w  ]
*
*  op(A) = A or A', A' denotes the conjugate transpose of
*  matrix A.
*
*  On input, X = [ c ].  On output, X = [ p ].
*                [ d ]                  [ q ]
*
*  This subroutine is designed for the condition number estimation
*  in routine DTRSNA.
*
*  Arguments
*  =========
*
*  LTRAN   (input) LOGICAL
*          On entry, LTRAN specifies the option of conjugate transpose:
*             = .FALSE.,    op(T+i*B) = T+i*B,
*             = .TRUE.,     op(T+i*B) = (T+i*B)'.
*
*  LREAL   (input) LOGICAL
*          On entry, LREAL specifies the input matrix structure:
*             = .FALSE.,    the input is complex
*             = .TRUE.,     the input is real
*
*  N       (input) INTEGER
*          On entry, N specifies the order of T+i*B. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, T contains a matrix in Schur canonical form.
*          If LREAL = .FALSE., then the first diagonal block of T mu
*          be 1 by 1.
*
*  LDT     (input) INTEGER
*          The leading dimension of the matrix T. LDT >= max(1,N).
*
*  B       (input) DOUBLE PRECISION array, dimension (N)
*          On entry, B contains the elements to form the matrix
*          B as described above.
*          If LREAL = .TRUE., B is not referenced.
*
*  W       (input) DOUBLE PRECISION
*          On entry, W is the diagonal element of the matrix B.
*          If LREAL = .TRUE., W is not referenced.
*
*  SCALE   (output) DOUBLE PRECISION
*          On exit, SCALE is the scale factor.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (2*N)
*          On entry, X contains the right hand side of the system.
*          On exit, X is overwritten by the solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          On exit, INFO is set to
*             0: successful exit.
*               1: the some diagonal 1 by 1 block has been perturbed by
*                  a small number SMIN to keep nonsingularity.
*               2: the some diagonal 2 by 2 block has been perturbed by
*                  a small number in DLALN2 to keep nonsingularity.
*          NOTE: In the interests of speed, this routine does not
*                check the inputs for errors.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            I, IERR, J, J1, J2, JNEXT, K, N1, N2
      DOUBLE PRECISION   BIGNUM, EPS, REC, SCALOC, SI, SMIN, SMINW,
     $                   SMLNUM, SR, TJJ, TMP, XJ, XMAX, XNORM, Z
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   D( 2, 2 ), V( 2, 2 )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH, DLANGE
      EXTERNAL           IDAMAX, DASUM, DDOT, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLADIV, DLALN2, DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Do not test the input parameters for errors
*
      NOTRAN = .NOT.LTRAN
      INFO = 0
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Set constants to control overflow
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
*
      XNORM = DLANGE( 'M', N, N, T, LDT, D )
      IF( .NOT.LREAL )
     $   XNORM = MAX( XNORM, ABS( W ), DLANGE( 'M', N, 1, B, N, D ) )
      SMIN = MAX( SMLNUM, EPS*XNORM )
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      WORK( 1 ) = ZERO
      DO 10 J = 2, N
         WORK( J ) = DASUM( J-1, T( 1, J ), 1 )
   10 CONTINUE
*
      IF( .NOT.LREAL ) THEN
         DO 20 I = 2, N
            WORK( I ) = WORK( I ) + ABS( B( I ) )
   20    CONTINUE
      END IF
*
      N2 = 2*N
      N1 = N
      IF( .NOT.LREAL )
     $   N1 = N2
      K = IDAMAX( N1, X, 1 )
      XMAX = ABS( X( K ) )
      SCALE = ONE
*
      IF( XMAX.GT.BIGNUM ) THEN
         SCALE = BIGNUM / XMAX
         CALL DSCAL( N1, SCALE, X, 1 )
         XMAX = BIGNUM
      END IF
*
      IF( LREAL ) THEN
*
         IF( NOTRAN ) THEN
*
*           Solve T*p = scale*c
*
            JNEXT = N
            DO 30 J = N, 1, -1
               IF( J.GT.JNEXT )
     $            GO TO 30
               J1 = J
               J2 = J
               JNEXT = J - 1
               IF( J.GT.1 ) THEN
                  IF( T( J, J-1 ).NE.ZERO ) THEN
                     J1 = J - 1
                     JNEXT = J - 2
                  END IF
               END IF
*
               IF( J1.EQ.J2 ) THEN
*
*                 Meet 1 by 1 diagonal block
*
*                 Scale to avoid overflow when computing
*                     x(j) = b(j)/T(j,j)
*
                  XJ = ABS( X( J1 ) )
                  TJJ = ABS( T( J1, J1 ) )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMIN ) THEN
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  END IF
*
                  IF( XJ.EQ.ZERO )
     $               GO TO 30
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J1 ) = X( J1 ) / TMP
                  XJ = ABS( X( J1 ) )
*
*                 Scale x if necessary to avoid overflow when adding a
*                 multiple of column j1 of T.
*
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( WORK( J1 ).GT.( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     K = IDAMAX( J1-1, X, 1 )
                     XMAX = ABS( X( K ) )
                  END IF
*
               ELSE
*
*                 Meet 2 by 2 diagonal block
*
*                 Call 2 by 2 linear system solve, to take
*                 care of possible overflow by scaling factor.
*
                  D( 1, 1 ) = X( J1 )
                  D( 2, 1 ) = X( J2 )
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, T( J1, J1 ),
     $                         LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2,
     $                         SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 2
*
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( N, SCALOC, X, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
*
*                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
*                 to avoid overflow in updating right-hand side.
*
                  XJ = MAX( ABS( V( 1, 1 ) ), ABS( V( 2, 1 ) ) )
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( MAX( WORK( J1 ), WORK( J2 ) ).GT.
     $                   ( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
*
*                 Update right-hand side
*
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     CALL DAXPY( J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 )
                     K = IDAMAX( J1-1, X, 1 )
                     XMAX = ABS( X( K ) )
                  END IF
*
               END IF
*
   30       CONTINUE
*
         ELSE
*
*           Solve T'*p = scale*c
*
            JNEXT = 1
            DO 40 J = 1, N
               IF( J.LT.JNEXT )
     $            GO TO 40
               J1 = J
               J2 = J
               JNEXT = J + 1
               IF( J.LT.N ) THEN
                  IF( T( J+1, J ).NE.ZERO ) THEN
                     J2 = J + 1
                     JNEXT = J + 2
                  END IF
               END IF
*
               IF( J1.EQ.J2 ) THEN
*
*                 1 by 1 diagonal block
*
*                 Scale if necessary to avoid overflow in forming the
*                 right-hand side element by inner product.
*
                  XJ = ABS( X( J1 ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( WORK( J1 ).GT.( BIGNUM-XJ )*REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
*
                  X( J1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )
*
                  XJ = ABS( X( J1 ) )
                  TJJ = ABS( T( J1, J1 ) )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMIN ) THEN
                     TMP = SMIN
                     TJJ = SMIN
                     INFO = 1
                  END IF
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J1 ) = X( J1 ) / TMP
                  XMAX = MAX( XMAX, ABS( X( J1 ) ) )
*
               ELSE
*
*                 2 by 2 diagonal block
*
*                 Scale if necessary to avoid overflow in forming the
*                 right-hand side elements by inner product.
*
                  XJ = MAX( ABS( X( J1 ) ), ABS( X( J2 ) ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( MAX( WORK( J2 ), WORK( J1 ) ).GT.( BIGNUM-XJ )*
     $                   REC ) THEN
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
*
                  D( 1, 1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X,
     $                        1 )
                  D( 2, 1 ) = X( J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X,
     $                        1 )
*
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, T( J1, J1 ),
     $                         LDT, ONE, ONE, D, 2, ZERO, ZERO, V, 2,
     $                         SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 2
*
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( N, SCALOC, X, 1 )
                     SCALE = SCALE*SCALOC
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  XMAX = MAX( ABS( X( J1 ) ), ABS( X( J2 ) ), XMAX )
*
               END IF
   40       CONTINUE
         END IF
*
      ELSE
*
         SMINW = MAX( EPS*ABS( W ), SMIN )
         IF( NOTRAN ) THEN
*
*           Solve (T + iB)*(p+iq) = c+id
*
            JNEXT = N
            DO 70 J = N, 1, -1
               IF( J.GT.JNEXT )
     $            GO TO 70
               J1 = J
               J2 = J
               JNEXT = J - 1
               IF( J.GT.1 ) THEN
                  IF( T( J, J-1 ).NE.ZERO ) THEN
                     J1 = J - 1
                     JNEXT = J - 2
                  END IF
               END IF
*
               IF( J1.EQ.J2 ) THEN
*
*                 1 by 1 diagonal block
*
*                 Scale if necessary to avoid overflow in division
*
                  Z = W
                  IF( J1.EQ.1 )
     $               Z = B( 1 )
                  XJ = ABS( X( J1 ) ) + ABS( X( N+J1 ) )
                  TJJ = ABS( T( J1, J1 ) ) + ABS( Z )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMINW ) THEN
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  END IF
*
                  IF( XJ.EQ.ZERO )
     $               GO TO 70
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  CALL DLADIV( X( J1 ), X( N+J1 ), TMP, Z, SR, SI )
                  X( J1 ) = SR
                  X( N+J1 ) = SI
                  XJ = ABS( X( J1 ) ) + ABS( X( N+J1 ) )
*
*                 Scale x if necessary to avoid overflow when adding a
*                 multiple of column j1 of T.
*
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( WORK( J1 ).GT.( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
*
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     CALL DAXPY( J1-1, -X( N+J1 ), T( 1, J1 ), 1,
     $                           X( N+1 ), 1 )
*
                     X( 1 ) = X( 1 ) + B( J1 )*X( N+J1 )
                     X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 )
*
                     XMAX = ZERO
                     DO 50 K = 1, J1 - 1
                        XMAX = MAX( XMAX, ABS( X( K ) )+
     $                         ABS( X( K+N ) ) )
   50                CONTINUE
                  END IF
*
               ELSE
*
*                 Meet 2 by 2 diagonal block
*
                  D( 1, 1 ) = X( J1 )
                  D( 2, 1 ) = X( J2 )
                  D( 1, 2 ) = X( N+J1 )
                  D( 2, 2 ) = X( N+J2 )
                  CALL DLALN2( .FALSE., 2, 2, SMINW, ONE, T( J1, J1 ),
     $                         LDT, ONE, ONE, D, 2, ZERO, -W, V, 2,
     $                         SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 2
*
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( 2*N, SCALOC, X, 1 )
                     SCALE = SCALOC*SCALE
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  X( N+J1 ) = V( 1, 2 )
                  X( N+J2 ) = V( 2, 2 )
*
*                 Scale X(J1), .... to avoid overflow in
*                 updating right hand side.
*
                  XJ = MAX( ABS( V( 1, 1 ) )+ABS( V( 1, 2 ) ),
     $                 ABS( V( 2, 1 ) )+ABS( V( 2, 2 ) ) )
                  IF( XJ.GT.ONE ) THEN
                     REC = ONE / XJ
                     IF( MAX( WORK( J1 ), WORK( J2 ) ).GT.
     $                   ( BIGNUM-XMAX )*REC ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                     END IF
                  END IF
*
*                 Update the right-hand side.
*
                  IF( J1.GT.1 ) THEN
                     CALL DAXPY( J1-1, -X( J1 ), T( 1, J1 ), 1, X, 1 )
                     CALL DAXPY( J1-1, -X( J2 ), T( 1, J2 ), 1, X, 1 )
*
                     CALL DAXPY( J1-1, -X( N+J1 ), T( 1, J1 ), 1,
     $                           X( N+1 ), 1 )
                     CALL DAXPY( J1-1, -X( N+J2 ), T( 1, J2 ), 1,
     $                           X( N+1 ), 1 )
*
                     X( 1 ) = X( 1 ) + B( J1 )*X( N+J1 ) +
     $                        B( J2 )*X( N+J2 )
                     X( N+1 ) = X( N+1 ) - B( J1 )*X( J1 ) -
     $                          B( J2 )*X( J2 )
*
                     XMAX = ZERO
                     DO 60 K = 1, J1 - 1
                        XMAX = MAX( ABS( X( K ) )+ABS( X( K+N ) ),
     $                         XMAX )
   60                CONTINUE
                  END IF
*
               END IF
   70       CONTINUE
*
         ELSE
*
*           Solve (T + iB)'*(p+iq) = c+id
*
            JNEXT = 1
            DO 80 J = 1, N
               IF( J.LT.JNEXT )
     $            GO TO 80
               J1 = J
               J2 = J
               JNEXT = J + 1
               IF( J.LT.N ) THEN
                  IF( T( J+1, J ).NE.ZERO ) THEN
                     J2 = J + 1
                     JNEXT = J + 2
                  END IF
               END IF
*
               IF( J1.EQ.J2 ) THEN
*
*                 1 by 1 diagonal block
*
*                 Scale if necessary to avoid overflow in forming the
*                 right-hand side element by inner product.
*
                  XJ = ABS( X( J1 ) ) + ABS( X( J1+N ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( WORK( J1 ).GT.( BIGNUM-XJ )*REC ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
*
                  X( J1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X, 1 )
                  X( N+J1 ) = X( N+J1 ) - DDOT( J1-1, T( 1, J1 ), 1,
     $                        X( N+1 ), 1 )
                  IF( J1.GT.1 ) THEN
                     X( J1 ) = X( J1 ) - B( J1 )*X( N+1 )
                     X( N+J1 ) = X( N+J1 ) + B( J1 )*X( 1 )
                  END IF
                  XJ = ABS( X( J1 ) ) + ABS( X( J1+N ) )
*
                  Z = W
                  IF( J1.EQ.1 )
     $               Z = B( 1 )
*
*                 Scale if necessary to avoid overflow in
*                 complex division
*
                  TJJ = ABS( T( J1, J1 ) ) + ABS( Z )
                  TMP = T( J1, J1 )
                  IF( TJJ.LT.SMINW ) THEN
                     TMP = SMINW
                     TJJ = SMINW
                     INFO = 1
                  END IF
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.BIGNUM*TJJ ) THEN
                        REC = ONE / XJ
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  CALL DLADIV( X( J1 ), X( N+J1 ), TMP, -Z, SR, SI )
                  X( J1 ) = SR
                  X( J1+N ) = SI
                  XMAX = MAX( ABS( X( J1 ) )+ABS( X( J1+N ) ), XMAX )
*
               ELSE
*
*                 2 by 2 diagonal block
*
*                 Scale if necessary to avoid overflow in forming the
*                 right-hand side element by inner product.
*
                  XJ = MAX( ABS( X( J1 ) )+ABS( X( N+J1 ) ),
     $                 ABS( X( J2 ) )+ABS( X( N+J2 ) ) )
                  IF( XMAX.GT.ONE ) THEN
                     REC = ONE / XMAX
                     IF( MAX( WORK( J1 ), WORK( J2 ) ).GT.
     $                   ( BIGNUM-XJ ) / XMAX ) THEN
                        CALL DSCAL( N2, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
*
                  D( 1, 1 ) = X( J1 ) - DDOT( J1-1, T( 1, J1 ), 1, X,
     $                        1 )
                  D( 2, 1 ) = X( J2 ) - DDOT( J1-1, T( 1, J2 ), 1, X,
     $                        1 )
                  D( 1, 2 ) = X( N+J1 ) - DDOT( J1-1, T( 1, J1 ), 1,
     $                        X( N+1 ), 1 )
                  D( 2, 2 ) = X( N+J2 ) - DDOT( J1-1, T( 1, J2 ), 1,
     $                        X( N+1 ), 1 )
                  D( 1, 1 ) = D( 1, 1 ) - B( J1 )*X( N+1 )
                  D( 2, 1 ) = D( 2, 1 ) - B( J2 )*X( N+1 )
                  D( 1, 2 ) = D( 1, 2 ) + B( J1 )*X( 1 )
                  D( 2, 2 ) = D( 2, 2 ) + B( J2 )*X( 1 )
*
                  CALL DLALN2( .TRUE., 2, 2, SMINW, ONE, T( J1, J1 ),
     $                         LDT, ONE, ONE, D, 2, ZERO, W, V, 2,
     $                         SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 2
*
                  IF( SCALOC.NE.ONE ) THEN
                     CALL DSCAL( N2, SCALOC, X, 1 )
                     SCALE = SCALOC*SCALE
                  END IF
                  X( J1 ) = V( 1, 1 )
                  X( J2 ) = V( 2, 1 )
                  X( N+J1 ) = V( 1, 2 )
                  X( N+J2 ) = V( 2, 2 )
                  XMAX = MAX( ABS( X( J1 ) )+ABS( X( N+J1 ) ),
     $                   ABS( X( J2 ) )+ABS( X( N+J2 ) ), XMAX )
*
               END IF
*
   80       CONTINUE
*
         END IF
*
      END IF
*
      RETURN
*
*     End of DLAQTR
*
      END
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILADLC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILADLR(M, LASTV, C, LDC)
         END IF
      END IF
!     Note that lastc.eq.0 renders the BLAS operations null; no special
!     case is needed at this level.
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastv,1:lastc)' * v(1:lastv,1)
*
            CALL DGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,
     $           ZERO, WORK, 1 )
*
*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)'
*
            CALL DGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
            CALL DGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC,
     $           V, INCV, ZERO, WORK, 1 )
*
*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)'
*
            CALL DGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of DLARF
*
      END
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFB applies a real block reflector H or its transpose H' to a
*  real m by n matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'T': apply H' (Transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) DOUBLE PRECISION array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
*          The triangular k by k matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J, LASTV, LASTC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILADLR, ILADLC
      EXTERNAL           LSAME, ILADLR, ILADLC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
               LASTV = MAX( K, ILADLR( M, K, V, LDV ) )
               LASTC = ILADLC( LASTV, N, C, LDC )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL DCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL DGEMM( 'Transpose', 'No transpose',
     $                 LASTC, K, LASTV-K,
     $                 ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( LASTV.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL DGEMM( 'No transpose', 'Transpose',
     $                 LASTV-K, LASTC, K,
     $                 -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                 C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, LASTC
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
               LASTV = MAX( K, ILADLR( N, K, V, LDV ) )
               LASTC = ILADLR( M, LASTV, C, LDC )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL DCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose',
     $                 LASTC, K, LASTV-K,
     $                 ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( LASTV.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL DGEMM( 'No transpose', 'Transpose',
     $                 LASTC, LASTV-K, K,
     $                 -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                 C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, LASTC
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
               LASTV = MAX( K, ILADLR( M, K, V, LDV ) )
               LASTC = ILADLC( LASTV, N, C, LDC )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL DCOPY( LASTC, C( LASTV-K+J, 1 ), LDC,
     $                 WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV,
     $              WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL DGEMM( 'Transpose', 'No transpose',
     $                 LASTC, K, LASTV-K, ONE, C, LDC, V, LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( LASTV.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL DGEMM( 'No transpose', 'Transpose',
     $                 LASTV-K, LASTC, K, -ONE, V, LDV, WORK, LDWORK,
     $                 ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV,
     $              WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, LASTC
                     C( LASTV-K+J, I ) = C( LASTV-K+J, I ) - WORK(I, J)
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
               LASTV = MAX( K, ILADLR( N, K, V, LDV ) )
               LASTC = ILADLR( M, LASTV, C, LDC )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL DCOPY( LASTC, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV,
     $              WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose',
     $                 LASTC, K, LASTV-K, ONE, C, LDC, V, LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( LASTV.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL DGEMM( 'No transpose', 'Transpose',
     $                 LASTC, LASTV-K, K, -ONE, WORK, LDWORK, V, LDV,
     $                 ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V( LASTV-K+1, 1 ), LDV,
     $              WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, LASTC
                     C( I, LASTV-K+J ) = C( I, LASTV-K+J ) - WORK(I, J)
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
               LASTV = MAX( K, ILADLC( K, M, V, LDV ) )
               LASTC = ILADLC( LASTV, N, C, LDC )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL DCOPY( LASTC, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL DGEMM( 'Transpose', 'Transpose',
     $                 LASTC, K, LASTV-K,
     $                 ONE, C( K+1, 1 ), LDC, V( 1, K+1 ), LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( LASTV.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL DGEMM( 'Transpose', 'Transpose',
     $                 LASTV-K, LASTC, K,
     $                 -ONE, V( 1, K+1 ), LDV, WORK, LDWORK,
     $                 ONE, C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, LASTC
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
               LASTV = MAX( K, ILADLC( K, N, V, LDV ) )
               LASTC = ILADLR( M, LASTV, C, LDC )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL DCOPY( LASTC, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL DGEMM( 'No transpose', 'Transpose',
     $                 LASTC, K, LASTV-K,
     $                 ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( LASTV.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL DGEMM( 'No transpose', 'No transpose',
     $                 LASTC, LASTV-K, K,
     $                 -ONE, WORK, LDWORK, V( 1, K+1 ), LDV,
     $                 ONE, C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, LASTC
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
               LASTV = MAX( K, ILADLC( K, M, V, LDV ) )
               LASTC = ILADLC( LASTV, N, C, LDC )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL DCOPY( LASTC, C( LASTV-K+J, 1 ), LDC,
     $                 WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV,
     $              WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL DGEMM( 'Transpose', 'Transpose',
     $                 LASTC, K, LASTV-K, ONE, C, LDC, V, LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( LASTV.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL DGEMM( 'Transpose', 'Transpose',
     $                 LASTV-K, LASTC, K, -ONE, V, LDV, WORK, LDWORK,
     $                 ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV,
     $              WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, LASTC
                     C( LASTV-K+J, I ) = C( LASTV-K+J, I ) - WORK(I, J)
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
               LASTV = MAX( K, ILADLC( K, N, V, LDV ) )
               LASTC = ILADLR( M, LASTV, C, LDC )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL DCOPY( LASTC, C( 1, LASTV-K+J ), 1,
     $                 WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit',
     $              LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV,
     $              WORK, LDWORK )
               IF( LASTV.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL DGEMM( 'No transpose', 'Transpose',
     $                 LASTC, K, LASTV-K, ONE, C, LDC, V, LDV,
     $                 ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit',
     $              LASTC, K, ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( LASTV.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL DGEMM( 'No transpose', 'No transpose',
     $                 LASTC, LASTV-K, K, -ONE, WORK, LDWORK, V, LDV,
     $                 ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit',
     $              LASTC, K, ONE, V( 1, LASTV-K+1 ), LDV,
     $              WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, LASTC
                     C( I, LASTV-K+J ) = C( I, LASTV-K+J ) - WORK(I, J)
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of DLARFB
*
      END
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         TAU = ( BETA-ALPHA ) / BETA
         CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of DLARFG
*
      END
      SUBROUTINE DLARFP( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFP generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, beta is non-negative, and x is
*  an (n-1)-element real vector.  H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0
*
         IF( ALPHA.GE.ZERO ) THEN
!           When TAU.eq.ZERO, the vector is special-cased to be
!           all zeros in the application routines.  We do not need
!           to clear it.
            TAU = ZERO
         ELSE
!           However, the application routines rely on explicit
!           zero checks when TAU.ne.ZERO, and we must clear X.
            TAU = TWO
            DO J = 1, N-1
               X( 1 + (J-1)*INCX ) = 0
            END DO
            ALPHA = -ALPHA
         END IF
      ELSE
*
*        general case
*
         BETA = SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         KNT = 0
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            RSAFMN = ONE / SAFMIN
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         END IF
         ALPHA = ALPHA + BETA
         IF( BETA.LT.ZERO ) THEN
            BETA = -BETA
            TAU = -ALPHA / BETA
         ELSE
            ALPHA = XNORM * (XNORM/ALPHA)
            TAU = ALPHA / BETA
            ALPHA = -ALPHA
         END IF
         CALL DSCAL( N-1, ONE / ALPHA, X, INCX )
*
*        If BETA is subnormal, it may lose relative accuracy
*
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
*
      RETURN
*
*     End of DLARFP
*
      END

      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFT forms the triangular factor T of a real block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) DOUBLE PRECISION array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
      DOUBLE PRECISION   VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         PREVLASTV = N
         DO 20 I = 1, K
            PREVLASTV = MAX( I, PREVLASTV )
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( LASTV, I ).NE.ZERO ) EXIT
                  END DO
                  J = MIN( LASTV, PREVLASTV )
*
*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)' * V(i:j,i)
*
                  CALL DGEMV( 'Transpose', J-I+1, I-1, -TAU( I ),
     $                        V( I, 1 ), LDV, V( I, I ), 1, ZERO,
     $                        T( 1, I ), 1 )
               ELSE
!                 Skip any trailing zeros.
                  DO LASTV = N, I+1, -1
                     IF( V( I, LASTV ).NE.ZERO ) EXIT
                  END DO
                  J = MIN( LASTV, PREVLASTV )
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)'
*
                  CALL DGEMV( 'No transpose', I-1, J-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
               IF( I.GT.1 ) THEN
                  PREVLASTV = MAX( PREVLASTV, LASTV )
               ELSE
                  PREVLASTV = LASTV
               END IF
            END IF
   20    CONTINUE
      ELSE
         PREVLASTV = 1
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( LASTV, I ).NE.ZERO ) EXIT
                     END DO
                     J = MAX( LASTV, PREVLASTV )
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(j:n-k+i,i+1:k)' * V(j:n-k+i,i)
*
                     CALL DGEMV( 'Transpose', N-K+I-J+1, K-I, -TAU( I ),
     $                           V( J, I+1 ), LDV, V( J, I ), 1, ZERO,
     $                           T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
!                    Skip any leading zeros.
                     DO LASTV = 1, I-1
                        IF( V( I, LASTV ).NE.ZERO ) EXIT
                     END DO
                     J = MAX( LASTV, PREVLASTV )
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)'
*
                     CALL DGEMV( 'No transpose', K-I, N-K+I-J+1,
     $                    -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV,
     $                    ZERO, T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                  IF( I.GT.1 ) THEN
                     PREVLASTV = MIN( PREVLASTV, LASTV )
                  ELSE
                     PREVLASTV = LASTV
                  END IF
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of DLARFT
*
      END
      SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARFX applies a real elementary reflector H to a real m by n
*  matrix C, from either the left or the right. H is represented in the
*  form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix
*
*  This version uses inline code if H has order < 11.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension (M) if SIDE = 'L'
*                                     or (N) if SIDE = 'R'
*          The vector v in the representation of H.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= (1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                      (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*          WORK is not referenced if H has order < 11.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9,
     $                   V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF
*     ..
*     .. Executable Statements ..
*
      IF( TAU.EQ.ZERO )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C, where H has order m.
*
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150,
     $           170, 190 )M
*
*        Code for general M
*
         CALL DLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
         GO TO 410
   10    CONTINUE
*
*        Special code for 1 x 1 Householder
*
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 20 J = 1, N
            C( 1, J ) = T1*C( 1, J )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
*
*        Special code for 2 x 2 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 40 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
   40    CONTINUE
         GO TO 410
   50    CONTINUE
*
*        Special code for 3 x 3 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 60 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
   60    CONTINUE
         GO TO 410
   70    CONTINUE
*
*        Special code for 4 x 4 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 80 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
   80    CONTINUE
         GO TO 410
   90    CONTINUE
*
*        Special code for 5 x 5 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 100 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J ) + V5*C( 5, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
  100    CONTINUE
         GO TO 410
  110    CONTINUE
*
*        Special code for 6 x 6 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 120 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
  120    CONTINUE
         GO TO 410
  130    CONTINUE
*
*        Special code for 7 x 7 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 140 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
     $            V7*C( 7, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
  140    CONTINUE
         GO TO 410
  150    CONTINUE
*
*        Special code for 8 x 8 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 160 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
     $            V7*C( 7, J ) + V8*C( 8, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
  160    CONTINUE
         GO TO 410
  170    CONTINUE
*
*        Special code for 9 x 9 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 180 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
     $            V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
  180    CONTINUE
         GO TO 410
  190    CONTINUE
*
*        Special code for 10 x 10 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 200 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) +
     $            V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) +
     $            V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) +
     $            V10*C( 10, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
            C( 10, J ) = C( 10, J ) - SUM*T10
  200    CONTINUE
         GO TO 410
      ELSE
*
*        Form  C * H, where H has order n.
*
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350,
     $           370, 390 )N
*
*        Code for general N
*
         CALL DLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
         GO TO 410
  210    CONTINUE
*
*        Special code for 1 x 1 Householder
*
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 220 J = 1, M
            C( J, 1 ) = T1*C( J, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
*
*        Special code for 2 x 2 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 240 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
  240    CONTINUE
         GO TO 410
  250    CONTINUE
*
*        Special code for 3 x 3 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 260 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
  260    CONTINUE
         GO TO 410
  270    CONTINUE
*
*        Special code for 4 x 4 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 280 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
  280    CONTINUE
         GO TO 410
  290    CONTINUE
*
*        Special code for 5 x 5 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 300 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 ) + V5*C( J, 5 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
  300    CONTINUE
         GO TO 410
  310    CONTINUE
*
*        Special code for 6 x 6 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 320 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
  320    CONTINUE
         GO TO 410
  330    CONTINUE
*
*        Special code for 7 x 7 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 340 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
     $            V7*C( J, 7 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
  340    CONTINUE
         GO TO 410
  350    CONTINUE
*
*        Special code for 8 x 8 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 360 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
     $            V7*C( J, 7 ) + V8*C( J, 8 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
  360    CONTINUE
         GO TO 410
  370    CONTINUE
*
*        Special code for 9 x 9 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 380 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
     $            V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
  380    CONTINUE
         GO TO 410
  390    CONTINUE
*
*        Special code for 10 x 10 Householder
*
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 400 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) +
     $            V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) +
     $            V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) +
     $            V10*C( J, 10 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
            C( J, 10 ) = C( J, 10 ) - SUM*T10
  400    CONTINUE
         GO TO 410
      END IF
  410 CONTINUE
      RETURN
*
*     End of DLARFX
*
      END
      SUBROUTINE DLARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLARNV returns a vector of n random real numbers from a uniform or
*  normal distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  uniform (0,1)
*          = 2:  uniform (-1,1)
*          = 3:  normal (0,1)
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  N       (input) INTEGER
*          The number of random numbers to be generated.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The generated random numbers.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine DLARUV to generate random
*  real numbers from a uniform (0,1) distribution, in batches of up to
*  128 using vectorisable code. The Box-Muller method is used to
*  transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   U( LV )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARUV
*     ..
*     .. Executable Statements ..
*
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
*
*        Call DLARUV to generate IL2 numbers from a uniform (0,1)
*        distribution (IL2 <= LV)
*
         CALL DLARUV( ISEED, IL2, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           Copy generated numbers
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           Convert generated numbers to uniform (-1,1) distribution
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           Convert generated numbers to normal (0,1) distribution
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
*
*     End of DLARNV
*
      END
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*        FIRST = .FALSE.
*     END IF
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE.GE.SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
*
*     End of DLARTG
*
      END
      SUBROUTINE DLARUV( ISEED, N, X )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( N )
*     ..
*
*  Purpose
*  =======
*
*  DLARUV returns a vector of n random real numbers from a uniform (0,1)
*  distribution (n <= 128).
*
*  This is an auxiliary routine called by DLARNV and ZLARNV.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  N       (input) INTEGER
*          The number of random numbers to be generated. N <= 128.
*
*  X       (output) DOUBLE PRECISION array, dimension (N)
*          The generated random numbers.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      INTEGER            LV, IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
*     ..
*     .. Local Arrays ..
      INTEGER            MM( LV, 4 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN, MOD
*     ..
*     .. Data statements ..
      DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508,
     $                   2549 /
      DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754,
     $                   1145 /
      DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766,
     $                   2253 /
      DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572,
     $                   305 /
      DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893,
     $                   3301 /
      DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307,
     $                   1065 /
      DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297,
     $                   3133 /
      DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966,
     $                   2913 /
      DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758,
     $                   3285 /
      DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598,
     $                   1241 /
      DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406,
     $                   1197 /
      DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922,
     $                   3729 /
      DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038,
     $                   2501 /
      DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934,
     $                   1673 /
      DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091,
     $                   541 /
      DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451,
     $                   2753 /
      DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580,
     $                   949 /
      DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958,
     $                   2361 /
      DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055,
     $                   1165 /
      DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507,
     $                   4081 /
      DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078,
     $                   2725 /
      DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273,
     $                   3305 /
      DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17,
     $                   3069 /
      DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854,
     $                   3617 /
      DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916,
     $                   3733 /
      DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971,
     $                   409 /
      DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889,
     $                   2157 /
      DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831,
     $                   1361 /
      DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621,
     $                   3973 /
      DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541,
     $                   1865 /
      DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893,
     $                   2525 /
      DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736,
     $                   1409 /
      DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992,
     $                   3445 /
      DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787,
     $                   3577 /
      DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125,
     $                   77 /
      DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364,
     $                   3761 /
      DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460,
     $                   2149 /
      DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257,
     $                   1449 /
      DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574,
     $                   3005 /
      DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912,
     $                   225 /
      DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216,
     $                   85 /
      DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248,
     $                   3673 /
      DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401,
     $                   3117 /
      DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124,
     $                   3089 /
      DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762,
     $                   1349 /
      DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149,
     $                   2057 /
      DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245,
     $                   413 /
      DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166,
     $                   65 /
      DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466,
     $                   1845 /
      DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018,
     $                   697 /
      DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399,
     $                   3085 /
      DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190,
     $                   3441 /
      DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879,
     $                   1573 /
      DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153,
     $                   3689 /
      DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320,
     $                   2941 /
      DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18,
     $                   929 /
      DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712,
     $                   533 /
      DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159,
     $                   2841 /
      DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318,
     $                   4077 /
      DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091,
     $                   721 /
      DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443,
     $                   2821 /
      DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510,
     $                   2249 /
      DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449,
     $                   2397 /
      DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956,
     $                   2817 /
      DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201,
     $                   245 /
      DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137,
     $                   1913 /
      DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399,
     $                   1997 /
      DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321,
     $                   3121 /
      DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271,
     $                   997 /
      DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667,
     $                   1833 /
      DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703,
     $                   2877 /
      DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629,
     $                   1633 /
      DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365,
     $                   981 /
      DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431,
     $                   2009 /
      DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113,
     $                   941 /
      DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922,
     $                   2449 /
      DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554,
     $                   197 /
      DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184,
     $                   2441 /
      DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099,
     $                   285 /
      DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228,
     $                   1473 /
      DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012,
     $                   2741 /
      DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921,
     $                   3129 /
      DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452,
     $                   909 /
      DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901,
     $                   2801 /
      DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572,
     $                   421 /
      DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309,
     $                   4073 /
      DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171,
     $                   2813 /
      DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817,
     $                   2337 /
      DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039,
     $                   1429 /
      DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696,
     $                   1177 /
      DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256,
     $                   1901 /
      DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715,
     $                   81 /
      DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077,
     $                   1669 /
      DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019,
     $                   2633 /
      DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497,
     $                   2269 /
      DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101,
     $                   129 /
      DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717,
     $                   1141 /
      DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51,
     $                   249 /
      DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981,
     $                   3917 /
      DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978,
     $                   2481 /
      DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813,
     $                   3941 /
      DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881,
     $                   2217 /
      DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76,
     $                   2749 /
      DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846,
     $                   3041 /
      DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694,
     $                   1877 /
      DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682,
     $                   345 /
      DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124,
     $                   2861 /
      DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660,
     $                   1809 /
      DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997,
     $                   3141 /
      DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479,
     $                   2825 /
      DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141,
     $                   157 /
      DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886,
     $                   2881 /
      DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514,
     $                   3637 /
      DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301,
     $                   1465 /
      DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604,
     $                   2829 /
      DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888,
     $                   2161 /
      DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836,
     $                   3365 /
      DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990,
     $                   361 /
      DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058,
     $                   2685 /
      DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692,
     $                   3745 /
      DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194,
     $                   2325 /
      DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20,
     $                   3609 /
      DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285,
     $                   3821 /
      DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046,
     $                   3537 /
      DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107,
     $                   517 /
      DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508,
     $                   3017 /
      DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525,
     $                   2141 /
      DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801,
     $                   1537 /
*     ..
*     .. Executable Statements ..
*
      I1 = ISEED( 1 )
      I2 = ISEED( 2 )
      I3 = ISEED( 3 )
      I4 = ISEED( 4 )
*
      DO 10 I = 1, MIN( N, LV )
*	  
  20     CONTINUE
*
*        Multiply the seed by i-th power of the multiplier modulo 2**48
*
         IT4 = I4*MM( I, 4 )
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM( I, 4 ) + I4*MM( I, 3 )
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM( I, 4 ) + I3*MM( I, 3 ) + I4*MM( I, 2 )
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM( I, 4 ) + I2*MM( I, 3 ) + I3*MM( I, 2 ) +
     $         I4*MM( I, 1 )
         IT1 = MOD( IT1, IPW2 )
*
*        Convert 48-bit integer to a real number in the interval (0,1)
*
         X( I ) = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $            DBLE( IT4 ) ) ) )
*
         IF (X( I ).EQ.1.0D0) THEN
*           If a real number has n bits of precision, and the first
*           n bits of the 48-bit integer above happen to be all 1 (which
*           will occur about once every 2**n calls), then X( I ) will
*           be rounded to exactly 1.0. 
*           Since X( I ) is not supposed to return exactly 0.0 or 1.0,
*           the statistically correct thing to do in this situation is
*           simply to iterate again.
*           N.B. the case X( I ) = 0.0 should not be possible.	
            I1 = I1 + 2
            I2 = I2 + 2
            I3 = I3 + 2
            I4 = I4 + 2
            GOTO 20
         END IF
*
   10 CONTINUE
*
*     Return final value of seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RETURN
*
*     End of DLARUV
*
      END
      SUBROUTINE DLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH, DISNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) THEN
         INFO = -4
      ELSE IF( DISNAN(CTO) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      IF( CFROM1.EQ.CFROMC ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = .TRUE.
         CTO1 = CTOC
      ELSE
         CTO1 = CTOC / BIGNUM
         IF( CTO1.EQ.CTOC ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            MUL = CTOC
            DONE = .TRUE.
            CFROMC = ONE
         ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
            MUL = SMLNUM
            DONE = .FALSE.
            CFROMC = CFROM1
         ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
            MUL = BIGNUM
            DONE = .FALSE.
            CTOC = CTO1
         ELSE
            MUL = CTOC / CFROMC
            DONE = .TRUE.
         END IF
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DLASCL
*
      END
      SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of A is not changed.
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of A is not changed.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) DOUBLE PRECISION
*          The constant to which the diagonal elements are to be set.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the leading m-by-n submatrix of A is set as follows:
*
*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*
*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the strictly upper triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the strictly lower triangular or trapezoidal part of the
*        array to ALPHA.
*
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
*
      ELSE
*
*        Set the leading m-by-n submatrix to ALPHA.
*
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
*
*     Set the first min(M,N) diagonal elements to BETA.
*
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
*
      RETURN
*
*     End of DLASET
*
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE DLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR,
     $                   LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      DOUBLE PRECISION   SCALE, XNORM
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ),
     $                   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
*
*         op(TL)*X + ISGN*X*op(TR) = SCALE*B,
*
*  where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
*  -1.  op(T) = T or T', where T' denotes the transpose of T.
*
*  Arguments
*  =========
*
*  LTRANL  (input) LOGICAL
*          On entry, LTRANL specifies the op(TL):
*             = .FALSE., op(TL) = TL,
*             = .TRUE., op(TL) = TL'.
*
*  LTRANR  (input) LOGICAL
*          On entry, LTRANR specifies the op(TR):
*            = .FALSE., op(TR) = TR,
*            = .TRUE., op(TR) = TR'.
*
*  ISGN    (input) INTEGER
*          On entry, ISGN specifies the sign of the equation
*          as described before. ISGN may only be 1 or -1.
*
*  N1      (input) INTEGER
*          On entry, N1 specifies the order of matrix TL.
*          N1 may only be 0, 1 or 2.
*
*  N2      (input) INTEGER
*          On entry, N2 specifies the order of matrix TR.
*          N2 may only be 0, 1 or 2.
*
*  TL      (input) DOUBLE PRECISION array, dimension (LDTL,2)
*          On entry, TL contains an N1 by N1 matrix.
*
*  LDTL    (input) INTEGER
*          The leading dimension of the matrix TL. LDTL >= max(1,N1).
*
*  TR      (input) DOUBLE PRECISION array, dimension (LDTR,2)
*          On entry, TR contains an N2 by N2 matrix.
*
*  LDTR    (input) INTEGER
*          The leading dimension of the matrix TR. LDTR >= max(1,N2).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,2)
*          On entry, the N1 by N2 matrix B contains the right-hand
*          side of the equation.
*
*  LDB     (input) INTEGER
*          The leading dimension of the matrix B. LDB >= max(1,N1).
*
*  SCALE   (output) DOUBLE PRECISION
*          On exit, SCALE contains the scale factor. SCALE is chosen
*          less than or equal to 1 to prevent the solution overflowing.
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,2)
*          On exit, X contains the N1 by N2 solution.
*
*  LDX     (input) INTEGER
*          The leading dimension of the matrix X. LDX >= max(1,N1).
*
*  XNORM   (output) DOUBLE PRECISION
*          On exit, XNORM is the infinity-norm of the solution.
*
*  INFO    (output) INTEGER
*          On exit, INFO is set to
*             0: successful exit.
*             1: TL and TR have too close eigenvalues, so TL or
*                TR is perturbed to get a nonsingular equation.
*          NOTE: In the interests of speed, this routine does not
*                check the inputs for errors.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO, HALF, EIGHT
      PARAMETER          ( TWO = 2.0D+0, HALF = 0.5D+0, EIGHT = 8.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      DOUBLE PRECISION   BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1,
     $                   TEMP, U11, U12, U22, XMAX
*     ..
*     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ),
     $                   LOCU22( 4 )
      DOUBLE PRECISION   BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Data statements ..
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / ,
     $                   LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     Do not check the input parameters for errors
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N1.EQ.0 .OR. N2.EQ.0 )
     $   RETURN
*
*     Set constants to control overflow
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SGN = ISGN
*
      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K
*
*     1 by 1: TL11*X + SGN*X*TR11 = B11
*
   10 CONTINUE
      TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
      BET = ABS( TAU1 )
      IF( BET.LE.SMLNUM ) THEN
         TAU1 = SMLNUM
         BET = SMLNUM
         INFO = 1
      END IF
*
      SCALE = ONE
      GAM = ABS( B( 1, 1 ) )
      IF( SMLNUM*GAM.GT.BET )
     $   SCALE = ONE / GAM
*
      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN
*
*     1 by 2:
*     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
*                                       [TR21 TR22]
*
   20 CONTINUE
*
      SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ),
     $       ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ),
     $       SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      IF( LTRANR ) THEN
         TMP( 2 ) = SGN*TR( 2, 1 )
         TMP( 3 ) = SGN*TR( 1, 2 )
      ELSE
         TMP( 2 ) = SGN*TR( 1, 2 )
         TMP( 3 ) = SGN*TR( 2, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40
*
*     2 by 1:
*          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
*            [TL21 TL22] [X21]         [X21]         [B21]
*
   30 CONTINUE
      SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ),
     $       ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ),
     $       SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      IF( LTRANL ) THEN
         TMP( 2 ) = TL( 1, 2 )
         TMP( 3 ) = TL( 2, 1 )
      ELSE
         TMP( 2 ) = TL( 2, 1 )
         TMP( 3 ) = TL( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
   40 CONTINUE
*
*     Solve 2 by 2 system using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*
      IPIV = IDAMAX( 4, TMP, 1 )
      U11 = TMP( IPIV )
      IF( ABS( U11 ).LE.SMIN ) THEN
         INFO = 1
         U11 = SMIN
      END IF
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      IF( ABS( U22 ).LE.SMIN ) THEN
         INFO = 1
         U22 = SMIN
      END IF
      IF( BSWAP ) THEN
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      ELSE
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      END IF
      SCALE = ONE
      IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR.
     $    ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      END IF
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      IF( XSWAP ) THEN
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      END IF
      X( 1, 1 ) = X2( 1 )
      IF( N1.EQ.1 ) THEN
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
      ELSE
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
      END IF
      RETURN
*
*     2 by 2:
*     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
*       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
*
*     Solve equivalent 4 by 4 system using complete pivoting.
*     Set pivots less than SMIN to SMIN.
*
   50 CONTINUE
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ),
     $       ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
      SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ),
     $       ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      BTMP( 1 ) = ZERO
      CALL DCOPY( 16, BTMP, 0, T16, 1 )
      T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
      IF( LTRANL ) THEN
         T16( 1, 2 ) = TL( 2, 1 )
         T16( 2, 1 ) = TL( 1, 2 )
         T16( 3, 4 ) = TL( 2, 1 )
         T16( 4, 3 ) = TL( 1, 2 )
      ELSE
         T16( 1, 2 ) = TL( 1, 2 )
         T16( 2, 1 ) = TL( 2, 1 )
         T16( 3, 4 ) = TL( 1, 2 )
         T16( 4, 3 ) = TL( 2, 1 )
      END IF
      IF( LTRANR ) THEN
         T16( 1, 3 ) = SGN*TR( 1, 2 )
         T16( 2, 4 ) = SGN*TR( 1, 2 )
         T16( 3, 1 ) = SGN*TR( 2, 1 )
         T16( 4, 2 ) = SGN*TR( 2, 1 )
      ELSE
         T16( 1, 3 ) = SGN*TR( 2, 1 )
         T16( 2, 4 ) = SGN*TR( 2, 1 )
         T16( 3, 1 ) = SGN*TR( 1, 2 )
         T16( 4, 2 ) = SGN*TR( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )
*
*     Perform elimination
*
      DO 100 I = 1, 3
         XMAX = ZERO
         DO 70 IP = I, 4
            DO 60 JP = I, 4
               IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   60       CONTINUE
   70    CONTINUE
         IF( IPSV.NE.I ) THEN
            CALL DSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV.NE.I )
     $      CALL DSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
            INFO = 1
            T16( I, I ) = SMIN
         END IF
         DO 90 J = I + 1, 4
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
            DO 80 K = I + 1, 4
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      IF( ABS( T16( 4, 4 ) ).LT.SMIN )
     $   T16( 4, 4 ) = SMIN
      SCALE = ONE
      IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR.
     $    ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
         SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ),
     $           ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
         BTMP( 4 ) = BTMP( 4 )*SCALE
      END IF
      DO 120 I = 1, 4
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 110 J = K + 1, 4
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
  110    CONTINUE
  120 CONTINUE
      DO 130 I = 1, 3
         IF( JPIV( 4-I ).NE.4-I ) THEN
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         END IF
  130 CONTINUE
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ),
     $        ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
      RETURN
*
*     End of DLASY2
*
      END
      SUBROUTINE DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
     $                   CNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), CNORM( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLATRS solves one of the triangular systems
*
*     A *x = s*b  or  A'*x = s*b
*
*  with scaling to prevent overflow.  Here A is an upper or lower
*  triangular matrix, A' denotes the transpose of A, x and b are
*  n-element vectors, and s is a scaling factor, usually less than
*  or equal to 1, chosen so that the components of x will be less than
*  the overflow threshold.  If the unscaled problem will not cause
*  overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A
*  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
*  non-trivial solution to A*x = 0 is returned.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  TRANS   (input) CHARACTER*1
*          Specifies the operation applied to A.
*          = 'N':  Solve A * x = s*b  (No transpose)
*          = 'T':  Solve A'* x = s*b  (Transpose)
*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  NORMIN  (input) CHARACTER*1
*          Specifies whether CNORM has been set or not.
*          = 'Y':  CNORM contains the column norms on entry
*          = 'N':  CNORM is not set on entry.  On exit, the norms will
*                  be computed and stored in CNORM.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of the array A contains the upper
*          triangular matrix, and the strictly lower triangular part of
*          A is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of the array A contains the lower triangular
*          matrix, and the strictly upper triangular part of A is not
*          referenced.  If DIAG = 'U', the diagonal elements of A are
*          also not referenced and are assumed to be 1.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max (1,N).
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the right hand side b of the triangular system.
*          On exit, X is overwritten by the solution vector x.
*
*  SCALE   (output) DOUBLE PRECISION
*          The scaling factor s for the triangular system
*             A * x = s*b  or  A'* x = s*b.
*          If SCALE = 0, the matrix A is singular or badly scaled, and
*          the vector x is an exact or approximate solution to A*x = 0.
*
*  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
*
*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*          contains the norm of the off-diagonal part of the j-th column
*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*          must be greater than or equal to the 1-norm.
*
*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*          returns the 1-norm of the offdiagonal part of the j-th column
*          of A.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*
*  Further Details
*  ======= =======
*
*  A rough bound on x is computed; if that is less than overflow, DTRSV
*  is called, otherwise, specific code is used which checks for possible
*  overflow or divide-by-zero at every operation.
*
*  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*  if A is lower triangular is
*
*       x[1:n] := b[1:n]
*       for j = 1, ..., n
*            x(j) := x(j) / A(j,j)
*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*       end
*
*  Define bounds on the components of x after j iterations of the loop:
*     M(j) = bound on x[1:j]
*     G(j) = bound on x[j+1:n]
*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*
*  Then for iteration j+1 we have
*     M(j+1) <= G(j) / | A(j+1,j+1) |
*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*
*  where CNORM(j+1) is greater than or equal to the infinity-norm of
*  column j+1 of A, not counting the diagonal.  Hence
*
*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*                  1<=i<=j
*  and
*
*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*                                   1<=i< j
*
*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the
*  reciprocal of the largest M(j), j=1,..,n, is larger than
*  max(underflow, 1/overflow).
*
*  The bound on x(j) is also used to determine when a step in the
*  columnwise method can be performed without fear of overflow.  If
*  the computed bound is greater than a large constant, x is scaled to
*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*
*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
*  algorithm for A upper triangular is
*
*       for j = 1, ..., n
*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*       end
*
*  We simultaneously compute two bounds
*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*       M(j) = bound on x(i), 1<=i<=j
*
*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*  Then the bound on x(j) is
*
*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*
*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*                      1<=i<=j
*
*  and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater
*  than max(underflow, 1/overflow).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS,
     $                   TMAX, TSCAL, USCAL, XBND, XJ, XMAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DASUM, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DTRSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT.
     $         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
*
      IF( LSAME( NORMIN, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( UPPER ) THEN
*
*           A is upper triangular.
*
            DO 10 J = 1, N
               CNORM( J ) = DASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N - 1
               CNORM( J ) = DASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM.
*
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine DTRSV can be used.
*
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              M(j) = G(j-1) / abs(A(j,j))
*
               TJJ = ABS( A( J, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
*
      ELSE
*
*        Compute the growth in A' * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
               TJJ = ABS( A( J, J ) )
               IF( XJ.GT.TJJ )
     $            XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
*
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL DTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 110 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = A( J, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE )
     $               GO TO 100
               END IF
               TJJ = ABS( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by 1/b(j).
*
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        REC = REC / CNORM( J )
                     END IF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  100          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     CALL DAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X,
     $                           1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     CALL DAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1,
     $                           X( J+1 ), 1 )
                     I = J + IDAMAX( N-J, X( J+1 ), 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               END IF
  110       CONTINUE
*
         ELSE
*
*           Solve A' * x = b
*
            DO 160 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               SUMJ = ZERO
               IF( USCAL.EQ.ONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call DDOT to perform the dot product.
*
                  IF( UPPER ) THEN
                     SUMJ = DDOT( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     SUMJ = DDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     DO 120 I = 1, J - 1
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  120                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 130 I = J + 1, N
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.TSCAL ) THEN
*
*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE )
     $                  GO TO 150
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0, and compute a solution to A'*x = 0.
*
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  150             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
*
      RETURN
*
*     End of DLATRS
*
      END
*-----------------------------------------------------------------------
*  Routine:    DMOUT
*
*  Purpose:    Real matrix output routine.
*
*  Usage:      CALL DMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
*
*  Arguments
*     M      - Number of rows of A.  (Input)
*     N      - Number of columns of A.  (Input)
*     A      - Real M by N matrix to be printed.  (Input)
*     LDA    - Leading dimension of A exactly as specified in the
*              dimension statement of the calling program.  (Input)
*     IFMT   - Format to be used in printing matrix A.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*-----------------------------------------------------------------------
*
      SUBROUTINE DMOUT( LOUT, M, N, A, LDA, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
*     .. Scalar Arguments ..
      CHARACTER*( * )    IFMT
      INTEGER            IDIGIT, LDA, LOUT, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*     .. Local Scalars ..
      CHARACTER*80       LINE
      INTEGER            I, J, K1, K2, LLL, NDIGIT
*     ..
*     .. Local Arrays ..
      CHARACTER          ICOL( 3 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MIN, MIN0
*     ..
*     .. Data statements ..
      DATA               ICOL( 1 ), ICOL( 2 ), ICOL( 3 ) / 'C', 'o',
     $                   'l' /
*     ..
*     .. Executable Statements ..
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
      LLL = MIN( LEN( IFMT ), 80 )
      DO 10 I = 1, LLL
         LINE( I: I ) = '-'
   10 CONTINUE
*
      DO 20 I = LLL + 1, 80
         LINE( I: I ) = ' '
   20 CONTINUE
*
      WRITE( LOUT, FMT = 9999 )IFMT, LINE( 1: LLL )
 9999 FORMAT( / 1X, A, / 1X, A )
*
      IF( M.LE.0 .OR. N.LE.0 .OR. LDA.LE.0 )
     $   RETURN
      NDIGIT = IDIGIT
      IF( IDIGIT.EQ.0 )
     $   NDIGIT = 4
*
*=======================================================================
*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
*=======================================================================
*
      IF( IDIGIT.LT.0 ) THEN
         NDIGIT = -IDIGIT
         IF( NDIGIT.LE.4 ) THEN
            DO 40 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, FMT = 9998 )( ICOL, I, I = K1, K2 )
               DO 30 I = 1, M
                  WRITE( LOUT, FMT = 9994 )I, ( A( I, J ), J = K1, K2 )
   30          CONTINUE
   40       CONTINUE
*
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 60 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               WRITE( LOUT, FMT = 9997 )( ICOL, I, I = K1, K2 )
               DO 50 I = 1, M
                  WRITE( LOUT, FMT = 9993 )I, ( A( I, J ), J = K1, K2 )
   50          CONTINUE
   60       CONTINUE
*
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 80 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               WRITE( LOUT, FMT = 9996 )( ICOL, I, I = K1, K2 )
               DO 70 I = 1, M
                  WRITE( LOUT, FMT = 9992 )I, ( A( I, J ), J = K1, K2 )
   70          CONTINUE
   80       CONTINUE
*
         ELSE
            DO 100 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, FMT = 9995 )( ICOL, I, I = K1, K2 )
               DO 90 I = 1, M
                  WRITE( LOUT, FMT = 9991 )I, ( A( I, J ), J = K1, K2 )
   90          CONTINUE
  100       CONTINUE
         END IF
*
*=======================================================================
*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
*=======================================================================
*
      ELSE
         IF( NDIGIT.LE.4 ) THEN
            DO 120 K1 = 1, N, 10
               K2 = MIN0( N, K1+9 )
               WRITE( LOUT, FMT = 9998 )( ICOL, I, I = K1, K2 )
               DO 110 I = 1, M
                  WRITE( LOUT, FMT = 9994 )I, ( A( I, J ), J = K1, K2 )
  110          CONTINUE
  120       CONTINUE
*
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 140 K1 = 1, N, 8
               K2 = MIN0( N, K1+7 )
               WRITE( LOUT, FMT = 9997 )( ICOL, I, I = K1, K2 )
               DO 130 I = 1, M
                  WRITE( LOUT, FMT = 9993 )I, ( A( I, J ), J = K1, K2 )
  130          CONTINUE
  140       CONTINUE
*
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 160 K1 = 1, N, 6
               K2 = MIN0( N, K1+5 )
               WRITE( LOUT, FMT = 9996 )( ICOL, I, I = K1, K2 )
               DO 150 I = 1, M
                  WRITE( LOUT, FMT = 9992 )I, ( A( I, J ), J = K1, K2 )
  150          CONTINUE
  160       CONTINUE
*
         ELSE
            DO 180 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, FMT = 9995 )( ICOL, I, I = K1, K2 )
               DO 170 I = 1, M
                  WRITE( LOUT, FMT = 9991 )I, ( A( I, J ), J = K1, K2 )
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
      WRITE( LOUT, FMT = 9990 )
*
 9998 FORMAT( 10X, 10( 4X, 3A1, I4, 1X ) )
 9997 FORMAT( 10X, 8( 5X, 3A1, I4, 2X ) )
 9996 FORMAT( 10X, 6( 7X, 3A1, I4, 4X ) )
 9995 FORMAT( 10X, 5( 9X, 3A1, I4, 6X ) )
 9994 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 10D12.3 )
 9993 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 8D14.5 )
 9992 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 6D18.9 )
 9991 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 5D22.13 )
 9990 FORMAT( 1X, ' ' )
*
      RETURN
      END
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnaitr
c
c\Description: 
c  Reverse communication interface for applying NP additional steps to 
c  a K step nonsymmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dnaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dnaitr
c     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH, 
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and do not need to be
c          recompute in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.  See dnaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current size of V and H.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  NB      Integer.  (INPUT)
c          Blocksize to be used in the recurrence.          
c          Only work for NB = 1 right now.  The goal is to have a 
c          program that implement both the block and non-block method.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          B-norm of the starting residual on input.
c          B-norm of the updated residual r_{k+p} on output.
c
c  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K 
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
c          H is used to store the generated upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not 
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On input, WORKD(1:N) = B*RESID and is used to save some 
c          computation at the first step.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of the spanning invariant subspace of OP found.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlabad  LAPACK routine that computes machine constants.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c 
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c  
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        ( At present tol is zero )
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dnaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        H(:,j) = w_{j};
c        H(j,j-1) = rnorm
c        rnorm = || r_(j) ||
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};   
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf 
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaitr
     &   (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, 
     &    ipntr, workd, info)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, nb, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      logical    first, orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, msglvl,
     &           jj
      Double precision
     &           betaj, ovfl, temp1, rnorm1, smlnum, tst1, ulp, unfl, 
     &           wnorm
      save       first, orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl, ovfl,
     &           betaj, rnorm1, smlnum, ulp, unfl, wnorm
c
c     %-----------------------%
c     | Local Array Arguments | 
c     %-----------------------%
c
      Double precision
     &           xtemp(2)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dgemv, dgetv0, dlabad, 
     &           dvout, dmout, ivout, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlanhs, dlamch
      external   ddot, dnrm2, dlanhs, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, sqrt
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------%
c        | Set machine-dependent constants for the |
c        | the splitting and deflation criterion.  |
c        | If norm(H) <= sqrt(OVFL),               |
c        | overflow should not occur.              |
c        | REFERENCE: LAPACK subroutine dlahqr     |
c        %-----------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = one / unfl
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call second (t0)
         msglvl = mnaitr
c 
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
         j      = k + 1
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c 
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true. when ....                        |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %-----------------------------%
c     | Else this is the first step |
c     %-----------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
 
 1000 continue
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, j, ndigit, 
     &                  '_naitr: generating Arnoldi vector number')
            call dvout (logfil, 1, rnorm, ndigit, 
     &                  '_naitr: B-norm of the current residual is')
         end if
c 
c        %---------------------------------------------------%
c        | STEP 1: Check if the B norm of j-th residual      |
c        | vector is zero. Equivalent to determing whether   |
c        | an exact j-step Arnoldi factorization is present. |
c        %---------------------------------------------------%
c
         betaj = rnorm
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
            if (msglvl .gt. 0) then
               call ivout (logfil, 1, j, ndigit,
     &                     '_naitr: ****** RESTART AT STEP ******')
            end if
c 
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c 
            betaj  = zero
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call dgetv0 (ido, bmat, itry, .false., n, j, v, ldv, 
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               call second (t1)
               tnaitr = tnaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
c 
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. unfl) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), 1)
             call dscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call dlascl ('General', i, i, rnorm, one, n, 1, 
     &                    v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1, 
     &                    workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call second (t2)
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c 
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c 
         go to 9000 
   50    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
c        | if step3 = .true.                |
c        %----------------------------------%
c
         call second (t3)
         tmvopx = tmvopx + (t3 - t2)
 
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call dcopy (n, workd(irj), 1, resid, 1)
c 
c        %---------------------------------------%
c        | STEP 4:  Finish extending the Arnoldi |
c        |          factorization to length j.   |
c        %---------------------------------------%
c
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
c        | if step4 = .true.                |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
         if (bmat .eq. 'G') then  
             wnorm = ddot (n, resid, 1, workd(ipj), 1)
             wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c 
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, h(1,j), 1)
c
c        %--------------------------------------%
c        | Orthogonalize r_{j} against V_{j}.   |
c        | RESID contains OP*v_{j}. See STEP 3. | 
c        %--------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, h(1,j), 1,
     &               one, resid, 1)
c
         if (j .gt. 1) h(j,j-1) = betaj
c
         call second (t4)
c 
         orth1 = .true.
c
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if 
   70    continue
c 
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
c 
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        | The following test determines whether the sine of the     |
c        | angle between  OP*x and the computed residual is less     |
c        | than or equal to 0.717.                                   |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         iter  = 0
         nrorth = nrorth + 1
c 
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c 
   80    continue
c
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call dvout (logfil, 2, xtemp, ndigit, 
     &           '_naitr: re-orthonalization; wnorm and rnorm are')
            call dvout (logfil, j, h(1,j), ndigit,
     &                  '_naitr: j-th column of H')
         end if
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1, 
     &               zero, workd(irj), 1)
c
c        %---------------------------------------------%
c        | Compute the correction to the residual:     |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
c        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
c        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
c        %---------------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1, 
     &               one, resid, 1)
         call daxpy (j, one, workd(irj), 1, h(1,j), 1)
c 
         orth2 = .true.
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c 
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if 
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c 
         if (bmat .eq. 'G') then         
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
c
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            call ivout (logfil, 1, j, ndigit,
     &           '_naitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call dvout (logfil, 2, xtemp, ndigit,
     &           '_naitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %---------------------------------------%
c           | No need for further refinement.       |
c           | The cosine of the angle between the   |
c           | corrected residual vector and the old |
c           | residual vector is greater than 0.717 |
c           | In other words the corrected residual |
c           | and the old residual vector share an  |
c           | angle of less than arcCOS(0.717)      |
c           %---------------------------------------%
c
            rnorm = rnorm1
c 
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c 
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c 
  100    continue
c 
         rstart = .false.
         orth2  = .false.
c 
         call second (t5)
         titref = titref + (t5 - t4)
c 
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            call second (t1)
            tnaitr = tnaitr + (t1 - t0)
            ido = 99
            do 110 i = max(1,k), k+np-1
c     
c              %--------------------------------------------%
c              | Check for splitting and deflation.         |
c              | Use a standard test as in the QR algorithm |
c              | REFERENCE: LAPACK subroutine dlahqr        |
c              %--------------------------------------------%
c     
               tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
               if( tst1.eq.zero )
     &              tst1 = dlanhs( '1', k+np, h, ldh, workd(n+1) )
               if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) 
     &              h(i+1,i) = zero
 110        continue
c     
            if (msglvl .gt. 2) then
               call dmout (logfil, k+np, k+np, h, ldh, ndigit, 
     &          '_naitr: Final upper Hessenberg matrix H of order K+NP')
            end if
c     
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c 
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of dnaitr |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge chage sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call dnapps
c     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ, 
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW.  KEV is only 
c          updated on ouput when fewer than NP shifts are applied in
c          order to keep the conjugate pair together.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFTR, Double precision array of length NP.  (INPUT)
c  SHIFTI  Real and imaginary part of the shifts to be applied.
c          Upon, entry to dnapps, the shifts must be sorted so that the 
c          conjugate pairs are in consecutive locations.
c
c  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper 
c          Hessenber matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k} 
c          in the first KEV locations.
c
c  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Double precision work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices.
c     dvout   ARPACK utility routine that prints vectors.
c     dlabad  LAPACK routine that computes machine constants.
c     dlacpy  LAPACK matrix copy routine.
c     dlamch  LAPACK routine that determines machine constants. 
c     dlanhs  LAPACK routine that computes various norms of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlarf   LAPACK routine that applies Householder reflection to
c             a matrix.
c     dlarfg  LAPACK Householder reflection construction routine.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlaset  LAPACK matrix initialization routine.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dscal   Level 1 BLAS that scales a vector.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine dlahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     The subdiagonals of H are enforced to be non-negative.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnapps
     &   ( n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq, 
     &     workl, workd )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           h(ldh,kev+np), resid(n), shifti(np), shiftr(np), 
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, ir, istart, j, jj, kplusp, msglvl, nr
      logical    cconj, first
      Double precision
     &           c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai, 
     &           sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
      save       first, ovfl, smlnum, ulp, unfl 
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dlacpy, dlarfg, dlarf,
     &           dlaset, dlabad, second, dlartg
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch, dlanhs, dlapy2
      external   dlamch, dlanhs, dlapy2
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, max, min
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine dlahqr           |
c        %-----------------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = one / unfl
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call second (t0)
      msglvl = mnapps
      kplusp = kev + np 
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call dlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Quick return if there are no shifts to apply |
c     %----------------------------------------------%
c
      if (np .eq. 0) go to 9000
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      cconj = .false.
      do 110 jj = 1, np
         sigmar = shiftr(jj)
         sigmai = shifti(jj)
c
         if (msglvl .gt. 2 ) then
            call ivout (logfil, 1, jj, ndigit, 
     &               '_napps: shift number.')
            call dvout (logfil, 1, sigmar, ndigit, 
     &               '_napps: The real part of the shift ')
            call dvout (logfil, 1, sigmai, ndigit, 
     &               '_napps: The imaginary part of the shift ')
         end if
c
c        %-------------------------------------------------%
c        | The following set of conditionals is necessary  |
c        | in order that complex conjugate pairs of shifts |
c        | are applied together or not at all.             |
c        %-------------------------------------------------%
c
         if ( cconj ) then
c
c           %-----------------------------------------%
c           | cconj = .true. means the previous shift |
c           | had non-zero imaginary part.            |
c           %-----------------------------------------%
c
            cconj = .false.
            go to 110
         else if ( jj .lt. np .and. abs( sigmai ) .gt. zero ) then
c
c           %------------------------------------%
c           | Start of a complex conjugate pair. |
c           %------------------------------------%
c
            cconj = .true.
         else if ( jj .eq. np .and. abs( sigmai ) .gt. zero ) then
c
c           %----------------------------------------------%
c           | The last shift has a nonzero imaginary part. |
c           | Don't apply it; thus the order of the        |
c           | compressed H is order KEV+1 since only np-1  |
c           | were applied.                                |
c           %----------------------------------------------%
c
            kev = kev + 1
            go to 110
         end if
         istart = 1
   20    continue
c
c        %--------------------------------------------------%
c        | if sigmai = 0 then                               |
c        |    Apply the jj-th shift ...                     |
c        | else                                             |
c        |    Apply the jj-th and (jj+1)-th together ...    |
c        |    (Note that jj < np at this point in the code) |
c        | end                                              |
c        | to the current block of H. The next do loop      |
c        | determines the current block ;                   |
c        %--------------------------------------------------%
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine dlahqr    |
c           %----------------------------------------%
c
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
            if( tst1.eq.zero )
     &         tst1 = dlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) then
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, i, ndigit, 
     &                 '_napps: matrix splitting at row/column no.')
                  call ivout (logfil, 1, jj, ndigit, 
     &                 '_napps: matrix splitting with shift number.')
                  call dvout (logfil, 1, h(i+1,i), ndigit, 
     &                 '_napps: off diagonal element.')
               end if
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (msglvl .gt. 2) then
             call ivout (logfil, 1, istart, ndigit, 
     &                   '_napps: Start of current block ')
             call ivout (logfil, 1, iend, ndigit, 
     &                   '_napps: End of current block ')
         end if
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        %------------------------------------------------%
c
         if ( istart .eq. iend ) go to 100
c
c        %------------------------------------------------------%
c        | If istart + 1 = iend then no reason to apply a       |
c        | complex conjugate pair of shifts on a 2 by 2 matrix. |
c        %------------------------------------------------------%
c
         if ( istart + 1 .eq. iend .and. abs( sigmai ) .gt. zero ) 
     &      go to 100
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         if ( abs( sigmai ) .le. zero ) then
c
c           %---------------------------------------------%
c           | Real-valued shift ==> apply single shift QR |
c           %---------------------------------------------%
c
            f = h11 - sigmar
            g = h21
c 
            do 80 i = istart, iend-1
c
c              %-----------------------------------------------------%
c              | Contruct the plane rotation G to zero out the bulge |
c              %-----------------------------------------------------%
c
               call dlartg (f, g, c, s, r)
               if (i .gt. istart) then
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
                  if (r .lt. zero) then
                     r = -r
                     c = -c
                     s = -s
                  end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
c
c              %---------------------------------------------%
c              | Apply rotation to the left of H;  H <- G'*H |
c              %---------------------------------------------%
c
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t   
   50          continue
c
c              %---------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G |
c              %---------------------------------------------%
c
               do 60 j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + s*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t   
   60          continue
c
c              %----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
c              %----------------------------------------------------%
c
               do 70 j = 1, min( i+jj, kplusp ) 
                  t        =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t   
   70          continue
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
               if (i .lt. iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
c
c           %-----------------------------------%
c           | Finished applying the real shift. |
c           %-----------------------------------%
c 
         else
c
c           %----------------------------------------------------%
c           | Complex conjugate shifts ==> apply double shift QR |
c           %----------------------------------------------------%
c
            h12 = h(istart,istart+1)
            h22 = h(istart+1,istart+1)
            h32 = h(istart+2,istart+1)
c
c           %---------------------------------------------------------%
c           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
c           %---------------------------------------------------------%
c
            s    = 2.0*sigmar
            t = dlapy2 ( sigmar, sigmai ) 
            u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
            u(2) = h11 + h22 - s 
            u(3) = h32
c
            do 90 i = istart, iend-1
c
               nr = min ( 3, iend-i+1 )
c
c              %-----------------------------------------------------%
c              | Construct Householder reflector G to zero out u(1). |
c              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
c              %-----------------------------------------------------%
c
               call dlarfg ( nr, u(1), u(2), 1, tau )
c
               if (i .gt. istart) then
                  h(i,i-1)   = u(1)
                  h(i+1,i-1) = zero
                  if (i .lt. iend-1) h(i+2,i-1) = zero
               end if
               u(1) = one
c
c              %--------------------------------------%
c              | Apply the reflector to the left of H |
c              %--------------------------------------%
c
               call dlarf ('Left', nr, kplusp-i+1, u, 1, tau,
     &                     h(i,i), ldh, workl)
c
c              %---------------------------------------%
c              | Apply the reflector to the right of H |
c              %---------------------------------------%
c
               ir = min ( i+3, iend )
               call dlarf ('Right', ir, nr, u, 1, tau,
     &                     h(1,i), ldh, workl)
c
c              %-----------------------------------------------------%
c              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
c              %-----------------------------------------------------%
c
               call dlarf ('Right', kplusp, nr, u, 1, tau, 
     &                     q(1,i), ldq, workl)
c
c              %----------------------------%
c              | Prepare for next reflector |
c              %----------------------------%
c
               if (i .lt. iend-1) then
                  u(1) = h(i+1,i)
                  u(2) = h(i+2,i)
                  if (i .lt. iend-2) u(3) = h(i+3,i)
               end if
c
   90       continue
c
c           %--------------------------------------------%
c           | Finished applying a complex pair of shifts |
c           | to the current block                       |
c           %--------------------------------------------%
c 
         end if
c
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non negative sub diagonals |
c     %--------------------------------------------------%
c
      do 120 j=1,kev
         if ( h(j+1,j) .lt. zero ) then
              call dscal( kplusp-j+1, -one, h(j+1,j), ldh )
              call dscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
              call dscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
         end if
 120  continue
c
      do 130 i = 1, kev
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine dlahqr        |
c        %--------------------------------------------%
c
         tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
         if( tst1.eq.zero )
     &       tst1 = dlanhs( '1', kev, h, ldh, workl )
         if( h( i+1,i ) .le. max( ulp*tst1, smlnum ) ) 
     &       h(i+1,i) = zero
 130  continue
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &    call dgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero, 
     &                workd(n+1), 1)
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call dgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call dcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call dlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      if (h(kev+1,kev) .gt. zero)
     &   call dcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call dscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,kev) .gt. zero)
     &   call daxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
      if (msglvl .gt. 1) then
         call dvout (logfil, 1, q(kplusp,kev), ndigit,
     &        '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call dvout (logfil, 1, h(kev+1,kev), ndigit,
     &        '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call ivout (logfil, 1, kev, ndigit, 
     &               '_napps: Order of the final Hessenberg matrix ')
         if (msglvl .gt. 2) then
            call dmout (logfil, kev, kev, h, ldh, ndigit,
     &      '_napps: updated Hessenberg matrix H for next iteration')
         end if
c
      end if
c 
 9000 continue
      call second (t1)
      tnapps = tnapps + (t1 - t0)
c 
      return
c
c     %---------------%
c     | End of dnapps |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dnaup2
c
c\Description: 
c  Intermediate level interface called by dnaupd.
c
c\Usage:
c  call dnaup2
c     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
c       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS, 
c       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
c
c\Arguments
c
c  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd.
c  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd.
c
c  NP      Integer.  (INPUT/OUTPUT)
c          Contains the number of implicit shifts to apply during 
c          each Arnoldi iteration.  
c          If ISHIFT=1, NP is adjusted dynamically at each iteration 
c          to accelerate convergence and prevent stagnation.
c          This is also roughly equal to the number of matrix-vector 
c          products (involving the operator OP) per Arnoldi iteration.
c          The logic for adjusting is contained within the current
c          subroutine.
c          If ISHIFT=0, NP is the number of shifts the user needs
c          to provide via reverse comunication. 0 < NP < NCV-NEV.
c          NP may be less than NCV-NEV for two reasons. The first, is
c          to keep complex conjugate pairs of "wanted" Ritz values 
c          together. The second, is that a leading block of the current
c          upper Hessenberg matrix has split off and contains "unwanted"
c          Ritz values.
c          Upon termination of the IRA iteration, NP contains the number 
c          of "converged" wanted Ritz values.
c
c  IUPD    Integer.  (INPUT)
c          IUPD .EQ. 0: use explicit restart instead implicit update.
c          IUPD .NE. 0: use implicit update.
c
c  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT)
c          The Arnoldi basis vectors are returned in the first NEV 
c          columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling 
c          program.
c
c  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT)
c          H is used to store the generated upper Hessenberg matrix
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT)
c  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
c          imaginary) part of the computed Ritz values of OP.
c
c  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT)
c          BOUNDS(1:NEV) contain the error bounds corresponding to 
c          the computed Ritz values.
c          
c  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
c          Private (replicated) work array used to accumulate the
c          rotation in the shift application step.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length at least 
c          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  It is used in shifts calculation, shifts
c          application and convergence checking.
c
c          On exit, the last 3*(NEV+NP) locations of WORKL contain
c          the Ritz values (real,imaginary) and associated Ritz
c          estimates of the current Hessenberg matrix.  They are
c          listed in the same order as returned from dneigh.
c
c          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
c          of WORKL are used in reverse communication to hold the user 
c          supplied shifts.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD for 
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the 
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (WORKSPACE)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note in DNAUPD.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =     0: Normal return.
c          =     1: Maximum number of iterations taken.
c                   All possible eigenvalues of OP has been found.  
c                   NP returns the number of converged Ritz values.
c          =     2: No shifts could be applied.
c          =    -8: Error return from LAPACK eigenvalue calculation;
c                   This should never happen.
c          =    -9: Starting vector is zero.
c          = -9999: Could not build an Arnoldi factorization.
c                   Size that was built in returned in NP.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c
c\Routines called:
c     dgetv0  ARPACK initial vector generation routine. 
c     dnaitr  ARPACK Arnoldi factorization routine.
c     dnapps  ARPACK application of implicit shifts routine.
c     dnconv  ARPACK convergence of Ritz values routine.
c     dneigh  ARPACK compute Ritz values and error bounds routine.
c     dngets  ARPACK reorder Ritz values and error bounds routine.
c     dsortc  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors. 
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dswap   Level 1 BLAS that swaps two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c 
c\SCCS Information: @(#) 
c FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
c
c\Remarks
c     1. None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaup2
     &   ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, 
     &     ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds, 
     &     q, ldq, workl, ipntr, workd, info )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter,
     &           n, nev, np
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(13)
      Double precision
     &           bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n),
     &           ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np), 
     &           workd(3*n), workl( (nev+np)*(nev+np+3) )
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  wprime*2
      logical    cnorm , getv0, initv, update, ushift
      integer    ierr  , iter , j    , kplusp, msglvl, nconv, 
     &           nevbef, nev0 , np0  , nptemp, numcnv
      Double precision
     &           rnorm , temp , eps23
      save       cnorm , getv0, initv, update, ushift,
     &           rnorm , iter , eps23, kplusp, msglvl, nconv , 
     &           nevbef, nev0 , np0  , numcnv
c
c     %-----------------------%
c     | Local array arguments |
c     %-----------------------%
c
      integer    kp(4)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy , dgetv0, dnaitr, dnconv, dneigh, 
     &           dngets, dnapps, dvout , ivout , second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlapy2, dlamch
      external   ddot, dnrm2, dlapy2, dlamch
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    min, max, abs, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (ido .eq. 0) then
c 
         call second (t0)
c 
         msglvl = mnaup2
c 
c        %-------------------------------------%
c        | Get the machine dependent constant. |
c        %-------------------------------------%
c
         eps23 = dlamch('Epsilon-Machine')
         eps23 = eps23**(2.0D+0 / 3.0D+0)
c
         nev0   = nev
         np0    = np
c
c        %-------------------------------------%
c        | kplusp is the bound on the largest  |
c        |        Lanczos factorization built. |
c        | nconv is the current number of      |
c        |        "converged" eigenvlues.      |
c        | iter is the counter on the current  |
c        |      iteration step.                |
c        %-------------------------------------%
c
         kplusp = nev + np
         nconv  = 0
         iter   = 0
c 
c        %---------------------------------------%
c        | Set flags for computing the first NEV |
c        | steps of the Arnoldi factorization.   |
c        %---------------------------------------%
c
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
c
         if (info .ne. 0) then
c
c           %--------------------------------------------%
c           | User provides the initial residual vector. |
c           %--------------------------------------------%
c
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
c 
c     %---------------------------------------------%
c     | Get a possibly random starting vector and   |
c     | force it into the range of the operator OP. |
c     %---------------------------------------------%
c
   10 continue
c
      if (getv0) then
         call dgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm,
     &                ipntr, workd, info)
c
         if (ido .ne. 99) go to 9000
c
         if (rnorm .eq. zero) then
c
c           %-----------------------------------------%
c           | The initial vector is zero. Error exit. | 
c           %-----------------------------------------%
c
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
c 
c     %-----------------------------------%
c     | Back from reverse communication : |
c     | continue with update step         |
c     %-----------------------------------%
c
      if (update) go to 20
c
c     %-------------------------------------------%
c     | Back from computing user specified shifts |
c     %-------------------------------------------%
c
      if (ushift) go to 50
c
c     %-------------------------------------%
c     | Back from computing residual norm   |
c     | at the end of the current iteration |
c     %-------------------------------------%
c
      if (cnorm)  go to 100
c 
c     %----------------------------------------------------------%
c     | Compute the first NEV steps of the Arnoldi factorization |
c     %----------------------------------------------------------%
c
      call dnaitr (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, 
     &             h, ldh, ipntr, workd, info)
c 
c     %---------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication  |
c     | to compute operations involving OP and possibly B |
c     %---------------------------------------------------%
c
      if (ido .ne. 99) go to 9000
c
      if (info .gt. 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
c 
c     %--------------------------------------------------------------%
c     |                                                              |
c     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
c     |           Each iteration implicitly restarts the Arnoldi     |
c     |           factorization in place.                            |
c     |                                                              |
c     %--------------------------------------------------------------%
c 
 1000 continue
c
         iter = iter + 1
c
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, iter, ndigit, 
     &           '_naup2: **** Start of major iteration number ****')
         end if
c 
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        | Adjust NP since NEV might have been updated by last call  |
c        | to the shift application routine dnapps.                  |
c        %-----------------------------------------------------------%
c
         np  = kplusp - nev
c
         if (msglvl .gt. 1) then
            call ivout (logfil, 1, nev, ndigit, 
     &     '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, np, ndigit, 
     &           '_naup2: Extend the Arnoldi factorization by')
         end if
c
c        %-----------------------------------------------------------%
c        | Compute NP additional steps of the Arnoldi factorization. |
c        %-----------------------------------------------------------%
c
         ido = 0
   20    continue
         update = .true.
c
         call dnaitr (ido  , bmat, n  , nev, np , mode , resid, 
     &                rnorm, v   , ldv, h  , ldh, ipntr, workd,
     &                info)
c 
c        %---------------------------------------------------%
c        | ido .ne. 99 implies use of reverse communication  |
c        | to compute operations involving OP and possibly B |
c        %---------------------------------------------------%
c
         if (ido .ne. 99) go to 9000
c
         if (info .gt. 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.
c
         if (msglvl .gt. 1) then
            call dvout (logfil, 1, rnorm, ndigit, 
     &           '_naup2: Corresponding B-norm of the residual')
         end if
c 
c        %--------------------------------------------------------%
c        | Compute the eigenvalues and corresponding error bounds |
c        | of the current upper Hessenberg matrix.                |
c        %--------------------------------------------------------%
c
         call dneigh (rnorm, kplusp, h, ldh, ritzr, ritzi, bounds,
     &                q, ldq, workl, ierr)
c
         if (ierr .ne. 0) then
            info = -8
            go to 1200
         end if
c
c        %----------------------------------------------------%
c        | Make a copy of eigenvalues and corresponding error |
c        | bounds obtained from dneigh.                       |
c        %----------------------------------------------------%
c
         call dcopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call dcopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call dcopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
c
c        %---------------------------------------------------%
c        | Select the wanted Ritz values and their bounds    |
c        | to be used in the convergence test.               |
c        | The wanted part of the spectrum and corresponding |
c        | error bounds are in the last NEV loc. of RITZR,   |
c        | RITZI and BOUNDS respectively. The variables NEV  |
c        | and NP may be updated if the NEV-th wanted Ritz   |
c        | value has a non zero imaginary part. In this case |
c        | NEV is increased by one and NP decreased by one.  |
c        | NOTE: The last two arguments of dngets are no     |
c        | longer used as of version 2.1.                    |
c        %---------------------------------------------------%
c
         nev = nev0
         np = np0
         numcnv = nev
         call dngets (ishift, which, nev, np, ritzr, ritzi, 
     &                bounds, workl, workl(np+1))
         if (nev .eq. nev0+1) numcnv = nev0+1
c 
c        %-------------------%
c        | Convergence test. | 
c        %-------------------%
c
         call dcopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
         call dnconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1), 
     &        tol, nconv)
c 
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit, 
     &                  '_naup2: NEV, NP, NUMCNV, NCONV are')
            call dvout (logfil, kplusp, ritzr, ndigit,
     &           '_naup2: Real part of the eigenvalues of H')
            call dvout (logfil, kplusp, ritzi, ndigit,
     &           '_naup2: Imaginary part of the eigenvalues of H')
            call dvout (logfil, kplusp, bounds, ndigit, 
     &          '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
c
c        %---------------------------------------------------------%
c        | Count the number of unwanted Ritz values that have zero |
c        | Ritz estimates. If any Ritz estimates are equal to zero |
c        | then a leading block of H of order equal to at least    |
c        | the number of Ritz values with zero Ritz estimates has  |
c        | split off. None of these Ritz values may be removed by  |
c        | shifting. Decrease NP the number of shifts to apply. If |
c        | no shifts may be applied, then prepare to exit          |
c        %---------------------------------------------------------%
c
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
c     
         if ( (nconv .ge. numcnv) .or. 
     &        (iter .gt. mxiter) .or.
     &        (np .eq. 0) ) then
c
            if (msglvl .gt. 4) then
               call dvout(logfil, kplusp, workl(kplusp**2+1), ndigit,
     &             '_naup2: Real part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp+1),
     &                     ndigit,
     &             '_naup2: Imag part of the eig computed by _neigh:')
               call dvout(logfil, kplusp, workl(kplusp**2+kplusp*2+1),
     &                     ndigit,
     &             '_naup2: Ritz eistmates computed by _neigh:')
            end if
c     
c           %------------------------------------------------%
c           | Prepare to exit. Put the converged Ritz values |
c           | and corresponding bounds in RITZ(1:NCONV) and  |
c           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
c           | careful when NCONV > NP                        |
c           %------------------------------------------------%
c
c           %------------------------------------------%
c           |  Use h( 3,1 ) as storage to communicate  |
c           |  rnorm to _neupd if needed               |
c           %------------------------------------------%

            h(3,1) = rnorm
c
c           %----------------------------------------------%
c           | To be consistent with dngets, we first do a  |
c           | pre-processing sort in order to keep complex |
c           | conjugate pairs together.  This is similar   |
c           | to the pre-processing sort used in dngets    |
c           | except that the sort is done in the opposite |
c           | order.                                       |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SR'
            if (which .eq. 'SM') wprime = 'LR'
            if (which .eq. 'LR') wprime = 'SM'
            if (which .eq. 'SR') wprime = 'LM'
            if (which .eq. 'LI') wprime = 'SM'
            if (which .eq. 'SI') wprime = 'LM'
c
            call dsortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %----------------------------------------------%
c           | Now sort Ritz values so that converged Ritz  |
c           | values appear within the first NEV locations |
c           | of ritzr, ritzi and bounds, and the most     |
c           | desired one appears at the front.            |
c           %----------------------------------------------%
c
            if (which .eq. 'LM') wprime = 'SM'
            if (which .eq. 'SM') wprime = 'LM'
            if (which .eq. 'LR') wprime = 'SR'
            if (which .eq. 'SR') wprime = 'LR'
            if (which .eq. 'LI') wprime = 'SI'
            if (which .eq. 'SI') wprime = 'LI'
c
            call dsortc(wprime, .true., kplusp, ritzr, ritzi, bounds)
c
c           %--------------------------------------------------%
c           | Scale the Ritz estimate of each Ritz value       |
c           | by 1 / max(eps23,magnitude of the Ritz value).   |
c           %--------------------------------------------------%
c
            do 35 j = 1, numcnv
                temp = max(eps23,dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
c
c           %----------------------------------------------------%
c           | Sort the Ritz values according to the scaled Ritz  |
c           | esitmates.  This will push all the converged ones  |
c           | towards the front of ritzr, ritzi, bounds          |
c           | (in the case when NCONV < NEV.)                    |
c           %----------------------------------------------------%
c
            wprime = 'LR'
            call dsortc(wprime, .true., numcnv, bounds, ritzr, ritzi)
c
c           %----------------------------------------------%
c           | Scale the Ritz estimate back to its original |
c           | value.                                       |
c           %----------------------------------------------%
c
            do 40 j = 1, numcnv
                temp = max(eps23, dlapy2(ritzr(j),
     &                                   ritzi(j)))
                bounds(j) = bounds(j)*temp
 40         continue
c
c           %------------------------------------------------%
c           | Sort the converged Ritz values again so that   |
c           | the "threshold" value appears at the front of  |
c           | ritzr, ritzi and bound.                        |
c           %------------------------------------------------%
c
            call dsortc(which, .true., nconv, ritzr, ritzi, bounds)
c
            if (msglvl .gt. 1) then
               call dvout (logfil, kplusp, ritzr, ndigit,
     &            '_naup2: Sorted real part of the eigenvalues')
               call dvout (logfil, kplusp, ritzi, ndigit,
     &            '_naup2: Sorted imaginary part of the eigenvalues')
               call dvout (logfil, kplusp, bounds, ndigit,
     &            '_naup2: Sorted ritz estimates.')
            end if
c
c           %------------------------------------%
c           | Max iterations have been exceeded. | 
c           %------------------------------------%
c
            if (iter .gt. mxiter .and. nconv .lt. numcnv) info = 1
c
c           %---------------------%
c           | No shifts to apply. | 
c           %---------------------%
c
            if (np .eq. 0 .and. nconv .lt. numcnv) info = 2
c
            np = nconv
            go to 1100
c
         else if ( (nconv .lt. numcnv) .and. (ishift .eq. 1) ) then
c     
c           %-------------------------------------------------%
c           | Do not have all the requested eigenvalues yet.  |
c           | To prevent possible stagnation, adjust the size |
c           | of NEV.                                         |
c           %-------------------------------------------------%
c
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 3) then
               nev = 2
            end if
            np = kplusp - nev
c     
c           %---------------------------------------%
c           | If the size of NEV was just increased |
c           | resort the eigenvalues.               |
c           %---------------------------------------%
c     
            if (nevbef .lt. nev) 
     &         call dngets (ishift, which, nev, np, ritzr, ritzi, 
     &              bounds, workl, workl(np+1))
c
         end if              
c     
         if (msglvl .gt. 0) then
            call ivout (logfil, 1, nconv, ndigit, 
     &           '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, 
     &              '_naup2: NEV and NP are')
               call dvout (logfil, nev, ritzr(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- real part')
               call dvout (logfil, nev, ritzi(np+1), ndigit,
     &              '_naup2: "wanted" Ritz values -- imag part')
               call dvout (logfil, nev, bounds(np+1), ndigit,
     &              '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
c
         if (ishift .eq. 0) then
c
c           %-------------------------------------------------------%
c           | User specified shifts: reverse comminucation to       |
c           | compute the shifts. They are returned in the first    |
c           | 2*NP locations of WORKL.                              |
c           %-------------------------------------------------------%
c
            ushift = .true.
            ido = 3
            go to 9000
         end if
c 
   50    continue
c
c        %------------------------------------%
c        | Back from reverse communication;   |
c        | User specified shifts are returned |
c        | in WORKL(1:2*NP)                   |
c        %------------------------------------%
c
         ushift = .false.
c
         if ( ishift .eq. 0 ) then
c 
c            %----------------------------------%
c            | Move the NP shifts from WORKL to |
c            | RITZR, RITZI to free up WORKL    |
c            | for non-exact shift case.        |
c            %----------------------------------%
c
             call dcopy (np, workl,       1, ritzr, 1)
             call dcopy (np, workl(np+1), 1, ritzi, 1)
         end if
c
         if (msglvl .gt. 2) then 
            call ivout (logfil, 1, np, ndigit, 
     &                  '_naup2: The number of shifts to apply ')
            call dvout (logfil, np, ritzr, ndigit,
     &                  '_naup2: Real part of the shifts')
            call dvout (logfil, np, ritzi, ndigit,
     &                  '_naup2: Imaginary part of the shifts')
            if ( ishift .eq. 1 ) 
     &          call dvout (logfil, np, bounds, ndigit,
     &                  '_naup2: Ritz estimates of the shifts')
         end if
c
c        %---------------------------------------------------------%
c        | Apply the NP implicit shifts by QR bulge chasing.       |
c        | Each shift is applied to the whole upper Hessenberg     |
c        | matrix H.                                               |
c        | The first 2*N locations of WORKD are used as workspace. |
c        %---------------------------------------------------------%
c
         call dnapps (n, nev, np, ritzr, ritzi, v, ldv, 
     &                h, ldh, resid, q, ldq, workl, workd)
c
c        %---------------------------------------------%
c        | Compute the B-norm of the updated residual. |
c        | Keep B*RESID in WORKD(1:N) to be used in    |
c        | the first step of the next call to dnaitr.  |
c        %---------------------------------------------%
c
         cnorm = .true.
         call second (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
c 
c           %----------------------------------%
c           | Exit in order to compute B*RESID |
c           %----------------------------------%
c 
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd, 1)
         end if
c 
  100    continue
c 
c        %----------------------------------%
c        | Back from reverse communication; |
c        | WORKD(1:N) := B*RESID            |
c        %----------------------------------%
c
         if (bmat .eq. 'G') then
            call second (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c 
         if (bmat .eq. 'G') then         
            rnorm = ddot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
         cnorm = .false.
c
         if (msglvl .gt. 2) then
            call dvout (logfil, 1, rnorm, ndigit, 
     &      '_naup2: B-norm of residual for compressed factorization')
            call dmout (logfil, nev, nev, h, ldh, ndigit,
     &        '_naup2: Compressed upper Hessenberg matrix H')
         end if
c 
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 1100 continue
c
      mxiter = iter
      nev = numcnv
c     
 1200 continue
      ido = 99
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      call second (t1)
      tnaup2 = t1 - t0
c     
 9000 continue
c
c     %---------------%
c     | End of dnaup2 |
c     %---------------%
c
      return
      end
c\BeginDoc
c
c\Name: dnaupd
c
c\Description: 
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This subroutine computes approximations to a few eigenpairs 
c  of a linear operator "OP" with respect to a semi-inner product defined by 
c  a symmetric positive semi-definite real matrix B. B may be the identity 
c  matrix. NOTE: If the linear operator "OP" is real and symmetric 
c  with respect to the real positive semi-definite symmetric matrix B, 
c  i.e. B*OP = (OP`)*B, then subroutine dsaupd should be used instead.
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dnaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
c           Note: If sigma is real, i.e. imaginary part of sigma is zero;
c                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M 
c                 amu == 1/(lambda-sigma). 
c  
c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
c           ===> shift-and-invert mode (in real arithmetic)
c           If OP*x = amu*x, then 
c           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
c
c  Both mode 3 and 4 give the same enhancement to eigenvalues close to
c  the (complex) shift sigma.  However, as lambda goes to infinity,
c  the operator OP in mode 4 dampens the eigenvalues more strongly than
c  does OP defined in mode 3.
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dnaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to dnaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dnaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3 and 4, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) real and imaginary parts 
c                    of the shifts where INPTR(14) is the pointer
c                    into WORKL for placing the shifts. See Remark
c                    5 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT/OUTPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 2 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Arnoldi vectors are generated, the algorithm generates 
c          approximately NCV-NEV Arnoldi vectors at each subsequent update 
c          iteration. Most of the cost in generating each Arnoldi vector is 
c          in the matrix-vector operation OP*x. 
c          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz 
c          values are kept together. (See remark 4 below)
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The real and imaginary
c                      parts of the NCV eigenvalues of the Hessenberg
c                      matrix H are returned in the part of the WORKL 
c                      array corresponding to RITZR and RITZI. See remark 
c                      5 below.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to 
c                      restarting the iteration with a starting vector
c                      that is a linear combination of approximate Schur
c                      vectors associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.  
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4; See under \Description of dnaupd for the 
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dnaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          5 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
c                    H in WORKL.
c          IPNTR(6): pointer to the real part of the ritz value array 
c                    RITZR in WORKL.
c          IPNTR(7): pointer to the imaginary part of the ritz value array
c                    RITZI in WORKL.
c          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZR and RITZI in WORKL.
c
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
c
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the 
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of 
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
c          associated with the converged Ritz values is desired, see remark
c          2 below, subroutine dneupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 6*NCV.
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when
c     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
c     original problem may be obtained with the ARPACK subroutine dneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz 
c     values is needed, the user must call dneupd immediately following 
c     completion of dnaupd. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular 
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
c     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically. 
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
c     NP = IPARAM(8) real and imaginary parts of the shifts in locations 
c         real part                  imaginary part
c         -----------------------    --------------
c     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
c     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
c                        .                          .
c                        .                          .
c                        .                          .
c     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
c
c     Only complex conjugate pairs of shifts may be applied and the pairs 
c     must be placed in consecutive locations. The real part of the 
c     eigenvalues of the current upper Hessenberg matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part 
c     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
c     according to the order defined by WHICH. The complex conjugate
c     pairs are kept together and the associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note: 
c
c  Fortran-D syntax:
c  ================
c  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c  
c  CM2/CM5 syntax:
c  ==============
c  
c-----------------------------------------------------------------------
c
c     include   'ex-nonsym.doc'
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     dnaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/16/93: Version '1.1'
c
c\SCCS Information: @(#) 
c FILE: naupd.F   SID: 2.10   DATE OF SID: 08/23/02   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnaupd
     &   ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, 
     &     ipntr, workd, workl, lworkl, info )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      Double precision
     &           resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, 
     &           ldh, ldq, levec, mode, msglvl, mxiter, nb,
     &           nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq,
     &           levec, mode, msglvl, mxiter, nb, nev0, next,
     &           np, ritzi, ritzr
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dnaup2, dvout, ivout, second, dstatn
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlamch
      external   dlamch
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
      if (ido .eq. 0) then
c 
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call dstatn
         call second (t0)
         msglvl = mnaupd
c
c        %----------------%
c        | Error checking |
c        %----------------%
c
         ierr   = 0
         ishift = iparam(1)
c         levec  = iparam(2)
         mxiter = iparam(3)
c         nb     = iparam(4)
         nb     = 1
c
c        %--------------------------------------------%
c        | Revision 2 performs only implicit restart. |
c        %--------------------------------------------%
c
         iupd   = 1
         mode   = iparam(7)
c
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev+1 .or.  ncv .gt. n) then
            ierr = -3
         else if (mxiter .le.          0) then
            ierr = 4
         else if (which .ne. 'LM' .and.
     &       which .ne. 'SM' .and.
     &       which .ne. 'LR' .and.
     &       which .ne. 'SR' .and.
     &       which .ne. 'LI' .and.
     &       which .ne. 'SI') then
            ierr = -5
         else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode .lt. 1 .or. mode .gt. 4) then
            ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
            ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
            ierr = -12
         end if
c 
c        %------------%
c        | Error Exit |
c        %------------%
c
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
c 
c        %------------------------%
c        | Set default parameters |
c        %------------------------%
c
         if (nb .le. 0)				nb = 1
         if (tol .le. zero)			tol = dlamch('EpsMach')
c
c        %----------------------------------------------%
c        | NP is the number of additional steps to      |
c        | extend the length NEV Lanczos factorization. |
c        | NEV0 is the local variable designating the   |
c        | size of the invariant subspace desired.      |
c        %----------------------------------------------%
c
         np     = ncv - nev
         nev0   = nev 
c 
c        %-----------------------------%
c        | Zero out internal workspace |
c        %-----------------------------%
c
         do 10 j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
  10     continue
c 
c        %-------------------------------------------------------------%
c        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
c        | etc... and the remaining workspace.                         |
c        | Also update pointer to be used on output.                   |
c        | Memory is laid out as follows:                              |
c        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
c        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
c        |                                   parts of ritz values      |
c        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
c        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
c        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
c        | The final workspace is needed by subroutine dneigh called   |
c        | by dnaup2. Subroutine dneigh calls LAPACK routines for      |
c        | calculating eigenvalues and the last row of the eigenvector |
c        | matrix.                                                     |
c        %-------------------------------------------------------------%
c
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
c
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw 
c
      end if
c
c     %-------------------------------------------------------%
c     | Carry out the Implicitly restarted Arnoldi Iteration. |
c     %-------------------------------------------------------%
c
      call dnaup2 
     &   ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     &     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), 
     &     workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), 
     &     ipntr, workd, info )
c 
c     %--------------------------------------------------%
c     | ido .ne. 99 implies use of reverse communication |
c     | to compute operations involving OP or shifts.    |
c     %--------------------------------------------------%
c
      if (ido .eq. 3) iparam(8) = np
      if (ido .ne. 99) go to 9000
c 
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
c
c     %------------------------------------%
c     | Exit if there was an informational |
c     | error within dnaup2.               |
c     %------------------------------------%
c
      if (info .lt. 0) go to 9000
      if (info .eq. 2) info = 3
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, mxiter, ndigit,
     &               '_naupd: Number of update iterations taken')
         call ivout (logfil, 1, np, ndigit,
     &               '_naupd: Number of wanted "converged" Ritz values')
         call dvout (logfil, np, workl(ritzr), ndigit, 
     &               '_naupd: Real part of the final Ritz values')
         call dvout (logfil, np, workl(ritzi), ndigit, 
     &               '_naupd: Imaginary part of the final Ritz values')
         call dvout (logfil, np, workl(bounds), ndigit, 
     &               '_naupd: Associated Ritz estimates')
      end if
c
      call second (t1)
      tnaupd = t1 - t0
c
      if (msglvl .gt. 0) then
c
c        %--------------------------------------------------------%
c        | Version Number & Version Date are defined in version.h |
c        %--------------------------------------------------------%
c
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     &                  tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref,
     &                  tgetv0, tneigh, tngets, tnapps, tnconv, trvec
 1000    format (//,
     &      5x, '=============================================',/
     &      5x, '= Nonsymmetric implicit Arnoldi update code =',/
     &      5x, '= Version Number: ', ' 2.4', 21x, ' =',/
     &      5x, '= Version Date:   ', ' 07/31/96', 16x,   ' =',/
     &      5x, '=============================================',/
     &      5x, '= Summary of timing statistics              =',/
     &      5x, '=============================================',//)
 1100    format (
     &      5x, 'Total number update iterations             = ', i5,/
     &      5x, 'Total number of OP*x operations            = ', i5,/
     &      5x, 'Total number of B*x operations             = ', i5,/
     &      5x, 'Total number of reorthogonalization steps  = ', i5,/
     &      5x, 'Total number of iterative refinement steps = ', i5,/
     &      5x, 'Total number of restart steps              = ', i5,/
     &      5x, 'Total time in user OP*x operation          = ', f12.6,/
     &      5x, 'Total time in user B*x operation           = ', f12.6,/
     &      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     &      5x, 'Total time in naup2 routine                = ', f12.6,/
     &      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     &      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     &      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     &      5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/
     &      5x, 'Total time in getting the shifts           = ', f12.6,/
     &      5x, 'Total time in applying the shifts          = ', f12.6,/
     &      5x, 'Total time in convergence testing          = ', f12.6,/
     &      5x, 'Total time in computing final Ritz vectors = ', f12.6/)
      end if
c
 9000 continue
c
      return
c
c     %---------------%
c     | End of dnaupd |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dnconv
c
c\Description: 
c  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dnconv
c     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZR,  Double precision arrays of length N.  (INPUT)
c  RITZI   Real and imaginary parts of the Ritz values to be checked
c          for convergence.

c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates for the Ritz values in RITZR and RITZI.
c
c  TOL     Double precision scalar.  (INPUT)
c          Desired backward error for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     second  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dnconv (n, ritzr, ritzi, bounds, tol, nconv)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    n, nconv
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%

      Double precision
     &           ritzr(n), ritzi(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      Double precision
     &           temp, eps23
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dlamch
      external   dlapy2, dlamch

c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %-------------------------------------------------------------%
c     | Convergence test: unlike in the symmetric code, I am not    |
c     | using things like refined error bounds and gap condition    |
c     | because I don't know the exact equivalent concept.          |
c     |                                                             |
c     | Instead the i-th Ritz value is considered "converged" when: |
c     |                                                             |
c     |     bounds(i) .le. ( TOL * | ritz | )                       |
c     |                                                             |
c     | for some appropriate choice of norm.                        |
c     %-------------------------------------------------------------%
c
      call second (t0)
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
      nconv  = 0
      do 20 i = 1, n
         temp = max( eps23, dlapy2( ritzr(i), ritzi(i) ) )
         if (bounds(i) .le. tol*temp)   nconv = nconv + 1
   20 continue
c 
      call second (t1)
      tnconv = tnconv + (t1 - t0)
c 
      return
c
c     %---------------%
c     | End of dnconv |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dneigh
c
c\Description:
c  Compute the eigenvalues of the current upper Hessenberg matrix
c  and the corresponding Ritz estimates given the current residual norm.
c
c\Usage:
c  call dneigh
c     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          Residual norm corresponding to the current upper Hessenberg 
c          matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the matrix H.
c
c  H       Double precision N by N array.  (INPUT)
c          H contains the current upper Hessenberg matrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RITZR,  Double precision arrays of length N.  (OUTPUT)
c  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real 
c          (respectively imaginary) parts of the eigenvalues of H.
c
c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the Ritz estimates associated with
c          the eigenvalues RITZR and RITZI.  This is equal to RNORM 
c          times the last components of the eigenvectors corresponding 
c          to the eigenvalues in RITZR and RITZI.
c
c  Q       Double precision N by N array.  (WORKSPACE)
c          Workspace needed to store the eigenvectors of H.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  This is needed to keep the full Schur form
c          of H and also in the calculation of the eigenvectors of H.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dlaqrb or dtrevc.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dlaqrb  ARPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix and last row of the Schur vectors.
c     second  ARPACK utility routine for timing.
c     dmout   ARPACK utility routine that prints matrices
c     dvout   ARPACK utility routine that prints vectors.
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dneigh (rnorm, n, h, ldh, ritzr, ritzi, bounds, 
     &                   q, ldq, workl, ierr)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, n, ldh, ldq
      Double precision     
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision     
     &           bounds(n), h(ldh,n), q(ldq,n), ritzi(n), ritzr(n),
     &           workl(n*(n+3))
c 
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision     
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c 
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      logical    select(1)
      integer    i, iconj, msglvl
      Double precision     
     &           temp, vl(1)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dlacpy, dlaqrb, dtrevc, dvout, second
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           dlapy2, dnrm2
      external   dlapy2, dnrm2
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic  abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call second (t0)
      msglvl = mneigh
c 
      if (msglvl .gt. 2) then
          call dmout (logfil, n, n, h, ldh, ndigit, 
     &         '_neigh: Entering upper Hessenberg matrix H ')
      end if
c 
c     %-----------------------------------------------------------%
c     | 1. Compute the eigenvalues, the last components of the    |
c     |    corresponding Schur vectors and the full Schur form T  |
c     |    of the current upper Hessenberg matrix H.              |
c     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  |
c     | and the last components of the Schur vectors in BOUNDS.   |
c     %-----------------------------------------------------------%
c
      call dlacpy ('All', n, n, h, ldh, workl, n)
      call dlaqrb (.true., n, 1, n, workl, n, ritzr, ritzi, bounds,
     &             ierr)
      if (ierr .ne. 0) go to 9000
c
      if (msglvl .gt. 1) then
         call dvout (logfil, n, bounds, ndigit,
     &              '_neigh: last row of the Schur matrix for H')
      end if
c
c     %-----------------------------------------------------------%
c     | 2. Compute the eigenvectors of the full Schur form T and  |
c     |    apply the last components of the Schur vectors to get  |
c     |    the last components of the corresponding eigenvectors. |
c     | Remember that if the i-th and (i+1)-st eigenvalues are    |
c     | complex conjugate pairs, then the real & imaginary part   |
c     | of the eigenvector components are split across adjacent   |
c     | columns of Q.                                             |
c     %-----------------------------------------------------------%
c
      call dtrevc ('R', 'A', select, n, workl, n, vl, n, q, ldq,
     &             n, n, workl(n*n+1), ierr)
c
      if (ierr .ne. 0) go to 9000
c
c     %------------------------------------------------%
c     | Scale the returning eigenvectors so that their |
c     | euclidean norms are all one. LAPACK subroutine |
c     | dtrevc returns each eigenvector normalized so  |
c     | that the element of largest magnitude has      |
c     | magnitude 1; here the magnitude of a complex   |
c     | number (x,y) is taken to be |x| + |y|.         |
c     %------------------------------------------------%
c
      iconj = 0
      do 10 i=1, n
         if ( abs( ritzi(i) ) .le. zero ) then
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c    
            temp = dnrm2( n, q(1,i), 1 )
            call dscal ( n, one / temp, q(1,i), 1 )
         else
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we further normalize by the      |
c           | square root of two.                       |
c           %-------------------------------------------%
c
            if (iconj .eq. 0) then
               temp = dlapy2( dnrm2( n, q(1,i), 1 ), 
     &                        dnrm2( n, q(1,i+1), 1 ) )
               call dscal ( n, one / temp, q(1,i), 1 )
               call dscal ( n, one / temp, q(1,i+1), 1 )
               iconj = 1
            else
               iconj = 0
            end if
         end if         
   10 continue
c
      call dgemv ('T', n, n, one, q, ldq, bounds, 1, zero, workl, 1)
c
      if (msglvl .gt. 1) then
         call dvout (logfil, n, workl, ndigit,
     &              '_neigh: Last row of the eigenvector matrix for H')
      end if
c
c     %----------------------------%
c     | Compute the Ritz estimates |
c     %----------------------------%
c
      iconj = 0
      do 20 i = 1, n
         if ( abs( ritzi(i) ) .le. zero ) then
c
c           %----------------------%
c           | Real eigenvalue case |
c           %----------------------%
c    
            bounds(i) = rnorm * abs( workl(i) )
         else
c
c           %-------------------------------------------%
c           | Complex conjugate pair case. Note that    |
c           | since the real and imaginary part of      |
c           | the eigenvector are stored in consecutive |
c           | columns, we need to take the magnitude    |
c           | of the last components of the two vectors |
c           %-------------------------------------------%
c
            if (iconj .eq. 0) then
               bounds(i) = rnorm * dlapy2( workl(i), workl(i+1) )
               bounds(i+1) = bounds(i)
               iconj = 1
            else
               iconj = 0
            end if
         end if
   20 continue
c
      if (msglvl .gt. 2) then
         call dvout (logfil, n, ritzr, ndigit,
     &              '_neigh: Real part of the eigenvalues of H')
         call dvout (logfil, n, ritzi, ndigit,
     &              '_neigh: Imaginary part of the eigenvalues of H')
         call dvout (logfil, n, bounds, ndigit,
     &              '_neigh: Ritz estimates for the eigenvalues of H')
      end if
c
      call second (t1)
      tneigh = tneigh + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of dneigh |
c     %---------------%
c
      end
c\BeginDoc
c
c\Name: dneupd 
c
c\Description: 
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to DNAUPD .  DNAUPD  must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.  The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  See documentation in the header of the subroutine DNAUPD  for 
c  definition of OP as well as other terms and the relation of computed
c  Ritz values and Ritz vectors of OP with respect to the given problem
c  A*z = lambda*B*z.  For a brief description, see definitions of 
c  IPARAM(7), MODE and WHICH in the documentation of DNAUPD .
c
c\Usage:
c  call dneupd  
c     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, 
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, 
c       LWORKL, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether a basis for the invariant subspace corresponding 
c          to the converged Ritz value approximations for the eigenproblem 
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
c                                See Remarks below. 
c 
c  HOWMNY  Character*1  (INPUT) 
c          Specifies the form of the basis for the invariant subspace 
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors; 
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
c
c  DR      Double precision  array of dimension NEV+1.  (OUTPUT)
c          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains 
c          the real part of the Ritz  approximations to the eigenvalues of 
c          A*z = lambda*B*z. 
c          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
c          DR contains the real part of the Ritz values of OP computed by 
c          DNAUPD . A further computation must be performed by the user
c          to transform the Ritz values computed for OP by DNAUPD  to those
c          of the original system A*z = lambda*B*z. See remark 3 below.
c
c  DI      Double precision  array of dimension NEV+1.  (OUTPUT)
c          On exit, DI contains the imaginary part of the Ritz value 
c          approximations to the eigenvalues of A*z = lambda*B*z associated
c          with DR.
c
c          NOTE: When Ritz values are complex, they will come in complex 
c                conjugate pairs.  If eigenvectors are requested, the 
c                corresponding Ritz vectors will also come in conjugate 
c                pairs and the real and imaginary parts of these are 
c                represented in two consecutive columns of the array Z 
c                (see below).
c
c  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
c          Z represent approximate eigenvectors (Ritz vectors) corresponding 
c          to the NCONV=IPARAM(5) Ritz values for eigensystem 
c          A*z = lambda*B*z. 
c 
c          The complex Ritz vector associated with the Ritz value 
c          with positive imaginary part is stored in two consecutive 
c          columns.  The first column holds the real part of the Ritz 
c          vector and the second column holds the imaginary part.  The 
c          Ritz vector associated with the Ritz value with negative 
c          imaginary part is simply the complex conjugate of the Ritz vector 
c          associated with the positive imaginary part.
c
c          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by DNAUPD .  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
c
c  SIGMAR  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the real part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  SIGMAI  Double precision   (INPUT)
c          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift. 
c          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
c
c  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPD  that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DNEUPD  following the last call
c         to DNAUPD .  These arguments MUST NOT BE MODIFIED between
c         the the last call to DNAUPD  and the call to DNEUPD .
c
c  Three of these parameters (V, WORKL, INFO) are also output parameters:
c
c  V       Double precision  N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by DNAUPD  .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.  See Remark 2 below.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
c          dnaupd .  They are not changed by dneupd .
c          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
c          real and imaginary part of the untransformed Ritz values,
c          the upper quasi-triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by dneupd .
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
c                     original system.
c          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
c                     the original system.
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     dneupd  if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine dlahqr 
c                could not be reordered by LAPACK routine dtrsen .
c                Re-enter subroutine dneupd  with IPARAM(5)=NCV and 
c                increase the size of the arrays DR and DI to have 
c                dimension at least dimension NCV and allocate at least NCV 
c                columns for Z. NOTE: Not necessary if Z and V share 
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from calculation of a real Schur form.
c                Informational error from LAPACK routine dlahqr .
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine dtrevc .
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: DNAUPD  did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: DNEUPD got a different count of the number of converged
c                 Ritz values than DNAUPD got.  This indicates the user
c                 probably made an error in passing data from DNAUPD to
c                 DNEUPD or that the data was modified before entering
c                 DNEUPD
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
c     Real Matrices", Linear Algebra and its Applications, vol 88/89,
c     pp 575-595, (1987).
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     dmout    ARPACK utility routine that prints matrices
c     dvout    ARPACK utility routine that prints vectors.
c     dgeqr2   LAPACK routine that computes the QR factorization of 
c             a matrix.
c     dlacpy   LAPACK matrix copy routine.
c     dlahqr   LAPACK routine to compute the real Schur form of an
c             upper Hessenberg matrix.
c     dlamch   LAPACK routine that determines machine constants.
c     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlaset   LAPACK matrix initialization routine.
c     dorm2r   LAPACK routine that applies an orthogonal matrix in 
c             factored form.
c     dtrevc   LAPACK routine to compute the eigenvectors of a matrix
c             in upper quasi-triangular form.
c     dtrsen   LAPACK routine that re-orders the Schur form.
c     dtrmm    Level 3 BLAS matrix times an upper triangular matrix.
c     dger     Level 2 BLAS rank one update to a matrix.
c     dcopy    Level 1 BLAS that copies one vector to another .
c     ddot     Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2    Level 1 BLAS that computes the norm of a vector.
c     dscal    Level 1 BLAS that scales a vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c     Let trans(X) denote the transpose of X.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .TRUE. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately 
c     satisfied. Here T is the leading submatrix of order IPARAM(5) of the 
c     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
c     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
c     each 2-by-2 diagonal block has its diagonal elements equal and its
c     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
c     diagonal block is a complex conjugate pair of Ritz values. The real
c     Ritz values are stored on the diagonal of T.
c
c  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
c     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
c     values computed by DNAUPD  for OP to those of A*z = lambda*B*z. 
c     Set RVEC = .true. and HOWMNY = 'A', and
c     compute 
c           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
c     If DI(I) is not equal to zero and DI(I+1) = - D(I), 
c     then the desired real and imaginary parts of the Ritz value are
c           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
c           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I), 
c     respectively.
c     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
c     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
c     quasi-triangular matrix of order IPARAM(5) is computed. See remark
c     2 above.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics          
c     Rice University           
c     Houston, Texas            
c 
c\SCCS Information: @(#) 
c FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2 
c
c\EndLib
c
c-----------------------------------------------------------------------
      subroutine dneupd (rvec , howmny, select, dr    , di,    
     &                   z    , ldz   , sigmar, sigmai, workev,
     &                   bmat , n     , which , nev   , tol,
     &                   resid, ncv   , v     , ldv   , iparam,
     &                   ipntr, workd , workl , lworkl, info)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision      
     &           sigmar, sigmai, tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Double precision 
     &           dr(nev+1)    , di(nev+1), resid(n)  , 
     &           v(ldv,ncv)   , z(ldz,*) , workd(3*n), 
     &           workl(lworkl), workev(3*ncv)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision 
     &           one, zero
      parameter (one = 1.0D+0 , zero = 0.0D+0 )
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character  type*6
      integer    bounds, ierr  , ih    , ihbds   , 
     &           iheigr, iheigi, iconj , nconv   , 
     &           invsub, iuptri, iwev  , iwork(1),
     &           j     , k     , ldh   , ldq     ,
     &           mode  , msglvl, outncv, ritzr   ,
     &           ritzi , wri   , wrr   , irr     ,
     &           iri   , ibd   , ishift, numcnv  ,
     &           np    , jj 
      logical    reord
      Double precision 
     &           conds  , rnorm, sep  , temp,
     &           vl(1,1), temp1, eps23
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy  , dger   , dgeqr2 , dlacpy , 
     &           dlahqr , dlaset , dmout  , dorm2r , 
     &           dtrevc , dtrmm  , dtrsen , dscal  , 
     &           dvout  , ivout
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision 
     &           dlapy2 , dnrm2 , dlamch , ddot 
      external   dlapy2 , dnrm2 , dlamch , ddot 
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs, min, sqrt
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 
c     %------------------------%
c     | Set default parameters |
c     %------------------------%
c
      msglvl = mneupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
c
c     %---------------------------------%
c     | Get machine dependent constant. |
c     %---------------------------------%
c
      eps23 = dlamch ('Epsilon-Machine')
      eps23 = eps23**(2.0D+0  / 3.0D+0 )
c
c     %--------------%
c     | Quick return |
c     %--------------%
c
      ierr = 0
c
      if (nconv .le. 0) then
         ierr = -14
      else if (n .le. 0) then
         ierr = -1
      else if (nev .le. 0) then
         ierr = -2
      else if (ncv .le. nev+1 .or.  ncv .gt. n) then
         ierr = -3
      else if (which .ne. 'LM' .and.
     &        which .ne. 'SM' .and.
     &        which .ne. 'LR' .and.
     &        which .ne. 'SR' .and.
     &        which .ne. 'LI' .and.
     &        which .ne. 'SI') then
         ierr = -5
      else if (bmat .ne. 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl .lt. 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny .ne. 'A' .and.
     &           howmny .ne. 'P' .and.
     &           howmny .ne. 'S') .and. rvec ) then
         ierr = -13
      else if (howmny .eq. 'S' ) then
         ierr = -12
      end if
c     
      if (mode .eq. 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode .eq. 3 .and. sigmai .eq. zero) then
         type = 'SHIFTI'
      else if (mode .eq. 3 ) then
         type = 'REALPT'
      else if (mode .eq. 4 ) then
         type = 'IMAGPT'
      else 
                                              ierr = -10
      end if
      if (mode .eq. 1 .and. bmat .eq. 'G')    ierr = -11
c
c     %------------%
c     | Error Exit |
c     %------------%
c
      if (ierr .ne. 0) then
         info = ierr
         go to 9000
      end if
c 
c     %--------------------------------------------------------%
c     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
c     | etc... and the remaining workspace.                    |
c     | Also update pointer to be used on output.              |
c     | Memory is laid out as follows:                         |
c     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
c     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
c     |                                   parts of ritz values |
c     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
c     %--------------------------------------------------------%
c
c     %-----------------------------------------------------------%
c     | The following is used and set by DNEUPD .                  |
c     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
c     |                             real part of the Ritz values. |
c     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
c     |                        imaginary part of the Ritz values. |
c     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
c     |                           error bounds of the Ritz values |
c     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
c     |                             quasi-triangular matrix for H |
c     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
c     |       associated matrix representation of the invariant   |
c     |       subspace for H.                                     |
c     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
c     %-----------------------------------------------------------%
c     
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
c
c     %-----------------------------------------%
c     | irr points to the REAL part of the Ritz |
c     |     values computed by _neigh before    |
c     |     exiting _naup2.                     |
c     | iri points to the IMAGINARY part of the |
c     |     Ritz values computed by _neigh      |
c     |     before exiting _naup2.              |
c     | ibd points to the Ritz estimates        |
c     |     computed by _neigh before exiting   |
c     |     _naup2.                             |
c     %-----------------------------------------%
c
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
c
c     %------------------------------------%
c     | RNORM is B-norm of the RESID(1:N). |
c     %------------------------------------%
c
      rnorm = workl(ih+2)
      workl(ih+2) = zero
c
      if (msglvl .gt. 2) then
         call dvout (logfil, ncv, workl(irr), ndigit,
     &   '_neupd: Real part of Ritz values passed in from _NAUPD.')
         call dvout (logfil, ncv, workl(iri), ndigit,
     &   '_neupd: Imag part of Ritz values passed in from _NAUPD.')
         call dvout (logfil, ncv, workl(ibd), ndigit,
     &   '_neupd: Ritz estimates passed in from _NAUPD.')
      end if
c
      if (rvec) then
c     
         reord = .false.
c
c        %---------------------------------------------------%
c        | Use the temporary bounds array to store indices   |
c        | These will be used to mark the select array later |
c        %---------------------------------------------------%
c
         do 10 j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
   10    continue
c
c        %-------------------------------------%
c        | Select the wanted Ritz values.      |
c        | Sort the Ritz values so that the    |
c        | wanted ones appear at the tailing   |
c        | NEV positions of workl(irr) and     |
c        | workl(iri).  Move the corresponding |
c        | error estimates in workl(bound)     |
c        | accordingly.                        |
c        %-------------------------------------%
c
         np     = ncv - nev
         ishift = 0
         call dngets (ishift       , which     , nev       , 
     &                np           , workl(irr), workl(iri),
     &                workl(bounds), workl     , workl(np+1))
c
         if (msglvl .gt. 2) then
            call dvout (logfil, ncv, workl(irr), ndigit,
     &      '_neupd: Real part of Ritz values after calling _NGETS.')
            call dvout (logfil, ncv, workl(iri), ndigit,
     &      '_neupd: Imag part of Ritz values after calling _NGETS.')
            call dvout (logfil, ncv, workl(bounds), ndigit,
     &      '_neupd: Ritz value indices after calling _NGETS.')
         end if
c
c        %-----------------------------------------------------%
c        | Record indices of the converged wanted Ritz values  |
c        | Mark the select array for possible reordering       |
c        %-----------------------------------------------------%
c
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23,
     &                 dlapy2 ( workl(irr+ncv-j), workl(iri+ncv-j) ))
            jj = workl(bounds + ncv - j)
            if (numcnv .lt. nconv .and.
     &          workl(ibd+jj-1) .le. tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj .gt. nev) reord = .true.
            endif
   11    continue
c
c        %-----------------------------------------------------------%
c        | Check the count (numcnv) of converged Ritz values with    |
c        | the number (nconv) reported by dnaupd.  If these two      |
c        | are different then there has probably been an error       |
c        | caused by incorrect passing of the dnaupd data.           |
c        %-----------------------------------------------------------%
c
         if (msglvl .gt. 2) then
             call ivout(logfil, 1, numcnv, ndigit,
     &            '_neupd: Number of specified eigenvalues')
             call ivout(logfil, 1, nconv, ndigit,
     &            '_neupd: Number of "converged" eigenvalues')
         end if
c
         if (numcnv .ne. nconv) then
            info = -15
            go to 9000
         end if
c
c        %-----------------------------------------------------------%
c        | Call LAPACK routine dlahqr  to compute the real Schur form |
c        | of the upper Hessenberg matrix returned by DNAUPD .        |
c        | Make a copy of the upper Hessenberg matrix.               |
c        | Initialize the Schur vector matrix Q to the identity.     |
c        %-----------------------------------------------------------%
c     
         call dcopy (ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call dlaset ('All', ncv, ncv, 
     &                zero , one, workl(invsub),
     &                ldq)
         call dlahqr (.true., .true.       , ncv, 
     &                1     , ncv          , workl(iuptri), 
     &                ldh   , workl(iheigr), workl(iheigi),
     &                1     , ncv          , workl(invsub), 
     &                ldq   , ierr)
         call dcopy (ncv         , workl(invsub+ncv-1), ldq, 
     &               workl(ihbds), 1)
c     
         if (ierr .ne. 0) then
            info = -8
            go to 9000
         end if
c     
         if (msglvl .gt. 1) then
            call dvout (logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H')
            call dvout (logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imaginary part of the Eigenvalues of H')
            call dvout (logfil, ncv, workl(ihbds), ndigit,
     &           '_neupd: Last row of the Schur vector matrix')
            if (msglvl .gt. 3) then
               call dmout (logfil       , ncv, ncv   , 
     &                     workl(iuptri), ldh, ndigit,
     &              '_neupd: The upper quasi-triangular matrix ')
            end if
         end if 
c
         if (reord) then
c     
c           %-----------------------------------------------------%
c           | Reorder the computed upper quasi-triangular matrix. | 
c           %-----------------------------------------------------%
c     
            call dtrsen ('None'       , 'V'          , 
     &                   select       , ncv          ,
     &                   workl(iuptri), ldh          , 
     &                   workl(invsub), ldq          , 
     &                   workl(iheigr), workl(iheigi), 
     &                   nconv        , conds        ,
     &                   sep          , workl(ihbds) , 
     &                   ncv          , iwork        ,
     &                   1            , ierr)
c
            if (ierr .eq. 1) then
               info = 1
               go to 9000
            end if
c
            if (msglvl .gt. 2) then
                call dvout (logfil, ncv, workl(iheigr), ndigit,
     &           '_neupd: Real part of the eigenvalues of H--reordered')
                call dvout (logfil, ncv, workl(iheigi), ndigit,
     &           '_neupd: Imag part of the eigenvalues of H--reordered')
                if (msglvl .gt. 3) then
                   call dmout (logfil       , ncv, ncv   , 
     &                         workl(iuptri), ldq, ndigit,
     &             '_neupd: Quasi-triangular matrix after re-ordering')
                end if
            end if
c     
         end if
c
c        %---------------------------------------%
c        | Copy the last row of the Schur vector |
c        | into workl(ihbds).  This will be used |
c        | to compute the Ritz estimates of      |
c        | converged Ritz values.                |
c        %---------------------------------------%
c
         call dcopy (ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
c
c        %----------------------------------------------------%
c        | Place the computed eigenvalues of H into DR and DI |
c        | if a spectral transformation was not used.         |
c        %----------------------------------------------------%
c
         if (type .eq. 'REGULR') then 
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
         end if
c     
c        %----------------------------------------------------------%
c        | Compute the QR factorization of the matrix representing  |
c        | the wanted invariant subspace located in the first NCONV |
c        | columns of workl(invsub,ldq).                            |
c        %----------------------------------------------------------%
c     
         call dgeqr2 (ncv, nconv , workl(invsub), 
     &               ldq, workev, workev(ncv+1),
     &               ierr)
c
c        %---------------------------------------------------------%
c        | * Postmultiply V by Q using dorm2r .                     |   
c        | * Copy the first NCONV columns of VQ into Z.            |
c        | * Postmultiply Z by R.                                  |
c        | The N by NCONV matrix Z is now a matrix representation  |
c        | of the approximate invariant subspace associated with   |
c        | the Ritz values in workl(iheigr) and workl(iheigi)      |
c        | The first NCONV columns of V are now approximate Schur  |
c        | vectors associated with the real upper quasi-triangular |
c        | matrix of order NCONV in workl(iuptri)                  |
c        %---------------------------------------------------------%
c     
         call dorm2r ('Right', 'Notranspose', n            , 
     &                ncv   , nconv        , workl(invsub),
     &                ldq   , workev       , v            , 
     &                ldv   , workd(n+1)   , ierr)
         call dlacpy ('All', n, nconv, v, ldv, z, ldz)
c
         do 20 j=1, nconv
c     
c           %---------------------------------------------------%
c           | Perform both a column and row scaling if the      |
c           | diagonal element of workl(invsub,ldq) is negative |
c           | I'm lazy and don't take advantage of the upper    |
c           | quasi-triangular form of workl(iuptri,ldq)        |
c           | Note that since Q is orthogonal, R is a diagonal  |
c           | matrix consisting of plus or minus ones           |
c           %---------------------------------------------------%
c     
            if (workl(invsub+(j-1)*ldq+j-1) .lt. zero) then
               call dscal (nconv, -one, workl(iuptri+j-1), ldq)
               call dscal (nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if
c     
 20      continue
c     
         if (howmny .eq. 'A') then
c     
c           %--------------------------------------------%
c           | Compute the NCONV wanted eigenvectors of T | 
c           | located in workl(iuptri,ldq).              |
c           %--------------------------------------------%
c     
            do 30 j=1, ncv
               if (j .le. nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
c
            call dtrevc ('Right', 'Select'     , select       , 
     &                   ncv    , workl(iuptri), ldq          , 
     &                   vl     , 1            , workl(invsub),
     &                   ldq    , ncv          , outncv       ,
     &                   workev , ierr)
c
            if (ierr .ne. 0) then
                info = -9
                go to 9000
            end if
c     
c           %------------------------------------------------%
c           | Scale the returning eigenvectors so that their |
c           | Euclidean norms are all one. LAPACK subroutine |
c           | dtrevc  returns each eigenvector normalized so  |
c           | that the element of largest magnitude has      |
c           | magnitude 1;                                   |
c           %------------------------------------------------%
c     
            iconj = 0
            do 40 j=1, nconv
c
               if ( workl(iheigi+j-1) .eq. zero ) then
c     
c                 %----------------------%
c                 | real eigenvalue case |
c                 %----------------------%
c     
                  temp = dnrm2 ( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call dscal ( ncv, one / temp, 
     &                 workl(invsub+(j-1)*ldq), 1 )
c
               else
c     
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 | columns, we further normalize by the      |
c                 | square root of two.                       |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     temp = dlapy2 (dnrm2 (ncv, 
     &                                   workl(invsub+(j-1)*ldq), 
     &                                   1),
     &                             dnrm2 (ncv, 
     &                                   workl(invsub+j*ldq),
     &                                   1))  
                     call dscal (ncv, one/temp, 
     &                           workl(invsub+(j-1)*ldq), 1 )
                     call dscal (ncv, one/temp, 
     &                           workl(invsub+j*ldq), 1 )
                     iconj = 1
                  else
                     iconj = 0
                  end if
c
               end if
c
 40         continue
c
            call dgemv ('T', ncv, nconv, one, workl(invsub),
     &                 ldq, workl(ihbds), 1, zero,  workev, 1)
c
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) .ne. zero) then
c
c                 %-------------------------------------------%
c                 | Complex conjugate pair case. Note that    |
c                 | since the real and imaginary part of      |
c                 | the eigenvector are stored in consecutive |
c                 %-------------------------------------------%
c
                  if (iconj .eq. 0) then
                     workev(j) = dlapy2 (workev(j), workev(j+1))
                     workev(j+1) = workev(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
 45         continue
c
            if (msglvl .gt. 2) then
               call dcopy (ncv, workl(invsub+ncv-1), ldq,
     &                    workl(ihbds), 1)
               call dvout (logfil, ncv, workl(ihbds), ndigit,
     &              '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl .gt. 3) then
                  call dmout (logfil, ncv, ncv, workl(invsub), ldq, 
     &                 ndigit, '_neupd: The eigenvector matrix for T')
               end if
            end if
c
c           %---------------------------------------%
c           | Copy Ritz estimates into workl(ihbds) |
c           %---------------------------------------%
c
            call dcopy (nconv, workev, 1, workl(ihbds), 1)
c
c           %---------------------------------------------------------%
c           | Compute the QR factorization of the eigenvector matrix  |
c           | associated with leading portion of T in the first NCONV |
c           | columns of workl(invsub,ldq).                           |
c           %---------------------------------------------------------%
c     
            call dgeqr2 (ncv, nconv , workl(invsub), 
     &                   ldq, workev, workev(ncv+1),
     &                   ierr)
c     
c           %----------------------------------------------%
c           | * Postmultiply Z by Q.                       |   
c           | * Postmultiply Z by R.                       |
c           | The N by NCONV matrix Z is now contains the  | 
c           | Ritz vectors associated with the Ritz values |
c           | in workl(iheigr) and workl(iheigi).          |
c           %----------------------------------------------%
c     
            call dorm2r ('Right', 'Notranspose', n            ,
     &                   ncv  , nconv        , workl(invsub),
     &                   ldq  , workev       , z            ,
     &                   ldz  , workd(n+1)   , ierr)
c     
            call dtrmm ('Right'   , 'Upper'       , 'No transpose',
     &                  'Non-unit', n            , nconv         ,
     &                  one       , workl(invsub), ldq           ,
     &                  z         , ldz)
c     
         end if
c     
      else 
c
c        %------------------------------------------------------%
c        | An approximate invariant subspace is not needed.     |
c        | Place the Ritz values computed DNAUPD  into DR and DI |
c        %------------------------------------------------------%
c
         call dcopy (nconv, workl(ritzr), 1, dr, 1)
         call dcopy (nconv, workl(ritzi), 1, di, 1)
         call dcopy (nconv, workl(ritzr), 1, workl(iheigr), 1)
         call dcopy (nconv, workl(ritzi), 1, workl(iheigi), 1)
         call dcopy (nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
c 
c     %------------------------------------------------%
c     | Transform the Ritz values and possibly vectors |
c     | and corresponding error bounds of OP to those  |
c     | of A*x = lambda*B*x.                           |
c     %------------------------------------------------%
c
      if (type .eq. 'REGULR') then
c
         if (rvec) 
     &      call dscal (ncv, rnorm, workl(ihbds), 1)     
c     
      else 
c     
c        %---------------------------------------%
c        |   A spectral transformation was used. |
c        | * Determine the Ritz estimates of the |
c        |   Ritz values in the original system. |
c        %---------------------------------------%
c     
         if (type .eq. 'SHIFTI') then
c
            if (rvec) 
     &         call dscal (ncv, rnorm, workl(ihbds), 1)
c
            do 50 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(ihbds+k-1) = abs( workl(ihbds+k-1) ) 
     &                          / temp / temp
 50         continue
c
         else if (type .eq. 'REALPT') then
c
            do 60 k=1, ncv
 60         continue
c
         else if (type .eq. 'IMAGPT') then
c
            do 70 k=1, ncv
 70         continue
c
         end if
c     
c        %-----------------------------------------------------------%
c        | *  Transform the Ritz values back to the original system. |
c        |    For TYPE = 'SHIFTI' the transformation is              |
c        |             lambda = 1/theta + sigma                      |
c        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
c        |    Rayleigh quotients or a projection. See remark 3 above.| 
c        | NOTES:                                                    |
c        | *The Ritz vectors are not affected by the transformation. |
c        %-----------------------------------------------------------%
c     
         if (type .eq. 'SHIFTI') then 
c
            do 80 k=1, ncv
               temp = dlapy2 ( workl(iheigr+k-1), 
     &                        workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1)/temp/temp 
     &                           + sigmar   
               workl(iheigi+k-1) = -workl(iheigi+k-1)/temp/temp
     &                           + sigmai   
 80         continue
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         else if (type .eq. 'REALPT' .or. type .eq. 'IMAGPT') then
c
            call dcopy (nconv, workl(iheigr), 1, dr, 1)
            call dcopy (nconv, workl(iheigi), 1, di, 1)
c
         end if
c
      end if
c
      if (type .eq. 'SHIFTI' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, dr, ndigit,
     &   '_neupd: Untransformed real part of the Ritz valuess.')
         call dvout  (logfil, nconv, di, ndigit,
     &   '_neupd: Untransformed imag part of the Ritz valuess.')
         call dvout (logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Ritz estimates of untransformed Ritz values.')
      else if (type .eq. 'REGULR' .and. msglvl .gt. 1) then
         call dvout (logfil, nconv, dr, ndigit,
     &   '_neupd: Real parts of converged Ritz values.')
         call dvout  (logfil, nconv, di, ndigit,
     &   '_neupd: Imag parts of converged Ritz values.')
         call dvout (logfil, nconv, workl(ihbds), ndigit,
     &   '_neupd: Associated Ritz estimates.')
      end if
c 
c     %-------------------------------------------------%
c     | Eigenvector Purification step. Formally perform |
c     | one of inverse subspace iteration. Only used    |
c     | for MODE = 2.                                   |
c     %-------------------------------------------------%
c
      if (rvec .and. howmny .eq. 'A' .and. type .eq. 'SHIFTI') then
c
c        %------------------------------------------------%
c        | Purify the computed Ritz vectors by adding a   |
c        | little bit of the residual vector:             |
c        |                      T                         |
c        |          resid(:)*( e    s ) / theta           |
c        |                      NCV                       |
c        | where H s = s theta. Remember that when theta  |
c        | has nonzero imaginary part, the corresponding  |
c        | Ritz vector is stored across two columns of Z. |
c        %------------------------------------------------%
c
         iconj = 0
         do 110 j=1, nconv
            if (workl(iheigi+j-1) .eq. zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) /
     &                      workl(iheigr+j-1)
            else if (iconj .eq. 0) then
               temp = dlapy2 ( workl(iheigr+j-1), workl(iheigi+j-1) )
               workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) * 
     &                       workl(iheigr+j-1) +
     &                       workl(invsub+j*ldq+ncv-1) * 
     &                       workl(iheigi+j-1) ) / temp / temp
               workev(j+1) = ( workl(invsub+j*ldq+ncv-1) * 
     &                         workl(iheigr+j-1) -
     &                         workl(invsub+(j-1)*ldq+ncv-1) * 
     &                         workl(iheigi+j-1) ) / temp / temp
               iconj = 1
            else
               iconj = 0
            end if
 110     continue
c
c        %---------------------------------------%
c        | Perform a rank one update to Z and    |
c        | purify all the Ritz vectors together. |
c        %---------------------------------------%
c
         call dger (n, nconv, one, resid, 1, workev, 1, z, ldz)
c
      end if
c
 9000 continue
c
      return
c     
c     %---------------%
c     | End of DNEUPD  |
c     %---------------%
c
      end
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dngets
c
c\Description: 
c  Given the eigenvalues of the upper Hessenberg matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of 
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: call this even in the case of user specified shifts in order
c  to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dngets
c     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> want the KEV eigenvalues of largest magnitude.
c          'SM' -> want the KEV eigenvalues of smallest magnitude.
c          'LR' -> want the KEV eigenvalues of largest real part.
c          'SR' -> want the KEV eigenvalues of smallest real part.
c          'LI' -> want the KEV eigenvalues of largest imaginary part.
c          'SI' -> want the KEV eigenvalues of smallest imaginary part.
c
c  KEV      Integer.  (INPUT/OUTPUT)
c           INPUT: KEV+NP is the size of the matrix H.
c           OUTPUT: Possibly increases KEV by one to keep complex conjugate
c           pairs together.
c
c  NP       Integer.  (INPUT/OUTPUT)
c           Number of implicit shifts to be computed.
c           OUTPUT: Possibly decreases NP by one to keep complex conjugate
c           pairs together.
c
c  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary 
c          parts of the eigenvalues of H.
c          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
c          eigenvalues are in the first NP locations and the wanted
c          portion is in the last KEV locations.  When exact shifts are 
c          selected, the unwanted part corresponds to the shifts to 
c          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
c          are further sorted so that the ones with largest Ritz values
c          are first.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
c  
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortc  ARPACK sorting routine.
c     dcopy   Level 1 BLAS that copies one vector to another .
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\Remarks
c     1. xxxx
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dngets ( ishift, which, kev, np, ritzr, ritzi, bounds,
     &                    shiftr, shifti )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           bounds(kev+np), ritzr(kev+np), ritzi(kev+np), 
     &           shiftr(1), shifti(1)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0, zero = 0.0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dsortc, second
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c 
      call second (t0)
      msglvl = mngets
c 
c     %----------------------------------------------------%
c     | LM, SM, LR, SR, LI, SI case.                       |
c     | Sort the eigenvalues of H into the desired order   |
c     | and apply the resulting order to BOUNDS.           |
c     | The eigenvalues are sorted so that the wanted part |
c     | are always in the last KEV locations.              |
c     | We first do a pre-processing sort in order to keep |
c     | complex conjugate pairs together                   |
c     %----------------------------------------------------%
c
      if (which .eq. 'LM') then
         call dsortc ('LR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SM') then
         call dsortc ('SR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LR') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SR') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'LI') then
         call dsortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which .eq. 'SI') then
         call dsortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      end if
c      
      call dsortc (which, .true., kev+np, ritzr, ritzi, bounds)
c     
c     %-------------------------------------------------------%
c     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
c     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero |
c     | Accordingly decrease NP by one. In other words keep   |
c     | complex conjugate pairs together.                     |
c     %-------------------------------------------------------%
c     
      if (       ( ritzr(np+1) - ritzr(np) ) .eq. zero
     &     .and. ( ritzi(np+1) + ritzi(np) ) .eq. zero ) then
         np = np - 1
         kev = kev + 1
      end if
c
      if ( ishift .eq. 1 ) then
c     
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first        |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when they shifts |
c        | are applied in subroutine dnapps.                     |
c        | Be careful and use 'SR' since we want to sort BOUNDS! |
c        %-------------------------------------------------------%
c     
         call dsortc ( 'SR', .true., np, bounds, ritzr, ritzi )
      end if
c     
      call second (t1)
      tngets = tngets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, kev, ndigit, '_ngets: KEV is')
         call ivout (logfil, 1, np, ndigit, '_ngets: NP is')
         call dvout (logfil, kev+np, ritzr, ndigit,
     &        '_ngets: Eigenvalues of current H matrix -- real part')
         call dvout (logfil, kev+np, ritzi, ndigit,
     &        '_ngets: Eigenvalues of current H matrix -- imag part')
         call dvout (logfil, kev+np, bounds, ndigit, 
     &      '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
c     
      return
c     
c     %---------------%
c     | End of dngets |
c     %---------------%
c     
      end
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORG2R generates an m by n real matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of DORG2R
*
      END
      SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGHR generates a real orthogonal matrix Q which is defined as the
*  product of IHI-ILO elementary reflectors of order N, as returned by
*  DGEHRD:
*
*  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI must have the same values as in the previous call
*          of DGEHRD. Q is equal to the unit matrix except in the
*          submatrix Q(ilo+1:ihi,ilo+1:ihi).
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by DGEHRD.
*          On exit, the N-by-N orthogonal matrix Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEHRD.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= IHI-ILO.
*          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORGQR, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NH = IHI - ILO
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, NH ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DORGQR', ' ', NH, NH, NH, -1 )
         LWKOPT = MAX( 1, NH )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Shift the vectors which define the elementary reflectors one
*     column to the right, and set the first ilo and the last n-ihi
*     rows and columns to those of the unit matrix
*
      DO 40 J = IHI, ILO + 1, -1
         DO 10 I = 1, J - 1
            A( I, J ) = ZERO
   10    CONTINUE
         DO 20 I = J + 1, IHI
            A( I, J ) = A( I, J-1 )
   20    CONTINUE
         DO 30 I = IHI + 1, N
            A( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      DO 60 J = 1, ILO
         DO 50 I = 1, N
            A( I, J ) = ZERO
   50    CONTINUE
         A( J, J ) = ONE
   60 CONTINUE
      DO 80 J = IHI + 1, N
         DO 70 I = 1, N
            A( I, J ) = ZERO
   70    CONTINUE
         A( J, J ) = ONE
   80 CONTINUE
*
      IF( NH.GT.0 ) THEN
*
*        Generate Q(ilo+1:ihi,ilo+1:ihi)
*
         CALL DORGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ),
     $                WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORGHR
*
      END
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORM2R overwrites the general real m by n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'T', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'T',
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'T': apply Q' (Transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of DORM2R
*
      END
      SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORMHR overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix of order nq, with nq = m if
*  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
*  IHI-ILO elementary reflectors, as returned by DGEHRD:
*
*  Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI must have the same values as in the previous call
*          of DGEHRD. Q is equal to the unit matrix except in the
*          submatrix Q(ilo+1:ihi,ilo+1:ihi).
*          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
*          ILO = 1 and IHI = 0, if M = 0;
*          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
*          ILO = 1 and IHI = 0, if N = 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                               (LDA,M) if SIDE = 'L'
*                               (LDA,N) if SIDE = 'R'
*          The vectors which define the elementary reflectors, as
*          returned by DGEHRD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
*
*  TAU     (input) DOUBLE PRECISION array, dimension
*                               (M-1) if SIDE = 'L'
*                               (N-1) if SIDE = 'R'
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEHRD.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DORMQR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NH = IHI - ILO
      LEFT = LSAME( SIDE, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) )
     $          THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, NQ ) ) THEN
         INFO = -5
      ELSE IF( IHI.LT.MIN( ILO, NQ ) .OR. IHI.GT.NQ ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( LEFT ) THEN
            NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, NH, N, NH, -1 )
         ELSE
            NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, NH, NH, -1 )
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMHR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NH.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( LEFT ) THEN
         MI = NH
         NI = N
         I1 = ILO + 1
         I2 = 1
      ELSE
         MI = M
         NI = NH
         I1 = 1
         I2 = ILO + 1
      END IF
*
      CALL DORMQR( SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA,
     $             TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO )
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMHR
*
      END
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DORMQR overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,
     $                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of DORMQR
*
      END
      SUBROUTINE DRSCL( N, SA, SX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
*     ..
*
*  Purpose
*  =======
*
*  DRSCL multiplies an n-element real vector x by the real scalar 1/a.
*  This is done without overflow or underflow as long as
*  the final result x/a does not overflow or underflow.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of components of the vector x.
*
*  SA      (input) DOUBLE PRECISION
*          The scalar a which is used to divide each component of x.
*          SA must be >= 0, or the subroutine will divide by zero.
*
*  SX      (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-1)*abs(INCX))
*          The n-element vector x.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector SX.
*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      DOUBLE PRECISION   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Initialize the denominator to SA and the numerator to 1.
*
      CDEN = SA
      CNUM = ONE
*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
*
*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
*
*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
*
*        Multiply X by CNUM / CDEN and return.
*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
*
*     Scale the vector X by MUL
*
      CALL DSCAL( N, MUL, SX, INCX )
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DRSCL
*
      END
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortc
c
c\Description:
c  Sorts the complex array in XREAL and XIMAG into the order 
c  specified by WHICH and optionally applies the permutation to the
c  real array Y. It is assumed that if an element of XIMAG is
c  nonzero, then its negative is also an element. In other words,
c  both members of a complex conjugate pair are to be sorted and the
c  pairs are kept adjacent to each other.
c
c\Usage:
c  call dsortc
c     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
c          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
c          'LR' -> sort XREAL into increasing order of algebraic.
c          'SR' -> sort XREAL into decreasing order of algebraic.
c          'LI' -> sort XIMAG into increasing order of magnitude.
c          'SI' -> sort XIMAG into decreasing order of magnitude.
c          NOTE: If an element of XIMAG is non-zero, then its negative
c                is also an element.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to array Y.
c          APPLY = .FALSE. -> do not apply the sorted order to array Y.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  XREAL,  Double precision array of length N.  (INPUT/OUTPUT)
c  XIMAG   Real and imaginary part of the array to be sorted.
c
c  Y       Double precision array of length N.  (INPUT/OUTPUT)
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/92: Version ' 2.1'
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#) 
c FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortc (which, apply, n, xreal, ximag, y)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision     
     &           xreal(0:n-1), ximag(0:n-1), y(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision     
     &           temp, temp1, temp2
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision     
     &           dlapy2
      external   dlapy2
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c 
      if (which .eq. 'LM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into increasing order of magnitude. |
c        %------------------------------------------------------%
c
   10    continue
         if (igap .eq. 0) go to 9000
c
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.gt.temp2) then
                temp = xreal(j)
                xreal(j) = xreal(j+igap)
                xreal(j+igap) = temp
c
                temp = ximag(j)
                ximag(j) = ximag(j+igap)
                ximag(j+igap) = temp
c
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        %------------------------------------------------------%
c        | Sort XREAL,XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------------%
c
   40    continue
         if (igap .eq. 0) go to 9000
c
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j .lt. 0) go to 60
c
            temp1 = dlapy2(xreal(j),ximag(j))
            temp2 = dlapy2(xreal(j+igap),ximag(j+igap))
c
            if (temp1.lt.temp2) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c 
      else if (which .eq. 'LR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into increasing order of algebraic. |
c        %------------------------------------------------%
c
   70    continue
         if (igap .eq. 0) go to 9000
c
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (xreal(j).gt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c 
      else if (which .eq. 'SR') then
c
c        %------------------------------------------------%
c        | Sort XREAL into decreasing order of algebraic. |
c        %------------------------------------------------%
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (xreal(j).lt.xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
c 
      else if (which .eq. 'LI') then
c
c        %------------------------------------------------%
c        | Sort XIMAG into increasing order of magnitude. |
c        %------------------------------------------------%
c
  130    continue
         if (igap .eq. 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
c
            if (j.lt.0) go to 150
c
            if (abs(ximag(j)).gt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
c 
      else if (which .eq. 'SI') then
c
c        %------------------------------------------------%
c        | Sort XIMAG into decreasing order of magnitude. |
c        %------------------------------------------------%
c
  160    continue
         if (igap .eq. 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
c
            if (j.lt.0) go to 180
c
            if (abs(ximag(j)).lt.abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
c
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
c 
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
c 
 9000 continue
      return
c
c     %---------------%
c     | End of dsortc |
c     %---------------%
c
      end
c
c     %---------------------------------------------%
c     | Initialize statistic and timing information |
c     | for nonsymmetric Arnoldi code.              |
c     %---------------------------------------------%
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas    
c
c\SCCS Information: @(#) 
c FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2
c
      subroutine dstatn
c
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
      include   'stat.h'
c 
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
c 
      tnaupd = 0.0D+0
      tnaup2 = 0.0D+0
      tnaitr = 0.0D+0
      tneigh = 0.0D+0
      tngets = 0.0D+0
      tnapps = 0.0D+0
      tnconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
c 
c     %----------------------------------------------------%
c     | User time including reverse communication overhead |
c     %----------------------------------------------------%
c
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
c 
      return
c
c
c     %---------------%
c     | End of dstatn |
c     %---------------%
c
      end
      SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, MM, M, WORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DTREVC computes some or all of the right and/or left eigenvectors of
*  a real upper quasi-triangular matrix T.
*  Matrices of this type are produced by the Schur factorization of
*  a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.
*  
*  The right eigenvector x and the left eigenvector y of T corresponding
*  to an eigenvalue w are defined by:
*  
*     T*x = w*x,     (y**H)*T = w*(y**H)
*  
*  where y**H denotes the conjugate transpose of y.
*  The eigenvalues are not input to this routine, but are read directly
*  from the diagonal blocks of T.
*  
*  This routine returns the matrices X and/or Y of right and left
*  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
*  input matrix.  If Q is the orthogonal factor that reduces a matrix
*  A to Schur form T, then Q*X and Q*Y are the matrices of right and
*  left eigenvectors of A.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  compute right eigenvectors only;
*          = 'L':  compute left eigenvectors only;
*          = 'B':  compute both right and left eigenvectors.
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A':  compute all right and/or left eigenvectors;
*          = 'B':  compute all right and/or left eigenvectors,
*                  backtransformed by the matrices in VR and/or VL;
*          = 'S':  compute selected right and/or left eigenvectors,
*                  as indicated by the logical array SELECT.
*
*  SELECT  (input/output) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*          computed.
*          If w(j) is a real eigenvalue, the corresponding real
*          eigenvector is computed if SELECT(j) is .TRUE..
*          If w(j) and w(j+1) are the real and imaginary parts of a
*          complex eigenvalue, the corresponding complex eigenvector is
*          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
*          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
*          .FALSE..
*          Not referenced if HOWMNY = 'A' or 'B'.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          The upper quasi-triangular matrix T in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of Schur vectors returned by DHSEQR).
*          On exit, if SIDE = 'L' or 'B', VL contains:
*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*Y;
*          if HOWMNY = 'S', the left eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VL, in the same order as their
*                           eigenvalues.
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part, and the second the imaginary part.
*          Not referenced if SIDE = 'R'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.  LDVL >= 1, and if
*          SIDE = 'L' or 'B', LDVL >= N.
*
*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of Schur vectors returned by DHSEQR).
*          On exit, if SIDE = 'R' or 'B', VR contains:
*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*X;
*          if HOWMNY = 'S', the right eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VR, in the same order as their
*                           eigenvalues.
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part and the second the imaginary part.
*          Not referenced if SIDE = 'L'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= 1, and if
*          SIDE = 'R' or 'B', LDVR >= N.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR actually
*          used to store the eigenvectors.
*          If HOWMNY = 'A' or 'B', M is set to N.
*          Each selected real eigenvector occupies one column and each
*          selected complex eigenvector occupies two columns.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The algorithm used in this program is basically backward (forward)
*  substitution, with scaling to make the the code robust against
*  possible overflow.
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x| + |y|.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV
      INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, N2
      DOUBLE PRECISION   BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE,
     $                   SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR,
     $                   XNORM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DLALN2, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   X( 2, 2 )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
*
      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )
*
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE
*
*        Set M to the number of columns required to store the selected
*        eigenvectors, standardize the array SELECT if necessary, and
*        test MM.
*
         IF( SOMEV ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 J = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
                  SELECT( J ) = .FALSE.
               ELSE
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).EQ.ZERO ) THEN
                        IF( SELECT( J ) )
     $                     M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( J ) .OR. SELECT( J+1 ) ) THEN
                           SELECT( J ) = .TRUE.
                           M = M + 2
                        END IF
                     END IF
                  ELSE
                     IF( SELECT( N ) )
     $                  M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
*
         IF( MM.LT.M ) THEN
            INFO = -11
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Set the constants to control overflow.
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      WORK( 1 ) = ZERO
      DO 30 J = 2, N
         WORK( J ) = ZERO
         DO 20 I = 1, J - 1
            WORK( J ) = WORK( J ) + ABS( T( I, J ) )
   20    CONTINUE
   30 CONTINUE
*
*     Index IP is used to specify the real or complex eigenvalue:
*       IP = 0, real eigenvalue,
*            1, first of conjugate complex pair: (wr,wi)
*           -1, second of conjugate complex pair: (wr,wi)
*
      N2 = 2*N
*
      IF( RIGHTV ) THEN
*
*        Compute right eigenvectors.
*
         IP = 0
         IS = M
         DO 140 KI = N, 1, -1
*
            IF( IP.EQ.1 )
     $         GO TO 130
            IF( KI.EQ.1 )
     $         GO TO 40
            IF( T( KI, KI-1 ).EQ.ZERO )
     $         GO TO 40
            IP = -1
*
   40       CONTINUE
            IF( SOMEV ) THEN
               IF( IP.EQ.0 ) THEN
                  IF( .NOT.SELECT( KI ) )
     $               GO TO 130
               ELSE
                  IF( .NOT.SELECT( KI-1 ) )
     $               GO TO 130
               END IF
            END IF
*
*           Compute the KI-th eigenvalue (WR,WI).
*
            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 )
     $         WI = SQRT( ABS( T( KI, KI-1 ) ) )*
     $              SQRT( ABS( T( KI-1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
*
            IF( IP.EQ.0 ) THEN
*
*              Real right eigenvector
*
               WORK( KI+N ) = ONE
*
*              Form right-hand side
*
               DO 50 K = 1, KI - 1
                  WORK( K+N ) = -T( K, KI )
   50          CONTINUE
*
*              Solve the upper quasi-triangular system:
*                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
*
               JNXT = KI - 1
               DO 60 J = KI - 1, 1, -1
                  IF( J.GT.JNXT )
     $               GO TO 60
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J.GT.1 ) THEN
                     IF( T( J, J-1 ).NE.ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
*
                  IF( J1.EQ.J2 ) THEN
*
*                    1-by-1 diagonal block
*
                     CALL DLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ),
     $                            LDT, ONE, ONE, WORK( J+N ), N, WR,
     $                            ZERO, X, 2, SCALE, XNORM, IERR )
*
*                    Scale X(1,1) to avoid overflow when updating
*                    the right-hand side.
*
                     IF( XNORM.GT.ONE ) THEN
                        IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE )
     $                  CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
*
*                    Update right-hand side
*
                     CALL DAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1,
     $                           WORK( 1+N ), 1 )
*
                  ELSE
*
*                    2-by-2 diagonal block
*
                     CALL DLALN2( .FALSE., 2, 1, SMIN, ONE,
     $                            T( J-1, J-1 ), LDT, ONE, ONE,
     $                            WORK( J-1+N ), N, WR, ZERO, X, 2,
     $                            SCALE, XNORM, IERR )
*
*                    Scale X(1,1) and X(2,1) to avoid overflow when
*                    updating the right-hand side.
*
                     IF( XNORM.GT.ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA.GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 2, 1 ) = X( 2, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE )
     $                  CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
*
*                    Update right-hand side
*
                     CALL DAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1,
     $                           WORK( 1+N ), 1 )
                     CALL DAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1,
     $                           WORK( 1+N ), 1 )
                  END IF
   60          CONTINUE
*
*              Copy the vector x or Q*x to VR and normalize.
*
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( KI, WORK( 1+N ), 1, VR( 1, IS ), 1 )
*
                  II = IDAMAX( KI, VR( 1, IS ), 1 )
                  REMAX = ONE / ABS( VR( II, IS ) )
                  CALL DSCAL( KI, REMAX, VR( 1, IS ), 1 )
*
                  DO 70 K = KI + 1, N
                     VR( K, IS ) = ZERO
   70             CONTINUE
               ELSE
                  IF( KI.GT.1 )
     $               CALL DGEMV( 'N', N, KI-1, ONE, VR, LDVR,
     $                           WORK( 1+N ), 1, WORK( KI+N ),
     $                           VR( 1, KI ), 1 )
*
                  II = IDAMAX( N, VR( 1, KI ), 1 )
                  REMAX = ONE / ABS( VR( II, KI ) )
                  CALL DSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
*
            ELSE
*
*              Complex right eigenvector.
*
*              Initial solve
*                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
*                [ (T(KI,KI-1)   T(KI,KI)   )               ]
*
               IF( ABS( T( KI-1, KI ) ).GE.ABS( T( KI, KI-1 ) ) ) THEN
                  WORK( KI-1+N ) = ONE
                  WORK( KI+N2 ) = WI / T( KI-1, KI )
               ELSE
                  WORK( KI-1+N ) = -WI / T( KI, KI-1 )
                  WORK( KI+N2 ) = ONE
               END IF
               WORK( KI+N ) = ZERO
               WORK( KI-1+N2 ) = ZERO
*
*              Form right-hand side
*
               DO 80 K = 1, KI - 2
                  WORK( K+N ) = -WORK( KI-1+N )*T( K, KI-1 )
                  WORK( K+N2 ) = -WORK( KI+N2 )*T( K, KI )
   80          CONTINUE
*
*              Solve upper quasi-triangular system:
*              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
*
               JNXT = KI - 2
               DO 90 J = KI - 2, 1, -1
                  IF( J.GT.JNXT )
     $               GO TO 90
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J.GT.1 ) THEN
                     IF( T( J, J-1 ).NE.ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
*
                  IF( J1.EQ.J2 ) THEN
*
*                    1-by-1 diagonal block
*
                     CALL DLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ),
     $                            LDT, ONE, ONE, WORK( J+N ), N, WR, WI,
     $                            X, 2, SCALE, XNORM, IERR )
*
*                    Scale X(1,1) and X(1,2) to avoid overflow when
*                    updating the right-hand side.
*
                     IF( XNORM.GT.ONE ) THEN
                        IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 1, 2 ) = X( 1, 2 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL DSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
*
*                    Update the right-hand side
*
                     CALL DAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1,
     $                           WORK( 1+N ), 1 )
                     CALL DAXPY( J-1, -X( 1, 2 ), T( 1, J ), 1,
     $                           WORK( 1+N2 ), 1 )
*
                  ELSE
*
*                    2-by-2 diagonal block
*
                     CALL DLALN2( .FALSE., 2, 2, SMIN, ONE,
     $                            T( J-1, J-1 ), LDT, ONE, ONE,
     $                            WORK( J-1+N ), N, WR, WI, X, 2, SCALE,
     $                            XNORM, IERR )
*
*                    Scale X to avoid overflow when updating
*                    the right-hand side.
*
                     IF( XNORM.GT.ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA.GT.BIGNUM / XNORM ) THEN
                           REC = ONE / XNORM
                           X( 1, 1 ) = X( 1, 1 )*REC
                           X( 1, 2 ) = X( 1, 2 )*REC
                           X( 2, 1 ) = X( 2, 1 )*REC
                           X( 2, 2 ) = X( 2, 2 )*REC
                           SCALE = SCALE*REC
                        END IF
                     END IF
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL DSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
                     WORK( J-1+N2 ) = X( 1, 2 )
                     WORK( J+N2 ) = X( 2, 2 )
*
*                    Update the right-hand side
*
                     CALL DAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1,
     $                           WORK( 1+N ), 1 )
                     CALL DAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1,
     $                           WORK( 1+N ), 1 )
                     CALL DAXPY( J-2, -X( 1, 2 ), T( 1, J-1 ), 1,
     $                           WORK( 1+N2 ), 1 )
                     CALL DAXPY( J-2, -X( 2, 2 ), T( 1, J ), 1,
     $                           WORK( 1+N2 ), 1 )
                  END IF
   90          CONTINUE
*
*              Copy the vector x or Q*x to VR and normalize.
*
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( KI, WORK( 1+N ), 1, VR( 1, IS-1 ), 1 )
                  CALL DCOPY( KI, WORK( 1+N2 ), 1, VR( 1, IS ), 1 )
*
                  EMAX = ZERO
                  DO 100 K = 1, KI
                     EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+
     $                      ABS( VR( K, IS ) ) )
  100             CONTINUE
*
                  REMAX = ONE / EMAX
                  CALL DSCAL( KI, REMAX, VR( 1, IS-1 ), 1 )
                  CALL DSCAL( KI, REMAX, VR( 1, IS ), 1 )
*
                  DO 110 K = KI + 1, N
                     VR( K, IS-1 ) = ZERO
                     VR( K, IS ) = ZERO
  110             CONTINUE
*
               ELSE
*
                  IF( KI.GT.2 ) THEN
                     CALL DGEMV( 'N', N, KI-2, ONE, VR, LDVR,
     $                           WORK( 1+N ), 1, WORK( KI-1+N ),
     $                           VR( 1, KI-1 ), 1 )
                     CALL DGEMV( 'N', N, KI-2, ONE, VR, LDVR,
     $                           WORK( 1+N2 ), 1, WORK( KI+N2 ),
     $                           VR( 1, KI ), 1 )
                  ELSE
                     CALL DSCAL( N, WORK( KI-1+N ), VR( 1, KI-1 ), 1 )
                     CALL DSCAL( N, WORK( KI+N2 ), VR( 1, KI ), 1 )
                  END IF
*
                  EMAX = ZERO
                  DO 120 K = 1, N
                     EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+
     $                      ABS( VR( K, KI ) ) )
  120             CONTINUE
                  REMAX = ONE / EMAX
                  CALL DSCAL( N, REMAX, VR( 1, KI-1 ), 1 )
                  CALL DSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
            END IF
*
            IS = IS - 1
            IF( IP.NE.0 )
     $         IS = IS - 1
  130       CONTINUE
            IF( IP.EQ.1 )
     $         IP = 0
            IF( IP.EQ.-1 )
     $         IP = 1
  140    CONTINUE
      END IF
*
      IF( LEFTV ) THEN
*
*        Compute left eigenvectors.
*
         IP = 0
         IS = 1
         DO 260 KI = 1, N
*
            IF( IP.EQ.-1 )
     $         GO TO 250
            IF( KI.EQ.N )
     $         GO TO 150
            IF( T( KI+1, KI ).EQ.ZERO )
     $         GO TO 150
            IP = 1
*
  150       CONTINUE
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) )
     $            GO TO 250
            END IF
*
*           Compute the KI-th eigenvalue (WR,WI).
*
            WR = T( KI, KI )
            WI = ZERO
            IF( IP.NE.0 )
     $         WI = SQRT( ABS( T( KI, KI+1 ) ) )*
     $              SQRT( ABS( T( KI+1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
*
            IF( IP.EQ.0 ) THEN
*
*              Real left eigenvector.
*
               WORK( KI+N ) = ONE
*
*              Form right-hand side
*
               DO 160 K = KI + 1, N
                  WORK( K+N ) = -T( KI, K )
  160          CONTINUE
*
*              Solve the quasi-triangular system:
*                 (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
*
               VMAX = ONE
               VCRIT = BIGNUM
*
               JNXT = KI + 1
               DO 170 J = KI + 1, N
                  IF( J.LT.JNXT )
     $               GO TO 170
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).NE.ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
*
                  IF( J1.EQ.J2 ) THEN
*
*                    1-by-1 diagonal block
*
*                    Scale if necessary to avoid overflow when forming
*                    the right-hand side.
*
                     IF( WORK( J ).GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
*
                     WORK( J+N ) = WORK( J+N ) -
     $                             DDOT( J-KI-1, T( KI+1, J ), 1,
     $                             WORK( KI+1+N ), 1 )
*
*                    Solve (T(J,J)-WR)'*X = WORK
*
                     CALL DLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ),
     $                            LDT, ONE, ONE, WORK( J+N ), N, WR,
     $                            ZERO, X, 2, SCALE, XNORM, IERR )
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE )
     $                  CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     VMAX = MAX( ABS( WORK( J+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
*
                  ELSE
*
*                    2-by-2 diagonal block
*
*                    Scale if necessary to avoid overflow when forming
*                    the right-hand side.
*
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA.GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
*
                     WORK( J+N ) = WORK( J+N ) -
     $                             DDOT( J-KI-1, T( KI+1, J ), 1,
     $                             WORK( KI+1+N ), 1 )
*
                     WORK( J+1+N ) = WORK( J+1+N ) -
     $                               DDOT( J-KI-1, T( KI+1, J+1 ), 1,
     $                               WORK( KI+1+N ), 1 )
*
*                    Solve
*                      [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
*                      [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
*
                     CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, T( J, J ),
     $                            LDT, ONE, ONE, WORK( J+N ), N, WR,
     $                            ZERO, X, 2, SCALE, XNORM, IERR )
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE )
     $                  CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+1+N ) = X( 2, 1 )
*
                     VMAX = MAX( ABS( WORK( J+N ) ),
     $                      ABS( WORK( J+1+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
*
                  END IF
  170          CONTINUE
*
*              Copy the vector x or Q*x to VL and normalize.
*
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
*
                  II = IDAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                  REMAX = ONE / ABS( VL( II, IS ) )
                  CALL DSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
*
                  DO 180 K = 1, KI - 1
                     VL( K, IS ) = ZERO
  180             CONTINUE
*
               ELSE
*
                  IF( KI.LT.N )
     $               CALL DGEMV( 'N', N, N-KI, ONE, VL( 1, KI+1 ), LDVL,
     $                           WORK( KI+1+N ), 1, WORK( KI+N ),
     $                           VL( 1, KI ), 1 )
*
                  II = IDAMAX( N, VL( 1, KI ), 1 )
                  REMAX = ONE / ABS( VL( II, KI ) )
                  CALL DSCAL( N, REMAX, VL( 1, KI ), 1 )
*
               END IF
*
            ELSE
*
*              Complex left eigenvector.
*
*               Initial solve:
*                 ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
*                 ((T(KI+1,KI) T(KI+1,KI+1))                )
*
               IF( ABS( T( KI, KI+1 ) ).GE.ABS( T( KI+1, KI ) ) ) THEN
                  WORK( KI+N ) = WI / T( KI, KI+1 )
                  WORK( KI+1+N2 ) = ONE
               ELSE
                  WORK( KI+N ) = ONE
                  WORK( KI+1+N2 ) = -WI / T( KI+1, KI )
               END IF
               WORK( KI+1+N ) = ZERO
               WORK( KI+N2 ) = ZERO
*
*              Form right-hand side
*
               DO 190 K = KI + 2, N
                  WORK( K+N ) = -WORK( KI+N )*T( KI, K )
                  WORK( K+N2 ) = -WORK( KI+1+N2 )*T( KI+1, K )
  190          CONTINUE
*
*              Solve complex quasi-triangular system:
*              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
*
               VMAX = ONE
               VCRIT = BIGNUM
*
               JNXT = KI + 2
               DO 200 J = KI + 2, N
                  IF( J.LT.JNXT )
     $               GO TO 200
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J.LT.N ) THEN
                     IF( T( J+1, J ).NE.ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
*
                  IF( J1.EQ.J2 ) THEN
*
*                    1-by-1 diagonal block
*
*                    Scale if necessary to avoid overflow when
*                    forming the right-hand side elements.
*
                     IF( WORK( J ).GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
*
                     WORK( J+N ) = WORK( J+N ) -
     $                             DDOT( J-KI-2, T( KI+2, J ), 1,
     $                             WORK( KI+2+N ), 1 )
                     WORK( J+N2 ) = WORK( J+N2 ) -
     $                              DDOT( J-KI-2, T( KI+2, J ), 1,
     $                              WORK( KI+2+N2 ), 1 )
*
*                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
*
                     CALL DLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ),
     $                            LDT, ONE, ONE, WORK( J+N ), N, WR,
     $                            -WI, X, 2, SCALE, XNORM, IERR )
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     VMAX = MAX( ABS( WORK( J+N ) ),
     $                      ABS( WORK( J+N2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
*
                  ELSE
*
*                    2-by-2 diagonal block
*
*                    Scale if necessary to avoid overflow when forming
*                    the right-hand side elements.
*
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA.GT.VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
*
                     WORK( J+N ) = WORK( J+N ) -
     $                             DDOT( J-KI-2, T( KI+2, J ), 1,
     $                             WORK( KI+2+N ), 1 )
*
                     WORK( J+N2 ) = WORK( J+N2 ) -
     $                              DDOT( J-KI-2, T( KI+2, J ), 1,
     $                              WORK( KI+2+N2 ), 1 )
*
                     WORK( J+1+N ) = WORK( J+1+N ) -
     $                               DDOT( J-KI-2, T( KI+2, J+1 ), 1,
     $                               WORK( KI+2+N ), 1 )
*
                     WORK( J+1+N2 ) = WORK( J+1+N2 ) -
     $                                DDOT( J-KI-2, T( KI+2, J+1 ), 1,
     $                                WORK( KI+2+N2 ), 1 )
*
*                    Solve 2-by-2 complex linear equation
*                      ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
*                      ([T(j+1,j) T(j+1,j+1)]             )
*
                     CALL DLALN2( .TRUE., 2, 2, SMIN, ONE, T( J, J ),
     $                            LDT, ONE, ONE, WORK( J+N ), N, WR,
     $                            -WI, X, 2, SCALE, XNORM, IERR )
*
*                    Scale if necessary
*
                     IF( SCALE.NE.ONE ) THEN
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL DSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     WORK( J+1+N ) = X( 2, 1 )
                     WORK( J+1+N2 ) = X( 2, 2 )
                     VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ),
     $                      ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
*
                  END IF
  200          CONTINUE
*
*              Copy the vector x or Q*x to VL and normalize.
*
               IF( .NOT.OVER ) THEN
                  CALL DCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
                  CALL DCOPY( N-KI+1, WORK( KI+N2 ), 1, VL( KI, IS+1 ),
     $                        1 )
*
                  EMAX = ZERO
                  DO 220 K = KI, N
                     EMAX = MAX( EMAX, ABS( VL( K, IS ) )+
     $                      ABS( VL( K, IS+1 ) ) )
  220             CONTINUE
                  REMAX = ONE / EMAX
                  CALL DSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
                  CALL DSCAL( N-KI+1, REMAX, VL( KI, IS+1 ), 1 )
*
                  DO 230 K = 1, KI - 1
                     VL( K, IS ) = ZERO
                     VL( K, IS+1 ) = ZERO
  230             CONTINUE
               ELSE
                  IF( KI.LT.N-1 ) THEN
                     CALL DGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ),
     $                           LDVL, WORK( KI+2+N ), 1, WORK( KI+N ),
     $                           VL( 1, KI ), 1 )
                     CALL DGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ),
     $                           LDVL, WORK( KI+2+N2 ), 1,
     $                           WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  ELSE
                     CALL DSCAL( N, WORK( KI+N ), VL( 1, KI ), 1 )
                     CALL DSCAL( N, WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  END IF
*
                  EMAX = ZERO
                  DO 240 K = 1, N
                     EMAX = MAX( EMAX, ABS( VL( K, KI ) )+
     $                      ABS( VL( K, KI+1 ) ) )
  240             CONTINUE
                  REMAX = ONE / EMAX
                  CALL DSCAL( N, REMAX, VL( 1, KI ), 1 )
                  CALL DSCAL( N, REMAX, VL( 1, KI+1 ), 1 )
*
               END IF
*
            END IF
*
            IS = IS + 1
            IF( IP.NE.0 )
     $         IS = IS + 1
  250       CONTINUE
            IF( IP.EQ.-1 )
     $         IP = 0
            IF( IP.EQ.1 )
     $         IP = -1
*
  260    CONTINUE
*
      END IF
*
      RETURN
*
*     End of DTREVC
*
      END
      SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DTREXC reorders the real Schur factorization of a real matrix
*  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
*  moved to row ILST.
*
*  The real Schur form T is reordered by an orthogonal similarity
*  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
*  is updated by postmultiplying it with Z.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'V':  update the matrix Q of Schur vectors;
*          = 'N':  do not update Q.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          Schur canonical form.
*          On exit, the reordered upper quasi-triangular matrix, again
*          in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*          orthogonal transformation matrix Z which reorders T.
*          If COMPQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= max(1,N).
*
*  IFST    (input/output) INTEGER
*  ILST    (input/output) INTEGER
*          Specify the reordering of the diagonal blocks of T.
*          The block with row index IFST is moved to row ILST, by a
*          sequence of transpositions between adjacent blocks.
*          On exit, if IFST pointed on entry to the second row of a
*          2-by-2 block, it is changed to point to the first row; ILST
*          always points to the first row of the block in its final
*          position (which may differ from its input value by +1 or -1).
*          1 <= IFST <= N; 1 <= ILST <= N.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          = 1:  two adjacent blocks were too close to swap (the problem
*                is very ill-conditioned); T may have been partially
*                reordered, and ILST points to the first row of the
*                current position of the block being moved.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAEXC, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input arguments.
*
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -7
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREXC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Determine the first row of specified block
*     and find out it is 1 by 1 or 2 by 2.
*
      IF( IFST.GT.1 ) THEN
         IF( T( IFST, IFST-1 ).NE.ZERO )
     $      IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST.LT.N ) THEN
         IF( T( IFST+1, IFST ).NE.ZERO )
     $      NBF = 2
      END IF
*
*     Determine the first row of the final block
*     and find out it is 1 by 1 or 2 by 2.
*
      IF( ILST.GT.1 ) THEN
         IF( T( ILST, ILST-1 ).NE.ZERO )
     $      ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST.LT.N ) THEN
         IF( T( ILST+1, ILST ).NE.ZERO )
     $      NBL = 2
      END IF
*
      IF( IFST.EQ.ILST )
     $   RETURN
*
      IF( IFST.LT.ILST ) THEN
*
*        Update ILST
*
         IF( NBF.EQ.2 .AND. NBL.EQ.1 )
     $      ILST = ILST - 1
         IF( NBF.EQ.1 .AND. NBL.EQ.2 )
     $      ILST = ILST + 1
*
         HERE = IFST
*
   10    CONTINUE
*
*        Swap block with next one below
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1 by 1 or 2 by 2
*
            NBNEXT = 1
            IF( HERE+NBF+1.LE.N ) THEN
               IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT,
     $                   WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE + NBNEXT
*
*           Test if 2 by 2 block breaks into two 1 by 1 blocks
*
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1 by 1 blocks each of which
*           must be swapped individually
*
            NBNEXT = 1
            IF( HERE+3.LE.N ) THEN
               IF( T( HERE+3, HERE+2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT,
     $                   WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1 by 1 blocks, no problems possible
*
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT,
     $                      WORK, INFO )
               HERE = HERE + 1
            ELSE
*
*              Recompute NBNEXT in case 2 by 2 split
*
               IF( T( HERE+2, HERE+1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2 by 2 Block did not split
*
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1,
     $                         NBNEXT, WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 2
               ELSE
*
*                 2 by 2 Block did split
*
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1,
     $                         WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1,
     $                         WORK, INFO )
                  HERE = HERE + 2
               END IF
            END IF
         END IF
         IF( HERE.LT.ILST )
     $      GO TO 10
*
      ELSE
*
         HERE = IFST
   20    CONTINUE
*
*        Swap block with next one above
*
         IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
*
*           Current block either 1 by 1 or 2 by 2
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT,
     $                   NBF, WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE - NBNEXT
*
*           Test if 2 by 2 block breaks into two 1 by 1 blocks
*
            IF( NBF.EQ.2 ) THEN
               IF( T( HERE+1, HERE ).EQ.ZERO )
     $            NBF = 3
            END IF
*
         ELSE
*
*           Current block consists of two 1 by 1 blocks each of which
*           must be swapped individually
*
            NBNEXT = 1
            IF( HERE.GE.3 ) THEN
               IF( T( HERE-1, HERE-2 ).NE.ZERO )
     $            NBNEXT = 2
            END IF
            CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT,
     $                   1, WORK, INFO )
            IF( INFO.NE.0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT.EQ.1 ) THEN
*
*              Swap two 1 by 1 blocks, no problems possible
*
               CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1,
     $                      WORK, INFO )
               HERE = HERE - 1
            ELSE
*
*              Recompute NBNEXT in case 2 by 2 split
*
               IF( T( HERE, HERE-1 ).EQ.ZERO )
     $            NBNEXT = 1
               IF( NBNEXT.EQ.2 ) THEN
*
*                 2 by 2 Block did not split
*
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1,
     $                         WORK, INFO )
                  IF( INFO.NE.0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 2
               ELSE
*
*                 2 by 2 Block did split
*
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1,
     $                         WORK, INFO )
                  CALL DLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1,
     $                         WORK, INFO )
                  HERE = HERE - 2
               END IF
            END IF
         END IF
         IF( HERE.GT.ILST )
     $      GO TO 20
      END IF
      ILST = HERE
*
      RETURN
*
*     End of DTREXC
*
      END
      SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI,
     $                   M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
      DOUBLE PRECISION   S, SEP
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ),
     $                   WR( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSEN reorders the real Schur factorization of a real matrix
*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
*  the leading diagonal blocks of the upper quasi-triangular matrix T,
*  and the leading columns of Q form an orthonormal basis of the
*  corresponding right invariant subspace.
*
*  Optionally the routine computes the reciprocal condition numbers of
*  the cluster of eigenvalues and/or the invariant subspace.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elemnts equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for the
*          cluster of eigenvalues (S) or the invariant subspace (SEP):
*          = 'N': none;
*          = 'E': for eigenvalues only (S);
*          = 'V': for invariant subspace only (SEP);
*          = 'B': for both eigenvalues and invariant subspace (S and
*                 SEP).
*
*  COMPQ   (input) CHARACTER*1
*          = 'V': update the matrix Q of Schur vectors;
*          = 'N': do not update Q.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster. To
*          select a real eigenvalue w(j), SELECT(j) must be set to
*          .TRUE.. To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*          On exit, T is overwritten by the reordered matrix T, again in
*          Schur canonical form, with the selected eigenvalues in the
*          leading diagonal blocks.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*          orthogonal transformation matrix which reorders T; the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If COMPQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          The real and imaginary parts, respectively, of the reordered
*          eigenvalues of T. The eigenvalues are stored in the same
*          order as on the diagonal of T, with WR(i) = T(i,i) and, if
*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
*          sufficiently ill-conditioned, then its value may differ
*          significantly from its value before reordering.
*
*  M       (output) INTEGER
*          The dimension of the specified invariant subspace.
*          0 < = M <= N.
*
*  S       (output) DOUBLE PRECISION
*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
*          condition number for the selected cluster of eigenvalues.
*          S cannot underestimate the true reciprocal condition number
*          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
*          If JOB = 'N' or 'V', S is not referenced.
*
*  SEP     (output) DOUBLE PRECISION
*          If JOB = 'V' or 'B', SEP is the estimated reciprocal
*          condition number of the specified invariant subspace. If
*          M = 0 or N, SEP = norm(T).
*          If JOB = 'N' or 'E', SEP is not referenced.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If JOB = 'N', LWORK >= max(1,N);
*          if JOB = 'E', LWORK >= max(1,M*(N-M));
*          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If JOB = 'N' or 'E', LIWORK >= 1;
*          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: reordering of T failed because some eigenvalues are too
*               close to separate (the problem is very ill-conditioned);
*               T may have been partially reordered, and WR and WI
*               contain the eigenvalues in the same order as in T; S and
*               SEP (if requested) are set to zero.
*
*  Further Details
*  ===============
*
*  DTRSEN first collects the selected eigenvalues by computing an
*  orthogonal transformation Z to move them to the top left corner of T.
*  In other words, the selected eigenvalues are the eigenvalues of T11
*  in:
*
*                Z'*T*Z = ( T11 T12 ) n1
*                         (  0  T22 ) n2
*                            n1  n2
*
*  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns
*  of Z span the specified invariant subspace of T.
*
*  If T has been obtained from the real Schur factorization of a matrix
*  A = Q*T*Q', then the reordered real Schur factorization of A is given
*  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span
*  the corresponding invariant subspace of A.
*
*  The reciprocal condition number of the average of the eigenvalues of
*  T11 may be returned in S. S lies between 0 (very badly conditioned)
*  and 1 (very well conditioned). It is computed as follows. First we
*  compute R so that
*
*                         P = ( I  R ) n1
*                             ( 0  0 ) n2
*                               n1 n2
*
*  is the projector on the invariant subspace associated with T11.
*  R is the solution of the Sylvester equation:
*
*                        T11*R - R*T22 = T12.
*
*  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
*  the two-norm of M. Then S is computed as the lower bound
*
*                      (1 + F-norm(R)**2)**(-1/2)
*
*  on the reciprocal of 2-norm(P), the true reciprocal condition number.
*  S cannot underestimate 1 / 2-norm(P) by more than a factor of
*  sqrt(N).
*
*  An approximate error bound for the computed average of the
*  eigenvalues of T11 is
*
*                         EPS * norm(T) / S
*
*  where EPS is the machine precision.
*
*  The reciprocal condition number of the right invariant subspace
*  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
*  SEP is defined as the separation of T11 and T22:
*
*                     sep( T11, T22 ) = sigma-min( C )
*
*  where sigma-min(C) is the smallest singular value of the
*  n1*n2-by-n1*n2 matrix
*
*     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
*
*  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
*  product. We estimate sigma-min(C) by the reciprocal of an estimate of
*  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
*  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
*
*  When SEP is small, small changes in T can cause large changes in
*  the invariant subspace. An approximate bound on the maximum angular
*  error in the computed right invariant subspace is
*
*                      EPS * norm(T) / SEP
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, PAIR, SWAP, WANTBH, WANTQ, WANTS,
     $                   WANTSP
      INTEGER            IERR, K, KASE, KK, KS, LIWMIN, LWMIN, N1, N2,
     $                   NN
      DOUBLE PRECISION   EST, RNORM, SCALE
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLANGE
      EXTERNAL           LSAME, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DLACPY, DTREXC, DTRSYL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP )
     $     THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -8
      ELSE
*
*        Set M to the dimension of the specified invariant subspace,
*        and test LWORK and LIWORK.
*
         M = 0
         PAIR = .FALSE.
         DO 10 K = 1, N
            IF( PAIR ) THEN
               PAIR = .FALSE.
            ELSE
               IF( K.LT.N ) THEN
                  IF( T( K+1, K ).EQ.ZERO ) THEN
                     IF( SELECT( K ) )
     $                  M = M + 1
                  ELSE
                     PAIR = .TRUE.
                     IF( SELECT( K ) .OR. SELECT( K+1 ) )
     $                  M = M + 2
                  END IF
               ELSE
                  IF( SELECT( N ) )
     $               M = M + 1
               END IF
            END IF
   10    CONTINUE
*
         N1 = M
         N2 = N - M
         NN = N1*N2
*
         IF( WANTSP ) THEN
            LWMIN = MAX( 1, 2*NN )
            LIWMIN = MAX( 1, NN )
         ELSE IF( LSAME( JOB, 'N' ) ) THEN
            LWMIN = MAX( 1, N )
            LIWMIN = 1
         ELSE IF( LSAME( JOB, 'E' ) ) THEN
            LWMIN = MAX( 1, NN )
            LIWMIN = 1
         END IF
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -15
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -17
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTS )
     $      S = ONE
         IF( WANTSP )
     $      SEP = DLANGE( '1', N, N, T, LDT, WORK )
         GO TO 40
      END IF
*
*     Collect the selected blocks at the top-left corner of T.
*
      KS = 0
      PAIR = .FALSE.
      DO 20 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         ELSE
            SWAP = SELECT( K )
            IF( K.LT.N ) THEN
               IF( T( K+1, K ).NE.ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT( K+1 )
               END IF
            END IF
            IF( SWAP ) THEN
               KS = KS + 1
*
*              Swap the K-th block to position KS.
*
               IERR = 0
               KK = K
               IF( K.NE.KS )
     $            CALL DTREXC( COMPQ, N, T, LDT, Q, LDQ, KK, KS, WORK,
     $                         IERR )
               IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
*
*                 Blocks too close to swap: exit.
*
                  INFO = 1
                  IF( WANTS )
     $               S = ZERO
                  IF( WANTSP )
     $               SEP = ZERO
                  GO TO 40
               END IF
               IF( PAIR )
     $            KS = KS + 1
            END IF
         END IF
   20 CONTINUE
*
      IF( WANTS ) THEN
*
*        Solve Sylvester equation for R:
*
*           T11*R - R*T22 = scale*T12
*
         CALL DLACPY( 'F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 )
         CALL DTRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ),
     $                LDT, WORK, N1, SCALE, IERR )
*
*        Estimate the reciprocal of the condition number of the cluster
*        of eigenvalues.
*
         RNORM = DLANGE( 'F', N1, N2, WORK, N1, WORK )
         IF( RNORM.EQ.ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )*
     $          SQRT( RNORM ) )
         END IF
      END IF
*
      IF( WANTSP ) THEN
*
*        Estimate sep(T11,T22).
*
         EST = ZERO
         KASE = 0
   30    CONTINUE
         CALL DLACN2( NN, WORK( NN+1 ), WORK, IWORK, EST, KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Solve  T11*R - R*T22 = scale*X.
*
               CALL DTRSYL( 'N', 'N', -1, N1, N2, T, LDT,
     $                      T( N1+1, N1+1 ), LDT, WORK, N1, SCALE,
     $                      IERR )
            ELSE
*
*              Solve  T11'*R - R*T22' = scale*X.
*
               CALL DTRSYL( 'T', 'T', -1, N1, N2, T, LDT,
     $                      T( N1+1, N1+1 ), LDT, WORK, N1, SCALE,
     $                      IERR )
            END IF
            GO TO 30
         END IF
*
         SEP = SCALE / EST
      END IF
*
   40 CONTINUE
*
*     Store the output eigenvalues in WR and WI.
*
      DO 50 K = 1, N
         WR( K ) = T( K, K )
         WI( K ) = ZERO
   50 CONTINUE
      DO 60 K = 1, N - 1
         IF( T( K+1, K ).NE.ZERO ) THEN
            WI( K ) = SQRT( ABS( T( K, K+1 ) ) )*
     $                SQRT( ABS( T( K+1, K ) ) )
            WI( K+1 ) = -WI( K )
         END IF
   60 CONTINUE
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of DTRSEN
*
      END
      SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     Modified to call DLACN2 in place of DLACON, 5 Feb 03, SJH.
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, JOB
      INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSNA estimates reciprocal condition numbers for specified
*  eigenvalues and/or right eigenvectors of a real upper
*  quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
*  orthogonal).
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for
*          eigenvalues (S) or eigenvectors (SEP):
*          = 'E': for eigenvalues only (S);
*          = 'V': for eigenvectors only (SEP);
*          = 'B': for both eigenvalues and eigenvectors (S and SEP).
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute condition numbers for all eigenpairs;
*          = 'S': compute condition numbers for selected eigenpairs
*                 specified by the array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
*          condition numbers are required. To select condition numbers
*          for the eigenpair corresponding to a real eigenvalue w(j),
*          SELECT(j) must be set to .TRUE.. To select condition numbers
*          corresponding to a complex conjugate pair of eigenvalues w(j)
*          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
*          set to .TRUE..
*          If HOWMNY = 'A', SELECT is not referenced.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          The upper quasi-triangular matrix T, in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M)
*          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
*          must be stored in consecutive columns of VL, as returned by
*          DHSEIN or DTREVC.
*          If JOB = 'V', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.
*          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
*
*  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M)
*          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
*          must be stored in consecutive columns of VR, as returned by
*          DHSEIN or DTREVC.
*          If JOB = 'V', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
*
*  S       (output) DOUBLE PRECISION array, dimension (MM)
*          If JOB = 'E' or 'B', the reciprocal condition numbers of the
*          selected eigenvalues, stored in consecutive elements of the
*          array. For a complex conjugate pair of eigenvalues two
*          consecutive elements of S are set to the same value. Thus
*          S(j), SEP(j), and the j-th columns of VL and VR all
*          correspond to the same eigenpair (but not in general the
*          j-th eigenpair, unless all eigenpairs are selected).
*          If JOB = 'V', S is not referenced.
*
*  SEP     (output) DOUBLE PRECISION array, dimension (MM)
*          If JOB = 'V' or 'B', the estimated reciprocal condition
*          numbers of the selected eigenvectors, stored in consecutive
*          elements of the array. For a complex eigenvector two
*          consecutive elements of SEP are set to the same value. If
*          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)
*          is set to 0; this can only occur when the true value would be
*          very small anyway.
*          If JOB = 'E', SEP is not referenced.
*
*  MM      (input) INTEGER
*          The number of elements in the arrays S (if JOB = 'E' or 'B')
*           and/or SEP (if JOB = 'V' or 'B'). MM >= M.
*
*  M       (output) INTEGER
*          The number of elements of the arrays S and/or SEP actually
*          used to store the estimated condition numbers.
*          If HOWMNY = 'A', M is set to N.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,N+6)
*          If JOB = 'E', WORK is not referenced.
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
*
*  IWORK   (workspace) INTEGER array, dimension (2*(N-1))
*          If JOB = 'E', IWORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The reciprocal of the condition number of an eigenvalue lambda is
*  defined as
*
*          S(lambda) = |v'*u| / (norm(u)*norm(v))
*
*  where u and v are the right and left eigenvectors of T corresponding
*  to lambda; v' denotes the conjugate-transpose of v, and norm(u)
*  denotes the Euclidean norm. These reciprocal condition numbers always
*  lie between zero (very badly conditioned) and one (very well
*  conditioned). If n = 1, S(lambda) is defined to be 1.
*
*  An approximate error bound for a computed eigenvalue W(i) is given by
*
*                      EPS * norm(T) / S(i)
*
*  where EPS is the machine precision.
*
*  The reciprocal of the condition number of the right eigenvector u
*  corresponding to lambda is defined as follows. Suppose
*
*              T = ( lambda  c  )
*                  (   0    T22 )
*
*  Then the reciprocal condition number is
*
*          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
*
*  where sigma-min denotes the smallest singular value. We approximate
*  the smallest singular value by the reciprocal of an estimate of the
*  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
*  defined to be abs(T(1,1)).
*
*  An approximate error bound for a computed right eigenvector VR(i)
*  is given by
*
*                      EPS * norm(T) / SEP(i)
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            PAIR, SOMCON, WANTBH, WANTS, WANTSP
      INTEGER            I, IERR, IFST, ILST, J, K, KASE, KS, N2, NN
      DOUBLE PRECISION   BIGNUM, COND, CS, DELTA, DUMM, EPS, EST, LNRM,
     $                   MU, PROD, PROD1, PROD2, RNRM, SCALE, SMLNUM, SN
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
      DOUBLE PRECISION   DUMMY( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLAPY2, DNRM2
      EXTERNAL           LSAME, DDOT, DLAMCH, DLAPY2, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACN2, DLACPY, DLAQTR, DTREXC, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
*
      SOMCON = LSAME( HOWMNY, 'S' )
*
      INFO = 0
      IF( .NOT.WANTS .AND. .NOT.WANTSP ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' ) .AND. .NOT.SOMCON ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( WANTS .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTS .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE
*
*        Set M to the number of eigenpairs for which condition numbers
*        are required, and test MM.
*
         IF( SOMCON ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 K = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
               ELSE
                  IF( K.LT.N ) THEN
                     IF( T( K+1, K ).EQ.ZERO ) THEN
                        IF( SELECT( K ) )
     $                     M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( K ) .OR. SELECT( K+1 ) )
     $                     M = M + 2
                     END IF
                  ELSE
                     IF( SELECT( N ) )
     $                  M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
*
         IF( MM.LT.M ) THEN
            INFO = -13
         ELSE IF( LDWORK.LT.1 .OR. ( WANTSP .AND. LDWORK.LT.N ) ) THEN
            INFO = -16
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSNA', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( SOMCON ) THEN
            IF( .NOT.SELECT( 1 ) )
     $         RETURN
         END IF
         IF( WANTS )
     $      S( 1 ) = ONE
         IF( WANTSP )
     $      SEP( 1 ) = ABS( T( 1, 1 ) )
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
      KS = 0
      PAIR = .FALSE.
      DO 60 K = 1, N
*
*        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
*
         IF( PAIR ) THEN
            PAIR = .FALSE.
            GO TO 60
         ELSE
            IF( K.LT.N )
     $         PAIR = T( K+1, K ).NE.ZERO
         END IF
*
*        Determine whether condition numbers are required for the k-th
*        eigenpair.
*
         IF( SOMCON ) THEN
            IF( PAIR ) THEN
               IF( .NOT.SELECT( K ) .AND. .NOT.SELECT( K+1 ) )
     $            GO TO 60
            ELSE
               IF( .NOT.SELECT( K ) )
     $            GO TO 60
            END IF
         END IF
*
         KS = KS + 1
*
         IF( WANTS ) THEN
*
*           Compute the reciprocal condition number of the k-th
*           eigenvalue.
*
            IF( .NOT.PAIR ) THEN
*
*              Real eigenvalue.
*
               PROD = DDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               RNRM = DNRM2( N, VR( 1, KS ), 1 )
               LNRM = DNRM2( N, VL( 1, KS ), 1 )
               S( KS ) = ABS( PROD ) / ( RNRM*LNRM )
            ELSE
*
*              Complex eigenvalue.
*
               PROD1 = DDOT( N, VR( 1, KS ), 1, VL( 1, KS ), 1 )
               PROD1 = PROD1 + DDOT( N, VR( 1, KS+1 ), 1, VL( 1, KS+1 ),
     $                 1 )
               PROD2 = DDOT( N, VL( 1, KS ), 1, VR( 1, KS+1 ), 1 )
               PROD2 = PROD2 - DDOT( N, VL( 1, KS+1 ), 1, VR( 1, KS ),
     $                 1 )
               RNRM = DLAPY2( DNRM2( N, VR( 1, KS ), 1 ),
     $                DNRM2( N, VR( 1, KS+1 ), 1 ) )
               LNRM = DLAPY2( DNRM2( N, VL( 1, KS ), 1 ),
     $                DNRM2( N, VL( 1, KS+1 ), 1 ) )
               COND = DLAPY2( PROD1, PROD2 ) / ( RNRM*LNRM )
               S( KS ) = COND
               S( KS+1 ) = COND
            END IF
         END IF
*
         IF( WANTSP ) THEN
*
*           Estimate the reciprocal condition number of the k-th
*           eigenvector.
*
*           Copy the matrix T to the array WORK and swap the diagonal
*           block beginning at T(k,k) to the (1,1) position.
*
            CALL DLACPY( 'Full', N, N, T, LDT, WORK, LDWORK )
            IFST = K
            ILST = 1
            CALL DTREXC( 'No Q', N, WORK, LDWORK, DUMMY, 1, IFST, ILST,
     $                   WORK( 1, N+1 ), IERR )
*
            IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
*
*              Could not swap because blocks not well separated
*
               SCALE = ONE
               EST = BIGNUM
            ELSE
*
*              Reordering successful
*
               IF( WORK( 2, 1 ).EQ.ZERO ) THEN
*
*                 Form C = T22 - lambda*I in WORK(2:N,2:N).
*
                  DO 20 I = 2, N
                     WORK( I, I ) = WORK( I, I ) - WORK( 1, 1 )
   20             CONTINUE
                  N2 = 1
                  NN = N - 1
               ELSE
*
*                 Triangularize the 2 by 2 block by unitary
*                 transformation U = [  cs   i*ss ]
*                                    [ i*ss   cs  ].
*                 such that the (1,1) position of WORK is complex
*                 eigenvalue lambda with positive imaginary part. (2,2)
*                 position of WORK is the complex eigenvalue lambda
*                 with negative imaginary  part.
*
                  MU = SQRT( ABS( WORK( 1, 2 ) ) )*
     $                 SQRT( ABS( WORK( 2, 1 ) ) )
                  DELTA = DLAPY2( MU, WORK( 2, 1 ) )
                  CS = MU / DELTA
                  SN = -WORK( 2, 1 ) / DELTA
*
*                 Form
*
*                 C' = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
*                                        [   mu                     ]
*                                        [         ..               ]
*                                        [             ..           ]
*                                        [                  mu      ]
*                 where C' is conjugate transpose of complex matrix C,
*                 and RWORK is stored starting in the N+1-st column of
*                 WORK.
*
                  DO 30 J = 3, N
                     WORK( 2, J ) = CS*WORK( 2, J )
                     WORK( J, J ) = WORK( J, J ) - WORK( 1, 1 )
   30             CONTINUE
                  WORK( 2, 2 ) = ZERO
*
                  WORK( 1, N+1 ) = TWO*MU
                  DO 40 I = 2, N - 1
                     WORK( I, N+1 ) = SN*WORK( 1, I+1 )
   40             CONTINUE
                  N2 = 2
                  NN = 2*( N-1 )
               END IF
*
*              Estimate norm(inv(C'))
*
               EST = ZERO
               KASE = 0
   50          CONTINUE
               CALL DLACN2( NN, WORK( 1, N+2 ), WORK( 1, N+4 ), IWORK,
     $                      EST, KASE, ISAVE )
               IF( KASE.NE.0 ) THEN
                  IF( KASE.EQ.1 ) THEN
                     IF( N2.EQ.1 ) THEN
*
*                       Real eigenvalue: solve C'*x = scale*c.
*
                        CALL DLAQTR( .TRUE., .TRUE., N-1, WORK( 2, 2 ),
     $                               LDWORK, DUMMY, DUMM, SCALE,
     $                               WORK( 1, N+4 ), WORK( 1, N+6 ),
     $                               IERR )
                     ELSE
*
*                       Complex eigenvalue: solve
*                       C'*(p+iq) = scale*(c+id) in real arithmetic.
*
                        CALL DLAQTR( .TRUE., .FALSE., N-1, WORK( 2, 2 ),
     $                               LDWORK, WORK( 1, N+1 ), MU, SCALE,
     $                               WORK( 1, N+4 ), WORK( 1, N+6 ),
     $                               IERR )
                     END IF
                  ELSE
                     IF( N2.EQ.1 ) THEN
*
*                       Real eigenvalue: solve C*x = scale*c.
*
                        CALL DLAQTR( .FALSE., .TRUE., N-1, WORK( 2, 2 ),
     $                               LDWORK, DUMMY, DUMM, SCALE,
     $                               WORK( 1, N+4 ), WORK( 1, N+6 ),
     $                               IERR )
                     ELSE
*
*                       Complex eigenvalue: solve
*                       C*(p+iq) = scale*(c+id) in real arithmetic.
*
                        CALL DLAQTR( .FALSE., .FALSE., N-1,
     $                               WORK( 2, 2 ), LDWORK,
     $                               WORK( 1, N+1 ), MU, SCALE,
     $                               WORK( 1, N+4 ), WORK( 1, N+6 ),
     $                               IERR )
*
                     END IF
                  END IF
*
                  GO TO 50
               END IF
            END IF
*
            SEP( KS ) = SCALE / MAX( EST, SMLNUM )
            IF( PAIR )
     $         SEP( KS+1 ) = SEP( KS )
         END IF
*
         IF( PAIR )
     $      KS = KS + 1
*
   60 CONTINUE
      RETURN
*
*     End of DTRSNA
*
      END
      SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $                   LDC, SCALE, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSYL solves the real Sylvester matrix equation:
*
*     op(A)*X + X*op(B) = scale*C or
*     op(A)*X - X*op(B) = scale*C,
*
*  where op(A) = A or A**T, and  A and B are both upper quasi-
*  triangular. A is M-by-M and B is N-by-N; the right hand side C and
*  the solution X are M-by-N; and scale is an output scale factor, set
*  <= 1 to avoid overflow in X.
*
*  A and B must be in Schur canonical form (as returned by DHSEQR), that
*  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
*  each 2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  TRANA   (input) CHARACTER*1
*          Specifies the option op(A):
*          = 'N': op(A) = A    (No transpose)
*          = 'T': op(A) = A**T (Transpose)
*          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
*
*  TRANB   (input) CHARACTER*1
*          Specifies the option op(B):
*          = 'N': op(B) = B    (No transpose)
*          = 'T': op(B) = B**T (Transpose)
*          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
*
*  ISGN    (input) INTEGER
*          Specifies the sign in the equation:
*          = +1: solve op(A)*X + X*op(B) = scale*C
*          = -1: solve op(A)*X - X*op(B) = scale*C
*
*  M       (input) INTEGER
*          The order of the matrix A, and the number of rows in the
*          matrices X and C. M >= 0.
*
*  N       (input) INTEGER
*          The order of the matrix B, and the number of columns in the
*          matrices X and C. N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
*          The upper quasi-triangular matrix A, in Schur canonical form.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The upper quasi-triangular matrix B, in Schur canonical form.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N right hand side matrix C.
*          On exit, C is overwritten by the solution matrix X.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M)
*
*  SCALE   (output) DOUBLE PRECISION
*          The scale factor, scale, set <= 1 to avoid overflow in X.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: A and B have common or very close eigenvalues; perturbed
*               values were used to solve the equation (but the matrices
*               A and B are unchanged).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT
      DOUBLE PRECISION   A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,
     $                   SMLNUM, SUML, SUMR, XNORM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT, DLAMCH, DLANGE
      EXTERNAL           LSAME, DDOT, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, DLALN2, DLASY2, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and Test input parameters
*
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
*
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT.
     $    LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT.
     $         LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN.NE.1 .AND. ISGN.NE.-1 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRSYL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Set constants to control overflow
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
*
      SMIN = MAX( SMLNUM, EPS*DLANGE( 'M', M, M, A, LDA, DUM ),
     $       EPS*DLANGE( 'M', N, N, B, LDB, DUM ) )
*
      SGN = ISGN
*
      IF( NOTRNA .AND. NOTRNB ) THEN
*
*        Solve    A*X + ISGN*X*B = scale*C.
*
*        The (K,L)th block of X is determined starting from
*        bottom-left corner column by column by
*
*         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
*
*        Where
*                  M                         L-1
*        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
*                I=K+1                       J=1
*
*        Start column loop (index = L)
*        L1 (L2) : column index of the first (first) row of X(K,L).
*
         LNEXT = 1
         DO 60 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 60
            IF( L.EQ.N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L ).NE.ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
*
*           Start row loop (index = K)
*           K1 (K2): row index of the first (last) row of X(K,L).
*
            KNEXT = M
            DO 50 K = M, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 50
               IF( K.EQ.1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 ).NE.ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
*
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 20 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 30 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   30                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*
                  CALL DLASY2( .FALSE., .FALSE., ISGN, 2, 2,
     $                         A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC,
     $                         2, SCALOC, X, 2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 40 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*
   50       CONTINUE
*
   60    CONTINUE
*
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
*
*        Solve    A' *X + ISGN*X*B = scale*C.
*
*        The (K,L)th block of X is determined starting from
*        upper-left corner column by column by
*
*          A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
*
*        Where
*                   K-1                        L-1
*          R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
*                   I=1                        J=1
*
*        Start column loop (index = L)
*        L1 (L2): column index of the first (last) row of X(K,L)
*
         LNEXT = 1
         DO 120 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 120
            IF( L.EQ.N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L ).NE.ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
*
*           Start row loop (index = K)
*           K1 (K2): row index of the first (last) row of X(K,L)
*
            KNEXT = 1
            DO 110 K = 1, M
               IF( K.LT.KNEXT )
     $            GO TO 110
               IF( K.EQ.M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K ).NE.ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
*
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 70 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   70                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 90 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*
                  CALL DLASY2( .TRUE., .FALSE., ISGN, 2, 2, A( K1, K1 ),
     $                         LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                         2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 100 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*
  110       CONTINUE
  120    CONTINUE
*
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
*
*        Solve    A'*X + ISGN*X*B' = scale*C.
*
*        The (K,L)th block of X is determined starting from
*        top-right corner column by column by
*
*           A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
*
*        Where
*                     K-1                          N
*            R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
*                     I=1                        J=L+1
*
*        Start column loop (index = L)
*        L1 (L2): column index of the first (last) row of X(K,L)
*
         LNEXT = N
         DO 180 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 180
            IF( L.EQ.1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 ).NE.ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
*
*           Start row loop (index = K)
*           K1 (K2): row index of the first (last) row of X(K,L)
*
            KNEXT = 1
            DO 170 K = 1, M
               IF( K.LT.KNEXT )
     $            GO TO 170
               IF( K.EQ.M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K ).NE.ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
*
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC,
     $                   B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 130 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  130                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  CALL DLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 140 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  140                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 150 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  150                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                   B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*
                  CALL DLASY2( .TRUE., .TRUE., ISGN, 2, 2, A( K1, K1 ),
     $                         LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                         2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 160 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  160                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*
  170       CONTINUE
  180    CONTINUE
*
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
*
*        Solve    A*X + ISGN*X*B' = scale*C.
*
*        The (K,L)th block of X is determined starting from
*        bottom-right corner column by column by
*
*            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
*
*        Where
*                      M                          N
*            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
*                    I=K+1                      J=L+1
*
*        Start column loop (index = L)
*        L1 (L2): column index of the first (last) row of X(K,L)
*
         LNEXT = N
         DO 240 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 240
            IF( L.EQ.1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 ).NE.ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
*
*           Start row loop (index = K)
*           K1 (K2): row index of the first (last) row of X(K,L)
*
            KNEXT = M
            DO 230 K = M, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 230
               IF( K.EQ.1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 ).NE.ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
*
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC,
     $                   B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
*
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 190 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  190                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
*
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 200 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  200                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
*
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
*
                  SUML = DDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
*
                  CALL DLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 210 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  210                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
*
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
*
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                   B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                   B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
*
                  SUML = DDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                   C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = DDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                   B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
*
                  CALL DLASY2( .FALSE., .TRUE., ISGN, 2, 2, A( K1, K1 ),
     $                         LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                         2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
*
                  IF( SCALOC.NE.ONE ) THEN
                     DO 220 J = 1, N
                        CALL DSCAL( M, SCALOC, C( 1, J ), 1 )
  220                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
*
  230       CONTINUE
  240    CONTINUE
*
      END IF
*
      RETURN
*
*     End of DTRSYL
*
      END
*-----------------------------------------------------------------------
*  Routine:    DVOUT
*
*  Purpose:    Real vector output routine.
*
*  Usage:      CALL DVOUT (LOUT, N, SX, IDIGIT, IFMT)
*
*  Arguments
*     N      - Length of array SX.  (Input)
*     SX     - Real array to be printed.  (Input)
*     IFMT   - Format to be used in printing array SX.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*-----------------------------------------------------------------------
*
      SUBROUTINE DVOUT( LOUT, N, SX, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
*     .. Scalar Arguments ..
      CHARACTER*( * )    IFMT
      INTEGER            IDIGIT, LOUT, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
*     ..
*     .. Local Scalars ..
      CHARACTER*80       LINE
      INTEGER            I, K1, K2, LLL, NDIGIT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MIN, MIN0
*     ..
*     .. Executable Statements ..
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
*
      LLL = MIN( LEN( IFMT ), 80 )
      DO 10 I = 1, LLL
         LINE( I: I ) = '-'
   10 CONTINUE
*
      DO 20 I = LLL + 1, 80
         LINE( I: I ) = ' '
   20 CONTINUE
*
      WRITE( LOUT, FMT = 9999 )IFMT, LINE( 1: LLL )
 9999 FORMAT( / 1X, A, / 1X, A )
*
      IF( N.LE.0 )
     $   RETURN
      NDIGIT = IDIGIT
      IF( IDIGIT.EQ.0 )
     $   NDIGIT = 4
*
*=======================================================================
*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
*=======================================================================
*
      IF( IDIGIT.LT.0 ) THEN
         NDIGIT = -IDIGIT
         IF( NDIGIT.LE.4 ) THEN
            DO 30 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, FMT = 9998 )K1, K2, ( SX( I ), I = K1, K2 )
   30       CONTINUE
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 40 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               WRITE( LOUT, FMT = 9997 )K1, K2, ( SX( I ), I = K1, K2 )
   40       CONTINUE
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 50 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               WRITE( LOUT, FMT = 9996 )K1, K2, ( SX( I ), I = K1, K2 )
   50       CONTINUE
         ELSE
            DO 60 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, FMT = 9995 )K1, K2, ( SX( I ), I = K1, K2 )
   60       CONTINUE
         END IF
*
*=======================================================================
*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
*=======================================================================
*
      ELSE
         IF( NDIGIT.LE.4 ) THEN
            DO 70 K1 = 1, N, 10
               K2 = MIN0( N, K1+9 )
               WRITE( LOUT, FMT = 9998 )K1, K2, ( SX( I ), I = K1, K2 )
   70       CONTINUE
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 80 K1 = 1, N, 8
               K2 = MIN0( N, K1+7 )
               WRITE( LOUT, FMT = 9997 )K1, K2, ( SX( I ), I = K1, K2 )
   80       CONTINUE
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 90 K1 = 1, N, 6
               K2 = MIN0( N, K1+5 )
               WRITE( LOUT, FMT = 9996 )K1, K2, ( SX( I ), I = K1, K2 )
   90       CONTINUE
         ELSE
            DO 100 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, FMT = 9995 )K1, K2, ( SX( I ), I = K1, K2 )
  100       CONTINUE
         END IF
      END IF
      WRITE( LOUT, FMT = 9994 )
      RETURN
 9998 FORMAT( 1X, I4, ' - ', I4, ':', 1P, 10D12.3 )
 9997 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P, 8D14.5 )
 9996 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P, 6D18.9 )
 9995 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P, 5D24.13 )
 9994 FORMAT( 1X, ' ' )
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      INTEGER FUNCTION ILADLC(M, N, A, LDA)
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ILADLC scans A for its last non-zero column.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( N.EQ.0 ) THEN
         ILADLC = N
      ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLC = N
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ILADLC = N, 1, -1
            DO I = 1, M
               IF( A(I, ILADLC).NE.ZERO ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END FUNCTION
      INTEGER FUNCTION ILADLR(M, N, A, LDA)
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ILADLR scans A for its last non-zero row.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER I, J
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILADLR = M
      ELSE IF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILADLR = M
      ELSE
*     Scan up each column tracking the last zero row seen.
         ILADLR = 0
         DO J = 1, N
            DO I = M, 1, -1
               IF( A(I, J).NE.ZERO ) EXIT
            END DO
            ILADLR = MAX( ILADLR, I )
         END DO
      END IF
      RETURN
      END FUNCTION
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.2.1)                        --
*
*  -- April 2009                                                      --
*
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  ILAENV returns an INTEGER
*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (DEPRECATED)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR method
*               for nonsymmetric eigenvalue problems (DEPRECATED)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ISPEC <= 16:
*               xHSEQR or one of its subroutines,
*               see IPARMQ for detailed explanation
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC+64 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1: 1 ) = CHAR( IC-32 )
            DO 40 I = 2, 6
               IC = ICHAR( SUBNAM( I: I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1: 1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. 
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*     
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  Purpose
*  =======
*
*       This program sets problem and machine dependent parameters
*       useful for xHSEQR and its subroutines. It is called whenever 
*       ILAENV is called with 12 <= ISPEC <= 16
*
*  Arguments
*  =========
*
*       ISPEC  (input) integer scalar
*              ISPEC specifies which tunable parameter IPARMQ should
*              return.
*
*              ISPEC=12: (INMIN)  Matrices of order nmin or less
*                        are sent directly to xLAHQR, the implicit
*                        double shift QR algorithm.  NMIN must be
*                        at least 11.
*
*              ISPEC=13: (INWIN)  Size of the deflation window.
*                        This is best set greater than or equal to
*                        the number of simultaneous shifts NS.
*                        Larger matrices benefit from larger deflation
*                        windows.
*
*              ISPEC=14: (INIBL) Determines when to stop nibbling and
*                        invest in an (expensive) multi-shift QR sweep.
*                        If the aggressive early deflation subroutine
*                        finds LD converged eigenvalues from an order
*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*                        then the next QR sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  Setting
*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*                        multi-shift QR sweep whenever early deflation
*                        finds a converged eigenvalue.  Setting
*                        IPARMQ(ISPEC=14) greater than or equal to 100
*                        prevents TTQRE from skipping a multi-shift
*                        QR sweep.
*
*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*                        a multi-shift QR iteration.
*
*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*                        following meanings.
*                        0:  During the multi-shift QR sweep,
*                            xLAQR5 does not accumulate reflections and
*                            does not use matrix-matrix multiply to
*                            update the far-from-diagonal matrix
*                            entries.
*                        1:  During the multi-shift QR sweep,
*                            xLAQR5 and/or xLAQRaccumulates reflections and uses
*                            matrix-matrix multiply to update the
*                            far-from-diagonal matrix entries.
*                        2:  During the multi-shift QR sweep.
*                            xLAQR5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*                        (If xTRMM is slower than xGEMM, then
*                        IPARMQ(ISPEC=16)=1 may be more efficient than
*                        IPARMQ(ISPEC=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice.)
*
*       NAME    (input) character string
*               Name of the calling subroutine
*
*       OPTS    (input) character string
*               This is a concatenation of the string arguments to
*               TTQRE.
*
*       N       (input) integer scalar
*               N is the order of the Hessenberg matrix H.
*
*       ILO     (input) INTEGER
*       IHI     (input) INTEGER
*               It is assumed that H is already upper triangular
*               in rows and columns 1:ILO-1 and IHI+1:N.
*
*       LWORK   (input) integer scalar
*               The amount of workspace available.
*
*  Further Details
*  ===============
*
*       Little is known about how best to choose these parameters.
*       It is possible to use different values of the parameters
*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*
*       It is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been
*       implemented --- yet.
*
*
*       The best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of xLAQR3 and xLAQR5 and on the nature of each
*       particular eigenvalue problem.  Experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       Following is a list of default values supplied by IPARMQ.
*       These defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*                        Default: 75. (Must be at least 11.)
*
*       IPARMQ(ISPEC=13) Recommended deflation window size.
*                        This depends on ILO, IHI and NS, the
*                        number of simultaneous shifts returned
*                        by IPARMQ(ISPEC=15).  The default for
*                        (IHI-ILO+1).LE.500 is NS.  The default
*                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*
*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*                        a multi-shift QR iteration.
*
*                        If IHI-ILO+1 is ...
*
*                        greater than      ...but less    ... the
*                        or equal to ...      than        default is
*
*                                0               30       NS =   2+
*                               30               60       NS =   4+
*                               60              150       NS =  10
*                              150              590       NS =  **
*                              590             3000       NS =  64
*                             3000             6000       NS = 128
*                             6000             infinity   NS = 256
*
*                    (+)  By default matrices of this order are
*                         passed to the implicit double shift routine
*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*                         values of NS are used only in case of a rare
*                         xLAHQR failure.
*
*                    (**) The asterisks (**) indicate an ad-hoc
*                         function increasing from 10 to 64.
*
*       IPARMQ(ISPEC=16) Select structured matrix multiply.
*                        (See ISPEC=16 above for details.)
*                        Default: 3.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $      NS = 4
         IF( NH.GE.60 )
     $      NS = 10
         IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
         IF( NH.GE.590 )
     $      NS = 64
         IF( NH.GE.3000 )
     $      NS = 128
         IF( NH.GE.6000 )
     $      NS = 256
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            IPARMQ = NS
         ELSE
            IPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
         IPARMQ = 0
         IF( NS.GE.KACMIN )
     $      IPARMQ = 1
         IF( NS.GE.K22MIN )
     $      IPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
         IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
C-----------------------------------------------------------------------
C  Routine:    IVOUT
C
C  Purpose:    Integer vector output routine.
C
C  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT)
C
C  Arguments
C     N      - Length of array IX. (Input)
C     IX     - Integer array to be printed. (Input)
C     IFMT   - Format to be used in printing array IX. (Input)
C     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
C              If IDIGIT .LT. 0, printing is done with 72 columns.
C              If IDIGIT .GT. 0, printing is done with 132 columns.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IVOUT (LOUT, N, IX, IDIGIT, IFMT)
C     ...
C     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IX(*), N, IDIGIT, LOUT
      CHARACTER  IFMT*(*)
C     ...
C     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, NDIGIT, K1, K2, LLL
      CHARACTER*80 LINE
*     ...
*     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
*
C
      LLL = MIN ( LEN ( IFMT ), 80 )
      DO 1 I = 1, LLL
          LINE(I:I) = '-'
    1 CONTINUE
C
      DO 2 I = LLL+1, 80
          LINE(I:I) = ' '
    2 CONTINUE
C
      WRITE ( LOUT, 2000 ) IFMT, LINE(1:LLL)
 2000 FORMAT ( /1X, A  /1X, A )
C
      IF (N .LE. 0) RETURN
      NDIGIT = IDIGIT
      IF (IDIGIT .EQ. 0) NDIGIT = 4
C
C=======================================================================
C             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
C=======================================================================
C
      IF (IDIGIT .LT. 0) THEN
C
      NDIGIT = -IDIGIT
      IF (NDIGIT .LE. 4) THEN
         DO 10 K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   10    CONTINUE
C
      ELSE IF (NDIGIT .LE. 6) THEN
         DO 30 K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
   30    CONTINUE
C
      ELSE IF (NDIGIT .LE. 10) THEN
         DO 50 K1 = 1, N, 5
            K2 = MIN0(N,K1+4)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
   50    CONTINUE
C
      ELSE
         DO 70 K1 = 1, N, 3
            K2 = MIN0(N,K1+2)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
   70    CONTINUE
      END IF
C
C=======================================================================
C             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
C=======================================================================
C
      ELSE
C
      IF (NDIGIT .LE. 4) THEN
         DO 90 K1 = 1, N, 20
            K2 = MIN0(N,K1+19)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   90    CONTINUE
C
      ELSE IF (NDIGIT .LE. 6) THEN
         DO 110 K1 = 1, N, 15
            K2 = MIN0(N,K1+14)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
  110    CONTINUE
C
      ELSE IF (NDIGIT .LE. 10) THEN
         DO 130 K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
  130    CONTINUE
C
      ELSE
         DO 150 K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
  150    CONTINUE
      END IF
      END IF
      WRITE (LOUT,1004)
C
 1000 FORMAT(1X,I4,' - ',I4,':',20(1X,I5))
 1001 FORMAT(1X,I4,' - ',I4,':',15(1X,I7))
 1002 FORMAT(1X,I4,' - ',I4,':',10(1X,I11))
 1003 FORMAT(1X,I4,' - ',I4,':',7(1X,I15))
 1004 FORMAT(1X,' ')
C
      RETURN
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      REAL             FUNCTION SECOND( )
*
*  -- LAPACK auxiliary routine (version 3.2) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*  Purpose
*  =======
*
*  SECOND returns nothing instead of returning the user time for a process in seconds.
*  If you are using that routine, it means that neither EXTERNAL ETIME,
*  EXTERNAL ETIME_, INTERNAL ETIME, INTERNAL CPU_TIME is available  on
*  your machine.
*
* =====================================================================
*
      SECOND = 0.0E+0
      RETURN
*
*     End of SECOND
*
      END

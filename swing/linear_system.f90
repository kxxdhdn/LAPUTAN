!******************************************************************************
!*
!*                          LINEAR SYSTEM SOLUTION
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007 
  !    - LINSYST_GAUSS & the LU routines updated on 11/2012
  ! 
  ! 3) DESCRIPTION: Solve a linear system of equations.
  !==========================================================================

MODULE linear_system

  USE utilities, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: linsyst_gauss, tridag, LU_decomp, linsyst_LU, cholesky_decomp, &
            linsyst_cholesky

  INTERFACE linsyst_gauss
    MODULE PROCEDURE linsyst_gauss_1D, linsyst_gauss_2D
  END INTERFACE linsyst_gauss


CONTAINS


  !==========================================================================
  ! X[N,M] = LINSYST_GAUSS_2D(A[N,N],B[N,M][,INVERSE[N,N],ONLY_INVERSE])
  !
  !   Solves the linear systems A[N,N].X[N,M] = B[N,M], using the Gauss-Jordan
  ! elimination with full pivot's algorithm (Num. Rec. 2007, p. 44).
  !   In other words, it simultaneously solves the M linear systems defined by
  ! A(:,:).X(:,1) = B(:,1), ..., A(:,:).X(:,M) = B(:,M). A natural product of 
  ! the operation is the inverse matrix of A (INVERSE). If ONLY_INVERSE is
  ! true, then only the inverse matrix is solved, the solution is not.
  !==========================================================================

  FUNCTION linsyst_gauss_2D (A,B,inverse,only_inverse)

    USE utilities, ONLY: DP, warning, swap, outerprod, outerand
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: A(:,:), B(:,:)
    REAL(DP), INTENT(OUT), OPTIONAL, DIMENSION(SIZE(A,1),SIZE(A,2)) :: inverse
    LOGICAL, INTENT(IN), OPTIONAL :: only_inverse
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: invA
    REAL(DP), DIMENSION(SIZE(B,1),SIZE(B,2)) :: linsyst_gauss_2D

    LOGICAL :: calcsol
    INTEGER, DIMENSION(SIZE(A,1)) :: indxc, indxr, ipiv
    REAL(DP), DIMENSION(SIZE(A,1)) :: dumc
    REAL(DP) :: pivinv
    INTEGER, TARGET, DIMENSION(2) :: irc
    INTEGER, POINTER :: irow, icol
    INTEGER :: N, M, i, l
    LOGICAL, DIMENSION(SIZE(A,1)) :: lpiv

    !-----------------------------------------------------------------------

    ! Initializing
    N = SIZE(B,1)
    M = SIZE(B,2)
    irow => irc(1)
    icol => irc(2)
    ipiv(:) = 0
    invA = A
    IF (PRESENT(only_inverse)) THEN
      calcsol = only_inverse
    ELSE
      calcsol = .True.
    END IF
    linsyst_gauss_2D = B

    ! Gauss-Jordan elimination with pivoting algorithm
    gauss_elim: DO i=1,N
      lpiv = (ipiv == 0)
      irc = MAXLOC(ABS(invA),OUTERAND(lpiv,lpiv))
      ipiv(icol) = ipiv(icol) + 1
      IF (ipiv(icol) > 1) CALL WARNING("LINSYST_GAUSS","Singular matrix")
      pivot: IF (irow /= icol) THEN
        CALL SWAP(invA(irow,:),invA(icol,:))
        IF (calcsol) &
          CALL SWAP(linsyst_gauss_2D(irow,:),linsyst_gauss_2D(icol,:))
      END IF pivot
      indxr(i) = irow
      indxc(i) = icol
      IF (invA(icol,icol) == 0._DP) &
        CALL WARNING("LINSYST_GAUSS","Singular matrix")
      pivinv = 1._DP / invA(icol,icol)
      invA(icol,icol) = 1._DP
      invA(icol,:) = invA(icol,:) * pivinv
      IF (calcsol) &
        linsyst_gauss_2D(icol,:) = linsyst_gauss_2D(icol,:) * pivinv
      dumc = invA(:,icol)
      invA(:,icol) = 0._DP
      invA(icol,icol) = pivinv
      invA(1:icol-1,:) = invA(1:icol-1,:) &
                         - OUTERPROD(dumc(1:icol-1),invA(icol,:))
      invA(icol+1:,:) = invA(icol+1:,:) &
                         - OUTERPROD(dumc(icol+1:),invA(icol,:)) 
      IF (calcsol) THEN
        linsyst_gauss_2D(1:icol-1,:) = linsyst_gauss_2D(1:icol-1,:) &
                            - OUTERPROD(dumc(1:icol-1),linsyst_gauss_2D(icol,:))
        linsyst_gauss_2D(icol+1:,:) = linsyst_gauss_2D(icol+1:,:) &
                            - OUTERPROD(dumc(icol+1:),linsyst_gauss_2D(icol,:))
      END IF
    END DO gauss_elim

    ! Unscramble the inverse matrix
    unscramble: IF (PRESENT(inverse)) THEN
      DO l=N,1,-1
        IF (indxr(l) /= indxc(l)) CALL SWAP(invA(:,indxr(l)),invA(:,indxc(l)))
      END DO
      inverse = invA
    END IF unscramble

    !-----------------------------------------------------------------------

  END FUNCTION linsyst_gauss_2D


  !==========================================================================
  ! Overload of LINSYST_GAUSS_2D in the case where X and B are of rank 1.
  !==========================================================================

  FUNCTION linsyst_gauss_1D (A,B,inverse)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: A(:,:), B(:)
    REAL(DP), INTENT(OUT), OPTIONAL, DIMENSION(SIZE(A,1),SIZE(A,2)) :: inverse
    REAL(DP), DIMENSION(SIZE(B)) :: linsyst_gauss_1D

    INTEGER :: N

    !-----------------------------------------------------------------------

    N = SIZE(B)

    IF (PRESENT(inverse)) THEN
      linsyst_gauss_1D = RESHAPE(LINSYST_GAUSS_2D(A,RESHAPE(B,[N,1]),inverse), &
                                 [N])
    ELSE
      linsyst_gauss_1D = RESHAPE(LINSYST_GAUSS_2D(A,RESHAPE(B,[N,1])),[N])
    END IF

    !-----------------------------------------------------------------------

  END FUNCTION linsyst_gauss_1D


  !==========================================================================
  ! X[N] = TRIDAG(A[N-1],B[N],C[N-1],R[N])
  !
  !   Solves a tridiagonal matrix system of diagonal B[N], off-diaognal 
  ! elements A[N-1] and C[N-1], and right-hand side R[N].
  !==========================================================================

  FUNCTION tridag (a,b,c,r)

    USE utilities, ONLY: DP, strike
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:), INTENT(IN) :: a, b, c, r
    REAL(DP), DIMENSION(SIZE(b)) :: tridag

    REAL(DP), DIMENSION(SIZE(b)) :: gam
    INTEGER :: N, j
    REAL(DP) :: bet

    !-----------------------------------------------------------------------

    ! Check input
    N = SIZE(b)
    IF (SIZE(r) /= N .OR. SIZE(a) /= N-1 .OR. SIZE(c) /= N-1) &
      CALL STRIKE("TRIDAG","wrong input.")
    bet = b(1)
    IF (bet == 0.0) CALL STRIKE("TRIDAG","system of order N-1.")

    ! Decomposition & substitution
    tridag(1) = r(1)/bet
    decomp: DO j=2,N
      gam(j) = c(j-1)/bet
      bet = b(j) - a(j-1)*gam(j)
      IF (bet == 0.0) CALL STRIKE("TRIDAG","system not solvable.")
      tridag(j) = (r(j)-a(j-1)*tridag(j-1))/bet
    END DO decomp
    backsubst: DO j=n-1,1,-1
      tridag(j) = tridag(j) - gam(j+1)*tridag(j+1)
    END DO backsubst

    !-----------------------------------------------------------------------

  END FUNCTION tridag


  !==========================================================================
  ! LU[N,N] = LU_DECOMP(A[N,N],indx[N],d)
  !
  !   L.U. decomposition of square matrix A, following Num. Rec. 2007, 
  ! Chap. 2.3.
  !==========================================================================

  FUNCTION LU_decomp (A,indx,d)

    USE utilities, ONLY: DP, warning, tinyDP, swap, outerprod
    IMPLICIT NONE
 
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: A
    INTEGER, INTENT(OUT), DIMENSION(SIZE(A,1)) :: indx
    INTEGER, INTENT(OUT) :: d
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: LU_decomp

    REAL(DP), DIMENSION(SIZE(A,1)) :: vv
    INTEGER :: j, imax, N

    !-----------------------------------------------------------------------

    d = 1
    N = SIZE(A,1)
    LU_decomp(:,:) = A(:,:)
    vv = MAXVAL(ABS(LU_decomp),DIM=2)
    IF (ANY(vv == 0._DP)) CALL WARNING("LU_DECOMP","Singular matrix")
    vv = 1._DP / vv
    DO j=1,N
      imax = (j-1) + MAXVAL(MAXLOC(vv(j:N)*ABS(LU_decomp(j:N,j))))
      IF (j /= imax) THEN
        CALL SWAP (LU_decomp(imax,:),LU_decomp(j,:))
        d = -d
        vv(imax) = vv(j)
      END IF
      indx(j) = imax
      IF (LU_decomp(j,j) == 0._DP) LU_decomp(j,j) = tinyDP
      LU_decomp(j+1:N,j) = LU_decomp(j+1:N,j) / LU_decomp(j,j)
      LU_decomp(j+1:N,j+1:N) = LU_decomp(j+1:N,j+1:N) &
                              - OUTERPROD(LU_decomp(j+1:N,j),LU_decomp(j,j+1:N))
    END DO

    !-----------------------------------------------------------------------

  END FUNCTION LU_decomp  

 
  !==========================================================================
  ! X[N] = LU_SOLVE(LU[N,N],B[N],indx[N])
  !
  !   L.U. decomposition for solving of linear system x=LU.b, following 
  ! Num. Rec. 2007, Chap. 2.3.
  !==========================================================================

  FUNCTION linsyst_LU (LU,b,indx)

    USE utilities, ONLY: DP
    IMPLICIT NONE
 
    REAL(DP), INTENT(IN), DIMENSION(:,:) :: LU
    REAL(DP), INTENT(IN), DIMENSION(SIZE(LU,1)) :: b
    INTEGER, INTENT(IN), DIMENSION(SIZE(LU,1)) :: indx
    REAL(DP), DIMENSION(SIZE(b)) :: linsyst_LU

    REAL(DP) :: summ
    INTEGER :: i, ii, ll, N

    !-----------------------------------------------------------------------

    N = SIZE(b)
    linsyst_LU(:) = b(:)
    ii = 0

    DO i=1,N
      ll = indx(i)
      summ = linsyst_LU(ll)
      linsyst_LU(ll) = linsyst_LU(i)
      IF (ii /= 0) THEN
        summ = summ - DOT_PRODUCT(LU(i,ii:i-1),linsyst_LU(ii:i-1))
      ELSE IF (summ /= 0._DP) THEN
        ii = i
      END IF
      linsyst_LU(i) = summ
    END DO
    DO i=N,1,-1
      linsyst_LU(i) = linsyst_LU(i) - DOT_PRODUCT(LU(i,i+1:N),linsyst_LU(i+1:N))
      linsyst_LU(i) = linsyst_LU(i) / LU(i,i)
    END DO

    !-----------------------------------------------------------------------

  END FUNCTION linsyst_LU
 
  !==========================================================================
  ! L[N,N] = CHOLESKY_DECOMP(A[N,N],NOPOSDEF)
  !
  !   Compute the Cholesky decomposition of the symmetric positive definite 
  ! matrix A[N,N]. L[N,N] is the lower triangle of the decomposition. In the
  ! end, A = L.L^T.
  !==========================================================================

  FUNCTION cholesky_decomp (A,noposdef)
  
    USE utilities, ONLY: DP, strike
    IMPLICIT NONE
  
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A
    LOGICAL, INTENT(OUT), OPTIONAL :: noposdef
    REAL(DP), DIMENSION(SIZE(A,1),SIZE(A,2)) :: cholesky_decomp

    REAL(DP), DIMENSION(SIZE(A,1)) :: diag
    REAL(DP) :: summ
    INTEGER :: i, N

    !-----------------------------------------------------------------------

    N = SIZE(A,1)
    cholesky_decomp = A
    IF (PRESENT(noposdef)) noposdef = .False.
    row: DO i=1,N
      summ = cholesky_decomp(i,i) &
             - DOT_PRODUCT(cholesky_decomp(i,1:i-1),cholesky_decomp(i,1:i-1))
      IF (summ <= 0._DP) THEN
        IF (PRESENT(noposdef)) THEN
          noposdef = .True.
          RETURN
        ELSE
          CALL STRIKE("CHOLESKY_DECOMP","Matrix is not positive definite")
        END IF
      END IF
      diag(i) = SQRT(summ)
      cholesky_decomp(i+1:N,i) = (cholesky_decomp(i,i+1:N) &
        - MATMUL(cholesky_decomp(i+1:N,1:i-1),cholesky_decomp(i,1:i-1)))/diag(i)
    END DO row

    ! Build the matrix
    lower_triangle: FORALL (i=1:N)
      cholesky_decomp(i,i+1:N) = 0._DP
      cholesky_decomp(i,i) = diag(i)
    END FORALL lower_triangle

    !-----------------------------------------------------------------------

  END FUNCTION cholesky_decomp


  !==========================================================================
  ! X[N] = LINSYST_CHOLESKY(L[N,N],b[N])
  !
  !   Solves the equation A[N,N].x[N]=b[N], given A and b, using the Cholesky de
  ! composition of A=L.L^T.
  !==========================================================================

  FUNCTION linsyst_cholesky (L,b)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), INTENT(IN) :: L
    REAL(DP), DIMENSION(SIZE(L,1)), INTENT(IN) :: b
    REAL(DP), DIMENSION(SIZE(L,1)) :: linsyst_cholesky 

    INTEGER :: i, N

    !-----------------------------------------------------------------------

    N = SIZE(L,1)
    DO i=1,N 
      linsyst_cholesky(i) = (b(i) &
                       - DOT_PRODUCT(L(i,1:i-1),linsyst_cholesky(1:i-1)))/L(i,i)
    END DO
    DO i=N,1,-1
      linsyst_cholesky(i) = (linsyst_cholesky(i) &
                       - DOT_PRODUCT(L(i+1:N,i),linsyst_cholesky(i+1:N)))/L(i,i)
    END DO

    !-----------------------------------------------------------------------

  END FUNCTION linsyst_cholesky


END MODULE linear_system

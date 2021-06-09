!******************************************************************************
!*
!*                    FUNCTIONS FOR NUMERICAL INTERPOLATION
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007 
  !    - Add interp_lin_sorted on 12/2012
  !    - Correct the starting indices of LOCATE_SORTED: jl=0->1, ju=N+1->N 
  !      (02/2014).
  !    - Allow zero values in the log case.
  ! 
  ! 3) DESCRIPTION: Set of functions performing numerical integration using
  !                 various methods.
  !==========================================================================

MODULE interpolation

  USE utilities, ONLY:
  USE arrays, ONLY:
  USE linear_system, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: interp_poly, interp_spline, interp_lin_sorted, locate_sorted

  INTERFACE interp_poly
    MODULE PROCEDURE interp_poly_v, interp_poly_s
  END INTERFACE interp_poly

  INTERFACE interp_lin_sorted
    MODULE PROCEDURE interp_lin_sorted_v, interp_lin_sorted_s
  END INTERFACE interp_lin_sorted


CONTAINS


  !==========================================================================
  ! index = LOCATE_SORTED(xx[N],x)
  !
  !   Finds the index position of the x value in the xx array, by bisection,
  ! provided that the array has been previously sorted.
  !==========================================================================

  PURE FUNCTION locate_sorted (xx,x)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: xx
    REAL(DP), INTENT(IN) :: x
    INTEGER :: locate_sorted

    INTEGER :: N, jl, jm, ju
    LOGICAL :: ascnd

    !-----------------------------------------------------------------------
    
    N = SIZE(xx)
    ascnd = (xx(N) >= xx(1))
    jl = 1
    ju = N
    midpoint: DO
      IF (ju-jl <= 1) EXIT
      jm = (ju+jl)/2
      IF (ascnd .EQV. (x >= xx(jm))) THEN ; jl = jm
                                     ELSE ; ju = jm ; END IF
    END DO midpoint  
    output: IF (x == xx(1)) THEN
      locate_sorted = 1
    ELSE IF (x == xx(N)) THEN
      locate_sorted = N-1
    ELSE
      locate_sorted = jl
    END IF output

    !-----------------------------------------------------------------------

  END FUNCTION locate_sorted


  !==========================================================================
  ! y_new[M] = INTERP_POLY(y_old[N],x_old[N],x_new[M],dy_new[M],DEGREE,XYLOG)
  !
  !   Evaluates the tabulated function (x_old,y_old) into the new array x_new
  ! and returns its values into y_new and error into err_y, using the 
  ! Neville's algorithm.
  !==========================================================================

  FUNCTION interp_poly_v (y_old,x_old,x_new,err_new,degree,xlog,ylog)
  
    USE utilities, ONLY: DP, strike, flEQ
    USE arrays, ONLY: closest
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: y_old, x_old
    REAL(DP), INTENT(IN), DIMENSION(:) :: x_new
    REAL(DP), INTENT(OUT), DIMENSION(SIZE(x_new)), OPTIONAL :: err_new
    INTEGER, INTENT(IN), OPTIONAL :: degree
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog
    REAL(DP), DIMENSION(SIZE(x_new)) :: interp_poly_v

    INTEGER :: i, i_new, m, Npoly, N, Na, jl
    INTEGER, DIMENSION(SIZE(x_new)) :: khi, klo
    REAL(DP), DIMENSION(SIZE(x_old)) :: xa, ya
    REAL(DP), DIMENSION(:), ALLOCATABLE :: C, D, A, dx
    REAL(DP), DIMENSION(SIZE(x_new)) :: x, y_new, dy_new
    LOGICAL :: xlog_flag, ylog_flag

    !-----------------------------------------------------------------------

    ! Check input
    !------------
    ! Sizes
    N = SIZE(x_new)
    Na = SIZE(x_old)
    IF (SIZE(y_old) /= Na) CALL STRIKE("INTERP_POLY","bad input.")
    IF (PRESENT(err_new)) dy_new = err_new

    ! Degree of the polynomial, default N=3
    IF (PRESENT(degree)) THEN ; Npoly = degree+1 ; ELSE ; Npoly = 3 ; END IF
    IF (Na < Npoly) CALL STRIKE("INTERP_POLY","there is not enough points.")
    ALLOCATE(C(Npoly),D(Npoly),A(Npoly),dx(Npoly))

    ! Type of axes
    IF (PRESENT(xlog)) THEN ; xlog_flag = xlog 
                       ELSE ; xlog_flag = .False. ; END IF
    IF (PRESENT(ylog)) THEN ; ylog_flag = ylog
                       ELSE ; ylog_flag = .False. ; END IF
    IF (xlog_flag .AND. (ANY(x_new <= 0.0) .OR. ANY(x_old <= 0.0))) &
      CALL STRIKE("INTERP_POLY","XLOG incompatible with X<=0.")
    IF (ylog_flag .AND. ANY(y_old <= 0.0)) &
      CALL STRIKE("INTERP_POLY","YLOG incompatible with Y<=0.")
    x = MERGE(LOG(x_new),x_new,xlog_flag)
    xa = MERGE(LOG(x_old),x_old,xlog_flag)
    ya = MERGE(LOG(y_old),y_old,ylog_flag)

    ! Interpolation grid
    !-------------------
    ! Locate the adjacent elements and reduce the search zone at each step,
    ! since the x array is sorted.
    DO i=1,N
      jl = MERGE(klo(i-1),1,i > 1)
      klo(i) = MAX(MIN(LOCATE_SORTED(xa(jl:Na),x(i))+jl-1,Na-1),1)
    END DO

    ! Distribute the points as evenly as possible around the value to be
    ! interpolated.
    klo(:) = klo(:) - (Npoly-1)/2
    khi(:) = klo(:) + (Npoly-1)
    ! Shift the values that are out of range, back inside the grid.
    edges: WHERE (klo <= 0)
      khi = khi - klo + 1
      klo = 1
    ELSEWHERE (khi > Na)
      klo = klo - (khi-Na)
      khi = Na
    END WHERE edges

    ! Polynomial interpolation
    !-------------------------
    bigloopx: DO i=1,N

      ! Degree 0
      C = ya(klo(i):khi(i))
      D = ya(klo(i):khi(i))
      dx = xa(klo(i):khi(i)) - x(i)
      i_new = CLOSEST(xa(klo(i):khi(i)),x(i))
      y_new(i) = ya(klo(i)+i_new-1)
      i_new = i_new-1

      ! Build the higher orders
      recursiv: DO m=1,Npoly-1
        A(1:Npoly-m) = dx(1:Npoly-m) - dx(1+m:Npoly)
        IF (ANY(A(1:Npoly-m) == 0._DP)) &
          CALL STRIKE("INTERP_POLY","negative denominator.")
        A(1:Npoly-m) = (C(2:Npoly-m+1)-D(1:Npoly-m)) / A(1:Npoly-m)
        D(1:Npoly-m) = dx(1+m:Npoly) * A(1:Npoly-m)
        C(1:Npoly-m) = dx(1:Npoly-m) * A(1:Npoly-m)
        middle: IF (2*i_new < Npoly-m) THEN
          dy_new(i) = C(i_new+1)
        ELSE
          dy_new(i) = D(i_new)
          i_new = i_new-1
        END IF middle
        y_new(i) = y_new(i) + dy_new(i)
      END DO recursiv

    END DO bigloopx

    ! Final unlogging
    interp_poly_v = MERGE(EXP(y_new),y_new,ylog_flag)
    IF (PRESENT(err_new)) err_new = dy_new
    DEALLOCATE (C,D,A,dx)

    !-----------------------------------------------------------------------

  END FUNCTION interp_poly_v

  !==========================================================================
  ! Same function when x_new is a scalar.
  !==========================================================================

  FUNCTION interp_poly_s (y_old,x_old,x_new,err_new,degree,xlog,ylog)
  
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: y_old, x_old
    REAL(DP), INTENT(IN) :: x_new
    REAL(DP), INTENT(OUT), OPTIONAL :: err_new
    INTEGER, INTENT(IN), OPTIONAL :: degree
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog
    REAL(DP) :: interp_poly_s

    INTEGER :: Npoly
    REAL(DP), DIMENSION(1) :: polint, err
    LOGICAL :: xlog_flag, ylog_flag

    !-----------------------------------------------------------------------

    IF (PRESENT(degree)) THEN ; Npoly = degree+1 ; ELSE ; Npoly = 3 ; END IF

    IF (PRESENT(xlog)) THEN ; xlog_flag = xlog
                       ELSE ; xlog_flag = .False. ; END IF
    IF (PRESENT(ylog)) THEN ; ylog_flag = ylog
                       ELSE ; ylog_flag = .False. ; END IF

    polint = INTERP_POLY_V(y_old,x_old,(/x_new/),err,Npoly,xlog_flag,ylog_flag)
    IF (PRESENT(err_new)) err_new = err(1)
    interp_poly_s = polint(1)

    !-----------------------------------------------------------------------

  END FUNCTION interp_poly_s


  !==========================================================================
  ! y_new[N] = SPLINE(y_old[N],x_old[N],yp1,ypn)
  !
  !   Returns the array that contains the second derivatives of the 
  ! interpolating function at the tabulated points x_old.
  !==========================================================================

  FUNCTION spline (x,y,yp1,ypn)

    USE utilities, ONLY: DP, strike
    USE linear_system, ONLY: tridag
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x, y
    REAL(DP), INTENT(IN) :: yp1, ypn
    REAL(DP), DIMENSION(SIZE(x)) :: spline

    INTEGER :: N
    REAL(DP), DIMENSION(SIZE(x)) :: a, b, c, r

    !-----------------------------------------------------------------------

    ! Check input
    N = SIZE(x)
    IF (SIZE(y) /= N) CALL STRIKE("SPLINE","bad input.")

    ! Set up the tridiagonal equations
    c(1:n-1) = x(2:n) - x(1:n-1)
    r(1:n-1) = 6.0_DP*((y(2:n)-y(1:n-1)) / c(1:n-1))
    r(2:n-1) = r(2:n-1) - r(1:n-2)
    a(2:n-1) = c(1:n-2)
    b(2:n-1) = 2.0_DP*(c(2:n-1)+a(2:n-1))
    b(1) = 1.0
    b(n) = 1.0
    lowbound: IF (yp1 > 0.99E30_DP) THEN
      r(1) = 0.0
      c(1) = 0.0
    ELSE
      r(1) = (3.0_DP/(x(2)-x(1))) * ((y(2)-y(1))/(x(2)-x(1)) - yp1)
      c(1) = 0.5
    END IF lowbound
    upbound: IF (ypn > 0.99E30_DP) THEN
      r(n) = 0.0
      a(n) = 0.0
    ELSE
      r(n) = (-3.0_DP/(x(n)-x(n-1))) * ((y(n)-y(n-1))/(x(n)-x(n-1)) - ypn)
      a(n) = 0.5
    END IF upbound

    ! Solve the system
    spline(1:N) = TRIDAG(a(2:n),b(1:n),c(1:n-1),r(1:n))

    !-----------------------------------------------------------------------

  END FUNCTION spline


  !==========================================================================
  ! y_new[M] = INTERP_SPLINE(y_old[N],x_old[N],x_new[M],XLOG=T/F,YLOG=T/F)
  !
  !   Evaluates the tabulated function (x_old,y_old) into the new array x_new
  ! and returns its values into y_new, using cubic spline interpolation. 
  ! x_new must be sorted and in ascending order.
  !==========================================================================

  FUNCTION interp_spline (y_old,x_old,x_new,xlog,ylog)

    USE utilities, ONLY: DP, strike, flEQ
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x_old, y_old
    REAL(DP), DIMENSION(:), INTENT(IN) :: x_new
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog
    REAL(DP), DIMENSION(SIZE(x_new)) :: interp_spline

    INTEGER, DIMENSION(SIZE(x_new)) :: khi, klo
    INTEGER :: Na, N, i, jl
    REAL(DP), DIMENSION(SIZE(x_old)) :: y2a, xa, ya
    REAL(DP), DIMENSION(SIZE(x_new)) :: h, a ,b, x
    LOGICAL :: xlog_flag, ylog_flag

    !-----------------------------------------------------------------------

    ! Check input
    !------------
    ! Sizes
    N = SIZE(x_new)
    Na = SIZE(x_old)
    IF (SIZE(y_old) /= Na) CALL STRIKE("INTERP_SPLINE","bad input.")

    ! Type of axes
    IF (PRESENT(xlog)) THEN
      xlog_flag = xlog
    ELSE 
      xlog_flag = .False.
    END IF
    IF (PRESENT(ylog)) THEN
      ylog_flag = ylog
    ELSE 
      ylog_flag = .False.
    END IF
    IF (xlog_flag .AND. (ANY(x_new <= 0.0) .OR. ANY(x_old <= 0.0))) &
      CALL STRIKE("INTERP_SPLINE","XLOG incompatible with X<=0.")
    IF (ylog_flag .AND. ANY(y_old <= 0.0)) &
      CALL STRIKE("INTERP_SPLINE","YLOG incompatible with Y<=0.")
    x = MERGE(LOG(x_new),x_new,xlog_flag)
    xa = MERGE(LOG(x_old),x_old,xlog_flag)
    ya = MERGE(LOG(y_old),y_old,ylog_flag)

    ! Interpolation grid
    !-------------------
    ! Locate the adjacent elements and reduce the search zone at each step,
    ! since the x array is sorted.
    DO i=1,N
      jl = MERGE(klo(i-1),1,i > 1)
      klo(i) = MAX(MIN(LOCATE_SORTED(xa(jl:Na),x(i))+jl-1,Na-1),1)
    END DO
    khi(:) = klo(:)+1
    
    ! Compute the spline
    !-------------------
    ! Finite difference
    h(:) = xa(khi(:))-xa(klo(:))
    IF (ANY(h == 0.0)) CALL STRIKE("INTERP_SPLINE","bad xa input")

    ! Spline calculation
    y2a(:) = SPLINE(xa,ya,1.E30_DP,1.E30_DP)
    a(:) = (xa(khi(:))-x(:)) / h(:)
    b(:) = (x(:)-xa(klo(:))) / h(:)
    interp_spline(:) = a(:)*ya(klo(:)) + b(:)*ya(khi(:)) &
                    + ( (a(:)**3-a(:))*y2a(klo(:)) &
                        + (b(:)**3-b(:))*y2a(khi(:)) )*h(:)**2 / 6.0_DP

    ! Final unlogging
    IF (ylog_flag) interp_spline(:) = EXP(interp_spline(:))

    !-----------------------------------------------------------------------

  END FUNCTION interp_spline 


  !==========================================================================
  ! y_new[M] = INTERP_LIN_SORTED(y_old[N],x_old[N],x_new[M],XLOG=T/F,YLOG=T/F,
  !                              FORCE=T/F)
  !   
  !   Performs fast log or linear interpolation on a previously sorted array.
  ! If log is used, the corresponding array must be positive. No control is 
  ! done. FORCE forces to zero the values that are NaN (usually LOG(0)).
  !==========================================================================

  PURE FUNCTION interp_lin_sorted_v (y_old,x_old,x_new,xlog,ylog,force)

    USE utilities, ONLY: DP, isNaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x_old, y_old
    REAL(DP), DIMENSION(:), INTENT(IN) :: x_new
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog, force
    REAL(DP), DIMENSION(SIZE(x_new)) :: interp_lin_sorted_v

    LOGICAL :: xl, yl, forced
    INTEGER :: i, Nnew, Nold
    INTEGER, DIMENSION(SIZE(x_new)) :: ju, jl
    REAL(DP), DIMENSION(SIZE(x_new)) :: x1
    REAL(DP), DIMENSION(SIZE(x_old)) :: x0, y0

    !-----------------------------------------------------------------------

    ! Preamble
    Nnew = SIZE(x_new)
    Nold = SIZE(x_old)
    xl = .False.
    yl = .False.
    IF (PRESENT(xlog)) xl = xlog
    IF (PRESENT(ylog)) yl = ylog    
    x0(:) = MERGE(LOG(x_old),x_old,xl)
    x1(:) = MERGE(LOG(x_new),x_new,xl)
    y0(:) = MERGE(LOG(y_old),y_old,yl)

    ! Find the right points to compute the slope
    FORALL (i=1:Nnew) jl(i) = LOCATE_SORTED(x0,x1(i))
    ju(:) = jl(:) + 1

    ! Interpolate
    interp_lin_sorted_v(:) = (y0(ju(:))-y0(jl(:)))/(x0(ju(:))-x0(jl(:))) &
                             * (x1(:)-x0(jl(:))) + y0(jl(:))

    ! Finalize
    interp_lin_sorted_v(:) = MERGE( EXP(interp_lin_sorted_v(:)), &
                                    interp_lin_sorted_v(:), yl )

    ! Deal with the potential zero/negative values if FORCE
    forced = .False.
    IF (PRESENT(force)) forced = force
    IF (forced .AND. xl) THEN
      WHERE (ISNAN(interp_lin_sorted_v(:))) interp_lin_sorted_v(:) = 0._DP
    END IF

    !-----------------------------------------------------------------------    

  END FUNCTION interp_lin_sorted_v

    !-----------------------------------------------------------------------    

  PURE FUNCTION interp_lin_sorted_s (y_old,x_old,x_new,xlog,ylog,force)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: x_old, y_old
    REAL(DP), INTENT(IN) :: x_new
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog, force
    REAL(DP) :: interp_lin_sorted_s
    LOGICAL :: xl, yl, forced

    !-----------------------------------------------------------------------

    xl = .False.
    yl = .False.
    forced = .False.
    IF (PRESENT(xlog)) xl = xlog
    IF (PRESENT(ylog)) yl = ylog    
    IF (PRESENT(force)) forced = force
    interp_lin_sorted_s = MINVAL(INTERP_LIN_SORTED_V(y_old,x_old,[x_new],xl,yl, &
                                                     forced))

    !-----------------------------------------------------------------------

  END FUNCTION interp_lin_sorted_s


END MODULE interpolation


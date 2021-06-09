!******************************************************************************
!*
!*                        INTEGRATION OF FUNCTIONS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007 
  !    - Correct bug of INTEG_TABULATED, when there is not point inside the
  ! interval to integrate over.
  !    - On 05/2014: add the log-lin integration.
  !    - On 05/2014: allow zero values in the function in the log case.
  ! 
  ! 3) DESCRIPTION: Various methods to integrate numerically a function.
  !==========================================================================

MODULE integration

  USE utilities, ONLY:
  USE arrays, ONLY:
  USE interpolation, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: dPrimitive, integ_tabulated, integ_romb


CONTAINS


  !==========================================================================
  ! dF[N] = DPRIMITIVE(x,f,XLOG,YLOG,lnx,lnf)
  !
  !   Returns the differential primitive of a function on a given grid, for 
  ! various function assumptions (lin/log).
  !   Log(x) and LOG(f) can eventually be passed as argument to avoid having
  ! to comput them twice.
  !==========================================================================

  PURE FUNCTION dPrimitive (x,f,xlog,ylog,lnx,lnf)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x, f
    LOGICAL, INTENT(IN) :: xlog, ylog
    REAL(DP), INTENT(IN), DIMENSION(:), OPTIONAL :: lnx, lnf
    REAL(DP), DIMENSION(SIZE(x)) :: dPrimitive
    
    INTEGER :: N
    REAL(DP), DIMENSION(SIZE(x)) :: logx, logf
    
    !------------------------------------------------------------------------

    N = SIZE(x(:))
    dPrimitive(1) = 0._DP
    linlog: IF (.NOT. xlog .AND. .NOT. ylog) THEN

      ! Linear: assumes f(x) ~ a*x + b 
      dPrimitive(2:N) = DPRIMITIVE_LINLIN(x(:),f(:),N)

    ELSE IF (xlog .AND. .NOT. ylog) THEN

      ! Logarithm: assumes f(x) ~ ln(a*x+b)
      dPrimitive(2:N) = DPRIMITIVE_LOGLIN(x(:),f(:),N)

    ELSE IF (.NOT. xlog .AND. ylog) THEN

      ! Exponential: assumes f(x) ~ exp(a*x+b)
      IF (PRESENT(lnf)) THEN 
        logf(:) = lnf(:) 
      ELSE 
        logf(:) = LOG(f(:))
      END IF
      dPrimitive(2:N) = DPRIMITIVE_LINLOG(x(:),f(:),N,lnf(:))
        
    ELSE IF (xlog .AND. ylog) THEN 

      ! Power-law: assumes f(x) ~ a*x^b
      IF (PRESENT(lnx)) THEN 
        logx(:) = lnx(:) 
      ELSE 
        logx(:) = LOG(x(:))
      END IF
      IF (PRESENT(lnf)) THEN 
        logf(:) = lnf(:) 
      ELSE 
        logf(:) = LOG(f(:))
      END IF
      dPrimitive(2:N) = DPRIMITIVE_LOGLOG(x(:),f(:),N,lnx(:),lnf(:))
      
    END IF linlog

    !------------------------------------------------------------------------

  END FUNCTION dPrimitive

  
  ! Formula for each method
  !------------------------
  ! Linear: assumes f(x) ~ a*x + b 
  !  => F(x1)-F(x0) = (f(x1)+f(x0)) / 2 * (x1-x0)
  PURE FUNCTION dPrimitive_linlin (x,f,N)
    
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(N), INTENT(IN) :: x, f
    INTEGER, INTENT(IN) :: N
    REAL(DP), DIMENSION(N-1) :: dPrimitive_linlin

    !------------------------------------------------------------------------

    dPrimitive_linlin(:) = 0.5_DP * ( f(2:N) + f(1:N-1) ) * ( x(2:N) - x(1:N-1) )

    !------------------------------------------------------------------------

  END FUNCTION dPrimitive_linlin
  
  ! Logarithm: assumes f(x) ~ ln(a*x+b)
  !  => F(x1)-F(x0) = ( (f(x1)-1)*EXP(f(x1)) - (f(x0)-1)*EXP(f(x0)) ) 
  !                 * (x1-x0) / ( EXP(f(x1)) - EXP(f(x0)) )
  PURE FUNCTION dPrimitive_loglin (x,f,N)
    
    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(N), INTENT(IN) :: x, f
    INTEGER, INTENT(IN) :: N
    REAL(DP), DIMENSION(N-1) :: dPrimitive_loglin
    REAL(DP), DIMENSION(N) :: expf, fm1
    REAL(DP), DIMENSION(N-1) :: dexpf
    
    !------------------------------------------------------------------------

    fm1 = f(:) - 1._DP
    expf = EXP(f(:))
    dexpf(:) = expf(2:N) - expf(1:N-1)
    WHERE (ABS(dexpf(:)) > epsDP*expf(2:N))
      dPrimitive_loglin(:) = ( fm1(2:N)*expf(2:N) - fm1(1:N-1)* expf(1:N-1) ) &
                           * ( x(2:N) - x(1:N-1) ) / dexpf(:)
    ELSEWHERE
      dPrimitive_loglin(:) = 0.5_DP * ( f(2:N) + f(1:N-1) ) * (x(2:N)-x(1:N-1))
    END WHERE

    !------------------------------------------------------------------------

  END FUNCTION dPrimitive_loglin
  
  ! Exponential: assumes f(x) ~ exp(a*x+b)
  !  => F(x1)-F(x0) = ( f(x1) - f(x0) ) * (x1-x0) / ( ln(f(x1)) - ln(f(x0)) )
  PURE FUNCTION dPrimitive_linlog (x,f,N,lnf)
    
    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(N), INTENT(IN) :: x, f, lnf
    INTEGER, INTENT(IN) :: N
    REAL(DP), DIMENSION(N-1) :: dPrimitive_linlog
    REAL(DP), DIMENSION(N-1) :: dlnf
    
    !------------------------------------------------------------------------

    dlnf(:) = lnf(2:N) - lnf(1:N-1)
    WHERE (f(1:N-1) > 0._DP .AND. f(2:N) > 0._DP &
           .AND. ABS(dlnf(:)) > epsDP*ABS(lnf(2:N)) )
      dPrimitive_linlog(:) = ( f(2:N) - f(1:N-1) ) * ( x(2:N) - x(1:N-1) ) &
                           / dlnf(:)                     
    ELSEWHERE
      dPrimitive_linlog(:) = 0.5_DP * ( f(2:N) + f(1:N-1) ) * (x(2:N)-x(1:N-1))
    END WHERE
      
    !------------------------------------------------------------------------

  END FUNCTION dPrimitive_linlog
  
  ! Power-law: assumes f(x) ~ a*x^b
  ! F(x1)-F(x0) = ( x1*f(x1) - x0*f(x0) ) / ( 1 + ( ln(f(x1)) - ln(f(x0) ) )
  !                                               / ( ln(x1) - ln(x0) ) )
  PURE FUNCTION dPrimitive_loglog (x,f,N,lnx,lnf)
    
    USE utilities, ONLY: DP, epsDP
    IMPLICIT NONE

    REAL(DP), DIMENSION(N), INTENT(IN) :: x, f, lnx, lnf
    INTEGER, INTENT(IN) :: N
    REAL(DP), DIMENSION(N-1) :: dPrimitive_loglog
    REAL(DP), DIMENSION(N-1) :: dlnx

    !------------------------------------------------------------------------

    dlnx(:) = lnx(2:N) - lnx(1:N-1)
    WHERE (f(1:N-1) > 0._DP .AND. f(2:N) > 0._DP &
           .AND. ABS(dlnx(:)) > epsDP*ABS(lnx(2:N)))
      dPrimitive_loglog(:) = ( x(2:N)*f(2:N) - x(1:N-1)*f(1:N-1) ) &
                           / ( 1._DP + ( lnf(2:N) - lnf(1:N-1) ) / dlnx(:) )
    ELSEWHERE
      dPrimitive_loglog(:) = 0.5_DP * ( f(2:N) + f(1:N-1) ) * (x(2:N)-x(1:N-1))
    END WHERE
      
    !------------------------------------------------------------------------

  END FUNCTION dPrimitive_loglog
  

  !==========================================================================
  ! I = INTEG_TABULATED(x[N],f[N],[x0,x1],PRIMITIVE[N],LINEAR,LOGARITHM, &
  !                     EXPONENTIAL,POWERLAW,XLOG,YLOG,FORCE,RESCALE)
  !
  !   Integration of tabulated functions. Returns the integration of f(x) from
  ! x0 to x1>x0, using the trapezium method. X must be SORTED without any 
  ! duplicates. The boolean keywords select the integration method (lin-lin, 
  ! log-lin, lin-log or log-log).
  !==========================================================================

  FUNCTION integ_tabulated (x,f,xrange,primitive,linear,logarithm,exponential, &
                            powerlaw,xlog,ylog,force,rescale)

    USE utilities, ONLY: DP, strike, warning, isNaN
    USE arrays, ONLY: iwhere

    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x, f
    REAL(DP), INTENT(IN), DIMENSION(2), OPTIONAL :: xrange
    REAL(DP), INTENT(OUT), DIMENSION(SIZE(x)), OPTIONAL :: primitive
    LOGICAL, INTENT(IN), OPTIONAL :: linear, logarithm, exponential, powerlaw
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog, force, rescale
    REAL(DP) :: integ_tabulated

    REAL(DP) :: x0, x1, min_x, max_x, x_inf, x_sup, f_inf, f_sup, logx2, logf2
    REAL(DP) :: scaling
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xi, fi, dPrim, logx, logf
    INTEGER :: Nx, Nin, i, i0, i1, i2
    LOGICAL, DIMENSION(SIZE(x)) :: boolin
    LOGICAL, DIMENSION(2) :: interpol, extrapol
    LOGICAL :: xl, yl, forced, rescl
    
    !-----------------------------------------------------------------------

    ! Check the arguments
    !--------------------
    ! Size
    Nx = SIZE(x(:))
    IF (Nx <= 1) CALL STRIKE("INTEG_TABULATED","x should have >=2 elements")

    ! Integration range
    min_x = MINVAL(x(:))
    max_x = MAXVAL(x(:))
    range: IF (PRESENT(xrange)) THEN
      x0 = xrange(1)
      x1 = xrange(2)
      IF ((x0 > max_x) .OR. (x1 < min_x)) &
        CALL STRIKE("INTEG_TABULATED","bad integration range.")
    ELSE
      x0 = min_x
      x1 = max_x
    END IF range
    IF (x0 >= x1) CALL STRIKE("INTEG_TABULATED","bad XRANGE.")

    ! Type of integration (xlog and ylog have priority over powerlaw, 
    ! exponential, linear, etc.)
    xl = .False.
    IF (PRESENT(xlog)) xl = xlog
    yl = .False.
    IF (PRESENT(ylog)) yl = ylog
    IF (PRESENT(linear)) THEN
      IF (linear) THEN
        xl = .False.
        yl = .False.
      END IF
    END IF
    IF (PRESENT(logarithm)) THEN
      IF (logarithm) THEN
        xl = .True.
        yl = .False.
      END IF
    END IF
    IF (PRESENT(exponential)) THEN
      IF (exponential) THEN
        xl = .False.
        yl = .True.
      END IF
    END IF
    IF (PRESENT(powerlaw)) THEN
      IF (powerlaw) THEN
        xl = .True.
        yl = .True.
      END IF
    END IF
    forced = .False.
    IF (PRESENT(force)) forced = force

    ! Rescale to avoid numerical problems
    rescl = .True.
    IF (PRESENT(rescale)) rescl = rescale


    ! Interpolation and extrapolation at the edges
    !---------------------------------------------
    ! Determine if extrapolation is necessary
    interpol(1) = ((.NOT. ANY(x(:) == x0)) .AND. (min_x < x0))
    interpol(2) = ((.NOT. ANY(x(:) == x1)) .AND. (max_x > x1))
    extrapol(1) = (min_x > x0)
    extrapol(2) = (max_x < x1)
    
    ! Logical mask of the points considered for integration
    boolin = (x(:) > x0 .AND. x(:) < x1)
    Nin = COUNT(boolin) + 2
    IF (rescl) THEN
      scaling = MERGE(MAXVAL(ABS(f(:)),MASK=( boolin(:) .AND. f(:) /= 0._DP &
                                              .AND. .NOT. isNaN(f(:)) )), &
                      1._DP,Nin > 2)
    ELSE
      scaling = 1._DP
    END IF
    ALLOCATE(xi(Nin),fi(Nin),dPrim(Nin),logx(Nin),logf(Nin))
    xi = [    x0, PACK(x(:),boolin),    x1 ]
    IF (.NOT. rescl) THEN
      fi = [ 1._DP, PACK(f(:),boolin), 1._DP ] ! edges will be replaced later
    ELSE
      fi = [ 1._DP, PACK(f(:),boolin)/scaling, 1._DP ]
    END IF
    logx = MERGE(LOG(xi(:)),xi(:),xl)
    logf = MERGE(LOG(fi(:)),fi(:),yl)
    
    ! Inter/extrapolates at the boundaries
    lower_boundary: IF (extrapol(1)) THEN
      logf(1) = logf(2) &
              - (logf(3)-logf(2)) / (logx(3)-logx(2)) * (logx(2)-logx(1))
    ELSE IF (interpol(1)) THEN
      x_inf = MAXVAL(PACK(x(:),x(:)<x0))
      IF (.NOT. rescl) THEN
        f_inf = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x_inf))), &
                      MINVAL(PACK(f(:),x(:)==x_inf)),yl)
      ELSE
        f_inf = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x_inf))/scaling), &
                      MINVAL(PACK(f(:),x(:)==x_inf))/scaling,yl)
      END IF
      x_inf = MERGE(LOG(x_inf),x_inf,xl)
      IF (Nin > 2) THEN
        logf(1) = logf(2) - (logf(2)-f_inf)/(logx(2)-x_inf) * (logx(2)-logx(1))
      ELSE
        CALL IWHERE(MERGE(LOG(x(:)),x(:),xl) == x_inf,i2)
        logx2 = MERGE(LOG(x(i2+1)),x(i2+1),xl)
        IF (.NOT. rescl) THEN
          logf2 = MERGE(LOG(f(i2+1)),f(i2+1),yl)
        ELSE
          logf2 = MERGE(LOG(f(i2+1)/scaling),f(i2+1)/scaling,yl)
        END IF
        logf(1) = logf2 - (logf2-f_inf)/(logx2-x_inf) * (logx2-logx(1))
      END IF
    ELSE 
      IF (.NOT. rescl) THEN
        logf(1) = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x0))), &
                        MINVAL(PACK(f(:),x(:)==x0)),yl)
      ELSE
        logf(1) = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x0))/scaling), &
                        MINVAL(PACK(f(:),x(:)==x0))/scaling,yl)
      END IF
    END IF lower_boundary
    fi(1) = MERGE(EXP(logf(1)),logf(1),yl)
    upper_boundary: IF (extrapol(2)) THEN
      logf(Nin) = logf(Nin-1) &
                + (logf(Nin-1)-logf(Nin-2)) / (logx(Nin-1)-logx(Nin-2)) &
                  * (logx(Nin)-logx(Nin-1))
    ELSE IF (interpol(2)) THEN
      x_sup = MINVAL(PACK(x(:),x(:)>x1))
      IF (.NOT. rescl) THEN
        f_sup = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x_sup))), &
                      MINVAL(PACK(f(:),x(:)==x_sup)),yl)
      ELSE
        f_sup = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x_sup))/scaling), &
                      MINVAL(PACK(f(:),x(:)==x_sup))/scaling,yl)
      END IF
      x_sup = MERGE(LOG(x_sup),x_sup,xl)
      IF (Nin > 2) THEN
        logf(Nin) = logf(Nin-1) &
                  + (f_sup-logf(Nin-1)) / (x_sup-logx(Nin-1)) &
                    * (logx(Nin)-logx(Nin-1))
      ELSE
        CALL IWHERE(MERGE(LOG(x(:)),x(:),xl) == x_sup,i2)
        logx2 = MERGE(LOG(x(i2-1)),x(i2-1),xl)
        IF (.NOT. rescl) THEN
          logf2 = MERGE(LOG(f(i2-1)),f(i2-1),yl)
        ELSE
          logf2 = MERGE(LOG(f(i2-1)/scaling),f(i2-1)/scaling,yl)
        END IF
        logf(Nin) = logf2 &
                  + (f_sup-logf2) / (x_sup-logx2) * (logx(Nin)-logx2)
      END IF
    ELSE 
      IF (.NOT. rescl) THEN
        logf(Nin) = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x1))), &
                          MINVAL(PACK(f(:),x(:)==x1)),yl)
      ELSE
        logf(Nin) = MERGE(LOG(MINVAL(PACK(f(:),x(:)==x1))/scaling), &
                          MINVAL(PACK(f(:),x(:)==x1))/scaling,yl)
      END IF
    END IF upper_boundary
    fi(Nin) = MERGE(EXP(logf(Nin)),logf(Nin),yl)


    ! Actual lin/log trapezium integration
    !-------------------------------------
    dPrim(:) = DPRIMITIVE(xi(:),fi(:),xl,yl,logx(:),logf(:))
    IF (rescl) dPrim(:) = dPrim(:) * scaling

    ! Remove the NaN
    IF (forced .AND. yl) THEN
      WHERE (ISNAN(dPrim(:))) dPrim(:) = 0._DP
    END IF

    ! Integral and primitive function if requested
    prim: IF (PRESENT(primitive)) THEN
      DO i=2,Nin ; dPrim(i) = dPrim(i) + dPrim(i-1) ; END DO
      integ_tabulated = dPrim(Nin)
      i0 = MERGE(2,1,(extrapol(1) .OR. interpol(1)))
      i1 = MERGE(Nin-1,Nin,(extrapol(2) .OR. interpol(2)))
      primitive(:) = UNPACK(dPrim(i0:i1),(x(:)>=x0 .AND. x(:)<=x1),FIELD=0._DP)
      primitive(:) = MERGE(integ_tabulated,primitive(:),x(:)>x1)
    ELSE 
      integ_tabulated = SUM(dPrim(:))
    END IF prim

    DEALLOCATE(xi,fi,dPrim,logx,logf)

    !-----------------------------------------------------------------------

  END FUNCTION integ_tabulated 


  !==========================================================================
  ! CALL TRAPEZ (func,a,b,s,n)
  !
  !   Returns the integration of f(x) from a to b>a, using the trapezium rule
  ! of order n.
  !==========================================================================

  SUBROUTINE trapez (func,a,b,s,n)
    
    USE utilities, ONLY: DP
    IMPLICIT NONE

    INTERFACE
      FUNCTION func (x)
        USE utilities, ONLY: DP
        REAL(DP), INTENT(IN), DIMENSION(:) :: x
        REAL(DP), DIMENSION(SIZE(x)) :: func        
      END FUNCTION func
    END INTERFACE
    REAL(DP), INTENT(IN) :: a, b
    REAL(DP), INTENT(INOUT) :: s
    INTEGER, INTENT(IN) :: n

    REAL(DP) :: del, fsum
    REAL(DP), DIMENSION(:), ALLOCATABLE :: arth
    INTEGER :: it, i

    !-----------------------------------------------------------------------

    IF (n == 1) THEN
      s = 0.5_DP*(b-a) * SUM(func((/a,b/)))
    ELSE
      it = 2**(n-2)
      del = (b-a)/it
      ALLOCATE(arth(it))
      FORALL (i=1:it) arth(i) = a+0.5_DP*del + (i-1)*del
      fsum = SUM(func(arth))
      s = 0.5_DP*(s+del*fsum)
      DEALLOCATE (arth)
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE trapez


  !==========================================================================
  ! I = INTEG_ROMB(func(x),a,b)
  !
  !   Returns the integration of f(x) from a to b>a, using the Romberg
  ! method of order k. 
  !==========================================================================

  FUNCTION integ_romb (func,a,b)

    USE utilities, ONLY: DP, strike, NaN
    USE arrays, ONLY: reverse
    USE interpolation, ONLY: interp_poly
    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: a, b
    INTERFACE
      FUNCTION func (x)
        USE utilities, ONLY: DP
        REAL(DP), INTENT(IN), DIMENSION(:) :: x
        REAL(DP), DIMENSION(SIZE(x)) :: func        
      END FUNCTION func
    END INTERFACE
    REAL(DP) :: integ_romb
    
    INTEGER, PARAMETER :: jmax=20, jmaxp=jmax+1, k=5, km=k-1
    REAL(DP), PARAMETER :: eps=1.E-6_DP

    REAL(DP), DIMENSION(jmaxp) :: h, s
    REAL(DP) :: dqromb
    INTEGER :: j

    !-----------------------------------------------------------------------

    integ_romb = NaN()
    h(1) = 1.0_DP
    DO j=1,jmax
      CALL TRAPEZ (func,a,b,s(j),j)
      IF (j >= k) THEN
        integ_romb = INTERP_POLY(REVERSE(s(j-km:j)),REVERSE(h(j-km:j)), &
                                 0._DP,ERR_NEW=dqromb,DEGREE=km-1)
        IF (ABS(dqromb) <= eps*ABS(integ_romb)) RETURN
      END IF
      s(j+1) = s(j)
      h(j+1) = 0.25_DP*h(j)
    END DO
    CALL STRIKE ("INTEG_ROMB","too many steps.")

    !-----------------------------------------------------------------------

  END FUNCTION integ_romb


END MODULE integration

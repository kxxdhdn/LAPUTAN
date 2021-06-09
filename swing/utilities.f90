!******************************************************************************
!*
!*                             Basic Utilities
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 07/2007 
  !    - Reallocate added 12/2012
  ! 
  ! 3) DESCRIPTION: Defines useful paths, variables and macros
  !==========================================================================

MODULE utilities

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: libF
  PUBLIC :: trimLR, trimeq, strike, warning, pring, swap, flEQ, incr, scl, NaN
  PUBLIC :: outerprod, outerand, arth, cumprod, cumsum, isNaN, isInf
  PUBLIC :: timinfo, timestring, zroots_unity, reallocate_pt, banner_program
  PUBLIC :: adjustc, strupcase, strlowcase, today, initiate_clock, strreplace

  ! General variables
  INTEGER, PARAMETER, PUBLIC :: LI = KIND(8)
  INTEGER, PARAMETER, PUBLIC :: DP = KIND(0D0)
  INTEGER, PARAMETER, PUBLIC :: CDP = KIND(CMPLX(0._DP,0._DP,DP))
  INTEGER, PARAMETER, PUBLIC :: lenstrnum = 20
  INTEGER, PARAMETER, PUBLIC :: ustd = 6
  INTEGER, PARAMETER, PUBLIC :: lenmax = 200
  REAL(DP), PARAMETER, PUBLIC :: tinyDP = TINY(0._DP), epsDP = EPSILON(0._DP)
  REAL(DP), PARAMETER, PUBLIC :: hugeDP = HUGE(0._DP)
  LOGICAL, SAVE, PUBLIC :: verbatim = .True.
  LOGICAL, SAVE, PUBLIC :: warnings = .True.
  CHARACTER(*), PARAMETER, PUBLIC :: bibi = "F. Galliano"
  CHARACTER(*), PARAMETER, PUBLIC :: programrunner = "D. Hu"

  ! Private parameters
  INTEGER, PARAMETER :: npar_arth = 16, npar2_arth = 8
  CHARACTER(*), PRIVATE, PARAMETER :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER(*), PRIVATE, PARAMETER :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  ! Procedure overload and interfaces
  INTERFACE pring
    MODULE PROCEDURE pring_int, pring_DP
  END INTERFACE pring

  INTERFACE trimeq
    MODULE PROCEDURE trimeq_scl, trimeq_1D
  END INTERFACE trimeq

  INTERFACE swap
    MODULE PROCEDURE swap_int, swap_dbl, swap_cdp 
    MODULE PROCEDURE swap_intarr, swap_dblarr, swap_cdparr
    MODULE PROCEDURE swap_intarr2, swap_dblarr2, swap_cdparr2
  END INTERFACE swap

  INTERFACE incr
    MODULE PROCEDURE incr_int, incr_dbl, incr_int1D, incr_dbl1D
    MODULE PROCEDURE incr_int2D, incr_dbl2D, incr_int3D, incr_dbl3D
  END INTERFACE incr

  INTERFACE scl
    MODULE PROCEDURE scl_int, scl_dbl, scl_int1D, scl_dbl1D
    MODULE PROCEDURE scl_int2D, scl_dbl2D, scl_int3D, scl_dbl3D
  END INTERFACE scl

  INTERFACE reallocate_pt
    MODULE PROCEDURE reallocate_ptint1D, reallocate_ptint2D
    MODULE PROCEDURE reallocate_ptDP1D, reallocate_ptDP2D
  END INTERFACE reallocate_pt

  INTERFACE flEQ
    MODULE PROCEDURE flEQ_scl, flEQ_1D, flEQ_1D_1D, flEQ_2D
  END INTERFACE flEQ

  INTERFACE cumsum
    MODULE PROCEDURE cumsum_DP, cumsum_int
  END INTERFACE

  INTERFACE cumprod
    MODULE PROCEDURE cumprod_DP, cumprod_int
  END INTERFACE

  INTERFACE arth
    MODULE PROCEDURE arth_dbl, arth_int
  END INTERFACE

  INTERFACE isNaN
    MODULE PROCEDURE isNaN_4D, isNaN_3D, isNaN_2D, isNaN_1D, isNaN_0D
  END INTERFACE isNaN

  INTERFACE isInf
    MODULE PROCEDURE isInf_v, isInf_s
  END INTERFACE isInf

  INTERFACE adjustc
    MODULE PROCEDURE adjustc_simple, adjustc_fixed
  END INTERFACE adjustc

  INTERFACE strreplace
    MODULE PROCEDURE strreplace_0_0, strreplace_0_1, strreplace_1_0
    MODULE PROCEDURE strreplace_1_1, strreplace_2_0, strreplace_2_1
    MODULE PROCEDURE strreplace_3_0, strreplace_3_1
  END INTERFACE strreplace

  ! Time structure
  TYPE, PUBLIC :: time_type
    INTEGER :: Ncycles
    REAL(DP) :: time0, rate, timemax
    LOGICAL :: cpu
  END TYPE time_type

  
CONTAINS


  !==========================================================================
  ! Specific directories of the library
  !
  ! dir = libF()      ! gives the machine dependent directory where library is.
  !==========================================================================

  ! Fortran library
  FUNCTION libF ()

    IMPLICIT NONE
    CHARACTER(lenmax) :: libF

    LibF = "/path/of/this/file/"

  END FUNCTION libF

  !==========================================================================
  ! string = TRIMLR(char)
  !
  !   Removes the leading and trailing blanks of a character string.
  !==========================================================================

  PURE FUNCTION trimLR (char)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: char
    CHARACTER(LEN_TRIM(ADJUSTL(char))) :: trimLR

    trimLR = TRIM(ADJUSTL(char))

  END FUNCTION trimLR


  !==========================================================================
  ! bool[N] = TRIMEQ(charray[N],char)
  !
  !   Tells where an array of strings is equal to a particular string, not
  ! accounting for leading and ending blanks.
  !==========================================================================

  ELEMENTAL FUNCTION trimeq_scl (charray,char)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: charray
    CHARACTER(*), INTENT(IN) :: char
    LOGICAL :: trimeq_scl

    trimeq_scl = ( TRIMLR(charray) == TRIMLR(char) )
    
  END FUNCTION trimeq_scl

    !----------------------------------------------------------------------

  PURE FUNCTION trimeq_1D (charray,char)

    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: charray
    CHARACTER(*), INTENT(IN) :: char
    LOGICAL, DIMENSION(SIZE(charray)) :: trimeq_1D
    INTEGER :: N, M

    N = LEN(charray)
    M = LEN_TRIM(char)
    IF (N-M >= 0) THEN
      trimeq_1D(:) = ( ADJUSTL(charray) == TRIMLR(char)//REPEAT(" ",N-M) )
    ELSE 
      trimeq_1D(:) = .False.
    END IF
    
  END FUNCTION trimeq_1D


  !==========================================================================
  ! string = PRING(number,Ndec)
  !
  !   Creates a string from a given number. NDEC is the number of decimal
  ! digits. 
  !==========================================================================

  FUNCTION pring_int (num)
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: num
    CHARACTER(lenstrnum) :: pring_int, form

    WRITE(form,"('(I',I2.2,')')") lenstrnum
    WRITE(pring_int,"(I20)") num

  END FUNCTION pring_int

    !----------------------------------------------------------------------

  FUNCTION pring_DP (num,Ndec)
  
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: num
    INTEGER, INTENT(IN), OPTIONAL :: Ndec
    CHARACTER(lenstrnum) :: pring_DP, form
    INTEGER :: N2

    decimal: IF (PRESENT(Ndec)) THEN
      N2 = Ndec
    ELSE 
      N2 = 2 
    END IF decimal

    format: IF ( (ABS(num) < 1.E4_DP .AND. ABS(num) > 1.E-4_DP) &
                 .OR. flEQ(num,0._DP,TOL=10**(-REAL(N2,DP))) ) THEN
      WRITE(form,"('(F',I2.2,'.',I2.2,')')") lenstrnum, N2
    ELSE
      WRITE(form,"('(ES',I2.2,'.',I2.2,')')") lenstrnum, N2      
    ENDIF format

    WRITE(pring_DP,form) num

  END FUNCTION pring_DP


  !==========================================================================
  ! CALL STRIKE(proc_name,comment)
  !
  !   Prints an error message and stops the code.
  !==========================================================================

  SUBROUTINE strike (proc,why)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: proc, why

    PRINT*, "!!! "//proc//" went on strike: "//why
    PRINT*
  PRINT*, "    /---\   /---\   |\  /|   |----       /---\   |    |   |----   |---\ "
  PRINT*, "    |       |   |   | \/ |   |           |   |   |    |   |       |   | "
  PRINT*, "    |  __   |---|   |    |   |---        |   |   \    /   |---    |---/ "
  PRINT*, "    |   |   |   |   |    |   |           |   |    \  /    |       |  \ "
  PRINT*, "    \---/   |   |   |    |   |----       \---/     \/     |----   |   \ "
  PRINT*
  PRINT*, "                         ***   INSERT COIN   ***"
  PRINT*
    STOP

  END SUBROUTINE strike


  !==========================================================================
  ! CALL WARNING(proc_name,comment)
  !
  !   Prints an error message but does not stop the code.
  !==========================================================================

  SUBROUTINE warning (proc,what)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: proc, what

    IF (warnings) PRINT*, "... "//proc//": "//what

  END SUBROUTINE warning


  !==========================================================================
  ! CALL SWAP(v1,v2)
  !
  !   Exchange the values of two variables.
  !==========================================================================

  ELEMENTAL SUBROUTINE swap_int (v1,v2)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: v1, v2
    INTEGER :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_int

    !----------------------------------------------------------------------

  ELEMENTAL SUBROUTINE swap_dbl (v1,v2)

    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: v1, v2
    REAL(DP) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_dbl

    !----------------------------------------------------------------------

  ELEMENTAL SUBROUTINE swap_cdp (v1,v2)

    IMPLICIT NONE
    COMPLEX(CDP), INTENT(INOUT) :: v1, v2
    COMPLEX(CDP) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_cdp

    !----------------------------------------------------------------------

  PURE SUBROUTINE swap_intarr (v1,v2)

    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: v1, v2
    INTEGER, DIMENSION(SIZE(v1)) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_intarr

    !----------------------------------------------------------------------

  PURE SUBROUTINE swap_dblarr (v1,v2)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: v1, v2
    REAL(DP), DIMENSION(SIZE(v1)) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_dblarr

    !----------------------------------------------------------------------

  PURE SUBROUTINE swap_cdparr (v1,v2)

    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:), INTENT(INOUT) :: v1, v2
    COMPLEX(CDP), DIMENSION(SIZE(v1)) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_cdparr

    !----------------------------------------------------------------------

  PURE SUBROUTINE swap_intarr2 (v1,v2)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: v1, v2
    INTEGER, DIMENSION(SIZE(v1,DIM=1),SIZE(v1,DIM=2)) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_intarr2

    !----------------------------------------------------------------------

  PURE SUBROUTINE swap_dblarr2 (v1,v2)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: v1, v2
    REAL(DP), DIMENSION(SIZE(v1,DIM=1),SIZE(v1,DIM=2)) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_dblarr2

    !----------------------------------------------------------------------

  PURE SUBROUTINE swap_cdparr2 (v1,v2)

    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:,:), INTENT(INOUT) :: v1, v2
    COMPLEX(CDP), DIMENSION(SIZE(v1,DIM=1),SIZE(v1,DIM=2)) :: temp
    
    temp = v1
    v1 = v2
    v2 = temp

  END SUBROUTINE swap_cdparr2


  !==========================================================================
  ! bool[N] = flEQ(v1[N,M],v2[N],TOL=1.E-5)
  !
  !   Determines equality of two real numbers within a tolerance interval.
  !==========================================================================

  PURE FUNCTION flEQ_2D (v1,v2,tol)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: v1
    REAL(DP), INTENT(IN) :: v2
    REAL(DP), INTENT(IN), OPTIONAL :: tol
    LOGICAL, DIMENSION(SIZE(v1,1),SIZE(v1,2)) :: flEQ_2D
    
    INTEGER, PARAMETER :: Nthresh = 5
    REAL(DP), PARAMETER :: epsil = Nthresh * epsDP
    REAL(DP), PARAMETER :: zero_min = Nthresh * tinyDP / epsDP
    REAL(DP), DIMENSION(SIZE(v1,1),SIZE(v1,2)) :: absv, threshold

    !----------------------------------------------------------------------

    absv(:,:) = MIN( ABS(v1(:,:)), ABS(v2) )
    precis: IF (.NOT. PRESENT(tol)) THEN
      threshold(:,:) = MERGE( epsil, epsil*absv(:,:), (absv(:,:) <= zero_min) )
    ELSE
      threshold(:,:) = tol * absv(:,:)
    END IF precis

    flEQ_2D(:,:) = ( ABS(v1(:,:)-v2) <= threshold(:,:) )

    !----------------------------------------------------------------------

  END FUNCTION flEQ_2D

  !------------------------------------------------------------------------

  PURE FUNCTION flEQ_1D_1D (v1,v2,tol)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: v1
    REAL(DP), DIMENSION(:),INTENT(IN) :: v2
    REAL(DP), INTENT(IN), OPTIONAL :: tol
    LOGICAL, DIMENSION(SIZE(v1)) :: flEQ_1D_1D
    
    INTEGER, PARAMETER :: Nthresh = 5
    REAL(DP), PARAMETER :: epsil = Nthresh * epsDP
    REAL(DP), PARAMETER :: zero_min = Nthresh * tinyDP / epsDP
    REAL(DP), DIMENSION(SIZE(v1)) :: absv, threshold

    !----------------------------------------------------------------------

    absv(:) = MIN( ABS(v1(:)), ABS(v2(:)) )
    precis: IF (.NOT. PRESENT(tol)) THEN
      threshold(:) = MERGE( epsil, epsil*absv(:), (absv(:) <= zero_min) )
    ELSE
      threshold(:) = tol * absv(:)
    END IF precis

    flEQ_1D_1D(:) = ( ABS(v1(:)-v2(:)) <= threshold(:) )

    !----------------------------------------------------------------------

  END FUNCTION flEQ_1D_1D

  ! Interface shells
  PURE FUNCTION flEQ_1D (v1,v2,tol)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: v1
    REAL(DP), INTENT(IN) :: v2
    REAL(DP), INTENT(IN), OPTIONAL :: tol
    LOGICAL, DIMENSION(SIZE(v1)) :: flEQ_1D
    flEQ_1D(:) = FLEQ_1D_1D(v1(:),SPREAD(v2,1,SIZE(v1(:))),tol)
  END FUNCTION flEQ_1D
  PURE FUNCTION flEQ_scl (v1,v2,tol)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: v1, v2
    REAL(DP), INTENT(IN), OPTIONAL :: tol
    LOGICAL :: flEQ_scl
    flEQ_scl = ANY(FLEQ_1D_1D([v1],[v2],tol))
  END FUNCTION flEQ_scl


  !==========================================================================
  ! CALL INCR(x[,val])
  !
  !   Increments x of 1 or val if set.
  !==========================================================================

  SUBROUTINE incr_int (x,val)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: val

    INTEGER :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1 ; END IF
    IF (warnings) THEN
      IF (ABS(x) > ABS(HUGE(x)-addval)) CALL WARNING("INCR","X is too high.")
    END IF
    x = x + addval

  END SUBROUTINE incr_int

    !----------------------------------------------------------------------

  SUBROUTINE incr_dbl (x,val)

    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: val

    REAL(DP) :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1._DP ; END IF
    IF (warnings) THEN
      IF (ABS(addval) < ABS(x*EPSILON(x))) CALL WARNING("INCR","X is too high.")
    END IF
    x = x + addval

  END SUBROUTINE incr_dbl

    !----------------------------------------------------------------------

  SUBROUTINE incr_int1D (x,val)

    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: val

    INTEGER :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1 ; END IF
    IF (warnings) THEN
      IF (ANY(ABS(x) > ABS(HUGE(x)-addval))) &
        CALL WARNING("INCR","X is too high.")
    END IF
    x(:) = x(:) + addval

  END SUBROUTINE incr_int1D

    !----------------------------------------------------------------------

  SUBROUTINE incr_dbl1D (x,val)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: val

    REAL(DP) :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1._DP ; END IF
    IF (warnings) THEN
      IF (ANY(ABS(addval) < ABS(x*EPSILON(x)))) &
        CALL WARNING("INCR","X is too high.")
    END IF
    x(:) = x(:) + addval

  END SUBROUTINE incr_dbl1D

    !----------------------------------------------------------------------

  SUBROUTINE incr_int2D (x,val)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: val

    INTEGER :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1 ; END IF
    IF (warnings) THEN
      IF (ANY(ABS(x) > ABS(HUGE(x)-addval))) &
        CALL WARNING("INCR","X is too high.")
    END IF
    x(:,:) = x(:,:) + addval

  END SUBROUTINE incr_int2D

    !----------------------------------------------------------------------

  SUBROUTINE incr_dbl2D (x,val)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: val

    REAL(DP) :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1._DP ; END IF
    IF (warnings) THEN
      IF (ANY(ABS(addval) < ABS(x*EPSILON(x)))) &
        CALL WARNING("INCR","X is too high.")
    END IF
    x(:,:) = x(:,:) + addval

  END SUBROUTINE incr_dbl2D

    !----------------------------------------------------------------------

  SUBROUTINE incr_int3D (x,val)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: val

    INTEGER :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1 ; END IF
    IF (warnings) THEN
      IF (ANY(ABS(x) > ABS(HUGE(x)-addval))) &
        CALL WARNING("INCR","X is too high.")
    END IF
    x(:,:,:) = x(:,:,:) + addval

  END SUBROUTINE incr_int3D

    !----------------------------------------------------------------------

  SUBROUTINE incr_dbl3D (x,val)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: val

    REAL(DP) :: addval

    IF (PRESENT(val)) THEN ; addval = val ; ELSE ; addval = 1._DP ; END IF
    IF (warnings) THEN
      IF (ANY(ABS(addval) < ABS(x*EPSILON(x)))) &
        CALL WARNING("INCR","X is too high.")
    END IF
    x(:,:,:) = x(:,:,:) + addval

  END SUBROUTINE incr_dbl3D


  !==========================================================================
  ! CALL SCL(x,fact)
  !
  !   Scale x by a factor FACT.
  !==========================================================================

  SUBROUTINE scl_int (x,fact)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: fact

    IF (warnings) THEN
      IF (ABS(x) > ABS(HUGE(x)/fact)) CALL WARNING("SCL","X is too high.")
    END IF
    x = x * fact

  END SUBROUTINE scl_int

    !----------------------------------------------------------------------

  SUBROUTINE scl_dbl (x,fact)

    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: fact

    x = x * fact

  END SUBROUTINE scl_dbl

    !----------------------------------------------------------------------

  SUBROUTINE scl_int1D (x,fact)

    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: fact

    IF (warnings) THEN
      IF (ANY(ABS(x) > ABS(HUGE(x)/fact))) CALL WARNING("SCL","X is too high.")
    END IF
    x(:) = x(:) * fact

  END SUBROUTINE scl_int1D

    !----------------------------------------------------------------------

  SUBROUTINE scl_dbl1D (x,fact)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: fact

    x(:) = x(:) * fact

  END SUBROUTINE scl_dbl1D

    !----------------------------------------------------------------------

  SUBROUTINE scl_int2D (x,fact)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: fact

    IF (warnings) THEN
      IF (ANY(ABS(x) > ABS(HUGE(x)/fact))) CALL WARNING("SCL","X is too high.")
    END IF
    x(:,:) = x(:,:) * fact

  END SUBROUTINE scl_int2D

    !----------------------------------------------------------------------

  SUBROUTINE scl_dbl2D (x,fact)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: fact

    x(:,:) = x(:,:) * fact

  END SUBROUTINE scl_dbl2D

    !----------------------------------------------------------------------

  SUBROUTINE scl_int3D (x,fact)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: fact

    IF (warnings) THEN
      IF (ANY(ABS(x) > ABS(HUGE(x)/fact))) CALL WARNING("SCL","X is too high.")
    END IF
    x(:,:,:) = x(:,:,:) * fact

  END SUBROUTINE scl_int3D

    !----------------------------------------------------------------------

  SUBROUTINE scl_dbl3D (x,fact)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN), OPTIONAL :: fact

    x(:,:,:) = x(:,:,:) * fact

  END SUBROUTINE scl_dbl3D


  !==========================================================================
  ! x[N,M] = OUTERPROD(a[N],b[M])
  !
  !   Computes the matrix X[N,M], the outer product of vectors A[N], B[M].
  !==========================================================================

  PURE FUNCTION outerprod (a,b)
  
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:) :: a, b
    REAL(DP), DIMENSION(SIZE(a,1),SIZE(b,1)) :: outerprod

    outerprod = SPREAD(a,DIM=2,NCOPIES=SIZE(b))*SPREAD(b,DIM=1,NCOPIES=SIZE(a))

  END FUNCTION outerprod

  !==========================================================================
  ! x[N,M] = OUTERAND(a[N],b[M])
  !
  !   Computes the matrix X[N,M], the outer logical of vectors A[N], B[M].
  !==========================================================================

  PURE FUNCTION outerand (a,b)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:), INTENT(IN) :: a, b
    LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand

    outerand = SPREAD(a,DIM=2,NCOPIES=SIZE(b)) .AND. & 
               SPREAD(b,DIM=1,NCOPIES=SIZE(a))

  END FUNCTION outerand


  !============================================================================
  ! pointer2 = REALLOCATE_PT(pointer1,N1,N2)
  !
  !   reallocate pointers (Num. Rec.).
  !============================================================================

  FUNCTION reallocate_ptint1D (p,n) 
    
    IMPLICIT NONE
    INTEGER, DIMENSION(:), POINTER :: p, reallocate_ptint1D
    INTEGER, INTENT(IN) :: N
    INTEGER :: Nold, ierr 

    !----------------------------------------------------------------------

    ALLOCATE (reallocate_ptint1D(N),STAT=ierr) 
    IF (ierr /= 0) &
      CALL STRIKE("REALLOCATE_PTINT1D","problem in attempt to allocate memory") 
    IF (.NOT. ASSOCIATED(p)) RETURN
    Nold = SIZE(p)
    reallocate_ptint1D(1:MIN(Nold,N)) = p(1:MIN(Nold,N))
    DEALLOCATE(p)

    !----------------------------------------------------------------------

  END FUNCTION reallocate_ptint1D

    !----------------------------------------------------------------------

  FUNCTION reallocate_ptint2D(p,n,m)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), POINTER :: p, reallocate_ptint2D 
    INTEGER, INTENT(IN) :: N, M
    INTEGER :: Nold, Mold, ierr 

    !----------------------------------------------------------------------

    ALLOCATE (reallocate_ptint2D(N,M),STAT=ierr)
    IF (ierr /= 0) CALL &
      STRIKE("REALLOCATE_PTINT2D","problem in attempt to allocate memory") 
    IF (.NOT. ASSOCIATED(p)) RETURN
    Nold = SIZE(p,1)
    Mold = SIZE(p,2)
    reallocate_ptint2D(1:MIN(Nold,N),1:MIN(Mold,M)) & 
      = p(1:MIN(Nold,N),1:MIN(Mold,M))
    DEALLOCATE(p)

    !----------------------------------------------------------------------
  
  END FUNCTION reallocate_ptint2D

    !----------------------------------------------------------------------

  FUNCTION reallocate_ptDP1D (p,N)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:), POINTER :: p, reallocate_ptDP1D 
    INTEGER, INTENT(IN) :: N
    INTEGER :: Nold, ierr 

    !----------------------------------------------------------------------
  
    ALLOCATE (reallocate_ptDP1D(N),STAT=ierr)
    IF (ierr /= 0) CALL &
      STRIKE("REALLOCATE_ptDP1D","problem in attempt to allocate memory")
    IF (.NOT. ASSOCIATED(p)) RETURN
    Nold = SIZE(p) 
    reallocate_ptDP1D(1:MIN(Nold,N)) = p(1:MIN(Nold,N)) 
    DEALLOCATE (p)

    !----------------------------------------------------------------------

  END FUNCTION reallocate_ptDP1D

    !----------------------------------------------------------------------

  FUNCTION reallocate_ptDP2D (p,N,M)

    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), POINTER :: p, reallocate_ptDP2D 
    INTEGER, INTENT(IN) :: N, M
    INTEGER :: Nold, Mold, ierr 
  
    !----------------------------------------------------------------------

    ALLOCATE (reallocate_ptDP2D(N,M),STAT=ierr)
    IF (ierr /= 0) CALL &
      STRIKE("REALLOCATE_PTDP2D","problem in attempt to allocate memory") 
    IF (.NOT. ASSOCIATED(p)) RETURN
    Nold = SIZE(p,1)
    Mold = SIZE(p,2)
    reallocate_ptDP2D(1:MIN(Nold,N),1:MIN(Mold,M)) &
      = p(1:MIN(Nold,N),1:MIN(Mold,M))
    DEALLOCATE(p)

    !----------------------------------------------------------------------

  END FUNCTION reallocate_ptDP2D


  !==========================================================================
  ! x[N] = ARTH(first,increment,N)
  !
  !   Returns an array containing an arithmetic progression.
  !==========================================================================

  PURE FUNCTION arth_dbl (first,increment,n)

    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: first, increment
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_dbl

    INTEGER :: k, k2
    REAL(DP) :: temp

    !----------------------------------------------------------------------

    IF (n > 0) arth_dbl(1) = first
    IF (n <= npar_arth) THEN
      DO k=2,n ; arth_dbl(k) = arth_dbl(k-1) + increment ; END DO
    ELSE
      DO k=2,npar2_arth ; arth_dbl(k) = arth_dbl(k-1) + increment ; END DO
      temp = increment * npar2_arth
      k = npar2_arth
      DO
        IF (k >= n) EXIT
        k2 = k + k
        arth_dbl(k+1:MIN(k2,n)) = temp + arth_dbl(1:MIN(k,n-k))
        temp = temp + temp
        k = k2
      END DO
    END IF

    !----------------------------------------------------------------------

  END FUNCTION arth_dbl

    !----------------------------------------------------------------------

  PURE FUNCTION arth_int (first,increment,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: first, increment, n
    INTEGER, DIMENSION(n) :: arth_int
    INTEGER :: k, k2, temp

    !----------------------------------------------------------------------

    IF (n > 0) arth_int(1) = first
    IF (n <= npar_arth) THEN
      DO k=2,n ; arth_int(k) = arth_int(k-1) + increment ; END DO
    ELSE
      DO k=2,npar2_arth ; arth_int(k) = arth_int(k-1) + increment ; END DO
      temp = increment * npar2_arth
      k = npar2_arth
      DO
        IF (k >= n) EXIT
        k2 = k + k
        arth_int(k+1:MIN(k2,n)) = temp + arth_int(1:MIN(k,n-k))
        temp = temp + temp
        k = k2
      END DO
    END IF

    !----------------------------------------------------------------------

  END FUNCTION arth_int


  !==========================================================================
  ! y[N] = CUMSUM(x[N],seed)
  !
  !   Cumulative sum of an array.
  !==========================================================================

  RECURSIVE FUNCTION cumsum_DP(arr,seed) RESULT(ans)
  
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    REAL(DP), OPTIONAL, INTENT(IN) :: seed
    REAL(DP), DIMENSION(SIZE(arr)) :: ans

    INTEGER, PARAMETER :: npar_cumsum = 16
    INTEGER :: n, j
    REAL(DP) :: sd

    !----------------------------------------------------------------------

    n = SIZE(arr)
    IF (n == 0) RETURN
      sd = 0._DP
    IF (PRESENT(seed)) sd = seed
    ans(1) = arr(1) + sd
    IF (n < npar_cumsum) THEN
      DO j=2,n ; ans(j) = ans(j-1) + arr(j) ; END DO
    ELSE
      ans(2:n:2) = CUMSUM_DP(arr(2:n:2)+arr(1:n-1:2),sd)
      ans(3:n:2) = ans(2:n-1:2) + arr(3:n:2)
    END IF

    !----------------------------------------------------------------------

  END FUNCTION cumsum_DP

  !------------------------------------------------------------------------

  RECURSIVE FUNCTION cumsum_int(arr,seed) RESULT(ans)
  
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, OPTIONAL, INTENT(IN) :: seed
    INTEGER, DIMENSION(SIZE(arr)) :: ans

    INTEGER, PARAMETER :: npar_cumsum = 16
    INTEGER :: n, j
    INTEGER :: sd

    !----------------------------------------------------------------------

    n = SIZE(arr)
    IF (n == 0) RETURN
      sd = 0
    IF (PRESENT(seed)) sd = seed
    ans(1) = arr(1) + sd
    IF (n < npar_cumsum) THEN
      DO j=2,n ; ans(j) = ans(j-1) + arr(j) ; END DO
    ELSE
      ans(2:n:2) = CUMSUM_INT(arr(2:n:2)+arr(1:n-1:2),sd)
      ans(3:n:2) = ans(2:n-1:2) + arr(3:n:2)
    END IF

    !----------------------------------------------------------------------

  END FUNCTION cumsum_int


  !==========================================================================
  ! y[N] = CUMPROD(x[N],seed)
  !
  !   Cumulative product of an array.
  !==========================================================================

  RECURSIVE FUNCTION cumprod_DP(arr,seed) RESULT(ans)
  
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    REAL(DP), OPTIONAL, INTENT(IN) :: seed
    REAL(DP), DIMENSION(SIZE(arr)) :: ans

    INTEGER, PARAMETER :: npar_cumprod = 8
    INTEGER :: n, j
    REAL(DP) :: sd

    !----------------------------------------------------------------------

    n = SIZE(arr)
    IF (n == 0) RETURN
      sd = 1.0_DP
    IF (PRESENT(seed)) sd = seed
    ans(1) = arr(1) * sd
    IF (n < npar_cumprod) THEN
      DO j=2,n ; ans(j) = ans(j-1) * arr(j) ; END DO
    ELSE
      ans(2:n:2) = CUMPROD_DP(arr(2:n:2)*arr(1:n-1:2),sd)
      ans(3:n:2) = ans(2:n-1:2) * arr(3:n:2)
    END IF

    !----------------------------------------------------------------------

  END FUNCTION cumprod_DP

  !------------------------------------------------------------------------

  RECURSIVE FUNCTION cumprod_int(arr,seed) RESULT(ans)
  
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, OPTIONAL, INTENT(IN) :: seed
    INTEGER, DIMENSION(SIZE(arr)) :: ans

    INTEGER, PARAMETER :: npar_cumprod = 8
    INTEGER :: n, j
    INTEGER :: sd

    !----------------------------------------------------------------------

    n = SIZE(arr)
    IF (n == 0) RETURN
      sd = 1
    IF (PRESENT(seed)) sd = seed
    ans(1) = arr(1) * sd
    IF (n < npar_cumprod) THEN
      DO j=2,n ; ans(j) = ans(j-1) * arr(j) ; END DO
    ELSE
      ans(2:n:2) = CUMPROD_int(arr(2:n:2)*arr(1:n-1:2),sd)
      ans(3:n:2) = ans(2:n-1:2) * arr(3:n:2)
    END IF

    !----------------------------------------------------------------------

  END FUNCTION cumprod_int


  !==========================================================================
  ! "" = TIMESTRING(time)
  !
  !   Print the time in h, m, s, etc.
  !==========================================================================

  FUNCTION timestring (time)

    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: time
    CHARACTER(15) :: timestring

    !----------------------------------------------------------------------

    IF (time < 60._DP) timestring = TRIMLR(PRING(NINT(time)))//"s"
    IF (time >= 60._DP .AND. time < 3600._DP) &
      timestring = TRIMLR(PRING(INT(time/60._DP)))//"m " &
                   //TRIMLR(PRING(INT(time-60._DP*INT(time/60._DP))))//"s"
    IF (time >= 3600._DP) &
      timestring = TRIMLR(PRING(INT(time/3600._DP)))//"h " &
                //TRIMLR(PRING(INT(time/60._DP-60._DP*INT(time/3600._DP))))//"m"

    !----------------------------------------------------------------------

  END FUNCTION timestring


  !==========================================================================
  ! "" = TIMINFO(time0,rate)
  !
  !   Print the time since time0, either CPU or wall-clock (if rate is entered).
  !==========================================================================

  FUNCTION timinfo (timestr)

    IMPLICIT NONE
    TYPE(time_type), INTENT(INOUT) :: timestr
    CHARACTER(15) :: timinfo

    INTEGER(LI) :: t1
    REAL(DP) :: dt, time1
    
    !----------------------------------------------------------------------

    IF (timestr%cpu) THEN
      CALL CPU_TIME(time1)
      dt = time1 - timestr%time0
    ELSE
      CALL SYSTEM_CLOCK(t1)
      time1 = REAL(t1,DP) + timestr%Ncycles * timestr%timemax
      IF (time1 < timestr%time0) THEN
        timestr%Ncycles = timestr%Ncycles + 1
        time1 = time1 + timestr%timemax
      END IF
      dt = ( time1 - timestr%time0 ) / timestr%rate
    END IF
    timinfo = TIMESTRING(dt)

    !----------------------------------------------------------------------

  END FUNCTION timinfo

    !----------------------------------------------------------------------

  SUBROUTINE initiate_clock (timestr,cpu)

    IMPLICIT NONE
    TYPE(time_type), INTENT(INOUT) :: timestr
    LOGICAL, INTENT(IN), OPTIONAL :: cpu

    INTEGER(LI) :: t0, r, cmax
    
    !----------------------------------------------------------------------

    IF (PRESENT(cpu)) THEN ; timestr%cpu = cpu
                      ELSE ; timestr%cpu = .False. ; END IF
    IF (timestr%cpu) THEN
      CALL CPU_TIME(timestr%time0)
    ELSE
      CALL SYSTEM_CLOCK(t0,COUNT_RATE=r,COUNT_MAX=cmax)
      timestr%time0 = REAL(t0,DP)
      timestr%rate = REAL(r,DP)
      timestr%timemax = REAL(cmax,DP) / REAL(r,DP)
      timestr%Ncycles = 0
    END IF
    
    !----------------------------------------------------------------------

  END SUBROUTINE initiate_clock

  
  !==========================================================================
  ! "" = TODAY()
  ! 
  !  Print the date
  !==========================================================================

  FUNCTION today()

    IMPLICIT NONE
    CHARACTER(80) :: today
    
    INTEGER, DIMENSION(8) :: values
    CHARACTER(2) :: minutes, hours
    CHARACTER(9), DIMENSION(12), PARAMETER :: &
      month = [ "January  ", "February ", "March    ", "April    ", &
                "May      ", "June     ", "July     ", "August   ", &
                "September", "October  ", "November ", "December " ]

    !----------------------------------------------------------------------

    CALL DATE_AND_TIME (VALUES=values)
    WRITE(hours,"(I2.2)") values(5)
    WRITE(minutes,"(I2.2)") values(6)
    today = TRIMLR(month(values(2)))//" " &
            //TRIMLR(PRING(values(3))) &
            //", "//TRIMLR(PRING(values(1))) &
            //", at "//hours//":"//minutes

    !----------------------------------------------------------------------

  END FUNCTION today


  !==========================================================================
  ! r = ZROOTS_UNITY(N,NN)
  !
  !   Complex function returning NN powers of the Nth root of unity.
  !==========================================================================

  FUNCTION zroots_unity(N,NN)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N, NN
    COMPLEX(CDP), DIMENSION(NN) :: zroots_unity
    INTEGER :: k
    REAL(DP) :: theta
    REAL(DP), PARAMETER :: twopi = 6.283185307179586476925286766559005768394_DP

    !----------------------------------------------------------------------

    zroots_unity(1) = 1._DP
    theta = twopi / N
    k = 1
    DO
      IF (k >= NN) EXIT
      zroots_unity(k+1) = CMPLX(COS(k*theta),SIN(k*theta),CDP)
      zroots_unity(k+2:MIN(2*k,NN)) = zroots_unity(k+1) &
                                    * zroots_unity(2:MIN(k,NN-k))
      k = 2 * k
    END DO

    !----------------------------------------------------------------------

  END FUNCTION zroots_unity


  !==========================================================================
  ! centered_text = ADJUSTC(text,textlen)
  !
  !   Center a string, containing leading or trailing blanks.
  !==========================================================================

  FUNCTION adjustc_simple (text)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: text
    CHARACTER(LEN(text)) :: adjustc_simple

    INTEGER :: strlen
    CHARACTER(12) :: form

    !----------------------------------------------------------------------

    strlen = LEN(text)
    WRITE(form,"('(A',I03,')')") strlen/2 + LEN(TRIMLR(text))/2
    WRITE(adjustc_simple,form) TRIMLR(text)

    !----------------------------------------------------------------------

  END FUNCTION adjustc_simple

  !------------------------------------------------------------------------

  FUNCTION adjustc_fixed (text,strlen)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: text
    INTEGER, INTENT(IN) :: strlen
    CHARACTER(strlen) :: adjustc_fixed

    CHARACTER(12) :: form

    !----------------------------------------------------------------------

    WRITE(form,"('(A',I03,')')") strlen/2 + LEN(TRIMLR(text))/2
    WRITE(adjustc_fixed,form) TRIMLR(text)

    !----------------------------------------------------------------------

  END FUNCTION adjustc_fixed


  !==========================================================================
  ! CALL BANNER_PROGRAM ('name',unit,SWING=T/F)
  !
  !   Print a generic banner at the start of a program.
  !==========================================================================

  SUBROUTINE banner_program (name,unit,swing)

    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: name
    INTEGER, INTENT(IN), OPTIONAL :: unit
    LOGICAL, INTENT(IN), OPTIONAL :: swing

    INTEGER, PARAMETER :: textwid = 80
    INTEGER :: u
    CHARACTER(12) :: form, form2
    CHARACTER(textwid) :: text
    LOGICAL :: sw

    !----------------------------------------------------------------------

    u = ustd
    sw = .False.
    IF (PRESENT(unit)) u = unit
    IF (PRESENT(swing)) sw = .True.

    WRITE(form2,"('(A',I2,')')") textwid
    WRITE(u,form2) REPEAT("=",textwid)
    text = "Running "//TRIMLR(name)
      WRITE(form,"('(A',I2,')')") textwid/2 + LEN(TRIMLR(text))/2
      WRITE(u,form) TRIMLR(text)
    IF (sw) THEN
      text = "from the SwING (SoftWares for Investigating Nebulae & Galaxies)"
        WRITE(form,"('(A',I2,')')") textwid/2 + LEN(TRIMLR(text))/2
        WRITE(u,form) TRIMLR(text)
    END IF
    text = "written by "//bibi
      WRITE(form,"('(A',I2,')')") textwid/2 + LEN(TRIMLR(text))/2
      WRITE(u,form) TRIMLR(text)
    WRITE(u,form2) REPEAT("=",textwid)
    WRITE(u,*)
    text = "("//TRIMLR(TODAY())//")"
      WRITE(form,"('(A',I2,')')") textwid
      WRITE(u,form) ADJUSTR(text)
    WRITE(u,*)

    !----------------------------------------------------------------------

  END SUBROUTINE banner_program


  !==========================================================================
  ! bool[N] = ISNAN(val[N])
  !
  !   Decide if a number is a NaN.
  !==========================================================================

  PURE FUNCTION isNaN_4D (val)
  
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: val
    LOGICAL, DIMENSION(SIZE(val,1),SIZE(val,2),SIZE(val,3),SIZE(val,4)) :: &
      isNaN_4D

    !----------------------------------------------------------------------

    isNaN_4D(:,:,:,:) = ( .NOT. ( ( val(:,:,:,:) <= 0._DP &
                                  .OR. val(:,:,:,:) >= 0._DP ) ) )

    !----------------------------------------------------------------------

  END FUNCTION isNaN_4D

  ! 3D
  PURE FUNCTION isNaN_3D (val)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: val
    LOGICAL, DIMENSION(SIZE(val,1),SIZE(val,2),SIZE(val,3)) :: isNaN_3D
    INTEGER :: N1, N2, N3
    N1 = SIZE(val(:,:,:),1)
    N2 = SIZE(val(:,:,:),2)
    N3 = SIZE(val(:,:,:),3)
    isNaN_3D(:,:,:) = RESHAPE(ISNAN_4D(RESHAPE(val(:,:,:),[N1,N2,N3,1])), &
                              [N1,N2,N3])
  END FUNCTION isNaN_3D

  ! 2D
  PURE FUNCTION isNaN_2D (val)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: val
    LOGICAL, DIMENSION(SIZE(val,1),SIZE(val,2)) :: isNaN_2D
    INTEGER :: N1, N2
    N1 = SIZE(val(:,:),1)
    N2 = SIZE(val(:,:),2)
    isNaN_2D(:,:) = RESHAPE(ISNAN_4D(RESHAPE(val(:,:),[N1,N2,1,1])),[N1,N2])
  END FUNCTION isNaN_2D

  ! 1D
  PURE FUNCTION isNaN_1D (val)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: val
    LOGICAL, DIMENSION(SIZE(val)) :: isNaN_1D
    INTEGER :: N1
    N1 = SIZE(val(:),1)
    isNaN_1D(:) = RESHAPE(ISNAN_4D(RESHAPE(val(:),[N1,1,1,1])),[N1])
  END FUNCTION isNaN_1D

  ! 0D
  PURE FUNCTION isNaN_0D (val)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: val
    LOGICAL :: isNaN_0D
    isNaN_0D = ALL(ISNAN_4D(RESHAPE([val],[1,1,1,1])))
  END FUNCTION isNaN_0D


  !==========================================================================
  ! NaN = NAN(val)
  ! 
  !  Returns a NaN with the same type as val
  !==========================================================================

  PURE FUNCTION NaN ()
  
    IMPLICIT NONE
    REAL(DP) :: NaN

    REAL(DP) :: a, b

    !----------------------------------------------------------------------

    a = 0._DP
    b = 0._DP
    NaN = a / b

    !----------------------------------------------------------------------

  END FUNCTION NaN


  !==========================================================================
  ! bool[N] = ISINF(val[N])
  !
  !   Decide if a number is a Infinity.
  !==========================================================================

  PURE FUNCTION isInf_v (val)
  
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: val
    LOGICAL, DIMENSION(SIZE(val)) :: isInf_v

    !----------------------------------------------------------------------

    isInf_v(:) = ( ABS(val(:)) >= hugeDP )

    !----------------------------------------------------------------------

  END FUNCTION isInf_v

    !----------------------------------------------------------------------

  PURE FUNCTION isInf_s (val)
  
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: val
    LOGICAL :: isInf_s

    !----------------------------------------------------------------------

    isInf_s = ANY(ISINF_V([val]))

    !----------------------------------------------------------------------

  END FUNCTION isInf_s



  !==========================================================================
  ! Conversion of a string to upper/lower case.
  !==========================================================================

  ELEMENTAL FUNCTION strupcase (input_string) 
    
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: input_string
    CHARACTER(LEN(input_string)) :: strupcase

    INTEGER :: i, n

    !----------------------------------------------------------------------

    strupcase = input_string
    DO i=1,LEN(strupcase)
      n = INDEX(lower_case,strupcase(i:i))
      IF (n /= 0) strupcase(i:i) = upper_case(n:n)
    END DO

    !----------------------------------------------------------------------

   END FUNCTION strupcase

    !----------------------------------------------------------------------

  ELEMENTAL FUNCTION strlowcase (input_string) 
    
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: input_string
    CHARACTER(LEN(input_string)) :: strlowcase

    INTEGER :: i, n

    !----------------------------------------------------------------------

    strlowcase = input_string
    DO i=1,LEN(strlowcase)
      n = INDEX(upper_case,strlowcase(i:i))
      IF (n /= 0) strlowcase(i:i) = lower_case(n:n)
    END DO

    !----------------------------------------------------------------------

   END FUNCTION strlowcase


  !==========================================================================
  ! Replace part of a string
  !==========================================================================

  FUNCTION strreplace_0_0 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*) :: s, text, repl
    CHARACTER(LEN(s)+100) :: strreplace_0_0 ! make sure it is long enough
    INTEGER :: i, Nt, Nr

    !----------------------------------------------------------------------

    strreplace_0_0 = s
    Nt = LEN_TRIM(text)
    Nr = LEN_TRIM(repl)
    DO
      i = INDEX(strreplace_0_0,text(:Nt))
      IF (i == 0) EXIT
      strreplace_0_0 = strreplace_0_0(:i-1) // repl(:nr) // strreplace_0_0(i+nt:)
    END DO

    !----------------------------------------------------------------------

  END FUNCTION strreplace_0_0
   
    !----------------------------------------------------------------------

  FUNCTION strreplace_0_1 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*) :: s
    CHARACTER(*), DIMENSION(:) :: text
    CHARACTER(*), DIMENSION(SIZE(text)) :: repl
    CHARACTER(LEN(s)+100) :: strreplace_0_1 ! make sure it is long enough
    INTEGER :: i

    !----------------------------------------------------------------------

    strreplace_0_1 = s
    DO i=1,SIZE(text)
      strreplace_0_1 = STRREPLACE_0_0(strreplace_0_1,text(i),repl(i))
    END DO

    !----------------------------------------------------------------------

  END FUNCTION strreplace_0_1
   
    !----------------------------------------------------------------------

  FUNCTION strreplace_1_0 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:) :: s
    CHARACTER(*) :: text, repl
    CHARACTER(LEN(s)+100), DIMENSION(SIZE(s)) :: strreplace_1_0
    INTEGER :: i
    
    !----------------------------------------------------------------------

    DO i=1,SIZE(s)
      strreplace_1_0(i) = STRREPLACE_0_0(s(i),text,repl)
    END DO
    
    !----------------------------------------------------------------------    
    
  END FUNCTION strreplace_1_0
  
    !----------------------------------------------------------------------

  FUNCTION strreplace_1_1 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:) :: s
    CHARACTER(*), DIMENSION(:) :: text
    CHARACTER(*), DIMENSION(SIZE(text)) :: repl
    CHARACTER(LEN(s)+100), DIMENSION(SIZE(s)) :: strreplace_1_1
    INTEGER :: i
    
    !----------------------------------------------------------------------

    DO i=1,SIZE(s)
      strreplace_1_1(i) = STRREPLACE_0_1(s(i),text,repl)
    END DO
    
    !----------------------------------------------------------------------    
    
  END FUNCTION strreplace_1_1

    !----------------------------------------------------------------------

  FUNCTION strreplace_2_0 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:,:) :: s
    CHARACTER(*) :: text, repl
    CHARACTER(LEN(s)+100), DIMENSION(SIZE(s,1),SIZE(s,2)) :: strreplace_2_0
    INTEGER :: i, j
    
    !----------------------------------------------------------------------

    DO i=1,SIZE(s,1)
      DO j=1,SIZE(s,2)
        strreplace_2_0(i,j) = STRREPLACE_0_0(s(i,j),text,repl)
      END DO
    END DO
    
    !----------------------------------------------------------------------    
    
  END FUNCTION strreplace_2_0
  
    !----------------------------------------------------------------------

  FUNCTION strreplace_2_1 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:,:) :: s
    CHARACTER(*), DIMENSION(:) :: text
    CHARACTER(*), DIMENSION(SIZE(text)) :: repl
    CHARACTER(LEN(s)+100), DIMENSION(SIZE(s,1),SIZE(s,2)) :: strreplace_2_1
    INTEGER :: i, j
    
    !----------------------------------------------------------------------

    DO i=1,SIZE(s,1)
      DO j=1,SIZE(s,2)
        strreplace_2_1(i,j) = STRREPLACE_0_1(s(i,j),text,repl)
      END DO
    END DO
    
    !----------------------------------------------------------------------    
    
  END FUNCTION strreplace_2_1

    !----------------------------------------------------------------------

  FUNCTION strreplace_3_0 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:,:,:) :: s
    CHARACTER(*) :: text, repl
    CHARACTER(LEN(s)+100), &
      DIMENSION(SIZE(s,1),SIZE(s,2),SIZE(s,3)) :: strreplace_3_0
    INTEGER :: i, j, k
    
    !----------------------------------------------------------------------

    DO i=1,SIZE(s,1)
      DO j=1,SIZE(s,2)
        DO k=1,SIZE(s,3)    
          strreplace_3_0(i,j,k) = STRREPLACE_0_0(s(i,j,k),text,repl)
        END DO
      END DO
    END DO
    
    !----------------------------------------------------------------------    
    
  END FUNCTION strreplace_3_0
  
    !----------------------------------------------------------------------

  FUNCTION strreplace_3_1 (s,text,repl)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:,:,:) :: s
    CHARACTER(*), DIMENSION(:) :: text
    CHARACTER(*), DIMENSION(SIZE(text)) :: repl
    CHARACTER(LEN(s)+100), &
      DIMENSION(SIZE(s,1),SIZE(s,2),SIZE(s,3)) :: strreplace_3_1
    INTEGER :: i, j, k
    
    !----------------------------------------------------------------------

    DO i=1,SIZE(s,1)
      DO j=1,SIZE(s,2)
        DO k=1,SIZE(s,3)
          strreplace_3_1(i,j,k) = STRREPLACE_0_1(s(i,j,k),text,repl)
        END DO
      END DO
    END DO
    
    !----------------------------------------------------------------------    
    
  END FUNCTION strreplace_3_1

  
END MODULE utilities

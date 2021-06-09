!******************************************************************************
!*
!*                 SIMPLE FUNCTIONS TO MANIPULATE ARRAYS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007 
  ! 
  ! 3) DESCRIPTION: Package for array manipulation.
  !==========================================================================

MODULE arrays

  USE utilities, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reallocate, ramp, ramp_step, print_vec, closest, sort, reverse
  PUBLIC :: sort_masked, uniq_sorted, iwhere, incrarr

  INTERFACE reallocate
    MODULE PROCEDURE reallocate_int1D, reallocate_int2D, reallocate_int3D
    MODULE PROCEDURE reallocate_int4D
    MODULE PROCEDURE reallocate_dbl1D, reallocate_dbl2D, reallocate_dbl3D
    MODULE PROCEDURE reallocate_dbl4D
    MODULE PROCEDURE reallocate_cmp1D, reallocate_cmp2D, reallocate_cmp3D
    MODULE PROCEDURE reallocate_cmp4D
    MODULE PROCEDURE reallocate_log1D, reallocate_log2D, reallocate_log3D
    MODULE PROCEDURE reallocate_log4D
    MODULE PROCEDURE reallocate_char1D, reallocate_char2D, reallocate_char3D
    MODULE PROCEDURE reallocate_char4D
  END INTERFACE reallocate

  INTERFACE closest
    MODULE PROCEDURE closest_int, closest_dbl_scl, closest_dbl_1D
  END INTERFACE closest

  INTERFACE uniq_sorted
    MODULE PROCEDURE uniq_sorted_int, uniq_sorted_dbl, uniq_sorted_dbl_pure
  END INTERFACE uniq_sorted

  INTERFACE iwhere
    MODULE PROCEDURE iwhere_general, iwhere_uniq
  END INTERFACE iwhere

  INTERFACE incrarr
    MODULE PROCEDURE incrarr_int, incrarr_DP, incrarr_bool, incrarr_char
    MODULE PROCEDURE incrarr1D_int, incrarr1D_DP, incrarr1D_bool, incrarr1D_char
    MODULE PROCEDURE incrarr2D_int, incrarr2D_DP, incrarr2D_bool, incrarr2D_char
  END INTERFACE incrarr


CONTAINS


  !============================================================================
  ! CALL REALLOCATE(array,N1,N2,N3,N4)
  !
  !   deallocate and reallocate an array of any kind and any profile.
  !============================================================================

  PURE SUBROUTINE reallocate_int1D (array,N1)

    IMPLICIT NONE
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1))

  END SUBROUTINE reallocate_int1D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_int2D (array,N1,N2)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2))

  END SUBROUTINE reallocate_int2D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_int3D (array,N1,N2,N3)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3))

  END SUBROUTINE reallocate_int3D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_int4D (array,N1,N2,N3,N4)

    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3, N4
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3,N4))

  END SUBROUTINE reallocate_int4D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_dbl1D (array,N1)
  
    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1))

  END SUBROUTINE reallocate_dbl1D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_dbl2D (array,N1,N2)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2))

  END SUBROUTINE reallocate_dbl2D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_dbl3D (array,N1,N2,N3)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3))

  END SUBROUTINE reallocate_dbl3D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_dbl4D (array,N1,N2,N3,N4)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3, N4
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3,N4))

  END SUBROUTINE reallocate_dbl4D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_cmp1D (array,N1)

    USE utilities, ONLY: CDP
    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1))

  END SUBROUTINE reallocate_cmp1D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_cmp2D (array,N1,N2)

    USE utilities, ONLY: CDP
    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2))

  END SUBROUTINE reallocate_cmp2D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_cmp3D (array,N1,N2,N3)

    USE utilities, ONLY: CDP
    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3))

  END SUBROUTINE reallocate_cmp3D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_cmp4D (array,N1,N2,N3,N4)

    USE utilities, ONLY: CDP
    IMPLICIT NONE
    COMPLEX(CDP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3, N4
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3,N4))

  END SUBROUTINE reallocate_cmp4D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_log1D (array,N1)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1))

  END SUBROUTINE reallocate_log1D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_log2D (array,N1,N2)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2))

  END SUBROUTINE reallocate_log2D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_log3D (array,N1,N2,N3)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3))

  END SUBROUTINE reallocate_log3D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_log4D (array,N1,N2,N3,N4)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3, N4
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3,N4))

  END SUBROUTINE reallocate_log4D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_char1D (array,N1)

    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1))

  END SUBROUTINE reallocate_char1D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_char2D (array,N1,N2)

    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2))

  END SUBROUTINE reallocate_char2D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_char3D (array,N1,N2,N3)

    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3))

  END SUBROUTINE reallocate_char3D

    !----------------------------------------------------------------------

  PURE SUBROUTINE reallocate_char4D (array,N1,N2,N3,N4)

    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: array
    INTEGER, INTENT(IN) :: N1, N2, N3, N4
    
    IF (ALLOCATED(array)) DEALLOCATE (array)
    ALLOCATE (array(N1,N2,N3,N4))

  END SUBROUTINE reallocate_char4D


  !==========================================================================
  ! r[N] = RAMP(N,Rinf,Rsup,XLOG=T/F,POWIND=)
  !
  !   Generates an array of N values increasing regularly from Rinf to Rsup
  !   either linearly or logarithmically.
  !==========================================================================

  PURE FUNCTION ramp (N,Rinf,Rsup,xlog,powind)
    
    USE utilities, ONLY: DP
    IMPLICIT NONE
  
    INTEGER, INTENT(IN) :: N                ! Number of elements of the array
    REAL(DP), INTENT(IN) :: Rinf, Rsup      ! Lower & upper bounds of the array
    REAL(DP), INTENT(IN),OPTIONAL :: powind ! Index of the power-law
    LOGICAL, INTENT(IN), OPTIONAL :: xlog   ! Types of progression

    REAL(DP), DIMENSION(N) :: ramp        ! Output array
    CHARACTER(3) :: mode
    INTEGER :: i

    !-----------------------------------------------------------------------

    mode = "LIN"
    IF (PRESENT(xlog)) mode = MERGE("LOG","LIN",xlog)
    IF (PRESENT(powind)) mode = "POW"

    generate: SELECT CASE (mode)
      CASE("LIN")
        ramp(:) = Rinf + (/(i,i=0,N-1)/)/REAL(N-1,DP)*(Rsup-Rinf)
      CASE("LOG")
        ramp(:) = Rinf * EXP( (/(i,i=0,N-1)/)/REAL(N-1,DP)*LOG(Rsup/Rinf) )
      CASE("POW")
        ramp(:) = ( Rinf**powind &
                  + (/(i,i=0,N-1)/)/REAL(N-1,DP) &
                    * (Rsup**powind-Rinf**powind) )**(1._DP/powind)
    END SELECT generate
    ramp(N) = Rsup ! in case there are cumulated rounding errors

    !-----------------------------------------------------------------------

  END FUNCTION ramp


  !==========================================================================
  ! CALL RAMP_STEP(Xinf,Xsup,X[N],DX,DlnX,N)
  !
  !   Generates an array X of increasing values from Xinf to Xsup, with a step
  ! DX. If DlNX is selected instead of DX then the grid is logarithmic. The
  ! number of points N is computed to ensure the minimum number of points
  ! with an actual step smaller or equal to the required one.
  !==========================================================================

  SUBROUTINE ramp_step (Xinf,Xsup,X,DX,DlnX,N)
    
    USE utilities, ONLY: DP, strike
    IMPLICIT NONE
  
    REAL(DP), INTENT(IN) :: Xinf, Xsup
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: X
    REAL(DP), INTENT(IN), OPTIONAL :: DX, DlnX
    INTEGER, INTENT(OUT), OPTIONAL :: N

    INTEGER :: Nx
    REAL(DP) :: step
    LOGICAL :: logmode

    !-----------------------------------------------------------------------

    ! Step and mode
    IF (PRESENT(DX)) THEN
      step = DX
      logmode = .False.
    ELSE IF (PRESENT(DlnX)) THEN
      step = DlnX
      logmode = .True.
    ELSE
      CALL STRIKE ("RAMP_STEP","Either DX or DlnX must be selected") 
      RETURN
    END IF 

    ! Ramp
    Nx = MERGE(FLOOR((LOG(Xsup)-LOG(Xinf))/step), &
               FLOOR((Xsup-Xinf)/step),logmode) + 1
    ALLOCATE (X(Nx))
    X(:) = RAMP(Nx,Xinf,Xsup,XLOG=logmode)
    IF (PRESENT(N)) N = Nx

    !-----------------------------------------------------------------------

  END SUBROUTINE ramp_step


  !==========================================================================
  ! i = CLOSEST(x[N],val)
  !
  !   Returns the index i corresponding to the closest x[i] to val.
  !==========================================================================

  PURE FUNCTION closest_int (x,val)

    IMPLICIT NONE
    INTEGER, INTENT(IN), DIMENSION(:) :: x
    INTEGER, INTENT(IN) :: val
    INTEGER :: closest_int

    closest_int = MINVAL(MINLOC(ABS(x-val)))

  END FUNCTION closest_int

    !-----------------------------------------------------------------------

  PURE FUNCTION closest_dbl_scl (x,val)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), INTENT(IN) :: val
    INTEGER :: closest_dbl_scl

    closest_dbl_scl = MINVAL(MINLOC(ABS(x-val)))

  END FUNCTION closest_dbl_scl

    !-----------------------------------------------------------------------

  PURE FUNCTION closest_dbl_1D (x,val)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:) :: x, val
    INTEGER, DIMENSION(SIZE(val)) :: closest_dbl_1D
    INTEGER :: N, i

    N = SIZE(val)
    FORALL (i=1:N) closest_dbl_1D(i) = MINVAL(MINLOC(ABS(x-val(i))))

  END FUNCTION closest_dbl_1D


  !==========================================================================
  ! bool = UNIQ_SORTED(x[N],Nsort(OUT),FIRST=T/F,LAST=T/F)
  !
  !  Returns the mask of the non repeated elements x[N], previously sorted.
  ! Nsort is the number of unique elements. Flags FIRST/LAST decide if the
  ! first or last occurence is kept.
  !==========================================================================

  FUNCTION uniq_sorted_int (x,Nsort,first,last)

    IMPLICIT NONE
    INTEGER, INTENT(IN), DIMENSION(:) :: x
    INTEGER, INTENT(OUT), OPTIONAL :: Nsort
    LOGICAL, INTENT(IN), OPTIONAL :: first, last
    LOGICAL, DIMENSION(SIZE(x)) :: uniq_sorted_int
    INTEGER :: N
    LOGICAL :: left

    !-----------------------------------------------------------------------

    N = SIZE(x)
    left = .True.
    IF (PRESENT(first)) left = first
    IF (PRESENT(last)) left = (.NOT. last)
    IF (left) THEN
      uniq_sorted_int(1) = .True.
      uniq_sorted_int(2:N) = (x(2:N) /= x(1:N-1))
    ELSE
      uniq_sorted_int(N) = .True.
      uniq_sorted_int(1:N-1) = (x(2:N) /= x(1:N-1))
    END IF
    IF (PRESENT(Nsort)) Nsort = SIZE(PACK(uniq_sorted_int,uniq_sorted_int))

  END FUNCTION uniq_sorted_int

    !-----------------------------------------------------------------------

  FUNCTION uniq_sorted_dbl (x,Nsort,first,last)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    INTEGER, INTENT(OUT) :: Nsort
    LOGICAL, INTENT(IN), OPTIONAL :: first, last
    LOGICAL, DIMENSION(SIZE(x)) :: uniq_sorted_dbl
    INTEGER :: N
    LOGICAL :: left

    !-----------------------------------------------------------------------

    N = SIZE(x)
    left = .True.
    IF (PRESENT(first)) left = first
    IF (PRESENT(last)) left = (.NOT. last)
    IF (left) THEN
      uniq_sorted_dbl(1) = .True.
      uniq_sorted_dbl(2:N) = (x(2:N) /= x(1:N-1))
    ELSE
      uniq_sorted_dbl(N) = .True.
      uniq_sorted_dbl(1:N-1) = (x(2:N) /= x(1:N-1))
    END IF
    Nsort = SIZE(PACK(uniq_sorted_dbl,uniq_sorted_dbl))

  END FUNCTION uniq_sorted_dbl

   !-----------------------------------------------------------------------

  PURE FUNCTION uniq_sorted_dbl_pure (x,first,last)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: first, last
    LOGICAL, DIMENSION(SIZE(x)) :: uniq_sorted_dbl_pure
    INTEGER :: N
    LOGICAL :: left

    !-----------------------------------------------------------------------

    N = SIZE(x)
    left = .True.
    IF (PRESENT(first)) left = first
    IF (PRESENT(last)) left = (.NOT. last)
    IF (left) THEN
      uniq_sorted_dbl_pure(1) = .True.
      uniq_sorted_dbl_pure(2:N) = (x(2:N) /= x(1:N-1))
    ELSE
      uniq_sorted_dbl_pure(N) = .True.
      uniq_sorted_dbl_pure(1:N-1) = (x(2:N) /= x(1:N-1))
    END IF

  END FUNCTION uniq_sorted_dbl_pure


  !==========================================================================
  ! list_sorted = SORT(list,IND=indices)
  !
  !   Routine of quick sorting using the Haore algorithm, from: Brainerd, 
  ! W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150, modified
  ! by Alan Miller.
  !
  !   Warning: this routine crashes if there are NaNs in the list...
  !==========================================================================

  FUNCTION sort (list,ind)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: list
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: ind
    REAL(DP), DIMENSION(SIZE(list)) :: sort

    INTEGER :: i, N
    REAL(DP), DIMENSION(SIZE(list)) :: listmp
    INTEGER, DIMENSION(SIZE(list)) :: indtmp

    !-----------------------------------------------------------------------

    N = SIZE(list)
    indtmp = [( i, i=1,N )]
    listmp = list
    CALL HOARE_METHOD (listmp,indtmp,1,N)
    sort = listmp
    IF (PRESENT(ind)) ind = indtmp

  END FUNCTION sort

    !-----------------------------------------------------------------------

  RECURSIVE SUBROUTINE hoare_method (list,ind,left_end,right_end)

    USE utilities, ONLY: DP, swap
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: list
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ind
    INTEGER, INTENT(IN) :: left_end, right_end

    INTEGER :: i, j
    REAL(DP) :: reference
    INTEGER, PARAMETER :: max_simple_sort_size = 6

    !-----------------------------------------------------------------------

    method: IF (right_end < left_end+max_simple_sort_size) THEN
      
      ! 1) Use interchange sort for small lists
      CALL sort_short(list,ind,left_end,right_end)

    ELSE

      ! 2) Use partition ("quick") sort
      reference = list((left_end+right_end)/2)
      i = left_end-1
      j = right_end+1
      
      scan: DO 
        ! Scan list from left end until element >= reference is found
        DO
          i = i+1
          IF (list(i) >= reference) EXIT
        END DO
        ! Scan list from right end until element <= reference is found
        DO
          j = j-1
          IF (list(j) <= reference) EXIT
        END DO
        IF (i < j) THEN
          ! Swap two out-of-order elements
          CALL SWAP (list(i),list(j))
          CALL SWAP (ind(i),ind(j))
        ELSE IF (i == j) THEN
          i = i+1
          EXIT
        ELSE
          EXIT
        END IF
      END DO scan
      
      IF (left_end < j) CALL hoare_method (list,ind,left_end,j)
      IF (i < right_end) CALL hoare_method (list,ind,i,right_end)

    END IF method

  END SUBROUTINE hoare_method

    !-----------------------------------------------------------------------

  SUBROUTINE sort_short(list,ind,left_end,right_end)

    USE utilities, ONLY: DP, swap
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(INOUT) :: list
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ind
    INTEGER, INTENT(IN) :: left_end, right_end

    INTEGER :: i, j

    !-----------------------------------------------------------------------

    DO i=left_end,right_end-1
      DO j=i+1,right_end
        IF (list(i) > list(j)) THEN
          CALL SWAP(list(i),list(j))
          CALL SWAP(ind(i),ind(j))
        END IF
      END DO
    END DO

    !-----------------------------------------------------------------------

  END SUBROUTINE sort_short


  !==========================================================================
  ! CALL SORT_MASKED(array,mask,array_sorted,IND=indices,N=N)
  !
  !   Same as SORT, but adding a MASK. As a consequence, the size of the array
  ! changes (it is N<=SIZE(array)), and therefore, ARRAY_SORTED needs to be
  ! allocatable, and we need to use a subroutine.
  !==========================================================================

  SUBROUTINE sort_masked (array,mask,array_sorted,ind,N) ! TO BE TESTED...

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: array
    LOGICAL, DIMENSION(:), INTENT(IN) :: mask
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: array_sorted
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: ind
    INTEGER, INTENT(OUT), OPTIONAL :: N

    INTEGER :: i, N0, Nm
    REAL(DP), DIMENSION(:), ALLOCATABLE :: arraym
    INTEGER, DIMENSION(SIZE(array)) :: ind0
    INTEGER, DIMENSION(:), ALLOCATABLE :: indm

    !-----------------------------------------------------------------------

    ! Apply the mask
    N0 = SIZE(array(:))
    Nm = COUNT(mask(:))
    ALLOCATE (arraym(Nm))
    ALLOCATE (indm(Nm))
    ind0(:) = [( i, i=1,N0 )]
    indm(:) = PACK(ind0(:),mask(:))
    arraym(:) = PACK(array(:),mask(:))

    ! Sort the masked array
    CALL HOARE_METHOD (arraym(:),indm(:),1,Nm)
    
    ! Optional outputs
    IF (PRESENT(array_sorted)) THEN
      ALLOCATE (array_sorted(Nm))
      array_sorted(:) = array(indm(:))
    END IF
    IF (PRESENT(ind)) THEN
      ALLOCATE (ind(Nm))
      ind(:) = indm(:)
    END IF
    IF (PRESENT(N)) N = Nm
    
  END SUBROUTINE sort_masked

  
  !==========================================================================
  ! CALL PRINT_VEC(v1[N1],v2[N2],v3[N3],v4[N4],PRECISION=2)
  !
  !   Prints one or several vectors in column.
  !==========================================================================

  SUBROUTINE print_vec (vec1,vec2,vec3,vec4,precision)
    
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: vec1
    REAL(DP), INTENT(IN), DIMENSION(:), OPTIONAL :: vec2, vec3, vec4
    INTEGER, INTENT(IN), OPTIONAL :: precision     ! Number of decimal figures

    REAL(DP), ALLOCATABLE, DIMENSION(:) :: lvec1, lvec2, lvec3, lvec4
    INTEGER :: N1, N2, N3, N4, Nmax
    INTEGER :: i, prec
    CHARACTER(2), PARAMETER :: colwid='15'         ! Width of each column
    CHARACTER(80) :: formtot
    CHARACTER(19) :: formv, formb, formv1, formv2, formv3, formv4

    !-----------------------------------------------------------------------

    ! Homogenize the sizes
    N1 = SIZE(vec1)
    IF (PRESENT(vec2)) THEN ; N2 = SIZE(vec2) ; ELSE ; N2 = 0 ; END IF
    IF (PRESENT(vec3)) THEN ; N3 = SIZE(vec3) ; ELSE ; N3 = 0 ; END IF
    IF (PRESENT(vec4)) THEN ; N4 = SIZE(vec4) ; ELSE ; N4 = 0 ; END IF
    Nmax = MAXVAL((/N1,N2,N3,N4/))
    ALLOCATE(lvec1(Nmax),lvec2(Nmax),lvec3(Nmax),lvec4(Nmax))
    
    ! General format
    prec = MERGE(precision,2,PRESENT(precision))
    WRITE(formv,"(A5,I1)") "ES"//colwid//".", prec
    formb = "TL"//colwid//","//colwid//"(' ')"

    ! Printing loop
    printing: DO i=1,Nmax

      lvec1(i) = MERGE(vec1(i),-1._DP,i<=N1)
      lvec2(i) = MERGE(vec2(i),-1._DP,i<=N2)
      lvec3(i) = MERGE(vec3(i),-1._DP,i<=N3)
      lvec4(i) = MERGE(vec4(i),-1._DP,i<=N4)

      formv1 = MERGE(REPEAT(" ",LEN_TRIM(formb)+1),","//TRIM(formb),i<=N1)
      formv1 = TRIM(formv)//formv1
      formv2 = MERGE(REPEAT(" ",LEN_TRIM(formb)+1),","//TRIM(formb),i<=N2)
      formv2 = TRIM(formv)//formv2
      formv3 = MERGE(REPEAT(" ",LEN_TRIM(formb)+1),","//TRIM(formb),i<=N3)
      formv3 = TRIM(formv)//formv3
      formv4 = MERGE(REPEAT(" ",LEN_TRIM(formb)+1),","//TRIM(formb),i<=N4)
      formv4 = TRIM(formv)//formv4
      formtot = "(I5,"//TRIM(formv1)//","//TRIM(formv2)//","//TRIM(formv3)// &
                ","//TRIM(formv4)//")"

      WRITE(*,formtot) i, lvec1(i), lvec2(i), lvec3(i), lvec4(i)

    END DO printing

    DEALLOCATE(lvec1,lvec2,lvec3,lvec4)

    !-----------------------------------------------------------------------

  END SUBROUTINE print_vec


  !==========================================================================
  ! x[N] = REVERSE(x[N])
  !
  !   Returns the reversed array.
  !==========================================================================

  PURE FUNCTION reverse (x)
    
    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x
    REAL(DP), DIMENSION(SIZE(x)) :: reverse

    INTEGER :: N

    !-----------------------------------------------------------------------

    N = SIZE(x)
    reverse(:) = x(N:1:-1)

    !-----------------------------------------------------------------------

  END FUNCTION reverse


  !==========================================================================
  ! CALL IWHERE(bool[N],ind[M])
  !
  !   Return the indices where BOOL is true, in IND.
  !==========================================================================

  SUBROUTINE iwhere_general (bool,ind)

    USE utilities, ONLY: DP, strike
    IMPLICIT NONE

    LOGICAL, DIMENSION(:), INTENT(IN) :: bool
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ind

    INTEGER :: i, N, M

    !-----------------------------------------------------------------------

    N = SIZE(bool)
    M = COUNT(bool)
    IF (M > 0) THEN
      ALLOCATE (ind(M))
      ind(:) = PACK([(i,i=1,N)],bool)
    ELSE
      CALL STRIKE ("IWHERE","no occurences found")
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE iwhere_general

    !-----------------------------------------------------------------------
    !   In case we are sure that there is only one true element of BOOL, 
    ! ind can be passed as a scalar (simplification).
    !-----------------------------------------------------------------------

  PURE SUBROUTINE iwhere_uniq (bool,ind)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    LOGICAL, DIMENSION(:), INTENT(IN) :: bool
    INTEGER, INTENT(OUT) :: ind

    INTEGER :: i, N, M

    !-----------------------------------------------------------------------

    N = SIZE(bool)
    M = COUNT(bool)
    IF (M == 1) THEN
      ind = MINVAL(PACK([(i,i=1,N)],bool))
    ELSE
      ind = 0
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE iwhere_uniq


  !==========================================================================
  ! CALL INCRARR (arr(1D or 2D),val(scal or 1D))
  !
  !   1D: ARR[N] -> ARR[N+1] & put VAL into ARR[N+1]
  !     : ARR[N] -> ARR[N+M] & put val[M] into ARR[N+1:N+M]
  !   2D: ARR[N,M] -> ARR[N+1,M] & put VAL[M] into ARR[N+1,M]
  !==========================================================================

  SUBROUTINE incrarr_int (arr,val)

    IMPLICIT NONE

    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
    INTEGER, INTENT(IN) :: val

    INTEGER :: N
    INTEGER, DIMENSION(:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (.NOT. ALLOCATED(arr)) THEN
      ALLOCATE (arr(1))
      arr(:) = val
    ELSE
      N = SIZE(arr)
      ALLOCATE (arr0(N))
      arr0(:) = arr(:)
      CALL REALLOCATE (arr,N+1)
      arr(1:N) = arr0(:)
      arr(N+1) = val
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr_int

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr_DP (arr,val)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
    REAL(DP), INTENT(IN) :: val

    INTEGER :: N
    REAL(DP), DIMENSION(:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (.NOT. ALLOCATED(arr)) THEN
      ALLOCATE (arr(1))
      arr(:) = val
    ELSE
      N = SIZE(arr)
      ALLOCATE (arr0(N))
      arr0(:) = arr(:)
      CALL REALLOCATE (arr,N+1)
      arr(1:N) = arr0(:)
      arr(N+1) = val
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr_DP

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr_bool (arr,val)

    IMPLICIT NONE

    LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
    LOGICAL, INTENT(IN) :: val

    INTEGER :: N
    LOGICAL, DIMENSION(:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (.NOT. ALLOCATED(arr)) THEN
      ALLOCATE (arr(1))
      arr(:) = val
    ELSE
      N = SIZE(arr)
      ALLOCATE (arr0(N))
      arr0(:) = arr(:)
      CALL REALLOCATE (arr,N+1)
      arr(1:N) = arr0(:)
      arr(N+1) = val
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr_bool

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr_char (arr,val)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
    CHARACTER(*), INTENT(IN) :: val

    INTEGER :: N
    CHARACTER(LEN(arr)), DIMENSION(:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (.NOT. ALLOCATED(arr)) THEN
      ALLOCATE (arr(1))
      arr(:) = val
    ELSE
      N = SIZE(arr(:))
      ALLOCATE (arr0(N))
      arr0(:) = arr(:)
      CALL REALLOCATE (arr,N+1)
      arr(1:N) = arr0(:)
      arr(N+1) = val
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr_char

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr1D_int (arr1D,vec)

    IMPLICIT NONE

    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr1D
    INTEGER, DIMENSION(:), INTENT(IN) :: vec

    INTEGER :: N, M
    INTEGER, DIMENSION(:), ALLOCATABLE :: arr0 
    
    !-----------------------------------------------------------------------

    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr1D)) THEN
      ALLOCATE (arr1D(M))
      arr1D(:) = vec(:)
    ELSE
      N = SIZE(arr1D)
      ALLOCATE (arr0(N))
      arr0(:) = arr1D(:)
      CALL REALLOCATE (arr1D,N+M)
      arr1D(1:N) = arr0(:)
      arr1D(N+1:N+M) = vec(:)
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr1D_int

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr1D_DP (arr1D,vec)

    USE utilities, ONLY: DP
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr1D
    REAL(DP), DIMENSION(:), INTENT(IN) :: vec

    INTEGER :: N, M
    REAL(DP), DIMENSION(:), ALLOCATABLE :: arr0 
    
    !-----------------------------------------------------------------------

    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr1D)) THEN
      ALLOCATE (arr1D(M))
      arr1D(:) = vec(:)
    ELSE
      N = SIZE(arr1D)
      ALLOCATE (arr0(N))
      arr0(:) = arr1D(:)
      CALL REALLOCATE (arr1D,N+M)
      arr1D(1:N) = arr0(:)
      arr1D(N+1:N+M) = vec(:)
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr1D_DP

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr1D_bool (arr1D,vec)

    IMPLICIT NONE

    LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr1D
    LOGICAL, DIMENSION(:), INTENT(IN) :: vec

    INTEGER :: N, M
    LOGICAL, DIMENSION(:), ALLOCATABLE :: arr0 
    
    !-----------------------------------------------------------------------

    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr1D)) THEN
      ALLOCATE (arr1D(M))
      arr1D(:) = vec(:)
    ELSE
      N = SIZE(arr1D)
      ALLOCATE (arr0(N))
      arr0(:) = arr1D(:)
      CALL REALLOCATE (arr1D,N+M)
      arr1D(1:N) = arr0(:)
      arr1D(N+1:N+M) = vec(:)
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr1D_bool

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr1D_char (arr1D,vec)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr1D
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: vec

    INTEGER :: N, M
    CHARACTER(LEN(arr1D)), DIMENSION(:), ALLOCATABLE :: arr0 
    
    !-----------------------------------------------------------------------

    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr1D)) THEN
      ALLOCATE (arr1D(M))
      arr1D(:) = vec(:)
    ELSE
      N = SIZE(arr1D(:))
      ALLOCATE (arr0(N))
      arr0(:) = arr1D(:)
      CALL REALLOCATE (arr1D,N+M)
      arr1D(1:N) = arr0(:)
      arr1D(N+1:N+M) = vec(:)
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr1D_char

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr2D_int (arr2D,vec,missing)

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: arr2D
    INTEGER, DIMENSION(:), INTENT(IN) :: vec
    INTEGER, INTENT(IN), OPTIONAL :: missing

    INTEGER :: N, M, Mold, Mmax, miss
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (PRESENT(missing)) THEN
      miss = missing
    ELSE
      miss = -1
    END IF
    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr2D)) THEN
      ALLOCATE (arr2D(1,M))
      arr2D(1,:) = vec(:)
    ELSE
      N = SIZE(arr2D,1)
      Mold = SIZE(arr2D,2)
      Mmax = MAX(M,Mold)
      ALLOCATE (arr0(N,Mold))
      arr0(:,:) = arr2D(:,:)
      CALL REALLOCATE (arr2D,N+1,Mmax)
      arr2D(1:N,1:Mold) = arr0(:,:)
      IF (Mold < Mmax) arr2D(1:N,Mold+1:Mmax) = miss
      arr2D(N+1,1:M) = vec(:)
      IF (M < Mmax) arr2D(N+1,M+1:Mmax) = miss
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr2D_int

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr2D_DP (arr2D,vec,missing)

    USE utilities, ONLY: DP, NaN
    IMPLICIT NONE

    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: arr2D
    REAL(DP), DIMENSION(:), INTENT(IN) :: vec
    REAL(DP), INTENT(IN), OPTIONAL :: missing

    INTEGER :: N, M, Mold, Mmax
    REAL(DP) :: miss
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (PRESENT(missing)) THEN
      miss = missing
    ELSE
      miss = NaN()
    END IF
    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr2D)) THEN
      ALLOCATE (arr2D(1,M))
      arr2D(1,:) = vec(:)
    ELSE
      N = SIZE(arr2D,1)
      Mold = SIZE(arr2D,2)
      Mmax = MAX(M,Mold)
      ALLOCATE (arr0(N,Mold))
      arr0(:,:) = arr2D(:,:)
      CALL REALLOCATE (arr2D,N+1,Mmax)
      arr2D(1:N,1:Mold) = arr0(:,:)
      IF (Mold < Mmax) arr2D(1:N,Mold+1:Mmax) = miss
      arr2D(N+1,1:M) = vec(:)
      IF (M < Mmax) arr2D(N+1,M+1:Mmax) = miss
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr2D_DP

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr2D_bool (arr2D,vec,missing)

    IMPLICIT NONE

    LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: arr2D
    LOGICAL, DIMENSION(:), INTENT(IN) :: vec
    LOGICAL, INTENT(IN), OPTIONAL :: missing

    INTEGER :: N, M, Mold, Mmax
    LOGICAL :: miss
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (PRESENT(missing)) THEN
      miss = missing
    ELSE
      miss = .False.
    END IF
    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr2D)) THEN
      ALLOCATE (arr2D(1,M))
      arr2D(1,:) = vec(:)
    ELSE
      N = SIZE(arr2D,1)
      Mold = SIZE(arr2D,2)
      Mmax = MAX(M,Mold)
      ALLOCATE (arr0(N,Mold))
      arr0(:,:) = arr2D(:,:)
      CALL REALLOCATE (arr2D,N+1,Mmax)
      arr2D(1:N,1:Mold) = arr0(:,:)
      IF (Mold < Mmax) arr2D(1:N,Mold+1:Mmax) = miss
      arr2D(N+1,1:M) = vec(:)
      IF (M < Mmax) arr2D(N+1,M+1:Mmax) = miss
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr2D_bool

    !-----------------------------------------------------------------------

  SUBROUTINE incrarr2D_char (arr2D,vec,missing)

    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: arr2D
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: vec
    CHARACTER(*), INTENT(IN), OPTIONAL :: missing

    INTEGER :: N, M, Mold, Mmax
    CHARACTER(LEN(arr2D)) :: miss
    CHARACTER(LEN(arr2D)), DIMENSION(:,:), ALLOCATABLE :: arr0    
    
    !-----------------------------------------------------------------------

    IF (PRESENT(missing)) THEN
      miss = missing
    ELSE
      miss = ""
    END IF
    M = SIZE(vec(:))
    IF (.NOT. ALLOCATED(arr2D)) THEN
      ALLOCATE (arr2D(1,M))
      arr2D(1,:) = vec(:)
    ELSE
      N = SIZE(arr2D,1)
      Mold = SIZE(arr2D,2)
      Mmax = MAX(M,Mold)
      ALLOCATE (arr0(N,Mold))
      arr0(:,:) = arr2D(:,:)
      CALL REALLOCATE (arr2D,N+1,Mmax)
      arr2D(1:N,1:Mold) = arr0(:,:)
      IF (Mold < Mmax) arr2D(1:N,Mold+1:Mmax) = miss
      arr2D(N+1,1:M) = vec(:)
      IF (M < Mmax) arr2D(N+1,M+1:Mmax) = miss
    END IF

    !-----------------------------------------------------------------------

  END SUBROUTINE incrarr2D_char


END MODULE arrays

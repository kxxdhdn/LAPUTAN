!******************************************************************************
!*
!*               GENERAL ADAPTATIVE GRID WITH ACCURACY CONTROL
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 04/2014 
  ! 
  ! 3) DESCRIPTION: function generating adaptative grid, with control of the 
  ! accuracy on the interpolation.
  !==========================================================================

MODULE adaptative_grid

  USE utilities, ONLY:
  USE arrays, ONLY:
  USE integration, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: gridadapt1D, grid_kernel

  INTERFACE gridadapt1D
    MODULE PROCEDURE gridadapt1D_simple, gridadapt1D_zgrid
  END INTERFACE gridadapt1D


CONTAINS


  ! Interface for DPRIMITIVE (removing the first element)
  !-------------------------
  PURE FUNCTION dPrim (x,f,xlog,ylog,lnx,lnf)

    USE utilities, ONLY: DP
    USE integration, ONLY: dPrimitive
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: x, f
    LOGICAL, INTENT(IN) :: xlog, ylog
    REAL(DP), INTENT(IN), DIMENSION(:), OPTIONAL :: lnx, lnf
    REAL(DP), DIMENSION(SIZE(x)-1) :: dPrim

    INTEGER :: N
    REAL(DP), DIMENSION(SIZE(x)) :: dF

    !------------------------------------------------------------------------

    N = SIZE(x)
    IF (.NOT. xlog .AND. ylog) THEN 
      dF(:) = DPRIMITIVE(x,f,XLOG=xlog,YLOG=ylog,LNF=lnf)
    ELSE IF (xlog .AND. ylog) THEN
      dF(:) = DPRIMITIVE(x,f,XLOG=xlog,YLOG=ylog,LNX=lnx,LNF=lnf)
    ELSE
      dF(:) = DPRIMITIVE(x,f,XLOG=xlog,YLOG=ylog)
    END IF
    dPrim(:) = dF(2:N)

    !------------------------------------------------------------------------  
 
  END FUNCTION dPrim


  !==========================================================================
  ! CALL GRIDADAPT1D_SIMPLE (xcoarse[N0],xfine[N],yfine[N],yfunc,accuracy, &
  !                          ABSACC=T/F,INTERP=T/F,INTEG=T/F,primitive[N], &
  !                          XLOG=T/F,YLOG=T/F,REEVALUATE=T/F,MINSTEP=, &
  !                          /LNFUNC,/RESCALE,/ADJUSTXLIM[2], &
  !                          SCALING=,NMAX=)
  !
  !   Returns an adaptative grid xfine[N], and function values yfine[N] of 
  ! function, YFUNC, from an initial grid xcoarse[N0] (between x0min and x0max).
  ! The grid adapted so that the accuracy on the linear interpretation between 
  ! two values is better than ACCURACY.
  !   If ABSACC is true, then it is an absolute accuracy (default is relative).
  !   If INTEG is true, the accuracy is on the integral and not on the 
  ! interpolation.
  !   If keyword REEVALUATE is set, then the function is evaluated on the entire
  ! grid at each iteration, while it is evaluated only at mid-points if it is 
  ! not set. REEVALUATE is used when the function needs the entire grid to be 
  ! evaluated properly (normalization, etc.).
  !   If keyword SLIM is set, then we keep only the points needed to achieve
  ! the required accuracy, and not necessary the last midpoints. This keyword
  ! is useful when one want use the computed grid later. This way the grid has
  ! the minimal size. 
  !  If /LNFUNC, then YFUNC is supposed to return ln(Y) instead of Y.
  !  If /RESCALE, then YFUNC is rescaled at each iteration so that its maximum 
  ! sampled is 1. In the end the results YFINE and PRIMITIVE are affected by 
  ! this rescaling. This scaling factor can be recovered through the SCALING
  ! argument. This function is useful when sampling a non normalized 
  ! probability distribution.
  !==========================================================================

  SUBROUTINE gridadapt1D_simple (xcoarse,xfine,yfine,yfunc,accuracy,Nmax,  &
                                 absacc,interp,integ,primitive,xlog,ylog, &
                                 reevaluate,minstep,slim,verbose,rescale, &
                                 lnfunc,adjustxlim,scaling)

    USE utilities, ONLY: DP, tinyDP, hugeDP, epsDP, trimlr, pring, isNaN, &
                         verbatim
    USE arrays, ONLY: reallocate, ramp, iwhere
    USE integration, ONLY: dPrimitive
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: xcoarse
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: xfine, yfine
    REAL(DP), INTENT(IN) :: accuracy
    INTEGER, INTENT(IN), OPTIONAL :: Nmax
    REAL(DP), INTENT(IN), OPTIONAL :: minstep
    LOGICAL, INTENT(IN), OPTIONAL :: absacc, interp, integ, xlog, ylog
    LOGICAL, INTENT(IN), OPTIONAL :: reevaluate, slim, verbose, rescale, lnfunc
    LOGICAL, DIMENSION(2), INTENT(IN), OPTIONAL :: adjustxlim
    REAL(DP), INTENT(OUT), OPTIONAL :: scaling
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:), OPTIONAL :: primitive
    INTERFACE
      FUNCTION yfunc (x)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: x
        REAL(DP), DIMENSION(SIZE(x)) :: yfunc
      END FUNCTION yfunc
    END INTERFACE

    INTEGER, PARAMETER :: Nrefmax = 50
    INTEGER, PARAMETER :: Nmax0 = HUGE(KIND(1))

    INTEGER :: i, k, N0, N1, Nref, Nmaxi, counter, Nslim, i1, i2, Nlim
    INTEGER, DIMENSION(:), ALLOCATABLE :: imid0, imid1, ind1
    REAL(DP), DIMENSION(2) :: xlim
    REAL(DP), DIMENSION(:), ALLOCATABLE :: x0, x1, f0, f1, xmid, fmid, fmidint
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnx0, lnx1, lnf0, lnf1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: lnxmid, lnfmid, lnfmidint
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dF0, dF1, prim1, threshint
    REAL(DP), DIMENSION(:), ALLOCATABLE :: deltadF0, deltadF1, relstep
    REAL(DP) :: threshprim, stepmin, scaling0, scaling1, ratscl, maxf
    LOGICAL :: integacc, relacc, xl, yl, evalall, dontkeep
    LOGICAL :: verb, rescl, lnfn
    LOGICAL, DIMENSION(2) :: adjxlim, limOK 
    LOGICAL, DIMENSION(:), ALLOCATABLE :: good0, good1, keep0, keep1
 
    !------------------------------------------------------------------------


    !------------------------------------------------------------------------
    !                             I. Preliminary
    !------------------------------------------------------------------------

    ! 1) Manage keywords and options
    !-------------------------------
    ! Method
    integacc = .False.
    IF (PRESENT(integ)) integacc = integ
    IF (PRESENT(interp)) integacc = (.NOT. interp)

    ! Accuracy
    relacc = .True.
    IF (PRESENT(absacc)) relacc = (.NOT. absacc)
   
    ! Maximum size if required
    IF (PRESENT(Nmax)) THEN ; Nmaxi = Nmax ; ELSE ; Nmaxi = Nmax0 ; END IF

    ! Minimum step
    stepmin = epsDP
    IF (PRESENT(minstep)) stepmin = minstep

    ! Verbose
    verb = verbatim
    IF (PRESENT(verbose)) verb = verbose

    ! Lin/log
    xl = .False.
    IF (PRESENT(xlog)) xl = xlog
    yl = .False.
    IF (PRESENT(ylog)) yl = ylog

    ! Function evaluation
    evalall = .False.
    IF (PRESENT(reevaluate)) evalall = reevaluate
    lnfn = .False.
    IF (PRESENT(lnfunc)) lnfn = lnfunc
    rescl = .False.
    IF (PRESENT(rescale)) rescl = rescale
    IF (lnfn) THEN
      scaling0 = 0._DP
      scaling1 = 0._DP
    ELSE
      scaling0 = 1._DP
      scaling1 = 1._DP
    END IF

    ! Optimal grid
    dontkeep = .True.
    IF (PRESENT(slim)) dontkeep = slim

    ! Refinement of the coarse grid
    adjxlim(:) = [.False.,.False.]
    IF (PRESENT(adjustxlim)) adjxlim(:) = adjustxlim(:)
    limOK(:) = ( .NOT. adjxlim(:) )
    xlim(:) = [ MINVAL(xcoarse(:)), MAXVAL(xcoarse(:)) ]


    ! 2) Design the coarse grid
    !--------------------------
    ! Initialize the coarse grid
    N0 = SIZE(xcoarse(:))
    ALLOCATE (x0(N0),f0(N0))
    IF (yl .OR. lnfn) ALLOCATE(lnf0(N0))
    x0(:) = xcoarse(:)
    IF (verb) PRINT*, "Adapatative grid: coarse grid; N = "//TRIMLR(PRING(N0))

    ! Push the limits away, if required
    counter = 0
    coarsegrid: DO 

      ! Swap
      N1 = N0
      CALL REALLOCATE(x1,N1)
      CALL REALLOCATE(f1,N1)
      IF (yl .OR. lnfn) CALL REALLOCATE(lnf1,N1)
      x1(:) = x0(:)
      scaling1 = scaling0
      
      ! Evalutate the function
      evalinit: IF (.NOT. evalall) THEN
        notallinit: IF (.NOT. lnfn) THEN
          IF (counter == 0) THEN
            f1(:) = YFUNC(x1(:))
          ELSE
            Nlim = COUNT(.NOT. limOK(:))
            IF (Nlim > 0) THEN
              CALL REALLOCATE(ind1,Nlim)
              ind1(:) = PACK( [1,N1], (.NOT. limOK(:)) )
              f1(ind1(:)) = YFUNC(x1(ind1(:))) / scaling1
            END IF
            f1(i1:i2) = f0(i1:i2)
          END IF
          IF (rescl) THEN
            maxf = MAXVAL(f1(:)) * scaling1
            scaling1 = MERGE( maxf, scaling0, &
                              ABS(maxf) < hugeDP .AND. ABS(maxf) > tinyDP )
            f1(:) = f1(:) * ( scaling0 / scaling1 )
          END IF
          IF (yl) THEN
            lnf1(:) = LOG(f1(:))
            WHERE (f1(:) < tinyDP) lnf1(:) = -hugeDP
          END IF
        ELSE
          IF (counter == 0) THEN
            lnf1(:) = YFUNC(x1(:))
          ELSE
            Nlim = COUNT(.NOT. limOK(:))
            IF (Nlim > 0) THEN
              CALL REALLOCATE(ind1,Nlim)
              ind1(:) = PACK( [1,N1], (.NOT. limOK(:)) )
              lnf1(ind1(:)) = YFUNC(x1(ind1(:))) - scaling1
            END IF
            lnf1(i1:i2) = lnf0(i1:i2) 
          END IF
          IF (rescl) THEN
            maxf = MAXVAL(lnf1(:)) + scaling0
            scaling1 = MERGE( maxf, scaling0, ABS(maxf) < hugeDP )
            lnf1(:) = lnf1(:) + ( scaling0 - scaling1 )
          END IF
          f1(:) = EXP(lnf1(:))
        END IF notallinit
      ELSE
        allinit: IF (.NOT. lnfn) THEN
          f1(:) = YFUNC(x1(:))
          IF (rescl) THEN
            maxf = MAXVAL(f1(:))
            scaling1 = MERGE( maxf, 1._DP, &
                              ABS(maxf) < hugeDP .AND. ABS(maxf) > tinyDP )
            f1(:) = f1(:) / scaling1
          END IF
          IF (yl) THEN
            lnf1(:) = LOG(f1(:))
            WHERE (f1(:) < tinyDP) lnf1(:) = -hugeDP
          END IF
        ELSE
          lnf1(:) = YFUNC(x1(:))
          IF (rescl) THEN
            maxf = MAXVAL(lnf1(:))
            scaling1 = MERGE( maxf, 0._DP, ABS(maxf) < hugeDP )
            lnf1(:) = lnf1(:) - scaling1
          END IF
          f1(:) = EXP(lnf1(:))
        END IF allinit
      END IF evalinit

      ! Iterate until we find zeros on each side. This method implicitly assumes
      ! that some non-zero values of the function are included in the original
      ! interval. The iteration looks for the edges (below machine precision).
      ! This method is not appropriate if the function has several narrow peaks 
      ! with zeros between them. 
      counter = counter + 1
      adjustlim: IF ( ALL(limOK(:)) .OR. counter == Nrefmax ) THEN
        EXIT
      ELSE
        IF (.NOT. limOK(1)) THEN
          IF (f1(1) > tinyDP .AND. f1(2) /= f1(1)) THEN
            xlim(1) = MERGE( x1(1)**2/x1(N1), 2._DP*x1(1)-x1(N1), xl )
          ELSE
            xlim(1) = x1(1)
            limOK(1) = .True.
          END IF
        END IF
        IF (.NOT. limOK(2)) THEN
          IF (f1(N1) > tinyDP .AND. f1(N1-1) /= f1(N1)) THEN
            xlim(2) = MERGE( x1(N1)**2/x1(1), 2._DP*x1(N1)-x1(1), xl )
          ELSE
            xlim(2) = x1(N1)
            limOK(2) = .True.
          END IF
        END IF
        N0 = N1 + COUNT(.NOT. limOK(:))
        IF (verb .AND. ANY(.NOT. limOK(:))) &
          PRINT*, "Adapatative grid: new limits "//TRIMLR(PRING(counter)) &
                  //"; N="//TRIMLR(PRING(N0))
        CALL REALLOCATE(x0,N0) 
        CALL REALLOCATE(f0,N0) 
        IF (yl .OR. lnfn) CALL REALLOCATE(lnf0,N0) 
        i1 = MERGE( 1, 2, limOK(1) )
        i2 = MERGE( N0, N0-1, limOK(2) )
        x0(i1:i2) = x1(:)
        f0(i1:i2) = f1(:)
        scaling0 = scaling1
        IF (yl .OR. lnfn) lnf0(i1:i2) = lnf1(:)
        IF (.NOT. limOK(1)) x0(1) = xlim(1)
        IF (.NOT. limOK(2)) x0(N0) = xlim(2)
      END IF adjustlim
  
    END DO coarsegrid

    ! Primitive on the coarse grid
    IF (xl) THEN
      ALLOCATE (lnx1(N1))
      lnx1(:) = LOG(x1(:))
    END IF
    IF (integacc) THEN
      ALLOCATE (dF1(N1))
      IF (.NOT. xl .AND. yl) THEN 
        dF1(:) = DPRIMITIVE(x1(:),f1(:),XLOG=xl,YLOG=yl,LNF=lnf1(:))
      ELSE IF (xl .AND. yl) THEN 
        dF1(:) = DPRIMITIVE(x1(:),f1(:),XLOG=xl,YLOG=yl, &
                            LNX=lnx1(:),LNF=lnf1(:))
      ELSE
        dF1(:) = DPRIMITIVE(x1(:),f1(:),XLOG=xl,YLOG=yl)
      END IF
      ALLOCATE(deltadF1(N1-1))
      deltadF1(:) = SUM(dF1(:)) ! large value for initialization
    END IF
    Nref = N1 - 1
    ALLOCATE (good1(N1-1))
    good1(:) = .False.
    IF (dontkeep) THEN
      ALLOCATE (keep1(N1))
      keep1(:) = .True.
    END IF

    
    !------------------------------------------------------------------------
    !                         II. Refining the grid
    !------------------------------------------------------------------------

    counter = 0
    refinement: DO
        

      ! 1) Preamble
      !------------
      ! Update the various arrays
      N0 = N1
      CALL REALLOCATE(x0,N0)
      CALL REALLOCATE(f0,N0)
      CALL REALLOCATE(good0,N0-1)
      IF (dontkeep) CALL REALLOCATE(keep0,N0)
      x0(:) = x1(:)
      f0(:) = f1(:)
      IF (xl) THEN
        CALL REALLOCATE(lnx0,N0)
        lnx0(:) = lnx1(:)
      END IF
      IF (yl .OR. lnfn) THEN
        CALL REALLOCATE(lnf0,N0)
        lnf0(:) = lnf1(:)
      END IF
      good0(:) = good1(:)
      IF (dontkeep) keep0(:) = keep1(:)
      IF (integacc) THEN 
        CALL REALLOCATE(dF0,N0)
        dF0(:) = dF1(:)
        CALL REALLOCATE(deltadF0,N0-1)
        deltadF0(:) = deltadF1(:)
      END IF
      IF (rescl) scaling0 = scaling1

      ! Indices where to refine      
      N1 = N0 + Nref
      IF (verb) &
        PRINT*, "Adapatative grid: iteration "//TRIMLR(PRING(counter+1)) &
                //"; N = "//TRIMLR(PRING(N1))
      CALL REALLOCATE(imid0,Nref) ! position of mid-points on grid 0 (x0 below)
      CALL REALLOCATE(imid1,Nref) ! position of mid-points on grid 1
      CALL REALLOCATE(ind1,N0)    ! position of grid 0 points on grid 1
      imid0(:) = PACK([(i,i=1,N0-1)],.NOT. good0(:))
      imid1(:) = imid0(:) + [(i,i=1,Nref)]
      k = 0
      ind1(1) = 1
      DO i=2,N0
        IF (.NOT. good0(i-1)) k = k + 1
        ind1(i) = i + k
      END DO

      
      ! 2) Refine the grid
      !-------------------
      ! Add mid-points where needed
      CALL REALLOCATE(xmid,Nref)
      IF (xl) THEN
        CALL REALLOCATE(lnxmid,Nref)
        lnxmid(:) = 0.5_DP * ( lnx0(imid0(:)+1) + lnx0(imid0(:)) )
        xmid(:) = EXP(lnxmid(:))
      ELSE 
        xmid(:) = 0.5_DP * ( x0(imid0(:)+1) + x0(imid0(:)) )      
      END IF

      ! Update the refinement of the grid
      CALL REALLOCATE(x1,N1)
      x1(ind1(:)) = x0(:)
      x1(imid1(:)) = xmid(:)
      IF (xl) THEN
        CALL REALLOCATE(lnx1,N1)
        lnx1(ind1(:)) = lnx0(:)
        lnx1(imid1(:)) = lnxmid(:)
      END IF

      ! Evaluate function on the new grid
      evalfunc: IF (.NOT. evalall) THEN
        CALL REALLOCATE(fmid,Nref)
        mideval: IF (.NOT. lnfn) THEN
          fmid(:) = YFUNC(xmid(:))
          midrescl: IF (rescl) THEN
            maxf = MAXVAL(fmid(:))
            scaling1 = MERGE( maxf, scaling0, &
                              ( ABS(maxf) < hugeDP .AND. ABS(maxf) > tinyDP ) )
            IF (scaling1 >= scaling0) THEN 
              ratscl = scaling0 / scaling1
              f0(:) = f0(:) * ratscl
              IF (yl) lnf0(:) = lnf0(:) + LOG(ratscl)
              IF (integacc) THEN
                dF0(:) = dF0(:) * ratscl
                deltadF0(:) = deltadF0(:) * ratscl
              END IF
            ELSE
              scaling1 = scaling0
            END IF
            fmid(:) = fmid(:) / scaling1
          END IF midrescl
          IF (yl) THEN
            CALL REALLOCATE(lnfmid,Nref)
            lnfmid(:) = LOG(fmid(:))
            WHERE (fmid(:) < tinyDP) lnfmid(:) = -hugeDP
          END IF
        ELSE
          CALL REALLOCATE(lnfmid,Nref)
          lnfmid(:) = YFUNC(xmid(:))
          midlnrescl: IF (rescl) THEN
            maxf = MAXVAL(lnfmid(:))
            scaling1 = MERGE( maxf, scaling0, ABS(maxf) < hugeDP )
            IF (scaling1 >= scaling0) THEN
              ratscl = EXP( scaling0 - scaling1 )
              f0(:) = f0(:) * ratscl
              lnf0(:) = lnf0(:) + ( scaling0 - scaling1 )
              IF (integacc) THEN
                dF0(:) = dF0(:) * ratscl
                deltadF0(:) = deltadF0(:) * ratscl            
              END IF
            ELSE
              scaling1 = scaling0
            END IF
          END IF midlnrescl
          lnfmid(:) = lnfmid(:) - scaling1
          fmid(:) = EXP(lnfmid(:))
        END IF mideval
        CALL REALLOCATE(f1,N1)
        f1(ind1(:)) = f0(:)
        f1(imid1(:)) = fmid(:)
        IF (yl .OR. lnfn) THEN
          CALL REALLOCATE(lnf1,N1)
          lnf1(ind1(:)) = lnf0(:)
          lnf1(imid1(:)) = lnfmid(:)
        END IF
      ELSE
        CALL REALLOCATE(f1,N1)
        CALL REALLOCATE(fmid,Nref)
        alleval: IF (.NOT. lnfn) THEN
          f1(:) = YFUNC(x1(:))
          allrescl: IF (rescl) THEN
            maxf = MAXVAL(f1(:))
            scaling1 = MERGE( maxf, 1._DP, &
                              ( ABS(maxf) < hugeDP .AND. ABS(maxf) > tinyDP ) )
            ratscl = scaling0 / scaling1
            f0(:) = f0(:) * ratscl
            IF (yl) lnf0(:) = lnf0(:) + LOG(ratscl)
            IF (integacc) THEN
              dF0(:) = dF0(:) * ratscl
              deltadF0(:) = deltadF0(:) * ratscl
            END IF
          END IF allrescl
          fmid(:) = f1(imid1(:))
          IF (yl) THEN
            CALL REALLOCATE(lnf1,N1)
            lnf1(:) = LOG(f1(:))
            WHERE (f1(:) < tinyDP) lnf1(:) = -hugeDP
            CALL REALLOCATE(lnfmid,Nref)
            lnfmid(:) = lnf1(imid1(:))
          END IF
        ELSE
          CALL REALLOCATE(lnf1,N1)
          lnf1(:) = YFUNC(x1(:))
          alllnrescl: IF (rescl) THEN
            maxf = MAXVAL(lnf1(:))
            scaling1 = MERGE( maxf, 0._DP, ABS(maxf) < hugeDP )
            ratscl = EXP( scaling0 - scaling1 )
            f0(:) = f0(:) * ratscl
            lnf0(:) = lnf0(:) + ( scaling0 - scaling1 )
            IF (integacc) THEN
              dF0(:) = dF0(:) * ratscl
              deltadF0(:) = deltadF0(:) * ratscl            
            END IF
          END IF alllnrescl
          CALL REALLOCATE(lnfmid,Nref)
          lnfmid(:) = lnf1(imid1(:))
          f1(:) = EXP(lnf1(:))
          fmid(:) = f1(imid1(:))
        END IF alleval
      END IF evalfunc
      
      ! Update the mask    
      CALL REALLOCATE(good1,N1-1)
      good1(ind1(1:N0-1)) = good0(:)


      ! 3) Interpolation method
      !------------------------
      interpol: IF (.NOT. integacc) THEN

        ! Interpolation at mid-point
        CALL REALLOCATE(fmidint,Nref)
        IF (yl) THEN
          CALL REALLOCATE(lnfmidint,Nref)
          lnfmidint(:) = lnf0(imid0(:)) &
                       + 0.5_DP * ( lnf0(imid0(:)+1) - lnf0(imid0(:)) )
          fmidint(:) = EXP(lnfmidint(:))
        ELSE
          fmidint(:) = f0(imid0(:)) + 0.5_DP * ( f0(imid0(:)+1) - f0(imid0(:)) )
        END IF

        ! Update the mask
        CALL REALLOCATE(threshint,Nref)
        IF (relacc) THEN
          WHERE( ABS(fmid(:)) > tinyDP/accuracy )
            threshint(:) = ABS(fmid(:)) * accuracy
          ELSEWHERE ( ABS(fmidint(:)) > tinyDP/accuracy )
            threshint(:) = ABS(fmidint(:)) * accuracy
          ELSEWHERE
            threshint(:) = ABS(fmid(:)-fmidint(:)) * 2._DP
          END WHERE
        ELSE
          threshint(:) = accuracy
        END IF
        FORALL (i=1:Nref) &
          good1(imid1(i)-1:imid1(i)) &
            = ( ABS(fmid(i) - fmidint(i)) <= threshint(i) .OR. isNaN(fmid(i)) &
                .OR. threshint(i) < tinyDP )
        IF (yl) THEN
          FORALL (i=1:Nref,fmid(i) < tinyDP) good1(imid1(i)-1:imid1(i)) = .True.
        END IF
      END IF interpol


      ! 4) Integration method
      !----------------------
      ! Compute the primitive
      integr: IF (integacc) THEN
        CALL REALLOCATE(dF1,N1)
        evalprim: IF (.NOT. evalall) THEN

          dF1(ind1(:)) = dF0(:)
          linlog: IF (.NOT. xl .AND. yl) THEN 
            FORALL (i=1:Nref) &
              dF1(imid1(i):imid1(i)+1) = &
                DPRIM( [x0(imid0(i)),xmid(i),x0(imid0(i)+1)], &
                       [f0(imid0(i)),fmid(i),f0(imid0(i)+1)], &
                       XLOG=xl,YLOG=(lnf0(imid0(i))/=lnf0(imid0(i)+1)), &
                       LNF=[lnf0(imid0(i)),lnfmid(i),lnf0(imid0(i)+1)] )
          ELSE IF (xl .AND. yl) THEN 
            FORALL (i=1:Nref) &
              dF1(imid1(i):imid1(i)+1) = &
                DPRIM( [x0(imid0(i)),xmid(i),x0(imid0(i)+1)], &
                       [f0(imid0(i)),fmid(i),f0(imid0(i)+1)], XLOG=xl, YLOG=yl,&
                       LNX=[lnx0(imid0(i)),lnxmid(i),lnx0(imid0(i)+1)], &
                       LNF=[lnf0(imid0(i)),lnfmid(i),lnf0(imid0(i)+1)] )
          ELSE
            FORALL (i=1:Nref) &
              dF1(imid1(i):imid1(i)+1) = &
                DPRIM( [x0(imid0(i)),xmid(i),x0(imid0(i)+1)], &
                       [f0(imid0(i)),fmid(i),f0(imid0(i)+1)], XLOG=xl, YLOG=yl )
          END IF linlog

        ELSE

          linlog2: IF (.NOT. xl .AND. yl) THEN 
            dF1(:) = DPRIMITIVE(x1(:),f1(:),XLOG=xl,YLOG=yl,LNF=lnf1(:))
          ELSE IF (xl .AND. yl) THEN 
            dF1(:) = DPRIMITIVE(x1(:),f1(:),XLOG=xl,YLOG=yl, &
                                LNX=lnx1(:),LNF=lnf1(:))
          ELSE
            dF1(:) = DPRIMITIVE(x1(:),f1(:),XLOG=xl,YLOG=yl)
          END IF linlog2
          
        END IF evalprim
        WHERE (ISNAN(dF1(:))) dF1(:) = 0._DP
        IF (yl) &
          WHERE (dF1(:) <= 0._DP) dF1(:) = 0._DP
        
        ! Decide where to refine
        CALL REALLOCATE(deltadF1,N1-1)
        deltadF1(ind1(1:N0-1)) = deltadF0(:)
        FORALL (i=1:Nref) &
          deltadF1(imid1(i)-1:imid1(i)) = ABS( SUM(dF1(imid1(i):imid1(i)+1)) &
                                             - dF0(imid0(i)+1) )
        threshprim = MERGE(accuracy * SUM(dF1(:)) / SQRT(REAL(N1,DP)), &
                           accuracy,relacc)
        good1(:) = ( deltadF1(:) <= threshprim )
        
      END IF integr


      ! 5) Adjustments and exit conditions
      !-----------------------------------
      ! Case of narrow support function
      IF (yl .OR. lnfn) THEN
        IF (COUNT(f1(:) > 0._DP) <= 1) good1(:) = .False.
      END IF

      ! Slim mode
      IF (dontkeep) THEN
        CALL REALLOCATE(keep1,N1)
        keep1(ind1(:)) = keep0(:)
        keep1(imid1(:)) = (.NOT. good1(imid1(:)) )
      END IF

      ! Enforce minimum step
      CALL REALLOCATE (relstep,Nref)
      relstep(:) = (x0(imid0(:)+1)-x0(imid0(:))) &
                 / MERGE(x0(imid0(:)),x0(imid0(:)+1),x0(imid0(:)) /= 0._DP)
      FORALL (i=1:Nref,ABS(relstep(i)) <= stepmin) &
        good1(imid1(i)-1:imid1(i)) = .True.

      ! Exit when no more refinement is needed
      Nref = COUNT(.NOT. good1)
      counter = counter + 1
      IF (counter == Nrefmax .OR. Nref == 0 .OR. N1 >= Nmaxi) EXIT

    END DO refinement


    !------------------------------------------------------------------------
    !                          III. Final grid
    !------------------------------------------------------------------------

    ! Slim mode
    nokeep: IF (.NOT. dontkeep) THEN
      ALLOCATE (xfine(N1),yfine(N1))
      xfine(:) = x1(:)
      yfine(:) = f1(:)
    ELSE
      Nslim = COUNT(keep1(:))
      ALLOCATE (xfine(Nslim),yfine(Nslim))
      xfine(:) = PACK(x1(:),keep1(:))
      yfine(:) = PACK(f1(:),keep1(:))
    END IF nokeep
    
    ! Primitive
    doprim: IF (integacc .AND. PRESENT(primitive)) THEN
      IF (.NOT. dontkeep) THEN
        ALLOCATE (primitive(N1))
        primitive(1) = 0._DP
        DO i=2,N1 ; primitive(i) = dF1(i) + primitive(i-1) ; END DO
      ELSE
        ALLOCATE (prim1(N1),primitive(Nslim))
        prim1(1) = 0._DP
        DO i=2,N1 ; prim1(i) = dF1(i) + prim1(i-1) ; END DO
        primitive(:) = PACK(prim1(:),keep1(:))
        DEALLOCATE (prim1)
      END IF
    END IF doprim

    ! Scaling
    IF (rescl .AND. PRESENT (scaling)) &
      scaling = MERGE(EXP(scaling1),scaling1,lnfn)

    ! Clean memory
    IF (verb) PRINT*, "Adaptative grid: packing and tagging..."
    DEALLOCATE(x0,x1,f0,f1,ind1,good0,good1,imid0,imid1,xmid,fmid,relstep)
    IF (yl .OR. lnfn) DEALLOCATE(lnf0,lnf1,lnfmid)
    IF (xl) DEALLOCATE(lnx0,lnx1,lnxmid)
    IF (integacc) DEALLOCATE(dF0,dF1,deltadF0,deltadF1)
    IF (dontkeep) DEALLOCATE(keep0,keep1)
    IF (.NOT. integacc) DEALLOCATE(fmidint)
    IF (.NOT. integacc .AND. yl) DEALLOCATE(lnfmidint)
    IF (.NOT. integacc) DEALLOCATE(threshint)

    !------------------------------------------------------------------------

  END SUBROUTINE gridadapt1D_simple


  !==========================================================================
  ! CALL GRIDADAPT1D_ZGRID (xcoarse[N0],zgrid[Nz],xfine[N],yfine[N,Nz],yfunc, &
  !                         accuracy,ABSACC=T/F,INTERP=T/F,INTEG=T/F, &
  !                         REEVALUATE=T/F,primitive[N],slim)
  !
  !   Returns an adaptative grid xfine[N], and function values yfine[N] of 
  ! function, YFUNC, from an initial grid xcoarse[N0] (between x0min and x0max).
  ! The grid adapted so that the accuracy on the linear interpretation between 
  ! two values is better than ACCURACY. The difference between the _SIMPLE and
  ! _ZGRID procedures, is that the interpolation of the data is simple in the 
  ! first case, and is done at one value of X, but several value of an 
  ! additional parameter Z, for _ZGRID. The latter case, can be the sampling
  ! in wavelength (X) of a black body spectrum, that has to be accurate for
  ! all temperatures (Z).
  !   IF ABSACC is true, then it is an absolute accuracy (default is relative).
  !   If INTEG is true, the accuracy is on the integral and not on the 
  ! interpolation.
  !   IF keyword REEVALUATE is set, then the function is evaluated on the entire
  ! grid at each iteration, while it is evaluated only at mid-points if it is 
  ! not set. REEVALUATE is used when the function needs the entire grid to be 
  ! evaluated properly (normalization, etc.).
  !   IF keyword SLIM is set, then we keep only the points needed to achieve
  ! the required accuracy, and not necessary the last midpoints. This keyword
  ! is useful when one want use the computed grid later. This way the grid has
  ! the minimal size. 
  !==========================================================================

  SUBROUTINE gridadapt1D_zgrid (xcoarse,zgrid,xfine,yfine,yfunc,accuracy, &
                                absacc,interp,integ,primitive,xlog,ylog, &
                                reevaluate,minstep,slim,verbose)

    USE utilities, ONLY: DP, tinyDP, hugeDP, epsDP, trimlr, pring, isNaN, &
                         verbatim
    USE arrays, ONLY: reallocate
    USE integration, ONLY: dPrimitive
    IMPLICIT NONE

    REAL(DP), INTENT(IN), DIMENSION(:) :: xcoarse, zgrid
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: xfine
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: yfine
    REAL(DP), INTENT(IN) :: accuracy
    REAL(DP), INTENT(IN), OPTIONAL :: minstep
    LOGICAL, INTENT(IN), OPTIONAL :: absacc, interp, integ, xlog, ylog
    LOGICAL, INTENT(IN), OPTIONAL :: reevaluate, slim, verbose
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:), OPTIONAL :: primitive
    INTERFACE
      FUNCTION yfunc (x,z)
        USE utilities, ONLY: DP
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: x, z
        REAL(DP), DIMENSION(SIZE(x),SIZE(z)) :: yfunc
      END FUNCTION yfunc
    END INTERFACE

    INTEGER, PARAMETER :: Nrefmax = 50

    INTEGER :: i, k, iz, N0, N1, Nz, Nref, counter, Nslim
    INTEGER, DIMENSION(:), ALLOCATABLE :: imid0, imid1, ind1
    REAL(DP), DIMENSION(SIZE(zgrid)) :: threshprim
    REAL(DP), DIMENSION(:), ALLOCATABLE :: x0, x1, xmid, lnx0, lnx1, lnxmid
    REAL(DP), DIMENSION(:), ALLOCATABLE :: relstep
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lnf0, lnf1, f0, f1, fmid, fmidint
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lnfmid, lnfmidint, threshint
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dF0, dF1, prim1
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: deltadF0, deltadF1
    REAL(DP) :: stepmin
    LOGICAL :: integacc, relacc, xl, yl, evalall, dontkeep, verb
    LOGICAL, DIMENSION(:), ALLOCATABLE :: good0, good1, keep0, keep1
 
    !------------------------------------------------------------------------

    ! 1. Preliminary
    !---------------
    ! Method
    integacc = .False.
    IF (PRESENT(integ)) integacc = integ
    IF (PRESENT(interp)) integacc = (.NOT. interp)

    ! Accuracy
    relacc = .True.
    IF (PRESENT(absacc)) relacc = (.NOT. absacc)

    ! Minimum step
    stepmin = epsDP
    IF (PRESENT(minstep)) stepmin = minstep

    ! Verbose
    verb = verbatim
    IF (PRESENT(verbose)) verb = verbose

    ! Lin/log
    xl = .False.
    IF (PRESENT(xlog)) xl = xlog
    yl = .False.
    IF (PRESENT(ylog)) yl = ylog

    ! Function evaluation
    evalall = .False.
    IF (PRESENT(reevaluate)) evalall = reevaluate

    ! Optimal grid
    dontkeep = .True.
    IF (PRESENT(slim)) dontkeep = slim

    ! Coarse grid
    N1 = SIZE(xcoarse)
    IF (verb) PRINT*, "Adapatative grid: coarse grid; N = "//TRIMLR(PRING(N1))
    Nz = SIZE(zgrid)
    ALLOCATE (x1(N1),f1(N1,Nz))
    x1(:) = xcoarse(:)
    f1(:,:) = YFUNC(x1(:),zgrid(:))
    IF (xl) THEN
      ALLOCATE (lnx1(N1))
      lnx1(:) = LOG(x1(:))
    END IF
    IF (yl) THEN
      ALLOCATE (lnf1(N1,Nz))
      lnf1(:,:) = LOG(f1(:,:))
      WHERE (f1(:,:) < tinyDP) lnf1(:,:) = -hugeDP
    END IF
    IF (integacc) THEN
      ALLOCATE (dF1(N1,Nz))
      IF (.NOT. xl .AND. yl) THEN 
        FORALL(iz=1:Nz) &
          dF1(:,iz) = dPRIMITIVE(x1,f1(:,iz),XLOG=xl,YLOG=yl,LNF=lnf1(:,iz))
      ELSE IF (xl .AND. yl) THEN 
        FORALL(iz=1:Nz) &
          dF1(:,iz) = DPRIMITIVE(x1(:),f1(:,iz),XLOG=xl,YLOG=yl, &
                            LNX=lnx1(:),LNF=lnf1(:,iz))
      ELSE
        FORALL(iz=1:Nz) &
          dF1(:,iz) = DPRIMITIVE(x1(:),f1(:,iz),XLOG=xl,YLOG=yl)
      END IF
      ALLOCATE(deltadF1(N1-1,Nz))
      FORALL (iz=1:Nz) deltadF1(:,iz) = SUM(dF1(:,iz)) ! large value for start
    END IF
    Nref = N1 - 1
    ALLOCATE (good1(N1-1))
    good1(:) = .False.
    IF (dontkeep) THEN
      ALLOCATE (keep1(N1))
      keep1(:) = .True.
    END IF

    
    ! 2. Refining the grid
    !---------------------
    counter = 0
    refinement: DO

      ! Preamble
      N0 = N1
      CALL REALLOCATE(x0,N0)
      CALL REALLOCATE(f0,N0,Nz)
      CALL REALLOCATE(good0,N0-1)
      IF (dontkeep) CALL REALLOCATE(keep0,N0)
      x0(:) = x1(:)
      f0(:,:) = f1(:,:)
      IF (xl) THEN
        CALL REALLOCATE(lnx0,N0)
        lnx0(:) = lnx1(:)
      END IF
      IF (yl) THEN
        CALL REALLOCATE(lnf0,N0,Nz)
        lnf0(:,:) = lnf1(:,:)
      END IF
      good0(:) = good1(:)
      IF (dontkeep) keep0(:) = keep1(:)
      IF (integacc) THEN 
        CALL REALLOCATE(dF0,N0,Nz)
        dF0(:,:) = dF1(:,:)
        CALL REALLOCATE(deltadF0,N0-1,Nz)
        deltadF0(:,:) = deltadF1(:,:)
      END IF

      ! Indices where to refine      
      N1 = N0 + Nref
      IF (verb) &
        PRINT*, "Adapatative grid: iteration "//TRIMLR(PRING(counter+1)) &
                //"; N = "//TRIMLR(PRING(N1))
      CALL REALLOCATE(imid0,Nref) ! position of mid-points on grid 0 (x0 below)
      CALL REALLOCATE(imid1,Nref) ! position of mid-points on grid 1
      CALL REALLOCATE(ind1,N0)    ! position of grid 0 points on grid 1
      imid0(:) = PACK([(i,i=1,N0-1)],.NOT. good0)
      imid1(:) = imid0(:) + [(i,i=1,Nref)]
      k = 0
      ind1(1) = 1
      DO i=2,N0
        IF (.NOT. good0(i-1)) k = k + 1
        ind1(i) = i + k
      END DO

      ! Add mid-points where needed
      CALL REALLOCATE(xmid,Nref)
      IF (xl) THEN
        CALL REALLOCATE(lnxmid,Nref)
        lnxmid(:) = 0.5_DP * ( lnx0(imid0(:)+1) + lnx0(imid0(:)) )
        xmid(:) = EXP(lnxmid(:))
      ELSE 
        xmid(:) = 0.5_DP * ( x0(imid0(:)+1) + x0(imid0(:)) )      
      END IF

      ! Update the refinement
      CALL REALLOCATE(x1,N1)
      x1(ind1(:)) = x0(:)
      x1(imid1(:)) = xmid(:)
      IF (xl) THEN
        CALL REALLOCATE(lnx1,N1)
        lnx1(ind1(:)) = lnx0(:)
        lnx1(imid1(:)) = lnxmid(:)
      END IF

      ! Evaluate the function on the new grid
      evalfunc: IF (.NOT. evalall) THEN
        CALL REALLOCATE(fmid,Nref,Nz)
        fmid(:,:) = YFUNC(xmid(:),zgrid(:))
        IF (yl) THEN
          CALL REALLOCATE(lnfmid,Nref,Nz)
          lnfmid(:,:) = LOG(fmid(:,:))
          WHERE (fmid(:,:) < tinyDP) lnfmid(:,:) = -hugeDP
        END IF
        CALL REALLOCATE(f1,N1,Nz)
        f1(ind1(:),:) = f0(:,:)
        f1(imid1(:),:) = fmid(:,:)
        IF (yl) THEN
          CALL REALLOCATE(lnf1,N1,Nz)
          lnf1(ind1(:),:) = lnf0(:,:)
          lnf1(imid1(:),:) = lnfmid(:,:)
        END IF
      ELSE 
        CALL REALLOCATE(f1,N1,Nz)
        CALL REALLOCATE(fmid,Nref,Nz)
        f1(:,:) = YFUNC(x1(:),zgrid(:))
        fmid(:,:) = f1(imid1(:),:)
        IF (yl) THEN
          CALL REALLOCATE(lnf1,N1,Nz)
          lnf1(:,:) = LOG(f1(:,:))
          WHERE (f1(:,:) < tinyDP) lnf1(:,:) = -hugeDP
          CALL REALLOCATE(lnfmid,Nref,Nz)
          lnfmid(:,:) = lnf1(imid1(:),:)
        END IF
      END IF evalfunc

      ! Interpolation at mid-point
      CALL REALLOCATE(fmidint,Nref,Nz)
      IF (yl) THEN
        CALL REALLOCATE(lnfmidint,Nref,Nz)
        lnfmidint(:,:) = lnf0(imid0(:),:) &
                       + 0.5_DP * ( lnf0(imid0(:)+1,:) - lnf0(imid0(:),:) )
        fmidint(:,:) = EXP(lnfmidint(:,:))
      ELSE
        fmidint(:,:) = f0(imid0(:),:) &
                     + 0.5_DP * ( f0(imid0(:)+1,:) - f0(imid0(:),:) )
      END IF

      ! Update the mask
      CALL REALLOCATE(good1,N1-1)
      good1(ind1(1:N0-1)) = good0(:)
      interpol: IF (.NOT. integacc) THEN 
        CALL REALLOCATE(threshint,Nref,Nz)
        IF (relacc) THEN
          FORALL(iz=1:Nz)
            WHERE( ABS(fmid(:,iz)) > tinyDP/accuracy )
              threshint(:,iz) = ABS(fmid(:,iz)) * accuracy
            ELSEWHERE ( ABS(fmidint(:,iz)) > tinyDP/accuracy )
              threshint(:,iz) = ABS(fmidint(:,iz)) * accuracy
            ELSEWHERE
              threshint(:,iz) = ABS(fmid(:,iz)-fmidint(:,iz)) * 2._DP
            END WHERE
          END FORALL
        ELSE
          threshint(:,:) = accuracy
        END IF
        FORALL (i=1:Nref) &
          good1(imid1(i)-1:imid1(i)) &
            = ALL( ABS(fmid(i,:) - fmidint(i,:)) <= threshint(i,:) &
                   .OR. isNaN(fmid(i,:)) .OR. threshint(i,:) < tinyDP)
        IF (yl) THEN
          FORALL (i=1:Nref,ALL(fmid(i,:) < tinyDP)) &
            good1(imid1(i)-1:imid1(i)) = .True.
        END IF
      END IF interpol

      ! Compute the primitive
      integr: IF (integacc) THEN
        CALL REALLOCATE(dF1,N1,Nz)
        evalprim: IF (.NOT. evalall) THEN

          dF1(ind1(:),:) = dF0(:,:)
          linlog: IF (.NOT. xl .AND. yl) THEN 
            FORALL (i=1:Nref,iz=1:Nz) &
              dF1(imid1(i):imid1(i)+1,iz) = &
                DPRIM( [x0(imid0(i)),xmid(i),x0(imid0(i)+1)], &
                       [f0(imid0(i),iz),fmid(i,iz),f0(imid0(i)+1,iz)], &
                       XLOG=xl, YLOG=yl, &
                      LNF=[lnf0(imid0(i),iz),lnfmid(i,iz),lnf0(imid0(i)+1,iz)] )
          ELSE IF (xl .AND. yl) THEN 
            FORALL (i=1:Nref,iz=1:Nz) &
              dF1(imid1(i):imid1(i)+1,iz) = &
                DPRIM( [x0(imid0(i)),xmid(i),x0(imid0(i)+1)], &
                       [f0(imid0(i),iz),fmid(i,iz),f0(imid0(i)+1,iz)], &
                       XLOG=xl, YLOG=yl, &
                       LNX=[lnx0(imid0(i)),lnxmid(i),lnx0(imid0(i)+1)], &
                      LNF=[lnf0(imid0(i),iz),lnfmid(i,iz),lnf0(imid0(i)+1,iz)] )
          ELSE
            FORALL (i=1:Nref,iz=1:Nz) &
              dF1(imid1(i):imid1(i)+1,iz) = &
                DPRIM( [x0(imid0(i)),xmid(i),x0(imid0(i)+1)], &
                       [f0(imid0(i),iz),fmid(i,iz),f0(imid0(i)+1,iz)], &
                       XLOG=xl, YLOG=yl )
          END IF linlog

        ELSE

          linlog2: IF (.NOT. xl .AND. yl) THEN 
            FORALL (iz=1:Nz) &
              dF1(:,iz) = DPRIMITIVE(x1(:),f1(:,iz),XLOG=xl,YLOG=yl, &
                                     LNF=lnf1(:,iz))
          ELSE IF (xl .AND. yl) THEN 
            FORALL (iz=1:Nz) &
              dF1(:,iz) = DPRIMITIVE(x1(:),f1(:,iz),XLOG=xl,YLOG=yl, &
                                     LNX=lnx1(:),LNF=lnf1(:,iz))
          ELSE
            FORALL (iz=1:Nz) &
              dF1(:,iz) = DPRIMITIVE(x1(:),f1(:,iz),XLOG=xl,YLOG=yl)
          END IF linlog2

        END IF evalprim
        WHERE(ISNAN(dF1(:,:))) dF1(:,:) = 0._DP

        ! Decide where to refine
        CALL REALLOCATE(deltadF1,N1-1,Nz)
        deltadF1(ind1(1:N0-1),:) = deltadF0(:,:)
        FORALL (i=1:Nref,iz=1:Nz) &
          deltadF1(imid1(i)-1:imid1(i),iz) &
            = ABS( SUM(dF1(imid1(i):imid1(i)+1,iz)) - dF0(imid0(i)+1,iz) )
        threshprim(:) = MERGE(accuracy*SUM(ABS(dF1(:,:)),1)/SQRT(REAL(N1,DP)), &
                              accuracy,relacc)
        FORALL (i=1:N1-1) good1(i) = ALL( deltadF1(i,:) <= threshprim(:) )
        IF (yl) THEN
          FORALL (i=1:Nref,ALL(fmid(i,:) < tinyDP)) &
            good1(imid1(i)-1:imid1(i)) = .True.
        END IF

      END IF integr

      ! Slim mode
      IF (dontkeep) THEN
        CALL REALLOCATE(keep1,N1)
        keep1(ind1(:)) = keep0(:)
        keep1(imid1(:)) = (.NOT. good1(imid1(:)) )
      END IF

      ! Enforce minimum step
      CALL REALLOCATE (relstep,Nref)
      relstep(:) = (x0(imid0(:)+1)-x0(imid0(:))) &
                 / MERGE(x0(imid0(:)),x0(imid0(:)+1),x0(imid0(:)) /= 0._DP)
      FORALL (i=1:Nref,ABS(relstep(i)) <= stepmin) &
        good1(imid1(i)-1:imid1(i)) = .True.

      ! Exit when no more refinement is needed
      Nref = COUNT(.NOT. good1)
      counter = counter + 1
      IF (counter == Nrefmax .OR. Nref == 0) EXIT

    END DO refinement


    ! 3. Fine grid
    !-------------
    nokeep: IF (.NOT. dontkeep) THEN
      ALLOCATE (xfine(N1),yfine(N1,Nz))
      xfine(:) = x1(:)
      yfine(:,:) = f1(:,:)
    ELSE
      Nslim = COUNT(keep1(:))
      ALLOCATE (xfine(Nslim),yfine(Nslim,Nz))
      xfine(:) = PACK(x1(:),keep1(:))
      FORALL (iz=1:Nz) yfine(:,iz) = PACK(f1(:,iz),keep1(:))
    END IF nokeep

    ! Primitive
    doprim: IF (integacc .AND. PRESENT(primitive)) THEN
      IF (.NOT. dontkeep) THEN
        ALLOCATE (primitive(N1,Nz))
        primitive(1,:) = 0._DP
        DO i=2,N1 ; primitive(i,:) = dF1(i,:) + primitive(i-1,:) ; END DO
      ELSE
        ALLOCATE (prim1(N1,Nz),primitive(Nslim,Nz))
        prim1(1,:) = 0._DP
        DO i=2,N1 ; prim1(i,:) = dF1(i,:) + prim1(i-1,:) ; END DO
        FORALL (iz=1:Nz) primitive(:,iz) = PACK(prim1(:,iz),keep1(:))
        DEALLOCATE (prim1)
      END IF
    END IF doprim

    ! Clean memory
    DEALLOCATE(x0,x1,f0,f1,ind1,good0,good1,imid0,imid1,xmid,fmid,relstep)
    IF (yl) DEALLOCATE(lnf0,lnf1,lnfmid)
    IF (xl) DEALLOCATE(lnx0,lnx1,lnxmid)
    IF (integacc) DEALLOCATE(dF0,dF1,deltadF0,deltadF1)
    IF (dontkeep) DEALLOCATE(keep0,keep1)
    IF (.NOT. integacc) DEALLOCATE(fmidint)
    IF (.NOT. integacc .AND. yl) DEALLOCATE(lnfmidint)
    IF (.NOT. integacc) DEALLOCATE(threshint)

    !------------------------------------------------------------------------

  END SUBROUTINE gridadapt1D_zgrid


  !===========================================================================
  !  x[N] = GRID_KERNEL(N,xinf,xsup,xgrid,fgrid,xlog,ylog)
  !
  !  Creates a grid by inversing a cumulative distribution function of f, 
  ! which is a tabulated function.
  !===========================================================================

  FUNCTION grid_kernel (N,xinf,xsup,xgrid,fgrid,xlog,ylog)
! TO BE TESTED!!!...

    USE utilities, ONLY: DP
    USE interpolation, ONLY: interp_lin_sorted
    USE integration, ONLY: integ_tabulated
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: xinf, xsup
    REAL(DP), DIMENSION(:), INTENT(IN) :: xgrid, fgrid
    REAL(DP), DIMENSION(N) :: grid_kernel
    LOGICAL, INTENT(IN), OPTIONAL :: xlog, ylog

    INTEGER :: i, Ngrid
    REAL(DP) :: intgrid
    REAL(DP), DIMENSION(N) :: ind
    REAL(DP), DIMENSION(SIZE(xgrid)) :: primgrid
    LOGICAL :: xl, yl

    !------------------------------------------------------------------------

    ! Options
    xl = .False.
    IF (PRESENT(xlog)) xl = xlog
    yl = .False.
    IF (PRESENT(ylog)) yl = ylog

    ! Primitive of the kernel
    Ngrid = SIZE(xgrid(:))
    intgrid = INTEG_TABULATED(xgrid(:),fgrid(:),XRANGE=[xinf,xsup], &
                              PRIMITIVE=primgrid,XLOG=xl,YLOG=yl)
    primgrid(:) = primgrid(:) / primgrid(Ngrid)
    
    ! Inversion
    ind(:) = [(i,i=0,N-1)] / REAL(N-1,DP)
    grid_kernel(:) = INTERP_LIN_SORTED(xgrid(:),primgrid(:),ind(:),XLOG=yl)

    !------------------------------------------------------------------------

  END FUNCTION grid_kernel


END MODULE adaptative_grid

!******************************************************************************
!*
!*                            SYNTHETIC PHOTMETRY
!*
!******************************************************************************


  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) DESCRIPTION: Least-squares near-IR-to-radio SED fits.
  !
  ! 3) HISTORY: 
  !    - 10/2016: created.
  !    - 02/2017: add calibration uncertainties.
  !==========================================================================

! TO DO:
!  - mask
!  - interpolation method
!  - uncertainty propagation


PROGRAM synthetic_photometry

  USE utilities, ONLY: DP
  USE inout, ONLY: h5ext, read_hdf5, write_hdf5
  USE instrument_filters, ONLY: lenfilter, read_filters, filter_type, &
                                calib_type, read_caliberr, wcen_filter, &
                                color_correction => synthetic_photometry
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: filIN = "./synthetic_photometry_input"//h5ext
  CHARACTER(*), PARAMETER :: filOUT = "./synthetic_photometry_output"//h5ext

  INTEGER :: Nx, Ny, Nfilt, x, y
  INTEGER, DIMENSION(:), ALLOCATABLE :: flags
  REAL(DP), DIMENSION(:), ALLOCATABLE :: w_spec, wcen
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Fnu_spec, Fnu_filt
  CHARACTER(lenfilter), DIMENSION(:), ALLOCATABLE :: filt
  LOGICAL :: docalib, dophot
  TYPE(filter_type) :: filters
  TYPE(calib_type) :: calibstr
  
  !---------------------------------------------------------------------------

  ! Read the input file
  CALL READ_HDF5(INTARR1D=flags,FILE=filIN,NAME="(docalib,dophot)")
  docalib = ( flags(1) == 1 )
  dophot = ( flags(2) == 1 )
  CALL READ_HDF5(STRARR1D=filt,FILE=filIN,NAME="Filter label",N1=Nfilt)
  IF (dophot) THEN
    CALL READ_HDF5(DBLARR1D=w_spec,FILE=filIN,NAME="Wavelength (microns)")  
    CALL READ_HDF5(DBLARR3D=Fnu_spec,FILE=filIN,NAME="Flux (x.Hz-1)", &
                    N1=Nx,N2=Ny)
  END IF

  ! Read the filters and the calibration matrices
  CALL READ_FILTERS(filters,filt(:))
  IF (docalib) CALL READ_CALIBERR(calibstr,filt(:))

  ! Compute the synthetic photometry of each pixel and each filter
  IF (dophot) THEN
    ALLOCATE (Fnu_filt(Nx,Ny,Nfilt),wcen(Nfilt))
    DO x=1,Nx
      DO y=1,Ny
        Fnu_filt(x,y,:) = COLOR_CORRECTION(w_spec(:),Fnu_spec(x,y,:), &
                                           filt(:),filters,wcen(:))
      END DO
    END DO
  ELSE
    ALLOCATE (wcen(Nfilt))
    wcen(:) = WCEN_FILTER(filt(:))
  END IF
  
  ! Write the output file
  IF (dophot) &
    CALL WRITE_HDF5(DBLARR3D=Fnu_filt(:,:,:),FILE=filOUT,NAME="Flux (x.Hz-1)",&
                    VERBOSE=.False.)
  CALL WRITE_HDF5(STRARR1D=filt(:),FILE=filOUT,NAME="Filter label", &
                  VERBOSE=.False.,APPEND=dophot)
  CALL WRITE_HDF5(DBLARR1D=wcen(:),FILE=filOUT,VERBOSE=.False., &
                  NAME="Central wavelength (microns)",APPEND=.True.)
  IF (docalib) THEN
    CALL WRITE_HDF5(DBLARR2D=calibstr%matcov(:,:),FILE=filOUT, &
                    NAME="Covariance matrix",VERBOSE=.False.,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR2D=calibstr%matR(:,:),FILE=filOUT, &
                    NAME="Correlation matrix",VERBOSE=.False.,APPEND=.True.)
    CALL WRITE_HDF5(DBLARR2D=calibstr%matS(:,:),FILE=filOUT, &
                    NAME="Standard deviation matrix",VERBOSE=.False., &
                    APPEND=.True.)
  END IF
  
  !---------------------------------------------------------------------------
  
END PROGRAM synthetic_photometry

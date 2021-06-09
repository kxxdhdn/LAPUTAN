!******************************************************************************
!*
!*                  EVERYTHING ABOUT ASTRONOMICAL FILTERS
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 09/2013.
  !    - 11/2014: adaptative grid for the integration.
  !    - 02/2015: update with HDF5 files.
  !    - 04/2016: add to WCEN_FILTER the possibility of having individual
  !      wavelenght elements.
  !    - 05/2016: add to READ_CALIBERR the possibility of having individual
  !      wavelength elements.
  !    - 07/2016: add DIRBE filters.
  !    - 10/2018: add LFI filters.
  ! 
  ! 3) DESCRIPTION: Defines the various filter quantities and color 
  !                 corrections.
  !==========================================================================


MODULE instrument_filters_ext

  USE utilities, ONLY: DP
  USE constants, ONLY:
  USE interpolation, ONLY:
  USE integration, ONLY:
  IMPLICIT NONE
  PRIVATE

  REAL(DP), SAVE, PUBLIC :: a1, a2, nu0
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wfilt, transfilt
  REAL(DP), DIMENSION(:), ALLOCATABLE, SAVE, PUBLIC :: wSED, FnuSED

  PUBLIC :: colcorr


CONTAINS


  FUNCTION colcorr (wad)

    USE utilities, ONLY: DP, isNaN
    USE constants, ONLY: MKS
    USE interpolation, ONLY: interp_lin_sorted
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: wad
    REAL(DP), DIMENSION(SIZE(wad)) :: colcorr

    REAL(DP), DIMENSION(SIZE(wad)) :: transad, Fnuad, nuad, nu0ad, nuwad

    !----------------------------------------------------------------------

    ! Fequency grid
    nuad(:) = MKS%clight / MKS%micron / wad(:)
    nuwad(:) = nuad(:) / wad(:)

    ! Interpolate the filter transmission on the adaptative grid
    nu0ad(:) = nu0 / nuad(:)
    transad(:) = INTERP_LIN_SORTED(transfilt(:),wfilt(:),wad(:), &
                                   XLOG=.False.,YLOG=.False.)
  
    ! Interpolate the SED on the adaptative grid (the templates are designed 
    ! to be interpolated in log-log)
    Fnuad(:) = INTERP_LIN_SORTED(FnusED(:),wSED(:),wad(:), &
                                 XLOG=.True.,YLOG=.True.)

    ! Output
    colcorr(:) = Fnuad(:) * nu0ad(:)**a1 * transad(:) * nuwad(:)

    !----------------------------------------------------------------------

  END FUNCTION colcorr

END MODULE instrument_filters_ext


  !==========================================================================


MODULE instrument_filters

  USE utilities, ONLY: DP
  USE constants, ONLY: 
  USE inout, ONLY:
  USE integration, ONLY: integ_tabulated
  USE adaptative_grid, ONLY:
  USE instrument_filters_ext, ONLY:
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: wcen_filter, read_filters, read_caliberr, synthetic_photometry

  INTERFACE wcen_filter
    MODULE PROCEDURE wcen_filter_scl, wcen_filter_1D
  END INTERFACE wcen_filter

  ! Length of filter string labels
  INTEGER, PARAMETER, PUBLIC :: lenfilter = 20

  ! Complete list
  INTEGER, PARAMETER, PUBLIC :: Nphotall = 61
  REAL(DP), DIMENSION(Nphotall), PARAMETER, PUBLIC :: &
    wcenall = [ 1.235_DP,   1.662_DP,   2.159_DP, &                 ! 2MASS
                 2.40_DP,    3.20_DP,    4.10_DP, &                 ! IRC N
                  7.0_DP,     9.0_DP,    11.0_DP, &                 ! IRC S
                 15.0_DP,    18.0_DP,    24.0_DP, &                 ! IRC L
                 65.0_DP,     90._DP,    140._DP,    160._DP, &     ! FIS
               3.5443_DP,  4.4870_DP,  5.7104_DP,  7.8413_DP, &     ! IRAC
               23.675_DP,   71.44_DP, 155.899_DP, &                 ! MIPS
               11.988_DP,  24.975_DP,   59.94_DP,   99.90_DP, &     ! IRAS
               3.3526_DP,  4.6028_DP, 11.5608_DP, 22.0883_DP, &     ! WISE
               22.800_DP, &                                         ! WISE old
                 8.28_DP,   12.13_DP,   14.65_DP,   21.34_DP, &     ! MSXA,C,D,E
                 1.25_DP,     2.2_DP,     3.5_DP,     4.9_DP, &     ! DIRBE1-4
                  12._DP,     25._DP,     60._DP,    100._DP, &     ! DIRBE5-8
                 140._DP,    240._DP, &                             ! DIRBE9-10
                  70._DP,    100._DP,    160._DP, &                 ! PACS
                 250._DP,    350._DP,    500._DP, &                 ! SPIRE
              349.708_DP, 549.908_DP, 849.008_DP, 1381.11_DP, &     ! HFI
              2095.80_DP, 2997.00_DP, &                             ! HFI
              4257.10_DP, 6795.92_DP, 10552.8_DP ]                  ! LFI
  CHARACTER(*), DIMENSION(Nphotall), PARAMETER, PUBLIC :: &
    filtall = [ "2MASS1      ", "2MASS2      ", "2MASS3      ", &
                "AKARI_IRC1  ", "AKARI_IRC2  ", "AKARI_IRC3  ", "AKARI_IRC4  ",&
                "AKARI_IRC5  ", "AKARI_IRC6  ", "AKARI_IRC7  ", "AKARI_IRC8  ",&
                "AKARI_IRC9  ", &
                "AKARI_FIS1  ", "AKARI_FIS2  ", "AKARI_FIS3  ", "AKARI_FIS4  ",&
                "IRAC1       ", "IRAC2       ", "IRAC3       ", "IRAC4       ",&
                "MIPS1       ", "MIPS2       ", "MIPS3       ", &
                "IRAS1       ", "IRAS2       ", "IRAS3       ", "IRAS4       ",&
                "WISE1       ", "WISE2       ", "WISE3       ", "WISE4       ",&
                "WISE4_old   ", &
                "MSX1        ", "MSX2        ", "MSX3        ", "MSX4        ",&
                "DIRBE1      ", "DIRBE2      ", "DIRBE3      ", "DIRBE4      ",&
                "DIRBE5      ", "DIRBE6      ", "DIRBE7      ", "DIRBE8      ",&
                "DIRBE9      ", "DIRBE10     ", &
                "PACS1       ", "PACS2       ", "PACS3       ", &
                "SPIRE1      ", "SPIRE2      ", "SPIRE3      ", &
                "HFI1        ", "HFI2        ", "HFI3        ", "HFI4        ",&
                "HFI5        ", "HFI6        ", &
                "LFI1        ", "LFI2        ", "LFI3        " ]
  INTEGER, PARAMETER, PUBLIC :: Nspecall = 7
  CHARACTER(*), DIMENSION(Nspecall), PARAMETER, PUBLIC :: &
    specall = [ "IRS_SL ", "IRS_LL ", "IRS_SH ", "IRS_LH ", &
                "IRC_NG ", "IRC_SG1", "IRC_SG2" ]
  REAL(DP), DIMENSION(Nspecall), PARAMETER, PUBLIC :: &
    caliberr_specall = [ 0.15_DP, 0.15_DP, 0.2_DP, 0.2_DP, & ! (Decin+2004)
                         0.2_DP, 0.2_DP, 0.2_DP ]            ! (Ohyama+2007)
    
  ! Filter structure
  TYPE, PUBLIC :: filter_type
    INTEGER :: Nfilt, Nwmax
    INTEGER, DIMENSION(:), ALLOCATABLE :: Nw
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: wave, nu, trans
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wcen, nucen
    CHARACTER(lenfilter), DIMENSION(:), ALLOCATABLE :: namefilt
  END TYPE filter_type
  
  ! Calibration error structure
  TYPE, PUBLIC :: calib_type
    INTEGER :: Nfilt    ! number of non fully correlated filters
    INTEGER :: Nfiltall ! number of input filters (Nfiltall >= Nfilt)
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matR, matS, matcov
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matRall, matSall, matcovall
    REAL(DP), DIMENSION(:), ALLOCATABLE :: sigma, sigmall ! standard-deviation
    LOGICAL, DIMENSION(:), ALLOCATABLE :: bool ! F for fully correlated filters
    CHARACTER(lenfilter), DIMENSION(:), ALLOCATABLE :: namefilt, instrument
    CHARACTER(lenfilter), DIMENSION(:), ALLOCATABLE :: namefiltall,instrumentall
    INTEGER, DIMENSION(:), ALLOCATABLE :: labspecall ! 0 if broad band; index of
                                                     ! the SPECALL instrument
                                                     ! else.
    INTEGER, DIMENSION(Nspecall) :: ispecmod ! index of the wavelength
                                             ! for which the calibration
                                             ! error is assumed non fully
                                             ! correlated.
  END TYPE calib_type


CONTAINS


  !==========================================================================
  ! wcen[N] = WCEN_FILTER(filter[N],spec,Nspec,maskspec,instrument)
  !
  !  If the filter belongs to the FILTALL list, then the nominal wavelength
  ! wavelength in microns is returned. If the filter has the form
  ! (spectrograph label)//(wavelength in microns) then the number is returned.
  !==========================================================================

  FUNCTION wcen_filter_scl (filter,spec,instrument)

    USE utilities, ONLY: DP, trimeq
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: filter
    LOGICAL, INTENT(OUT), OPTIONAL :: spec
    CHARACTER(lenfilter), INTENT(OUT), OPTIONAL :: instrument
    REAL(DP) :: wcen_filter_scl

    INTEGER :: i
    CHARACTER(lenfilter) :: instr
    LOGICAL :: spc
    
    !----------------------------------------------------------------------

    instr = "Broad Band"
    DO i=1,Nspecall
      spc = (INDEX(filter(:),TRIM(specall(i))) > 0)
      IF (spc) THEN
        instr = specall(i)
        EXIT
      END IF
    END DO
    IF (spc) THEN
      READ(filter(LEN_TRIM(instr)+1:),*) wcen_filter_scl
    ELSE
      wcen_filter_scl = MINVAL(PACK(wcenall(:),TRIMEQ(filtall(:),filter)))
    END IF
    IF (PRESENT(spec)) spec = spc
    IF (PRESENT(instrument)) instrument = instr

    !----------------------------------------------------------------------

  END FUNCTION wcen_filter_scl

  !------------------------------------------------------------------------

  FUNCTION wcen_filter_1D (filter,spec,Nspec,maskspec,instrument)

    USE utilities, ONLY: DP, trimeq
    IMPLICIT NONE

    CHARACTER(*), DIMENSION(:), INTENT(IN) :: filter
    LOGICAL, INTENT(OUT), OPTIONAL :: spec
    INTEGER, INTENT(OUT), OPTIONAL :: Nspec
    LOGICAL, DIMENSION(SIZE(filter)), INTENT(OUT), OPTIONAL :: maskspec
    CHARACTER(lenfilter), DIMENSION(SIZE(filter)), INTENT(OUT), OPTIONAL :: &
      instrument
    REAL(DP), DIMENSION(SIZE(filter)) :: wcen_filter_1D

    INTEGER :: i, Nfilt
    CHARACTER(lenfilter), DIMENSION(SIZE(filter)) :: instr
    LOGICAL, DIMENSION(SIZE(filter)) :: spc
    
    !----------------------------------------------------------------------

    Nfilt = SIZE(filter)
    DO i=1,Nfilt 
      wcen_filter_1D(i) = WCEN_FILTER_SCL(filter(i),SPEC=spc(i), &
                                          INSTRUMENT=instr(i))
    END DO 
    IF (PRESENT(spec)) spec = ANY(spc(:))
    IF (PRESENT(Nspec)) Nspec = COUNT(spc(:))
    IF (PRESENT(maskspec)) maskspec(:) = spc(:)
    IF (PRESENT(instrument)) instrument(:) = instr(:)
    
    !----------------------------------------------------------------------

  END FUNCTION wcen_filter_1D


  !==========================================================================
  ! CALL READ_FILTERS(filter_type,namefilt[Nfilt])
  !
  !   Read the filter profiles of each waveband.
  !==========================================================================

  SUBROUTINE read_filters (filters,list0)

    USE utilities, ONLY: DP, libF, trimlr, trimeq, strike
    USE inout, ONLY: read_hdf5, h5ext
    USE constants, ONLY: MKS
    USE arrays, ONLY: iwhere, reallocate
    IMPLICIT NONE
  
    TYPE(filter_type), INTENT(OUT) :: filters
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: list0

    INTEGER, PARAMETER :: unit = 1

    INTEGER :: i, j, Nwmax, Nlist
    INTEGER, DIMENSION(:), ALLOCATABLE :: Nw, ind
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wlist, wave_r, trans_r
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: wave, trans, wave2, trans2
    CHARACTER(lenfilter), DIMENSION(:), ALLOCATABLE :: list

    !----------------------------------------------------------------------

    ! 1) Preliminaries
    !-----------------
    ! Input
    IF (PRESENT(list0)) THEN 
      Nlist = SIZE(list0)
      ALLOCATE(list(Nlist))
      list(:) = ADJUSTL(list0(:))
    ELSE 
      Nlist = Nphotall
      ALLOCATE (list(Nlist))
      list(:) = ADJUSTL(filtall(:))
    END IF

    ! Corresponding nominal wavelength
    ALLOCATE (wlist(Nlist))
    sublist: DO i=1,Nlist
      CALL IWHERE(TRIMEQ(filtall(:),list(i)),ind)
      wlist(i) = wcenall(ind(1))
      IF (wlist(i) == 0._DP) CALL STRIKE ("READ_FILTERS","Wrong filter name")
    END DO sublist


    ! 2) Read the filter profiles
    !----------------------------
    ALLOCATE (Nw(Nlist))
    filtread: DO i=1,Nlist

      ! Read the file
      CALL READ_HDF5 (DBLARR1D=wave_r, NAME="Filter wavelength (microns)", &
                      FILE=TRIMLR(libF())//"Instruments/Data/filt_" &
                           //TRIMLR(list(i))//h5ext)
      CALL READ_HDF5 (DBLARR1D=trans_r, NAME="Filter transmission", &
                      FILE=TRIMLR(libF())//"Instruments/Data/filt_" &
                           //TRIMLR(list(i))//h5ext)

      ! Put filters in a single array. Resize the array, if needed.
      Nw(i) = SIZE(wave_r)
      adjust: IF (i == 1) THEN

        Nwmax = MAXVAL(Nw(1:i))
        ALLOCATE (wave(Nlist,Nwmax),trans(Nlist,Nwmax))
        wave(i,1:Nw(i)) = wave_r(:)
        trans(i,1:Nw(i)) = trans_r(:)

      ELSE IF (Nw(i) > Nwmax) THEN 

        Nwmax = MAXVAL(Nw(1:i))
        CALL REALLOCATE (wave2,Nlist,Nwmax)
        CALL REALLOCATE (trans2,Nlist,Nwmax)
        DO j=1,i-1
          wave2(j,1:Nw(j)) = wave(j,1:Nw(j))
          trans2(j,1:Nw(j)) = trans(j,1:Nw(j))
        END DO
        wave2(i,1:Nw(i)) = wave_r(:)
        trans2(i,1:Nw(i)) = trans_r(:)
        CALL REALLOCATE (wave,Nlist,Nwmax)
        CALL REALLOCATE (trans,Nlist,Nwmax)
        wave(:,:) = 0._DP
        trans(:,:) = 0._DP
        wave(:,:) = wave2(:,:)
        trans(:,:) = trans2(:,:)
        DEALLOCATE(wave2,trans2)

      ELSE

        wave(i,1:Nw(i)) = wave_r(1:Nw(i))
        trans(i,1:Nw(i)) = trans_r(1:Nw(i))

      END IF adjust

    END DO filtread


    ! 3) Fill the structure
    !----------------------
    filters%Nwmax = Nwmax
    filters%Nfilt = Nlist
    ALLOCATE (filters%Nw(Nlist), filters%wcen(Nlist), filters%nucen(Nlist), &
              filters%wave(Nlist,Nwmax), filters%trans(Nlist,Nwmax), &
              filters%nu(Nlist,Nwmax), filters%namefilt(Nlist))
    filters%Nw(:) = Nw(:)
    filters%wcen(:) = wlist(:)
    filters%nucen(:) = MKS%clight / MKS%micron / wlist(:)
    filters%wave(:,:) = wave(:,:)
    WHERE (filters%wave(:,:) > 0._DP) 
      filters%nu(:,:) = MKS%clight / MKS%micron / wave(:,:)
    END WHERE
    filters%trans(:,:) = trans(:,:)
    filters%namefilt(:) = list(:)

    ! Free memory space
    DEALLOCATE (Nw,ind,wlist,wave_r,trans_r,wave,trans,list)
 
    !----------------------------------------------------------------------

  END SUBROUTINE read_filters


  !==========================================================================
  ! CALL READ_CALIBERR (calib,list[Nfilt])
  !
  !   Read the covariance matrix of instrumental calibration uncertainties.
  ! The boolean array SUBDIM returns True if correlation is not unity with
  ! band.
  !==========================================================================

  SUBROUTINE read_caliberr (calib,list)

    USE utilities, ONLY: DP, libF, trimlr, trimeq, strike
    USE arrays, ONLY: iwhere
    USE inout, ONLY: read_hdf5, h5ext
    IMPLICIT NONE

    TYPE(calib_type), INTENT(OUT) :: calib
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: list

    INTEGER, PARAMETER :: unit = 1

    INTEGER :: i, j, Nfilt, Nlist, Nbroad, Nspec, Nspecmod, Ncalib
    INTEGER, DIMENSION(:), ALLOCATABLE :: ind, ibroad, ind1, ind2
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wband, wcen, sigbroad
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matS0, matR0, matcov0
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: matS, matR, matcov, matI, matIall
    CHARACTER(lenfilter), DIMENSION(:), ALLOCATABLE :: list1, namefilt, instr
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: boolcal, boolcalall

    !----------------------------------------------------------------------

    ! 1) Input
    !---------
    IF (PRESENT(list)) THEN 
      Nlist = SIZE(list)
      ALLOCATE(list1(Nlist))
      list1(:) = ADJUSTL(list(:))
    ELSE 
      Nlist = Nphotall
      ALLOCATE (list1(Nlist))
      list1(:) = ADJUSTL(filtall(:))
    END IF

    ! Spectrum
    ALLOCATE (wcen(Nlist),instr(Nlist))
    wcen(:) = WCEN_FILTER(list1(:),NSPEC=Nspec,INSTRUMENT=instr(:))
    Nbroad = Nlist - Nspec

    ! Fill-in the structure
    ALLOCATE (calib%namefiltall(Nlist),calib%instrumentall(Nlist), &
              calib%labspecall(Nlist))
    calib%namefiltall(:) = list(:)
    calib%instrumentall(:) = instr(:)
    calib%labspecall(:) = 0
    DO i=1,Nspecall
      WHERE (instr(:) == specall(i)) calib%labspecall(:) = i
    END DO
    

    ! 2) Spectrum uncertainties
    !--------------------------
    ! We currently consider that the calibration uncertainty of a given spectral
    ! module is perfectly correlated between wavelengths, and independent of
    ! the other modules and instruments. Therefore, we use only one independent
    ! calibration error variable per module.
    Nspecmod = 0
    calib%ispecmod(:) = 0 
    DO i=1,Nspecall
      IF (ANY(TRIMEQ(instr(:),specall(i)))) THEN
        CALL IWHERE(TRIMEQ(instr(:),specall(i)),ind)
        Nspecmod = Nspecmod + 1
        calib%ispecmod(i) = ind(1)
      END IF
    END DO

    ! Subdimension
    ALLOCATE (calib%bool(Nlist))
    calib%bool(:) = .False.
    WHERE (TRIMEQ(instr(:),"Broad Band")) calib%bool(:) = .True.
    IF (ANY(.NOT. calib%bool(:))) &
      calib%bool(PACK(calib%ispecmod(:),calib%ispecmod(:) > 0)) = .True.

    
    ! 3) Broadband uncertainties
    !---------------------------
    CALL READ_HDF5(STRARR1D=namefilt,NAME="Filter label", &
                    FILE=TRIMLR(libF())//"Instruments/Data/caliberr"//h5ext, &
                    N1=Nfilt)
    CALL READ_HDF5(DBLARR1D=wband,NAME="Filter wavelength (microns)", &
                    FILE=TRIMLR(libF())//"Instruments/Data/caliberr"//h5ext)
    CALL READ_HDF5(DBLARR2D=matS0,NAME="Standard-deviation matrix", &
                    FILE=TRIMLR(libF())//"Instruments/Data/caliberr"//h5ext)
    CALL READ_HDF5(DBLARR2D=matR0,NAME="Correlation matrix", &
                    FILE=TRIMLR(libF())//"Instruments/Data/caliberr"//h5ext)
    CALL READ_HDF5(DBLARR2D=matcov0,NAME="Covariance matrix", &
                    FILE=TRIMLR(libF())//"Instruments/Data/caliberr"//h5ext)
    
    ! Select the sublist of selected filters
    ALLOCATE (matS(Nbroad,Nbroad),matR(Nbroad,Nbroad),matcov(Nbroad,Nbroad))
    CALL IWHERE(TRIMEQ(instr(:),"Broad Band"),ibroad)
    DO i=1,Nbroad
      CALL IWHERE(TRIMEQ(namefilt(:),list1(ibroad(i))),ind1)
      DO j=1,Nbroad
        CALL IWHERE(TRIMEQ(namefilt(:),list1(ibroad(j))),ind2)
        matS(i,j) = matS0(ind1(1),ind2(1))
        matR(i,j) = matS0(ind1(1),ind2(1))         
      END DO
    END DO

    
    ! 4) Fill the structure
    !----------------------
    Ncalib = Nbroad + Nspecmod

    ! Initialize the structure
    calib%Nfilt = Ncalib
    ALLOCATE (calib%namefilt(Ncalib),calib%instrument(Ncalib), &
              calib%matS(Ncalib,Ncalib),calib%matR(Ncalib,Ncalib), &
              calib%matcov(Ncalib,Ncalib),calib%sigma(Ncalib), &
              calib%matSall(Nlist,Nlist),calib%matRall(Nlist,Nlist), &
              calib%matcovall(Nlist,Nlist),calib%sigmall(Nlist))

    ! List
    calib%namefilt(:) = PACK(list1(:),calib%bool(:))
    calib%instrument(:) = PACK(instr(:),calib%bool(:))

    ! Mask
    ALLOCATE (boolcal(Ncalib,Ncalib),boolcalall(Nlist,Nlist))
    FORALL (i=1:Ncalib,j=1:Ncalib) &
      boolcal(i,j) = ( TRIMEQ(calib%instrument(i),"Broad Band") &
                       .AND. TRIMEQ(calib%instrument(j),"Broad Band") )
    FORALL (i=1:Nlist,j=1:Nlist) &
      boolcalall(i,j) = ( TRIMEQ(calib%instrumentall(i),"Broad Band") &
                          .AND. TRIMEQ(calib%instrumentall(j),"Broad Band") )
    
    ! Identity matrix
    ALLOCATE (matI(Ncalib,Ncalib),matIall(Nlist,Nlist))
    matI(:,:) = 0._DP
    FORALL (i=1:Ncalib) matI(i,i) = 1._DP
    matIall(:,:) = 0._DP
    FORALL (i=1:Nlist) matIall(i,i) = 1._DP
    
    ! Broad bands
    calib%matS(:,:) = UNPACK(RESHAPE(matS(:,:),[Nbroad*Nbroad]), &
                             MASK=boolcal(:,:),FIELD=0._DP)
    calib%matR(:,:) = UNPACK(RESHAPE(matR(:,:),[Nbroad*Nbroad]), &
                             MASK=boolcal(:,:),FIELD=matI(:,:))
    ALLOCATE (sigbroad(Nbroad))
    sigbroad(:) = [(matS(i,i),i=1,Nbroad)]
    calib%sigmall(:) = UNPACK(sigbroad(:),TRIMEQ(instr(:),"Broad Band"), &
                              FIELD=0._DP)
    calib%matSall(:,:) = UNPACK(RESHAPE(matS(:,:),[Nbroad*Nbroad]), &
                                MASK=boolcalall(:,:),FIELD=0._DP)
    calib%matRall(:,:) = UNPACK(RESHAPE(matR(:,:),[Nbroad*Nbroad]), &
                                MASK=boolcalall(:,:),FIELD=matIall(:,:))
    
    ! Spectral modules
    IF (Nspecmod > 0) THEN
      DO i=1,Nspecall
        IF (ANY(TRIMEQ(calib%instrument(:),specall(i)))) THEN
          CALL IWHERE(TRIMEQ(calib%instrument(:),specall(i)),ind)
          calib%matS(ind(1),ind(1)) = caliberr_specall(i)
          CALL IWHERE(TRIMEQ(instr(:),specall(i)),ind)
          calib%sigmall(ind(:)) = caliberr_specall(i)
          FORALL (j=1:SIZE(ind)) &
            calib%matSall(ind(j),ind(j)) = caliberr_specall(i)
        END IF
      END DO
    END IF

    ! Final covariance matrix
    FORALL (i=1:Ncalib) calib%sigma(i) = calib%matS(i,i)
    calib%matcov(:,:) = MATMUL(MATMUL(calib%matS(:,:),calib%matR(:,:)), &
                               calib%matS(:,:))
    calib%matcovall(:,:) = MATMUL(MATMUL(calib%matSall(:,:), &
                                         calib%matRall(:,:)), &
                                  calib%matSall(:,:))
    
    ! Free memory space
    DEALLOCATE (wband,wcen,matS0,matR0,matcov0,matS,matR,matcov,list1,namefilt,&
                boolcal,boolcalall,matI,matIall,sigbroad)
    IF (ALLOCATED(ind)) DEALLOCATE(ind)
    
    !----------------------------------------------------------------------

  END SUBROUTINE read_caliberr


  !==========================================================================
  ! Fnu0[N] = SYNTHETIC_PHOTOMETRY (wave[Nw],Fnu[Nw],namefilt[N],filters, &
  !                                 [wnom[N]],/POINT,/EXTENDED)
  !
  !   Perform synthetic photometry from a given SED in several instrumental
  ! filters. WAVE must be sorted. By default, we assume extended SPIRE 
  ! calibration. 
  !==========================================================================

  FUNCTION synthetic_photometry (wave,Fnu,namefilt,filters,wnom,point, &
                                  extended,accuracy)

    USE utilities, ONLY: DP, trimlr, trimeq, strike
    USE arrays, ONLY: reallocate
    USE constants, ONLY: MKS
    USE integration, ONLY: integ_tabulated
    USE interpolation, ONLY: interp_lin_sorted
    USE adaptative_grid, ONLY: gridadapt1D
    USE instrument_filters_ext, ONLY: a1, a2, nu0, wfilt, transfilt, &
                                      wSED, FnuSED, colcorr

    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN) :: wave, Fnu
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: namefilt
    TYPE(filter_type), INTENT(IN) :: filters
    REAL(DP), DIMENSION(SIZE(namefilt)), INTENT(OUT), OPTIONAL :: wnom
    REAL(DP), DIMENSION(SIZE(namefilt)) :: synthetic_photometry
    LOGICAL, INTENT(IN), OPTIONAL :: point, extended
    REAL(DP), INTENT(IN), OPTIONAL :: accuracy

    INTEGER :: i, j, Nfilt, ifilt, Nw, Nad
    REAL(DP) :: acc, norm, sclSED
    REAL(DP), DIMENSION(:), ALLOCATABLE :: wad, yfine, prim, nuad, nuwad, nu0ad
    REAL(DP), DIMENSION(:), ALLOCATABLE :: transad
    LOGICAL :: extend

    !----------------------------------------------------------------------

    ! Accuracy
    IF (PRESENT(accuracy)) THEN 
      acc = accuracy
    ELSE
      acc = 1.E-3_DP
    END IF

    ! Array sizes
    Nw = SIZE(wave(:))
    Nfilt = SIZE(namefilt(:))

    ! SPIRE point source vs extended source calibration
    extend = .True.
    IF (PRESENT(point)) extend = (.NOT. point)
    IF (PRESENT(extended)) extend = extended

    ! Integrate over the filter profile
    filtinteg: DO i=1,Nfilt
      ifilt = MAXVAL(PACK([(j,j=1,filters%Nfilt)], &
                          TRIMEQ(filters%namefilt(:),namefilt(i))))

      ! Implement each flux convention
      bandpass: SELECT CASE (TRIMLR(namefilt(i)))
        CASE ("IRAS1","IRAS2","IRAS3","IRAS4", &
              "IRAC1","IRAC2","IRAC3","IRAC4") 
          a1 = 1._DP
          a2 = 2._DP
        CASE ("2MASS1","2MASS2","2MASS3", &
              "WISE1","WISE2","WISE3","WISE4","WISE4_old")
          a1 = 1._DP
          a2 = 3._DP
        CASE ("MIPS1","MIPS2","MIPS3")
          a1 = 0._DP
          a2 = -2._DP
        CASE ("AKARI_IRC1","AKARI_IRC2","AKARI_IRC3", &
              "AKARI_IRC4","AKARI_IRC5","AKARI_IRC6", &
              "AKARI_IRC7","AKARI_IRC8","AKARI_IRC9", &
              "AKARI_FIS1","AKARI_FIS2","AKARI_FIS3", "AKARI_FIS4", &
              "PACS1","PACS2","PACS3", &
              "HFI1","HFI2","HFI3","HFI4","HFI5","HFI6", &
              "LFI1", "LFI2", "LFI3", &
              "DIRBE1","DIRBE2","DIRBE3","DIRBE4","DIRBE5", &
              "DIRBE6","DIRBE7","DIRBE8","DIRBE9","DIRBE10")
          a1 = 0._DP
          a2 = 1._DP
        CASE ("SPIRE1","SPIRE2","SPIRE3") 
          a1 = 2._DP
          a2 = 3._DP
          IF (.NOT. extend) THEN
            a1 = a1 - 2._DP
            a2 = a2 - 2._DP
         END IF
       CASE ("MSX1","MSX2","MSX3","MSX4")
         a1 = 0._DP
         a2 = 0._DP 
       CASE DEFAULT
          CALL STRIKE("SYNTHETIC_PHOTOMETRY","Wrong filter")
      END SELECT bandpass

      ! Filter grid
      nu0 = filters%nucen(ifilt)
      CALL REALLOCATE(wfilt,filters%Nw(ifilt))
      CALL REALLOCATE(transfilt,filters%Nw(ifilt))
      wfilt(:) = filters%wave(ifilt,1:filters%Nw(ifilt))
      transfilt(:) = filters%trans(ifilt,1:filters%Nw(ifilt))

      ! SED grid
      CALL REALLOCATE(wSED,Nw)
      CALL REALLOCATE(FnuSED,Nw)
      wSED(:) = wave(:)
      FnuSED(:) = Fnu(:)
      sclSED = MAXVAL(FnuSED(:))
      FnuSED(:) = FnuSED(:) / sclSED ! rescale to avoid numerical problems

      ! Integration with adaptative grid
      CALL GRIDADAPT1D(wfilt(:),wad,yfine,COLCORR,acc,INTEG=.True., &
                       XLOG=.False.,YLOG=.False.,PRIMITIVE=prim,VERBOSE=.False.)
      Nad = SIZE(wad(:))

      ! Normalization
      CALL REALLOCATE(nuad,Nad)
      CALL REALLOCATE(nuwad,Nad)
      CALL REALLOCATE(nu0ad,Nad)
      CALL REALLOCATE(transad,Nad)
      nuad(:) = MKS%clight / MKS%micron / wad(:)
      nuwad(:) = nuad(:) / wad(:)
      nu0ad(:) = nu0 / nuad(:)
      transad(:) = INTERP_LIN_SORTED(transfilt(:),wfilt(:),wad(:), &
                                     XLOG=.False.,YLOG=.False.)
      norm = INTEG_TABULATED( wad(:), nu0ad(:)**a2 * transad(:) * nuwad(:), &
                              XLOG=.False., YLOG=.False. )

      ! Output
      synthetic_photometry(i) = prim(Nad) * sclSED / norm
      IF (PRESENT(wnom)) wnom(i) = filters%wcen(ifilt)

    END DO filtinteg

    ! Free memory space
    DEALLOCATE (wad,yfine,prim,nuad,nuwad,nu0ad,transad,wfilt,transfilt, &
                wSED,FnuSED)

    !----------------------------------------------------------------------

  END FUNCTION synthetic_photometry


END MODULE instrument_filters

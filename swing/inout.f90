!******************************************************************************
!*
!*              CUSTOMISED INPUT/OUPUT AND INTERACTION WITH IDL
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 08/2007.
  !    - 02/2015: clarify the header of WRITE_DATA/READ_DATA and parametrize
  !      the format of the double precision writing.
  !    - 02/2015: add the HDF5 routines, and rename the ASCII and BINARY 
  !      subroutines.
  !    - 07/2016: add the fractionation routines.
  ! 
  ! 3) DESCRIPTION: Package to read and write customised data files.
  !==========================================================================


MODULE inout

  USE utilities, ONLY:
  USE hdf5, ONLY: hsize_t
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: write_ascii, read_ascii, write_binary, read_binary
  PUBLIC :: write_hdf5, read_hdf5, getdim_hdf5, write_hdf5_frac, read_hdf5_frac
  PUBLIC :: ind_hdf5_frac, check_hdf5, read_input_line, edit_input_line

  INTERFACE write_binary
    MODULE PROCEDURE write_binary_1D, write_binary_2D, write_binary_3D, &
                     write_binary_4D
  END INTERFACE write_binary

  INTERFACE read_binary
    MODULE PROCEDURE read_binary_1D, read_binary_2D, read_binary_3D, &
                     read_binary_4D
  END INTERFACE read_binary

  ! Default extensions
  CHARACTER(*), PARAMETER, PUBLIC :: ascext = ".txt" ! default ASCII extension
  CHARACTER(*), PARAMETER, PUBLIC :: binext = ".dat" ! default binary extension
  CHARACTER(*), PARAMETER, PUBLIC :: h5ext = ".h5"   ! default HDF5 extension

  ! Default parameters for the HDF5 format
  INTEGER, PARAMETER :: rank0 = 6 ! rank of dataset
  INTEGER, PARAMETER :: arank = 1
  INTEGER(hsize_t), DIMENSION(arank), PARAMETER :: adims = [3]
  INTEGER(hsize_t), PARAMETER :: chunkdims = 20
  LOGICAL, PARAMETER :: compress_def = .True.
  CHARACTER(*), PARAMETER :: dsetname_def = "Array stored with WRITE_HDF5"
  CHARACTER(*), PARAMETER :: aname = "WRITE_HDF5 information"

  ! Format for the generic ASCII input/output routines
  INTEGER, PARAMETER, PUBLIC :: textwid = 80    ! Lenght of line for these files
  INTEGER, PARAMETER, PUBLIC :: unitdata = 1    ! Default unit
  INTEGER, PARAMETER, PUBLIC :: Ndecim_def =  8 ! Default number of decimals
  INTEGER, PARAMETER :: lencom = 48             ! Length of parameter comment
  CHARACTER(*), PARAMETER :: formline = "(A80)" ! String on the whole line
  CHARACTER(*), PARAMETER :: formparn = "(A12)" ! Parameter name
  CHARACTER(*), PARAMETER :: formint = "(I16)"  ! Integer parameter
  CHARACTER(*), PARAMETER :: formparc = "(A48)" ! Comment for parameter

  ! Length of parameter label and input line for READ_INPUT_LINE
  INTEGER, PARAMETER, PUBLIC :: lenpar = 30       ! length of parameter value
  INTEGER, PARAMETER, PUBLIC :: lenline = 150     ! length of line
  INTEGER, PARAMETER, PUBLIC :: lenpath = 200     ! length of path


CONTAINS


  !=========================================================================
  ! CALL WRITE_HDF5(DBLARR{1-6}D or STRARR1D or INTARR{1-3}D,FILE="", &
  !                 COMPRESS=T/F,INITDBLARR=[],INITINTARR=[], &
  !                 APPEND=T/F,NAME="",IND1=[idim1_inf,idim1_sup], &
  !                 IND2=[idim2_inf,idm2_sup],IND3=[idim3_inf,idm3_sup], &
  !                 IND4=[idim4_inf,idm4_sup],IND5=[idim5_inf,idm5_sup], &
  !                 IND6=[idim5_inf,idim6_sup],UNIT=)
  ! 
  !   Write a 1D to 6D double precision array or a 1D string array to an 
  ! HDF5 file, with or without compression.
  !   NAME should not contain any "/", "[" or "]".
  !   If INIT* is set, then the array is just initialized, supposed to be filled
  ! later with the IND* keyword. The value of INIT* is the list of dimensions.
  !=========================================================================

  SUBROUTINE write_hdf5 (dblarr1D,dblarr2D,dblarr3D,dblarr4D,dblarr5D,dblarr6D, &
                         strarr1D,intarr1D,intarr2D,intarr3D,initdblarr, &
                         initintarr, &
                         file,name,compress,append,ind1,ind2,ind3,ind4,ind5, &
                         ind6,unit,verbose)

    USE utilities, ONLY: DP, trimlr, verbatim, programrunner, today, strike, &
                         ustd, strreplace
    USE hdf5, ONLY: hid_t, hsize_t, size_t, h5f_acc_trunc_f, h5t_native_double,&
                    h5t_native_integer, h5t_native_character, h5s_unlimited_f, &
                    h5f_acc_rdwr_f, h5p_dataset_create_f, &
                    h5open_f, h5close_f, h5fopen_f, h5fcreate_f, h5fclose_f, &
                    h5sselect_hyperslab_f, h5s_select_set_f, h5dset_extent_f, &
                    h5screate_simple_f, h5sclose_f, h5dcreate_f, h5dclose_f, &
                    h5dwrite_f, h5pcreate_f, h5pset_chunk_f, h5pset_deflate_f, &
                    h5tcopy_f, h5tset_size_f, h5dopen_f, h5dget_space_f, &
                    h5acreate_f, h5awrite_f, h5aclose_f
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: dblarr1D
    REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: dblarr2D
    REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: dblarr3D
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL :: dblarr4D
    REAL(DP), DIMENSION(:,:,:,:,:), INTENT(IN), OPTIONAL :: dblarr5D
    REAL(DP), DIMENSION(:,:,:,:,:,:), INTENT(IN), OPTIONAL :: dblarr6D
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: strarr1D
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: intarr1D
    INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: intarr2D
    INTEGER, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: intarr3D
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: initdblarr, initintarr
    CHARACTER(*), INTENT(IN) :: file
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    LOGICAL, INTENT(IN), OPTIONAL :: compress, append
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: ind1, ind2, ind3, ind4, ind5
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: ind6
    INTEGER, INTENT(IN), OPTIONAL :: unit
    LOGICAL, INTENT(IN), OPTIONAL :: verbose

    INTEGER(hsize_t), DIMENSION(rank0) :: dims    ! dimension of dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: dims0   ! dimension of initial dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: maxdims ! max. dimension of dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: cdims   ! dimension of chunked dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: hyperdims ! dimension of dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: offset             ! hyperslab offset
    INTEGER(hsize_t), DIMENSION(2) :: i1, i2, i3, i4, i5, i6 ! minmax indices
    INTEGER, DIMENSION(2) :: u ! unit for notification
    INTEGER(hid_t) :: file_ID    ! file identifier
    INTEGER(hid_t) :: plist_ID   ! property list identifier
    INTEGER(hid_t) :: dspace_ID  ! dataspace identifier
    INTEGER(hid_t) :: dset_ID    ! dataset identifier
    INTEGER(hid_t) :: mspace_ID  ! memory space identifier
    INTEGER(hid_t) :: dtype_ID   ! dataset type identifier
    INTEGER(hid_t) :: atype_ID   ! attribute type identifier
    INTEGER(hid_t) :: aspace_ID  ! attribute's dataspace identifier
    INTEGER(hid_t) :: attr_ID    ! attribute identifier
    INTEGER :: i, strlen
    INTEGER :: rank ! actual rank
    INTEGER :: error ! error flag
    LOGICAL :: cmprss, newfile, subset, appnd, string, integ, verb, init
    CHARACTER(textwid), DIMENSION(adims(1)) :: adata
    CHARACTER(textwid) :: dsetname

    !-------------------------------------------------------------------------

    ! 1) Keywords and Preliminary settings
    !-------------------------------------
    ! Type of writing
    subset = ( PRESENT(ind1) .OR. PRESENT(ind2) .OR. PRESENT(ind3) &
               .OR. PRESENT(ind4) .OR. PRESENT(ind5) .OR. PRESENT(ind6) )
    IF (PRESENT(append)) THEN 
      appnd = append
    ELSE
      appnd = .False.
    END IF
    newfile = ( .NOT. appnd .AND. .NOT. subset )
    init = .False.
    IF (PRESENT(initdblarr) .OR. PRESENT(initintarr)) THEN
      init = .True.
      subset = .False.
      dims0(:) = 1
      maxdims(:) = h5s_unlimited_f
    END IF
    
    ! Dataset name (used to identify the struture in the HDF5 file)
    IF (PRESENT(name)) THEN
      dsetname = TRIMLR(STRREPLACE(name,["/"],["ov"]))
    ELSE
      dsetname = dsetname_def
    END IF

    ! Actual rank and dimensions
    rank = -1
    string = .False.
    integ = .False.
    IF (PRESENT(dblarr1D)) THEN
      rank = 1
      dims(1:rank) = [ SIZE(dblarr1D(:),1) ]
    END IF
    IF (PRESENT(dblarr2D)) THEN
      IF (rank == -1) THEN
        rank = 2
        dims(1:rank) = [ SIZE(dblarr2D(:,:),1), SIZE(dblarr2D(:,:),2) ]
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr3D)) THEN
      IF (rank == -1) THEN
        rank = 3
        dims(1:rank) = [ SIZE(dblarr3D(:,:,:),1), SIZE(dblarr3D(:,:,:),2), &
                         SIZE(dblarr3D(:,:,:),3) ]
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr4D)) THEN
      IF (rank == -1) THEN
        rank = 4
        dims(1:rank) = [ SIZE(dblarr4D(:,:,:,:),1), SIZE(dblarr4D(:,:,:,:),2), &
                         SIZE(dblarr4D(:,:,:,:),3), SIZE(dblarr4D(:,:,:,:),4) ]
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr5D)) THEN
      IF (rank == -1) THEN
        rank = 5
        dims(1:rank) = [ SIZE(dblarr5D(:,:,:,:,:),1), &
                         SIZE(dblarr5D(:,:,:,:,:),2), &
                         SIZE(dblarr5D(:,:,:,:,:),3), &
                         SIZE(dblarr5D(:,:,:,:,:),4), &
                         SIZE(dblarr5D(:,:,:,:,:),5) ]
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr6D)) THEN
      IF (rank == -1) THEN
        rank = 6
        dims(1:rank) = [ SIZE(dblarr6D(:,:,:,:,:,:),1), &
                         SIZE(dblarr6D(:,:,:,:,:,:),2), &
                         SIZE(dblarr6D(:,:,:,:,:,:),3), &
                         SIZE(dblarr6D(:,:,:,:,:,:),4), &
                         SIZE(dblarr6D(:,:,:,:,:,:),5), &
                         SIZE(dblarr6D(:,:,:,:,:,:),6) ]
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    strorDP: IF (PRESENT(strarr1D)) THEN
      IF (rank == -1) THEN
        rank = 1
        dims(1:rank) = [ SIZE(strarr1D(:),1) ]
        string = .True.
        strlen = LEN(strarr1D(:))
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF strorDP
    IF (PRESENT(intarr1D)) THEN
      IF (rank == -1) THEN
        rank = 1
        dims(1:rank) = [ SIZE(intarr1D(:),1) ]
        integ = .True.
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(intarr2D)) THEN
      IF (rank == -1) THEN
        rank = 2
        dims(1:rank) = [ SIZE(intarr2D(:,:),1), SIZE(intarr2D(:,:),2) ]
        integ = .True.
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(intarr3D)) THEN
      IF (rank == -1) THEN
        rank = 3
        dims(1:rank) = [ SIZE(intarr3D(:,:,:),1), SIZE(intarr3D(:,:,:),2), &
                         SIZE(intarr3D(:,:,:),3) ]
        integ = .True.
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(initdblarr)) THEN
      IF (rank == -1) THEN
        rank = SIZE(initdblarr(:))
        dims(1:rank) = initdblarr(:)
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(initintarr)) THEN
      IF (rank == -1) THEN
        rank = SIZE(initintarr(:))
        dims(1:rank) = initintarr(:)
        integ = .True.
      ELSE
        CALL STRIKE("WRITE_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (rank == -1) &
      CALL STRIKE("WRITE_HDF5","you have to select at least one ARRAY")

    ! Compression settings
    cmprss = compress_def
    IF (PRESENT(compress)) cmprss = compress
    cdims(:) = chunkdims
    IF (ALL(dims(1:rank) < cdims(1:rank))) THEN
      cmprss = .False.
    ELSE
      WHERE (dims(1:rank) < cdims(1:rank)) cdims(1:rank) = dims(1:rank)
    END IF

    ! Consistency check
    IF (PRESENT(ind2) .AND. rank < 2) &
      CALL STRIKE("WRITE_HDF5", &
                  "you cannot select IND2 if ARRAY is less than 2D")
    IF (PRESENT(ind3) .AND. rank < 3) &
      CALL STRIKE("WRITE_HDF5", &
                  "you cannot select IND3 if ARRAY is less than 3D")
    IF (PRESENT(ind4) .AND. rank < 4) &
      CALL STRIKE("WRITE_HDF5", &
                  "you cannot select IND4 if ARRAY is less than 4D")
    IF (PRESENT(ind5) .AND. rank < 5) &
      CALL STRIKE("WRITE_HDF5", &
                  "you cannot select IND5 if ARRAY is less than 5D")
    IF (appnd .AND. subset) &
      CALL STRIKE("WRITE_HDF5","you can not write a subset in a new ARRAY")


    ! 2) File preparation
    !--------------------
    ! Initialize the HDF5 interface
    CALL H5OPEN_F(ERROR=error)

    ! Open the file
    createoropen: IF (newfile) THEN
      CALL H5FCREATE_F(NAME=file,ACCESS_FLAGS=h5f_acc_trunc_f,FILE_ID=file_ID, &
                       HDFERR=error) ! create new file
    ELSE
      CALL H5FOPEN_F(NAME=file,ACCESS_FLAGS=h5f_acc_rdwr_f,FILE_ID=file_ID, &
                     HDFERR=error) ! append existing file
    END IF createoropen


    ! 3) Write the dataset with or without compression
    !-------------------------------------------------
    ! Define the type of the dataset (DP or string)
    strfullcompr: IF (.NOT. string) THEN
      dtype_ID = MERGE(h5t_native_integer,h5t_native_double,integ)
    ELSE
      CALL H5TCOPY_F(TYPE_ID=h5t_native_character,NEW_TYPE_ID=dtype_ID, &
                     HDFERR=error) ! define the type (character)
      CALL H5TSET_SIZE_F(TYPE_ID=dtype_ID,SIZE=INT(strlen,hsize_t), &
                         HDFERR=error) ! string
    END IF strfullcompr

    ! Actual writing
    fullorpartial: IF (.NOT. subset) THEN

      ! New dataset
      IF (.NOT. init) THEN
        CALL H5SCREATE_SIMPLE_F(RANK=rank,DIMS=dims(1:rank),SPACE_ID=dspace_ID,&
                                HDFERR=error) ! create dataspace
      ELSE
        CALL H5SCREATE_SIMPLE_F(RANK=rank,DIMS=dims0(1:rank), &
                                SPACE_ID=dspace_ID,HDFERR=error, &
                                MAXDIMS=maxdims(1:rank))
      END IF
      fullcompression: IF (cmprss .OR. init) THEN
        CALL H5PCREATE_F(CLASS=h5p_dataset_create_f,PRP_ID=plist_ID, &
                         HDFERR=error) ! create property list
        CALL H5PSET_CHUNK_F(PRP_ID=plist_ID,NDIMS=rank,DIMS=cdims(1:rank), &
                            HDFERR=error) ! chunk dataset
        IF (.NOT. init) &
          CALL H5PSET_DEFLATE_F(PRP_ID=plist_ID,LEVEL=6,HDFERR=error) ! compress
        CALL H5DCREATE_F(LOC_ID=file_ID,NAME=dsetname,TYPE_ID=dtype_ID, &
                         SPACE_ID=dspace_ID,DSET_ID=dset_ID,HDFERR=error, &
                         DCPL_ID=plist_ID) ! create dataset with default props
      ELSE
        CALL H5DCREATE_F(LOC_ID=file_ID,NAME=dsetname,TYPE_ID=dtype_ID, & 
                         SPACE_ID=dspace_ID,DSET_ID=dset_ID, &
                         HDFERR=error) ! create dataset with default props
      END IF fullcompression
     
      ! Write the data (minimum number of dimensions)
      IF (.NOT. init) THEN
        writefull: SELECT CASE (rank)
          CASE (1)
            IF (.NOT. string) THEN
              IF (.NOT. integ) THEN
                CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                                BUF=dblarr1D(:),DIMS=dims(1:rank),HDFERR=error)
              ELSE
                CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                                BUF=intarr1D(:),DIMS=dims(1:rank),HDFERR=error)
              END IF
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=strarr1D(:),DIMS=dims(1:rank),HDFERR=error)
            END IF
          CASE (2)
            IF (.NOT. integ) THEN
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=dblarr2D(:,:),DIMS=dims(1:rank),HDFERR=error)
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=intarr2D(:,:),DIMS=dims(1:rank),HDFERR=error)
            END IF
          CASE (3)
            IF (.NOT. integ) THEN
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                             BUF=dblarr3D(:,:,:),DIMS=dims(1:rank),HDFERR=error)
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                             BUF=intarr3D(:,:,:),DIMS=dims(1:rank),HDFERR=error)
            END IF
          CASE (4)
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=dblarr4D(:,:,:,:),DIMS=dims(1:rank),HDFERR=error)
          CASE (5)
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=dblarr5D(:,:,:,:,:),DIMS=dims(1:rank), &
                            HDFERR=error)
          CASE (6)
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=dblarr6D(:,:,:,:,:,:),DIMS=dims(1:rank), &
                            HDFERR=error)
        END SELECT writefull
      ELSE
        writeinit: SELECT CASE (rank)
          CASE (1)
            IF (.NOT. integ) THEN
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=[0._DP],DIMS=dims0(1:rank),HDFERR=error)
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=[0],DIMS=dims0(1:rank),HDFERR=error)
            END IF
          CASE (2)
            IF (.NOT. integ) THEN
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=[0._DP,0._DP],DIMS=dims0(1:rank),HDFERR=error)
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=[0,0],DIMS=dims0(1:rank),HDFERR=error)
            END IF
          CASE (3)
            IF (.NOT. integ) THEN
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=[0._DP,0._DP,0._DP],DIMS=dims0(1:rank), &
                              HDFERR=error)
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=[0,0,0],DIMS=dims0(1:rank),HDFERR=error)
            END IF
          CASE (4)
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=[0._DP,0._DP,0._DP,0._DP],DIMS=dims0(1:rank), &
                            HDFERR=error)
          CASE (5)
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=[0._DP,0._DP,0._DP,0._DP,0._DP], &
                            DIMS=dims0(1:rank),HDFERR=error)
          CASE (6)
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=[0._DP,0._DP,0._DP,0._DP,0._DP,0._DP], &
                            DIMS=dims0(1:rank),HDFERR=error)
        END SELECT writeinit        
        CALL H5DSET_EXTENT_F(DATASET_ID=dset_ID,SIZE=dims(1:rank),HDFERR=error)
      END IF

    ELSE

      ! Already existing dataset
      IF (PRESENT(ind1)) THEN 
        i1(:) = ind1(:) 
      ELSE 
        i1(:) = [ INT(1,hsize_t), dims(1) ]
      END IF
      IF (PRESENT(ind2)) THEN 
        i2(:) = ind2(:) 
      ELSE
        i2(:) = [ INT(1,hsize_t), dims(2) ]
      END IF
      IF (PRESENT(ind3)) THEN 
        i3(:) = ind3(:) 
      ELSE 
        i3(:) = [ INT(1,hsize_t), dims(3) ]
      END IF
      IF (PRESENT(ind4)) THEN 
        i4(:) = ind4(:) 
      ELSE 
        i4(:) = [ INT(1,hsize_t), dims(4) ]
      END IF
      IF (PRESENT(ind5)) THEN 
        i5(:) = ind5(:) 
      ELSE 
        i5(:) = [ INT(1,hsize_t), dims(5) ]
      END IF
      IF (PRESENT(ind6)) THEN 
        i6(:) = ind6(:) 
      ELSE 
        i6(:) = [ INT(1,hsize_t), dims(6) ]
      END IF
      offset(:) = [ i1(1), i2(1), i3(1), i4(1), i5(1), i6(1) ] - 1 ! offset
                                                                   ! starts at 0
      hyperdims(:) = [ i1(2)-i1(1)+1, i2(2)-i2(1)+1, i3(2)-i3(1)+1, &
                       i4(2)-i4(1)+1, i5(2)-i5(1)+1, i6(2)-i6(1)+1 ]
      CALL H5DOPEN_F(LOC_ID=file_ID,NAME=dsetname,DSET_ID=dset_ID, &
                     HDFERR=error) ! open dataset
      CALL H5DGET_SPACE_F(DATASET_ID=dset_ID,DATASPACE_ID=dspace_ID, &
                          HDFERR=error) ! load the data space
      subcompression: IF (cmprss) THEN
        CALL H5PCREATE_F(CLASS=h5p_dataset_create_f,PRP_ID=plist_ID, &
                         HDFERR=error) ! create property list
        CALL H5PSET_CHUNK_F(PRP_ID=plist_ID,NDIMS=rank,DIMS=cdims(1:rank), &
                            HDFERR=error) ! chunk dataset
        CALL H5PSET_DEFLATE_F(PRP_ID=plist_ID,LEVEL=6,HDFERR=error) ! compress
        CALL H5SSELECT_HYPERSLAB_F(SPACE_ID=dspace_ID, &
                                   OPERATOR=h5s_select_set_f, &
                                   START=offset(1:rank), &
                                   COUNT=hyperdims(1:rank), &
                                   HDFERR=error) ! select the subset
      ELSE
        CALL H5SSELECT_HYPERSLAB_F(SPACE_ID=dspace_ID, &
                                   OPERATOR=h5s_select_set_f, &
                                   START=offset(1:rank), &
                                   COUNT=hyperdims(1:rank), &
                                   HDFERR=error) ! select the subset
      END IF subcompression
      CALL H5SCREATE_SIMPLE_F(RANK=rank,DIMS=hyperdims(1:rank), &
                              SPACE_ID=mspace_ID, &
                              HDFERR=error) ! create sub-dataspace
 
      ! Write the data (minimum number of dimensions)
      writesub: SELECT CASE (rank)
        CASE (1)
          IF (.NOT. string) THEN
            IF (.NOT. integ) THEN
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=dblarr1D(:),DIMS=hyperdims(1:rank), &
                              HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                              FILE_SPACE_ID=dspace_ID)
            ELSE
              CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                              BUF=intarr1D(:),DIMS=hyperdims(1:rank), &
                              HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                              FILE_SPACE_ID=dspace_ID)
            END IF
          ELSE
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=strarr1D(:),DIMS=hyperdims(1:rank), &
                            HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                            FILE_SPACE_ID=dspace_ID)
          END IF
        CASE (2)
          IF (.NOT. integ) THEN
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=dblarr2D(:,:),DIMS=hyperdims(1:rank), &
                            HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                            FILE_SPACE_ID=dspace_ID)
          ELSE
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=intarr2D(:,:),DIMS=hyperdims(1:rank), &
                            HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                            FILE_SPACE_ID=dspace_ID)
          END IF
        CASE (3)
          IF (.NOT. integ) THEN
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=dblarr3D(:,:,:),DIMS=hyperdims(1:rank), &
                            HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                            FILE_SPACE_ID=dspace_ID)
          ELSE
            CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                            BUF=intarr3D(:,:,:),DIMS=hyperdims(1:rank), &
                            HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                            FILE_SPACE_ID=dspace_ID)
          END IF
        CASE (4)
          CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                          BUF=dblarr4D(:,:,:,:),DIMS=hyperdims(1:rank), &
                          HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                          FILE_SPACE_ID=dspace_ID)
        CASE (5)
          CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                          BUF=dblarr5D(:,:,:,:,:),DIMS=hyperdims(1:rank), &
                          HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                          FILE_SPACE_ID=dspace_ID)
        CASE (6)
          CALL H5DWRITE_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                          BUF=dblarr6D(:,:,:,:,:,:),DIMS=hyperdims(1:rank), &
                          HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                          FILE_SPACE_ID=dspace_ID)
      END SELECT writesub

    END IF fullorpartial


    ! 4) Add a signature attribute
    !-----------------------------
    ! Write an informative attribute for each data array. READ_HDF5 does not 
    ! read it. It is only here for archival purposes. It can be read in Unix
    ! with h5dump -A file.h5
    writribute: IF (.NOT. subset) THEN
      CALL H5TCOPY_F(TYPE_ID=h5t_native_character,NEW_TYPE_ID=atype_ID, &
                     HDFERR=error) ! define the type (character)
      CALL H5TSET_SIZE_F(TYPE_ID=atype_ID,SIZE=INT(textwid,hsize_t), &
                         HDFERR=error) ! string
      CALL H5SCREATE_SIMPLE_F(RANK=arank,DIMS=adims,SPACE_ID=aspace_ID, &
                              HDFERR=error) ! create attribute dataspace
      CALL H5ACREATE_F(LOC_ID=dset_ID,NAME=aname,TYPE_ID=atype_ID, &
                       SPACE_ID=aspace_ID,ATTR_ID=attr_ID,HDFERR=error) ! dims
      adata(1) = "Array written in the customized format of routines" &
               // " WRITE_HDF5/READ_HDF5"
      adata(2) = MERGE("Compressed with ZLIB DEFLATE level 6", &
                       "Uncompressed                        ",cmprss)
      adata(3) = "Generated in Fortran by "//programrunner//", on " &
                     //TRIMLR(TODAY())//"."
      CALL H5AWRITE_F(ATTR_ID=attr_ID,MEMTYPE_ID=atype_ID,BUF=adata(:), &
                      DIMS=adims(:),HDFERR=error) ! write attribute
    END IF writribute

    
    ! 5) Close the environments and the file
    !---------------------------------------
    IF (.NOT. subset) &
      CALL H5ACLOSE_F(ATTR_ID=attr_ID,HDFERR=error)  ! close attribute 
    IF (subset) CALL H5SCLOSE_F(SPACE_ID=mspace_ID,HDFERR=error) ! close memmory
    CALL H5DCLOSE_F(DSET_ID=dset_ID,HDFERR=error)    ! close dataset
    CALL H5SCLOSE_F(SPACE_ID=dspace_ID,HDFERR=error) ! close dataspace
    CALL H5FCLOSE_F(FILE_ID=file_ID,HDFERR=error)    ! close file access
    CALL H5CLOSE_F(ERROR=error)                      ! close Fortran interface

    ! Courtesy notification
    u(1) = ustd
    IF (PRESENT(unit)) u(2) = unit
    IF (PRESENT(verbose)) THEN ; verb = verbose ; ELSE ; verb = verbatim ; ENDIF
    courtesy: IF (verb) THEN
      DO i=MERGE(1,2,verb),MERGE(2,1,PRESENT(unit))
        IF (newfile) THEN
          WRITE(u(i),*) " - Array '"//TRIMLR(dsetname) &
                      //"' has been written in the new file:"
          WRITE(u(i),*) "   "//TRIMLR(file)//"."
        ELSE IF (subset) THEN
          WRITE(u(i),*) " - Array '"//TRIMLR(dsetname)//"' in the file:"
          WRITE(u(i),*) "   "//TRIMLR(file)//" has been modified."
        ELSE IF (appnd) THEN
          WRITE(u(i),*) " - Array '"//TRIMLR(dsetname) &
                        //"' has been added to the file:"
          WRITE(u(i),*) "   "//TRIMLR(file)//"."
        END IF
      END DO
    END IF courtesy

    !-------------------------------------------------------------------------

  END SUBROUTINE write_hdf5


  !=========================================================================
  ! CALL READ_HDF5(DBLARR{1-6}D or STRARR1D or INTARR{1-3}D,FILE="",NAME="", 
  !                IND1=[idim1_inf,idim1_sup], &
  !                IND2=[idim2_inf,idim2_sup],IND3=[idim3_inf,idim3_sup], &
  !                IND4=[idim4_inf,idim4_sup],IND5=[idim5_inf,idim5_sup], &
  !                IND6=[idim6_inf,idim6_sup],N1=,N2=,N3=,N4=,N5=,N6=)
  ! 
  !   Read the array to an HDF5 file.
  !=========================================================================

  SUBROUTINE read_hdf5 (dblarr1D,dblarr2D,dblarr3D,dblarr4D,dblarr5D,dblarr6D, &
                        strarr1D,intarr1D,intarr2D,intarr3D,file,name, &
                        ind1,ind2,ind3,ind4,ind5,ind6,N1,N2,N3,N4,N5,N6)

    USE utilities, ONLY: DP, trimlr, strike
    USE hdf5, ONLY: hid_t, hsize_t, h5f_acc_rdonly_f, h5t_native_double, &
                    h5t_native_integer, h5t_native_character, &
                    h5s_select_set_f, h5tcopy_f, &
                    h5tset_size_f, h5open_f, h5close_f, h5fopen_f, h5fclose_f, &
                    h5dopen_f, h5dread_f, h5dclose_f, h5dget_space_f, &
                    h5screate_f, h5sget_simple_extent_dims_f, &
                    h5sselect_hyperslab_f, h5screate_simple_f, h5sclose_f
    IMPLICIT NONE

    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr1D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr2D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr3D
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr4D
    REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(OUT), &
      OPTIONAL :: dblarr5D
    REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, INTENT(OUT), &
      OPTIONAL :: dblarr6D
    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: strarr1D
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: intarr1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: intarr2D
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: intarr3D
    CHARACTER(*), INTENT(IN) :: file
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: ind1, ind2, ind3, ind4, ind5
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: ind6
    INTEGER, INTENT(OUT), OPTIONAL :: N1, N2, N3, N4, N5, N6

    INTEGER(hsize_t), DIMENSION(rank0) :: dims      ! dimension of dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: maxdims   ! max. dimension of dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: hyperdims ! dimension of hyperslab
    INTEGER(hid_t) :: file_ID    ! file identifier
    INTEGER(hid_t) :: dspace_ID  ! dataspace identifier
    INTEGER(hid_t) :: dset_ID    ! dataset identifier
    INTEGER(hid_t) :: dtype_ID   ! data type identifier
    INTEGER(hid_t) :: mspace_ID  ! memmory identifier
    INTEGER(hsize_t), DIMENSION(rank0) :: offset   ! hyperslab offset
    INTEGER(hsize_t), DIMENSION(2) :: i1, i2, i3, i4 ,i5, i6 ! [min,max] indices
    INTEGER :: rank, strlen
    INTEGER :: error ! error flag
    LOGICAL :: subset, integ, string
    CHARACTER(textwid) :: dsetname

    !-------------------------------------------------------------------------

    ! 1) Keywords and preliminray settings
    !-------------------------------------
    ! Type of reading
    subset = ( PRESENT(ind1) .OR. PRESENT(ind2) .OR. PRESENT(ind3) &
               .OR. PRESENT(ind4) .OR. PRESENT(ind5) .OR. PRESENT(ind6) )

    ! Dataset name (used to identify the struture in the HDF5 file)
    IF (PRESENT(name)) THEN
      dsetname = name
    ELSE
      dsetname = dsetname_def
    END IF

    ! Actual rank
    rank = -1
    string = .False.
    integ = .False.
    IF (PRESENT(dblarr1D)) rank = 1
    IF (PRESENT(dblarr2D)) THEN
      IF (rank == -1) THEN
        rank = 2
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr3D)) THEN
      IF (rank == -1) THEN
        rank = 3
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr4D)) THEN
      IF (rank == -1) THEN
        rank = 4
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr5D)) THEN
      IF (rank == -1) THEN
        rank = 5
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(dblarr6D)) THEN
      IF (rank == -1) THEN
        rank = 6
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    strorDP: IF (PRESENT(strarr1D)) THEN
      IF (rank == -1) THEN
        rank = 1
        string = .True.
        strlen = LEN(strarr1D)
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF strorDP
    IF (PRESENT(intarr1D)) THEN
      IF (rank == -1) THEN
        rank = 1
        integ = .True.
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(intarr2D)) THEN
      IF (rank == -1) THEN
        rank = 2
        integ = .True.
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (PRESENT(intarr3D)) THEN
      IF (rank == -1) THEN
        rank = 3
        integ = .True.
      ELSE
        CALL STRIKE("READ_HDF5","you can select only one ARRAY shape")
      END IF
    END IF
    IF (rank == -1) &
      CALL STRIKE("READ_HDF5","you have to select at least one ARRAY")

    ! Consistency check
    IF (PRESENT(ind2) .AND. rank < 2) &
      CALL STRIKE("READ_HDF5","you cannot select IND2 if ARRAY is less than 2D")
    IF (PRESENT(N2) .AND. rank < 2) &
      CALL STRIKE("READ_HDF5","you cannot select N2 if ARRAY is less than 2D")
    IF (PRESENT(ind3) .AND. rank < 3) &
      CALL STRIKE("READ_HDF5","you cannot select IND3 if ARRAY is less than 3D")
    IF (PRESENT(N3) .AND. rank < 3) &
      CALL STRIKE("READ_HDF5","you cannot select N3 if ARRAY is less than 3D")
    IF (PRESENT(ind4) .AND. rank < 4) &
      CALL STRIKE("READ_HDF5","you cannot select IND4 if ARRAY is less than 4D")
    IF (PRESENT(N4) .AND. rank < 4) &
      CALL STRIKE("READ_HDF5","you cannot select N4 if ARRAY is less than 4D")
    IF (PRESENT(ind5) .AND. rank < 5) &
      CALL STRIKE("READ_HDF5","you cannot select IND5 if ARRAY is less than 5D")
    IF (PRESENT(N5) .AND. rank < 5) &
      CALL STRIKE("READ_HDF5","you cannot select N5 if ARRAY is less than 5D")

  
    ! 2) Open the file and the HDF5 environments
    !-------------------------------------------
    ! Initialize the HDF5 interface
    CALL H5OPEN_F(ERROR=error)

    ! Create new file
    CALL H5FOPEN_F(NAME=file,ACCESS_FLAGS=h5f_acc_rdonly_f,FILE_ID=file_ID, &
                   HDFERR=error)


    ! 3) Reading whole or partial datset
    !----------------------------------- 
    ! Define the type of the dataset (DP or string)
    strfullcompr: IF (.NOT. string) THEN
      dtype_ID = MERGE(h5t_native_integer,h5t_native_double,integ)
    ELSE
      CALL H5TCOPY_F(TYPE_ID=h5t_native_character,NEW_TYPE_ID=dtype_ID, &
                     HDFERR=error) ! define the type (character)
      CALL H5TSET_SIZE_F(TYPE_ID=dtype_ID,SIZE=INT(strlen,hsize_t), &
                         HDFERR=error) ! string
    END IF strfullcompr

    ! Get the dimenions
    CALL H5DOPEN_F(LOC_ID=file_ID,NAME=dsetname,DSET_ID=dset_id,HDFERR=error)
    CALL H5DGET_SPACE_F(DATASET_ID=dset_id,DATASPACE_ID=dspace_ID,HDFERR=error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(SPACE_ID=dspace_ID,DIMS=dims(1:rank), &
                                     MAXDIMS=maxdims(1:rank),HDFERR=error)

    ! Actual reading
    fullorpartial: IF (.NOT. subset) THEN

      ! Read the whole array
      fullread: SELECT CASE (rank)
        CASE (1)
          IF (.NOT. string) THEN
            IF (.NOT. integ) THEN
              ALLOCATE (dblarr1D(dims(1)))
              CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                             BUF=dblarr1D(:),DIMS=dims(1:rank),HDFERR=error)
            ELSE
              ALLOCATE (intarr1D(dims(1)))
              CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                             BUF=intarr1D(:),DIMS=dims(1:rank),HDFERR=error)
            END IF
          ELSE
            ALLOCATE (strarr1D(dims(1)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=strarr1D(:),DIMS=dims(1:rank),HDFERR=error)
          END IF
        CASE (2)
          IF (.NOT. integ) THEN
            ALLOCATE (dblarr2D(dims(1),dims(2)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=dblarr2D(:,:),DIMS=dims(1:rank),HDFERR=error)
          ELSE
            ALLOCATE (intarr2D(dims(1),dims(2)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=intarr2D(:,:),DIMS=dims(1:rank),HDFERR=error)
          END IF
        CASE (3)
          IF (.NOT. integ) THEN
            ALLOCATE (dblarr3D(dims(1),dims(2),dims(3)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=dblarr3D(:,:,:),DIMS=dims(1:rank),HDFERR=error)
          ELSE
            ALLOCATE (intarr3D(dims(1),dims(2),dims(3)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=intarr3D(:,:,:),DIMS=dims(1:rank),HDFERR=error)
          END IF
        CASE (4)
          ALLOCATE (dblarr4D(dims(1),dims(2),dims(3),dims(4)))
          CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                         BUF=dblarr4D(:,:,:,:),DIMS=dims(1:rank),HDFERR=error)
        CASE (5)
          ALLOCATE (dblarr5D(dims(1),dims(2),dims(3),dims(4),dims(5)))
          CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                         BUF=dblarr5D(:,:,:,:,:),DIMS=dims(1:rank),HDFERR=error)
        CASE (6)
          ALLOCATE (dblarr6D(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
          CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                        BUF=dblarr6D(:,:,:,:,:,:),DIMS=dims(1:rank),HDFERR=error)
      END SELECT fullread

    ELSE

      ! Read only a subset of an already existing dataset
      IF (PRESENT(ind1)) THEN 
        i1(:) = ind1(:) 
      ELSE 
        i1(:) = [ INT(1,hsize_t), dims(1) ]
      END IF
      IF (PRESENT(ind2)) THEN 
        i2(:) = ind2(:) 
      ELSE
        i2(:) = [ INT(1,hsize_t), dims(2) ]
      END IF
      IF (PRESENT(ind3)) THEN 
        i3(:) = ind3(:) 
      ELSE 
        i3(:) = [ INT(1,hsize_t), dims(3) ]
      END IF
      IF (PRESENT(ind4)) THEN 
        i4(:) = ind4(:) 
      ELSE 
        i4(:) = [ INT(1,hsize_t), dims(4) ]
      END IF
      IF (PRESENT(ind5)) THEN 
        i5(:) = ind5(:) 
      ELSE 
        i5(:) = [ INT(1,hsize_t), dims(5) ]
      END IF
      IF (PRESENT(ind6)) THEN 
        i6(:) = ind6(:) 
      ELSE 
        i6(:) = [ INT(1,hsize_t), dims(6) ]
      END IF
      offset(:) = [ i1(1), i2(1), i3(1), i4(1), i5(1), i6(1) ] - 1 ! offset
                                                                   ! starts at 0
      hyperdims(:) = [ i1(2)-i1(1)+1, i2(2)-i2(1)+1, i3(2)-i3(1)+1, &
                       i4(2)-i4(1)+1, i5(2)-i5(1)+1, i6(2)-i6(1)+1 ]
      CALL H5SSELECT_HYPERSLAB_F(SPACE_ID=dspace_ID,OPERATOR=h5s_select_set_f, &
                                 START=offset(1:rank),COUNT=hyperdims(1:rank), &
                                 HDFERR=error) ! select the subset
      CALL H5SCREATE_SIMPLE_F(RANK=rank,DIMS=hyperdims(1:rank), &
                              SPACE_ID=mspace_ID, &
                              HDFERR=error) ! create sub-dataspace
      subread: SELECT CASE (rank)
        CASE (1)
          IF (.NOT. string) THEN
            IF (.NOT. integ) THEN
              ALLOCATE (dblarr1D(hyperdims(1)))
              CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                             BUF=dblarr1D(:),DIMS=hyperdims(1:rank), &
                             HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                             FILE_SPACE_ID=dspace_ID)
            ELSE
              ALLOCATE (intarr1D(hyperdims(1)))
              CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                             BUF=intarr1D(:),DIMS=hyperdims(1:rank), &
                             HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                             FILE_SPACE_ID=dspace_ID)
            END IF
          ELSE
            ALLOCATE (strarr1D(hyperdims(1)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=strarr1D(:),DIMS=hyperdims(1:rank),HDFERR=error,&
                           MEM_SPACE_ID=mspace_ID,FILE_SPACE_ID=dspace_ID)
          END IF
        CASE (2)
          IF (.NOT. integ) THEN
            ALLOCATE (dblarr2D(hyperdims(1),hyperdims(2)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=dblarr2D(:,:),DIMS=hyperdims(1:rank), &
                           HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                           FILE_SPACE_ID=dspace_ID)
          ELSE
            ALLOCATE (intarr2D(hyperdims(1),hyperdims(2)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=intarr2D(:,:),DIMS=hyperdims(1:rank), &
                           HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                           FILE_SPACE_ID=dspace_ID)
          END IF
        CASE (3)
          IF (.NOT. integ) THEN
            ALLOCATE (dblarr3D(hyperdims(1),hyperdims(2),hyperdims(3)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=dblarr3D(:,:,:),DIMS=hyperdims(1:rank), &
                           HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                           FILE_SPACE_ID=dspace_ID)
          ELSE
            ALLOCATE (intarr3D(hyperdims(1),hyperdims(2),hyperdims(3)))
            CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                           BUF=intarr3D(:,:,:),DIMS=hyperdims(1:rank), &
                           HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                           FILE_SPACE_ID=dspace_ID)
          END IF
        CASE (4)
          ALLOCATE (dblarr4D(hyperdims(1),hyperdims(2),hyperdims(3), &
                    hyperdims(4)))
          CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                         BUF=dblarr4D(:,:,:,:),DIMS=hyperdims(1:rank), &
                         HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                         FILE_SPACE_ID=dspace_ID)
        CASE (5)
          ALLOCATE (dblarr5D(hyperdims(1),hyperdims(2),hyperdims(3), &
                             hyperdims(4),hyperdims(5)))
          CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                         BUF=dblarr5D(:,:,:,:,:),DIMS=hyperdims(1:rank), &
                         HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                         FILE_SPACE_ID=dspace_ID)
        CASE (6)
          ALLOCATE (dblarr6D(hyperdims(1),hyperdims(2),hyperdims(3), &
                             hyperdims(4),hyperdims(5),hyperdims(6)))
          CALL H5DREAD_F(DSET_ID=dset_ID,MEM_TYPE_ID=dtype_ID, &
                         BUF=dblarr6D(:,:,:,:,:,:),DIMS=hyperdims(1:rank), &
                         HDFERR=error,MEM_SPACE_ID=mspace_ID, &
                         FILE_SPACE_ID=dspace_ID)
      END SELECT subread

    END IF fullorpartial


    ! 4) Close the environments and the file
    !---------------------------------------
    IF (subset) CALL H5SCLOSE_F(SPACE_ID=mspace_ID,HDFERR=error) ! close memmory
    CALL H5SCLOSE_F(SPACE_ID=dspace_ID,HDFERR=error) ! close dataspace
    CALL H5DCLOSE_F(DSET_ID=dset_ID,HDFERR=error)    ! close dataset
    CALL H5FCLOSE_F(FILE_ID=file_ID,HDFERR=error)    ! close file access
    CALL H5CLOSE_F(ERROR=error)                      ! close Fortran interface

    ! Optional output
    IF (PRESENT(N1)) N1 = INT(dims(1),KIND(N1))
    IF (PRESENT(N2)) N2 = INT(dims(2),KIND(N2))
    IF (PRESENT(N3)) N3 = INT(dims(3),KIND(N3))
    IF (PRESENT(N4)) N4 = INT(dims(4),KIND(N4))
    IF (PRESENT(N5)) N5 = INT(dims(5),KIND(N5))
    IF (PRESENT(N6)) N6 = INT(dims(6),KIND(N6))
    
    !-------------------------------------------------------------------------

  END SUBROUTINE read_hdf5


  !=========================================================================
  ! CALL GETDIM_HDF5(FILE="",NAME="",N1=,N2=,N3=,N4=,N5=)
  ! 
  !   Read the dimensions of an array in an HDF5 file.
  !=========================================================================

  SUBROUTINE getdim_hdf5 (file,name,N1,N2,N3,N4,N5,N6)

    USE utilities, ONLY: DP, trimlr, strike
    USE hdf5, ONLY: hid_t, hsize_t, h5f_acc_rdonly_f, h5open_f, h5close_f, &
                    h5fopen_f, h5fclose_f, h5dopen_f, h5dread_f, h5dclose_f, &
                    h5dget_space_f, h5sget_simple_extent_dims_f, &
                    h5sget_simple_extent_ndims_f, h5sclose_f
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: file
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    INTEGER, INTENT(OUT), OPTIONAL :: N1, N2, N3, N4, N5, N6

    INTEGER(hsize_t), DIMENSION(rank0) :: dims      ! dimension of dataset
    INTEGER(hsize_t), DIMENSION(rank0) :: maxdims   ! max. dimension of dataset
    INTEGER(hid_t) :: file_ID    ! file identifier
    INTEGER(hid_t) :: dspace_ID  ! dataspace identifier
    INTEGER(hid_t) :: dset_ID    ! dataset identifier
    INTEGER :: rank, error

    !-------------------------------------------------------------------------

    ! Initialize the HDF5 interface
    CALL H5OPEN_F(ERROR=error)

    ! Open the file
    CALL H5FOPEN_F(NAME=file,ACCESS_FLAGS=h5f_acc_rdonly_f,FILE_ID=file_ID, &
                   HDFERR=error)

    ! Get the dimenions
    CALL H5DOPEN_F(LOC_ID=file_ID,NAME=name,DSET_ID=dset_id,HDFERR=error)
    CALL H5DGET_SPACE_F(DATASET_ID=dset_id,DATASPACE_ID=dspace_ID,HDFERR=error)
    CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(SPACE_ID=dspace_ID,RANK=rank,HDFERR=error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(SPACE_ID=dspace_ID,DIMS=dims(1:rank), &
                                     MAXDIMS=maxdims(1:rank),HDFERR=error)

    ! Closing
    CALL H5SCLOSE_F(SPACE_ID=dspace_ID,HDFERR=error) ! close dataspace
    CALL H5DCLOSE_F(DSET_ID=dset_ID,HDFERR=error)    ! close dataset
    CALL H5FCLOSE_F(FILE_ID=file_ID,HDFERR=error)    ! close file access
    CALL H5CLOSE_F(ERROR=error)                      ! close Fortran interface

    ! Optional output
    IF (PRESENT(N1)) N1 = INT(dims(1),KIND(N1))
    IF (PRESENT(N2)) N2 = INT(dims(2),KIND(N2))
    IF (PRESENT(N3)) N3 = INT(dims(3),KIND(N3))
    IF (PRESENT(N4)) N4 = INT(dims(4),KIND(N4))
    IF (PRESENT(N5)) N5 = INT(dims(5),KIND(N5))
    IF (PRESENT(N6)) N6 = INT(dims(6),KIND(N6))

    !-------------------------------------------------------------------------

  END SUBROUTINE getdim_hdf5

  
  !==========================================================================
  ! bool = CHECK_HDF5(file)
  !
  !    Returns F is the file is corrupted.
  !==========================================================================

  FUNCTION check_hdf5 (file)

    USE hdf5, ONLY: h5fis_hdf5_f
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: file
    LOGICAL :: check_hdf5

    INTEGER :: hdferr
    LOGICAL :: exist
    
    !-------------------------------------------------------------------------

    INQUIRE (FILE=file,EXIST=exist)
    IF (exist) THEN
      CALL H5FIS_HDF5_F(NAME=file,STATUS=check_hdf5,HDFERR=hdferr)
    ELSE
      check_hdf5 = .False.
    END IF
    
    !-------------------------------------------------------------------------
    
  END FUNCTION check_hdf5
  
  
  !==========================================================================
  ! CALL IND_HDF5_FRAC(N,Nfile,IFIRST=[Nfile],ILAST=[Nfile])
  !  
  !   Simple function to set the list of first and last indices of fractionated
  ! files, to be used by WRITE_HDF5_FRAC and READ_HDF5_FRAC.
  !==========================================================================

  SUBROUTINE ind_hdf5_frac (N,Nfile,ifirst,ilast)

    USE utilities, ONLY: strike
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: N, Nfile
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ifirst, ilast
    
    INTEGER :: i
    
    !-------------------------------------------------------------------------

    IF (Nfile > 0) THEN
      ALLOCATE (ifirst(Nfile),ilast(Nfile))
      IF (N < Nfile) & 
        CALL STRIKE("IND_HDF5_FRAC", &
                "the number of files must be smaller than the size of the array")
      DO i=1,Nfile
        ifirst(i) = NINT( (REAL(N)/Nfile) * (i-1) + 1 )
        ilast(i) = NINT( (REAL(N)/Nfile) * i )
        IF (i == Nfile) ilast(i) = N
      END DO
    ELSE
      ALLOCATE (ifirst(1),ilast(1))
      ifirst(1) = 1
      ilast(1) = N
    END IF
      
    !-------------------------------------------------------------------------

  END SUBROUTINE ind_hdf5_frac

  
  !==========================================================================
  ! CALL WRITE_HDF5_FRAC(IFIRST=[Nfile],ILAST=[Nfile], &
  !                      DBLARR{1-5}D or STRARR1D or INTARR{1-3}D,FILE=[Nfile], &
  !                      COMPRESS=T/F,INITDBLARR=[],INITINTARR=[], &
  !                      APPEND=T/F,NAME=[""],IND1=[idim1_inf,idim1_sup], &
  !                      IND2=[idim2_inf,idm2_sup],IND3=[idim3_inf,idm3_sup], &
  !                      IND4=[idim4_inf,idm4_sup],IND5=[idim5_inf,idm5_sup], &
  !                      UNIT=)
  !
  !    Interface of WRITE_HDF5, but when the array is spread over several files.
  ! This is useful when long runs are concerned as a power cut may corrupt the
  ! file.
  !    FILE has to be a list of files. IFIRST and ILAST are the first and last
  ! indices of the array in each file (should be ordered: IFIRST(i)<=ILAST(i) and
  ! IFIRST(i)=ILAST(i-1)+1; use IND_HDF5_FRAC to set them). It is assumed that
  ! the dimension to be fractionated is the last dimension the array.
  !    An array written with this function must have been initialized first.
  !==========================================================================

  SUBROUTINE write_hdf5_frac (ifirst,ilast,dblarr1D,dblarr2D,dblarr3D,dblarr4D, &
                              dblarr5D,strarr1D,intarr1D,intarr2D,intarr3D, &
                              initdblarr,initintarr,file,name,compress,append, &
                              ind1,ind2,ind3,ind4,ind5,unit,verbose)

    USE utilities, ONLY: DP
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) :: ifirst, ilast
    REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: dblarr1D
    REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: dblarr2D
    REAL(DP), DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: dblarr3D
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN), OPTIONAL :: dblarr4D
    REAL(DP), DIMENSION(:,:,:,:,:), INTENT(IN), OPTIONAL :: dblarr5D
    CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: strarr1D
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: intarr1D
    INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: intarr2D
    INTEGER, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: intarr3D
    INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: initdblarr, initintarr
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: file
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    LOGICAL, INTENT(IN), OPTIONAL :: compress, append
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: ind1, ind2, ind3, ind4, ind5
    INTEGER, INTENT(IN), OPTIONAL :: unit
    LOGICAL, INTENT(IN), OPTIONAL :: verbose

    INTEGER :: i, ifile, Ndim, Nfile, i0arr
    INTEGER, DIMENSION(2) :: ranfile
    INTEGER, DIMENSION(rank0) :: dimen
    INTEGER, DIMENSION(rank0,2) :: ind2D
    INTEGER, DIMENSION(SIZE(file),rank0,2) :: ind
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwh
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: initarr
    LOGICAL :: init
    
    !------------------------------------------------------------------------- 
    
    ! 1) Prepare the variables
    !-------------------------
    ! Set the dimensions of the array
    Ndim = 0
    Nfile = SIZE(file(:))
    ind2D(:,:) = 0
    IF (PRESENT(initdblarr)) THEN
      Ndim = SIZE(initdblarr(:))
      FORALL (i=1:Ndim)
        dimen(i) = initdblarr(i)
        ind2D(i,:) = [1,dimen(i)]
      END FORALL
      init = .True.
    ENDIF
    IF (PRESENT(initintarr)) THEN
      Ndim = SIZE(initintarr(:))
      FORALL (i=1:Ndim)
        dimen(i) = initintarr(i)
        ind2D(i,:) = [1,dimen(i)]
      END FORALL
      init = .True.
    ENDIF
    IF (PRESENT(dblarr1D) .OR. PRESENT(intarr1D) .OR. PRESENT(strarr1D)) THEN
      Ndim = 1
      IF (PRESENT(dblarr1D)) dimen(1) = SIZE(dblarr1D(:))
      IF (PRESENT(intarr1D)) dimen(1) = SIZE(intarr1D(:))
      IF (PRESENT(strarr1D)) dimen(1) = SIZE(strarr1D(:))
      dimen(2:) = 0
      IF (PRESENT(ind1)) THEN ; ind2D(1,:) = ind1(:)
                         ELSE ; ind2D(1,:) = [1,dimen(1)] ; END IF
      init = .False.
    END IF
    IF (PRESENT(dblarr2D) .OR. PRESENT(intarr2D)) THEN
      Ndim = 2
      IF (PRESENT(dblarr2D)) &
        FORALL (i=1:Ndim) dimen(i) = SIZE(dblarr2D(:,:),i)
      IF (PRESENT(intarr2D)) &
        FORALL (i=1:Ndim) dimen(i) = SIZE(intarr2D(:,:),i)
      dimen(3:) = 0
      IF (PRESENT(ind1)) THEN ; ind2D(1,:) = ind1(:)
                         ELSE ; ind2D(1,:) = [1,dimen(1)] ; END IF
      IF (PRESENT(ind2)) THEN ; ind2D(2,:) = ind2(:)
                         ELSE ; ind2D(2,:) = [1,dimen(2)] ; END IF
      init = .False.
    END IF
    IF (PRESENT(dblarr3D) .OR. PRESENT(intarr3D)) THEN
      Ndim = 3
      IF (PRESENT(dblarr3D)) &
        FORALL (i=1:Ndim) dimen(i) = SIZE(dblarr3D(:,:,:),i)
      IF (PRESENT(intarr3D)) &
        FORALL (i=1:Ndim) dimen(i) = SIZE(intarr3D(:,:,:),i)
      dimen(4:) = 0
      IF (PRESENT(ind1)) THEN ; ind2D(1,:) = ind1(:)
                         ELSE ; ind2D(1,:) = [1,dimen(1)] ; END IF
      IF (PRESENT(ind2)) THEN ; ind2D(2,:) = ind2(:)
                         ELSE ; ind2D(2,:) = [1,dimen(2)] ; END IF
      IF (PRESENT(ind3)) THEN ; ind2D(3,:) = ind3(:)
                         ELSE ; ind2D(3,:) = [1,dimen(3)] ; END IF
      init = .False.
    END IF
    IF (PRESENT(dblarr4D)) THEN
      Ndim = 4
      IF (PRESENT(dblarr4D)) &
        FORALL (i=1:Ndim) dimen(i) = SIZE(dblarr4D(:,:,:,:),i)
      dimen(5:) = 0
      IF (PRESENT(ind1)) THEN ; ind2D(1,:) = ind1(:)
                         ELSE ; ind2D(1,:) = [1,dimen(1)] ; END IF
      IF (PRESENT(ind2)) THEN ; ind2D(2,:) = ind2(:)
                         ELSE ; ind2D(2,:) = [1,dimen(2)] ; END IF
      IF (PRESENT(ind3)) THEN ; ind2D(3,:) = ind3(:)
                         ELSE ; ind2D(3,:) = [1,dimen(3)] ; END IF
      IF (PRESENT(ind4)) THEN ; ind2D(4,:) = ind4(:)
                         ELSE ; ind2D(4,:) = [1,dimen(4)] ; END IF
      init = .False.
    END IF
    IF (PRESENT(dblarr5D)) THEN
      Ndim = 5
      IF (PRESENT(dblarr5D)) &
        FORALL (i=1:Ndim) dimen(i) = SIZE(dblarr5D(:,:,:,:,:),i)
      IF (PRESENT(ind1)) THEN ; ind2D(1,:) = ind1(:)
                         ELSE ; ind2D(1,:) = [1,dimen(1)] ; END IF
      IF (PRESENT(ind2)) THEN ; ind2D(2,:) = ind2(:)
                         ELSE ; ind2D(2,:) = [1,dimen(2)] ; END IF
      IF (PRESENT(ind3)) THEN ; ind2D(3,:) = ind3(:)
                         ELSE ; ind2D(3,:) = [1,dimen(3)] ; END IF
      IF (PRESENT(ind4)) THEN ; ind2D(4,:) = ind4(:)
                         ELSE ; ind2D(4,:) = [1,dimen(4)] ; END IF
      IF (PRESENT(ind5)) THEN ; ind2D(5,:) = ind5(:)
                         ELSE ; ind2D(5,:) = [1,dimen(5)] ; END IF
      init = .False.
    END IF

    ! Set the indices for writing the different files
    ind(:,:,:) = 0
    CALL IWHERE(ind2D(Ndim,1) >= ifirst(:),iwh)
    ranfile(1) = MAXVAL(iwh(:))
    CALL IWHERE(ind2D(Ndim,2) <= ilast(:),iwh)
    ranfile(2) = MINVAL(iwh(:))
    DO ifile=ranfile(1),ranfile(2)
      IF (Ndim > 1) ind(ifile,1:Ndim-1,:) = ind2D(1:Ndim-1,:)
      ind(ifile,Ndim,1) = MERGE(ind2D(Ndim,1),ifirst(ifile),ifile == ranfile(1))
      ind(ifile,Ndim,2) = MERGE(ind2D(Ndim,2),ilast(ifile),ifile == ranfile(2))
    END DO


    ! 2) Write the files
    !-------------------
    writarr: IF (init) THEN

      ! Size of the array
      ALLOCATE (initarr(Nfile,Ndim))
      IF (Ndim > 1) &
        FORALL (ifile=1:Nfile) initarr(ifile,1:Ndim-1) = dimen(1:Ndim-1)
      initarr(:,Ndim) = ilast(:) - ifirst(:) + 1
      
      ! Initialize
      IF (PRESENT(initdblarr)) THEN
        DO ifile=1,Nfile
          CALL WRITE_HDF5(INITDBLARR=initarr(ifile,:), &
                          FILE=file(ifile),NAME=name,COMPRESS=compress, &
                          APPEND=append,UNIT=unit,VERBOSE=verbose)
        END DO
      ELSE IF (PRESENT(initintarr)) THEN
        DO ifile=1,Nfile
          CALL WRITE_HDF5(INITINTARR=initarr(ifile,:), &
                          FILE=file(ifile),NAME=name,COMPRESS=compress, &
                          APPEND=append,UNIT=unit,VERBOSE=verbose)          
        END DO        
      END IF
      
    ELSE

      ! Write the array
      i0arr = ind(ranfile(1),Ndim,1) - 1
      IF (PRESENT(dblarr1D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(DBLARR1D=dblarr1D(ind(ifile,1,1)-i0arr: &
                                            ind(ifile,1,2)-i0arr), &
                          IND1=ind(ifile,1,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(dblarr2D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(DBLARR2D=dblarr2D(ind(ifile,1,1):ind(ifile,1,2), &
                                            ind(ifile,2,1)-i0arr: &
                                                         ind(ifile,2,2)-i0arr), &
                          IND1=ind(ifile,1,:), &
                          IND2=ind(ifile,2,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(dblarr3D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(DBLARR3D=dblarr3D(ind(ifile,1,1):ind(ifile,1,2), &
                                            ind(ifile,2,1):ind(ifile,2,2), &
                                            ind(ifile,3,1)-i0arr: &
                                                         ind(ifile,3,2)-i0arr), &
                          IND1=ind(ifile,1,:), &
                          IND2=ind(ifile,2,:), &
                          IND3=ind(ifile,3,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(dblarr4D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(DBLARR4D=dblarr4D(ind(ifile,1,1):ind(ifile,1,2), &
                                            ind(ifile,2,1):ind(ifile,2,2), &
                                            ind(ifile,3,1):ind(ifile,3,2), &
                                            ind(ifile,4,1)-i0arr: &
                                                         ind(ifile,4,2)-i0arr), &
                          IND1=ind(ifile,1,:), &
                          IND2=ind(ifile,2,:), &
                          IND3=ind(ifile,3,:), &
                          IND4=ind(ifile,4,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(dblarr5D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(DBLARR5D=dblarr5D(ind(ifile,1,1):ind(ifile,1,2), &
                                            ind(ifile,2,1):ind(ifile,2,2), &
                                            ind(ifile,3,1):ind(ifile,3,2), &
                                            ind(ifile,4,1):ind(ifile,4,2), &
                                            ind(ifile,5,1)-i0arr: &
                                                         ind(ifile,5,2)-i0arr), &
                          IND1=ind(ifile,1,:), &
                          IND2=ind(ifile,2,:), &
                          IND3=ind(ifile,3,:), &
                          IND4=ind(ifile,4,:), &
                          IND5=ind(ifile,5,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(intarr1D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(INTARR1D=intarr1D(ind(ifile,1,1)-i0arr: &
                                            ind(ifile,1,2)-i0arr), &
                          IND1=ind(ifile,1,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(intarr2D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(INTARR2D=intarr2D(ind(ifile,1,1):ind(ifile,1,2), &
                                            ind(ifile,2,1)-i0arr: &
                                                         ind(ifile,2,2)-i0arr), &
                          IND1=ind(ifile,1,:), &
                          IND2=ind(ifile,2,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(intarr3D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(INTARR3D=intarr3D(ind(ifile,1,1):ind(ifile,1,2), &
                                            ind(ifile,2,1):ind(ifile,2,2), &
                                            ind(ifile,3,1)-i0arr: &
                                                         ind(ifile,3,2)-i0arr), &
                          IND1=ind(ifile,1,:), &
                          IND2=ind(ifile,2,:), &
                          IND3=ind(ifile,3,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      IF (PRESENT(strarr1D)) THEN
        DO ifile=ranfile(1),ranfile(2)
          CALL WRITE_HDF5(STRARR1D=strarr1D(ind(ifile,1,1)-i0arr: &
                                            ind(ifile,1,2)-i0arr), &
                          IND1=ind(ifile,1,:)-ifirst(ifile)+1, &
                          FILE=file(ifile),NAME=name,UNIT=unit, &
                          COMPRESS=compress,VERBOSE=verbose)
        END DO
      END IF
      
    END IF writarr
  
    !-------------------------------------------------------------------------

  END SUBROUTINE write_hdf5_frac
    

  !=========================================================================
  ! CALL READ_HDF5_FRAC(DBLARR{1-5}D or STRARR1D or INTARR{1-3}D,FILE=[Nfile], &
  !                     IFIRST=[Nfile],ILAST=[Nfile],NAME="", &
  !                     IND1=[idim1_inf,idim1_sup], &
  !                     IND2=[idim2_inf,idm2_sup],IND3=[idim3_inf,idm3_sup], &
  !                     IND4=[idim4_inf,idm4_sup],IND5=[idim5_inf,idm5_sup], &
  !                     N1=,N2=,N3=,N4=,N5=)
  ! 
  !   
  !=========================================================================

  SUBROUTINE read_hdf5_frac (ifirst,ilast,dblarr1D,dblarr2D,dblarr3D,dblarr4D, &
                             dblarr5D,strarr1D,intarr1D,intarr2D,intarr3D,file, &
                             name,ind1,ind2,ind3,ind4,ind5,N1,N2,N3,N4,N5)

    USE utilities, ONLY: DP
    USE arrays, ONLY: iwhere
    IMPLICIT NONE

    INTEGER, DIMENSION(:), INTENT(IN) :: ifirst, ilast
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr1D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr2D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr3D
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: dblarr4D
    REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL::dblarr5D
    CHARACTER(*), DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: strarr1D
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: intarr1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: intarr2D
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: intarr3D
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: file
    CHARACTER(*), INTENT(IN), OPTIONAL :: name
    INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: ind1, ind2, ind3, ind4, ind5
    INTEGER, INTENT(OUT), OPTIONAL :: N1, N2, N3, N4, N5

    INTEGER :: ifile, Ndim, Nfile, i0arr, Nfrac, N10, N20, N30, N40, N50
    INTEGER :: Nfracsub
    INTEGER, DIMENSION(2) :: ranfile
    INTEGER, DIMENSION(rank0,2) :: ind2D
    INTEGER, DIMENSION(SIZE(file),rank0,2) :: ind
    INTEGER, DIMENSION(:), ALLOCATABLE :: iwh
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dbltmp1D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dbltmp2D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dbltmp3D
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: dbltmp4D
    REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: dbltmp5D
    CHARACTER(lenline), DIMENSION(:), ALLOCATABLE :: strtmp1D
    INTEGER, DIMENSION(:), ALLOCATABLE :: inttmp1D
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: inttmp2D
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: inttmp3D
    
    !------------------------------------------------------------------------- 
    
    ! 1) Prepare the variables
    !-------------------------
    ! Dimensions of the fields
    CALL GETDIM_HDF5(FILE=file(1),NAME=name,N1=N10,N2=N20,N3=N30,N4=N40)
    N50 = 0
    
    ! Set the dimensions of the array
    Nfile = SIZE(file(:))
    Nfrac = ilast(Nfile) - ifirst(1) + 1
    ind2D(:,:) = 0
    IF (PRESENT(dblarr1D) .OR. PRESENT(intarr1D) .OR. PRESENT(strarr1D)) THEN
      Ndim = 1
      IF (PRESENT(ind2)) THEN ; ind2D(Ndim,:) = ind1(:)
                         ELSE ; ind2D(Ndim,:) = [1,Nfrac] ; END IF
    END IF
    IF (PRESENT(dblarr2D) .OR. PRESENT(intarr2D)) THEN
      Ndim = 2
      IF (PRESENT(ind1)) THEN
        ind2D(1,:) = ind1(:)
        N10 = ind1(2) - ind1(1) + 1
      ELSE
        ind2D(1,:) = [1,N10]
      END IF
      IF (PRESENT(ind2)) THEN ; ind2D(Ndim,:) = ind2(:)
                         ELSE ; ind2D(Ndim,:) = [1,Nfrac] ; END IF
    END IF
    IF (PRESENT(dblarr3D) .OR. PRESENT(intarr3D)) THEN
      Ndim = 3
      IF (PRESENT(ind1)) THEN
        ind2D(1,:) = ind1(:)
        N10 = ind1(2) - ind1(1) + 1
      ELSE
        ind2D(1,:) = [1,N10]
      END IF
      IF (PRESENT(ind2)) THEN
        ind2D(2,:) = ind2(:)
        N20 = ind2(2) - ind2(1) + 1
      ELSE
        ind2D(2,:) = [1,N20]
      END IF
      IF (PRESENT(ind3)) THEN ; ind2D(Ndim,:) = ind3(:)
                         ELSE ; ind2D(Ndim,:) = [1,Nfrac] ; END IF
    END IF
    IF (PRESENT(dblarr4D)) THEN
      Ndim = 4
      IF (PRESENT(ind1)) THEN
        ind2D(1,:) = ind1(:)
        N10 = ind1(2) - ind1(1) + 1
      ELSE
        ind2D(1,:) = [1,N10]
      END IF
      IF (PRESENT(ind2)) THEN
        ind2D(2,:) = ind2(:)
        N20 = ind2(2) - ind2(1) + 1
      ELSE
        ind2D(2,:) = [1,N20]
      END IF
      IF (PRESENT(ind3)) THEN
        ind2D(3,:) = ind3(:)
        N30 = ind3(2) - ind3(1) + 1
      ELSE
        ind2D(3,:) = [1,N30]
      END IF
      IF (PRESENT(ind4)) THEN ; ind2D(Ndim,:) = ind4(:)
                         ELSE ; ind2D(Ndim,:) = [1,Nfrac] ; END IF
    END IF
    IF (PRESENT(dblarr5D)) THEN
      Ndim = 5
      IF (PRESENT(ind1)) THEN
        ind2D(1,:) = ind1(:)
        N10 = ind1(2) - ind1(1) + 1
      ELSE
        ind2D(1,:) = [1,N10]
      END IF
      IF (PRESENT(ind2)) THEN
        ind2D(2,:) = ind2(:)
        N20 = ind2(2) - ind2(1) + 1
      ELSE
        ind2D(2,:) = [1,N20]
      END IF
      IF (PRESENT(ind3)) THEN
        ind2D(3,:) = ind3(:)
        N30 = ind3(2) - ind3(1) + 1
      ELSE
        ind2D(3,:) = [1,N30]
      END IF
      IF (PRESENT(ind4)) THEN
        ind2D(4,:) = ind4(:)
        N40 = ind4(2) - ind4(1) + 1
      ELSE
        ind2D(4,:) = [1,N40]
      END IF
      IF (PRESENT(ind5)) THEN ; ind2D(Ndim,:) = ind5(:)
                         ELSE ; ind2D(Ndim,:) = [1,Nfrac] ; END IF
    END IF
    Nfracsub = ind2D(Ndim,2) - ind2D(Ndim,1) + 1
    CALL IWHERE(ind2D(Ndim,1) >= ifirst(:),iwh)
    ranfile(1) = MAXVAL(iwh(:))
    CALL IWHERE(ind2D(Ndim,2) <= ilast(:),iwh)
    ranfile(2) = MINVAL(iwh(:))
    ind(:,:,:) = 0
    DO ifile=ranfile(1),ranfile(2)
      IF (Ndim > 1) ind(ifile,1:Ndim-1,:) = ind2D(1:Ndim-1,:)
      ind(ifile,Ndim,1) = MERGE(ind2D(Ndim,1),ifirst(ifile),ifile == ranfile(1))
      ind(ifile,Ndim,2) = MERGE(ind2D(Ndim,2),ilast(ifile),ifile == ranfile(2))
    END DO


    ! 2) Read the files
    !-------------------
    i0arr = ind(ranfile(1),Ndim,1) - 1
    IF (PRESENT(dblarr1D)) THEN
      ALLOCATE (dblarr1D(Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(DBLARR1D=dbltmp1D, &
                       IND1=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        dblarr1D(ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) = dbltmp1D(:)
      END DO
    END IF
    IF (PRESENT(dblarr2D)) THEN
      ALLOCATE (dblarr2D(N10,Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(DBLARR2D=dbltmp2D, &
                       IND1=ind(ifile,1,:), &
                       IND2=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        dblarr2D(:,ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) &
          = dbltmp2D(:,:)
      END DO
    END IF
    IF (PRESENT(dblarr3D)) THEN
      ALLOCATE (dblarr3D(N10,N20,Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(DBLARR3D=dbltmp3D, &
                       IND1=ind(ifile,1,:), &
                       IND2=ind(ifile,2,:), &
                       IND3=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        dblarr3D(:,:,ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) &
          = dbltmp3D(:,:,:)
      END DO
    END IF
    IF (PRESENT(dblarr4D)) THEN
      ALLOCATE (dblarr4D(N10,N20,N30,Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(DBLARR4D=dbltmp4D, &
                       IND1=ind(ifile,1,:), &
                       IND2=ind(ifile,2,:), &
                       IND3=ind(ifile,3,:), &
                       IND4=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        dblarr4D(:,:,:,ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) &
          = dbltmp4D(:,:,:,:)
      END DO
    END IF
    IF (PRESENT(dblarr5D)) THEN
      ALLOCATE (dblarr5D(N10,N20,N30,N40,Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(DBLARR5D=dbltmp5D, &
                       IND1=ind(ifile,1,:), &
                       IND2=ind(ifile,2,:), &
                       IND3=ind(ifile,3,:), &
                       IND4=ind(ifile,4,:), &
                       IND5=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        dblarr5D(:,:,:,:,ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) &
          = dbltmp5D(:,:,:,:,:)
      END DO
    END IF
    IF (PRESENT(intarr1D)) THEN
      ALLOCATE (intarr1D(Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(INTARR1D=inttmp1D, &
                       IND1=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        intarr1D(ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) = inttmp1D(:)
      END DO
    END IF
    IF (PRESENT(intarr2D)) THEN
      ALLOCATE (intarr2D(N10,Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(INTARR2D=inttmp2D, &
                       IND1=ind(ifile,1,:), &
                       IND2=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        intarr2D(:,ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) &
          = inttmp2D(:,:)
      END DO
    END IF
    IF (PRESENT(intarr3D)) THEN
      ALLOCATE (intarr3D(N10,N20,Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(INTARR3D=inttmp3D, &
                       IND1=ind(ifile,1,:), &
                       IND2=ind(ifile,2,:), &
                       IND3=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        intarr3D(:,:,ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) &
          = inttmp3D(:,:,:)
      END DO
    END IF
    IF (PRESENT(strarr1D)) THEN
      ALLOCATE (strarr1D(Nfracsub))
      DO ifile=ranfile(1),ranfile(2)
        CALL READ_HDF5(STRARR1D=strtmp1D, &
                       IND1=ind(ifile,Ndim,:)-ifirst(ifile)+1, &
                       FILE=file(ifile),NAME=name)
        strarr1D(ind(ifile,Ndim,1)-i0arr:ind(ifile,Ndim,2)-i0arr) = strtmp1D(:)
      END DO
    END IF


    ! 3) Eventually returns the sizes
    !--------------------------------
    SELECT CASE (Ndim)
      CASE (1)
        N10 = Nfrac
      CASE (2)
        N20 = Nfrac
      CASE (3)
        N30 = Nfrac
      CASE (4)
        N40 = Nfrac
      CASE (5)
        N50 = Nfrac
    END SELECT
    IF (PRESENT(N1)) N1 = N10
    IF (PRESENT(N2)) N2 = N20
    IF (PRESENT(N3)) N3 = N30
    IF (PRESENT(N4)) N4 = N40
    IF (PRESENT(N5)) N5 = N50

    ! Deallocate the arrays
    IF (ALLOCATED(dbltmp1D)) DEALLOCATE(dbltmp1D)
    IF (ALLOCATED(dbltmp2D)) DEALLOCATE(dbltmp2D)
    IF (ALLOCATED(dbltmp3D)) DEALLOCATE(dbltmp3D)
    IF (ALLOCATED(dbltmp4D)) DEALLOCATE(dbltmp4D)
    IF (ALLOCATED(dbltmp5D)) DEALLOCATE(dbltmp5D)
    IF (ALLOCATED(inttmp1D)) DEALLOCATE(inttmp1D)
    IF (ALLOCATED(inttmp2D)) DEALLOCATE(inttmp2D)
    IF (ALLOCATED(inttmp3D)) DEALLOCATE(inttmp3D)
    IF (ALLOCATED(strtmp1D)) DEALLOCATE(strtmp1D)
    
    !-------------------------------------------------------------------------

  END SUBROUTINE read_hdf5_frac


  !==========================================================================
  ! CALL WRITE_ASCII (file,vec1[N],vec2[N],vec3[N],vec4[N],vec5[N], &
  !                  mat1[N,M],mat2[N,M],mat3[N,M],mat4[N,M], &
  !                  comgen="",comvec1="",comvec2="",comvec3="",comvec4="", &
  !                  comvec5="",commat1="",commat2="",commat3="",commat4="", &
  !                  IDL=T/F)
  !
  !   Writes my customised files with a header and data blocks.
  !==========================================================================

  SUBROUTINE write_ascii (file,vec1,vec2,vec3,vec4,vec5,mat1,mat2,mat3,mat4, &
                         comgen,comvec1,comvec2,comvec3,comvec4,comvec5, &
                         commat1,commat2,commat3,commat4,idl,verbose,silent, &
                         Ndecimal)
    
    USE utilities, ONLY: DP, trimlr, adjustc, NaN, today, programrunner, &
                         verbatim, strike
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: file
    REAL(DP), INTENT(IN), DIMENSION(:) :: vec1
    REAL(DP), INTENT(IN), OPTIONAL, DIMENSION(:) :: vec2, vec3, vec4, vec5
    REAL(DP), INTENT(IN), OPTIONAL, DIMENSION(:,:) :: mat1, mat2, mat3, mat4
    CHARACTER(*), INTENT(IN), OPTIONAL :: comgen, comvec1, comvec2, comvec3
    CHARACTER(*), INTENT(IN), OPTIONAL :: comvec4, comvec5, commat1, commat2
    CHARACTER(*), INTENT(IN), OPTIONAL :: commat3, commat4
    LOGICAL, INTENT(IN), OPTIONAL :: idl, verbose, silent
    INTEGER, INTENT(IN), OPTIONAL :: Ndecimal
 
    INTEGER :: i, j, Nline, Nvec, Nblock, Nmat, lenreal, Ndecim
    REAL(DP), DIMENSION(SIZE(vec1)) :: v1, v2, v3, v4, v5
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: m1, m2, m3, m4
    CHARACTER(12) :: formDP
    CHARACTER(lencom) :: comment
    CHARACTER(textwid) :: comg, comv1, comv2, comv3, comv4, comv5
    CHARACTER(textwid) :: comm1, comm2, comm3, comm4
    LOGICAL :: v2_flag, v3_flag, v4_flag, v5_flag
    LOGICAL :: m1_flag, m2_flag, m3_flag, m4_flag, idl_flag, doprint

    !-----------------------------------------------------------------------

    ! Parameter checking
    !-------------------
    IF (PRESENT(vec5) .AND. .NOT. PRESENT(vec4)) &
      CALL STRIKE("WRITE_ASCII","if VEC5 is set, then VEC4 should be set, too")
    IF (PRESENT(vec4) .AND. .NOT. PRESENT(vec3)) &
      CALL STRIKE("WRITE_ASCII","if VEC4 is set, then VEC3 should be set, too")
    IF (PRESENT(vec3) .AND. .NOT. PRESENT(vec2)) &
      CALL STRIKE("WRITE_ASCII","if VEC3 is set, then VEC2 should be set, too")
    IF (PRESENT(mat4) .AND. .NOT. PRESENT(mat3)) &
      CALL STRIKE("WRITE_ASCII","if MAT4 is set, then MAT3 should be set, too")
    IF (PRESENT(mat3) .AND. .NOT. PRESENT(mat2)) &
      CALL STRIKE("WRITE_ASCII","if MAT3 is set, then MAT2 should be set, too")
    IF (PRESENT(mat2) .AND. .NOT. PRESENT(mat1)) &
      CALL STRIKE("WRITE_ASCII","if MAT2 is set, then MAT1 should be set, too")


    ! Header parameters
    !------------------
    ! Vector sizes
    Nblock = SIZE(vec1(:))
    v1(:) = vec1(:)
    v2pres: IF (PRESENT(vec2)) THEN 
      v2(:) = vec2(:)
      v2_flag = .True.
    ELSE 
      v2_flag = .False.
    END IF v2pres
    v3pres: IF (PRESENT(vec3)) THEN 
      v3(:) = vec3(:)
      v3_flag = .True.
    ELSE 
      v3_flag = .False.
    END IF v3pres
    v4pres: IF (PRESENT(vec4)) THEN 
      v4(:) = vec4(:)
      v4_flag = .True.
    ELSE 
      v4_flag = .False.
    END IF v4pres
    v5pres: IF (PRESENT(vec5)) THEN 
      v5(:) = vec5(:)
      v5_flag = .True.
    ELSE 
      v5_flag = .False.
    END IF v5pres

    ! Matrix sizes
    Nvec = 1 + MERGE(1,0,v2_flag) + MERGE(1,0,v3_flag) &
           + MERGE(1,0,v4_flag) + MERGE(1,0,v5_flag)
    m1pres: IF (PRESENT(mat1)) THEN 
      Nline = SIZE(mat1(:,:),DIM=2)
      ALLOCATE(m1(Nblock,Nline))
      m1(:,:) = mat1(:,:)
      m1_flag = .True.
      m2pres: IF (PRESENT(mat2)) THEN 
        ALLOCATE(m2(Nblock,Nline))
        m2(:,:) = mat2(:,:)
        m2_flag = .True.
      ELSE 
        ALLOCATE(m2(1,1))
        m2(:,:) = 0._DP
        m2_flag = .False.
      END IF m2pres
      m3pres: IF (PRESENT(mat3)) THEN 
        ALLOCATE(m3(Nblock,Nline))
        m3(:,:) = mat3(:,:)
        m3_flag = .True.
      ELSE 
        ALLOCATE(m3(1,1))
        m3(:,:) = 0._DP
        m3_flag = .False.
      END IF m3pres
      m4pres: IF (PRESENT(mat4)) THEN 
        ALLOCATE(m4(Nblock,Nline))
        m4(:,:) = mat4(:,:)
        m4_flag = .True.
      ELSE 
        ALLOCATE(m4(1,1))
        m4(:,:) = 0._DP
        m4_flag = .False.
      END IF m4pres
    ELSE
      Nline = 0
      m1_flag = .False.
      m2_flag = .False.
      m3_flag = .False.
      m4_flag = .False.
      ALLOCATE(m1(1,1),m2(1,1),m3(1,1),m4(1,1))
      m1(:,:) = 0._DP
      m2(:,:) = 0._DP
      m3(:,:) = 0._DP
      m4(:,:) = 0._DP
    END IF m1pres
    Nmat = MERGE(1,0,m1_flag) + MERGE(1,0,m2_flag) + MERGE(1,0,m3_flag) &
         + MERGE(1,0,m4_flag)

    ! Optional keywords
    IF (PRESENT(comgen)) THEN ; comg = comgen ; ELSE ; comg = "" ; END IF
    IF (PRESENT(comvec1)) THEN ; comv1 = comvec1 ; ELSE ; comv1 = "" ; END IF
    IF (PRESENT(comvec2)) THEN ; comv2 = comvec2 ; ELSE ; comv2 = "" ; END IF
    IF (PRESENT(comvec3)) THEN ; comv3 = comvec3 ; ELSE ; comv3 = "" ; END IF
    IF (PRESENT(comvec4)) THEN ; comv4 = comvec4 ; ELSE ; comv4 = "" ; END IF
    IF (PRESENT(comvec5)) THEN ; comv5 = comvec5 ; ELSE ; comv5 = "" ; END IF
    IF (PRESENT(commat1)) THEN ; comm1 = commat1 ; ELSE ; comm1 = "" ; END IF
    IF (PRESENT(commat2)) THEN ; comm2 = commat2 ; ELSE ; comm2 = "" ; END IF
    IF (PRESENT(commat3)) THEN ; comm3 = commat3 ; ELSE ; comm3 = "" ; END IF
    IF (PRESENT(commat4)) THEN ; comm4 = commat4 ; ELSE ; comm4 = "" ; END IF


    ! Precision
    !----------
    ! Format
    Ndecim = Ndecim_def
    IF (PRESENT(Ndecimal)) Ndecim = Ndecimal
    lenreal = Ndecim + 8
    WRITE(formDP,"('ES',I2,'.',I1)") lenreal, Ndecim

    ! (IDL does not read exponents larger than 99)
    IF (PRESENT(idl)) THEN ; idl_flag = idl 
      ELSE ; idl_flag = .False. ; END IF
    IF (idl_flag) THEN
      WHERE (ABS(v1(:)) >= 1.E100_DP) 
        v1(:) = NaN()
      ELSEWHERE (ABS(v1(:)) < 1.E-99_DP) 
        v1(:) = 0._DP
      END WHERE
      IF (v2_flag) THEN
        WHERE (ABS(v2(:)) >= 1.E100_DP) 
          v2(:) = NaN()
        ELSEWHERE (ABS(v2(:)) < 1.E-99_DP) 
          v2(:) = 0._DP
        END WHERE
      END IF
      IF (v3_flag) THEN
        WHERE (ABS(v3(:)) >= 1.E100_DP) 
          v3(:) = NaN()
        ELSEWHERE (ABS(v3(:)) < 1.E-99_DP) 
          v3(:) = 0._DP
        END WHERE
      END IF
      IF (v4_flag) THEN
        WHERE (ABS(v4(:)) >= 1.E100_DP) 
          v4(:) = NaN()
        ELSEWHERE (ABS(v4(:)) < 1.E-99_DP) 
          v4(:) = 0._DP
        END WHERE
      END IF
      IF (v5_flag) THEN
        WHERE (ABS(v5(:)) >= 1.E100_DP) 
          v5(:) = NaN()
        ELSEWHERE (ABS(v5(:)) < 1.E-99_DP) 
          v5(:) = 0._DP
        END WHERE
      END IF
      IF (m1_flag) THEN
        WHERE (ABS(m1(:,:)) >= 1.E100_DP) 
          m1(:,:) = NaN()
        ELSEWHERE (ABS(m1(:,:)) < 1.E-99_DP) 
          m1(:,:) = 0._DP
        END WHERE
      END IF
      IF (m2_flag) THEN
        WHERE (ABS(m2(:,:)) >= 1.E100_DP) 
          m2(:,:) = NaN()
        ELSEWHERE (ABS(m2(:,:)) < 1.E-99_DP) 
          m2(:,:) = 0._DP
        END WHERE
      END IF
      IF (m3_flag) THEN
        WHERE (ABS(m3(:,:)) >= 1.E100_DP) 
          m3(:,:) = NaN()
        ELSEWHERE (ABS(m3(:,:)) < 1.E-99_DP) 
          m3(:,:) = 0._DP
        END WHERE
      END IF
      IF (m4_flag) THEN
        WHERE (ABS(m4(:,:)) >= 1.E100_DP) 
          m4(:,:) = NaN()
        ELSEWHERE (ABS(m4(:,:)) < 1.E-99_DP) 
          m4(:,:) = 0._DP
        END WHERE
      END IF
    END IF


    ! File writing
    !-------------
    OPEN(unitdata,FILE=TRIMLR(file),STATUS="REPLACE",ACTION="WRITE")

      ! File Header
      WRITE(unitdata,"("//formline//")") REPEAT('*',textwid)
      WRITE(unitdata,"(A75)") &
        "CUSTOMIZED FORMAT DATA FILE USED BY SUBROUTINES WRITE_ASCII/READ_ASCII"
      WRITE(unitdata,"("//formline//")") REPEAT('*',textwid)
      WRITE(unitdata,*)
      WRITE(unitdata,"(A8)") "FORMAT"
      WRITE(unitdata,*)
      comment = " ! Number of vectors"
      WRITE(unitdata,"("//formparn//","//formint//","//formparc//")") &
        "Nvec = ", Nvec, comment
      comment = " ! Length of vectors"
      WRITE(unitdata,"("//formparn//","//formint//","//formparc//")") &
        "Nblock = ", Nblock, comment
      comment = " ! Number of matrices"
      WRITE(unitdata,"("//formparn//","//formint//","//formparc//")") &
        "Nmat = ", Nmat, comment
      comment = " ! Width of matrices"
      WRITE(unitdata,"("//formparn//","//formint//","//formparc//")") &
        "Nline = ", Nline, comment
      comment = " ! Length of real numbers"
      WRITE(unitdata,"("//formparn//","//formint//","//formparc//")") &
        "lenreal = ", lenreal, comment
      comment = " ! Number of decimal digits"
      WRITE(unitdata,"("//formparn//","//formint//","//formparc//")") &
        "Ndecim = ", Ndecim, comment
      WRITE(unitdata,*)
      WRITE(unitdata,"("//formline//")") REPEAT('-',textwid)
      WRITE(unitdata,*)
      WRITE(unitdata,"(A10)") "COMMENTS"
      WRITE(unitdata,*)
      WRITE(unitdata,"(A12,"//formline//")") "File: ", comg
      WRITE(unitdata,"(A12,"//formline//")") "Vector 1: ", comv1
      IF (v2_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Vector 2: ", comv2
      IF (v3_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Vector 3: ", comv3
      IF (v4_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Vector 4: ", comv4
      IF (v5_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Vector 5: ", comv5
      IF (m1_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Matrix 1: ", comm1
      IF (m2_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Matrix 2: ", comm2
      IF (m3_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Matrix 3: ", comm3
      IF (m4_flag) &
        WRITE(unitdata,"(A12,"//formline//")") "Matrix 4: ", comm4
      WRITE(unitdata,*)
      WRITE(unitdata,"("//formline//")") REPEAT('*',textwid)
      WRITE(unitdata,*)
      WRITE(unitdata,"(A6)") "DATA"
      WRITE(unitdata,*)
      WRITE(unitdata,"("//formline//")") REPEAT('-',textwid)

      ! Data blocks
      block: DO i=1,Nblock

        ifblock: IF (.NOT. v2_flag) THEN 
          WRITE(unitdata,"(1("//formDP//"))") v1(i)
        ELSE IF (v2_flag .AND. .NOT. v3_flag) THEN
          WRITE(unitdata,"(2("//formDP//"))") v1(i), v2(i)
        ELSE IF (v3_flag .AND. .NOT. v4_flag) THEN
          WRITE(unitdata,"(3("//formDP//"))") v1(i), v2(i), v3(i)
        ELSE IF (v4_flag .AND. .NOT. v5_flag) THEN
          WRITE(unitdata,"(4("//formDP//"))") v1(i), v2(i), v3(i), v4(i)
        ELSE
          WRITE(unitdata,"(5("//formDP//"))") v1(i), v2(i), v3(i), v4(i), v5(i)
        END IF ifblock
  
        line: DO j=1,Nline 
          ifline: IF (m1_flag .AND. .NOT. m2_flag) THEN
            WRITE(unitdata,"(1("//formDP//"))") m1(i,j)
          ELSE IF (m2_flag .AND. .NOT. m3_flag) THEN
            WRITE(unitdata,"(2("//formDP//"))") m1(i,j), m2(i,j)
          ELSE IF (m3_flag .AND. .NOT. m4_flag) THEN
            WRITE(unitdata,"(3("//formDP//"))") m1(i,j), m2(i,j), m3(i,j)
          ELSE
            WRITE(unitdata,"(4("//formDP//"))") m1(i,j), m2(i,j), m3(i,j), &
                                                m4(i,j)
          END IF ifline
        END DO line
        IF (Nmat>0 .AND. i<Nblock) WRITE(unitdata,"("//formline//")") &
          REPEAT('-',textwid)

      END DO block
      
      WRITE(unitdata,"("//formline//")") REPEAT('-',textwid)
      WRITE(unitdata,*)
      WRITE(unitdata,*) "Generated by "//programrunner//", on " &
                        //TRIMLR(TODAY())//"."

    CLOSE(unitdata)

    ! Courtesy notification
    doprint = verbatim
    IF (PRESENT(verbose)) doprint = verbose
    IF (PRESENT(silent)) doprint = ( .NOT. silent )
    IF (doprint) PRINT*, " - File "//TRIMLR(file)//" has been written."

    !-----------------------------------------------------------------------

  END SUBROUTINE write_ascii


  !==========================================================================
  ! CALL READ_ASCII (file,vec1,vec2,vec3,vec4,vec5,mat1,mat2,mat3,mat4)
  !
  !   Reads my customised files with a header and data blocks.
  !==========================================================================

  SUBROUTINE read_ascii (file,vec1,vec2,vec3,vec4,vec5,mat1,mat2,mat3,mat4)
    
    USE utilities, ONLY: DP
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: file
    REAL(DP), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: vec1
    REAL(DP), INTENT(OUT), ALLOCATABLE, OPTIONAL, DIMENSION(:) :: vec2, vec3
    REAL(DP), INTENT(OUT), ALLOCATABLE, OPTIONAL, DIMENSION(:) :: vec4, vec5
    REAL(DP), INTENT(OUT), ALLOCATABLE, OPTIONAL, DIMENSION(:,:) :: mat1, mat2
    REAL(DP), INTENT(OUT), ALLOCATABLE, OPTIONAL, DIMENSION(:,:) :: mat3, mat4

    INTEGER :: i, j, Nline, Nvec, Nblock, Nmat, lenreal, Ndecim
    CHARACTER(12) :: formDP
    CHARACTER(lenpar) :: parn

    !-----------------------------------------------------------------------

    ! File reading
    !-------------
    OPEN(unitdata,FILE=file,STATUS="OLD",ACTION="READ")

      ! File Header
      DO i=1,6 ; READ(unitdata,*) ; END DO
      READ(unitdata,"("//formparn//","//formint//")") parn, Nvec
      READ(unitdata,"("//formparn//","//formint//")") parn, Nblock
      READ(unitdata,"("//formparn//","//formint//")") parn, Nmat
      READ(unitdata,"("//formparn//","//formint//")") parn, Nline
      READ(unitdata,"("//formparn//","//formint//")") parn, lenreal
      READ(unitdata,"("//formparn//","//formint//")") parn, Ndecim
      DO i=1,6+Nvec+Nmat+6 ; READ(unitdata,*) ; END DO

      ! Format
      WRITE(formDP,"('ES',I2,'.',I1)") lenreal, Ndecim

      ! Data blocks
      block: DO i=1,Nblock

        ifblock: IF (.NOT. PRESENT(vec2) .AND. Nvec<2) THEN 
          IF (i==1) ALLOCATE(vec1(Nblock))
          READ(unitdata,"(1("//formDP//"))") vec1(i)
        ELSE IF (PRESENT(vec2) .AND. .NOT. PRESENT(vec3) .AND. Nvec<3) THEN
          IF (i==1) ALLOCATE(vec1(Nblock),vec2(Nblock))
          READ(unitdata,"(2("//formDP//"))") vec1(i), vec2(i)
        ELSE IF (PRESENT(vec3) .AND. .NOT. PRESENT(vec4) .AND. Nvec<4) THEN
          IF (i==1) ALLOCATE(vec1(Nblock),vec2(Nblock),vec3(Nblock))
          READ(unitdata,"(3("//formDP//"))") vec1(i), vec2(i), vec3(i)
        ELSE IF (PRESENT(vec4) .AND. .NOT. PRESENT(vec5) .AND. Nvec<5) THEN
          IF (i==1) ALLOCATE(vec1(Nblock),vec2(Nblock),vec3(Nblock), &
                             vec4(Nblock))
          READ(unitdata,"(4("//formDP//"))") vec1(i), vec2(i), vec3(i), vec4(i)
        ELSE
          IF (i==1) ALLOCATE(vec1(Nblock),vec2(Nblock),vec3(Nblock), &
                             vec4(Nblock),vec5(Nblock))
          READ(unitdata,"(5("//formDP//"))") vec1(i), vec2(i), vec3(i), &
            vec4(i), vec5(i)
        END IF ifblock

        line: DO j=1,Nline 
          ifline: IF (PRESENT(mat1) .AND. .NOT. PRESENT(mat2)) THEN
            IF (i==1 .AND. j==1) ALLOCATE(mat1(Nblock,Nline))
            READ(unitdata,"(1("//formDP//"))") mat1(i,j)
          ELSE IF (PRESENT(mat2) .AND. .NOT. PRESENT(mat3)) THEN
            IF (i==1 .AND. j==1) ALLOCATE(mat1(Nblock,Nline),mat2(Nblock,Nline))
            READ(unitdata,"(2("//formDP//"))") mat1(i,j), mat2(i,j)
          ELSE IF (PRESENT(mat3) .AND. .NOT. PRESENT(mat4)) THEN
            IF (i==1 .AND. j==1) &
              ALLOCATE(mat1(Nblock,Nline),mat2(Nblock,Nline),mat3(Nblock,Nline))
            READ(unitdata,"(3("//formDP//"))") mat1(i,j), mat2(i,j), mat3(i,j)
          ELSE
            IF (i==1 .AND. j==1) &
              ALLOCATE(mat1(Nblock,Nline),mat2(Nblock,Nline), &
                       mat3(Nblock,Nline),mat4(Nblock,Nline))
            READ(unitdata,"(4("//formDP//"))") mat1(i,j), mat2(i,j), &
                                               mat3(i,j), mat4(i,j)
          END IF ifline
        END DO line
        IF (Nmat>0) READ(unitdata,*)

      END DO block

    CLOSE(unitdata)

    !-----------------------------------------------------------------------

  END SUBROUTINE read_ascii


  !==========================================================================
  ! CALL WRITE_BINARY (array,FILE=file)
  !
  !   Writes binary file for any array.
  !==========================================================================

  SUBROUTINE write_binary_1D (array,file,verbose,silent)

    USE utilities, ONLY: DP, trimlr, verbatim
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: file
    LOGICAL, INTENT(IN), OPTIONAL :: verbose, silent

    INTEGER :: N
    LOGICAL :: doprint

    !-----------------------------------------------------------------------

    N = SIZE(array)  
    OPEN (unitdata,FILE=file,STATUS="REPLACE",ACTION="WRITE",FORM="UNFORMATTED")
      WRITE (unitdata) N
      WRITE (unitdata) array(:)
    CLOSE (unitdata)

    ! Courtesy notification
    doprint = verbatim
    IF (PRESENT(verbose)) doprint = verbose
    IF (PRESENT(silent)) doprint = ( .NOT. silent )
    IF (doprint) PRINT*, " - File "//TRIMLR(file)//" has been written."

    !-----------------------------------------------------------------------

  END SUBROUTINE write_binary_1D

  !==========================================================================

  SUBROUTINE write_binary_2D (array,file,verbose,silent)

    USE utilities, ONLY: DP, trimlr, verbatim
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: file
    LOGICAL, INTENT(IN), OPTIONAL :: verbose, silent

    INTEGER :: N1, N2
    LOGICAL :: doprint

    !-----------------------------------------------------------------------

    N1 = SIZE(array,1)  
    N2 = SIZE(array,2)  
    OPEN (unitdata,FILE=file,STATUS="REPLACE",ACTION="WRITE",FORM="UNFORMATTED")
      WRITE (unitdata) N1, N2
      WRITE (unitdata) array(:,:)
    CLOSE (unitdata)

    ! Courtesy notification
    doprint = verbatim
    IF (PRESENT(verbose)) doprint = verbose
    IF (PRESENT(silent)) doprint = ( .NOT. silent )
    IF (doprint) PRINT*, " - File "//TRIMLR(file)//" has been written."

    !-----------------------------------------------------------------------

  END SUBROUTINE write_binary_2D

  !==========================================================================

  SUBROUTINE write_binary_3D (array,file,verbose,silent)

    USE utilities, ONLY: DP, trimlr, verbatim
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: file
    LOGICAL, INTENT(IN), OPTIONAL :: verbose, silent

    INTEGER :: N1, N2, N3
    LOGICAL :: doprint

    !-----------------------------------------------------------------------

    N1 = SIZE(array,1)  
    N2 = SIZE(array,2)  
    N3 = SIZE(array,3)  
    OPEN (unitdata,FILE=file,STATUS="REPLACE",ACTION="WRITE",FORM="UNFORMATTED")
      WRITE (unitdata) N1, N2, N3
      WRITE (unitdata) array(:,:,:)
    CLOSE (unitdata)

    ! Courtesy notification
    doprint = verbatim
    IF (PRESENT(verbose)) doprint = verbose
    IF (PRESENT(silent)) doprint = ( .NOT. silent )
    IF (doprint) PRINT*, " - File "//TRIMLR(file)//" has been written."

    !-----------------------------------------------------------------------

  END SUBROUTINE write_binary_3D

  !==========================================================================

  SUBROUTINE write_binary_4D (array,file,verbose,silent)

    USE utilities, ONLY: DP, trimlr, verbatim
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: array
    CHARACTER(*), INTENT(IN) :: file
    LOGICAL, INTENT(IN), OPTIONAL :: verbose, silent

    INTEGER :: N1, N2, N3, N4
    LOGICAL :: doprint

    !-----------------------------------------------------------------------

    N1 = SIZE(array,1)  
    N2 = SIZE(array,2)  
    N3 = SIZE(array,3)  
    N4 = SIZE(array,4)  
    OPEN (unitdata,FILE=file,STATUS="REPLACE",ACTION="WRITE",FORM="UNFORMATTED")
      WRITE (unitdata) N1, N2, N3, N4
      WRITE(unitdata) array(:,:,:,:)
    CLOSE (unitdata)

    ! Courtesy notification
    doprint = verbatim
    IF (PRESENT(verbose)) doprint = verbose
    IF (PRESENT(silent)) doprint = ( .NOT. silent )
    IF (doprint) PRINT*, " - File "//TRIMLR(file)//" has been written."

    !-----------------------------------------------------------------------

  END SUBROUTINE write_binary_4D


  !==========================================================================
  ! CALL READ_BINARY (array,FILE=file)
  !
  !   Reads binary file for any array.
  !==========================================================================

  SUBROUTINE read_binary_1D (array,file)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: array
    CHARACTER(*), INTENT(IN) :: file

    INTEGER :: N

    !-----------------------------------------------------------------------

    OPEN (unitdata,FILE=file,STATUS="OLD",ACTION="READ",FORM="UNFORMATTED")
      READ (unitdata) N
      ALLOCATE (array(N))
      READ (unitdata) array
    CLOSE (unitdata)

    !-----------------------------------------------------------------------

  END SUBROUTINE read_binary_1D

  !==========================================================================

  SUBROUTINE read_binary_2D (array,file)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: array
    CHARACTER(*), INTENT(IN) :: file

    INTEGER :: N1, N2

    !-----------------------------------------------------------------------

    OPEN (unitdata,FILE=file,STATUS="OLD",ACTION="READ",FORM="UNFORMATTED")
      READ (unitdata) N1, N2
      ALLOCATE (array(N1,N2))
      READ (unitdata) array
    CLOSE (unitdata)

    !-----------------------------------------------------------------------

  END SUBROUTINE read_binary_2D

  !==========================================================================

  SUBROUTINE read_binary_3D (array,file)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: array
    CHARACTER(*), INTENT(IN) :: file

    INTEGER :: N1, N2, N3

    !-----------------------------------------------------------------------

    OPEN (unitdata,FILE=file,STATUS="OLD",ACTION="READ",FORM="UNFORMATTED")
      READ (unitdata) N1, N2, N3
      ALLOCATE (array(N1,N2,N3))
      READ (unitdata) array
    CLOSE (unitdata)

    !-----------------------------------------------------------------------

  END SUBROUTINE read_binary_3D

  !==========================================================================

  SUBROUTINE read_binary_4D (array,file)

    USE utilities, ONLY: DP
    IMPLICIT NONE
    
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(OUT) :: array
    CHARACTER(*), INTENT(IN) :: file

    INTEGER :: N1, N2, N3, N4

    !-----------------------------------------------------------------------

    OPEN (unitdata,FILE=file,STATUS="OLD",ACTION="READ",FORM="UNFORMATTED")
      READ (unitdata) N1, N2, N3, N4
      ALLOCATE (array(N1,N2,N3,N4))
      READ(unitdata) array(:,:,:,:)
    CLOSE (unitdata)

    !-----------------------------------------------------------------------

  END SUBROUTINE read_binary_4D


  !==========================================================================
  !  Subroutine to read individual lines of my customized format input files
  !==========================================================================

  SUBROUTINE read_input_line (unit,parname,parval,iostat)

    USE utilities, ONLY: DP, trimLR, pring
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(INOUT) :: iostat
    CHARACTER(lenline), INTENT(OUT) :: parname, parval

    INTEGER :: ieq, icom
    CHARACTER(lenline) :: line
    
    !------------------------------------------------------------------------

    READ(unit,"(A"//TRIMLR(PRING(lenline))//")",IOSTAT=iostat) line
    ieq = INDEX(line,"=")
    icom = INDEX(line,"!")
    IF (icom == 0) icom = lenline
    IF (ieq == 0) THEN
      parname = ""
      parval = ""
    ELSE
      parname = TRIMLR(line(1:ieq-1))
      parval = TRIMLR(line(ieq+1:icom-1))
    END IF
 
    !------------------------------------------------------------------------

  END SUBROUTINE read_input_line


  !==========================================================================
  !  Subroutine to read individual a line in my customized format input files
  !==========================================================================

  SUBROUTINE edit_input_line (parname,inputfile,outputfile, &
                              parval_DP,parval_Int,parval_Bool,parval_Str, &
                              verbose)

    USE utilities, ONLY: DP, trimlr, trimeq, pring, verbatim
    USE arrays, ONLY: incrarr
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: parname, inputfile
    CHARACTER(*), INTENT(IN), OPTIONAL :: outputfile
    REAL(DP), INTENT(IN), OPTIONAL :: parval_DP
    INTEGER, INTENT(IN), OPTIONAL :: parval_Int
    LOGICAL, INTENT(IN), OPTIONAL :: parval_Bool
    CHARACTER(*), INTENT(IN), OPTIONAL :: parval_Str
    LOGICAL, INTENT(IN), OPTIONAL :: verbose

    INTEGER, PARAMETER :: unit = 1

    INTEGER :: i, Nline, iostat, ieq, icom, lenl
    LOGICAL :: verb
    CHARACTER(lenpath) :: outfile
    CHARACTER(lenline) :: line, parname_read, parval_read, parval_write
    CHARACTER(lenline), DIMENSION(:), ALLOCATABLE :: filarr
    
    !------------------------------------------------------------------------
    
    ! 1) Arguments
    !-------------
    IF (PRESENT(outputfile)) THEN ; outfile = outputfile
                             ELSE ; outfile = inputfile ; END IF
    IF (PRESENT(verbose)) THEN ; verb = verbose ; ELSE ; verb = verbatim ; ENDIF

                               
    ! 2) Read the file
    !-----------------
    OPEN (unit,FILE=inputfile,STATUS="OLD",ACTION="READ",IOSTAT=iostat)
      DO 
        IF (iostat /= 0) EXIT
        READ(unit,"(A"//TRIMLR(PRING(lenline))//")",IOSTAT=iostat) line

        ! Read the parameter value
        ieq = INDEX(line,"=")
        icom = INDEX(line,"!")
        IF (icom == 0) icom = lenline
        IF (ieq == 0) THEN
          parname_read = ""
          parval_read = ""
        ELSE
          parname_read = TRIMLR(line(1:ieq-1))
          parval_read = TRIMLR(line(ieq+1:icom-1))
        END IF

        ! Update the parameter if relevant
        IF (TRIMEQ(parname,parname_read)) THEN
          IF (PRESENT(parval_DP)) WRITE(parval_write,*) parval_DP 
          IF (PRESENT(parval_Int)) WRITE(parval_write,*) parval_Int 
          IF (PRESENT(parval_Bool)) parval_write = MERGE("T","F",parval_Bool)
          IF (PRESENT(parval_Str)) parval_write = parval_Str
          WRITE(line(ieq+1:icom-1),"(' ',A"//TRIMLR(PRING(icom-ieq-4))//")") &
            TRIMLR(parval_write)
        END IF

        ! Make an array
        CALL INCRARR(filarr,line)
        
      END DO
    CLOSE (unit)
    Nline = SIZE(filarr)-1
    
    
    ! 3) Write the edited file
    !-------------------------
    OPEN (unit,FILE=outfile,STATUS="REPLACE",ACTION="WRITE")
      DO i=1,Nline
        lenl = LEN_TRIM(filarr(i))
        IF (lenl == 0) lenl = 1
        WRITE(unit,"(A"//TRIMLR(PRING(lenl))//")") TRIM(filarr(i))
      END DO
    CLOSE (unit)
    IF (verb) THEN  
      IF (PRESENT(outputfile)) THEN
        PRINT*, " - File "//TRIMLR(outfile)//" has been written."
      ELSE
        PRINT*, " - File "//TRIMLR(outfile)//" has been edited."
      END IF
    END IF
      
    !------------------------------------------------------------------------

  END SUBROUTINE edit_input_line

  
END MODULE inout

!******************************************************************************
!*
!*                     Physical and Geometrical Constants
!*
!******************************************************************************

  !==========================================================================
  ! 1) AUTHOR: F. Galliano
  !
  ! 2) HISTORY: 
  !    - Created 07/2007 
  !    - Add the CGS constants (02/2014)
  ! 
  ! 3) DESCRIPTION: Defines the main constants
  !==========================================================================

MODULE constants

  USE utilities, ONLY: DP
  IMPLICIT NONE
  PRIVATE

  
  ! Geometrical constants
  !----------------------
  REAL(DP), PARAMETER, PUBLIC :: pi = 3.141592653589793238_DP
  REAL(DP), PARAMETER, PUBLIC :: twopi = 6.283185307179586476925286766559005_DP
  REAL(DP), PARAMETER, PUBLIC :: oneoversqrt2pi = 0.39894227485064260_DP


  ! MKS constants
  !--------------
  TYPE const_MKS
    ! Fundamental constants
    REAL(DP) :: hplanck = 6.6255E-34_DP ! Planck constant [J.s]
    REAL(DP) :: clight = 2.997E8_DP     ! Speed of light in vacuum [m.s-1]
    REAL(DP) :: grav = 6.67E-11_DP      ! Gravitation const [N.m-2.kg-1]
    REAL(DP) :: kboltz = 1.381E-23_DP   ! Boltzman constant [J.K-1]
    REAL(DP) :: stefan = 5.67E-8_DP     ! Stefan-Boltzman [W.m-2.K-4]
    ! Atomic constants
    REAL(DP) :: u = 1.6605E-27_DP       ! Atomic mass unit [kg]
    REAL(DP) :: mp = 1.6726E-27_DP      ! Proton mass [kg]
    REAL(DP) :: mn = 1.6749E-27_DP      ! Neutron mass [kg]
    REAL(DP) :: me = 9.1096E-31_DP      ! Electron mass [kg]
    REAL(DP) :: mH = 1.6738E-27_DP      ! Hydrogen atom mass [kg]
    REAL(DP) :: a0 = 5.292E-11_DP       ! Bohr radius [m]
    REAL(DP) :: Ryd = 1.097373E7_DP     ! Rydberg constant [m-1]
    REAL(DP) :: cRyd = 3.289842E15_DP   ! c x Rydberg constant [Hz]
    REAL(DP) :: txsec = 6.65E-29_DP     ! Thomson cross section [m2]
    REAL(DP) :: eV = 1.602E-19_DP       ! Electron volt [J.eV-1]
    REAL(DP) :: eV2wave = 1.239E-6_DP   ! Wavelength <=> 1 eV [m]
    REAL(DP) :: eV2nu = 2.418E14_DP     ! Frequency <=> 1 eV [Hz]
    REAL(DP) :: eV2temp = 11604.8_DP    ! Temperature <=> 1 eV [K]
    ! Astronomical constants
    REAL(DP) :: Lsun = 3.845E26_DP      ! Solar luminosity [J.s-1]
    REAL(DP) :: Msun = 1.989E30_DP      ! Solar mass [kg]
    REAL(DP) :: Rsun = 6.596E8_DP       ! Solar radius [m]
    REAL(DP) :: pc = 3.086E16_DP        ! Parsec [m]
    REAL(DP) :: kpc = 3.086E19_DP       ! Kiloparsec [m]
    REAL(DP) :: Mpc = 3.086E22_DP       ! Megaparsec [m]
    REAL(DP) :: year = 3.1536E7_DP      ! Year [s]
    REAL(DP) :: AU = 1.496E11_DP        ! Sun-earth distance [m]
    REAL(DP) :: TCMB = 2.728_DP         ! CMB temperature [K]
    ! Conversion between units
    REAL(DP) :: micron = 1.E-6_DP       ! Micron [m]
    REAL(DP) :: angstrom = 1.E-10_DP    ! Angstrom [m]
    REAL(DP) :: erg = 1.E-7_DP          ! Erg [J]
    REAL(DP) :: mJy = 1.E-29_DP         ! mJy [W.m-2.Hz-1]
    REAL(DP) :: Jy = 1.E-26_DP          ! Jy [W.m-2.Hz-1]
  END TYPE const_MKS
  PUBLIC :: const_MKS
  TYPE(const_MKS), SAVE, PUBLIC :: MKS


  ! CGS constants
  !--------------
  TYPE const_CGS
    ! Fundamental constants
    REAL(DP) :: hplanck = 6.6255E-27_DP ! Planck constant [erg.s]
    REAL(DP) :: clight = 2.997E10_DP    ! Speed of light in vacuum [cm.s-1]
    REAL(DP) :: grav = 6.67E-8_DP       ! Gravitation const [dynes.cm-2.g-1]
    REAL(DP) :: kboltz = 1.381E-16_DP   ! Boltzman constant [erg.K-1]
    REAL(DP) :: stefan = 5.67E-5_DP     ! Stefan-Boltzman [erg.s-1.cm-2.K-4]
    ! Atomic constants
    REAL(DP) :: u = 1.6605E-24_DP       ! Atomic mass unit [g]
    REAL(DP) :: mp = 1.6726E-24_DP      ! Proton mass [g]
    REAL(DP) :: mn = 1.6749E-24_DP      ! Neutron mass [g]
    REAL(DP) :: me = 9.1096E-28_DP      ! Electron mass [g]
    REAL(DP) :: mH = 1.6738E-24_DP      ! Hydrogen atom mass [g]
    REAL(DP) :: a0 = 5.292E-9_DP        ! Bohr radius [cm]
    REAL(DP) :: Ryd = 1.097373E5_DP     ! Rydberg constant [cm-1]
    REAL(DP) :: cRyd = 3.289842E15_DP   ! c x Rydberg constant [Hz]
    REAL(DP) :: txsec = 6.65E-25_DP     ! Thomson cross section [cm2]
    REAL(DP) :: eV = 1.602E-12_DP       ! Electron volt [erg.eV-1]
    REAL(DP) :: eV2wave = 1.239E-4_DP   ! Wavelength <=> 1 eV [cm]
    REAL(DP) :: eV2nu = 2.418E14_DP     ! Frequency <=> 1 eV [Hz]
    REAL(DP) :: eV2temp = 11604.8_DP    ! Temperature <=> 1 eV [K]
    ! Astronomical constants
    REAL(DP) :: Lsun = 3.845E33_DP      ! Solar luminosity [erg.s-1]
    REAL(DP) :: Msun = 1.989E33_DP      ! Solar mass [g]
    REAL(DP) :: Rsun = 6.596E10_DP      ! Solar radius [cm]
    REAL(DP) :: pc = 3.086E18_DP        ! Parsec [cm]
    REAL(DP) :: kpc = 3.086E21_DP       ! Kiloparsec [cm]
    REAL(DP) :: Mpc = 3.086E24_DP       ! Megaparsec [cm]
    REAL(DP) :: year = 3.1536E7_DP      ! Year [s]
    REAL(DP) :: AU = 1.496E13_DP        ! Sun-earth distance [cm]
    REAL(DP) :: TCMB = 2.728_DP         ! CMB temperature [K]
    ! Conversion between units
    REAL(DP) :: micron = 1.E-4_DP       ! Micron [cm]
    REAL(DP) :: angstrom = 1.E-8_DP     ! Angstrom [cm]
    REAL(DP) :: joule = 1.E7_DP         ! J [erg]
    REAL(DP) :: mJy = 1.E-26_DP         ! mJy [erg.s-1.cm-2.Hz-1]
    REAL(DP) :: Jy = 1.E-23_DP          ! Jy [erg.s-1.cm-2.Hz-1]
  END TYPE const_CGS
  PUBLIC :: const_CGS
  TYPE(const_CGS), SAVE, PUBLIC :: CGS


  ! Atomic Quantities
  !------------------
  TYPE const_atom
    ! H
    CHARACTER(15) :: H_name = "Hydrogen"
    INTEGER :: H_number = 1
    REAL(DP) :: H_weight = 1.0079_DP
    ! He
    CHARACTER(15) :: He_name = "Helium"
    INTEGER :: He_number = 2
    REAL(DP) :: He_weight = 4.002602_DP
    ! Li
    CHARACTER(15) :: Li_name = "Lithium"
    INTEGER :: Li_number = 3
    REAL(DP) :: Li_weight = 6.941_DP
    ! Be
    CHARACTER(15) :: Be_name = "Beryllium"
    INTEGER :: Be_number = 4
    REAL(DP) :: Be_weight = 9.012182_DP
    ! B
    CHARACTER(15) :: B_name = "Boron"
    INTEGER :: B_number = 5
    REAL(DP) :: B_weight = 10.811_DP
    ! C
    CHARACTER(15) :: C_name = "Carbon"
    INTEGER :: C_number = 6
    REAL(DP) :: C_weight = 12.0107_DP
    ! N
    CHARACTER(15) :: N_name = "Nitrogen"
    INTEGER :: N_number = 7
    REAL(DP) :: N_weight = 14.0067_DP
    ! O
    CHARACTER(15) :: O_name = "Oxygen"
    INTEGER :: O_number = 8
    REAL(DP) :: O_weight = 15.9994_DP
    ! F
    CHARACTER(15) :: F_name = "Fluorine"
    INTEGER :: F_number = 9
    REAL(DP) :: F_weight = 18.998403_DP
    ! Ne
    CHARACTER(15) :: Ne_name = "Neon"
    INTEGER :: Ne_number = 10
    REAL(DP) :: Ne_weight = 20.179_DP
    ! Na
    CHARACTER(15) :: Na_name = "Sodium"
    INTEGER :: Na_number = 11
    REAL(DP) :: Na_weight = 22.98977_DP
    ! Mg
    CHARACTER(15) :: Mg_name = "Magnesium"
    INTEGER :: Mg_number = 12
    REAL(DP) :: Mg_weight = 24.305_DP
    ! Al
    CHARACTER(15) :: Al_name = "Aluminium"
    INTEGER :: Al_number = 13
    REAL(DP) :: Al_weight = 26.98154_DP
    ! Si
    CHARACTER(15) :: Si_name = "Silicon"
    INTEGER :: Si_number = 14
    REAL(DP) :: Si_weight = 28.086_DP
    ! P
    CHARACTER(15) :: P_name = "Phosphorous"
    INTEGER :: P_number = 15
    REAL(DP) :: P_weight = 30.97376_DP
    ! S
    CHARACTER(15) :: S_name = "Hydrogen"
    INTEGER :: S_number = 16
    REAL(DP) :: S_weight = 32.06_DP
    ! Cl
    CHARACTER(15) :: Cl_name = "Chlorine"
    INTEGER :: Cl_number = 17
    REAL(DP) :: Cl_weight = 35.453_DP
    ! Ar
    CHARACTER(15) :: Ar_name = "Argon"
    INTEGER :: Ar_number = 18
    REAL(DP) :: Ar_weight = 39.948_DP
    ! K
    CHARACTER(15) :: K_name = "Potassium"
    INTEGER :: K_number = 19
    REAL(DP) :: K_weight = 39.098_DP
    ! Ca
    CHARACTER(15) :: Ca_name = "Calcium"
    INTEGER :: Ca_number = 20
    REAL(DP) :: Ca_weight = 40.08_DP
    ! Sc
    CHARACTER(15) :: Sc_name = "Scandium"
    INTEGER :: Sc_number = 21
    REAL(DP) :: Sc_weight = 44.9559_DP
    ! Ti
    CHARACTER(15) :: Ti_name = "Titanium"
    INTEGER :: Ti_number = 22
    REAL(DP) :: Ti_weight = 47.90_DP
    ! V
    CHARACTER(15) :: V_name = "Vanadium"
    INTEGER :: V_number = 23
    REAL(DP) :: V_weight = 50.9414_DP
    ! Cr
    CHARACTER(15) :: Cr_name = "Chromium"
    INTEGER :: Cr_number = 24
    REAL(DP) :: Cr_weight = 51.996_DP
    ! Mn
    CHARACTER(15) :: Mn_name = "Manganese"
    INTEGER :: Mn_number = 25
    REAL(DP) :: Mn_weight = 54.9380_DP
    ! Fe
    CHARACTER(15) :: Fe_name = "Iron"
    INTEGER :: Fe_number = 26
    REAL(DP) :: Fe_weight = 55.847_DP
    ! Co
    CHARACTER(15) :: Co_name = "Cobalt"
    INTEGER :: Co_number = 27
    REAL(DP) :: Co_weight = 58.9332_DP
    ! Ni
    CHARACTER(15) :: Ni_name = "Nickel"
    INTEGER :: Ni_number = 28
    REAL(DP) :: Ni_weight = 58.70_DP
    ! Cu
    CHARACTER(15) :: Cu_name = "Copper"
    INTEGER :: Cu_number = 29
    REAL(DP) :: Cu_weight = 63.546_DP
    ! Zn
    CHARACTER(15) :: Zn_name = "Zinc"
    INTEGER :: Zn_number = 30
    REAL(DP) :: Zn_weight = 65.38_DP
    ! Ga
    CHARACTER(15) :: Ga_name = "Gallium"
    INTEGER :: Ga_number = 31
    REAL(DP) :: Ga_weight = 69.72_DP
    ! Ge
    CHARACTER(15) :: Ge_name = "Germanium"
    INTEGER :: Ge_number = 32
    REAL(DP) :: Ge_weight = 72.59_DP
    ! As
    CHARACTER(15) :: As_name = "Arsenic"
    INTEGER :: As_number = 33
    REAL(DP) :: As_weight = 74.9216_DP
    ! Se
    CHARACTER(15) :: Se_name = "Selenium"
    INTEGER :: Se_number = 34
    REAL(DP) :: Se_weight = 78.96_DP
    ! Br
    CHARACTER(15) :: Br_name = "Bromine"
    INTEGER :: Br_number = 35
    REAL(DP) :: Br_weight = 79.904_DP
    ! Kr
    CHARACTER(15) :: Kr_name = "Krypton"
    INTEGER :: Kr_number = 36
    REAL(DP) :: Kr_weight = 83.80_DP 
  END TYPE const_atom
  PUBLIC :: const_atom
  TYPE(const_atom), SAVE, PUBLIC :: atom


  ! Elemental abundances
  !---------------------
  ! From Asplund et al. (2009)
  TYPE const_abund
    REAL(DP) :: Xsun = 0.7381_DP ! M(H)/M(all)
    REAL(DP) :: Ysun = 0.2485_DP ! M(H)/M(all)
    REAL(DP) :: Zsun = 0.0134_DP ! M(H)/M(all)
  END TYPE const_abund
  PUBLIC :: const_abund
  TYPE(const_abund), SAVE, PUBLIC :: abund



END MODULE constants

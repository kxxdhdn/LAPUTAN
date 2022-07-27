Astrodust: Dielectric Function and Absorption and Scattering Cross Sections
===========================================================================

The dielectric functions provided here are from the paper "The
Dielectric Function of Astrodust, and Predictions for Polarization in
the 3.4um and 10um Features" (Draine & Hensley 2021a, ApJ, accepted;
arXiv:2009.11314).  Please cite that paper if you use these data.

The optical cross sections calculated for Astrodust spheroids provided
below are from the paper "Using the Starlight Polarization Efficiency
Integral to Constrain Shapes and Porosities of Interstellar Grains"
(Draine & Hensley 2021b, ApJ, submitted; arXiv:2101.07277).  Please
cite that paper if you use the tabulated cross sections.

Dielectric functions
====================

Draine & Hensley (2021a, ApJ, accepted; arXiv:2009.11314) present a
new estimate for the effective dielectric function of interstellar
grain material, from submm to X-ray energies.  This material consists
of a mixture of amorphous silicate, carbonaceous materials, and other
materials (likely metal oxides, but possibly also some metallic iron).
The files below present dielectric functions for Astrodust material,
for energies from E=2.485e-5 eV [lambda=5.0cm] to E=12.398 keV
[lambda=1.0 Angstrom]

The dielectric function is designed so that Astrodust material can
reproduce the observed infrared and sub-mm opacity of interstellar
dustl. The dielectric function therefore depends on the assumed

   P = porosity (vacuum fraction),

   f_Fe = fraction of the Fe in metallic form,

   axial ratio b/a (b/a < 1 for prolate spheroids, b/a > 1 for oblate
   spheroids).

These are plain ascii files.  Each file is organized as follows:

lines 1-2: comments
lines 3-6336:  E(eV)  [Re(m)-1]  Im(m)  [Re(eps)-1]  Im(eps)


Dielectric function for DH21 Astrodust with porosity=poro, a fraction
fFe of the iron in metallic for, and axial ratio b/a = ba

Files (50 files; each file is 355 kB) 

Filename: index_DH21Ad_Pporo_fFe_ba :

where poro =  0.00,  0.10,  0.20,  0.30,  0.40
      fFe  =  0.00,  0.10
      ba   = 0.333, 0.400, 0.500, 0.625, 1.400, 1.600, 2.000, 3.000

Note: for fFe = 0.10, we only provide data for ba = 0.500 and 2.000  

Scattering and Absorption Cross sections for Astrodust spheroids
================================================================

Draine & Hensley (2021b, ApJ, submitted; arXiv:2101.07277) have
calculated absorption and scattering cross sections for spheroids,
with size a_{eff} [=volume-equivalent radius], using small-spheroid
approximations when a_{eff}/lambda is small, and using the spheroid
code developed by Voshchinnikov and Farafanov (1993,
Astrophys.Sp.Sci. 204, 19) for finite a_{eff}/lambda.

The files contain cross sections for 
169 sizes a_{eff} [0.0003162um to 5.012um, uniform in log(a), 40 per decade] 
and 1129 wavelengths
[0.0912 micron to 3.981 cm, uniform in log(lambda), 200 per decade].

The size of a spheroid of volume V is characterized by 
a_eff = (3V/4pi)^{1/3}, the radius of an equal-volume sphere.

The 169 aeff values (microns) used are listed in the ascii file:
DH21_aeff

The 1129 wavelengths (microns) used are listed in the ascii file:
DH21_wave

Optical cross sections are calculated for oblate spheroids for selected
axial ratios b/a .
For each size and wavelength, three orientations (jori=1-3) are considered.
For a=spheroid symmetry axis,
k=direction of propagation, and E = incident electric field:

jori=1: k  ||  a, E perp a
jori=2: k perp a, E  ||  a
jori=3: k perp a, E perp a

These are (gzipped) ascii files.  Each file is organized as follows:

lines 1-12: comments

There are 169 columns in each row (jrad = 0 - 168)

lines 13-1141:   ((Qext(jrad, jwave, jori=1), jrad=0,168), jwave=0-1128)

lines 1142-2270: ((Qext(jrad, jwave, jori=2), jrad=0,168), jwave=0-1128)

lines 2271-3399: ((Qext(jrad, jwave, jori=3), jrad=0,168), jwave=0-1128)

lines 3400-4528: ((Qabs(jrad, jwave, jori=1), jrad=0,168), jwave=0-1128)

lines 4529-5657: ((Qabs(jrad, jwave, jori=2), jrad=0,168), jwave=0-1128)

lines 5658-6786: ((Qabs(jrad, jwave, jori=3), jrad=0,168), jwave=0-1128)

lines 6787-7915: ((Qsca(jrad, jwave, jori=1), jrad=0,168), jwave=0-1128)

lines 7916-9044: ((Qsca(jrad, jwave, jori=2), jrad=0,168), jwave=0-1128)

lines 9045-10173: ((Qsca(jrad, jwave, jori=3), jrad=0,168), jwave=0-1128)

Files: 50 files; each gzipped file is ~4.0 MB; 17.1 MB after unzipping.

Filenames = q_DH21Ad_Pporo_fFe_ba.dat.gz

where poro =  0.00,  0.10,  0.20,  0.30,  0.40
      fFe  =  0.00,  0.10
      ba   = 0.333, 0.400, 0.500, 0.625, 1.400, 1.600, 2.000, 3.000

Note: for fFe = 0.10, we only provide data for ba = 0.500 and 2.000  


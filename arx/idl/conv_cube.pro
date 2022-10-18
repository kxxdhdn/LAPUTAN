;+
; NAME: conv_cube
;
; PURPOSE: convolve a IRS spectral cube to a common resolution
;
; CATEGORY: data treatment
;
; CALLING SEQUENCE:
; conv_cube,
;          cubename,
;          psfrefname
;          [,psfsize=psfsize]
;          [,pixelsize=pixelsize]
;          [,instrument=instrument]
;          [,keeppsf=keeppsf]
;          [,besensible=besensible]
;          [,recreateconvolutionkernels=recreateconvolutionkernels]
;          [,createpsfref=createpsfref]
;          [,outfitsfilename=outfitsfilename]
;          [,kerneldirname=kerneldirname]
;          [,psfdirname=psfdirname]
;          [,nyquistthreshold=nyquistthreshold]
;          [,help=help]
;
; INPUTS:
;   cubename [string]:
;      name of fits file which contains the IRS spectralcube to be treated
;   psfrefname [string]: 
;      name of fits file which contains the psf to convolve to.
;
; OPTIONAL INPUTS:
;   psfsize [int]: 
;      number of pixels NxN of the psf to be created should match
;      refpsf size (default 1024)
;   pixelsize [float]:
;      size of one pixel in arcsec in the psf to be created should
;      match psfref pixelsize (default: 0.4 arcsec)
;   instrument [string]:
;      origin of the cube ('ll' or 'sl' implemented, default 'll')
;   besensible [float]:
;      if set then wavelengths which have a fwhm larger than besensible of
;      the target psf are not convolved since we will introduce more
;      artefacts than we will fix
;   createpsfref [string]:
;      if set create a reference psf from the longest wavelength found
;      in the cube and save it as createpsfref 
;   outfitsfilename [string]:
;      if set create a fitsfile which contains the resulting convolved
;      cube
;   kerneldirname [string]:
;      where the convolution kernels are to be sought/written (default './')
;   psfdirname [string]:
;      where the psf files are to be sought/written (default './')
;   nyquistthreshold [float]:
;      cut-off threshold in number of samplings per FWHM of the target
;      PSF (default 2*sqrt(2) = 2.83)
;   psfmode [string]:
;      'CUBISM': use psf appropiate for the output of CUBISM (default)
;      'STINYTIM': use psf appropiate simulated using stinytim
;
; KEYWORD PARAMETERS:
;   /keeppsf:
;      if set save the generated psffiles for each wavelength
;   /recreateconvolutionkernels:
;      if set recreate the kernels even if they already exist (useful
;      for example is the pixel size has changed)
;   /help: if set show this help
;
; OUTPUTS:
; a structure which contains the convolved cube, header, wavelengths
; and wavelength header
;
; SIDE EFFECTS:
;    fits files with convolution kernels are created in the current
;    directory
;    if keeppsf is set fits files for each stinytim psf are stored
;    aswell in the current directory
;
; RESTRICTIONS:
;   needs a working stiny2 program.
;
; EXAMPLE:
; cc = conv_cube('hii3_hii34_ll_LL1_cube.fits','psfref.fits')
; fo = conv_cube('hii34/hii3_hii34_sl_SL1_cube.fits','psfref.fits',/besensible,outfits='hii34/convolved_sl1.fits',kerneldir='sl1')
;
; DEPENDS ON:
; astrolib
; KG: change_image_scale conv_image ensure_psf_centered
; SH: default n2s
;
; CAVEATS and WORRIES:
;   I am not yet happy with the resulting kernels.
;   - Do we really need 0.4 arcsec pixels for maps which have 5.08
;   arcsec pixels.
;   - Do we really need 1024 pixels = a 400x400 arcsec map which
;     implies 40 fwhms even for the longest wavelength?
;   - When the psf widths are too close it does not make sense to convolve.
;   - How about stopping when there is less than 20% difference in the
;     FWHM (besensible option)
;   - Check the nyquist sampling frequency for a gaussian. If 2.8 correct?
;   - Why was there a 2*pi in the gordon et al paper?
;
;
; MODIFICATION HISTORY:
; (SH Apr 26 2010) initial version
; (SH Apr 26 2010) implemented besensible
; (SH Apr 29 2010) Allow to create a newfits file inside the routine
; (SH Apr 29 2010) Removed dependence on conv_image. I have more faith
;                  in making my own mistakes
;-

;; routine to take an image and convolve it with a psf. 
;; caveat: should only be used within conv_cube because there is
;; little proper checking on keywords etc
function conv_cube_convolve_image,img,himg_in,convolutionkernel,convolutionkernelpixelsize

;; The main problem is taking care of the pixels sizes. The real
;; problem is the regridding. The best would be to have integer
;; regridding. For now since we are inside the conv_cube routine and
;; we know that the kernel is finely sampled, resample the img to the
;; kernel resolution. This is more calcuation intensive but safer. The
;; best is perhaps to do a hastrom on the same header with a pixel
;; size corresponding to the kernel

  ;; determine pixel size of img from himg
  
  ;; check for existing astrometry
  ;; for the moment just do what is needed for the cubism cubes. Being
  ;; general is too much hassle.
  
  himg = himg_in
  sxaddpar,himg,'NAXIS',2
  sxdelpar,himg,'NAXIS3'
  
  sxaddpar,himg,'EQUINOX',2000.
  
  origcdelt1 = sxpar(himg,'cdelt1',comment=commentcdelt1)
  origcdelt2 = sxpar(himg,'cdelt2',comment=commentcdelt2)
  orignaxis1 = sxpar(himg,'naxis1',comment=commentnaxis1)
  orignaxis2 = sxpar(himg,'naxis2',comment=commentnaxis2)
  origcrpix1 = sxpar(himg,'crpix1',comment=commentcrpix1)
  origcrpix2 = sxpar(himg,'crpix2',comment=commentcrpix2)

  ;; extra border to be sure not to lose pixels
  extendborderfactor = 1.1d0
  scalefactor = ceil(3600d0*abs(origcdelt1)/convolutionkernelpixelsize*extendborderfactor)

  htargetimg = himg

  sxaddpar,htargetimg,'CDELT1',origcdelt1/abs(origcdelt1)*convolutionkernelpixelsize/3600d0,commentcdelt1
  sxaddpar,htargetimg,'CDELT2',origcdelt2/abs(origcdelt2)*convolutionkernelpixelsize/3600d0,commentcdelt2
  sxaddpar,htargetimg,'NAXIS1',orignaxis1*scalefactor,commentnaxis1
  sxaddpar,htargetimg,'NAXIS2',orignaxis2*scalefactor,commentnaxis2
  sxaddpar,htargetimg,'CRPIX1',origcrpix1*scalefactor,commentcrpix1
  sxaddpar,htargetimg,'CRPIX2',origcrpix2*scalefactor,commentcrpix2

  ;; here we have to do something about isolated NAN pixels in the
  ;; input image
  idx_NANs = where(finite(img) ne 0,cnt)
  if cnt ne 0 then begin
     fill_missing,img,id_NANs,2
  endif
  
;  ;; resample on the fine grid using hastrom
  hastrom,img,himg,fineimg,hfineimg,htargetimg,missing=!values.d_nan

  ;;(SH May 23 2013)
  ;; resample on the fine grid using hswarp
;  hswarp,img,himg,fineimg,hfineimg,htargetimg,/keepedge,/fixheader
  
  ;; (SH Apr 27 2010)
  ;; have verified that the reprojection is done correctly
  
  ; remove NANs
  NANidx = where(finite(fineimg) ne 1,NANcnt)
  if NANcnt ne 0 then fineimg[NANidx] = 0d0

  convolvedfineimg = convolve_for_conv_cube(fineimg,convolutionkernel)

;  ;; resample to the original grid
  hastrom,convolvedfineimg,hfineimg,convolvedimg,hconvolvedimg,himg,missing=!values.d_nan
  ;;(SH May 23 2013)
;  hswarp,convolvedfineimg,hfineimg,convolvedimg,hconvolvedimg,himg,/keepedge,/fixheader

  return,convolvedimg
end

;; calculate a convolution kernel based on a psf and a target
;; (reference) psf. Following Gordon 2008 to define a cutoff
;; frequency.
function conv_cube_create_convolutionkernel, $
   psfname, $
   psfrefname, $
   fftpsfref=fftpsfref, $ ;; allows to store the previous calculation and not redo the work
   hanning=hanning, $
   psfreffwhm=psfreffwhm, $
   nyquistthreshold=nyquistthreshold, $
   psfthreshold=psfthreshold
  
  default,nyquistthreshold,2d0*sqrt(2d0)

  psf = readfits(psfname,hpsf)
  psf = psf/total(psf)
  fftpsf = fft(psf,-1)
  
  if n_elements(fftpsfref) eq 0 then begin
     psfref = readfits(psfrefname,hpsfref)
     psfref = psfref/total(psfref)
     fftpsfref = fft(psfref,-1)
  endif

  ;; check inputs
  if n_elements(fftpsf) ne n_elements(fftpsfref) then begin
     message,'psf and psfref do not have the same size',/info
     message,'regridding of the psfs is not implemented'
     ;; could also check headers to make sure pixel sizes correspond.
     ;; for now assume that the user is not an imbecile :)
  endif
  
  nx = n_elements(fftpsfref[*,0])
  ny = n_elements(fftpsfref[0,*])
     
  if n_elements(hanning) eq 0 then begin
     if n_elements(psfreffwhm) ne 1 then begin
        if n_elements(psfref) eq 0 then begin
           psfref = readfits(psfrefname)
        endif
        gaussfit = gauss2dfit(psfref,gaussfitparams,/tilt)
        psfreffwhm = (mean(gaussfitparams[2:3]))*2d0*sqrt(2d0*alog(2d0))
     endif
     
     ;; create a map with distances in pixels to the center of the map
     idx = lindgen(nx,ny)
     xxx = (idx mod nx) - (nx-1)/2 ;; the int rounding is correct because of the way fft is defined
     yyy = (idx  /  nx) - (ny-1)/2
     
     dist = sqrt(xxx^2+yyy^2)
     
     ;; shift the thing to agree with the way frequencies are defined
     ;; in fft
     dist = shift(dist,-(nx-1)/2,-(ny-1)/2)

     ;; so this is wrong unless the image is square
     freq = dist/(double(nx)) ;; in pixel^-1

     ;; if the threshold equals 0 we dont want to filter at all 
     if double(nyquistthreshold) eq 0d0 then begin
        hanning = psf*0d0+1
     endif else begin
        hanning = psf*0d0

        ;; the nyquist frequency of the lowest resolution psf 
        ;; we use 2.8 here for 2*sqrt(2). 2 for a gaussian line and
        ;; sqrt(2) because we are dealing with a 2D image not a 1d line.
        nyquist_freq = (psfreffwhm/nyquistthreshold)^(-1d0) ;; in pixel^-1
        
        ;; from sandstrom
        idx = where(freq lt nyquist_freq,cnt)
        if cnt ne 0 then begin
           ;; like in gordon 2008 except that he has 2*pi while I believe
           ;; it should be pi
           hanning[idx] = 0.5*(1+cos(!dpi*freq[idx]/nyquist_freq))
        endif else begin
           message,'no points in the fft domain that fulfill the nyquist criterion'
        endelse
     endelse
  endif
  
  convolutionkernel = fft(hanning*fftpsfref/fftpsf,1)
  convolutionkernel = shift(convolutionkernel,nx/2,ny/2)
  
  if n_elements(psfthreshold) eq 1 then begin
     if finite(psfthreshold) then begin
        if n_elements(psfref) eq 0 then psfref = readfits(psfrefname)
        
        ;; give a weight of 1 to the inner part os the kernel
        centralregion = 1d0*(psfref gt max(psfref)*psfthreshold)
        
        ;; smoothly decrease to values of the outerparts (suppressed by the strength of the psfref wings
        ;; (SH Apr 12 2013) The may not work for real PSFs
        outerregion   = (psfref le max(psfref)*psfthreshold)*psfref
        outerregion = outerregion/max(outerregion)

        convolutionkernel = convolutionkernel*(centralregion+outerregion)
     endif
  endif

  convolutionkernel = convolutionkernel/total(convolutionkernel)

  return,convolutionkernel
end

;; calculate a convolution kernel based on a psf and a target
;; (reference) psf. Following Gordon 2008 to define a cutoff
;; frequency.
function conv_cube_create_gaussian_convolutionkernel, $
   psfrefname, $
   kernelfwhm=kernelfwhm
  
  psfref = readfits(psfrefname,hpsfref)
  
  nx = n_elements(psfref[*,0])
  ny = n_elements(psfref[0,*])
  
  ;; coordinates compared to the center of the image
  foo = lindgen(nx,ny)
  x = double(foo mod nx) - (nx-1.)/2.
  y = double(foo  /  nx) - (ny-1.)/2.
  r = sqrt(x^2+y^2)
  
  sigma = kernelfwhm/(2d0*sqrt(2d0*alog(2d0)))
  
  convolutionkernel = exp( -1d0* r^2/(2d0*sigma^2))
  convolutionkernel = convolutionkernel/total(convolutionkernel)

  return,convolutionkernel
end

;; write a .par file for stinytim
pro conv_cube_create_stinytim_parfile, $
   wave=wave, $
   instrument=instrument, $
   psfsize=psfsize, $
   pixelsize=pixelsize

  wantll = (strlowcase(strmid(instrument,0,2)) eq 'll')

  openw,lun,'stiny2.par',/get_lun
  
  printf,lun,'2.0    = Version of Tiny Tim that generated this file'
  printf,lun,'stiny2_psf = rootname of output image file(s)'

  if wantll then begin
     printf,lun,'10         = instrument code'
     printf,lun,'IRS_long_low        = instrument name'
  endif else begin
     printf,lun,'12         = instrument code'
     printf,lun,'IRS_short_low        = instrument name'
  endelse
  
  printf,lun,n2s(pixelsize)+'  = PSF model X-axis pixel size in arcsec (might be subsampled)'
  printf,lun,n2s(pixelsize)+'  = PSF model X-axis pixel size in arcsec (might be subsampled)'
  printf,lun,'1  = Is PSF subsampled? (0 = no, 1 = yes)'
  printf,lun,n2s(psfsize)+'    = model grid size (n by n pixels)'
  
  if wantll then begin
     if wave gt 25.86 then BEGIN
        printf,lun,'512    = Nyquist PSF grid size (n by n pixels)'
     endif else BEGIN
        printf,lun,'1024    = Nyquist PSF grid size (n by n pixels)'
     endelse
  endif else begin
     if wave gt 7.71 then BEGIN
        printf,lun,'512    = Nyquist PSF grid size (n by n pixels)'
     endif else BEGIN
        printf,lun,'1024    = Nyquist PSF grid size (n by n pixels)'
     endelse
  endelse
  
  printf,lun,'0.000000  = Jitter (arcsec) (ignored if < 0.001)'
  printf,lun,'(none)   = instrument throughput table file used'
  printf,lun,'(none)   = spectrum type or spectrum file used'
  printf,lun,'1        = Skip wavelengths with low weights? (1=yes)'
  printf,lun,'0.010000  = Weight skip low weight limit/max weight'
  printf,lun,'1    = number of wavelengths in model'
  printf,lun,n2s(wave)+'  1     = Wavelength 1 (microns) & weight'
  printf,lun,'1    = Use zonal error maps? (1=yes)'
  printf,lun,'1    = Use new secondary map? (1=yes)'
  printf,lun,'0    = Write pupil map? (1=yes)'
  printf,lun,'0    = Write OPD map? (1=yes)'
  printf,lun,'0    = Write Nyquist PSF? (1=yes)'
  printf,lun,'#---------------------------------------------------------------'

  if wantll then begin
     printf,lun,'# Zernike aberration coefficient file for IRS long-low spectrograph'
  endif else begin
     printf,lun,'# Zernike aberration coefficient file for IRS short-low spectrograph'
  endelse

  printf,lun,'# Code V predicted aberrations at the slit position'
  printf,lun,'# Microns RMS for 33% obscured Zernike polynomials'
  printf,lun,'# April 2000'
  printf,lun,'#---------------------------------------------------------------'
  printf,lun,'  0.       # Z1 = (Not used)'
  printf,lun,'  0.       # Z2 = X tilt'
  printf,lun,'  0.       # Z3 = Y tilt'
  printf,lun,'  0.       # Z4 = Focus for center of detector'

  if wantll then begin
     printf,lun,' -0.227    # Z5 = 0 degree astigmatism'
     printf,lun,' -0.105    # Z6 = 45 degree astigmatism'
     printf,lun,' -0.0066   # Z7 = X coma'
     printf,lun,'  0.1002   # Z8 = Y coma'
  endif else begin
     printf,lun,'  0.178    # Z5 = 0 degree astigmatism'
     printf,lun,' -0.080    # Z6 = 45 degree astigmatism'
     printf,lun,'  0.       # Z7 = X coma'
     printf,lun,'  0.1728   # Z8 = Y coma'
  endelse

  printf,lun,'  0.       # Z9 = X clover'
  printf,lun,'  0.       # Z10 = Y clover'
  printf,lun,'  0.       # Z11 = 3rd order spherical'
  printf,lun,'  0.       # Z12 = 0 degree Spherical astigmatism'
  printf,lun,'  0.       # Z13 = 45 degree Spherical astigmatism'
  printf,lun,'  0.       # Z14 = X Ashtray'
  printf,lun,'  0.       # Z15 = Y Ashtray'
  printf,lun,'  0.       # Z16'
  printf,lun,'  0.       # Z17'
  printf,lun,'  0.       # Z18'
  printf,lun,'  0.       # Z19'
  printf,lun,'  0.       # Z20'
  printf,lun,'  0.       # Z21'
  printf,lun,'  0.       # Z22 = 5th order spherical'
  printf,lun,'#-------------------------------------------------------------'
  printf,lun,'# SIRTF OTA parameter file'
  printf,lun,'# Created by John Krist, January 2000'
  printf,lun,'# These are preliminary settings'
  printf,lun,'#-------------------------------------------------------------'
  printf,lun,' 850.0      = OTA diameter in mm'
  printf,lun,' 0.3765     = Central obscuration radius normalized to pupil radius'
  printf,lun,' 0.0        = OTA spider rotation offset (degrees)'
  printf,lun,'#-------------------------------------------------------------'

  if wantll then begin
     printf,lun,'# IRS long-low instrument parameter file'
  endif else begin
     printf,lun,'# IRS short-low instrument parameter file'
  endelse

  printf,lun,'# Created by John Krist, June 2006'
  printf,lun,'#-------------------------------------------------------------'
  printf,lun,'# The pixel size given below is 1/3 the slit width'

  if wantll then begin
     printf,lun,' 14.0 40.0  = Instrument min, max wavelength (microns)'
     printf,lun,' 3.5        = X Pixel size in arcsec (1/3 slit width)'
     printf,lun,' 3.5        = Y Pixel size in arcsec (1/3 slit width)'
     printf,lun,' 4.41 -14.00 = Detector center offset from telescope axis (arcmin)'
     printf,lun,'  -1.2      = Degrees ccw from OTA z-y system to detector x-y'
  endif else begin
     printf,lun,'  4.0 15.0  = Instrument min, max wavelength (microns)'
     printf,lun,' 1.2        = X Pixel size in arcsec (1/3 slit width)'
     printf,lun,' 1.2        = Y Pixel size in arcsec (1/3 slit width)'
     printf,lun,'12.03 -2.79 = Detector center offset from telescope axis (arcmin)'
     printf,lun,' -84.7      = Degrees ccw from OTA z-y system to detector x-y'
  endelse

  printf,lun,'          0 = Flip -X to X? (1 = yes)'
  printf,lun,'          0 = Flip -Y to Y? (1 = yes)'
  
  close,lun
  free_lun,lun
  
end

pro conv_cube_create_stinytim_psf, $
   outname=outname, $
   wavelength=wavelength, $
   instrument=instrument, $
   psfsize=psfsize, $
   pixelsize=pixelsize

  ;;;; create the stinytim par file
  conv_cube_create_stinytim_parfile, $
     wave=max(wave), $
     instrument=instrument, $
     psfsize=psfsize, $
     pixelsize=pixelsize

  ;;;; run stinytim to create psf.fits
  spawn,'stiny2 stiny2.par'
  ;;;; move the created kernel to the psfrefname
  file_move,'stiny2_psf.fits',outname

end

;; create a psf which is applicable to the results from CUBISM.
;; The values are a very schematic capturing of the measurements from 

;; Pereira-Santaella, Miguel, Almudena Alonso-Herrero, George H.
;; Rieke, Luis Colina, Tanio Díaz-Santos, J.-D. T. Smith, Pablo G.
;; Pérez-González, and Charles W. Engelbracht. “Local Luminous
;; Infrared Galaxies. I. Spatially Resolved Observations with the
;; Spitzer Infrared Spectrograph.” The Astrophysical Journal
;; Supplement Series 188, no. 2 (June 1, 2010): 447.
;; doi:10.1088/0067-0049/188/2/447.

;; They only measure the short modules. We assume that after that the
;; FWHM follows the FWHM  ~ lambda
;; For spitzer we have FWHM = 4 * lambda

;; the schematic presciption below captures the behiviour measured
;; with a single straight line and a second straight line when the PSF
;; of spitzer become larger that the straight line segment

pro conv_cube_create_cubism_psf, $
   outname=outname, $
   wavelength=wavelength, $
   instrument=instrument, $
   psfsize=psfsize, $
   pixelsize=pixelsize
  
  sim_par_wave = [0,13.25,40]
  sim_par_fwhm = [2.8,3.26,10.1]
  
  sim_per_wave = [0,15.5,40]
  sim_per_fwhm = [3.8,3.8,10.1]
  
  psfref = make_array(psfsize,psfsize,value=0d0)
  
  nx = n_elements(psfref[*,0])
  ny = n_elements(psfref[0,*])
  
  ;; coordinates (arcsec) compared to the center of the image
  foo = lindgen(nx,ny)
  x = (double(foo mod nx) - (nx-1.)/2.)*pixelsize
  y = (double(foo  /  nx) - (ny-1.)/2.)*pixelsize

  ;; fwhm (arcsec)
  fwhm_par = interpol(sim_par_fwhm,sim_par_wave,wavelength)
  fwhm_per = interpol(sim_per_fwhm,sim_per_wave,wavelength)

  ;; sigma (arcsec)
  sigma_par = fwhm_par/(2d0*sqrt(2d0*alog(2d0)))
  sigma_per = fwhm_per/(2d0*sqrt(2d0*alog(2d0)))

  ;; the Parallel / Perpendicular directions of the PSF on the
  ;; generated CUBISM images
  ;;
  ;; X-axis is parallel
  ;; Y-axis is perpendicular

  ;; reduced coordinates in sigma_xxx
  
  r_reduced = sqrt((x/sigma_par)^2+(y/sigma_per)^2)

  psf = exp( -1d0/2d0 * r_reduced^2)
  psf = psf/total(psf)

  ;; create a suitable header
  mkhdr,hdr,psf,/image
  sxaddpar,hdr,'WAVELNTH',wavelength,'Wavelength of the psf in microns'
  sxaddpar,hdr,'INSTRUMENT',instrument,'Instrument name'
  sxaddpar,hdr,'PIXSCALE',pixelsize,'Size of pixels in arcsec'
  sxaddpar,hdr,'X_POS',(psfsize+1)/2,'X position of PSF center in detector pixels'
  sxaddpar,hdr,'Y_POS',(psfsize+1)/2,'Y position of PSF center in detector pixels'
  writefits,outname,float(psf),hdr
end

function conv_cube,cubename,psfrefname, $
                   psfsize=psfsize, $
                   pixelsize=pixelsize, $
                   instrument=instrument, $
                   keeppsf=keeppsf, $
                   recreateconvolutionkernel=recreateconvolutionkernels, $
                   createpsfref=createpsfref, $
                   outfitsfilename=outfitsfilename, $
                   kerneldirname=kerneldirname, $
                   psfdirname=psfdirname, $
                   psfmode=psfmode, $
                   nyquistthreshold=nyquistthreshold, $
                   besensible=besensible, $
                   psfthreshold=psfthreshold, $
                   help=help, $
                   _extra=_extra
  
  default,psfsize,512
  default,pixelsize,0.4 ;; arcsec
  default,psfrefname,'psfref.fits'
  default,keeppsf,0
  
  default,recreateconvolutionkernels,0
  default,createpsfref,0
  default,kerneldirname,'./'
  default,psfdirname,'./'

  default,psfmode,"CUBSIM" ;; can be CUBISM of STINYTIM

  default,besensible,0
  default,nyquistthreshold,2d0*sqrt(2d0)
  default,psfthreshold,!values.d_nan

  needhelp=0
  
  ;; check inputs
  if n_elements(cubename) eq 0 then begin
     message,'no cube name given',/info
     needhelp=1
  endif 

  foo = file_search(cubename,count=cnt)

  case cnt of
     0: begin
        message,'cannot find cube file',/info
        needhelp=1
     end
     1: begin
     end
     else: begin
        message,'cube file name is not unique',/info
        message,'did you use a * in the file name',/info
        needhelp=1
     end
  endcase

  case createpsfref of
     0: begin
        foo = file_search(psfrefname,count=cnt)
        case cnt of
           0: begin
              message,'cannot find psfref file',/info
              needhelp=1
           end
           1: begin
              docreatepsfref = 0
              ;; happy
           end
           else: begin
              message,'psfref file name is not unique',/info
              message,'did you use a * in the file name?',/info
              needhelp=1
           end
        endcase
     end
     1: begin
        ;; psfrefname to be created is psfrefname
        docreatepsfref = 1
     end
     else: begin
        psfrefname = (string(createpsfref))[0]
        docreatepsfref = 1
     end
  endcase
  
  case 1 of 
     strcmp(psfmode,"CUBISM",/fold_case): begin
        use_cubism_mode=1
        instrument='cubism'
     end
     strcmp(psfmode,"STINYTIM",/fold_case): begin
        use_cubism_mode=0
     end
     else: begin
        message,/info,'Could not understand the given PSFMode'
        needhelp=1
     end
  endcase
     
  if needhelp then begin
     doc_library,'conv_cube'
     return,-1
  endif

  ;; read the data
  struct = readfitscube(cubename)
  
  cube = struct.cube
  hcube = struct.header
  wave = struct.wave
  hwave = struct.hwave
  
  if n_elements(instrument) eq 0 then begin
     ;; determine the instrument
     fovname = strlowcase(sxpar(hcube,'FOVNAME'))
     case 1 of
        strpos(fovname,'long-lo') ne -1: begin
           instrument = 'll'
        end
        strpos(fovname,'short-lo') ne -1: begin
           instrument = 'sl'
        end
        else: begin
           message,'could not determine the instrument from which this cube came (ll/sl)'
           instrument=''
        end
     endcase
  endif
  
  ;; create the output cube
  convolved_cube = cube*!values.d_nan

  cube_pixelsize = (abs(sxpar(hcube,'CDELT1')))*3600d0 ;; arcsec

  ;; some convenience definitions
  nwave = n_elements(wave)
  swave = n2s(wave,format='(F10.3)')

  ;; do we want to self define the psfref?
  if docreatepsfref then begin
     if use_cubism_mode then begin
        instrument='cubism'
        conv_cube_create_cubism_psf, $
           wave=max(wave), $
           instrument=instrument, $
           psfsize=psfsize, $
           pixelsize=pixelsize, $
           outname=psfrefname
     endif else begin ;; use stinytim defined kernels
        conv_cube_create_stinytim_psf, $
           wave=max(wave), $
           instrument=instrument, $
           psfsize=psfsize, $
           pixelsize=pixelsize, $
           outname=psfrefname
     endelse
  endif

  ;; read the created or existing psfref file
  psfref = mrdfits(psfrefname,0,hpsfref)
  psfref_wavelength = sxpar(hpsfref,'WAVELNTH')
  
  ;; check the psfref size
  psfref_size = sxpar(hpsfref,'NAXIS1')
  if psfref_size ne psfsize then begin
     message,/info,string(format='("The existing PSF size (",I4,") does not match the requested size (",I4,")")',psfref_size,psfsize)
     message,/info,'Using size of the existing PSF'
     psfsize=psfref_size
  endif
  
  if (besensible ne 0) then begin
     if (besensible gt 0) and (besensible lt 1) then begin
        fwhmthreshold=besensible
     endif else begin
        fwhmthreshold=0.7d0 ;; default value from inspection of SL1
        ;; kernels > 13um for target 17 are too ugly
     endelse
     ;; determine the fwhm of the psfref in pixels
     gaussfit = gauss2dfit(psfref,gaussfitparams,/tilt)
     psfreffwhm = (mean(gaussfitparams[2:3]))*2d0*sqrt(2d0*alog(2d0))
  endif
  
  ;; create a header to stock the kernels based on the psfref
  tmp_psfname='psf.fits'
  houtputconvolutionkernel = hpsfref
  sxaddpar,houtputconvolutionkernel,'PIXSCALE',pixelsize,'Size of pixels in arcsec'
  sxaddpar,houtputconvolutionkernel,'REFWAVE',psfref_wavelength,'Wavelength of the target psf in microns'
  sxaddpar,houtputconvolutionkernel,'INSTRUMENT',instrument,'Instrument name'

  ;; a header for the notuseful files
  houtputconvolutionkernelnotuseful = houtputconvolutionkernel
  
  ;; for each wavelength
  for i=0,nwave-1 do begin
     
     ;; reset these toggles
     havekernel = 0
     wanttocreatekernel = 1
     wanttoconvolve = 1
     
  ;;; define the convolution kernel name
     convolutionkernelname = kerneldirname+'/'+'convolution_kernel_'+instrument+'_'+swave[i]+'.fits'
     
  ;;; define the convolution kernel notuseful name
     convolutionkernelnotusefulname = kerneldirname+'/'+'convolution_kernel_'+instrument+'_'+swave[i]+'_notuseful.fits'
     
     ;; there are many reasons for recreating the kernel:
     ;; 1) recreateconvolutionkernels = 1
     ;; 2) kernel fitsfile does not exist
     ;; 3) existing kernel does not match the request because
     ;;    a) target wave is not the same
     ;;    b) target size is not as requested
     ;;    c) target pixelsize is not as requested
     ;;
     ;; only if all these criteria are fullfilled do we reuse the
     ;; existing kernel
     
     if recreateconvolutionkernels eq 1 then begin
        GOTO,FINISHED_CHECKING_EXISTING_KERNEL
     endif
     
     foo = file_search(convolutionkernelname,count=cnt)
     if (cnt eq 0) then begin
        GOTO,FINISHED_CHECKING_EXISTING_KERNEL
     endif
     
     hconvolutionkernel = headfits(convolutionkernelname)
     convolutionkernel_refwave = sxpar(hconvolutionkernel,'REFWAVE')
     convolutionkernel_size = sxpar(hconvolutionkernel,'NAXIS2')
     convolutionkernel_pixscale = sxpar(hconvolutionkernel,'PIXSCALE')

     if float(convolutionkernel_refwave) ne float(psfref_wavelength) then begin
        message,/info,'convolution kernel exists but the reference wavelength does not match'
        message,/info,'recreating kernel for the correct reference wavelength'
        GOTO,FINISHED_CHECKING_EXISTING_KERNEL
     endif
     
     if float(convolutionkernel_size) ne float(psfsize) then begin
        message,/info,'convolution kernel exists but the image size does not match'
        message,/info,'recreating kernel with new size'
        GOTO,FINISHED_CHECKING_EXISTING_KERNEL
     endif
     
     if float(convolutionkernel_pixscale) ne float(pixelsize) then begin
        message,/info,'convolution kernel exists but the pixel scale does not match'
        message,/info,'recreating kernel with new pixel scale'
        GOTO,FINISHED_CHECKING_EXISTING_KERNEL
     endif
     
     ;; if we get here the existing kernel matches all criteria
     wanttocreatekernel = 0
     havekernel = 1
     
     FINISHED_CHECKING_EXISTING_KERNEL:

     ;; create wavelength psf
     if wanttocreatekernel then begin
        if use_cubism_mode then begin
           instrument='cubism'
           conv_cube_create_cubism_psf, $
              wave=wave[i], $
              instrument=instrument, $
              psfsize=psfsize, $
              pixelsize=pixelsize, $
              outname=tmp_psfname
        endif else begin ;; use stinytim defined kernels
           conv_cube_create_stinytim_psf, $
              wave=wave[i], $
              instrument=instrument, $
              psfsize=psfsize, $
              pixelsize=pixelsize, $
              outname=tmp_psfname
        endelse

        if (besensible ne 0) then begin
           ;; determine the fwhm of the psf in pixels
           psf = readfits(tmp_psfname,h)
           gaussfit = gauss2dfit(psf,gaussfitparams,/tilt)
           psffwhm = (mean(gaussfitparams[2:3]))*2d0*sqrt(2d0*alog(2d0))

           ;; the size of the kernel is over the threshold that we
           ;; have set make either a gaussian kernel of do nothing
           if (psffwhm ge fwhmthreshold*psfreffwhm) then begin
              ;;(SH May 28 2013)
              ;; old behaviour was to not convolve. Now we want to be
              ;; more conservative and stop making kernels at a shorter
              ;; wavelength. But to still have some convolution taking
              ;; place we calculate a pure gaussian convolution kernel
              
              ;; of course if the target PSF is very close to the
              ;; origin PSF (difference in fwhm < pixels) we still should
              ;; do nothing
              
              gaussian_kernel_fwhm = sqrt(psfreffwhm^2 - psffwhm^2)
           
              if gaussian_kernel_fwhm gt 2*pixelsize then begin
  ;;;; create the convolution kernel ala gordon 2008 + sandstrom
                 convolutionkernel = conv_cube_create_gaussian_convolutionkernel( $
                                     psfrefname, $
                                     kernelfwhm=gaussian_kernel_fwhm)
  ;;;; save the created kernel
                 sxaddpar,houtputconvolutionkernel,'WAVELNTH',wave[i],'Wavelength in microns'
                 sxaddpar,houtputconvolutionkernel,'X_POS',(psfsize+1)/2,'X position of PSF center in detector pixels'
                 sxaddpar,houtputconvolutionkernel,'Y_POS',(psfsize+1)/2,'Y position of PSF center in detector pixels'
                 writefits,convolutionkernelname,float(convolutionkernel),houtputconvolutionkernel
                 havekernel = 1
                 wanttoconvolve = 1
              endif else begin ;; really too broad to make anything useful
                 ;; write a file to say that we checked this psf and it is
                 ;; not useful
                 sxaddpar,houtputconvolutionkernelnotuseful,'WAVELNTH',wave[i],'Wavelength in microns'
                 sxaddpar,houtputconvolutionkernelnotuseful,'RFWHMMAX',fwhmthreshold,'Maximum ratio of size PSF/PSFREF'
                 writefits,convolutionkernelnotusefulname,[-1],houtputconvolutionkernelnotuseful
                 havekernel = 1
                 wanttoconvolve=0
              endelse ;; really too broad
           endif ;; broader then threshold
        endif ;; besensible
     endif ;; wanttocreatekernel

     if (wanttoconvolve eq 1) then begin
        if (havekernel eq 0) then begin
  ;;;; create the convolution kernel ala gordon 2008 + sandstrom
           convolutionkernel = conv_cube_create_convolutionkernel( $
                               tmp_psfname, $
                               psfrefname, $
                               fftpsfref=fftpsfref, $ ;; allows to store the previous calculation and not redo the work
                               hanning=hanning, $
                               psfreffwhm=psfreffwhm, $
                               nyquistthreshold=nyquistthreshold, $
                               psfthreshold=psfthreshold $
                                                                 )
  ;;;; save the created kernel
           sxaddpar,houtputconvolutionkernel,'WAVELNTH',wave[i],'Wavelength in microns'
           sxaddpar,houtputconvolutionkernel,'X_POS',(psfsize+1)/2,'X position of PSF center in detector pixels'
           sxaddpar,houtputconvolutionkernel,'Y_POS',(psfsize+1)/2,'Y position of PSF center in detector pixels'
           writefits,convolutionkernelname,float(convolutionkernel),houtputconvolutionkernel
           havekernel = 1
        endif

        convolutionkernel = readfits(convolutionkernelname,hconvolutionkernel)
        convolved_plane = conv_cube_convolve_image(reform(cube[*,*,i]),hcube,convolutionkernel,pixelsize)
        refwave = sxpar(hconvolutionkernel,'REFWAVE')
        sxaddpar,hcube,'HISTORY','CONV_CUBE: plane '+n2s(i)+' ('+swave[i]+' mum) convolved to '+n2s(refwave,format='(F10.3)')+' mum'
     endif else begin
        ;; do nothing psfs are too similar
        convolved_plane = cube[*,*,i]
     endelse
  ;;; store in output cube
     convolved_cube[*,*,i] = convolved_plane
     
  ;;;; some housekeeping
     if keeppsf then begin
        psfname = psfdirname+'/'+'psf_'+instrument+'_'+swave[i]+'.fits'
        foo = file_search(tmp_psfname,count=cnt)
        if cnt eq 1 then begin
           file_move,tmp_psfname,psfname
        endif
     endif
     file_delete,'stiny2.par',tmp_psfname,/allow_nonexistent
     
  endfor
  
  outstruct = {cube:convolved_cube,header:hcube,wave:wave,hwave:hwave}
  
  if n_elements(outfitsfilename) ne 0 then begin
     cleanoutfitsfilename = string(outfitsfilename[0])
     writefitscube,outstruct,cleanoutfitsfilename
  endif
  
  return,outstruct
  
end

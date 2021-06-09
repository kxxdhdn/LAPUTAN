; Hello!!!
;
; Please take 35 seconds to read carefully the following comments.
;
; The (main) routine convolve_image.pro will load an image, a convolution kernel, prepare the kernel, 
; prepare the image, convolve them, and write the result back into a fits file.
;
; It will load your images, kernels and write the result in the current directory, 
; althougt this can be adjusted by changing the path variabes (see below)
;
; If you have the image and kernel already loaded into the current IDL session, 
; you can convolve them by typing:
;
; .compile conv_image
; do_the_convolution,image,header,kernel_image,kernel_header,result_image,result_header,$
;                    result_kernel_image,result_kernel_header,do_we_write
;
; Where:
; original_image,original_header            are the the starting (original)    image and header
; kernel_image,kernel_header                are the the convolution kernel     image and header
; result_image,result_kernel                are the the final (result)         image and header
; result_kernel_image,result_kernel_header  are the the final (adapted) kernel image and header
; do_we_write should be either the number 0 if you do not want to have feedback in your screen
;   or the number 1 if you would like the routine to tell you the intermediate steeps... 
; (the resulting convolved image should be the same in either case, and error messages will always be displayed)
;
; At the begining of the main routine (pro conv_image) 4 variables are defined that you may want to change:
;
; do_we_write (idem as above)
;
; do_we_save_the_kernel should be either the number 0 if you do not want to store the adjusted kernel
;   or the number 1 if you want to save the adapted kernel (useful for sanity check and tests)
;
; images_path : set path to load/save the images relative to the crrent directory.
;  You may want to use something like:  images_path = './../Images/'
;
; kernels_path : set path to load the kernels relative to the crrent directory.
; You may want to do something like:  kernels_path = './../Kernels/'
;
; It should work with either 2D images or 3D cubes (like IRS data or Scanamorphos data).
; In the 3D case, the same kernel will be used to convolve every individual plane.
; (I have a slightly more complex version that load a family of kernels, useful
; when the PSF changes from frame to frame, or when the uncertanty image have correlations, 
; please let me know if you need it and I can send it to you)
;
; Please let me know if you find any bug or problem, if you would like to add some other feature
; to this routine, or you woud like to catch a beer next time we meet in a conference ;))
;
; Written by Gonzalo J. Aniano on Feb 2, 2011 (it coincided with my birthday...)
;
; Comments appreciated at ganiano@astro.princeton.edu
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro gaussian_kernel,kernel,sigma_kernel

; Dont worry, we will only use a gaussian kernel to smoothly interpolate over the 
; missing data (NAN) values, not for any kind of real convolution!!!

size_kernel  = 1 + 2 * fix(3.0*sigma_kernel)

size_kernel = size_kernel > 5

dist_sq = fltarr(size_kernel,size_kernel)

for index_x=0,size_kernel-1 do begin
   for index_y=0,size_kernel-1 do begin
      dist_sq[index_x,index_y] = ( float(index_x) - (float(size_kernel-1) / 2.0) ) ^2 + $
                                 ( float(index_y) - (float(size_kernel-1) / 2.0) ) ^2 
   endfor
endfor

kernel = exp (-dist_sq / ( 2*(sigma_kernel^2)))
kernel = kernel/total(kernel)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro remove_nan,image,header,places_with_original_data,do_we_write

if do_we_write eq 1 then print,'Checking for NaN values in the image.'
if do_we_write eq 1 then print,' '

places_with_original_data = image *0.0 +1.0

wh = where(finite(image) ne 1,cnt)
if cnt ne 0 then begin
  
   places_with_original_data(wh) = !Values.F_NAN
  
   if do_we_write eq 1 then print,'The image has NAN values, we will replace them during the convolution, and restore them later.'
   if do_we_write eq 1 then print,' '
   sigma_kernel = 2 ; pixels
   gaussian_kernel,kernel,sigma_kernel
   image(wh)=0
  
   check_FITS, image, header, dimen,/NOTYPE 
   if (N_elements(dimen) EQ 2) then dimen=[dimen[0], dimen[1],1]
  
   image_smooth = image * 0.0
   for nan_iter=1,5 do begin
      for frame=0,dimen[2]-1 do begin
         temp = image[*,*,frame]
         image_smooth [*,*,frame] = convolve (temp,kernel)
      endfor
     image(wh)=image_smooth(wh)
  endfor
  if do_we_write eq 1 then print,'The image has NAN values were replaced succesfully.'
  if do_we_write eq 1 then print,' '
  
endif else begin
  if do_we_write eq 1 then print,'The image do not have NAN values.'
  if do_we_write eq 1 then print,' '

endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro make_odd_square,image,header,do_we_write

; Pad the PSF with zeros into a square of odd number of pixels

check_FITS, image, header, dimen,/NOTYPE 
x_size_old = dimen[0] 
y_size_old = dimen[1]

size_new = x_size_old > y_size_old
if ((size_new mod 2) EQ 0) then size_new = size_new + 1

if (size_new gt x_size_old ) or (size_new gt y_size_old) then begin
   new_image = fltarr(size_new,size_new)
   new_image[0:x_size_old-1,0:y_size_old-1] = image
   image = new_image
   if do_we_write eq 1 then print,'The PSF was padded from a size ('+strtrim(string(x_size_old),2)+'x'+strtrim(string(y_size_old),2)+$
           ') into an odd square of size ('+strtrim(string(size_new),2)+'x'+strtrim(string(size_new),2)+') pixels.'
   if do_we_write eq 1 then print,' '
   sxaddpar,header,'NAXIS1',size_new
   sxaddpar,header,'NAXIS2',size_new

endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro get_maximun,image,x_max,y_max

rad_to_mean = 5

mean_im = image*0.0

for    i =-fix(rad_to_mean),fix(rad_to_mean) do begin
   for j =-fix(sqrt(rad_to_mean^2-i^2)),fix(sqrt(rad_to_mean^2-i^2)) do begin
      mean_im = mean_im + shift(image,i,j) 
   endfor
endfor

mx    = max(mean_im, location) 
index = ARRAY_INDICES(mean_im, location)  
x_max = index[0]
y_max = index[1]

wh = where (abs(mean_im - mx)/mx lt 5e-4, cnt)

if cnt gt 1 then begin
   print,'WARNING: The PSF has '+string(cnt)+' pixels with values similar to its maxximun...'
   print,'WARNING: we will take their centroid...'
   
   size_im = size(image)
   x_size = size_im[1]
   y_size = size_im[2]
   
   xpos = rebin(          (dindgen(x_size)),x_size,y_size)
   ypos = rebin( transpose(dindgen(y_size)),x_size,y_size)
   
   x_max = total(xpos(wh))/float(cnt)
   y_max = total(xpos(wh))/float(cnt)
endif


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro center_PSF,image,do_we_write

if do_we_write eq 1 then print,'Centering the PSF.'
if do_we_write eq 1 then print,' ' 

size_im      = (size(image))[1]
center_pixel = fix((size_im - 1) / 2)

get_maximun,image,x_max,y_max

; determine the needed shifts
shift_x = center_pixel - x_max
shift_y = center_pixel - y_max   

; make the shift if nonzero
if (shift_x ne 0) or (shift_y ne 0) then begin 
   if do_we_write eq 1 then print,'Shifting the PSF center by ('+ strtrim(string(shift_x),2) +','+ strtrim(string(shift_y),2) + ') pixels'
   image = shift(image,shift_x,shift_y)
   image[  0                     :abs(shift_x),*]=0.0
   image[  size_im-1-abs(shift_x):size_im-1   ,*]=0.0
   image[*,0                     :abs(shift_y)  ]=0.0
   image[*,size_im-1-abs(shift_y):size_im-1     ]=0.0
endif

; We check that the centering is OK:
get_maximun,image,x_max,y_max
shift_x = center_pixel - x_max
shift_y = center_pixel - y_max   
if (shift_x ne 0) or (shift_y ne 0) then if do_we_write eq 1 then print,'WARNING: Something went wrong in the image centering routine!!!'

if do_we_write eq 1 then print,'The PSF was centered successfully.'
if do_we_write eq 1 then print,' ' 

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro do_the_convolution,image,header,kernel_image,kernel_header,$
                       result_image,result_header,result_kernel_image,result_kernel_header,$
                       do_we_write

result_image  = image                  ; We make copies of the incoming images/headers so we can change
result_header = header                 ; then without changing the user original images/headers
result_kernel_image  = kernel_image    ; 
result_kernel_header = kernel_header   ;
                       
; First of all we replace the nan value with an interpolation of its neighbours, so 
; the convolution works better. After the convolution finishes, we will replace back
; the points that had NAN with NAN, but we will keep the image real during the convolution.

remove_nan,result_image,result_header,places_with_original_data,do_we_write

; Secondly we pad the images into a square with an even number of pixels per side.
; We add 100 Arcsec of black sky in each side to be able to include
; the boundaries contributions in the convolutions...

padding_arcseconds = 100.0

pixel_scale =                             fxpar(result_header,'PIXSCALE',count=count)
if (count EQ 0) then pixel_scale =        fxpar(result_header,'SECPIX'  ,count=count)
if (count EQ 0) then pixel_scale = sqrt( (fxpar(result_header,'CD1_1'   ,count=count)^2) + (fxpar(result_header,'CD1_2',count=count)^2) )*3600.0
if (count EQ 0) then pixel_scale = abs (  fxpar(result_header,'CDELT1'  ,count=count)   *3600.0)
if (count EQ 0) then begin
  print,'WARNING: we cannot get the pixel size in the image fits file header!!!!'
  print,' '
  print,'Please enter the pixel size of the image in arcsec.'
  print,'For example: 0.50'
  pixel_scale = 1.0
  read,': ',pixel_scale
  pixel_scale = pixel_scale > 0.001
  print,' '
endif else begin
   if do_we_write eq 1 then print,'The image has pixel of '+strtrim(string(pixel_scale,format='(F8.3)'),2)+' arcsec side.'
   if do_we_write eq 1 then print,' '
endelse
   
pixels_added_side = fix(padding_arcseconds / pixel_scale)
pixels_added      = 2 * pixels_added_side

check_FITS, result_image, result_header, dimen,/NOTYPE 

if (N_elements(dimen) EQ 2) then dimen=[dimen[0], dimen[1],1]
number_of_frames = dimen[2]

old_image_size_x = dimen[0] 
old_image_size_y = dimen[1]
new_image_size_x = old_image_size_x + pixels_added
new_image_size_y = old_image_size_y + pixels_added

if do_we_write eq 1 then print,'Padding the Original image for the convolution'
if do_we_write eq 1 then print,' '

padded_image = fltarr(new_image_size_x,new_image_size_y,dimen[2])

for index=0,number_of_frames-1 do begin
   padded_image [pixels_added_side:pixels_added_side+old_image_size_x-1,pixels_added_side:pixels_added_side+old_image_size_y-1,index] = result_image [*,*,index]
endfor

if do_we_write eq 1 then print,'The original image had size '+strtrim(old_image_size_x,2)+' x '+strtrim(old_image_size_y,2)+' pixels, and'
if do_we_write eq 1 then print,'was be padded to '+strtrim(new_image_size_x,2)+' x '+strtrim(new_image_size_y,2)+' pixels for the convolution.'
if do_we_write eq 1 then print,' '

; We need to put the psf and the original images in pixel of the same physical dimmensions.

if do_we_write eq 1 then print,'Adjusting the convolution kernel to the image resolution.'
if do_we_write eq 1 then print,' '

make_odd_square,result_kernel_image,result_kernel_header,do_we_write

size_ker = (size(result_kernel_image))[1]

if do_we_write eq 1 then print,'The kernel has '+strtrim(size_ker,2)+' x '+strtrim(size_ker,2)+' pixels.'
if do_we_write eq 1 then print,' '

pixel_scale_kernel = fxpar(result_kernel_header,'PIXSCALE',count=count)
if (count EQ 0) then pixel_scale_kernel = fxpar(result_kernel_header,'SECPIX',count=count)
if (count EQ 0) then pixel_scale_kernel = abs(fxpar(result_kernel_header,'CD1_1',count=count)*3600.0)
if (count EQ 0) then pixel_scale_kernel = abs(fxpar(result_kernel_header,'CDELT1',count=count)*3600.0)
if (count EQ 0) then begin
   print,'WARNING: we cannot get the pixel size in the image fits file header!!!!'
   print,' '
   print,'Please enter the pixel size of the image in arcsec.'
   print,'For example: 0.50'
   pixel_scale_kernel = 1.0
   read,': ',pixel_scale_kernel
   pixel_scale_kernel = pixel_scale_kernel > 0.001
   print,' '
endif 
    
if (abs(pixel_scale_kernel - pixel_scale)/pixel_scale) gt 0.05 then begin
   if do_we_write eq 1 then print,'The convolution kernel and the image are in grids of different pixel size,'
   if do_we_write eq 1 then print,'we will transform the kernel into the correct pixel size now.'
   if do_we_write eq 1 then print,' '
   size_ker = (size(kernel_image))[1]
   size_new = round( float(size_ker) * pixel_scale_kernel / pixel_scale ) 

  if ((size_new mod 2) EQ 0) then size_new = size_new + 1

  if do_we_write eq 1 then print,'Resampling the kernel from ' $ 
     +strtrim(size_ker,2)+' x '+strtrim(size_ker,2)+' to '$
     +strtrim(size_new,2)+' x '+strtrim(size_new,2)+' pixels.' 
  if do_we_write eq 1 then print,' '

  sxaddpar,result_kernel_header,'CD1_1' ,pixel_scale/3600.0
  sxaddpar,result_kernel_header,'CD1_2' ,0
  sxaddpar,result_kernel_header,'CD2_1' ,0
  sxaddpar,result_kernel_header,'CD2_2' ,pixel_scale/3600.0
  sxdelpar,result_kernel_header,'PIXSCALE' 
  sxdelpar,result_kernel_header,'SECPIX' 
  sxdelpar,result_kernel_header,'CDELT1' 
   
  result_kernel_image = congrid(result_kernel_image+0.0,size_new,size_new,cubic=-0.5,/center)
endif else begin
   if do_we_write eq 1 then print,'The convoution kernel and image are in the same pixel scale.'
   if do_we_write eq 1 then print,' '
endelse  
   
max_ker_size = new_image_size_x < new_image_size_y
if ((max_ker_size mod 2) EQ 0) then max_ker_size = max_ker_size -1

size_ker = (size(result_kernel_image))[1]

if max_ker_size lt size_ker then begin
      if do_we_write eq 1 then print,'Trimming the kernel from ' $ 
             +strtrim(size_ker,2)+' x '+strtrim(size_ker,2)+' to '$
             +strtrim(max_ker_size,2)+' x '+strtrim(max_ker_size,2)+' pixels.' 
   if do_we_write eq 1 then print,' '
   trim_side = fix(size_ker-max_ker_size)/2
   result_kernel_image = result_kernel_image (trim_side:trim_side+max_ker_size-1,trim_side:trim_side+max_ker_size-1)
endif

center_PSF,result_kernel_image,do_we_write

result_kernel_image = result_kernel_image/total(result_kernel_image)

size_ker = (size(result_kernel_image))[1]
sxaddpar,result_kernel_header,'NAXIS1',size_ker
sxaddpar,result_kernel_header,'NAXIS2',size_ker
    
if do_we_write eq 1 then print,'The convolution kernel has been adapted to the image resolution successfully.'
if do_we_write eq 1 then print,' '

check_FITS,result_kernel_image,result_kernel_header, dimen,/NOTYPE

if size_ker lt 3 then result_kernel_image = float([[0,0,0],[0,1,0],[0,0,0]])

if do_we_write eq 1 then print,'Convolving the image.'
if do_we_write eq 1 then print,' '

for index=0,number_of_frames-1 do begin
   if do_we_write eq 1 then print,'Convolving the frame '+strtrim(index+1,2)+' of '+strtrim(number_of_frames,2)+'.'
   padded_image [*,*,index] = convol_fft(padded_image [*,*,index],result_kernel_image)
endfor
if do_we_write eq 1 then print,'Convolution ready!!!'
if do_we_write eq 1 then print,' '

for index=0,number_of_frames-1 do begin
   result_image [*,*,index] = padded_image [pixels_added_side:pixels_added_side+old_image_size_x-1,pixels_added_side:pixels_added_side+old_image_size_y-1,index] 
endfor

result_image = result_image * places_with_original_data
check_FITS, result_image, result_header, dimen,/NOTYPE
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;'++++++++++++++++++++++++++++++++++++++++++++++++++++'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro convolve_image

;images_path  = './../Images/'
;kernels_path = './../Kernels/'
;convolutions_path = './../Convolved/'

images_path  = '/Users/dhu/Github/laputan/tests/out/conv/'
kernels_path = '/Users/dhu/Github/laputan/tests/lib/' 
convolutions_path = '/Users/dhu/Github/laputan/tests/out/conv/'
do_we_write           = 1
do_we_save_the_kernel = 1

if do_we_write eq 1 then print,'------------------------------------------------------------------------------------------'
if do_we_write eq 1 then print,' '



if do_we_write eq 1 then print,'------------------------------------------------------------------------------------------'
if do_we_write eq 1 then print,' '

do_the_convolution,image,header,kernel_image,kernel_header,$
                   result_image,result_header,result_kernel_image,result_kernel_header,do_we_write

fits_write,convolutions_path+Filename+'_conv.fits',result_image,result_header 
if do_we_save_the_kernel eq 1 then fits_write,convolutions_path+Filename+'_kernel.fits',result_kernel_image,result_kernel_header

if do_we_write eq 1 then print,'The image was convolved and saved successfully' 
if do_we_write eq 1 then print,' '
if do_we_write eq 1 then print,'------------------------------------------------------------------------------------------'

end


;; swarp drop-in replacement for hastrom. Being better flux
;; conservative hope fully

;; we use the full calling sequence:
;; hswarp,oldimage,oldheader,newimage,newheader,refheader

pro hswarp,oldimage,oldheader,newimage,newheader,refheader, $
           fixheader=fixheader, $
           debug=debug, $
           keepedge=keepedge

;; a temp directory to hold the files
  tmpdir = '/tmp/hswarp_'+strcompress(string(floor(randomu(time)*1d7)),/remove_all)+'/'
  file_mkdir,tmpdir
  
  pushd,tmpdir

  if keyword_set(fixheader) then begin
     refh=fitsheader_pc2cd(refheader)
     oldh=fitsheader_pc2cd(oldheader)
  endif else begin
     refh=refheader
     oldh=oldheader
  endelse
  
  ;;  writefits,'old.fits',float(oldimage),oldh
  writefits,'old.fits',float(oldimage),oldh

  openw,lun,'coadd.head',/get_lun
  printf,lun,refh
  close,lun
  free_lun,lun

  spawn,'swarp -d > swarp.cfg'
  swarp_options=' -c swarp.cfg -SUBTRACT_BACK N '
  if keyword_set(debug) then verbose='FULL' else verbose='QUIET'
  swarp_options=swarp_options+' -VERBOSE_TYPE '+verbose+' '
  
  spawn,'swarp '+swarp_options+'-RESAMPLING_TYPE LANCZOS3 '+' old.fits'
  newimage=readfits('coadd.fits',newheader,/silent)
  
  ;; add back in the edges because lanzcos3 kills the edges
  ;; do it in steps of less and less precision
  if keyword_set(keepedge) then begin

     oldweight=readfits('coadd.weight.fits',/silent)
     if total((oldweight eq 0d0)) ne 0 then begin
        spawn,'swarp '+swarp_options+' -RESAMPLING_TYPE LANCZOS2 '+' old.fits'
        edgeimage=readfits('coadd.fits',/silent)
        newweight=readfits('coadd.weight.fits',/silent)
        edgeidx=where((oldweight eq 0d0) and (newweight ne 0d0),cnt)
        if cnt ne 0 then newimage[edgeidx]=edgeimage[edgeidx]

        oldweight=readfits('coadd.weight.fits',/silent)
        if total((oldweight eq 0d0)) ne 0 then begin
           spawn,'swarp '+swarp_options+' -RESAMPLING_TYPE BILINEAR '+' old.fits'
           edgeimage=readfits('coadd.fits',/silent)
           newweight=readfits('coadd.weight.fits',/silent)
           edgeidx=where((oldweight eq 0d0) and (newweight ne 0d0),cnt)
           if cnt ne 0 then newimage[edgeidx]=edgeimage[edgeidx]

           oldweight=readfits('coadd.weight.fits',/silent)
           if total((oldweight eq 0d0)) ne 0 then begin
              spawn,'swarp '+swarp_options+' -RESAMPLING_TYPE NEAREST '+' old.fits'
              edgeimage=readfits('coadd.fits',/silent)
              newweight=readfits('coadd.weight.fits',/silent)
              edgeidx=where((oldweight eq 0d0) and (newweight ne 0d0),cnt)
              if cnt ne 0 then newimage[edgeidx]=edgeimage[edgeidx]
           endif

        endif

     endif

  endif

  popd

  ;; this is needed because swarp is conserving surface
  ;; brightness/pixel and the pixels size changes

  getrot,oldh,rot,oldcdelt
  getrot,refh,rot,refcdelt
  
  oldpixelfieldofview = abs(oldcdelt[0]*oldcdelt[1])
  newpixelfieldofview = abs(refcdelt[0]*refcdelt[1])

  newimage = newimage*oldpixelfieldofview/newpixelfieldofview

  if keyword_set(debug) then begin
     stop
  endif
  
  file_delete,/recursive,tmpdir

end

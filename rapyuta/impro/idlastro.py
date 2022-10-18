#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Python conversion of the IDL Astronomy Users Library
https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/

    hextract, hswarp

See also impro.irebin

"""

import os
from pathlib import Path
import subprocess as SP
import numpy as np

## Local
import rapyuta.utbox as UT
import rapyuta.inout as IO

def hextract(oldim, oldhd, x0, x1, y0, y1):
    '''
    Crop 2D image with pixel sequence numbers

    [REF] IDL lib hextract
    https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/hextract.pro
    '''
    newhd = oldhd.copy()
    # hdr['NAXIS1'] = x1 - x0 + 1
    # hdr['NAXIS2'] = y1 - y0 + 1
    newhd['CRPIX1'] += -x0
    newhd['CRPIX2'] += -y0
    newim = oldim[y0:y1+1, x0:x1+1]

    return newim, newhd

def hswarp(oldimage, oldheader, refheader,
           keepedge=False, tmpdir=None, verbose=True):
    '''
    Python version of hswarp (IDL), 
    a SWarp drop-in replacement for hastrom, 
    created by S. Hony

    ------ INPUT ------
    oldimage            ndarray
    oldheader           header object
    refheader           ref header
    keepedge            default: False
    tmpdir              default: None
    verbose             default: True
    ------ OUTPUT ------
    ds                  output object
      image               newimage
      header              newheader
    '''
    if verbose==False:
        devnull = open(os.devnull, 'w')
    else:
        devnull = None

    ## Initialize output object
    ds = type('', (), {})()

    ## Set path of tmp files
    UT.maketmp(os.getcwd()+'/tmp_hswarp')
    path_tmp = UT.fclean(tmpdir+'/coadd*')
    ## Make input
    IO.write_fits(path_tmp+'old', oldheader, oldimage)
    with open(path_tmp+'coadd.head', 'w') as f:
        f.write(str(refheader))

    ## Create config file
    SP.call('swarp -d > swarp.cfg',
            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
    ## Config param list
    swarp_opt = ' -c swarp.cfg -SUBTRACT_BACK N '
    if verbose=='quiet':
        swarp_opt += ' -VERBOSE_TYPE QUIET '
    ## Run SWarp
    SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS3 '+' old.fits',
            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
    coadd = IO.read_fits(path_tmp+'coadd')
    newimage = coadd.data
    newheader = coadd.header

    ## Add back in the edges because LANCZOS3 kills the edges
    ## Do it in steps of less and less precision
    if keepedge==True:
        oldweight = IO.read_fits(path_tmp+'coadd.weight').data
        if np.sum(oldweight==0)!=0:
            SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS2 '+' old.fits',
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
            edgeimage = IO.read_fits(path_tmp+'coadd').data
            newweight = IO.read_fits(path_tmp+'coadd.weight').data
            edgeidx = np.logical_and(oldweight==0, newweight!=0)
            if edgeidx.any():
                newimage[edgeidx] = edgeimage[edgeidx]

            oldweight = IO.read_fits(path_tmp+'coadd.weight').data
            if np.sum(oldweight==0)!=0:
                SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+' old.fits',
                        shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                edgeimage = IO.read_fits(path_tmp+'coadd').data
                newweight = IO.read_fits(path_tmp+'coadd.weight').data
                edgeidx = np.logical_and(oldweight==0, newweight!=0)
                if edgeidx.any():
                    newimage[edgeidx] = edgeimage[edgeidx]

                oldweight = IO.read_fits(path_tmp+'coadd.weight').data
                if np.sum(oldweight==0)!=0:
                    SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+' old.fits',
                            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                    edgeimage = IO.read_fits(path_tmp+'coadd').data
                    newweight = IO.read_fits(path_tmp+'coadd.weight').data
                    edgeidx = np.logical_and(oldweight==0, newweight!=0)
                    if edgeidx.any():
                        newimage[edgeidx] = edgeimage[edgeidx]

    ## SWarp is conserving surface brightness/pixel
    ## while the pixels size changes
    oldcdelt = IO.get_pc(wcs=IO.patch_wcs_3D(header=oldheader).wcs).cdelt
    refcdelt = IO.get_pc(wcs=IO.patch_wcs_3D(header=refheader).wcs).cdelt
    old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
    new_pixel_fov = abs(refcdelt[0]*refcdelt[1])
    newimage = newimage * old_pixel_fov/new_pixel_fov
    newimage[newimage==0] = np.nan
    # print('-------------------')
    # print(old_pixel_fov/new_pixel_fov)
    IO.write_fits(path_tmp+'new', newheader, newimage)
    # print('-------------------')
    
    ## Delete tmp file if tmpdir was not specified
    if tmpdir is None:
        UT.fclean(path_tmp)

    ds.data = newimage
    ds.header = newheader

    return ds

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Input & Output

    fclean, write_fits, read_fits,
    arr2tab, tab2arr, get_cd, get_pc, patch_wcs_3D,
    write_hdf5, read_hdf5,
    write_ascii, read_ascii, write_csv, read_csv

    fitsext, h5ext, ascext, csvext, savext

"""

import sys, os, logging
## Hide FITSFixedWarning:
## Removed redundant SCAMP distortion parameters
## because SIP parameters are also present [astropy.wcs.wcs]
logging.disable(sys.maxsize)

import subprocess as SP
from pathlib import Path
import warnings

import numpy as np
from astropy.io import fits#, registry
from astropy.wcs import WCS
from astropy import units as u
from specutils.io.registers import identify_spectrum_format
import h5py as H5
import csv

## Local
import rapyuta.utbox as UT

fitsext = '.fits'
h5ext = '.h5'
ascext = '.txt'
csvext = '.csv'
savext = '.sav'

## Quick file cleaning
def fclean(fname, *alert):
    '''
    Clean folder/file

    '''
    SP.call('rm -rf '+fname, shell=True)
    for text in alert:
        print(text)


##------------------------------------------------
##
##     Flexible Image Transport System (FITS)
##
##------------------------------------------------

def write_fits(fname, header, data,
               wave=None, wmod=0, header_w=None,
               ext=fitsext, **hdrl):
    '''
    Write fits file

    ------ INPUT ------
    fname               output FITS filename
    header              header of primary HDU
    data                data in primary HDU
    wave                data in table 1 (ndarray. Default: None)
    wmod                wave table format (0 - Image; 1 - BinTable. Default: 0)
    header_w            header of WAVE-TAB
    ------ OUTPUT ------
    '''
    for key, value in hdrl.items():
        header[key] = value
    primary_hdu = fits.PrimaryHDU(header=header, data=data)
    hdul = fits.HDUList(primary_hdu)
    
    ## Add table
    if wave is not None:
        ## Convert wave format
        if isinstance(wave, fits.fitsrec.FITS_rec):
            if wmod==0:
                wave = tab2arr(wave)
        else:
            Nw = len(wave)
            if wmod==1:
                wave = arr2tab(wave, header=header_w)
        ## Create table
        if wmod==0:
            hdu = fits.ImageHDU(data=wave, header=header_w, name='WAVE-TAB')
        elif wmod==1:
            hdu = fits.BinTableHDU(data=wave, header=header_w, name='WCS-TAB ')

        hdul.append(hdu)

    hdul.writeto(fname+ext, overwrite=True)
    
def read_fits(fname, wmod=0, instr=None, instr_auto=True,
              file_unc=None, ext=fitsext):
    '''
    Read fits file (auto detect dim)

    ------ INPUT ------
    fname               input FITS filename
    wmod                output wave mode (Default: 0)
                          0 - 1darray; 
                          1 - FITS_rec.
    instr               data format depends on instruments (Default: None)
                          None - primary array contains image/cube data
                                 table 1 contains wavelength data
                          'JWST s3d' - JWST spectral cube
    instr_auto          auto-identify FITS data format (Default: True)
                          - only works when instr is None
    file_unc            unc filename (Default: None)
                          None - auto detection (ext should be fitsext)
    ------ OUTPUT ------
    ds                  output dataset
      HDUL                header data unit list
      header              header of data HDU
      data                data in data HDU
      header_w            header of WAVE-TAB
      wave                data in table 1 (None if does not exist)
      unc                 uncertainty array
    '''
    ## Initialize output object
    ds = type('', (), {})()
    ds.header_w = None
    ds.wave = None
    ds.unc = None

    ## Check for valid string input
    filename = fname+ext
    if not isinstance(filename, (str, Path)) or not os.path.isfile(filename):
        UT.strike('read_fits', f'{filename} is not a valid string path to a file',
                  cat='InputError')

    ## Identify data format
    if instr_auto and (instr is None):
        instr = identify_spectrum_format(filename)
    
    with fits.open(filename) as hdul:
        ds.HDUL = hdul
        if instr=='JWST s3d':
            ## Read header & data
            hdr = hdul[1].header
            ds.data = hdul[1].data

            ## JWST cubes have no ready-made WAVE-TAB
            if hdr['NAXIS']==3:
               w0 = hdr['CRVAL3']
               dw = hdr['CDELT3']
               wave = []
               for iw in range(hdr['NAXIS3']):
                   wave.append(w0+iw*dw)
               Nw = len(wave)
               if wmod==0:
                   wave = np.array(wave)
               elif wmod==1:
                   wave = arr2tab(wave, header=ds.header_w)

            ## Uncertainties
            ds.unc = hdul[2].data

        else:
            ## Read header & data
            hdr = hdul[0].header
            ds.data = hdul[0].data

            ## Read wavelength
            if len(hdul)==2:
                ds.header_w = hdul[1].header
                wave = hdul[1].data

                if isinstance(hdul[1], fits.BinTableHDU):
                    if wmod==0:
                        wave = tab2arr(wave) ## Convert FITS_rec to 1darray
                elif isinstance(hdul[1], fits.ImageHDU):
                    Nw = len(wave)
                    if wmod==1:
                        wave = arr2tab(wave, header=ds.header_w)

            ## Uncertainties
            if (file_unc is None) and ext==fitsext:
                file_unc = fname+'_unc'+ext
            if Path(file_unc).exists():
                ## Read uncertainty data
                with fits.open(file_unc) as hdul:
                    ds.unc = hdul[0].data

    ## WCS
    try:
        ds.wcs = WCS(hdr)
    except:
        ds.wcs = patch_wcs_3D(header=hdr, ext=fitsext)
        warnings.warn('2D WCS extracted from 3D data! ')

    ds.header = hdr
    ds.wave = wave
            
    return ds

def arr2tab(arr, header=None, unit='um'):
    '''
    Convert 1darray to FITS_rec
    '''
    N = len(arr)
    arr = np.array(arr).reshape((N,1))
    col = fits.Column(array=[arr], format=str(N)+'E',
                      name='WAVE-TAB', unit=unit, dim='(1,{})'.format(N))
    tab = fits.BinTableHDU.from_columns([col], header=header, name='WCS-TAB ')

    return tab.data

def tab2arr(tab):
    '''
    Convert FITS_rec to 1darray
    '''
    return tab[0][0][:,0]

def get_cd(pc=None, cdelt=None, header=None, wcs=None):
    '''
    Convert CDELTia + PCi_ja to CDi_ja (2D only)
    (astropy.wcs use PC/CDELT by default)

    ------ INPUT ------
    pc                  PC matrix (priority if co-exist)
    cdelt               Coordinate increment at ref point
    header              header object (2nd priority if co-exist)
    wcs                 WCS object (3rd priority)
    ------ OUTPUT ------
    ds                  output object
      cd                  CD matrix
      pc                  PC matrix
      cdelt               CDELTia
    '''
    ## Initialize output object
    ds = type('', (), {})()
    ds.cd = np.zeros((2,2))
    ds.pc = np.zeros((2,2))
    ds.cdelt = np.zeros(2)

    if pc is not None and cdelt is not None:
        ds.pc = pc
        ds.cdelt = cdelt
        ## CDi_j = PCi_j * CDELTi
        ds.cd = ds.pc * ds.cdelt.reshape((2,1))
    else:
        if header is not None:
            # w = WCS(header)
            w = patch_wcs_3D(header=header).wcs # force 2D
        else:
            if wcs is not None:
                w = wcs
            else:
                UT.strike('get_cd', 'no input.', cat='InputError')

        ds.cd = w.pixel_scale_matrix

        if w.wcs.has_pc():
            ds.pc = w.wcs.get_pc()
            ds.cdelt = w.wcs.get_cdelt()
        else:
            ## See astropy.wcs.utils.proj_plane_pixel_scales
            ds.cdelt = np.sqrt((ds.cd**2).sum(axis=0, dtype=float))
            ## See Calabretta&Greisen paper sec-6.2 [A&A 395, 1077-1122 (2002)]
            ds.cdelt[0] = -ds.cdelt[0]
            ds.pc = ds.cd / ds.cdelt.reshape((np.size(ds.cdelt),1))

    return ds

def get_pc(cd=None, header=None, wcs=None):
    '''
    Convert CDi_ja to CDELTia + PCi_ja (2D only)

    ------ INPUT ------
    cd                  CD matrix (priority if co-exist)
    header              header object (2nd priority if co-exist)
    wcs                 WCS object (3rd priority)
    ------ OUTPUT ------
    ds                  output object
      cd                  CD matrix
      pc                  PC matrix
      cdelt               CDELTia
    '''
    ## Initialize output object
    ds = type('', (), {})()
    ds.cd = np.zeros((2,2))
    ds.pc = np.zeros((2,2))
    ds.cdelt = np.zeros(2)

    if cd is not None:
        ds.cd = cd
        ## See astropy.wcs.utils.proj_plane_pixel_scales
        ds.cdelt = np.sqrt((cd**2).sum(axis=0, dtype=float))
        ## See Calabretta&Greisen paper sec-6.2 [A&A 395, 1077-1122 (2002)]
        ds.cdelt[0] = -ds.cdelt[0]
        ## CDi_j = PCi_j * CDELTi
        ds.pc = ds.cd / ds.cdelt.reshape((2,1))
    else:
        if header is not None:
            # w = WCS(header)
            w = patch_wcs_3D(header=header).wcs # force 2D
        else:
            if wcs is not None:
                w = wcs
            else:
                UT.strike('get_pc', 'no input.', cat='InputError')

        ds.cd = w.pixel_scale_matrix

        if w.wcs.has_pc():
            ds.pc = w.wcs.get_pc()
            ds.cdelt = w.wcs.get_cdelt()
        else:
            ## See astropy.wcs.utils.proj_plane_pixel_scales
            ds.cdelt = np.sqrt((ds.cd**2).sum(axis=0, dtype=float))
            ## See Calabretta&Greisen paper sec-6.2 [A&A 395, 1077-1122 (2002)]
            ds.cdelt[0] = -ds.cdelt[0]
            ds.pc = ds.cd / ds.cdelt.reshape((np.size(ds.cdelt),1))

    return ds

def patch_wcs_3D(fname=None, header=None, del_kw=None, ext=fitsext):
    '''
    Auto-detect & reduce dim if WCS is 3D with distortion

    ------ INPUT ------
    fname               target FITS filename
    header              header object
    del_kw              kw to delete
                          None: reduce dimension with naxis kwarg (Default)
                          'all': All kw with the character "3"
                          'sip': SIP kw
    ------ OUTPUT ------
    ds                  output object
      header              header of primary HDU
      wcs                 2D WCS
      was3d               True: if input data is 3D
    '''
    ## Initialize output object
    ds = type('', (), {})()
    ds.wcs = WCS(None)
    ds.header = None
        
    ## Read file/header
    if fname is not None:
        hdr = fits.open(fname+ext)[0].header
        header = hdr.copy()
    else:
        if header is not None:
            hdr = header.copy()
        else:
            UT.strike('patch_wcs_3D', 'no input.', cat='InputError')

    if header['NAXIS']==3:
        ds.was3d = True
    else:
        ds.was3d = False

    ## Extract WCS from 3D data
    if del_kw=='all':
        ## Opt.1: delete keywords containing "3" (not recommanded)
        for kw in hdr.keys():
            if '3' in kw:
                del header[kw]
        header['NAXIS'] = 2
        header['COMMENT'] = '3D WCS extraction: All keywords with the character "3" are deleted. '
        
        ds.wcs = WCS(header)
    elif del_kw=='sip':
        ## Opt.2: delete SIP keywords (not recommanded)
        ## https://fits.gsfc.nasa.gov/registry/sip/SIP_distortion_v1_0.pdf
        for kw in hdr.keys():
            if ('A_' in kw) and (not 'PA' in kw) and (not 'RA' in kw):
                del header[kw]
            if ('B_' in kw) or ('AP_' in kw) or ('BP_' in kw):
                del header[kw]
        if 'CTYPE3' in hdr.keys():
            del header['CTYPE3']
        header['COMMENT'] = '3D WCS extraction: SIP keywords are deleted. '
        
        ds.wcs = WCS(header)
    else:
        ## Default: reduce dimension (default)
        if 'CTYPE3' in hdr.keys():
            del header['CTYPE3']
        if 'NAXIS3' in hdr.keys():
            del header['NAXIS3']
        header['COMMENT'] = '3D WCS extraction: Shrink to 2D with naxis kwarg. '

        ds.wcs = WCS(header, naxis=2)

    ds.header = header

    return ds


##------------------------------------------------
##
##                     HDF5
##
##------------------------------------------------

def write_hdf5(fname, name, data, group='/',
               ind1=None, ind2=None, ind3=None, ind4=None,
               ind5=None, ind6=None, ind7=None,
               append=False, verbose=False, ext=h5ext):
    '''
    Write dataset into a h5 file (a single name/data_array per time, dim < 7)
    Inspired by SwING inout library

    ------ INPUT ------
    fname               output h5 filename
    name                name of the dataset (len <= 80)
    data                dataset (dim < 7)
    group               name of the group (Default: '/')
    indx                array index ([idimx_inf,idimx_sup] or idimx if data is scalar, Default: None)
    append              True: if not overwrite (Default: False)
    verbose             courtesy notification (Default: False)
    ------ OUTPUT ------
    '''
    ## Preliminaries
    ##---------------
    if (append):
        hf = H5.File(fname+ext, 'a')
    else:
        hf = H5.File(fname+ext, 'w')
    
    h5 = hf.require_group(group) # Default: group = '/'

    ## Check if it is a subarray
    subarr = ( (ind1 is not None) or (ind2 is not None) or (ind3 is not None) or \
               (ind4 is not None) or (ind5 is not None) or (ind6 is not None) )

    ## Convert all data to array format
    darr = np.array(data)
    sharr = darr.shape
    if darr.dtype.kind == 'U':
        ## create_dataset does not support lists of UTF-8 yet
        ## "these strings are supposed to store only ASCII-encoded text"
        ## See http://docs.h5py.org/en/stable/strings.html
        asciiList = [n.encode('ascii','ignore') for n in darr.flatten()]
        darr = np.array(asciiList).reshape(sharr)

    ## Write dataset (dset)
    ##----------------------
    
    ## Case 0: check if the dataset already exists
    try:
        dset = h5[name]
        hasdset = True
    except:
        hasdset = False
    if (hasdset):
        createdset = False
        ## No indx input
        if (not subarr):
            del h5[name]
            createdset = True
    else:
        createdset = True
        if (subarr):
            h5.close()
            UT.strike('write_hdf5', 'Dataset does not exist.',
                      cat='InputError')

    ## Case 1: no dataset OR no subarr
    if (createdset):
        dset = h5.create_dataset(name, data=darr)

    ## Case 2: sub-array filling
    if (subarr):
        ## Dimension of the array
        Ndim = len(sharr)

        ## Indices indx[idimx_inf,idimx_sup]
        if (Ndim > 0):
            if ind1 is None:
                ind1 = [0,sharr[0]]
            elif (np.size(ind1) == 1):
                ind1 = [ind1, ind1+1]
        if (Ndim > 1):
            if ind2 is None:
                ind2 = [0,sharr[1]]
            elif (np.size(ind2) == 1):
                ind2 = [ind2, ind2+1]
        if (Ndim > 2):
            if ind3 is None:
                ind3 = [0,sharr[2]]
            elif (np.size(ind3) == 1):
                ind3 = [ind3, ind3+1]
        if (Ndim > 3):
            if ind4 is None:
                ind4 = [0,sharr[3]]
            elif (np.size(ind4) == 1):
                ind4 = [ind4, ind4+1]
        if (Ndim > 4):
            if ind5 is None:
                ind5 = [0,sharr[4]]
            elif (np.size(ind5) == 1):
                ind5 = [ind5, ind5+1]
        if (Ndim > 5):
            if ind6 is None:
                ind6 = [0,sharr[5]]
            elif (np.size(ind6) == 1):
                ind6 = [ind6, ind6+1]
        if (Ndim > 6):
            if ind7 is None:
                ind7 = [0,sharr[6]]
            elif (np.size(ind7) == 1):
                ind7 = [ind7, ind7+1]

        ## Write the sub-array if indx are set
        if Ndim == 0:
            dset = arr[0]
        elif Ndim == 1:
            dset[ind1[0]:ind1[1]] = np.reshape( darr, (ind1[1]-ind1[0],) )
        elif Ndim == 2:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0]) )
        elif Ndim == 3:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0],ind3[0]-ind3[1]) )
        elif Ndim == 4:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1],ind4[0]:ind4[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1]) )
        elif Ndim == 5:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                 ind4[0]:ind4[1],ind5[0]:ind5[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1], \
                                   ind5[0]-ind5[1]) )
        elif Ndim == 6:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                 ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1], \
                                   ind5[0]-ind5[1],ind6[0]-ind6[1]) )
        elif Ndim == 7:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                 ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1],ind7[0]:ind7[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1], \
                                   ind5[0]-ind5[1],ind6[0]-ind6[1],ind7[0]-ind7[1]) )

    hf.flush()
    hf.close()
    
    ## Courtesy notification
    ##-----------------------
    if (verbose):
        if (append):
            if (subarr):
                print('[write_hdf5] Dataset {}/{} in the file:'.format(group,name))
                print('    ', fname+ext, ' has been modified.')
            else:
                print('[write_hdf5] Dataset {}/{} has been added in the file:'.format(group,name))
                print('    ', fname+ext, '.')
        else:
            print('[write_hdf5] Dataset {}/{} has been written in the new file:'.format(group,name))
            print('    ', fname+ext, '.')

def read_hdf5(fname, name, group='/',
              ind1=None, ind2=None, ind3=None, ind4=None,
              ind5=None, ind6=None, ind7=None,
              ext=h5ext):
    '''
    Read h5 file (a single name/data_array per time, dim < 7)
    Inspired by SwING inout library

    ------ INPUT ------
    fname               input h5 filename
    name                name of the dataset (len <= 80)
    group               name of the group (Default: '/')
    indx                array index ([idimx_inf,idimx_sup] or idimx if data is scalar, Default: None)
    ------ OUTPUT ------
    dset                dataset
    '''
    ## Preliminaries
    ##---------------
    hf = H5.File(fname+ext, 'r')
    h5 = hf.require_group(group) # Default: group = '/'

    Ndim = h5[name].ndim
    shapIN = h5[name].shape

    ## Read dataset (dset)
    ##----------------------
    
    ## Indices indx[idimx_inf,idimx_sup]
    if (Ndim > 0):
        if ind1 is None:
            ind1 = [0,shapIN[0]]
        elif (np.size(ind1) == 1):
            ind1 = [ind1, ind1+1]
    if (Ndim > 1):
        if ind2 is None:
            ind2 = [0,shapIN[1]]
        elif (np.size(ind2) == 1):
            ind2 = [ind2, ind2+1]
    if (Ndim > 2):
        if ind3 is None:
            ind3 = [0,shapIN[2]]
        elif (np.size(ind3) == 1):
            ind3 = [ind3, ind3+1]
    if (Ndim > 3):
        if ind4 is None:
            ind4 = [0,shapIN[3]]
        elif (np.size(ind4) == 1):
            ind4 = [ind4, ind4+1]
    if (Ndim > 4):
        if ind5 is None:
            ind5 = [0,shapIN[4]]
        elif (np.size(ind5) == 1):
            ind5 = [ind5, ind5+1]
    if (Ndim > 5):
        if ind6 is None:
            ind6 = [0,shapIN[5]]
        elif (np.size(ind6) == 1):
            ind6 = [ind6, ind6+1]
    if (Ndim > 6):
        if ind7 is None:
            ind7 = [0,shapIN[6]]
        elif (np.size(ind7) == 1):
            ind7 = [ind7, ind7+1]
        
    # Read the array or the sub-array if some indx are set
    if Ndim==0:
        dset = h5[name]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==1:
        dset = h5[name][ind1[0]:ind1[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==2:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==3:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==4:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==5:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1],ind5[0]:ind5[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==6:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==7:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1], \
                        ind7[0]:ind7[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')

    hf.close()

    return dset


##------------------------------------------------
##
##                     ASCII
##
##------------------------------------------------

def write_ascii(fname, header=None, dset=None, trans=False,
                append=False, comment=None, ext=ascext):
    '''
    Write ASCII file
    Supported format: commented_header
    See also: astropy.io.ascii.write, numpy.savetxt
    
    ------ INPUT ------
    fname               input ASCII filename
    header              data header
    dset                dataset (1darray[nrow], 2darray[nrow,ncol] or list)
    trans               transpose dset (Default: False)
    append              True: if not overwrite (Default: False)
    comment             single line comment on top
    ------ OUTPUT ------
    '''
    if append==True:
        mod = 'a'
    else:
        mod = 'w'

    with open(fname+ext, mod) as f:
        pass
        ## Comment on the top
        ##--------------------
        if comment is not None:
            f.write('## '+comment+'\n')
        
        ## Number of columns/rows
        ncol = 0
        nrow = 0
        ## Dataset
        ##---------
        rows = ''
        if dset is not None:
            dset = np.array(dset)
            ## Transpose dset
            if trans==True:
                Ny, Nx = dset.shape
                darr = []
                for x in range(Nx):
                    for y in range(Ny):
                        darr.append(dset[y,x])
                darr = np.array(darr).reshape((Nx,Ny))
            else:
                darr = dset
            ## Update number of columns/rows
            Ndim = darr.ndim
            if Ndim==1:
                ncol = 1
                nrow = len(darr)
            else:
                nrow, ncol = darr.shape
            ## Convert dataset array to string
            if Ndim==1 or Ndim==2:
                if nrow<500:
                    rows += str(darr).replace('[','').replace(']','')+'\n'
                else:
                    ## Python 500 row limit
                    for i in range(nrow//500+1):
                        rows += str(darr[500*i:500*(i+1)]).replace('[','').replace(']','')+'\n'
            else:
                UT.strike('write_ascii', f'Expected 1D or 2D array, got {Ndim:.d}D array instead.',
                          cat='ValueError')
        ## Header
        ##--------
        hrow = ''
        if header is not None:
            if isinstance(header, str):
                header = [header]
            else:
                header = list(header)
                
            for i, h in enumerate(header):
                if i==0:
                    hrow += h
                else:
                    hrow += ' '+h
            hrow += '\n'
            ## Update number of columns/rows
            ncol = len(header)
            nrow += 1
        ## Column width control
        ##----------------------
        ## Table in a single string
        if nrow==0:
            fwrows = ''
        else:
            cells = hrow.split()
            cells.extend(rows.split())
            cells = np.array(cells).reshape((-1,ncol))
            col_width = []
            for x in range(ncol):
               col_width.append(len(max(cells[:,x], key=len)))
            ## Table with fixed column widths in a single string
            if hrow=='':
                fwrows = '' # no header
            else:
                fwrows = '# '
            for y in range(nrow):
                for x in range(ncol):
                    c = cells[y,x]
                    ## '{message:{fill}{align}{width}}'.format(
                    ##     message=c, fill=' ', align='>', width=col_width+1)
                    if y==0 and x==0 and header is not None:
                        ## Due to hashtag
                        fwrows += ''.ljust(col_width[x]-len(c))+c
                    else:
                        fwrows += ''.ljust(col_width[x]+2-len(c))+c
                    if x==ncol-1:
                        fwrows += '\n'

        f.write(fwrows)

def read_ascii(fname, dtype=str, start_header=-1, ext=ascext):
    '''
    Read ASCII file
    Supported format: commented_header
    See also: astropy.io.ascii.read, numpy.genfromtxt

    ------ INPUT ------
    fname               input ASCII filename
    dtype               data type (Default: 'str')
    start_header        line number of the header line (Default: -1)
    ------ OUTPUT ------
    dset                dataset (dict)
    '''
    with open(fname+ext, 'r') as f:
        ## f.read() -> str | f.readlines() -> list
        header = None
        darr = []
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            # print(line)
            if line[0]=='#':
                ## By default, header is taken from the last commented line before data
                if isinstance(start_header, int)==True:
                    if start_header==-1:
                        header = line.split()[1:]
                    elif start_header>0 and i==start_header-1:
                        header = line.split()[1:]
            else:
                line = list(map(dtype, line.split()))
                row = []
                for cell in line:
                    row.append(cell)
                darr.append(row)

    darr = np.array(darr)

    ## Number of columns/rows
    nrow, ncol = darr.shape

    ## Default header = ['col1','col2',...]
    header_default = []
    for x in range(ncol):
        header_default.append('col'+str(x+1))

    ## Dict out
    dset = {}
    for x, h in enumerate(header_default): 
        dset[h] = darr[:,x]
        if header is not None:
            dset[header[x]] = darr[:,x]

    return dset

def write_csv(fname, header, dset, trans=False, append=False, ext=csvext):
    '''
    Read fits file
    See also: astropy.io.ascii.write

    ------ INPUT ------
    fname               output csv filename
    header              data labels in list('label1', 'label2', ...)
    dset                dataset (1darray[nrow], 2darray[nrow,ncol] or list)
    trans               transpose dset (Default: False)
    append              True: if not overwrite (Default: False)
    ------ OUTPUT ------
    '''
    if append==True:
        mod = 'a'
    else:
        mod = 'w'

    with open(fname+ext, mod, newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)

        ## Header
        ##--------
        writer.writeheader()

        ## Dataset
        ##---------
        dset = np.array(dset)
        ## Transpose dset
        if trans==True:
            Ny, Nx = dset.shape
            darr = []
            for x in range(Nx):
                for y in range(Ny):
                    darr.append(dset[y,x])
            darr = np.array(darr).reshape((Nx,Ny))
        else:
            darr = dset
        ## Number of columns/rows
        Ndim = darr.ndim
        if Ndim==1:
            ncol = 1
            nrow = darr.shape[0]
        else:
            nrow, ncol = darr.shape
        ## Write
        for y in range(nrow):
            ## Init dict
            rows = {h: [] for h in header}
            ## Write dict
            for x in range(ncol):
                rows[header[x]] = darr[y,x]
            ## Write csv row
            writer.writerow(rows)
            
def read_csv(fname, ext=csvext, *header):
    '''
    Read csv file

    ------ INPUT ------
    fname               input csv filename
    header              labels of data to read
    ------ OUTPUT ------
    dset                dataset (dict)
    '''
    with open(fname+ext, 'r', newline='') as csvfile:
        dset = {}
        for h in header:
            csvfile.seek(0) # reset pointer
            reader = csv.DictReader(csvfile) # pointer at end
            col = []
            for row in reader:
                col.append(row[h])
            data = np.array(col)
            dset[h] = data

    return dset

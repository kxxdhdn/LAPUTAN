#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Astrometry

    pix2sr, sr2arcsec2, rad2arcsec, hour2deg, deg2hour,
    get_cd, get_pc, fixwcs

See also numpy.deg2rad, etc.
"""

import sys, logging
logging.disable(sys.maxsize)

import math
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

def pix2sr(X, CDELT):
    '''
    X pixel = Y sr
    ------ INPUT ------
    X                   float or ndarray in pixel
    CDELT               float or ndarray (same dim if X is also ndarray)
    '''
    PFOV = abs(CDELT)
    
    return X * (PFOV * 2. * math.pi / 360.)**2.

def sr2arcsec2(X):
    '''
    X sr = Y arcsec^2
    '''
    return X * (360. * 3600. / (2. * math.pi))**2.

def rad2arcsec(X):
    '''
    X rad = Y arcsec
    '''
    return X * 360. * 3600. / (2. * math.pi)

def hour2deg(h, m, s, deg, arcmin, arcsec):
    '''
    Hour angle to degree conversion
    '''
    ra = (h + m/60. + s/3600.) * 360./24.
    if deg<0:
        dec = -(-deg + arcmin/60. + arcsec/3600.)
    else:
        dec = deg + arcmin/60. + arcsec/3600.

    return ra, dec

def deg2hour(ra, dec):
    '''
    Degree to hour angle conversion
    '''
    h = math.floor(ra * 24./360.)
    m = math.floor((ra*24./360. - h) * 60.)
    s = ((ra*24./360. - h)*60 - m) * 60.
    if dec<0:
        deg = math.ceil(dec) # negtive
        arcmin = -math.ceil((dec - deg) * 60.) # positive
        arcsec = -((dec - deg)*60. + arcmin) * 60. # positive
    else:
        deg = math.floor(dec)
        arcmin = math.floor((dec - deg) * 60.)
        arcsec = ((dec - deg)*60. - arcmin) * 60.
        
    print('{:d}h{:d}m{:04.2f}s, {:d}d{:d}m{:04.2f}s'.format(h,m,s,deg,arcmin,arcsec))

    return h, m, s, deg, arcmin, arcsec

def get_cd(pc=None, cdelt=None, header=None, wcs=None):
    '''
    Convert CDELTia + PCi_ja to CDi_ja
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
            w = WCS(header)
        else:
            if wcs is not None:
                w = wcs
            else:
                raise ValueError('No input!')

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
    Convert CDi_ja to CDELTia + PCi_ja

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
            w = WCS(header)
        else:
            if wcs is not None:
                w = wcs
            else:
                raise ValueError('No input!')

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

def fixwcs(file=None, header=None, mode='red_dim'):
    '''
    Auto-detect & reduce dim if WCS is 3D with distortion

    ------ INPUT ------
    file                FITS file (priority if co-exist)
    header              header object
    mode                ways to fix 3D WCS with SIP distortion
                          'red_dim' (Default)
                          'del_sip'
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
    if file is not None:
        hdr = fits.open(file)[0].header
        header = hdr.copy()
    else:
        if header is not None:
            hdr = header.copy()
        else:
            raise ValueError('No input!')

    if header['NAXIS']==3:
        ds.was3d = True
    else:
        ds.was3d = False

    ## Fix 3D WCS with SIP distortion
    if mode=='red_dim':
        ## Opt.1: reduce dimension
        for kw in hdr.keys():
            if '3' in kw:
                del header[kw]
        header['NAXIS'] = 2
        header['COMMENT'] = "Removed 3D keywords. "
        
        ds.wcs = WCS(header, naxis=2)
    elif mode=='del_sip':
        ## Opt.2: delete SIP keywords
        for kw in hdr.keys():
            if ('A_' in kw) and (not 'PA' in kw) and (not 'RA' in kw):
                del header[kw]
            if ('B_' in kw) or ('AP_' in kw) or ('BP_' in kw):
                del header[kw]
        header['COMMENT'] = "Removed SIP keywords. "

        ds.wcs = WCS(header)
    
    ds.header = header # fixed header

    return ds


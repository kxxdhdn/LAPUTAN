#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Calibration

    intercalib:
        read_filter, synthetic_photometry, correct_spec

"""

import os
import math
import numpy as np
# from astropy.io import ascii
from matplotlib.ticker import ScalarFormatter, NullFormatter
import subprocess as SP
import warnings

## Local
from rapyuta import rapyroot
from inout import (
    ascext, fitsext, h5ext,
    read_fits, write_fits, patch_wcs_3D,
    read_hdf5, write_hdf5#, read_ascii
)
from latte import listize, pix2sup, sup2pix
from plots import pplot

##-----------------------------------------------
##
##            <intercalib> based tools
##
##-----------------------------------------------

class intercalib:
    '''
    Intercalibration

    ------ INPUT ------
    filIN               target FITS file (Default: None)
    '''
    def __init__(self, filIN=None):
        
        ## INPUTS
        self.filIN = filIN

        if filIN is not None:
            # self.hdr = fixwcs(filIN+fitsext).header
            w = fixwcs(filIN+fitsext).wcs
            ds = read_fits(filIN)
            self.hdr = ds.header
            self.im = ds.data
            self.wvl = ds.wave

    def read_filter(self, filt, w_spec=None):
        '''
        Return center wavelength of the filters
        The input offset corresponds to the integrated broad band value (bboff).
        Suppose a flat spectral offset to correct the spectral value (specoff).
        Return specoff/bboff

        ------ INPUT ------
        filt                photometry names (string, tuple or list)
        w_spec              wavelengths (Default: None - via filIN)
        ------ OUTPUT ------
        self
          wcen                center wavelength
          specoff_ov_bboff    spectral/broad band offset ratio
        '''
        ## Convert all format phot names to list
        filt = listize(filt)

        if self.filIN is not None:
            w_spec = self.wvl
        for phot in filt:
            w_grid = read_hdf5(rapyroot+'/lib/data/filt_'+phot,
                               'Filter wavelength (microns)')
            if w_spec[0]>w_grid[0] or w_spec[-1]<w_grid[-1]:
                warnings.warn('Synthetic photometry of {} can be underestimated' \
                              'due to uncovered wavelengths. '.format(phot))
        ## Insert 2 wvl (0.01 um & w_spec[0]-0.01 um) with 0 value
        wave = np.insert(w_spec, 0, (.01, w_spec[0]-.01))
        Fnu_uni = np.ones(len(wave))
        Fnu_uni[:2] = 0
        flux = Fnu_uni[:,np.newaxis,np.newaxis]

        ## Write input.h5
        ##----------------
        fortIN = os.getcwd()+'/synthetic_photometry_input'
        
        write_hdf5(fortIN, 'Filter label', filt)
        write_hdf5(fortIN, 'Wavelength (microns)', wave, append=True)
        write_hdf5(fortIN, 'Flux (x.Hz-1)', flux, append=True)
        write_hdf5(fortIN, '(docalib,dophot)', [1,1], append=True)

        ## Call the Fortran lib
        ##----------------------
        SP.call('synthetic_photometry', shell=True)

        ## Read output.h5
        ##----------------
        fortOUT = os.getcwd()+'/synthetic_photometry_output'

        self.wcen = read_hdf5(fortOUT, 'Central wavelength (microns)')
        Fnu_filt = read_hdf5(fortOUT, 'Flux (x.Hz-1)')[:,0,0]
        self.specoff_ov_bboff = 1. / Fnu_filt

        ## Clean temperary h5 files
        ##--------------------------
        SP.call('rm -rf '+fortIN+h5ext, shell=True, cwd=os.getcwd())
        SP.call('rm -rf '+fortOUT+h5ext, shell=True, cwd=os.getcwd())

    def synthetic_photometry(self, filt, w_spec=None, Fnu_spec=None,
                             xscale=1, yscale=1, 
                             extrapoff=True, verbose=False):
        '''
        External Fortran library (SwING) needed

        ------ INPUT ------
        filt                photometry names (string, tuple or list)
        w_spec              wavelengths (Default: None - via filIN)
        Fnu_spec            spectra (Default: None - via filIN)
        extrapoff           set zeros for uncovered wave grid (Default: True)
        verbose             keep tmp files (Default: False)
        ------ OUTPUT ------
        ds                  output dataset
          wcen              center wavelength
          Fnu_filt          synthetic photometry
          smat              standard deviation matrices
        '''
        ds = type('', (), {})()

        ## Convert all format phot names to list
        filt = listize(filt)

        ## Input is a FITS file
        if self.filIN is not None:
            w_spec = self.wvl
            Fnu_spec = self.im

        w_spec = np.array(w_spec)
        Fnu_spec = np.array(Fnu_spec)
        if len(Fnu_spec.shape)==1:
            Ndim = 1
            Fnu_spec = Fnu_spec[:,np.newaxis,np.newaxis]
        else:
            Ndim = 3

        ## Super pixels
        Nw, Ny, Nx = Fnu_spec.shape
        Nxs = math.ceil(Nx/xscale)
        spec_supx = np.zeros((Nw,Ny,Nxs))
        for xs in range(Nxs):
            x = sup2pix(xs, xscale, Npix=Nx, origin=0)
            spec_supx[:,:,xs] += np.nanmean(Fnu_spec[:,:,x[0]:x[-1]+1],axis=2)
        Nys = math.ceil(Ny/yscale)
        sup_spec = np.zeros((Nw,Nys,Nxs))
        for ys in range(Nys):
            y = sup2pix(ys, yscale, Npix=Ny, origin=0)
            sup_spec[:,ys,:] += np.nanmean(spec_supx[:,y[0]:y[-1]+1,:],axis=1)

        ## Do not extrapolate the wave grid that is not covered by input spectra
        ##-----------------------------------------------------------------------
        if extrapoff==True:
            for phot in filt:
                # w_grid = read_ascii(rapyroot+'/lib/filt_'+phot, dtype=float)[:,0]
                # w_grid = ascii.read(rapyroot+'/lib/filt_'+phot+ascext,
                #                     names=['Wave','Spectral Response'])['Wave']
                w_grid = read_hdf5(rapyroot+'/lib/data/filt_'+phot,
                                   'Filter wavelength (microns)')
                # print(w_spec[0], w_grid[0])
                # print(w_spec[-1], w_grid[-1])
                if w_spec[0]>w_grid[0] or w_spec[-1]<w_grid[-1]:
                    warnings.warn('Synthetic photometry of {} can be underestimated' \
                                  'due to uncovered wavelengths'.format(phot))
            ## Insert 2 wvl (0.01 um & w_spec[0]-0.01 um) with 0 value
            wave = np.insert(w_spec, 0, (.01, w_spec[0]-.01))
            flux = np.insert(sup_spec, 0, np.zeros(sup_spec.shape[-1]), axis=0)
            flux = np.insert(flux, 0, np.zeros(sup_spec.shape[-1]), axis=0)
        else:
            wave = w_spec
            flux = sup_spec

        ## Write input.h5
        ##----------------
        fortIN = os.getcwd()+'/synthetic_photometry_input'
        
        write_hdf5(fortIN, 'Filter label', filt)
        write_hdf5(fortIN, 'Wavelength (microns)', wave, append=True)
        write_hdf5(fortIN, 'Flux (x.Hz-1)', flux, append=True)
        write_hdf5(fortIN, '(docalib,dophot)', [1,1], append=True)

        ## Call the Fortran lib
        ##----------------------
        SP.call('synthetic_photometry', shell=True)

        ## Read output.h5
        ##----------------
        fortOUT = os.getcwd()+'/synthetic_photometry_output'

        ds.wcen = read_hdf5(fortOUT, 'Central wavelength (microns)')
        ds.sup_filt = read_hdf5(fortOUT, 'Flux (x.Hz-1)')
        ds.smat = read_hdf5(fortOUT, 'Standard deviation matrix')

        ## Original pixels
        ds.Fnu_filt = np.zeros((Nw,Ny,Nx))
        for x in range(Nx):
            for y in range(Ny):
                xs = pix2sup(x, xscale, origin=0)
                ys = pix2sup(y, yscale, origin=0)
                ds.Fnu_filt[:,y,x] = ds.sup_filt[:,ys,xs]

        ## Convert zeros to NaNs
        ds.sup_filt[ds.sup_filt==0] = np.nan
        ds.Fnu_filt[ds.Fnu_filt==0] = np.nan

        ## Reform outputs
        if Ndim==1:
            ds.Fnu_filt = ds.Fnu_filt[:,0,0]
            ds.sup_filt = ds.sup_filt[:,0,0]
        if len(ds.wcen)==1:
            ds.wcen = ds.wcen[0]
            ds.Fnu_filt = ds.Fnu_filt[0]
            ds.sup_filt = ds.sup_filt[0]
            ds.smat = ds.smat[0][0]
        
        ## Clean temperary h5 files
        ##--------------------------
        if not verbose:
            SP.call('rm -rf '+fortIN+h5ext, shell=True, cwd=os.getcwd())
            SP.call('rm -rf '+fortOUT+h5ext, shell=True, cwd=os.getcwd())

        return ds

    def correct_spec(self, gain=1., offset=0., w_spec=None, Fnu_spec=None,
                     wlim=(None,None), ylim=(None,None), xlim=(None,None),
                     header=None, filOUT=None):
        '''
        Correct spectra
        
        ------ INPUT ------
        gain                scalar or ndarray (Default: 1.)
        offset              scalar or ndarray (Default: 0.)
        w_spec              wavelengths (Default: None - via filIN)
        Fnu_spec            spectra (Default: None - via filIN)
        wlim                wave limits (Default: (None,None))
        ylim                y limits (Default: (None,None))
        xlim                x limits (Default: (None,None))
        filOUT              overwrite fits file (Default: None)
        ------ OUTPUT ------
        new_spec            new_spec = gain * Fnu_spec + offset
        '''
        ## Input is a FITS file
        if self.filIN is not None:
            w_spec = self.wvl
            Fnu_spec = self.im

        w_spec = np.array(w_spec)
        Fnu_spec = np.array(Fnu_spec)
        if len(Fnu_spec.shape)==1:
            Ndim = 1
            Fnu_spec = Fnu_spec[:,np.newaxis,np.newaxis]
        else:
            Ndim = 3
        Nw, Ny, Nx = Fnu_spec.shape

        if np.isscalar(gain):
            a = np.array(gain)
            a = np.array(a[np.newaxis,np.newaxis])
            a = np.repeat(a[:,:], Ny, axis=0)
            a = np.repeat(a[:,:], Nx, axis=1)
        else:
            a = gain
        if np.isscalar(offset):
            b = np.array(offset)
            b = np.array(b[np.newaxis,np.newaxis])
            b = np.repeat(b[:,:], Ny, axis=0)
            b = np.repeat(b[:,:], Nx, axis=1)
        else:
            b = offset

        ## Truncate wavelengths
        if wlim[0] is None:
            wmin = w_spec[0]
        else:
            wmin = wlim[0]
        if wlim[1] is None:
            wmax = w_spec[-1]
        else:
            wmax = wlim[1]

        ## Crop map if 3D
        xmin, xmax = xlim
        ymin, ymax = ylim

        ## Modify spectra
        new_spec = np.copy(Fnu_spec)
        for k, lam in enumerate(w_spec):
            if lam>=wmin and lam<=wmax:
                new_spec[k,ymin:ymax,xmin:xmax] = \
                    a[ymin:ymax,xmin:xmax] * Fnu_spec[k,ymin:ymax,xmin:xmax] \
                    + b[ymin:ymax,xmin:xmax]

        ## Reform outputs
        if Ndim==1:
            new_spec = new_spec[:,0,0]
                    
        if filOUT is not None:
            if header is None:
                header = self.hdr
            write_fits(filOUT, header, new_spec, wave=w_spec)
        
        return new_spec

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Calibration

    intercalib:
        read_filter, synthetic_photometry, correct_spec
    photometry_profile

"""

import os
import math
import numpy as np
# from astropy.io import ascii
from matplotlib.ticker import ScalarFormatter, NullFormatter
import subprocess as SP
import warnings
DEVNULL = open(os.devnull, 'w')

## Local
from arrays import listize, pix2sup, sup2pix
from inout import (ascext, fitsext, h5ext,
                   read_fits, write_fits,
                   read_hdf5, write_hdf5#, read_ascii
)
from astrom import fixwcs
from plots import pplot

## Path of current file
croot = os.path.dirname(os.path.abspath(__file__))

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
            w_grid = read_hdf5(croot+'/lib/data/filt_'+phot,
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
                # w_grid = read_ascii(croot+'/lib/filt_'+phot, dtype=float)[:,0]
                # w_grid = ascii.read(croot+'/lib/filt_'+phot+ascext,
                #                     names=['Wave','Spectral Response'])['Wave']
                w_grid = read_hdf5(croot+'/lib/data/filt_'+phot,
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

"""
class spec2phot(intercalib):
    '''
    Intercalibration between spectrometry and photometry (REF)

    --- INPUT ---
    filIN       to convolve
    filREF      convolution ref
    phot        photometry name (once a phot)
    filKER      convolution kernel(s) (Default: None)
    --- OUTPUT ---
    '''
    def __init__(self, filIN, filREF, phot, filKER=None, saveKER=None, \
        uncIN=None, Nmc=0, filOUT=None):
        super().__init__(filIN)
        self.phot = phot

        if self.wave is not None: # filIN is spec
            ## Convolve filIN (spec)
            if filKER is not None:
                conv = iconvolve(filIN, filKER, saveKER, \
                    uncIN, filOUT=filPRO)
            else:
                filPRO = filIN # filPRO is spec
            
            ## Reprojection to spec (filIN)
            pro = imontage(filREF, filPRO)
            F_phot = pro.reproject(filOUT=filOUT)

        else: # filIN is phot
            ## Reset header (should be spec)
            self.hdr = read_fits(filREF).header
            self.im = read_fits(filREF).data
            self.wvl = read_fits(filREF).wave
            
            ## Convolve filIN (phot)
            if filKER is not None:
                conv = iconvolve(filIN, filKER, saveKER, \
                    uncIN, filOUT=filPRO)
            else:
                filPRO = filIN # filPRO is phot
            
            ## Reprojection to spec (filREF)
            pro = imontage(filPRO, filREF)
            F_phot = pro.reproject(filOUT=filOUT)

        ## Synthetic photometry
        wcen, Fsyn, Fsig = self.synthetic_photometry((phot))
        self.wcen = wcen[0]
        self.Fsyn = Fsyn[0]
        self.Fsig = Fsig[0][0]

        self.gain = F_phot / self.Fsyn

    def calib_gain(self):
        return self.gain

    def image(self):
        return self.Fsyn
    
    def write_image(self, filSYN):
        comment = "Synthetic photometry with " + self.phot
        write_fits(filSYN, self.hdr, self.Fsyn, self.wvl, COMMENT=comment)

class phot2phot(intercalib):
    '''
    Intercalibration between two photometry
    '''
    def __init__(self, filIN, filREF, filKER=None, saveKER=None, \
        uncIN=None, Nmc=0, filOUT=None):

        ## Convolution (optional)
        if filKER is not None:
            conv = iconvolve(filIN, filKER, saveKER, \
                uncIN, filOUT=filPRO)
        else:
            filPRO = filIN

        ## Reprojection config
        pro = imontage(filPRO, filREF)
        self.im = pro.reproject(filOUT=filOUT)

    def image(self):
        return self.im
"""

def photometry_profile(datdir=None, *photometry):
    '''
    ------ INPUT ------
    datdir              profile data path (Default: ./lib/)
    photometry          photometry
    ------ OUTPUT ------
    '''
    ## Read data
    ##-----------
    lam = []
    val = []
    for phot in photometry:
        if datdir is None:
            # datdir = croot+'/lib/'
            datdir = croot+'/lib/data/'
        # dat = read_ascii(datdir+'filt_'+phot, dtype=float)
        # lam.append(dat[:,0])
        # val.append(dat[:,1])

        # dat = ascii.read(datdir+'filt_'+phot+ascext,
        #                  names=['Wave','Spectral Response'])
        # lam.append(dat['Wave'])
        # val.append(dat['Spectral Response'])

        lam.append(read_hdf5(datdir+'filt_'+phot,
            'Filter wavelength (microns)'))
        val.append(read_hdf5(datdir+'filt_'+phot,
            'Filter transmission'))
    # lam = np.array(lam)
    # val = np.array(val)

    ## Plotting setting
    ##------------------
    p = pplot(xlim=(1.9, 40.), ylim=(.01, 1.01),
              xlog=1, ylog=0,
              xlabel=r'$Wavelength,\,\,\lambda\,\,[\mu m]$',
              ylabel='Response',
              # ylabel='Spectral response\n(electrons / photon)',
              legend='upper left', legendalpha=.1,
              figsize=(12,3), title=None, clib='tableau')
    for i,w in enumerate(lam):
        p.add_plot(w, val[i], lw=1.8, label=photometry[i])

    p.set_fig(left=.05, bottom=.2, right=.99, top=.99)

    # sizeXL = 50
    # p.set_font(xticksize=sizeXL, yticksize=sizeXL, \
    #     axesize=sizeXL, legendsize=sizeXL)

    ## vlines (e.g. band markers)
    ##----------------------------
    greylines = []
    pinklines = []
    greylines.extend([2.3567863, 5.1532226]) # AKARI/IRC
    greylines.extend([5.242817, 7.597705]) # Spitzer/IRS-SL2
    pinklines.extend([7.3675313, 8.66892]) # Spitzer/IRS-SL3
    greylines.extend([7.5337057, 14.736635]) # Spitzer/IRS-SL1
    greylines.extend([14.266611, 21.051888]) # Spitzer/IRS-LL2
    pinklines.extend([19.483675, 21.50092]) # Spitzer/IRS-LL3
    greylines.extend([20.555237, 38.41488]) # Spitzer/IRS-LL1
    for gl in greylines:
        p.ax.axvline(gl, linestyle='dotted', color='grey')
    for pl in pinklines:
        p.ax.axvline(pl, linestyle='dotted', color='pink')

    ## tick setting
    ##-------------------- x --------------------------
    xtic = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 40]
    xtic_min = np.arange(2., 41., 1.)
    p.ax.set_xticks(xtic, minor=False) # major
    p.ax.set_xticks(xtic_min, minor=True) # minor
    # ScalarFormatter().set_scientific(False)
    p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major
    p.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
    # p.ax.minorticks_off()
    ##--------------------- y --------------------------
    ytic = np.arange(0, 1.01, .2)
    ytic_min = np.arange(0, 1., .1)
    p.ax.set_yticks(ytic, minor=False) # major
    p.ax.set_yticks(ytic_min, minor=True) # minor
    # ScalarFormatter().set_scientific(False)
    p.ax.yaxis.set_major_formatter(ScalarFormatter()) # major
    p.ax.yaxis.set_minor_formatter(NullFormatter()) # minor
    # p.ax.minorticks_off()
    ##-----------------------------------------------

    return p

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

    pass

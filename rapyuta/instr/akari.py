#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

AKARI data processing and analyzing

    cupid(IM.improve):
        spec_build, sav_build,
        header, image, wave

"""

import os
import math
import numpy as np
from scipy.io import readsav
from astropy.io import ascii
import subprocess as SP
import warnings
DEVNULL = open(os.devnull, 'w')

## Local
import rapyuta.utbox as UT
import rapyuta.inout as IO
import rapyuta.impro as IM
from inout import fitsext, h5ext, savext
from .utils import *


##------------------------------------------------
##
##     AKARI/IRC slit spectroscopic data cube
##
##------------------------------------------------

class cupid(IM.improve):
    '''
    
    AKARI/IRC (slit-)spectroscopy data cube builder
    via IRC pipeline extracted spectra or SAV file

    ------ INPUT ------
    filOUT              output FITS file
    ircdir              path of IRC dataset
    obsid               observation id
    slit                slit name ('Nh'/'Ns')
    spec                spec disperser ('NG'/'NP')
    imref               reference image (IRC N3 long exposure frame; 2MASS corrected; 90 deg rot)
    ------ OUTPUT ------
    '''
    def __init__(self, ircdir=None, obsid=None,
                 slit=None, spec=None, imref=None, verbose=False):
        self.path = ircdir + obsid + '/irc_specred_out_' + slit+'/'
        filref = self.path + imref
        super().__init__(filref) # use N3 header
        
        self.filsav = self.path + obsid + '.N3_' + spec + '.IRC_SPECRED_OUT'
        self.table = readsav(self.filsav+savext, python_dict=True)['source_table']

        ## Slit width will be corrected during reprojection
        if slit=='Ns':
            self.slit_width = 5 # 5" (Ns)
            # self.slit_width = 3 # 5"/1.446" = 3.458 pix (Ns)
        elif slit=='Nh':
            self.slit_width = 3 # 3" (Nh)
            # self.slit_width = 2 # 3"/1.446" = 2.075 pix (Nh)

        if verbose==True:
            print('\n----------')
            print('Slit extracted from ')
            print('obs_id: {} \n slit: {}'.format(obsid, slit))
            print('----------\n')
            self.verbose = verbose

        self.obsid = obsid
        self.slit = slit
        self.spec = spec
        self.imref = imref
        
    def spec_build(self, filOUT=None, filRAW=None, dist=None,
                   Nx=None, Ny=32, Nsub=1, pixscale=None,
                   wmin=None, wmax=None, tmpdir=None, fiLOG=None,
                   sig_pt=0, fill_pt='med', supix=False, swarp=False):
        '''
        Build the spectral cube/slit from spectra extracted by IDL pipeline
        (see IRC_SPEC_TOOL, plot_spec_with_image)

        ------ INPUT ------
        Nx                  number of (identical) pixels to fit slit width
                              Default: None (5 for Ns and 3 for Nh when pixscale=1)
        Ny                  number of pixels in spatial direction (Max=32)
                              Y axis in N3 frame (or X axis in focal plane arrays)
        Nsub                number of subslits
                              Default: 1 (exact division of Ny)
        pixscale            pixel size of final grid (Default: None)
        wmin,wmax           truncate wavelengths
        sig_pt              pointing accuracy in arcsec (Default: 0)
        supix               regrouped super pixel of size (xscale,yscale) (Default: False)
        fill_pt             fill value of no data regions after shift
                              'med': axis median (default)
                              'avg': axis average
                              'near': nearest non-NaN value on the same axis
                              float: constant
        swarp               use SWarp to perform position shifts
                              Default: False (not support supix)
        fiLOG               build info (Default: None)
        '''
        if Nx is None:
            Nx = self.slit_width
        ref_x = self.table['image_y'][0] # slit ref x
        ref_y = 512 - self.table['image_x'][0] # slit ref y

        ## Get slit coord from 2MASS corrected N3 frame
        ## Do NOT touch self.im (N3 frame, 2D) before this step
        self.crop(sizpix=(1, Ny), cenpix=(ref_x, ref_y))
        # self.hdr['CTYPE3'] = 'WAVE-TAB'
        self.hdr['CUNIT1'] = 'deg'
        self.hdr['CUNIT2'] = 'deg'
        self.hdr['BUNIT'] = 'MJy/sr'
        self.hdr['EQUINOX'] = 2000.0

        ## Read spec - Ny/Nsub should be integer, or there will be a shi(f)t
        yscale = math.ceil(Ny/Nsub)
        spec_arr = []
        for j in range(Ny):
            ## ATTENTION: the kw 'space_shift' in IRC pipeline follows the
            ## focal plane array coordinates, whose x axis corresponds to
            ## 'image_x'. That means if space_shift>0, x increases.
            ## When we work in N3 frame coordinates, which rotates 90 deg,
            ## if space_shift>0, y decreases.
            ispec = Nsub - 1 - math.floor(j / yscale)
            # ispec = math.floor(j / yscale) # inverse, see tests/test_build_slit
            readspec = ascii.read(self.path+'spec'+str(ispec)+'.spc')
            subslit = []
            for k in readspec.keys():
                subslit.append(readspec[k])
            subslit = np.array(subslit)
            
            spec_arr.append(subslit)
        ## spec_arr.shape = (Ny,4,Nw)
        spec_arr = np.array(spec_arr)
        ## Save spec in wave ascending order
        for i in range(spec_arr.shape[1]):
            for j in range(Ny):
                spec_arr[j,i,:] = spec_arr[j,i,::-1]
        wave = spec_arr[0,0,:]
        if wmin is None:
            wmin = wave[0]
        if wmax is None:
            wmax = wave[-1]
        iwi = LA.closest(wave, wmin, side='right')
        iws = LA.closest(wave, wmax, side='left')+1
        wave = wave[iwi:iws]
        Nw = len(wave)
        
        ## Build cube array
        cube = np.empty([Nw,Ny,1])
        unc = np.empty([Nw,Ny,1]) # Symmetric unc
        unc_N = np.empty([Nw,Ny,1]) # Asymmetric negtive
        unc_P = np.empty([Nw,Ny,1]) # Asymmetric positive
        for k in range(Nw):
            for j in range(Ny):
                cube[k][j] = spec_arr[j,1,k+iwi]
                unc[k][j] = (spec_arr[j,3,k+iwi]-spec_arr[j,2,k+iwi])/2
                unc_N[k][j] = (spec_arr[j,1,k+iwi]-spec_arr[j,2,k+iwi])
                unc_P[k][j] = (spec_arr[j,3,k+iwi]-spec_arr[j,1,k+iwi])

        ## Update self variables (for next steps)
        self.reinit(header=self.hdr, image=cube, wave=wave,
                    wmod=self.wmod, verbose=self.verbose)

        if filRAW is not None:
            comment = "Assembled AKARI/IRC slit spec cube. "
            IO.write_fits(filRAW, self.hdr, cube, self.wvl,
                          COMMENT=comment)
            comment = "Assembled AKARI/IRC slit spec uncertainty cube. "
            IO.write_fits(filRAW+'_unc', self.hdr, unc, self.wvl,
                          COMMENT=comment)
            comment = "Assembled AKARI/IRC slit spec uncertainty (N) cube. "
            IO.write_fits(filRAW+'_unc_N', self.hdr, unc_N, self.wvl,
                          COMMENT=comment)
            comment = "Assembled AKARI/IRC slit spec uncertainty (P) cube. "
            IO.write_fits(filRAW+'_unc_P', self.hdr, unc_P, self.wvl,
                          COMMENT=comment)

        ## Uncertainty propagation
        if dist=='norm':
            self.rand_norm(unc=unc)
        elif dist=='splitnorm':
            self.rand_splitnorm(unc=[unc_N,unc_P])
        
        ## Rescale pixel size
        if pixscale is not None:
            self.rebin(pixscale=pixscale)
            self.im = np.delete(self.im, [i+1 for i in range(self.Nx-1)], axis=2) # homo x
        ## Broaden slit width
        self.im = np.repeat(self.im, Nx, axis=2)
        self.hdr['CRPIX1'] = (Nx+1)/2
        ## Extrapolate discrepancy of rebinned slit length (up to pixscale*Nsub)
        newyscale =  math.ceil(self.Ny/Nsub)
        dNy = Nsub * newyscale - self.Ny
        # ystart = (Nsub - 1) * newyscale
        if newyscale==1:
            ystart = -1
        else:
            ystart = 1 - newyscale
        val = np.nanmean(self.im[:,ystart:,:], axis=1)
        arr = np.repeat(val[:,np.newaxis,:], dNy, axis=1)
        self.im = np.append(self.im, arr, axis=1)
        self.hdr['CRPIX2'] += (Nx+1)/2
        ## Update self variables (for next steps)
        self.reinit(header=self.hdr, image=self.im, wave=self.wvl,
                    wmod=self.wmod, verbose=self.verbose)
        # print(self.Nx, self.hdr['CRPIX1'])
        # print(self.Ny, self.hdr['CRPIX2'])
        self.xscale = self.Nx
        self.yscale =  math.ceil(self.Ny/Nsub) # A priori Ny/Nsub is integer...

        ## Add pointing unc
        if supix:
            self.rand_pointing(sig_pt, fill=fill_pt, tmpdir=tmpdir,
                               xscale=self.xscale, yscale=self.yscale, swarp=False)
        else:
            self.rand_pointing(sig_pt, fill=fill_pt, tmpdir=tmpdir,
                               xscale=1, yscale=1, swarp=swarp)
                
        ## Update self variables (for next steps)
        self.reinit(header=self.hdr, image=self.im, wave=self.wvl,
                    wmod=self.wmod, verbose=self.verbose)

        if filOUT is not None:
            comment = "<cupid> Assembled AKARI/IRC slit spectroscopy cube. "
            IO.write_fits(filOUT, self.hdr, self.im, self.wvl,
                          COMMENT=comment)
            
        if fiLOG is not None:
            IO.write_hdf5(fiLOG, 'Observation ID', [self.obsid])
            IO.write_hdf5(fiLOG, 'NIR Slit', [self.slit], append=True)
            IO.write_hdf5(fiLOG, 'Spectral disperser', [self.spec], append=True)
            IO.write_hdf5(fiLOG, 'N3 image',[self.imref], append=True)
            IO.write_hdf5(fiLOG, 'Pointing accuracy', [sig_pt], append=True)
            IO.write_hdf5(fiLOG, 'Spectral sampling size', [self.Nw], append=True)
            IO.write_hdf5(fiLOG, 'Slit length', [self.Ny], append=True)
            IO.write_hdf5(fiLOG, 'Slit width', [self.Nx], append=True)
            IO.write_hdf5(fiLOG, 'Subslit number', [Nsub], append=True)
            IO.write_hdf5(fiLOG, 'Pixel size', [pixscale], append=True)
            IO.write_hdf5(fiLOG, 'Super pixel size',
                          [self.xscale, self.yscale], append=True)
        
        return self.im

    def sav_build(self):
        '''
        Alternative extraction from SAV file
        Including wave calib, ?, etc. 
        (see IRC_SPEC_TOOL, plot_spec_with_image)
        '''
        print('NOT AVAILABLE')
        exit()
        
        filsav = self.filsav
        table = self.table
        ## Read SAV file
        image = readsav(filsav+savext, python_dict=True)['specimage_n_wc']
        image = image[::-1] # -> ascending order
        noise = readsav(filsav+savext, python_dict=True)['noisemap_n']
        noise = noise[::-1]
        wave = readsav(filsav+savext, python_dict=True)['wave_array']
        wave = wave[::-1] # -> ascending order
        Nw = image.shape[0] # num of wave
        Ny = image.shape[1] # slit length
        ref_x = table['image_y'][0] # slit ref x
        ref_y = 512-table['image_x'][0] # slit ref y
        spec_y = table['spec_y'][0] # ref pts of wavelength
        
        d_wave_offset_pix = -(spec_y-round(spec_y[0])) # Wave shift
        warr = np.arange(Nw)
        wave_shift = np.interp(warr+d_wave_offset_pix, warr, wave)
        
        for k in range(Nw):
            for j in range(Ny):
                for i in range(Nx):
                    cube[k][j][i] = image[k][j]
                    unc[k][j][i] = noise[k][j]

    def header(self):
        return self.hdr
    
    def image(self):
        return self.im

    def wave(self):
        return self.wvl

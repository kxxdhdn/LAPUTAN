#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Utilities for impro

    improve:
        reinit, BGunc, 
        rand_norm, rand_splitnorm, rand_pointing, 
        slice, slice_inv_sq, crop, rebin, groupixel
        smooth, artifact, mask

"""

import os
import math
import numpy as np
from scipy.interpolate import interp1d
from astropy import wcs
import warnings

## Local
import rapyuta.utbox as UT
import rapyuta.inout as IO
import rapyuta.latte as LA
import rapyuta.maths as MA
from rapyuta.inout import fitsext

##------------------------------------------------
##
##           IMage PROcessing VEssel
##
##------------------------------------------------

class improve:
    '''
    IMage PROcessing VEssel

    '''
    def __init__(self, filIN=None, header=None, images=None,
                 wave=None, wmod=0, whdr=None,
                 filUNC=None, verbose=False, filext=fitsext,
                 instr=None, instr_auto=True):
        '''
        filUNC - if not None, should be full name with ".fits"!

        self: filIN, header, images, wave, wmod, whdr,
              filUNC, verbose, filext, instr, instr_auto
              ## [derived]
              unc, wcs, cdelt, pc, cd, Ndim, Nx, Ny, Nw, 
        '''
        ## INPUTS
        self.filIN = filIN
        self.header = header
        self.images = images
        self.wave = wave
        self.wmod = wmod
        self.whdr = whdr
        self.filUNC = filUNC
        self.verbose = verbose
        self.filext = filext
        self.instr = instr
        self.instr_auto = self.instr_auto

        ## Default initialization
        self.unc = None
        self.wcs = None
        self.cdelt, self.pc, self.cd = None, None, None
        self.Ndim = 2
        self.Nx, self.Ny, self.Nw = 1, 1, 0

        ## filIN update: header, images, wave, whdr
        if filIN is not None:
            ds = IO.read_fits(filIN, wmod=wmod,
                              instr=instr, instr_auto=instr_auto,
                              uncfiles=filUNC, filext=filext)
            self.header = ds.header
            self.images = ds.data
            self.wave = ds.wave
            self.whdr = ds.whdr
            self.unc = ds.unc
            self.wcs = ds.wcs
            self.cdelt, self.pc, self.cd = ds.cdelt, ds.pc, ds.cd
            self.Ndim = ds.Ndim
            self.Nx, self.Ny, self.Nw = ds.Nx, ds.Ny, ds.Nw
        else:
            ## self.unc
            filUNC = LA.listize(filUNC)
            unc = []
            for uncf in filUNC:
                if Path(uncf).exists():
                    self.unc = IO.read_fits(uncf, filext='').data
                    # with fits.open(uncf) as hdul:
                    #     unc.append(hdul[0].data)
            if len(unc)==0:
                ## no unc
                self.unc = None
            elif len(unc)==1:
                ## symmetric unc
                self.unc = unc[0]
            else:
                ## asymmetric unc
                self.unc = unc
            
            ## self.wcs
            if header is not None:
                try:
                    self.wcs = wcs.WCS(header)
                except:
                    self.wcs = patch_wcs_3D(header=header).wcs
            
            ## self: cdelt, pc, cd
            if self.wcs is not None:
                pcdelt = get_pc(wcs=w.wcs)
                self.cdelt = pcdelt.cdelt
                self.pc = pcdelt.pc
                self.cd = pcdelt.cd

            ## self: Ndim, Nx, Ny, Nw
            if self.images is not None:
                self.Ndim = self.images.ndim
                if self.Ndim==3:
                    self.Nw, self.Ny, self.Nx = self.images.shape
            
                    ## Nw=1 patch
                    if self.images.shape[0]==1:
                        self.Ndim = 2
                elif self.Ndim==2:
                    self.Ny, self.Nx = self.images.shape
                    self.Nw = None

        if self.verbose:
            print(f'<improve> file: {filIN}')
            print(f'Raw size (pix): {self.Nx} * {self.Ny}')

    def reinit(self, **kwargs):
        '''
               instr=None, instr_auto=True

        Update init variables

        filUNC - if not None, should be full name with ".fits"!
        '''
        ## flags for updating derived instances (self)
        flag_wcs = False
        flag_img = False
        flag_unc = False
        
        ## INPUTS
        if 'filIN' in kwargs:
            self.filIN = kwargs['filIN']
            flag_wcs = True
        if 'header' in kwargs:
            self.header = kwargs['header']
            flag_wcs = True
        if 'images' in kwargs:
            self.images = kwargs['images']
            flag_img = True
        if 'wave' in kwargs:
            self.wave = kwargs['wave']
        if 'wmod' in kwargs:
            self.wmod = kwargs['wmod']
        if 'whdr' in kwargs:
            self.whdr = kwargs['whdr']
        if 'filUNC' in kwargs:
            self.filUNC = kwargs['filUNC']
            flag_unc = True
        if 'verbose' in kwargs:
            self.verbose = kwargs['verbose']
        if 'filext' in kwargs:
            self.filext = kwargs['filext']
        if 'instr' in kwargs:
            self.instr = kwargs['instr']
        if 'instr_auto' in kwargs:
            self.instr_auto = kwargs['instr_auto']

        ## filIN update: all
        if ('filIN' in kwargs) and (self.filIN is not None):
            ds = IO.read_fits(self.filIN, wmod=self.wmod,
                              instr=self.instr, instr_auto=self.instr_auto,
                              uncfiles=self.filUNC, filext=self.filext)
            self.header = ds.header
            self.images = ds.data
            self.wave = ds.wave
            self.whdr = ds.whdr
            self.unc = ds.unc
            self.wcs = ds.wcs
            self.cdelt, self.pc, self.cd = ds.cdelt, ds.pc, ds.cd
            self.Ndim = ds.Ndim
            self.Nx, self.Ny, self.Nw = ds.Nx, ds.Ny, ds.Nw
        else:
            ## filUNC update: self.unc
            if flag_unc:
                ## self.unc
                filUNC = LA.listize(self.filUNC)
                unc = []
                for uncf in filUNC:
                    if Path(uncf).exists():
                        self.unc = IO.read_fits(uncf, filext='').data
                        # with fits.open(uncf) as hdul:
                        #     unc.append(hdul[0].data)
                if len(unc)==0:
                    ## no unc
                    self.unc = None
                elif len(unc)==1:
                    ## symmetric unc
                    self.unc = unc[0]
                else:
                    ## asymmetric unc
                    self.unc = unc
            
            ## header update: self.wcs
            if flag_wcs and (self.header is not None):
                try:
                    self.wcs = wcs.WCS(self.header)
                except:
                    self.wcs = patch_wcs_3D(header=self.header).wcs

                ## wcs update: self: cdelt, pc, cd
                if self.wcs is not None:
                    pcdelt = get_pc(wcs=w.wcs)
                    self.cdelt = pcdelt.cdelt
                    self.pc = pcdelt.pc
                    self.cd = pcdelt.cd

            ## images update: self: Ndim, Nx, Ny, Nw
            if flag_img and (self.images is not None):
                self.Ndim = self.images.ndim
                if self.Ndim==3:
                    self.Nw, self.Ny, self.Nx = self.images.shape
            
                    ## Nw=1 patch
                    if self.images.shape[0]==1:
                        self.Ndim = 2
                elif self.Ndim==2:
                    self.Ny, self.Nx = self.images.shape
                    self.Nw = None

        if self.verbose:
            print(f'<improve> file: {filIN}')
            print(f'Raw size (pix): {self.Nx} * {self.Ny}')

    def BGunc(self, filOUT=None, filWGT=None, wfac=1.,
              BG_images=None, BG_weight=None, fill_zeros=np.nan
              filext=fitsext):
        '''
        Estimate uncertainties from the background map
        So made error map is uniform/weighted

        Updated: self.unc
        ------ INPUT ------
        filOUT              output uncertainty map (FITS)
        filWGT              input weight map (FITS)
                              if not None, should be full name with ".fits"!
        wfac                multiplication factor for filWGT (Default: 1)
        BG_images           background image array used to generate unc map
        BG_weight           background weight array
        fill_zeros          value used to replace zero value (Default: NaN)
        filext              suffix of filOUT (Default: ".fits")
        ------ OUTPUT ------
        self.unc                 estimated unc map
        '''
        # if self.unc is None:
        if BG_images is not None:
            images = BG_images
            Ny, Nx = BG_images.shape
        else:
            images = self.images
            Ny, Nx = self.Ny, self.Nx
        Nw = self.Nw

        ## sigma: std dev of (weighted) flux distribution of bg region
        if BG_weight is not None:
            if self.Ndim==3:
                sigma = np.nanstd(images * BG_weight, axis=(1,2))
            elif self.Ndim==2:
                sigma = np.nanstd(images * BG_weight)
        else:
            if self.Ndim==3:
                sigma = np.nanstd(images, axis=(1,2))
            elif self.Ndim==2:
                sigma = np.nanstd(images)

        ## wgt: weight map
        if filWGT is not None:
            wgt = IO.read_fits(filWGT, filext='').data * wfac
        else:
            wgt = np.ones(images.shape) * wfac

        ## unc: weighted rms = root of var/wgt
        if self.Ndim==3:
            unc = []
            for w in range(Nw):
                unc.append(np.sqrt(1./wgt[w,:,:]) * sigma(w))
            unc = np.array(unc)
        elif self.Ndim==2:
            unc = np.sqrt(1./wgt) * sigma

        ## Replace zero values
        unc[unc==0] = fill_zeros

        self.unc = unc
        
        if filOUT is not None:
            IO.write_fits(filOUT, header=self.header, data=self.unc,
                          wave=self.wave, wmod=self.wmod, filext=filext)
            
        return self.unc

    def rand_norm(self, mu=0., sigma=1.):
        '''
        Add random N(0,1) errors to images

        Updated: self.images
        ------ INPUT ------
        mu                  mean
        sigma               standard deviation
        ------ OUTPUT ------
        self.images         error added images
        '''
        if self.unc is not None:
            ## unc should have the same dimension with images
            theta = np.random.normal(mu, sigma, self.images.shape)
            self.images += theta * self.unc
        else:
            if self.verbose:
                warnings.warn('Uncertainty data unavailable. Nothing changed.')

        return self.images

    def rand_splitnorm(self, mu=0., sigma=1.):
        '''
        Add random SN(0,lam,lam*tau) errors to images

        Updated: self.images
        ------ INPUT ------
        mu                  mean
        sigma               standard deviation
        ------ OUTPUT ------
        self.images         error added images
        '''
        if len(self.unc)==2:
            unc = self.unc
            ## unc[i] should have the same dimension with images
            tau = unc[1]/unc[0]
            peak = 1/(1+tau)
            theta = np.random.normal(mu, sigma, self.images.shape) # ~N(0,1)
            flag = np.random.random(self.images.shape) # ~U(0,1)
            if self.Ndim==2:
                for x in range(self.Nx):
                    for y in range(self.Ny):
                        if flag[y,x]<peak[y,x]:
                            self.images[y,x] += -abs(theta[y,x]) * unc[0][y,x]
                        else:
                            self.images[y,x] += abs(theta[y,x]) * unc[1][y,x]
            elif self.Ndim==3:
                for x in range(self.Nx):
                    for y in range(self.Ny):
                        for k in range(self.Nw):
                            if flag[k,y,x]<peak[k,y,x]:
                                self.images[k,y,x] += -abs(
                                    theta[k,y,x]) * unc[0][k,y,x]
                            else:
                                self.images[k,y,x] += abs(
                                    theta[k,y,x]) * unc[1][k,y,x]
        else:
            if self.verbose:
                warnings.warn('Uncertainty data unavailable. Nothing changed.')

        return self.images

    def rand_pointing(self, accrand=0, filltype='near',
                      xscale=1, yscale=1):
        '''
        Add flux errors converted from pointing errors in coordinates to images (without changing header/WCS)

        Updated: self.images
        ------ INPUT ------
        accrand             pointing accuracy (in arcsec)
        filltype            filling value of no data regions after image shifting
                              'near': the nearest non-NaN value on the same axis (default)
                              'med': axis median
                              'avg': axis average
                              float: constant
        xscale,yscale       regrouped super pixel size
        swarp               use SWarp to perform position shifts
                              Default: False (not support supix)
        ------ OUTPUT ------
        '''
        ## Check inputs
        ##--------------
        if accrand<0:
            UT.strike('improve.rand_pointing', 'negtive pointing accuracy.',
                      cat='InputError')

        ## Effects of pointing shift on header: d_ro, d_phi (WCS) -> dx, dy (pix)
        ##------------------------------------------------------------------------
        # WCS increments: d_ro * cos(d_phi) and d_ro * sin(d_phi)
        accrand /= 3600.
        d_ro = abs(np.random.normal(0., accrand)) # N(0,accrand)
        d_phi = np.random.random() * 2. * np.pi # U(0,2*pi)
        # d_ro, d_phi = 0.0002, 4.5
        # print('d_ro, d_phi = ', d_ro, d_phi)
        
        ## Old header/WCS
        oldheader = self.header
        oldwcs = self.wcs
        Nx, Ny, Nw = self.Nx, self.Ny, self.Nw
        ## New header/WCS
        newheader = oldheader.copy()
        newheader['CRVAL1'] += d_ro * np.cos(d_phi)
        newheader['CRVAL2'] += d_ro * np.sin(d_phi)
        try:
            newcs = wcs.WCS(newheader)
        except:
            newcs = IO.patch_wcs_3D(header=newheader).wcs
        ## Convert WCS increments to image increments (consistent with DS9)
        newheader['CRPIX1'], newheader['CRPIX2'] = \
            oldwcs.all_world2pix(newheader['CRVAL1'], newheader['CRVAL2'], 1)
        dx = newheader['CRPIX1'] - oldheader['CRPIX1']
        dy = newheader['CRPIX2'] - oldheader['CRPIX2']
        ## Check increments at (1,1)
        # print('Near CRPIXn increments: ', dx, dy)
        # val1 = np.array(newcs.all_pix2world(0.5, 0.5, 1))
        # dx, dy = oldwcs.all_world2pix(val1[np.newaxis,:], 1)[0] - 0.5
        # print('Near (1,1) increments: ', dx, dy)

        ## Flux resampling (with rescaled pixels)
        ##----------------------------------------
        oldimage = self.images
        if self.Ndim==3:
            ## X axis rescaling (if not, set xscale=1)
            Nxs = math.ceil(Nx/xscale) # Number of x axis pixel after scaling
            IMGxs = np.zeros((Nw,Ny,Nxs)) # Image after x axis scaling
            frac2 = dx / xscale
            f2 = math.floor(frac2)
            frac1 = 1 - frac2
            for xs in range(Nxs):
                if frac2>=0:
                    x0 = LA.zoom(0, 1/xscale)
                else:
                    x0 = LA.zoom(Nxs-1, 1/xscale)
                ## fi
                if filltype=='med':
                    fill_value = np.nanmedian(self.images,axis=2)
                elif filltype=='avg':
                    fill_value = np.nanmean(self.images,axis=2)
                elif filltype=='near':
                    fill_value = np.nanmean(self.images[:,:,x0[0]:x0[-1]+1],axis=2)
                else:
                    fill_value = fill
                if frac2>=0:
                    if xs>=f2:
                        x1 = LA.zoom(xs-f2, 1/xscale)
                        IMGxs[:,:,xs] += (f2+frac1) * np.nanmean(self.images[:,:,x1[0]:x1[-1]+1],axis=2)
                        if xs>f2:
                            x2 = LA.zoom(xs-f2-1, 1/xscale)
                            IMGxs[:,:,xs] += (frac2-f2) * np.nanmean(self.images[:,:,x2[0]:x2[-1]+1],axis=2)
                        else:
                            IMGxs[:,:,xs] += (frac2-f2) * fill_value
                    else:
                        IMGxs[:,:,xs] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super x = {}'.format(xs))
                else:
                    if xs<=Nxs+f2:
                        x2 = LA.zoom(xs-f2-1, 1/xscale)
                        IMGxs[:,:,xs] += (frac2-f2) * np.nanmean(self.images[:,:,x2[0]:x2[-1]+1],axis=2)
                        if xs<Nxs+f2:
                            x1 = LA.zoom(xs-f2, 1/xscale)
                            IMGxs[:,:,xs] += (f2+frac1) * np.nanmean(self.images[:,:,x1[0]:x1[-1]+1],axis=2)
                        else:
                            IMGxs[:,:,xs] += (f2+frac1) * fill_value
                    else:
                        IMGxs[:,:,xs] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super x = {}'.format(xs))
            
            Nys = math.ceil(Ny/yscale) # Number of y axis pixel after scaling
            newimage = np.zeros((Nw,Nys,Nxs)) # Image after (x and) y axis scaling
            frac2 = dy / yscale
            f2 = math.floor(frac2)
            frac1 = 1 - frac2
            for ys in range(Nys):
                if frac2>=0:
                    y0 = LA.zoom(0, 1/yscale)
                else:
                    y0 = LA.zoom(Nys-1, 1/yscale)
                if filltype=='med':
                    fill_value = np.nanmedian(IMGxs,axis=1)
                elif filltype=='avg':
                    fill_value = np.nanmean(IMGxs,axis=1)
                elif filltype=='near':
                    fill_value = np.nanmean(IMGxs[:,y0[0]:y0[-1]+1,:],axis=1)
                else:
                    fill_value = fill
                if frac2>=0:
                    if ys>=f2:
                        y1 = LA.zoom(ys-f2, 1/yscale)
                        newimage[:,ys,:] += (f2+frac1) * np.nanmean(IMGxs[:,y1[0]:y1[-1]+1,:],axis=1)
                        if ys>f2:
                            y2 = LA.zoom(ys-f2-1, 1/yscale)
                            newimage[:,ys,:] += (frac2-f2) * np.nanmean(IMGxs[:,y2[0]:y2[-1]+1,:],axis=1)
                        else:
                            newimage[:,ys,:] += (frac2-f2) * fill_value
                    else:
                        newimage[:,ys,:] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super y = {}'.format(ys))
                else:
                    if ys<=Nys+f2:
                        y2 = LA.zoom(ys-f2-1, 1/yscale)
                        newimage[:,ys,:] += (frac2-f2) * np.nanmean(IMGxs[:,y2[0]:y2[-1]+1,:],axis=1)
                        if ys<Nys+f2:
                            y1 = LA.zoom(ys-f2, 1/yscale)
                            newimage[:,ys,:] += (f2+frac1) * np.nanmean(IMGxs[:,y1[0]-1:y1[-1],:],axis=1)
                        else:
                            newimage[:,ys,:] += (f2+frac1) * fill_value
                    else:
                        newimage[:,ys,:] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super y = {}'.format(ys))
            
            for x in range(Nx):
                for y in range(Ny):
                    xs = LA.zoom(x, xscale)
                    ys = LA.zoom(y, yscale)
                    self.images[:,y,x] = newimage[:,ys,xs]
            
        elif self.Ndim==2:
            Nxs = math.ceil(Nx/xscale)
            IMGxs = np.zeros((Ny,Nxs))
            frac2 = dx / xscale
            f2 = math.floor(frac2)
            frac1 = 1 - frac2
            for xs in range(Nxs):
                if frac2>=0:
                    x0 = LA.zoom(0, 1/xscale)
                else:
                    x0 = LA.zoom(Nxs-1, 1/xscale)
                if filltype=='med':
                    fill_value = np.nanmedian(self.images,axis=1)
                elif filltype=='avg':
                    fill_value = np.nanmean(self.images,axis=1)
                elif filltype=='near':
                    fill_value = np.nanmean(self.images[:,x0[0]:x0[-1]+1],axis=1)
                else:
                    fill_value = fill
                if frac2>=0:
                    if xs>=f2:
                        x1 = LA.zoom(xs-f2, 1/xscale)
                        IMGxs[:,xs] += (f2+frac1) * np.nanmean(self.images[:,x1[0]:x1[-1]+1],axis=1)
                        if xs>f2:
                            x2 = LA.zoom(xs-f2-1, 1/xscale)
                            IMGxs[:,xs] += (frac2-f2) * np.nanmean(self.images[:,x2[0]:x2[-1]+1],axis=1)
                        else:
                            IMGxs[:,xs] += (frac2-f2) * fill_value
                    else:
                        IMGxs[:,xs] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super x = {}'.format(xs))
                else:
                    if xs<=Nxs+f2:
                        x2 = LA.zoom(xs-f2-1, 1/xscale)
                        IMGxs[:,xs] += (frac2-f2) * np.nanmean(self.images[:,x2[0]:x2[-1]+1],axis=1)
                        if xs<Nxs+f2:
                            x1 = LA.zoom(xs-f2, 1/xscale)
                            IMGxs[:,xs] += (f2+frac1) * np.nanmean(self.images[:,x1[0]:x1[-1]+1],axis=1)
                        else:
                            IMGxs[:,xs] += (f2+frac1) * fill_value
                    else:
                        IMGxs[:,xs] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super x = {}'.format(xs))
            
            Nys = math.ceil(Ny/yscale)
            newimage = np.zeros((Nys,Nxs))
            frac2 = dy / yscale
            f2 = math.floor(frac2)
            frac1 = 1 - frac2
            for ys in range(Nys):
                if frac2>=0:
                    y0 = LA.zoom(0, 1/yscale)
                else:
                    y0 = LA.zoom(Nys-1, 1/yscale)
                if filltype=='med':
                    fill_value = np.nanmedian(IMGxs,axis=0)
                elif filltype=='avg':
                    fill_value = np.nanmean(IMGxs,axis=0)
                elif filltype=='near':
                    fill_value = np.nanmean(IMGxs[y0[0]:y0[-1]+1,:],axis=0)
                else:
                    fill_value = fill
                if frac2>=0:
                    if ys>=f2:
                        y1 = LA.zoom(ys-f2, 1/yscale)
                        newimage[ys,:] += (f2+frac1) * np.nanmean(IMGxs[y1[0]:y1[-1]+1,:],axis=0)
                        if ys>f2:
                            y2 = LA.zoom(ys-f2-1, 1/yscale)
                            newimage[ys,:] += (frac2-f2) * np.nanmean(IMGxs[y2[0]:y2[-1]+1,:],axis=0)
                        else:
                            newimage[ys,:] += (frac2-f2) * fill_value
                    else:
                        newimage[ys,:] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super y = {}'.format(ys))
                else:
                    if ys<=Nys+f2:
                        y2 = LA.zoom(ys-f2-1, 1/yscale)
                        newimage[ys,:] += (frac2-f2) * np.nanmean(IMGxs[y2[0]:y2[-1]+1,:],axis=0)
                        if ys<Nys+f2:
                            y1 = LA.zoom(ys-f2, 1/yscale)
                            newimage[ys,:] += (f2+frac1) * np.nanmean(IMGxs[y1[0]-1:y1[-1],:],axis=0)
                        else:
                            newimage[ys,:] += (f2+frac1) * fill_value
                    else:
                        newimage[ys,:] += fill_value
                        # if self.verbose:
                        #     warnings.warn('Zero appears at super y = {}'.format(ys))
            
            for x in range(Nx):
                for y in range(Ny):
                    xs = LA.zoom(x, xscale)
                    ys = LA.zoom(y, yscale)
                    self.images[y,x] = newimage[ys,xs]

        ## Original NaN mask
        mask_nan = np.isnan(oldimage)
        self.images[mask_nan] = np.nan
        ## Recover new NaN pixels with zeros
        mask_recover = np.logical_and(np.isnan(self.images), ~mask_nan)
        self.images[mask_recover] = 0
            
        return self.images

    def slice(self, filSLC, postfix='', filext=fitsext):
        '''
        Slice a cube

        ------ INPUT ------
        filSLC              filename root of sliced images
        postfix             postfix after filSLC+"_0000", before ".fits"
        filext              suffix of sliced image filename (Default: ".fits")
        ------ OUTPUT ------
        slcnames            list of sliced image filenames
        '''
        ## 3D cube slicing
        slcnames = []
        if self.Ndim==3:
            for k in range(self.Nw):
                ## output filename list
                fname = filSLC+'_'+'0'*(4-len(str(k+1)))+str(k+1)+postfix
                slcnames.append(fname+filext)
                hd2D = patch_wcs_3D(header=self.header).header
                IO.write_fits(fname, header=hd2D, data=self.images[k,:,:],
                              filext=filext)
        elif self.Ndim==2:
            if self.verbose:
                warnings.warn('2D image cannot be sliced. Nothing changed.')

        return slcnames

    def slice_inv_sq(self, filSLC, postfix='', filext=fitsext):
        '''
        Slice a cube with its inversed square

        ------ INPUT ------
        filSLC              filename root of sliced images
        postfix             postfix after filSLC+"_0000", before ".fits"
        filext              suffix of sliced image filename (Default: ".fits")
        ------ OUTPUT ------
        slcnames            list of sliced image filenames
        '''
        ## Inversed square cube slicing
        inv_sq = 1./self.images**2
        
        slcnames = []
        if self.Ndim==3:
            for k in range(self.Nw):
                ## output filename list
                fname = filSLC+'_'+'0'*(4-len(str(k+1)))+str(k+1)+postfix
                slcnames.append(fname+filext)
                hd2D = patch_wcs_3D(header=self.header).header
                IO.write_fits(fname, header=hd2D, data=inv_sq[k,:,:],
                              filext=filext)
        elif self.Ndim==2:
            if self.verbose:
                warnings.warn('2D image cannot be sliced. Nothing changed.')

        return slcnames
    
    def crop(self, filOUT=None, sizpix=None, cenpix=None,
             sizval=None, cenval=None):
        '''
        If pix and val co-exist, pix will be taken.

        ------ INPUT ------
        filOUT              output file
        sizpix              crop size in pix (dx, dy)
        cenpix              crop center in pix (x, y)
        sizval              crop size in deg (dRA, dDEC) -> (dx, dy)
        cenval              crop center in deg (RA, DEC) -> (x, y)
        ------ OUTPUT ------
        self.images         cropped image array
        '''
        oldimage = self.images
        newheader = self.header.copy()
        
        ## Crop center
        ##-------------
        if cenpix is None:
            if cenval is None:
                UT.strike('improve.crop', 'miss crop center', cat='InputError')
            else:
                ## Convert coord
                try:
                    cenpix = np.array(self.w.all_world2pix(cenval[0], cenval[1], 1))
                except wcs.wcs.NoConvergence as e:
                    cenpix = e.best_solution
                    print("Best solution:\n{0}".format(e.best_solution))
                    print("Achieved accuracy:\n{0}".format(e.accuracy))
                    print("Number of iterations:\n{0}".format(e.niter))
        else:
            cenval = self.w.all_pix2world(np.array([cenpix]), 1)[0]
        if not (0<cenpix[0]-0.5<self.Nx and 0<cenpix[1]-0.5<self.Ny):
            UT.strike('improve.crop', 'crop centre overpassed the image edge.',
                      cat='ValueError')

        ## Crop size
        ##-----------
        if sizpix is None:
            if sizval is None:
                UT.strike('improve.crop', 'miss crop size', cat='InputError')
            else:
                ## CDELTn needed (Physical increment at the reference pixel)
                sizpix = np.array(sizval) / abs(self.cdelt)
                sizpix = np.array([math.floor(n) for n in sizpix])
        else:
            sizval = np.array(sizpix) * abs(self.cdelt)

        if self.verbose==True:
            print('----------')
            print("Crop centre (RA, DEC): [{:.8}, {:.8}]".format(*cenval))
            print("Crop size (dRA, dDEC): [{}, {}]\n".format(*sizval))
            print("Crop centre (x, y): [{}, {}]".format(*cenpix))
            print("Crop size (dx, dy): [{}, {}]".format(*sizpix))
            print('----------')
        
        ## Lowerleft origin
        ##------------------
        xmin = math.floor(cenpix[0] - sizpix[0]/2.)
        ymin = math.floor(cenpix[1] - sizpix[1]/2.)
        xmax = xmin + sizpix[0]
        ymax = ymin + sizpix[1]

        if not (xmin>=0 and xmax<=self.Nx and ymin>=0 and ymax<=self.Ny):
            UT.strike('improve.crop', 'crop region overpassed the image edge.',
                      cat='ValueError')

        ## OUTPUTS
        ##---------
        ## New images
        if self.Ndim==3:
            newimage = oldimage[:, ymin:ymax, xmin:xmax]
        elif self.Ndim==2:
            newimage = oldimage[ymin:ymax, xmin:xmax]

        ## Modify header
        ##---------------
        newheader['CRPIX1'] = math.floor(sizpix[0]/2. + 0.5)
        newheader['CRPIX2'] = math.floor(sizpix[1]/2. + 0.5)
        newheader['CRVAL1'] = cenval[0]
        newheader['CRVAL2'] = cenval[1]
        
        self.header = newheader
        self.images = newimage
        
        ## Write cropped image/cube
        if filOUT is not None:
            # comment = "[ICROP]ped at centre: [{:.8}, {:.8}]. ".format(*cenval)
            # comment = "with size [{}, {}] (pix).".format(*sizpix)
            IO.write_fits(filOUT, self.header, self.images, self.wave, self.wmod)

        ## Update self variables
        self.reinit(header=self.header, images=self.images, wave=self.wave,
                    wmod=self.wmod, verbose=self.verbose)

        return self.images

    def rebin(self, pixscale=None, total=False, extrapol=False, filOUT=None):
        '''
        Shrinking (box averaging) or expanding (bilinear interpolation) astro images
        New/old images collimate on zero point.

        [REF] IDL lib frebin/hrebin
        https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/hrebin.pro
        https://github.com/wlandsman/IDLAstro/blob/master/pro/frebin.pro

        ------ INPUT ------
        pixscale            output pixel scale in arcsec/pixel
                              scalar - square pixel
                              tuple - same Ndim with images
        total               Default: False
                              True - sum the non-NaN pixels
                              False - mean
        extrapol            Default: False
                              True - value weighted by non NaN fractions
                              False - NaN if any fraction is NaN
        filOUT              output file
        ------ OUTPUT ------
        newimage            rebinned image array
        '''
        oldimage = self.images
        newheader = self.header
        oldheader = self.header.copy()
        oldw = self.w
        # cd = w.pixel_scale_matrix
        oldcd = self.cd
        oldcdelt = self.cdelt
        oldNx = self.Nx
        oldNy = self.Ny
        
        if pixscale is not None:
            pixscale = LA.listize(pixscale)
            if len(pixscale)==1:
                pixscale.extend(pixscale)
            else:
                warnings.warn('Non-square pixels present as square on DS9. '
                              'WCS will not work either.')
            ## convert arcsec to degree
            cdelt = np.array(pixscale) / 3600.
            ## Expansion (>1) or contraction (<1) in X/Y
            xratio = cdelt[0] / abs(oldcdelt[0])
            yratio = cdelt[1] / abs(oldcdelt[1])
        else:
            pixscale = LA.listize(abs(oldcdelt) * 3600.)
            xratio = 1.
            yratio = 1.

            if self.verbose==True:
                print('----------')
                print('The actual map size is {} * {}'.format(self.Nx, self.Ny))
                print('The actual pixel scale is {} * {} arcsec'.format(*pixscale))
                print('----------')

            UT.strike('improve.rebin', 'miss pixscale. Nothing has been done.',
                      cat='InputError')

        ## Modify header
        ##---------------
        ## Fix CRVALn
        crpix1 = newheader['CRPIX1']
        crpix2 = newheader['CRPIX2']
        newheader['CRPIX1'] = (crpix1 - 0.5) / xratio + 0.5
        newheader['CRPIX2'] = (crpix2 - 0.5) / yratio + 0.5
    
        cd = oldcd * [xratio,yratio]
        newheader['CD1_1'] = cd[0][0]
        newheader['CD2_1'] = cd[1][0]
        newheader['CD1_2'] = cd[0][1]
        newheader['CD2_2'] = cd[1][1]
    
        for kw in oldheader.keys():
            if 'PC' in kw:
                del newheader[kw]
            if 'CDELT' in kw:
                del newheader[kw]
            
        # lam = yratio/xratio
        # pix_ratio = xratio*yratio
        Nx = math.ceil(oldNx / xratio)
        Ny = math.ceil(oldNy / yratio)
        # Nx = int(oldNx/xratio + 0.5)
        # Ny = int(oldNy/yratio + 0.5)

        ## Rebin
        ##-------
        '''
        ## Ref: poppy(v0.3.4).utils.krebin
        ## Klaus P's fastrebin from web
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
        return a.reshape(sh).sum(-1).sum(1)
        '''

        if self.Ndim==3:
            image_newx = np.zeros((self.Nw,oldNy,Nx))
            newimage = np.zeros((self.Nw,Ny,Nx))
            nanbox = np.zeros((self.Nw,Ny,Nx))
        elif self.Ndim==2:
            image_newx = np.zeros((oldNy,Nx))
            newimage = np.zeros((Ny,Nx))
            nanbox = np.zeros((Ny,Nx))

        ## istart/old1, istop/old2, rstart/new1, rstop/new2 are old grid indices

        if not extrapol:
            
            ## Sample x axis
            ##---------------
            for x in range(Nx):
                rstart = x * xratio # float
                istart = int(rstart) # int
                frac1 = rstart - istart
                rstop = rstart + xratio # float
                if int(rstop)<oldNx:
                    ## Full covered new pixels
                    istop = int(rstop) # int
                    frac2 = 1. - (rstop - istop)
                else:
                    ## Upper edge (value 0 for uncovered frac: frac2)
                    istop = oldNx - 1 # int
                    frac2 = 0
            
                if istart==istop:
                    ## Shrinking case with old pix containing whole new pix (box averaging)
                    if self.Ndim==3:
                        image_newx[:,:,x] = (1.-frac1-frac2) * oldimage[:,:,istart]
                    elif self.Ndim==2:
                        image_newx[:,x] = (1.-frac1-frac2) * oldimage[:,istart]
                else:
                    ## Other cases (bilinear interpolation)
                    if self.Ndim==3:
                        edges = frac1*oldimage[:,:,istart] + frac2*oldimage[:,:,istop]
                        image_newx[:,:,x] = np.sum(oldimage[:,:,istart:istop+1],axis=2) - edges
                    elif self.Ndim==2:
                        edges = frac1*oldimage[:,istart] + frac2*oldimage[:,istop]
                        image_newx[:,x] = np.sum(oldimage[:,istart:istop+1],axis=1) - edges
                        
            ## Sample y axis
            ##---------------
            for y in range(Ny):
                rstart = y * yratio # float
                istart = int(rstart) # int
                frac1 = rstart - istart
                rstop = rstart + yratio # float
                if int(rstop)<oldNy:
                    ## Full covered new pixels
                    istop = int(rstop) # int
                    frac2 = 1. - (rstop - istop)
                else:
                    ## Upper edge (value 0 for uncovered frac: frac2)
                    istop = oldNy - 1 # int
                    frac2 = 0
            
                if istart==istop:
                    ## Shrinking case with old pix containing whole new pix (box averaging)
                    if self.Ndim==3:
                        newimage[:,y,:] = (1.-frac1-frac2) * image_newx[:,istart,:]
                    elif self.Ndim==2:
                        newimage[y,:] = (1.-frac1-frac2) * image_newx[istart,:]
                else:
                    ## Other cases (bilinear interpolation)
                    if self.Ndim==3:
                        edges = frac1*image_newx[:,istart,:] + frac2*image_newx[:,istop,:]
                        newimage[:,y,:] = np.sum(image_newx[:,istart:istop+1,:],axis=1) - edges
                    elif self.Ndim==2:
                        edges = frac1*image_newx[istart,:] + frac2*image_newx[istop,:]
                        newimage[y,:] = np.sum(image_newx[istart:istop+1,:],axis=0) - edges

            if not total:
                newimage = newimage / (xratio*yratio)

        else:
            
            ## Sample y axis
            ##---------------
            for y in range(Ny):
                rstart = y * yratio # float
                istart = int(rstart) # int
                frac1 = rstart - istart
                rstop = rstart + yratio # float
                if int(rstop)<oldNy:
                    ## Full covered new pixels
                    istop = int(rstop) # int
                    frac2 = 1. - (rstop - istop)
                else:
                    ## Upper edge (value 0 for uncovered frac: frac2)
                    istop = oldNy - 1 # int
                    frac2 = (rstop - istop) - 1.
    
                ## Sample x axis
                ##---------------
                for x in range(Nx):
                    new1 = x * xratio # float
                    old1 = int(new1) # int
                    f1 = new1 - old1
                    new2 = new1 + xratio # float
                    if int(new2)<oldNx:
                        ## Full covered new pixels
                        old2 = int(new2) # int
                        f2 = 1. - (new2 - old2)
                    else:
                        ## Upper edge (value 0 for uncovered frac: f2)
                        old2 = oldNx - 1 # int
                        f2 = (new2 - old2) - 1. # out frac

                    ## For each pixel (x,y) in new grid,
                    ## find NaNs in old grid and
                    ## recalculate nanbox[w,y,x] taking into account fractions
                    for j in range(istop+1-istart):
                        for i in range(old2+1-old1):
                                
                            ## old y grid
                            if j==0:
                                ybox = 1.-frac1
                            elif j==istop-istart:
                                if int(rstop)<oldNy:
                                    ybox = 1.-frac2
                                else:
                                    ybox = rstop-istop-1.
                            else:
                                ybox = 1.
                                
                            ## old x grid
                            if i==0:
                                xbox = 1.-f1
                            elif i==old2-old1:
                                if int(new2)<oldNx:
                                    xbox = 1.-f2
                                else:
                                    xbox = f2
                            else:
                                xbox = 1.
                                
                            ## old 2D grid
                            if self.Ndim==3:
                                for w in range(self.Nw):
                                    if ~np.isnan(oldimage[w,istart+j,old1+i]):
                                        newimage[w,y,x] += oldimage[w,istart+j,old1+i] * ybox * xbox
                                        nanbox[w,y,x] += ybox * xbox
                            elif self.Ndim==2:
                                if ~np.isnan(oldimage[istart+j,old1+i]):
                                    newimage[y,x] += oldimage[istart+j,old1+i] * ybox * xbox
                                    nanbox[y,x] += ybox * xbox

            if not total:
                newimage = np.where(nanbox==0, np.nan, newimage/nanbox)
                newimage[newimage==0] = np.nan
            
        if filOUT is not None:
            IO.write_fits(filOUT, header=newheader, data=newimage,
                          wave=self.wave, wmod=self.wmod)

        ## Update self variables
        self.reinit(header=newheader, images=newimage, wave=self.wave,
                    wmod=self.wmod, verbose=self.verbose)

        if self.verbose==True:
            print('----------')
            print(f'The actual map size is {self.Nx} * {self.Ny}')
            print('The actual pixel scale is {} * {} arcsec'.format(*pixscale))
            print('----------')
            
        return newimage

    def groupixel(self, xscale=1, yscale=1, filOUT=None):
        '''
        Group a cluster of pixels (with their mean value)

        ------ INPUT ------
        xscale,yscale       grouped super pixel size
        '''
        Nxs = math.ceil(self.Nx/xscale)
        Nys = math.ceil(self.Ny/yscale)
        if self.Ndim==3:
            ## Super pixels
            image_sup = np.zeros((self.Nw,Nys,Nxs))
            for xs in range(Nxs):
                xarr = LA.zoom(xs, 1/xscale)
                for ys in range(Nys):
                    yarr = LA.zoom(ys, 1/yscale)
                    im = self.images[:,yarr[0]:yarr[-1]+1,xarr[0]:xarr[-1]+1]
                    image_sup[:,ys,xs] += np.nanmean(im,axis=(1,2))
            ## Grouped pixels
            image_grp = np.zeros((self.Nw,self.Ny,self.Nx))
            for x in range(self.Nx):
                for y in range(self.Ny):
                    xs = LA.zoom(x, xscale)
                    ys = LA.zoom(y, yscale)
                    image_grp[:,y,x] = image_sup[:,ys,xs]
        elif self.Ndim==2:
            ## Super pixels
            image_sup = np.zeros((Nys,Nxs))
            for xs in range(Nxs):
                xarr = LA.zoom(xs, 1/xscale)
                for ys in range(Nys):
                    yarr = LA.zoom(ys, 1/yscale)
                    im = self.images[yarr[0]:yarr[-1]+1,xarr[0]:xarr[-1]+1]
                    image_sup[ys,xs] += np.nanmean(im)
            ## Grouped pixels
            image_grp = np.zeros((self.Ny,self.Nx))
            for x in range(self.Nx):
                for y in range(self.Ny):
                    xs = LA.zoom(x, xscale)
                    ys = LA.zoom(y, yscale)
                    image_grp[y,x] = image_sup[ys,xs]

        if filOUT is not None:
            IO.write_fits(filOUT, self.header, image_grp, self.wave, self.wmod)
            
        ## Update self variables
        self.reinit(header=self.header, images=image_grp, wave=self.wave,
                    wmod=self.wmod, verbose=self.verbose)

        return image_grp

    def smooth(self, smooth=1, wgrid=None, wstart=None, filOUT=None):
        '''
        Smooth wavelengths
        If shift, not compatible with unc which needs MC propagation

        See also specutils.smoothing
        ------ INPUT ------
        smooth              smooth wavelength grid by linear interpolation (Default: 1)
        wgrid               external wavelength grid (Default: None)
        wstart              shift wavelength grid to wstart origin (Default: None)
        ------ OUTPUT ------
        newimage            wavelength smoothed spectral cube
        '''
        ## Replace wavelength grid
        if wgrid is not None:
            wvl = wgrid
            Nw0 = len(wgrid)
        else:
            wvl = self.wave
            Nw0 = self.Nw
            
        ## Wavelength shift (within original w range)
        if wstart is not None:
            wshift = wstart - wvl[0]
        else:
            wshift = 0
            
        newave = []
        nan_left = 0
        nan_right = 0
        for k in range(Nw0):
            if k%smooth==0:
                w = wvl[k]+wshift
                ## New wgrid should be within the interpolation range (or give NaNs)
                newave.append(w)
                if w<self.wave[0]:
                    nan_left+=1
                elif w>self.wave[-1]:
                    nan_right+=1
        newave = np.array(newave)
        Nw = len(newave)
        newimage = np.empty([Nw,self.Ny,self.Nx])
        for x in range(self.Nx):
            for y in range(self.Ny):
                f = interp1d(self.wave, self.images[:,y,x], kind='linear')
                newimage[nan_left:Nw-nan_right,y,x] = f(newave[nan_left:Nw-nan_right])
                newimage[:nan_left,y,x] = np.nan
                newimage[Nw-nan_right:,y,x] = np.nan

        if filOUT is not None:
            IO.write_fits(filOUT, self.header, newimage, newave, self.wmod)
            
        ## Update self variables
        self.reinit(header=self.header, images=newimage, wave=newave,
                    wmod=self.wmod, verbose=self.verbose)

        return newimage

    def artifact(self, filUNC=None, BG_images=None, fill_zeros=np.nan,
                 wmin=None, wmax=None, lim_unc=1.e2, fltr_pn=None, cmin=5,
                 filOUT=None):
        '''
        Remove spectral artifacts (Interpolate aberrant wavelengths)
        Anormaly if:
          abs(v - v_med) / unc > lim_unc

        See also specutils.smoothing
        ------ INPUT ------
        filUNC              input uncertainty map (FITS)
        filOUT              output spectral map (FITS)
        BG_images           background images used to generate unc map
        fill_zeros          value used to replace zero value (Default:NaN)
        wmin,wmax           wavelength range to clean (float)
        lim_unc             uncertainty dependant factor limit (positive float)
        fltr_pn             positive/negtive filter (Default: None)
                              'p' - clean only positive aberrant
                              'n' - clean only negtive aberrant
        cmin                minimum neighboring artifacts
        ------ OUTPUT ------
        im                  cleaned spectral map
        '''
        im = self.images
        wvl = self.wave
        unc = self.BGunc(filUNC=filUNC,BG_images=BG_images,fill_zeros=fill_zeros)

        if wmin is None:
            wmin = wvl[0]
        iwi = LA.listize(wvl).index(wvl[LA.closest(wvl,wmin)])
        if wmax is None:
            wmax = wvl[-1]
        iws = LA.listize(wvl).index(wvl[LA.closest(wvl,wmax)])

        if lim_unc<0:
            raise ValueError('lim_unc must be positive!')

        ## Scan every pixel/spectrum at each wavelength
        for w in trange(self.Nw, leave=False,
                        desc='<improve> Cleaning spectral artifacts'):
            if w>=iwi and w<=iws:
                pix_x = []
                pix_y = []
                for y in range(self.Ny):
                    for x in range(self.Nx):
                        v_med = np.median(im[iwi:iws,y,x])
                        dv = (im[w,y,x] - v_med) / unc[w,y,x]
                        if fltr_pn is None or fltr_pn=='p':
                            if dv > lim_unc:
                                pix_x.append(x)
                                pix_y.append(y)
                        if fltr_pn is None or fltr_pn=='n':
                            if dv < -lim_unc:
                                pix_x.append(x)
                                pix_y.append(y)
                pix_x = np.array(pix_x)
                pix_y = np.array(pix_y)
                
                ## If the neighbors share the feature, not an artifact
                for ix, x in enumerate(pix_x):
                    counter = 0
                    for iy, y in enumerate(pix_y):
                        if abs(y-pix_y[ix]+pix_x[iy]-x)<=2:
                            counter += 1
                    ## max(counter) == 12
                    if counter<cmin:
                        if w==0:
                            im[w,pix_y[ix],x] = im[w+1,pix_y[ix],x]
                        elif w==self.Nw-1:
                            im[w,pix_y[ix],x] = im[w-1,pix_y[ix],x]
                        else:
                            im[w,pix_y[ix],x] = (im[w-1,pix_y[ix],x]+im[w+1,pix_y[ix],x])/2
                            # im[w,pix_y[ix],x] = np.median(im[iwi:iws,pix_y[ix],x])

        if filOUT is not None:
            comment = "Cleaned by <improve.artifact>"
            IO.write_fits(filOUT, self.header, im, wvl,
                          COMMENT=comment)

        return im
                
    def mask(self):
        '''
        '''
        pass

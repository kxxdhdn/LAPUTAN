#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Imaging

    improve:
        reinit, uncert, 
        rand_norm, rand_splitnorm, rand_pointing, 
        slice, slice_inv_sq, crop, rebin, groupixel
        smooth, artifact, mask
    Jy_per_pix_to_MJy_per_sr(improve):
        header, image, wave
    iuncert(improve):
        unc
    islice(improve):
        image, wave, filenames, clean
    icrop(improve):
        header, image, wave
    irebin(improve):
        header, image, wave
    igroupixel(improve):
        header, image, wave
    ismooth(improve):
        header, image, wave
    imontage(improve):
        reproject, reproject_mc, coadd, clean
    iswarp(improve):
        footprint, combine, reproject_mc, clean
    iconvolve(improve):
        spitzer_irs, choker, do_conv, image, wave,
        filenames, clean
    cupid(improve):
        spec_build, sav_build,
        header, image, wave
    wmask, wclean, interfill, hextract, hswarp, 
    concatenante

"""

from tqdm import tqdm, trange
import os
import math
import numpy as np
from scipy.io import readsav
from scipy.interpolate import interp1d
from astropy import wcs
from astropy.io import ascii
from astropy.table import Table
from reproject import reproject_interp, reproject_exact, reproject_adaptive
from reproject.mosaicking import reproject_and_coadd
import subprocess as SP
import warnings

## Local
from utilities import InputError
from inout import (fitsext, csvext, ascext, fclean,
                   read_fits, write_fits, savext, write_hdf5,
                   # read_csv, write_csv, read_ascii,
)
from arrays import listize, closest, pix2sup, sup2pix
from maths import nanavg, bsplinterpol
from astrom import fixwcs, get_pc, pix2sr

##-----------------------------------------------
##
##            <improve> based tools
##
##-----------------------------------------------

class improve:
    '''
    IMage PROcessing VEssel
    '''
    def __init__(self, filIN=None, header=None, image=None, wave=None,
                 wmod=0, verbose=False):
        '''
        self: filIN, wmod, hdr, w, cdelt, pc, cd, Ndim, Nx, Ny, Nw, im, wvl
        '''
        
        ## INPUTS
        self.filIN = filIN
        self.wmod = wmod
        self.verbose = verbose

        ## Read image/cube
        if filIN is not None:
            ds = read_fits(filIN)
            self.hdr = ds.header
            self.im = ds.data
            self.wvl = ds.wave
        else:
            self.hdr = header
            self.im = image
            self.wvl = wave
        if self.im is not None:
            self.Ndim = self.im.ndim
            if self.Ndim==3:
                self.Nw, self.Ny, self.Nx = self.im.shape
    
                ## Nw=1 patch
                if self.im.shape[0]==1:
                    self.Ndim = 2
            elif self.Ndim==2:
                self.Ny, self.Nx = self.im.shape
                self.Nw = None

        if self.hdr is not None:
            hdr = self.hdr.copy()
            ws = fixwcs(header=hdr, mode='red_dim')
            self.hdred = ws.header # reduced header
            self.w = ws.wcs
            pcdelt = get_pc(wcs=ws.wcs)
            self.cdelt = pcdelt.cdelt
            self.pc = pcdelt.pc
            self.cd = pcdelt.cd
        
        if verbose==True:
            print('<improve> file: ', filIN)
            print('Raw size (pix): {} * {}'.format(self.Nx, self.Ny))

    def reinit(self, filIN=None, header=None, image=None, wave=None,
               wmod=0, verbose=False):
        '''
        Update init variables
        '''
        ## INPUTS
        self.filIN = filIN
        self.wmod = wmod
        self.verbose = verbose

        ## Read image/cube
        if filIN is not None:
            ds = read_fits(filIN)
            self.hdr = ds.header
            self.im = ds.data
            self.wvl = ds.wave
        else:
            self.hdr = header
            self.im = image
            self.wvl = wave
        self.Ndim = self.im.ndim
        self.hdr['NAXIS'] = self.Ndim
        if self.Ndim==3:
            self.Nw, self.Ny, self.Nx = self.im.shape
            ## Nw=1 patch
            if self.im.shape[0]==1:
                self.Ndim = 2
                del hdr['NAXIS3']
            else:
                self.hdr['NAXIS3'] = self.Nw
        elif self.Ndim==2:
            self.Ny, self.Nx = self.im.shape
            self.Nw = None
        self.hdr['NAXIS2'] = self.Ny
        self.hdr['NAXIS1'] = self.Nx

        hdr = self.hdr.copy()
        ws = fixwcs(header=hdr, mode='red_dim')
        self.hdred = ws.header # reduced header
        self.w = ws.wcs
        pcdelt = get_pc(wcs=ws.wcs)
        self.cdelt = pcdelt.cdelt
        self.pc = pcdelt.pc
        self.cd = pcdelt.cd
        
        if verbose==True:
            print('<improve> file: ', filIN)
            print('Image size (pix): {} * {}'.format(self.Nx, self.Ny))

    def uncert(self, filOUT=None, filUNC=None, filWGT=None, wfac=1.,
               BG_image=None, BG_weight=None, zerovalue=np.nan):
        '''
        Estimate uncertainties from the background map
        So made error map is uniform/weighted

        ------ INPUT ------
        filOUT              output uncertainty map (FITS)
        filUNC              input uncertainty map (FITS)
        filWGT              input weight map (FITS)
        wfac                multiplication factor for filWGT (Default: 1)
        BG_image            background image array used to generate unc map
        BG_weight           background weight array
        zerovalue           value used to replace zero value (Default: NaN)
        ------ OUTPUT ------
        unc                 estimated unc map
        '''
        if filUNC is not None:
            unc = read_fits(filUNC).data
        else:
            if BG_image is not None:
                im = BG_image
                Ny, Nx = BG_image.shape
            else:
                im = self.im
                Ny = self.Ny
                Nx = self.Nx
            Nw = self.Nw

            ## sigma: std dev of (weighted) flux distribution of bg region
            if BG_weight is not None:
                if self.Ndim==3:
                    sigma = np.nanstd(im * BG_weight, axis=(1,2))
                else:
                    sigma = np.nanstd(im * BG_weight)
            else:
                if self.Ndim==3:
                    sigma = np.nanstd(im, axis=(1,2))
                else:
                    sigma = np.nanstd(im)

            ## wgt: weight map
            if filWGT is not None:
                wgt = read_fits(filWGT).data * wfac
            else:
                wgt = np.ones(self.im.shape) * wfac

            ## unc: weighted rms = root of var/wgt
            if self.Ndim==3:
                unc = []
                for w in range(Nw):
                    unc.append(np.sqrt(1./wgt[w,:,:]) * sigma(w))
                unc = np.array(unc)
            else:
                unc = np.sqrt(1./wgt) * sigma

            ## Replace zero values
            unc[unc==0] = zerovalue

        self.unc = unc
        
        if filOUT is not None:
            write_fits(filOUT, self.hdr, unc, self.wvl, self.wmod)
            
        return unc

    def rand_norm(self, filUNC=None, unc=None, sigma=1., mu=0.):
        '''
        Add random N(0,1) noise
        '''
        if filUNC is not None:
            unc = read_fits(filUNC).data

        if unc is not None:
            ## unc should have the same dimension with im
            theta = np.random.normal(mu, sigma, self.im.shape)
            self.im += theta * unc

        return self.im

    def rand_splitnorm(self, filUNC=None, unc=None, sigma=1., mu=0.):
        '''
        Add random SN(0,lam,lam*tau) noise

        ------ INPUT ------
        filUNC              2 FITS files for unc of left & right sides
        unc                 2 uncertainty ndarrays
        ------ OUTPUT ------
        '''
        if filUNC is not None:
            unc = []
            for f in filUNC:
                unc.append(read_fits(f).data)
            
        if unc is not None:
            ## unc[i] should have the same dimension with self.im
            tau = unc[1]/unc[0]
            peak = 1/(1+tau)
            theta = np.random.normal(mu, sigma, self.im.shape) # ~N(0,1)
            flag = np.random.random(self.im.shape) # ~U(0,1)
            if self.Ndim==2:
                for x in range(self.Nx):
                    for y in range(self.Ny):
                        if flag[y,x]<peak[y,x]:
                            self.im[y,x] += -abs(theta[y,x]) * unc[0][y,x]
                        else:
                            self.im[y,x] += abs(theta[y,x]) * unc[1][y,x]
            elif self.Ndim==3:
                for x in range(self.Nx):
                    for y in range(self.Ny):
                        for k in range(self.Nw):
                            if flag[k,y,x]<peak[k,y,x]:
                                self.im[k,y,x] += -abs(
                                    theta[k,y,x]) * unc[0][k,y,x]
                            else:
                                self.im[k,y,x] += abs(
                                    theta[k,y,x]) * unc[1][k,y,x]

        return self.im

    def rand_pointing(self, sigma=0, header=None, fill='med',
                      xscale=1, yscale=1, swarp=False, tmpdir=None):
        '''
        Add pointing uncertainty to WCS

        ------ INPUT ------
        sigma               pointing accuracy (arcsec)
        header              baseline
        fill                fill value of no data regions after shift
                              'med': axis median (default)
                              'avg': axis average
                              float: constant
        xscale,yscale       regrouped super pixel size
        swarp               use SWarp to perform position shifts
                              Default: False (not support supix)
        ------ OUTPUT ------
        '''
        ## Set path of tmp files
        if tmpdir is None:
            path_tmp = os.getcwd()+'/tmp_swp/'
        else:
            path_tmp = tmpdir
        if not os.path.exists(path_tmp):
            os.makedirs(path_tmp)
        
        if sigma!=0:
            sigma /= 3600.
            d_ro = abs(np.random.normal(0., sigma)) # N(0,sigma)
            d_phi = np.random.random() *2. * np.pi # U(0,2*pi)
            # d_ro, d_phi = 0.0002, 4.5
            # print('d_ro,d_phi = ', d_ro,d_phi)
            ## New header/WCS
            if header is None:
                header = self.hdr
            wcs = fixwcs(header=header, mode='red_dim').wcs
            Nx = header['NAXIS1']
            Ny = header['NAXIS2']
            newheader = header.copy()
            newheader['CRVAL1'] += d_ro * np.cos(d_phi)
            newheader['CRVAL2'] += d_ro * np.sin(d_phi)
            newcs = fixwcs(header=newheader, mode='red_dim').wcs
    
            ## Convert world increment to pix increment
            pix = wcs.all_world2pix(newheader['CRVAL1'], newheader['CRVAL2'], 1)
            d_x = pix[0] - header['CRPIX1']
            d_y = pix[1] - header['CRPIX2']
            # print('Near CRPIXn increments: ', d_x, d_y)
            # val1 = np.array(newcs.all_pix2world(0.5, 0.5, 1))
            # d_x, d_y = wcs.all_world2pix(val1[np.newaxis,:], 1)[0] - 0.5
            # print('Near (1,1) increments: ', d_x, d_y)

            ## Resampling
            if swarp:
                write_fits(path_tmp+'tmp_rand_shift',
                           newheader, self.im, self.wvl)
                swp = iswarp(refheader=self.hdr, tmpdir=path_tmp)
                rep = swp.combine(path_tmp+'tmp_rand_shift',
                                  combtype='avg', keepedge=True)
                self.im = rep.image
                self.im[np.isnan(rep.image)] = 0
            else:
                Nxs = math.ceil(Nx/xscale)
                cube_supx = np.zeros((self.Nw,Ny,Nxs))
                frac2 = d_x / xscale
                f2 = math.floor(frac2)
                frac1 = 1 - frac2
                for xs in range(Nxs):
                    if fill=='med':
                        fill_value = np.nanmedian(self.im,axis=2)
                    elif fill=='avg':
                        fill_value = np.nanmean(self.im,axis=2)
                    else:
                        fill_value = fill
                    if frac2>=0:
                        if xs>=f2:
                            x0 = sup2pix(xs-f2, xscale, Npix=Nx, origin=0)
                            cube_supx[:,:,xs] += (f2+frac1) * np.nanmean(self.im[:,:,x0[0]:x0[-1]+1],axis=2)
                            if xs>f2:
                                x1 = sup2pix(xs-f2-1, xscale, Npix=Nx, origin=0)
                                cube_supx[:,:,xs] += (frac2-f2) * np.nanmean(self.im[:,:,x1[0]:x1[-1]+1],axis=2)
                            else:
                                cube_supx[:,:,xs] += (frac2-f2) * fill_value
                        else:
                            cube_supx[:,:,xs] += fill_value
                            # if self.verbose:
                            #     warnings.warn('Zero appears at super x = {}'.format(xs))
                    else:
                        if xs<=Nxs+f2:
                            x1 = sup2pix(xs-f2-1, xscale, Npix=Nx, origin=0)
                            cube_supx[:,:,xs] += (frac2-f2) * np.nanmean(self.im[:,:,x1[0]:x1[-1]+1],axis=2)
                            if xs<Nxs+f2:
                                x0 = sup2pix(xs-f2, xscale, Npix=Nx, origin=0)
                                cube_supx[:,:,xs] += (f2+frac1) * np.nanmean(self.im[:,:,x0[0]:x0[-1]+1],axis=2)
                            else:
                                cube_supx[:,:,xs] += (f2+frac1) * fill_value
                        else:
                            cube_supx[:,:,xs] += fill_value
                            # if self.verbose:
                            #     warnings.warn('Zero appears at super x = {}'.format(xs))
                
                Nys = math.ceil(Ny/yscale)
                supcube = np.zeros((self.Nw,Nys,Nxs))
                frac2 = d_y / yscale
                f2 = math.floor(frac2)
                frac1 = 1 - frac2
                for ys in range(Nys):
                    if fill=='med':
                        fill_value = np.nanmedian(cube_supx,axis=1)
                    elif fill=='avg':
                        fill_value = np.nanmean(cube_supx,axis=1)
                    else:
                        fill_value = fill
                    if frac2>=0:
                        if ys>=f2:
                            y0 = sup2pix(ys-f2, yscale, Npix=Ny, origin=0)
                            supcube[:,ys,:] += (f2+frac1) * np.nanmean(cube_supx[:,y0[0]:y0[-1]+1,:],axis=1)
                            if ys>f2:
                                y1 = sup2pix(ys-f2-1, yscale, Npix=Ny, origin=0)
                                supcube[:,ys,:] += (frac2-f2) * np.nanmean(cube_supx[:,y1[0]:y1[-1]+1,:],axis=1)
                            else:
                                supcube[:,ys,:] += (frac2-f2) * fill_value
                        else:
                            supcube[:,ys,:] += fill_value
                            # if self.verbose:
                            #     warnings.warn('Zero appears at super y = {}'.format(ys))
                    else:
                        if ys<=Nys+f2:
                            y1 = sup2pix(ys-f2-1, yscale, Npix=Ny, origin=0)
                            supcube[:,ys,:] += (frac2-f2) * np.nanmean(cube_supx[:,y1[0]:y1[-1]+1,:],axis=1)
                            if ys<Nys+f2:
                                y0 = sup2pix(ys-f2, yscale, Npix=Ny, origin=0)
                                supcube[:,ys,:] += (f2+frac1) * np.nanmean(cube_supx[:,y0[0]-1:y0[-1],:],axis=1)
                            else:
                                supcube[:,ys,:] += (f2+frac1) * fill_value
                        else:
                            supcube[:,ys,:] += fill_value
                            # if self.verbose:
                            #     warnings.warn('Zero appears at super y = {}'.format(ys))
                
                for x in range(Nx):
                    for y in range(Ny):
                        xs = pix2sup(x, xscale, origin=0)
                        ys = pix2sup(y, yscale, origin=0)
                        self.im[:,y,x] = supcube[:,ys,xs]

        return self.im

    def slice(self, filSL, postfix='', ext=''):
        ## 3D cube slicing
        slist = []
        if self.Ndim==3:
            # hdr = self.hdr.copy()
            # for kw in self.hdr.keys():
            #     if '3' in kw:
            #         del hdr[kw]
            # hdr['NAXIS'] = 2
            for k in range(self.Nw):
                ## output filename list
                f = filSL+'_'+'0'*(4-len(str(k)))+str(k)+postfix
                slist.append(f+ext)
                write_fits(f, self.hdred, self.im[k,:,:]) # gauss_noise inclu
        else:
            f = filSL+'_0000'+postfix
            slist.append(f+ext)
            write_fits(f, self.hdred, self.im) # gauss_noise inclu
            if self.verbose==True:
                print('Input file is a 2D image which cannot be sliced! ')
                print('Rewritten with only random noise added (if provided).')

        return slist

    def slice_inv_sq(self, filSL, postfix=''):
        ## Inversed square cube slicing
        inv_sq = 1./self.im**2
        slist = []
        if self.Ndim==3:
            # hdr = self.hdr.copy()
            # for kw in self.hdr.keys():
            #     if '3' in kw:
            #         del hdr[kw]
            # hdr['NAXIS'] = 2
            for k in range(self.Nw):
                ## output filename list
                f = filSL+'_'+'0'*(4-len(str(k)))+str(k)+postfix
                slist.append(f)
                write_fits(f, self.hdred, inv_sq[k,:,:]) # gauss_noise inclu
        else:
            f = filSL+'_0000'+postfix
            slist.append(f)
            write_fits(f, self.hdred, inv_sq) # gauss_noise inclu

        return slist
    
    def crop(self, filOUT=None,
             sizpix=None, cenpix=None, sizval=None, cenval=None):
        '''
        If pix and val co-exist, pix will be taken.

        ------ INPUT ------
        filOUT              output file
        sizpix              crop size in pix (dx, dy)
        cenpix              crop center in pix (x, y)
        sizval              crop size in deg (dRA, dDEC) -> (dx, dy)
        cenval              crop center in deg (RA, DEC) -> (x, y)
        ------ OUTPUT ------
        self.im             cropped image array
        '''
        oldimage = self.im
        hdr = self.hdr
        
        ## Crop center
        ##-------------
        if cenpix is None:
            if cenval is None:
                raise ValueError('Crop center unavailable! ')
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
            raise ValueError('Crop centre overpassed image border! ')

        ## Crop size
        ##-----------
        if sizpix is None:
            if sizval is None:
                raise ValueError('Crop size unavailable! ')
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
            raise ValueError('Crop region overpassed image border! ')

        ## OUTPUTS
        ##---------
        ## New image
        if self.Ndim==3:
            newimage = oldimage[:, ymin:ymax, xmin:xmax] # gauss_noise inclu
            ## recover 3D non-reduced header
            # hdr = read_fits(self.filIN).header
        elif self.Ndim==2:
            newimage = oldimage[ymin:ymax, xmin:xmax] # gauss_noise inclu

        ## Modify header
        ##---------------
        hdr['CRPIX1'] = math.floor(sizpix[0]/2. + 0.5)
        hdr['CRPIX2'] = math.floor(sizpix[1]/2. + 0.5)
        hdr['CRVAL1'] = cenval[0]
        hdr['CRVAL2'] = cenval[1]
        
        self.hdr = hdr
        self.im = newimage
        
        ## Write cropped image/cube
        if filOUT is not None:
            # comment = "[ICROP]ped at centre: [{:.8}, {:.8}]. ".format(*cenval)
            # comment = "with size [{}, {}] (pix).".format(*sizpix)
            write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod)

        ## Update self variables
        self.reinit(header=self.hdr, image=self.im, wave=self.wvl,
                    wmod=self.wmod, verbose=self.verbose)

        return self.im

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
                              tuple - same Ndim with image
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
        oldimage = self.im
        hdr = self.hdr
        oldheader = hdr.copy()
        oldw = self.w
        # cd = w.pixel_scale_matrix
        oldcd = self.cd
        oldcdelt = self.cdelt
        oldNx = self.Nx
        oldNy = self.Ny
        
        if pixscale is not None:
            pixscale = listize(pixscale)
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
            pixscale = listize(abs(oldcdelt) * 3600.)
            xratio = 1.
            yratio = 1.

            if self.verbose==True:
                print('----------')
                print('The actual map size is {} * {}'.format(self.Nx, self.Ny))
                print('The actual pixel scale is {} * {} arcsec'.format(*pixscale))
                print('----------')
                
            raise InputError('<improve.rebin>',
                             'No pixscale, nothing has been done!')

        ## Modify header
        ##---------------

        ## Fix CRVALn
        crpix1 = hdr['CRPIX1']
        crpix2 = hdr['CRPIX2']
        hdr['CRPIX1'] = (crpix1 - 0.5) / xratio + 0.5
        hdr['CRPIX2'] = (crpix2 - 0.5) / yratio + 0.5
    
        cd = oldcd * [xratio,yratio]
        hdr['CD1_1'] = cd[0][0]
        hdr['CD2_1'] = cd[1][0]
        hdr['CD1_2'] = cd[0][1]
        hdr['CD2_2'] = cd[1][1]
    
        for kw in oldheader.keys():
            if 'PC' in kw:
                del hdr[kw]
            if 'CDELT' in kw:
                del hdr[kw]
            
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
        else:
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
                    else:
                        image_newx[:,x] = (1.-frac1-frac2) * oldimage[:,istart]
                else:
                    ## Other cases (bilinear interpolation)
                    if self.Ndim==3:
                        edges = frac1*oldimage[:,:,istart] + frac2*oldimage[:,:,istop]
                        image_newx[:,:,x] = np.sum(oldimage[:,:,istart:istop+1],axis=2) - edges
                    else:
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
                    else:
                        newimage[y,:] = (1.-frac1-frac2) * image_newx[istart,:]
                else:
                    ## Other cases (bilinear interpolation)
                    if self.Ndim==3:
                        edges = frac1*image_newx[:,istart,:] + frac2*image_newx[:,istop,:]
                        newimage[:,y,:] = np.sum(image_newx[:,istart:istop+1,:],axis=1) - edges
                    else:
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
                            else:
                                if ~np.isnan(oldimage[istart+j,old1+i]):
                                    newimage[y,x] += oldimage[istart+j,old1+i] * ybox * xbox
                                    nanbox[y,x] += ybox * xbox

            if not total:
                newimage = np.where(nanbox==0, np.nan, newimage/nanbox)
                newimage[newimage==0] = np.nan
            
        if filOUT is not None:
            write_fits(filOUT, hdr, newimage, self.wvl, self.wmod)

        ## Update self variables
        self.reinit(header=hdr, image=newimage, wave=self.wvl,
                    wmod=self.wmod, verbose=self.verbose)

        if self.verbose==True:
            print('----------')
            print('The actual map size is {} * {}'.format(self.Nx, self.Ny))
            print('The actual pixel scale is {} * {} arcsec'.format(*pixscale))
            print('\n <improve> Rebin [done]')
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
            image_supx = np.zeros((self.Nw,self.Ny,Nxs))
            for xs in range(Nxs):
                xarr = sup2pix(xs, xscale, Npix=self.Nx, origin=0)
                print('xarr',xarr)
                image_supx[:,:,xs] += np.nanmean(self.im[:,:,xarr[0]:xarr[-1]+1],axis=2)
            supimage = np.zeros((self.Nw,Nys,Nxs))
            for ys in range(Nys):
                yarr = sup2pix(ys, yscale, Npix=self.Ny, origin=0)
                print('yarr',yarr)
                supimage[:,ys,:] += np.nanmean(image_supx[:,yarr[0]:yarr[-1]+1,:],axis=1)
            ## Grouped pixels
            image_grp = np.zeros((self.Nw,self.Ny,self.Nx))
            for x in range(self.Nx):
                for y in range(self.Ny):
                    xs = pix2sup(x, xscale, origin=0)
                    ys = pix2sup(y, yscale, origin=0)
                    image_grp[:,y,x] = supimage[:,ys,xs]
        elif self.Ndim==2:
            ## Super pixels
            image_supx = np.zeros((self.Ny,Nxs))
            for xs in range(Nxs):
                xarr = sup2pix(xs, xscale, Npix=self.Nx, origin=0)
                image_supx[:,xs] += np.nanmean(self.im[:,xarr[0]:xarr[-1]+1],axis=1)
            supimage = np.zeros((Nys,Nxs))
            for ys in range(Nys):
                yarr = sup2pix(ys, yscale, Npix=self.Ny, origin=0)
                supimage[ys,:] += np.nanmean(image_supx[yarr[0]:yarr[-1]+1,:],axis=0)
            ## Grouped pixels
            image_grp = np.zeros((self.Ny,self.Nx))
            for x in range(self.Nx):
                for y in range(self.Ny):
                    xs = pix2sup(x, xscale, origin=0)
                    ys = pix2sup(y, yscale, origin=0)
                    image_grp[y,x] = supimage[ys,xs]

        if filOUT is not None:
            write_fits(filOUT, self.hdr, image_grp, self.wvl, self.wmod)
            
        ## Update self variables
        self.reinit(header=self.hdr, image=image_grp, wave=self.wvl,
                    wmod=self.wmod, verbose=self.verbose)

        return image_grp
    
    def smooth(self, smooth=1, wstart=None, filOUT=None):
        '''
        Smooth wavelengths
        If shift, not compatible with unc which needs MC propagation

        ------ INPUT ------
        smooth              smooth wavelength grid by linear interpolation (Default: 1)
        wstart              shift wavelength grid to wstart origin (Default: None)
        '''
        ## Wavelength shift (within original w range)
        if wstart is not None:
            wshift = wstart - self.wvl[0]
        else:
            wshift = 0
            
        newave = []
        for k in range(self.Nw):
            if k%smooth==0:
                w = self.wvl[k]+wshift
                if (w>=self.wvl[0] and w<=self.wvl[-1]):
                    newave.append(w)
        newave = np.array(newave)
        Nw = len(newave)
        newcube = np.empty([Nw,self.Ny,self.Nx])
        for x in range(self.Nx):
            for y in range(self.Ny):
                f = interp1d(self.wvl, self.im[:,y,x], kind='linear')
                newcube[:,y,x] = f(newave)

        if filOUT is not None:
            write_fits(filOUT, self.hdr, newcube, newave, self.wmod)
            
        ## Update self variables
        self.reinit(header=self.hdr, image=newcube, wave=newave,
                    wmod=self.wmod, verbose=self.verbose)

        return newcube

    def artifact(self, filUNC=None, BG_image=None, zerovalue=np.nan,
                 wmin=None, wmax=None, lim_unc=1.e2, fltr_pn=None, cmin=5,
                 filOUT=None):
        '''
        Remove spectral artifacts (Interpolate aberrant wavelengths)
        Anormaly if:
          abs(v - v_med) / unc > lim_unc

        ------ INPUT ------
        filUNC              input uncertainty map (FITS)
        filOUT              output spectral map (FITS)
        BG_image            background image used to generate unc map
        zerovalue           value used to replace zero value (Default:NaN)
        wmin,wmax           wavelength range to clean (float)
        lim_unc             uncertainty dependant factor limit (positive float)
        fltr_pn             positive/negtive filter (Default: None)
                              'p' - clean only positive aberrant
                              'n' - clean only negtive aberrant
        cmin                minimum neighboring artifacts
        ------ OUTPUT ------
        im                  cleaned spectral map
        '''
        im = self.im
        wvl = self.wvl
        unc = self.uncert(filUNC=filUNC,BG_image=BG_image,zerovalue=zerovalue)

        if wmin is None:
            wmin = wvl[0]
        iwi = listize(wvl).index(wvl[closest(wvl,wmin)])
        if wmax is None:
            wmax = wvl[-1]
        iws = listize(wvl).index(wvl[closest(wvl,wmax)])

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
            write_fits(filOUT, self.hdr, im, wvl,
                       COMMENT=comment)

        return im
                
    def mask(self):
        '''
        '''
        pass

class Jy_per_pix_to_MJy_per_sr(improve):
    '''
    Convert image unit from Jy/pix to MJy/sr

    ------ INPUT ------
    filIN               input FITS file
    filOUT              output FITS file
    ------ OUTPUT ------
    '''
    def __init__(self, filIN, filOUT=None, wmod=0, verbose=False):
        super().__init__(filIN, wmod=wmod, verbose=verbose)

        ## gmean( Jy/MJy / sr/pix )
        ufactor = np.sqrt(np.prod(1.e-6/pix2sr(1., self.cdelt)))
        self.im = self.im * ufactor
        self.hdr['BUNIT'] = 'MJy/sr'

        if filOUT is not None:
            write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod)
            
    def header(self):
        return self.hdr
            
    def image(self):
        return self.im

    def wave(self):
        return self.wvl

class iuncert(improve):
    '''
    Generate uncertainties

    ------ INPUT ------
    filIN               input map (FITS)
    filOUT              output weight map (FITS)
    filWGT              input weight map (FITS)
    wfac                multiplication factor for filWGT (Default: 1)
    BG_image            background image array
    BG_weight           background weight array
    zerovalue           value to replace zeros (Default: NaN)
    ------ OUTPUT ------
    '''
    def __init__(self, filIN, filOUT=None, filWGT=None, wfac=1,
                 BG_image=None, BG_weight=None, zerovalue=np.nan):
        super().__init__(filIN, wmod=0, verbose=False)

        self.uncert(filOUT=filOUT, BG_image=BG_image, zerovalue=zerovalue,
                    filWGT=filWGT, wfac=wfac, BG_weight=BG_weight)

    def unc(self):
        return self.unc

class islice(improve):
    '''
    Slice a cube

    ------ INPUT ------
    filIN               input FITS file
    filSL               ouput path+basename
    filUNC              input uncertainty FITS
    dist                unc pdf
    slicetype           Default: None
                          None - normal slices
                          'inv_sq' - inversed square slices
    postfix             postfix of output slice names
    ------ OUTPUT ------
    self: slist, path_tmp, 
          (filIN, wmod, hdr, w, cdelt, pc, cd, Ndim, Nx, Ny, Nw, im, wvl)
    '''
    def __init__(self, filIN, filSL=None, filUNC=None, dist=None,
                 slicetype=None, postfix=''):
        super().__init__(filIN)

        if filSL is None:
            path_tmp = os.getcwd()+'/tmp_proc/'
            if not os.path.exists(path_tmp):
                os.makedirs(path_tmp)

            filSL = path_tmp+'slice'
        self.filSL = filSL

        if dist=='norm':
            self.rand_norm(filUNC)
        elif dist=='splitnorm':
            self.rand_splitnorm(filUNC)

        if slicetype is None:
            self.slist = self.slice(filSL, postfix) # gauss_noise inclu
        elif slicetype=='inv_sq':
            self.slist = self.slice_inv_sq(filSL, postfix)

    def image(self):
        return self.im

    def wave(self):
        return self.wvl

    def filenames(self):
        return self.slist

    def clean(self, filIN=None):
        if filIN is not None:
            fclean(filIN)
        else:
            fclean(self.filSL+'*')

class icrop(improve):
    '''
    CROP 2D image or 3D cube
    '''
    def __init__(self, filIN, filOUT=None,
                 sizpix=None, cenpix=None, sizval=None, cenval=None,
                 filUNC=None, dist=None, wmod=0, verbose=False):
        ## slicrop: slice 
        super().__init__(filIN, wmod=wmod, verbose=verbose)
        
        if dist=='norm':
            self.rand_norm(filUNC)
        elif dist=='splitnorm':
            self.rand_splitnorm(filUNC)
        
        im_crop = self.crop(filOUT=filOUT, sizpix=sizpix, cenpix=cenpix,
                            sizval=sizval, cenval=cenval) # gauss_noise inclu

    def header(self):
        return self.hdr
    
    def image(self):
        return self.im

    def wave(self):
        return self.wvl

class irebin(improve):
    '''
    REBIN 2D image or 3D cube
    '''
    def __init__(self, filIN, filOUT=None,
                 pixscale=None, total=False, extrapol=False,
                 filUNC=None, dist=None, wmod=0, verbose=False):
        super().__init__(filIN, wmod=wmod, verbose=verbose)
        
        if dist=='norm':
            self.rand_norm(filUNC)
        elif dist=='splitnorm':
            self.rand_splitnorm(filUNC)

        im_rebin = self.rebin(filOUT=filOUT, pixscale=pixscale,
                              total=total, extrapol=extrapol)

    def header(self):
        return self.hdr
        
    def image(self):
        return self.im

    def wave(self):
        return self.wvl

class igroupixel(improve):
    '''
    GROUP a cluster of PIXELs (with their mean value)
    '''
    def __init__(self, filIN, filOUT=None,
                 xscale=1, yscale=1,
                 wmod=0, verbose=False):
        super().__init__(filIN, wmod=wmod, verbose=verbose)

        im_grp = self.groupixel(xscale=xscale, yscale=yscale, filOUT=filOUT)

    def header(self):
        return self.hdr
        
    def image(self):
        return self.im

    def wave(self):
        return self.wvl
    
class ismooth(improve):
    '''
    SMOOTH wavelengths
    '''
    def __init__(self, filIN, filOUT=None,
                 smooth=1, wstart=None,
                 wmod=0, verbose=False):
        super().__init__(filIN, wmod=wmod, verbose=verbose)

        im_smooth = self.smooth(smooth=smooth, wstart=wstart, filOUT=filOUT)

    def header(self):
        return self.hdr
        
    def image(self):
        return self.im

    def wave(self):
        return self.wvl

class imontage(improve):
    '''
    2D image or 3D cube montage toolkit
    Based on reproject v0.7.1 or later

    ------ INPUT ------
    reproject_function  resampling algorithms
                          'interp': fastest (Default)
                          'exact': slowest
                          'adaptive': DeForest2004
    tmpdir              tmp file path
    verbose             (Default: False)
    ------ OUTPUT ------
    '''
    def __init__(self, reproject_function='interp',
                 tmpdir=None, verbose=False):
        '''
        self: func, path_tmp, verbose
        '''
        if reproject_function=='interp':
            self.func = reproject_interp
        elif reproject_function=='exact':
            self.func = reproject_exact
        elif reproject_function=='adaptive':
            self.func = reproject_adaptive
        else:
            raise InputError('<imontage>',
                             'Unknown reprojection !')
        
        ## Set path of tmp files
        if tmpdir is None:
            path_tmp = os.getcwd()+'/tmp_mtg/'
        else:
            path_tmp = tmpdir
        if not os.path.exists(path_tmp):
            os.makedirs(path_tmp)
        self.path_tmp = path_tmp

        ## Verbose
        if verbose==False:
            devnull = open(os.devnull, 'w')
        else:
            devnull = None
        self.verbose = verbose
        self.devnull = devnull
    
    def reproject(self, flist, refheader,
                  filOUT=None, dist=None):
        '''
        Reproject 2D image or 3D cube

        ------ INPUT ------
        flist               FITS files to reproject
        refheader           reprojection header
        filOUT              output FITS file
        dist                uncertainty distribution
                              'norm' - N(0,1)
                              'splitnorm' - SN(0,lam,lam*tau)
        ------ OUTPUT ------
        self.images         reprojected images
        '''
        flist = listize(flist)
        
        # if refheader is None:
        #     raise InputError('<imontage>','No reprojection header!')
        
        images = []
        for f in flist:
            super().__init__(f)

            ## Set tmp and out
            filename = os.path.basename(f)
            if filOUT is None:
                filOUT = self.path_tmp+filename+'_rep'

            ## Uncertainty propagation
            if dist=='norm':
                self.rand_norm(f+'_unc')
            elif dist=='splitnorm':
                self.rand_splitnorm([f+'_unc_N', f+'_unc_P'])
            write_fits(filOUT, self.hdr, self.im, self.wvl, wmod=0)
            
            ## Do reprojection
            ##-----------------
            im = self.func(filOUT+fitsext, refheader)[0]
            images.append(im)
    
            comment = "Reprojected by <imontage>. "
            write_fits(filOUT, refheader, im, self.wvl, wmod=0,
                       COMMENT=comment)
        
        return images

    def reproject_mc(self, filIN, refheader,
                     filOUT=None, dist=None, Nmc=0):
        '''
        Generate Monte-Carlo uncertainties for reprojected input file
        '''
        dataset = type('', (), {})()

        hyperim = [] # [j,(w,)y,x]
        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Reprojection [MC]'):

            if j==0:
                im0 = self.reproject(filIN, refheader, filOUT, dist)[0]
            else:
                hyperim.append(self.reproject(filIN, refheader,
                                              filOUT+'_'+str(j), dist)[0])
        im0 = np.array(im0)
        hyperim = np.array(hyperim)
        unc = np.nanstd(hyperim, axis=0)
        comment = "Reprojected by <imontage>. "

        if Nmc>0:
            write_fits(filOUT+'_unc', refheader, unc, self.wvl,
                       COMMENT=comment)

        dataset.im0 = im0
        dataset.unc = unc
        dataset.hyperim = hyperim

        return dataset

    def coadd(self, flist, refheader,
              filOUT=None, dist=None, Nmc=0):
        '''
        Reproject and coadd
        '''
        flist = listize(flist)
        dataset = type('', (), {})()
        comment = "Created by <imontage>"

        slist = [] # slist[j,if,iw]
        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Slicing... [MC]'):
            sl = [] # sl[f,w]
            for f in flist:
                super().__init__(f)

                ## Set tmp and out
                filename = os.path.basename(f)
                if filOUT is None:
                    filOUT = self.path_tmp+filename+'_rep'
                        
                coadd_tmp = self.path_tmp+filename+'/'
                if not os.path.exists(coadd_tmp):
                    os.makedirs(coadd_tmp)
                        
                if j==0:
                    sl.append(self.slice(coadd_tmp+'slice', ext=fitsext))
                else:
                    if dist=='norm':
                        self.rand_norm(f+'_unc')
                    elif dist=='splitnorm':
                        self.rand_splitnorm([f+'_unc_N', f+'_unc_P'])
                        
                    sl.append(self.slice(coadd_tmp+'slice',
                                         postfix='_'+str(j), ext=fitsext))
            slist.append(np.array(sl))
        slist = np.array(slist)
        
        if self.Nw is None:
            Nw = 1
        else:
            Nw = self.Nw
        superim = []
        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Coadding... [MC]'):
            if j==0:
                im = []
                for iw in range(Nw):
                    im.append(reproject_and_coadd(slist[j,:,iw], refheader,
                                                  reproject_function=self.func)[0])
                im = np.array(im)

                write_fits(filOUT, refheader, im, self.wvl, wmod=0,
                           COMMENT=comment)
            else:
                hyperim = []
                for iw in range(Nw):
                    hyperim.append(reproject_and_coadd(slist[j,:,iw], refheader,
                                                       reproject_function=self.func)[0])
                superim.append(np.array(hyperim))

                write_fits(filOUT+'_'+str(j), refheader, hyperim, self.wvl, wmod=0,
                           COMMENT=comment)
        superim = np.array(superim)
        unc = np.nanstd(superim, axis=0)

        if Nmc>0:
            write_fits(filOUT+'_unc', refheader, unc, self.wvl, wmod=0,
                       COMMENT=comment)

        dataset.wvl = self.wvl
        dataset.im = im
        dataset.unc = unc
        dataset.superim = superim
        
        return dataset

    def clean(self, filIN=None):
        if filIN is not None:
            fclean(filIN)
        else:
            fclean(self.path_tmp)
            
class imontage_v0_4(improve):
    '''
    2D image or 3D cube montage toolkit
    (Archived reproject v0.4 version)

    ------ INPUT ------
    flist               FITS file (list, cf improve.filIN)
    filREF              ref file (priority if co-exist with input header)
    hdREF               ref header
    fmod                output image frame mode
                          'ref' - same as ref frame (Default)
                          'rec' - recenter back to input frame
                          'ext' - cover both input and ref frame
    ext_pix             number of pixels to extend to save edge
    tmpdir              tmp file path
    ------ OUTPUT ------
    '''
    def __init__(self, flist, filREF=None, hdREF=None,
                 fmod='ref', ext_pix=0, tmpdir=None):
        '''
        self: hdr_ref, path_tmp, 
        (filIN, wmod, hdr, w, Ndim, Nx, Ny, Nw, im, wvl)
        '''
        ## Set path of tmp files
        if tmpdir is None:
            path_tmp = os.getcwd()+'/tmp_mtg/'
        else:
            path_tmp = tmpdir
        if not os.path.exists(path_tmp):
            os.makedirs(path_tmp)
        self.path_tmp = path_tmp

        ## Inputs
        self.flist = listize(flist)
        self.filREF = filREF
        self.hdREF = hdREF
        self.fmod = fmod
        self.ext_pix = ext_pix
        
        ## Init ref header
        self.hdr_ref = None

    def make_header(self, filIN, filREF=None, hdREF=None, fmod='ref', ext_pix=0):
        '''
        Header maker

        ------ INPUT ------
        filIN               single FITS file
        '''
        super().__init__(filIN)

        ## Prepare reprojection header
        if filREF is not None:
            hdREF = read_fits(filREF).header
            # hdREF['EQUINOX'] = 2000.0

        if hdREF is not None:
            ## Frame mode (fmod) options
            ##---------------------------
            if fmod=='ref':
                pass
            else:
                ## Input WCS (old)
                pix_old = [[0, 0]]
                pix_old.append([0, self.Ny])
                pix_old.append([self.Nx, 0])
                pix_old.append([self.Nx, self.Ny])
                world_arr = self.w.all_pix2world(np.array(pix_old), 1)
                ## Ref WCS (new)
                w = fixwcs(header=hdREF).wcs
                try:
                    pix_new = w.all_world2pix(world_arr, 1)
                except wcs.wcs.NoConvergence as e:
                    pix_new = e.best_solution
                    print("Best solution:\n{0}".format(e.best_solution))
                    print("Achieved accuracy:\n{0}".format(e.accuracy))
                    print("Number of iterations:\n{0}".format(e.niter))
                xmin = min(pix_new[:,0])
                xmax = max(pix_new[:,0])
                ymin = min(pix_new[:,1])
                ymax = max(pix_new[:,1])

                ## Modify ref header
                if fmod=='rec': 
                    hdREF['CRPIX1'] += -xmin
                    hdREF['CRPIX2'] += -ymin
                    hdREF['NAXIS1'] = math.ceil(xmax - xmin)
                    hdREF['NAXIS2'] = math.ceil(ymax - ymin)
                elif fmod=='ext':
                    if xmin<0:
                        hdREF['CRPIX1'] += -xmin
                    if ymin<0:
                        hdREF['CRPIX2'] += -ymin
                    hdREF['NAXIS1'] = math.ceil(max(xmax, hdREF['NAXIS1']-xmin,
                                                    xmax-xmin, hdREF['NAXIS1'])) + ext_pix # save edges
                    hdREF['NAXIS2'] = math.ceil(max(ymax, hdREF['NAXIS2']-ymin,
                                                    ymax-ymin, hdREF['NAXIS2'])) + ext_pix
            ## Save hdREF
            self.hdr_ref = hdREF

            ## Test hdREF (Quick check: old=new or old<new)
            # w_new = fixwcs(header=hdREF).wcs
            # print('old: ', w.all_world2pix(
            #     self.hdr['CRVAL1'], self.hdr['CRVAL2'], 1))
            # print('new: ', w_new.all_world2pix(
            #     self.hdr['CRVAL1'], self.hdr['CRVAL2'], 1))
            # exit()
        else:
            raise ValueError('Cannot find reprojection reference! ')

    def make(self):
        '''
        Preparation (make header)
        '''
        flist = self.flist
        filREF = self.filREF
        hdREF = self.hdREF
        fmod = self.fmod
        ext_pix = self.ext_pix

        # if isinstance(flist, str):
        #     self.make_header(flist, filREF, hdREF, fmod, ext_pix)
        # elif isinstance(flist, list):
        self.make_header(flist[0], filREF, hdREF, fmod, ext_pix)
        if fmod=='ext':
            ## Refresh self.hdr_ref in every circle
            for f in flist:
                self.make_header(filIN=f, filREF=None,
                                 hdREF=self.hdr_ref, fmod='ext', ext_pix=ext_pix)
        
        tqdm.write('<imontage> Making ref header...[done]')

        return self.hdr_ref

    def footprint(self, filOUT=None):
        '''
        Save reprojection footprint
        '''
        if filOUT is None:
            filOUT = self.path_tmp+'footprint'
        
        Nx = self.hdr_ref['NAXIS1']
        Ny = self.hdr_ref['NAXIS2']
        im_fp = np.ones((Ny, Nx))
        
        comment = "<imontage> footprint"
        write_fits(filOUT, self.hdr_ref, im_fp, COMMENT=comment)

        return im_fp

    def reproject(self, filIN, filOUT=None,
                  dist=None, postfix=''):
        '''
        Reproject 2D image or 3D cube

        ------ INPUT ------
        filIN               single FITS file to reproject
        filOUT              output FITS file
        dist                uncertainty distribution
                              'norm' - N(0,1)
                              'splitnorm' - SN(0,lam,lam*tau)
        postfix              
        ------ OUTPUT ------

        '''
        super().__init__(filIN)
        
        if dist=='norm':
            self.rand_norm(filIN+'_unc')
        elif dist=='splitnorm':
            self.rand_splitnorm([filIN+'_unc_N', filIN+'_unc_P'])
        
        ## Set reprojection tmp path
        ##---------------------------
        filename = os.path.basename(filIN)
        rep_tmp = self.path_tmp+filename+postfix+'/'
        if not os.path.exists(rep_tmp):
            os.makedirs(rep_tmp)

        self.slist = self.slice(rep_tmp+'slice', '_') # gauss_noise inclu
        ## Do reprojection
        ##-----------------
        cube_rep = []
        # for k in range(self.Nw):
            # hdr = self.hdr.copy()
            # for kw in self.hdr.keys():
            #     if '3' in kw:
            #         del hdr[kw]
            # hdr['NAXIS'] = 2
            # phdu = fits.PrimaryHDU(header=hdr, data=self.im[k,:,:])
            # im_rep = reproject_interp(phdu, self.hdr_ref)[0]
        for s in self.slist:
            im_rep = reproject_interp(s+fitsext, self.hdr_ref)[0]
            cube_rep.append(im_rep)
            write_fits(s+'rep_', self.hdr_ref, im_rep)
            fclean(s+fitsext)
        self.im = np.array(cube_rep)

        comment = "Reprojected by <imontage>. "
        if filOUT is None:
            filOUT = self.path_tmp+filename+postfix+'_rep'
        self.file_rep = filOUT

        write_fits(filOUT, self.hdr_ref, self.im, self.wvl, wmod=0,
                   COMMENT=comment)
        
        return self.im

    def reproject_mc(self, filIN, filOUT=None, Nmc=0, dist=None):
        '''
        Generate Monte-Carlo uncertainties for reprojected input file
        '''
        dataset = type('', (), {})()

        hyperim = [] # [j,(w,)y,x]
        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Reprojection (MC level)'):

            if j==0:
                im0 = self.reproject(filIN, filOUT=filOUT, dist=dist)
                file_rep = self.file_rep
            else:
                hyperim.append(self.reproject(filIN, filOUT=filOUT,
                                              dist=dist, postfix='_'+str(j)))
        im0 = np.array(im0)
        hyperim = np.array(hyperim)
        unc = np.nanstd(hyperim, axis=0)
        comment = "Created by <imontage>"

        if Nmc>0:
            write_fits(file_rep+'_unc', self.hdr_ref, unc, self.wvl,
                       COMMENT=comment)

        dataset.im0 = im0
        dataset.unc = unc
        dataset.hyperim = hyperim

        return dataset

    def combine(self, flist, filOUT=None, method='avg',
                do_rep=True, Nmc=0, dist=None):
        '''
        Stitching input files (with the same wavelengths) to the ref WCS

        If Nmc==0, no MC
        '''
        flist = listize(flist)
        dataset = type('', (), {})()
        wvl = read_fits(flist[0]).wave
        dataset.wvl = wvl

        superim0 = [] # [i,(w,)y,x]
        superunc = [] # [i,(w,)y,x]
        superim = [] # [i,j,(w,)y,x]
        Nf = np.size(flist)
        for i in trange(Nf, leave=False,
                        desc='<imontage> Reprojection (file level)'):
            ## (Re)do reprojection
            ##---------------------
            if do_rep==True:
                ## With MC
                if Nmc>0:
                    rep = self.reproject_mc(flist[i], Nmc=Nmc, dist=dist)
                    im0 = rep.im0
                    
                    superunc.append(rep.unc)
                    superim.append(rep.hyperim)
                ## Without MC
                else:
                    im0 = self.reproject(flist[i])                    
                superim0.append(im0)

            ## Read archives
            ##---------------
            else:
                filename = os.path.basename(flist[i])
                file_rep = self.path_tmp+filename+'_rep'
                if Nmc>0:
                    hyperim = [] # [j,(w,)y,x]
                    for j in range(Nmc+1):
                        if j==0:
                            superunc.append(read_fits(file_rep+'_unc').data)
                        else:
                            file_rep = self.path_tmp+filename+'_'+str(j)+'_rep'
                            hyperim.append(read_fits(file_rep).data)
                    hyperim = np.array(hyperim)
                    superim.append(hyperim)
                superim0.append(read_fits(file_rep).data)

        superim0 = np.array(superim0)
        superunc = np.array(superunc)
        superim = np.array(superim)

        ## Combine images
        ##----------------
        hyperim_comb = []
        unc_comb = np.nanstd(hyperim_comb)
        ## Think about using 'try - except'
        if Nmc>0:
            inv_var = 1./superunc**2
            for j in trange(Nmc+1, leave=False,
                            desc='<imontage> Stitching'):
                if j==0:
                    if method=='avg':
                        im0_comb = nanavg(superim0, axis=0)
                    elif method=='wgt_avg':
                        im0_comb = nanavg(superim0, axis=0, weights=inv_var)
                else:
                    if method=='avg':
                        hyperim_comb.append(nanavg(superim[:,j-1], axis=0))
                    elif method=='wgt_avg':
                        hyperim_comb.append(
                            nanavg(superim[:,j-1], axis=0, weights=inv_var))
            hyperim_comb = np.array(hyperim_comb)
            unc_comb = np.nanstd(hyperim_comb)
        else:
            ## If no unc, inverse variance weighted mean not available
            im0_comb = nanavg(superim0, axis=0)

        if filOUT is not None:
            comment = "An <imontage> production"

            write_fits(filOUT, self.hdr_ref, im0_comb, wvl,
                       COMMENT=comment)
            write_fits(filOUT+'_unc', self.hdr_ref, unc_comb, wvl,
                       COMMENT=comment)
        
        dataset.im0_comb = im0_comb
        dataset.unc_comb = unc_comb
        dataset.hyperim_comb = hyperim_comb
        dataset.superim0 = superim0
        dataset.superunc = superunc
        dataset.superim = superim

        tqdm.write('<imontage> Combining images...[done]')
        
        return dataset

    def coadd(self, flist, filOUT=None,
              Nmc=0, dist=None):
        '''
        Same function with combine() using reproject.reproject_and_coadd()
        '''
        flist = listize(flist)
        dataset = type('', (), {})()
        comment = "Created by <imontage>"

        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Coadd [MC]'):
            sl = []
            if j==0:
                for f in flist:
                    super().__init__(f)

                    filename = os.path.basename(f)
                    coadd_tmp = self.path_tmp+filename+'/'
                    if not os.path.exists(coadd_tmp):
                        os.makedirs(coadd_tmp)
                    sl.append(self.slice(coadd_tmp+'slice', ext=fitsext))
                slist = [np.array(sl)]
            else:
                for f in flist:
                    super().__init__(f)

                    filename = os.path.basename(f)
                    coadd_tmp = self.path_tmp+filename+'/'
                    if not os.path.exists(coadd_tmp):
                        os.makedirs(coadd_tmp)
                    if dist=='norm':
                        self.rand_norm(f+'_unc')
                    elif dist=='splitnorm':
                        self.rand_splitnorm([f+'_unc_N', f+'_unc_P'])
                
                    sl.append(self.slice(coadd_tmp+'slice',
                                            postfix='_'+str(j), ext=fitsext))
                slist.append(np.array(sl))
        slist = np.array(slist)
        
        if self.Nw is None:
            Nw = 1
        else:
            Nw = self.Nw
        superim = []
        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Coadd [MC]'):
            if j==0:
                im = []
                for i in range(Nw):
                    im.append(reproject_and_coadd(slist[j,:,i], self.hdr_ref,
                                                  reproject_function=reproject_interp)[0])
                im = np.array(im)

                write_fits(filOUT, self.hdr_ref, im, self.wvl, wmod=0,
                           COMMENT=comment)
            else:
                hyperim = []
                for i in range(Nw):
                    hyperim.append(reproject_and_coadd(slist[j,:,i], self.hdr_ref,
                                                       reproject_function=reproject_interp)[0])
                superim.append(np.array(hyperim))
        superim = np.array(superim)
        unc = np.nanstd(superim, axis=0)

        if Nmc>0:
            write_fits(filOUT+'_unc', self.hdr_ref, unc, self.wvl, wmod=0,
                       COMMENT=comment)

        dataset.wvl = self.wvl
        dataset.im = im
        dataset.unc = unc
        dataset.superim = superim
        
        return dataset

    def clean(self, filIN=None):
        if filIN is not None:
            fclean(filIN)
        else:
            fclean(self.path_tmp)

class iswarp(improve):
    '''
    SWarp drop-in image montage toolkit
    i means <improve>-based
    Alternative to its fully Python-based twin <imontage>

    ------ INPUT ------
    flist               ref FITS files used to make header (footprint)
    refheader           scaling matrix adopted if co-exist with file
    center              center of output image frame
                          None - contains all input fields
                          str('hh:mm:ss,dd:mm:ss') - manual input RA,DEC
    pixscale            pixel scale (arcsec)
                          None - median of pixscale at center input frames
                          float() - in arcseconds
    verbose             default: True
    tmpdir              tmp file path
    ------ OUTPUT ------
    coadd.fits
    
    By default, SWarp reprojects all input to a WCS with diag CD matrix.
    "To implement the unusual output features required, 
     one must write a coadd.head ASCII file that contains 
     a custom anisotropic scaling matrix. "
    '''
    def __init__(self, flist=None, refheader=None,
                 center=None, pixscale=None, 
                 verbose=False, tmpdir=None):
        '''
        self: path_tmp, verbose
        (filIN, wmod, hdr, w, Ndim, Nx, Ny, Nw, im, wvl)
        '''
        if verbose==False:
            devnull = open(os.devnull, 'w')
        else:
            devnull = None
        self.verbose = verbose
        self.devnull = devnull
        
        ## Set path of tmp files
        if tmpdir is None:
            path_tmp = os.getcwd()+'/tmp_swp/'
        else:
            path_tmp = tmpdir
        if not os.path.exists(path_tmp):
            os.makedirs(path_tmp)
        
        self.path_tmp = path_tmp

        fclean(path_tmp+'coadd*') # remove previous coadd.fits/.head

        if flist is None:
            if refheader is None:
                raise InputError('<iswarp>','No input!')
            
            ## Define coadd frame via refheader
            else:
                if center is not None or pixscale is not None:
                    warnings.warn('The keywords center and pixscale are dumb. ')

                self.refheader = refheader
        else:
            ## Input files in list object
            flist = listize(flist)
                
            ## Images
            image_files = ' '
            list_ref = []
            for i in range(len(flist)):
                image = read_fits(flist[i]).data
                hdr = fixwcs(flist[i]+fitsext).header
                file_ref = flist[i]
                if image.ndim==3:
                    ## Extract 1st frame of the cube
                    file_ref = path_tmp+os.path.basename(flist[i])+'_ref'
                    write_fits(file_ref, hdr, image[0])
                
                image_files += file_ref+fitsext+' ' # SWarp input str
                list_ref.append(file_ref+fitsext) # reproject input

            ## Define coadd frame
            ##--------------------
            
            ## via SWarp without refheader (isotropic scaling matrix)
            
            ## Create config file
            SP.call('swarp -d > swarp.cfg',
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
            
            ## Config param list
            swarp_opt = ' -c swarp.cfg -SUBTRACT_BACK N -IMAGEOUT_NAME coadd.ref.fits '
            if center is not None:
                swarp_opt += ' -CENTER_TYPE MANUAL -CENTER '+center
            if pixscale is not None:
                swarp_opt += ' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE '+str(pixscale)
            if verbose=='quiet':
                swarp_opt += ' -VERBOSE_TYPE QUIET '
            
            ## Run SWarp
            SP.call('swarp '+swarp_opt+image_files,
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
            print('Running SWarp...')

            self.refheader = read_fits(path_tmp+'coadd.ref').header
            
            ## via reproject with refheader (custom anisotropic scaling matrix)
            if refheader is not None:
                if center is not None or pixscale is not None:
                    warnings.warn('The keywords center and pixscale are dumb. ')

                super().__init__(path_tmp+'coadd.ref')
                pix_old = [[0, 0]]
                pix_old.append([0, self.Ny])
                pix_old.append([self.Nx, 0])
                pix_old.append([self.Nx, self.Ny])
                world_arr = self.w.all_pix2world(np.array(pix_old), 1)
                
                w = fixwcs(header=refheader).wcs
                try:
                    pix_new = w.all_world2pix(world_arr, 1)
                except wcs.wcs.NoConvergence as e:
                    pix_new = e.best_solution
                    print("Best solution:\n{0}".format(e.best_solution))
                    print("Achieved accuracy:\n{0}".format(e.accuracy))
                    print("Number of iterations:\n{0}".format(e.niter))
                xmin = min(pix_new[:,0])
                xmax = max(pix_new[:,0])
                ymin = min(pix_new[:,1])
                ymax = max(pix_new[:,1])

                refheader['CRPIX1'] += -xmin
                refheader['CRPIX2'] += -ymin
                refheader['NAXIS1'] = math.ceil(xmax - xmin)
                refheader['NAXIS2'] = math.ceil(ymax - ymin)
                
                self.refheader = refheader

        # fclean(path_tmp+'*ref.fits')

    def footprint(self, filOUT=None):
        '''
        Save reprojection footprint
        '''
        if filOUT is None:
            filOUT = self.path_tmp+'footprint'
        
        Nx = self.refheader['NAXIS1']
        Ny = self.refheader['NAXIS2']
        im_fp = np.ones((Ny, Nx))
        
        comment = "<iswarp> footprint"
        write_fits(filOUT, self.refheader, im_fp, COMMENT=comment)

        return im_fp

    def combine(self, flist, combtype='med',
                keepedge=False, cropedge=False,
                dist=None, filOUT=None, tmpdir=None):
        '''
        Combine 

        ------ INPUT ------
        flist               input FITS files should have the same wvl
        combtype            combine type
                              med - median
                              avg - average
                              wgt_avg - inverse variance weighted average
        keepedge            default: False
        cropedge            crop the NaN edge of the frame
        dist                add uncertainties (filename+'_unc.fits' needed)
        filOUT              output FITS file
        ------ OUTPUT ------
        coadd.head          key for SWarp (inherit self.refheader)
        '''
        ds = type('', (), {})()
        
        verbose = self.verbose
        devnull = self.devnull
        path_tmp = self.path_tmp
        
        if tmpdir is None:
            path_comb = path_tmp+'comb/'
        else:
            path_comb = tmpdir
        if not os.path.exists(path_comb):
            os.makedirs(path_comb)

        ## Input files in list format
        flist = listize(flist)
        
        ## Header
        ##--------
        with open(path_tmp+'coadd.head', 'w') as f:
            f.write(str(self.refheader))

        ## Images and weights
        ##--------------------
        Nf = len(flist)
        
        imshape = read_fits(flist[0]).data.shape
        if len(imshape)==3:
            Nw = imshape[0]
            wvl = read_fits(flist[0]).wave
        else:
            Nw = 1
            wvl = None
        
        ## Build imlist & wgtlist (size=Nf)
        imlist = []
        wgtlist = []
        for i in range(Nf):
            filename = os.path.basename(flist[i])
            ## Set slice file
            file_slice = path_comb+filename
            
            ## Slice
            super().__init__(flist[i])
            if dist=='norm':
                self.rand_norm(flist[i]+'_unc')
            elif dist=='splitnorm':
                self.rand_splitnorm([flist[i]+'_unc_N', flist[i]+'_unc_P'])
            imlist.append(self.slice(file_slice, ''))
            
            if combtype=='wgt_avg':
                super().__init__(flist[i]+'_unc')
                wgtlist.append(self.slice_inv_sq(file_slice, '.weight'))

        ## Build image_files & weight_files (size=Nw)
        image_files = [' ']*Nw
        weight_files = [' ']*Nw

        ## Let's SWarp
        ##-------------
        hyperimage = []
        for k in trange(Nw, leave=False, 
            desc='<iswarp> Combining (by wvl)'):
            for i in range(Nf):
                image_files[k] += imlist[i][k]+fitsext+' '

                if combtype=='wgt_avg':
                    weight_files[k] += wgtlist[i][k]+fitsext+' '

            ## Create config file
            SP.call('swarp -d > swarp.cfg',
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
            ## Config param list
            swarp_opt = ' -c swarp.cfg -SUBTRACT_BACK N '
            if combtype=='med':
                pass
            elif combtype=='avg':
                swarp_opt += ' -COMBINE_TYPE AVERAGE '
            elif combtype=='wgt_avg':
                swarp_opt += ' -COMBINE_TYPE WEIGHTED '
                swarp_opt += ' -WEIGHT_TYPE MAP_WEIGHT '
                swarp_opt += ' -WEIGHT_SUFFIX .weight.fits '
                # swarp_opt += ' -WEIGHT_IMAGE '+weight_files[k] # not worked
            if verbose=='quiet':
                swarp_opt += ' -VERBOSE_TYPE QUIET '
            ## Run SWarp
            SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS3 '+image_files[k],
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
            coadd = read_fits(path_tmp+'coadd')
            newimage = coadd.data
            newheader = coadd.header

            ## Add back in the edges because LANCZOS3 kills the edges
            ## Do it in steps of less and less precision
            if keepedge==True:
                oldweight = read_fits(path_tmp+'coadd.weight').data
                if np.sum(oldweight==0)!=0:
                    SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS2 '+image_files[k],
                        shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                    edgeimage = read_fits(path_tmp+'coadd').data
                    newweight = read_fits(path_tmp+'coadd.weight').data
                    edgeidx = np.logical_and(oldweight==0, newweight!=0)
                    if edgeidx.any():
                        newimage[edgeidx] = edgeimage[edgeidx]

                    oldweight = read_fits(path_tmp+'coadd.weight').data
                    if np.sum(oldweight==0)!=0:
                        SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+image_files[k],
                            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                        edgeimage = read_fits(path_tmp+'coadd').data
                        newweight = read_fits(path_tmp+'coadd.weight').data
                        edgeidx = np.logical_and(oldweight==0, newweight!=0)
                        if edgeidx.any():
                            newimage[edgeidx] = edgeimage[edgeidx]

                        oldweight = read_fits(path_tmp+'coadd.weight').data
                        if np.sum(oldweight==0)!=0:
                            SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+image_files[k],
                                shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                            edgeimage = read_fits(path_tmp+'coadd').data
                            newweight = read_fits(path_tmp+'coadd.weight').data
                            edgeidx = np.logical_and(oldweight==0, newweight!=0)
                            if edgeidx.any():
                                newimage[edgeidx] = edgeimage[edgeidx]
            
            ## Astrometric flux-rescaling based on the local ratio of pixel scale
            ## Complementary for lack of FITS kw 'FLXSCALE'
            ## Because SWarp is conserving surface brightness/pixel
            oldcdelt = get_pc(wcs=fixwcs(flist[i]+fitsext).wcs).cdelt
            newcdelt = get_pc(wcs=fixwcs(path_tmp+'coadd'+fitsext).wcs).cdelt
            old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
            new_pixel_fov = abs(newcdelt[0]*newcdelt[1])
            newimage = newimage * old_pixel_fov/new_pixel_fov
            newimage[newimage==0] = np.nan
            # write_fits(path_comb+'coadd_'+str(k), newheader, newimage)
            # tqdm.write(str(old_pixel_fov))
            # tqdm.write(str(new_pixel_fov))
            # tqdm.write(str(abs(newheader['CD1_1']*newheader['CD2_2'])))

            if Nw==1:
                hyperimage = newimage
            else:
                hyperimage.append(newimage)

        hyperimage = np.array(hyperimage)

        if cropedge:
            reframe = improve(header=newheader, image=hyperimage, wave=wvl)
            xlist = []
            for x in range(reframe.Nx):
                if not np.isnan(reframe.im[:,:,x]).all():
                    xlist.append(x)
            ylist = []
            for y in range(reframe.Ny):
                if not np.isnan(reframe.im[:,y,:]).all():
                    ylist.append(y)
            xmin = min(xlist)
            xmax = max(xlist)+1
            ymin = min(ylist)
            ymax = max(ylist)+1
            dx = xmax-xmin
            dy = ymax-ymin
            x0 = xmin+dx/2
            y0 = ymin+dy/2

            reframe.crop(filOUT=path_tmp+'coadd.ref',
                         sizpix=(dx,dy), cenpix=(x0,y0))
            newheader = reframe.hdr
            hyperimage = reframe.im
            
        if filOUT is not None:
            write_fits(filOUT, newheader, hyperimage, wvl)

        ds.header = newheader
        ds.image = hyperimage
        ds.wvl = wvl

        return ds

    def reproject_mc(self, filIN, Nmc=0,
                     combtype='med', keepedge=False, cropedge=False,
                     dist=None, filOUT=None, tmpdir=None):
        '''
        Generate Monte-Carlo uncertainties for reprojected input file
        '''
        ds = type('', (), {})()

        hyperim = [] # [j,(w,)y,x]
        for j in trange(Nmc+1, leave=False,
                        desc='<iswarp> Reprojection (MC level)'):

            if j==0:
                comb = self.combine(filIN, filOUT=filOUT, tmpdir=tmpdir,
                                    keepedge=keepedge, cropedge=cropedge,
                                    combtype=combtype)
                im0 = comb.image
            else:
                hyperim.append( self.combine(filIN, filOUT=filOUT+'_'+str(j),
                                             dist=dist, tmpdir=tmpdir,
                                             keepedge=keepedge, cropedge=cropedge,
                                             combtype=combtype).image )
        im0 = np.array(im0)
        hyperim = np.array(hyperim)
        unc = np.nanstd(hyperim, axis=0)
        comment = "Created by <iswarp>"

        if Nmc>0:
            write_fits(filOUT+'_unc', comb.header, unc, comb.wvl,
                       COMMENT=comment)

        ds.im0 = im0
        ds.unc = unc
        ds.hyperim = hyperim

        return ds
    
    def clean(self, filIN=None):
        if filIN is not None:
            fclean(filIN)
        else:
            fclean(self.path_tmp)

class iconvolve(improve):
    '''
    Convolve 2D image or 3D cube with given kernels
    i means <improve>-based or IDL-based

    ------ INPUT ------
    filIN               input FITS file
    kfile               convolution kernel(s) (tuple or list)
    klist               CSV file storing kernel names
    filUNC              unc file (add gaussian noise)
    psf                 PSF list
    convdir             do_conv path (Default: None -> filIN path)
    filOUT              output file
    ------ OUTPUT ------
    '''
    def __init__(self, filIN, kfile, klist,
                 filUNC=None, dist=None, psf=None, convdir=None, filOUT=None):
        ## INPUTS
        super().__init__(filIN)
        
        if dist=='norm':
            self.rand_norm(filUNC)
        elif dist=='splitnorm':
            self.rand_splitnorm(filUNC)

        ## Input kernel file in list format
        self.kfile = listize(kfile)

        ## doc (csv) file of kernel list
        self.klist = klist
        self.path_conv = convdir
        self.filOUT = filOUT

        ## INIT
        if psf is None:
            self.psf = [1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.]
        else:
            self.psf = psf
        self.sigma_lam = None
        
    def spitzer_irs(self):
        '''
        Spitzer/IRS PSF profil
        [REF]
        Pereira-Santaella, Miguel, Almudena Alonso-Herrero, George H.
        Rieke, Luis Colina, Tanio Díaz-Santos, J.-D. T. Smith, Pablo G.
        Pérez-González, and Charles W. Engelbracht. “Local Luminous
        Infrared Galaxies. I. Spatially Resolved Observations with the
        Spitzer Infrared Spectrograph.” The Astrophysical Journal
        Supplement Series 188, no. 2 (June 1, 2010): 447.
        doi:10.1088/0067-0049/188/2/447.
        https://iopscience.iop.org/article/10.1088/0067-0049/188/2/447/pdf
        '''
        sim_par_wave = [0, 13.25, 40.]
        sim_par_fwhm = [2.8, 3.26, 10.1]
        sim_per_wave = [0, 15.5, 40.]
        sim_per_fwhm = [3.8, 3.8, 10.1]
        
        ## fwhm (arcsec)
        fwhm_par = np.interp(self.wvl, sim_par_wave, sim_par_fwhm)
        fwhm_per = np.interp(self.wvl, sim_per_wave, sim_per_fwhm)
        #fwhm_lam = np.sqrt(fwhm_par * fwhm_per)
        
        ## sigma (arcsec)
        sigma_par = fwhm_par / (2. * np.sqrt(2.*np.log(2.)))
        sigma_per = fwhm_per / (2. * np.sqrt(2.*np.log(2.)))
        self.sigma_lam = np.sqrt(sigma_par * sigma_per)
        
    # def choker(self, flist):
    #     '''
    #     ------ INPUT ------
    #     flist               FITS files to be convolved
    #     ------ OUTPUT ------
    #     '''
    #     ## Input files in list format
    #     flist = listize(flist)
        
    #     ## CHOose KERnel(s)
    #     lst = []
    #     for i, image in enumerate(flist):
    #         ## check PSF profil (or is not a cube)
    #         if self.sigma_lam is not None:
    #             image = flist[i]
    #             ind = closest(self.psf, self.sigma_lam[i])
    #             kernel = self.kfile[ind]
    #         else:
    #             image = flist[0]
    #             kernel = self.kfile[0]
    #         ## lst line elements: image, kernel
    #         k = [image, kernel]
    #         lst.append(k)

    #     ## write csv file
    #     write_csv(self.klist, header=['Images', 'Kernels'], dset=lst)

    def choker(self, flist):
        '''
        ------ INPUT ------
        flist               FITS files to be convolved
        ------ OUTPUT ------
        '''
        ## Input files in list format
        flist = listize(flist)
        
        ## CHOose KERnel(s)
        image = []
        kernel = []
        for i, filim in enumerate(flist):
            ## check PSF profil (or is not a cube)
            if self.sigma_lam is not None:
                image.append(filim)
                ind = closest(self.psf, self.sigma_lam[i])
                kernel.append(self.kfile[ind])
            else:
                image.append(flist[0])
                kernel.append(self.kfile[0])

        ## write csv file
        dataset = Table([image, kernel], names=['Images', 'Kernels'])
        ascii.write(dataset, self.klist+csvext, format='csv')

    def do_conv(self, idldir, verbose=False):
        '''
        ------ INPUT ------
        idldir              path of IDL routines
        ------ OUTPUT ------
        '''
        if verbose==False:
            devnull = open(os.devnull, 'w')
        else:
            devnull = None

        filename = os.path.basename(self.filIN)

        if self.Ndim==3:
            if self.path_conv is not None:
                f2conv = self.slice(self.path_conv+filename) # gauss_noise inclu
            else:
                f2conv = self.slice(self.filIN) # gauss_noise inclu
            
            self.spitzer_irs()

        else:
            f2conv = [self.filIN]
        
        self.choker(f2conv)

        SP.call('idl conv.pro',
                shell=True, cwd=idldir, stdout=devnull, stderr=SP.STDOUT)

        ## OUTPUTS
        ##---------
        if self.Ndim==3:
            im = []
            self.slist = []
            for f in f2conv:
                im.append(read_fits(f+'_conv').data)
                self.slist.append(f+'_conv')

            self.convim = np.array(im)
            ## recover 3D header cause the lost of WCS due to PS3_0='WCS-TAB'
            # self.hdr = read_fits(self.filIN).header

            fclean(f+'_conv'+fitsext)
        else:
            self.convim = read_fits(self.filIN+'_conv').data

            fclean(self.filIN+'_conv'+fitsext)
        
        if self.filOUT is not None:
            comment = "Convolved by G. Aniano's IDL routine."
            write_fits(self.filOUT, self.hdr, self.convim, self.wvl, 
                COMMENT=comment)

    def image(self):
        return self.convim

    def wave(self):
        return self.wvl

    def filenames(self):
        return self.slist

    def clean(self, filIN=None):
        if filIN is not None:
            fclean(filIN)
        else:
            if self.path_conv is not None:
                fclean(self.path_conv)
        
class cupid(improve):
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
                        
    def spec_build(self, filOUT=None, write_unc=False, dist=None,
                   Nx=None, Ny=32, Nsub=1, pixscale=None,
                   wmin=None, wmax=None, tmpdir=None, fiLOG=None,
                   sig_pt=0, fill='med', supix=False, swarp=False):
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
        fill                fill value of no data regions after shift
                              'med': axis median (default)
                              'avg': axis average
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
        iwi = closest(wave, wmin, side='right')
        iws = closest(wave, wmax, side='left')+1
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
        unc = np.repeat(unc, Nx, axis=2)
        unc_N = np.repeat(unc_N, Nx, axis=2)
        unc_P = np.repeat(unc_P, Nx, axis=2)
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
            self.rand_pointing(sig_pt, fill=fill, tmpdir=tmpdir,
                               xscale=self.xscale, yscale=self.yscale, swarp=False)
        else:
            self.rand_pointing(sig_pt, fill=fill, tmpdir=tmpdir,
                               xscale=1, yscale=1, swarp=swarp)
                
        ## Update self variables (for next steps)
        self.reinit(header=self.hdr, image=self.im, wave=self.wvl,
                    wmod=self.wmod, verbose=self.verbose)

        if filOUT is not None:
            comment = "Assembled AKARI/IRC slit spectroscopy cube. "
            write_fits(filOUT, self.hdr, self.im, self.wvl,
                       COMMENT=comment)

            if write_unc==True:
                uncom = "Assembled AKARI/IRC slit spec uncertainty cube. "
                write_fits(filOUT+'_unc', self.hdr, unc, self.wvl,
                           COMMENT=uncom)

                uncom_N = "Assembled AKARI/IRC slit spec uncertainty (N) cube. "
                write_fits(filOUT+'_unc_N', self.hdr, unc_N, self.wvl,
                           COMMENT=uncom)

                uncom_P = "Assembled AKARI/IRC slit spec uncertainty (P) cube. "
                write_fits(filOUT+'_unc_P', self.hdr, unc_P, self.wvl,
                           COMMENT=uncom)

        if fiLOG is not None:
            write_hdf5(fiLOG, 'Observation ID', [self.obsid])
            write_hdf5(fiLOG, 'NIR Slit', [self.slit], append=True)
            write_hdf5(fiLOG, 'Spectral disperser', [self.spec], append=True)
            write_hdf5(fiLOG, 'N3 image',[self.imref], append=True)
            write_hdf5(fiLOG, 'Pointing accuracy', sig_pt, append=True)
            write_hdf5(fiLOG, 'Spectral sampling size', self.Nw, append=True)
            write_hdf5(fiLOG, 'Slit length', self.Ny, append=True)
            write_hdf5(fiLOG, 'Slit width', self.Nx, append=True)
            write_hdf5(fiLOG, 'Subslit number', Nsub, append=True)
            write_hdf5(fiLOG, 'Pixel size', pixscale, append=True)
            write_hdf5(fiLOG, 'Super pixel size',
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

def wmask(filIN, filOUT=None):
    '''
    MASK Wavelengths

    --- INPUT ---
    filIN       input fits file 
    filOUT      overwrite fits file (Default: NO)
    --- OUTPUT ---
    data_new    new fits data
    wave_new    new fits wave
    '''
    pass

def wclean(filIN, cmod='eq', cfile=None,
           wmod=0, filOUT=None, verbose=False):
    '''
    CLEAN Wavelengths, alternative to the wrange option of concatenate()

    --- INPUT ---
    filIN       input fits file
    wmod        wave mode
    cmod        clean mode (Default: 'eq')
    cfile       input csv file (archived info)
    filOUT      overwrite fits file (Default: NO)
    verbose     display wclean info (Default: False)
    --- OUTPUT ---
    data_new    new fits data
    wave_new    new fits wave
    '''
    ds = read_fits(filIN)
    hdr = ds.header
    data = ds.data
    wave = ds.wave
    Nw = len(wave)
    
    ind = [] # list of indices of wvl to remove
    if cfile is not None:
        # indarxiv = read_csv(cfile, 'Ind')[0]
        indarxiv = ascii.read(cfile+csvext)['Ind']
        ind = []
        for i in indarxiv:
            ind.append(int(i))
    else:
        ## Detect crossing wvl
        ##---------------------
        for i in range(Nw-1):
            if wave[i]>=wave[i+1]: # found wave(i+1), i_max=Nw-2
                
                wmin = -1 # lower limit: closest wave smaller than wave[i+1]
                wmax = 0 # upper limit: closest wave larger than wave[i]
                
                for j in range(i+1):
                    dw = wave[i+1] - wave[i-j]
                    if dw>0: # found the closest smaller wave[i-j]
                        wmin = i-j
                        break # only the innermost loop
                if wmin==-1:
                    warnings.warn('Left side fully covered! ')
                
                for j in range(Nw-i-1):
                    dw = wave[i+1+j] - wave[i]
                    if dw>0: # found the closest larger wave[i+1+j]
                        wmax = i+1+j
                        break
                if wmax==0:
                    warnings.warn('Right side fully covered! ')

                Nw_seg = wmax-wmin-1 # number of crossing wvl in segment
                wave_seg = [] # a segment (every detect) of wave
                ind_seg = [] # corresponing segment for sort use
                for k in range(Nw_seg):
                    wave_seg.append(wave[wmin+1+k])
                    ind_seg.append(wmin+1+k)
                ## index list of sorted wave_seg
                ilist = sorted(range(len(wave_seg)), key=wave_seg.__getitem__)
                ## index of wave_seg center
                icen = math.floor((Nw_seg-1)/2)

                ## Visualisation (for test use)
                ##------------------------------
                # print('wave, i: ', wave[i], i)
                # print('wave_seg: ', wave_seg)
                # print('ind_seg: ', ind_seg)
                # print('ilist: ', ilist)
                # print('icen: ', icen)

                ## Remove all crossing wvl between two channels
                ##----------------------------------------------
                if cmod=='all': # most conservative but risk having holes
                    pass
                ## Remove (almost) equal wvl (NOT nb of wvl!) for both sides
                ##-----------------------------------------------------------
                elif cmod=='eq': # (default)
                    ## Select ascendant pair closest to segment center
                    for k in range(icen):
                        if ilist[icen]>ilist[0]: # large center
                            if ilist[icen-k]<ilist[0]:
                                for p in range(ilist[icen-k]+1):
                                    del ind_seg[0]
                                for q in range(Nw_seg-ilist[icen]):
                                    del ind_seg[-1]
                                break
                        else: # small center
                            if ilist[icen+k]>ilist[0]:
                                for p in range(ilist[icen]+1):
                                    del ind_seg[0]
                                for q in range(Nw_seg-ilist[icen+k]):
                                    del ind_seg[-1]
                                break
                ## Leave 2 closest wvl not crossing
                ##----------------------------------
                elif cmod=='closest_left':
                    for k in range(ilist[0]):
                        del ind_seg[0]
                elif cmod=='closest_right':
                    for k in range(Nw_seg-ilist[0]):
                        del ind_seg[-1]
                ## Others
                ##--------
                else:
                    raise ValueError('Non-supported clean mode! ')

                # print('ind_seg (final): ', ind_seg)
                ind.extend(ind_seg)

    ## Do clean
    ##----------
    data_new = np.delete(data, ind, axis=0)
    wave_new = list(np.delete(np.array(wave), ind))

    ## Display clean detail
    ##----------------------
    if verbose==True:
        print('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print('Number of wavelengths deleted: ', len(ind))
        print('Ind, wavelengths: ')
        for i in ind:
            print(i, wave[i])
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n')

    ## Overwrite fits file
    ##---------------------
    if filOUT is not None:
        # comment = 'Wavelength removal info in _wclean_info.csv'
        write_fits(filOUT, hdr, data_new, wave_new, wmod) # hdr auto changed
        
        ## Write csv file
        wlist = []
        for i in ind:
            wlist.append([i, wave[i]])
        write_csv(filOUT+'_wclean_info',
                  header=['Ind', 'Wavelengths'], dset=wlist)

    return data_new, wave_new

def interfill(arr, axis):
    '''
    FILL undersampling/artificial gap by (bspl)INTERpolation

    --- INPUT ---
    arr         array
    axis        axis along which interpolation
    --- OUTPUT ---
    newarr      new array
    '''
    print(">> fill gaps with b-splines <<")

    axsh = arr.shape
    NAXIS = np.size(axsh)
    newarr = np.copy(arr)
    if NAXIS==1: # 1D array
        x = np.arange(axsh[0])
        for i in range(axsh[0]):
            newarr = bsplinterpol(x, arr, x)
    if NAXIS==2: # no wavelength
        if axis==0: # col direction
            y = np.arange(axsh[0])
            for i in range(axsh[1]):
                col = bsplinterpol(y, arr[:,i], y)
                for j in range(axsh[0]):
                    newarr[j,i] = col[j]
        elif axis==1: # row direction
            x = np.arange(axsh[1])
            for j in range(axsh[0]):
                row = bsplinterpol(x, arr[j,:], x)
                for i in range(axsh[1]):
                    newarr[j,i] = row[i]
        else:
            raise ValueError('Unknown axis! ')
    elif NAXIS==3:
        if axis==0: # fill wavelength
            z = np.arange(axsh[0])
            for i in range(axsh[2]):
                for j in range(axsh[1]):
                    wvl = bsplinterpol(z, arr[:,j,i], z)
                    for k in range(axsh[0]):
                        newarr[k,j,i] = wvl[k]
        elif axis==1: # col direction
            y = np.arange(axsh[1])
            for k in range(axsh[0]):
                for i in range(axsh[2]):
                    col = bsplinterpol(y, arr[k,:,i], y)
                    for j in range(axsh[1]):
                        newarr[k,j,i] = col[j]
        elif axis==2: # row direction
            x = np.arange(axsh[2])
            for k in range(axsh[0]):
                for j in range(axsh[1]):
                    row = bsplinterpol(x, arr[k,j,:], x)
                    for i in range(axsh[2]):
                        newarr[k,j,i] = row[i]
        else:
            raise ValueError('Unknown axis! ')
    else:
        raise ValueError('Non-supported array shape! ')

    return newarr

def hextract(filIN, filOUT, x0, x1, y0, y1):
    '''
    Crop 2D image with pixel sequence numbers
    [REF] IDL lib hextract
    https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/hextract.pro
    '''
    ds = read_fits(filIN)
    oldimage = ds.data
    hdr = ds.header
    # hdr['NAXIS1'] = x1 - x0 + 1
    # hdr['NAXIS2'] = y1 - y0 + 1
    hdr['CRPIX1'] += -x0
    hdr['CRPIX2'] += -y0
    newimage = oldimage[y0:y1+1, x0:x1+1]

    write_fits(filOUT, hdr, newimage)

    return newimage

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
    if tmpdir is None:
        path_tmp = os.getcwd()+'/tmp_hswarp/'
    else:
        path_tmp = tmpdir
    if not os.path.exists(path_tmp):
        os.makedirs(path_tmp)

    fclean(path_tmp+'coadd*')
    ## Make input
    write_fits(path_tmp+'old', oldheader, oldimage)
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
    coadd = read_fits(path_tmp+'coadd')
    newimage = coadd.data
    newheader = coadd.header

    ## Add back in the edges because LANCZOS3 kills the edges
    ## Do it in steps of less and less precision
    if keepedge==True:
        oldweight = read_fits(path_tmp+'coadd.weight').data
        if np.sum(oldweight==0)!=0:
            SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE LANCZOS2 '+' old.fits',
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
            edgeimage = read_fits(path_tmp+'coadd').data
            newweight = read_fits(path_tmp+'coadd.weight').data
            edgeidx = np.logical_and(oldweight==0, newweight!=0)
            if edgeidx.any():
                newimage[edgeidx] = edgeimage[edgeidx]

            oldweight = read_fits(path_tmp+'coadd.weight').data
            if np.sum(oldweight==0)!=0:
                SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+' old.fits',
                        shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                edgeimage = read_fits(path_tmp+'coadd').data
                newweight = read_fits(path_tmp+'coadd.weight').data
                edgeidx = np.logical_and(oldweight==0, newweight!=0)
                if edgeidx.any():
                    newimage[edgeidx] = edgeimage[edgeidx]

                oldweight = read_fits(path_tmp+'coadd.weight').data
                if np.sum(oldweight==0)!=0:
                    SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+' old.fits',
                            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                    edgeimage = read_fits(path_tmp+'coadd').data
                    newweight = read_fits(path_tmp+'coadd.weight').data
                    edgeidx = np.logical_and(oldweight==0, newweight!=0)
                    if edgeidx.any():
                        newimage[edgeidx] = edgeimage[edgeidx]

    ## SWarp is conserving surface brightness/pixel
    ## while the pixels size changes
    oldcdelt = get_pc(wcs=fixwcs(header=oldheader).wcs).cdelt
    refcdelt = get_pc(wcs=fixwcs(header=refheader).wcs).cdelt
    old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
    new_pixel_fov = abs(refcdelt[0]*refcdelt[1])
    newimage = newimage * old_pixel_fov/new_pixel_fov
    newimage[newimage==0] = np.nan
    # print('-------------------')
    # print(old_pixel_fov/new_pixel_fov)
    write_fits(path_tmp+'new', newheader, newimage)
    # print('-------------------')
    
    ## Delete tmp file if tmpdir not given
    if tmpdir is None:
        fclean(path_tmp)

    ds.image = newimage
    ds.header = newheader

    return ds

def concatenate(flist, filOUT=None, comment=None,
                wsort=False, wrange=None,
                keepfrag=True, cropedge=False):
    '''
    wsort=True can be used with wclean
    When wsort=False, wrange is used to avoid wavelength overlapping

    '''
    dataset = type('', (), {})()

    if wrange is None:
        wrange = [ (2.50, 5.00), # irc
                   (5.21, 7.56), # sl2
                   (7.57, 14.28), # sl1
                   (14.29, 20.66), # ll2
                   (20.67, 38.00), ] # ll1
    wmin = []
    wmax = []
    for i in range(len(wrange)):
        wmin.append(wrange[i][0])
        wmax.append(wrange[i][1])
    
    ## Read data
    wave = []
    data = []

    ## Keep all wavelengths and sort them in ascending order
    if wsort==True:
        for f in flist:
            ds = read_fits(f)
            data.append(ds.data)
            wave.append(ds.wave)
    ## Keep wavelengths in the given ranges (wrange)
    else:
        for f in flist:
            ds = read_fits(f)
            imin = closest(wmin, ds.wave[0])
            imax = closest(wmax, ds.wave[-1])
            iwi = 0
            iws = -1
            for i, w in enumerate(ds.wave[:-2]):
                if w<wmin[imin] and ds.wave[i+1]>wmin[imin]:
                    iwi = i+1
                if w<wmax[imax] and ds.wave[i+1]>wmax[imax]:
                    iws = i+1
            data.append(ds.data[iwi:iws])
            wave.append(ds.wave[iwi:iws])

    data = np.concatenate(data, axis=0)
    wave = np.concatenate(wave)
    hdr = ds.header
    ## Sort
    ind = sorted(range(len(wave)), key=wave.__getitem__)
    # wave = np.sort(wave)
    wave = wave[ind]
    data = data[ind]
    ## NaN mask
    if not keepfrag:
        ma_any = np.isnan(data).any(axis=0)
        for k in range(len(wave)):
            data[k][ma_any] = np.nan

    if cropedge:
        reframe = improve(header=hdr, image=data, wave=wave)
        xlist = []
        for x in range(reframe.Nx):
            if not np.isnan(reframe.im[:,:,x]).all():
                xlist.append(x)
        ylist = []
        for y in range(reframe.Ny):
            if not np.isnan(reframe.im[:,y,:]).all():
                ylist.append(y)
        xmin = min(xlist)
        xmax = max(xlist)+1
        ymin = min(ylist)
        ymax = max(ylist)+1
        dx = xmax-xmin
        dy = ymax-ymin
        x0 = xmin+dx/2
        y0 = ymin+dy/2

        reframe.crop(sizpix=(dx,dy), cenpix=(x0,y0))
        data = reframe.im
        hdr = reframe.hdr

    dataset.wave = wave
    dataset.data = data
    
    ## Write FITS file
    if filOUT is not None:
        write_fits(filOUT, hdr, data, wave, COMMENT=comment)

    return dataset

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

    pass

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Imaging

    improve:
        uncert, rand_norm, rand_splitnorm, 
        slice, slice_inv_sq, crop, rebin
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
    imontage(improve):
        make_header, make, footprint, reproject, reproject_mc,
        combine, coadd, clean
    iswarp(improve):
        footprint, combine, clean
    iconvolve(improve):
        spitzer_irs, choker, do_conv, image, wave,
        filenames, clean
    respect(improve):
        concat, smooth, mask
    sextract(improve):
        rand_pointing, spec_build, sav_build,
        header, image, wave
    wmask, wclean, interfill, hextract, hswarp, concatenante, 

"""

from tqdm import tqdm, trange
import os
import math
import numpy as np
from scipy.io import readsav
from astropy import wcs
from astropy.io import ascii
from astropy.table import Table
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
import subprocess as SP

## Local
from utilities import InputError
from inout import (fitsext, csvext, ascext, fclean,
                   read_fits, write_fits, savext
                   # read_csv, write_csv, read_ascii,
)
from arrays import allist, closest
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
    def __init__(self, filIN, wmod=0, verbose=False):
        '''
        self: filIN, wmod, hdr, w, cdelt, pc, cd, Ndim, Nx, Ny, Nw, im, wvl
        '''
        
        ## INPUTS
        self.filIN = filIN
        self.wmod = wmod
        self.verbose = verbose

        ## read image/cube
        ## self.hdr is a 2D (reduced) header
        ws = fixwcs(filIN+fitsext)
        self.hdr = ws.header
        self.w = ws.wcs
        pcdelt = get_pc(wcs=ws.wcs)
        self.cdelt = pcdelt.cdelt
        self.pc = pcdelt.pc
        self.cd = pcdelt.cd
        self.Nx = self.hdr['NAXIS1']
        self.Ny = self.hdr['NAXIS2']
        self.Nw = None
        if verbose==True:
            print('<improve> file: ', filIN)
            print('Raw size (pix): {} * {}'.format(self.Nx, self.Ny))
        ## 3D cube slicing
        ds = read_fits(filIN)
        self.im = ds.data
        self.wvl = ds.wave
        self.Ndim = self.im.ndim
        if self.Ndim==3:
            if self.im.shape[0]==1:
                self.Ndim = 2 # Nw=1 patch
            else:
                self.Nw = len(self.wvl)

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
        ## 
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

    def rand_norm(self, filIN=None, unc=None, sigma=1., mu=0.):
        '''
        Add random N(0,1) noise
        '''
        if filIN is not None:
            unc = read_fits(filIN).data

        if unc is not None:
            ## unc should have the same dimension with im
            theta = np.random.normal(mu, sigma, self.im.shape)
            self.im += theta * unc

        return self.im

    def rand_splitnorm(self, filIN=None, unc=None, sigma=1., mu=0.):
        '''
        Add random SN(0,lam,lam*tau) noise

        ------ INPUT ------
        filIN               2 FITS files for unc of left & right sides
        unc                 2 uncertainty ndarrays
        ------ OUTPUT ------
        '''
        if filIN is not None:
            unc = []
            for f in filIN:
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

    def slice(self, filSL, postfix='', ext=''):
        ## 3D cube slicing
        slist = []
        if self.Ndim==3:
            hdr = self.hdr.copy()
            for kw in self.hdr.keys():
                if '3' in kw:
                    del hdr[kw]
            hdr['NAXIS'] = 2
            for k in range(self.Nw):
                ## output filename list
                f = filSL+'_'+'0'*(4-len(str(k)))+str(k)+postfix
                slist.append(f+ext)
                write_fits(f, hdr, self.im[k,:,:]) # gauss_noise inclu
        else:
            f = filSL+'_0000'+postfix
            slist.append(f+ext)
            write_fits(f, self.hdr, self.im) # gauss_noise inclu
            if self.verbose==True:
                print('Input file is a 2D image which cannot be sliced! ')
                print('Rewritten with only random noise added (if provided).')

        return slist

    def slice_inv_sq(self, filSL, postfix=''):
        ## Inversed square cube slicing
        inv_sq = 1./self.im**2
        slist = []
        if self.Ndim==3:
            hdr = self.hdr.copy()
            for kw in self.hdr.keys():
                if '3' in kw:
                    del hdr[kw]
            hdr['NAXIS'] = 2
            for k in range(self.Nw):
                ## output filename list
                f = filSL+'_'+'0'*(4-len(str(k)))+str(k)+postfix
                slist.append(f)
                write_fits(f, hdr, inv_sq[k,:,:]) # gauss_noise inclu
        else:
            f = filSL+'_0000'+postfix
            slist.append(f)
            write_fits(f, self.hdr, inv_sq) # gauss_noise inclu

        return slist
    
    def crop(self, filOUT=None,
             sizpix=None, cenpix=None, sizval=None, cenval=None):
        '''
        If pix and val co-exist, pix will be taken.
        Note that some self variables such as cdelt are
        not renewed...idem. for rebin

        ------ INPUT ------
        filOUT              output file
        sizpix              crop size in pix (dx, dy)
        cenpix              crop center in pix (x, y)
        sizval              crop size in deg (dRA, dDEC) -> (dx, dy)
        cenval              crop center in deg (RA, DEC) -> (x, y)
        ------ OUTPUT ------
        self.im             cropped image array
        '''
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
        if not (0<cenpix[0]<self.Nx and 0<cenpix[1]<self.Ny):
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
            self.im = self.im[:, ymin:ymax, xmin:xmax] # gauss_noise inclu
            ## recover 3D non-reduced header
            self.hdr = read_fits(self.filIN).header
        elif self.Ndim==2:
            self.im = self.im[ymin:ymax, xmin:xmax] # gauss_noise inclu
        ## Modify header
        ## Suppose no non-linear distortion
        self.hdr['CRPIX1'] = sizpix[0] / 2.
        self.hdr['CRPIX2'] = sizpix[1] / 2.
        self.hdr['CRVAL1'] = cenval[0]
        self.hdr['CRVAL2'] = cenval[1]
        ## Write cropped image/cube
        if filOUT is not None:
            # comment = "[ICROP]ped at centre: [{:.8}, {:.8}]. ".format(*cenval)
            # comment = "with size [{}, {}] (pix).".format(*sizpix)
            write_fits(filOUT, self.hdr, self.im, self.wvl, self.wmod)

        return self.im

    def rebin(self, filOUT=None, pixscale=None):
        '''
        Shrink or expand the size of an array an arbitrary amount 
        using interpolation
        Note that the actual version collimates on zero point,
        it does not conserve flux on upper right edges
        if not integer multiplication!
        [REF] IDL lib hrebin/frebin
        https://idlastro.gsfc.nasa.gov/ftp/pro/astrom/hrebin.pro
        https://github.com/wlandsman/IDLAstro/blob/master/pro/frebin.pro

        ------ INPUT ------
        filOUT              output file
        pixscale            output pixel scale in arcsec
        ------ OUTPUT ------
        newimage            rebinned image array
        '''
        oldimage = self.im
        hdr = self.hdr
        header = hdr.copy()
        oldw = self.w
        # cd = w.pixel_scale_matrix
        oldcd = self.cd
        cdelt = self.cdelt
        
        if pixscale is not None:
            pixscale = pixscale / 3600. # convert arcsec to degree
            xratio = math.ceil(pixscale / abs(cdelt[0])) # Expansion or contraction in X
            yratio = math.ceil(pixscale / abs(cdelt[1])) # Expansion or contraction in Y
        else:
            xratio = 1
            yratio = 1
        if self.verbose==True:
            print('----------')
            print('The actual (rounded) pixel scale of x axis is {} arcsec'.format(
                abs(cdelt[0])*xratio*3600))
            print('The actual (rounded) pixel scale of y axis is {} arcsec'.format(
                abs(cdelt[1])*yratio*3600))
            print('----------')
            
        # lam = yratio/xratio
        # pix_ratio = xratio*yratio
        Nx = math.ceil(self.Nx / xratio)
        Ny = math.ceil(self.Ny / yratio)
    
        ## Rebin
        ##-------
        if xratio<1:
            newimage = np.repeat(oldimage, Nx/self.Nx, axis=2)
        else:
            x_edge = Nx*xratio-self.Nx
            # x_after = math.ceil(x_edge / 2)
            # x_before = x_edge - x_after
            if self.Ndim==3:
                # npad is a tuple of (n_before, n_after) for each dimension
                xpad = ((0,0), (0,0), (0, x_edge))
                newimage = np.pad(oldimage, pad_width=xpad,
                                  mode='constant', constant_values=np.nan)
                newimage = newimage.reshape((self.Nw,self.Ny,Nx,xratio)).mean(3)
            else:
                xpad = ((0,0), (0, x_edge))
                newimage = np.pad(oldimage, pad_width=xpad,
                                  mode='constant', constant_values=np.nan)
                newimage = newimage.reshape((self.Ny,Nx,xratio)).mean(2)
        if yratio<1:
            newimage = np.repeat(newimage, Ny/self.Ny, axis=1)
        else:
            y_edge = Ny*yratio-self.Ny
            # y_after = math.ceil(y_edge / 2)
            # y_before = y_edge - y_after
            if self.Ndim==3:
                ypad = ((0,0), (0, y_edge), (0,0))
                newimage = np.pad(newimage, pad_width=ypad,
                                  mode='constant', constant_values=np.nan)
                newimage = newimage.reshape((self.Nw,Ny,yratio,Nx)).mean(2)
            else:
                ypad = ((0, y_edge), (0,0))
                newimage = np.pad(newimage, pad_width=ypad,
                                  mode='constant', constant_values=np.nan)
                newimage = newimage.reshape((Ny,yratio,Nx)).mean(1)
                
        ## Modify header
        ##---------------
        zval = oldw.all_pix2world(np.array([[0.5,0.5]]), 1)[0]
        # hdr['NAXIS1'] = Nx
        # hdr['NAXIS2'] = Ny
        hdr['CRPIX1'] = 0.5
        hdr['CRPIX2'] = 0.5
        hdr['CRVAL1'] = zval[0]
        hdr['CRVAL2'] = zval[1]
    
        cd = oldcd * [xratio,yratio]
        hdr['CD1_1'] = cd[0][0]
        hdr['CD2_1'] = cd[1][0]
        hdr['CD1_2'] = cd[0][1]
        hdr['CD2_2'] = cd[1][1]
    
        for kw in header.keys():
            if 'PC' in kw:
                del hdr[kw]
            if 'CDELT' in kw:
                del hdr[kw]

        if filOUT is not None:
            write_fits(filOUT, hdr, newimage, self.wvl, self.wmod)

        self.hdr = hdr
        self.im = newimage
    
        return newimage

class Jy_per_pix_to_MJy_per_sr(improve):
    '''
    Convert image unit from Jy/pix to MJy/sr

    ------ INPUT ------
    filIN               input FITS file
    filOUT              output FITS file
    ------ OUTPUT ------
    '''
    def __init__(self, filIN, filOUT=None, wmod=0, verbose=False):
        super().__init__(filIN, wmod, verbose)

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
        super().__init__(filIN, wmod, verbose)
        
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
                 pixscale=None,
                 filUNC=None, dist=None, wmod=0, verbose=False):
        super().__init__(filIN, wmod, verbose)
        
        if dist=='norm':
            self.rand_norm(filUNC)
        elif dist=='splitnorm':
            self.rand_splitnorm(filUNC)

        im_rebin = self.rebin(filOUT=filOUT, pixscale=pixscale)

    def header(self):
        return self.hdr
        
    def image(self):
        return self.im

    def wave(self):
        return self.wvl
        
class imontage(improve):
    '''
    2D image or 3D cube montage toolkit
    i means <improve>-based or initialize

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
            path_tmp = os.getcwd()+'/tmp_proc/'
        else:
            path_tmp = tmpdir
        if not os.path.exists(path_tmp):
            os.makedirs(path_tmp)
        self.path_tmp = path_tmp

        ## Inputs
        self.flist = allist(flist)
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
            self.rand_splitnorm([filIN[i]+'_unc_N', filIN[i]+'_unc_P'])
        
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
        flist = allist(flist)
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
        flist = allist(flist)
        dataset = type('', (), {})()
        comment = "Created by <imontage>"

        for j in trange(Nmc+1, leave=False,
                        desc='<imontage> Coadd (MC)'):
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
                        desc='<imontage> Coadd (MC)'):
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
                    print('Warning: The keywords center and pixscale are dumb. ')

                self.refheader = refheader
        else:
            ## Input files in list object
            flist = allist(flist)
                
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
                    print('Warning: The keywords center and pixscale are dumb. ')

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
                uncpdf=None, filOUT=None, tmpdir=None):
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
        uncpdf              add uncertainties (filename+'_unc.fits' needed)
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
        flist = allist(flist)
        
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
            if uncpdf=='norm':
                self.rand_norm(flist[i]+'_unc')
            elif uncpdf=='splitnorm':
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
                    edgeidx = np.ma.array(oldweight, 
                        mask=np.logical_and(oldweight==0, newweight!=0)).mask
                    if edgeidx.any():
                        newimage[edgeidx] = edgeimage[edgeidx]

                    oldweight = read_fits(path_tmp+'coadd.weight').data
                    if np.sum(oldweight==0)!=0:
                        SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+image_files[k],
                            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                        edgeimage = read_fits(path_tmp+'coadd').data
                        newweight = read_fits(path_tmp+'coadd.weight').data
                        edgeidx = np.ma.array(oldweight, 
                            mask=np.logical_and(oldweight==0, newweight!=0)).mask
                        if edgeidx.any():
                            newimage[edgeidx] = edgeimage[edgeidx]

                        oldweight = read_fits(path_tmp+'coadd.weight').data
                        if np.sum(oldweight==0)!=0:
                            SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+image_files[k],
                                shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                            edgeimage = read_fits(path_tmp+'coadd').data
                            newweight = read_fits(path_tmp+'coadd.weight').data
                            edgeidx = np.ma.array(oldweight, 
                                mask=np.logical_and(oldweight==0, newweight!=0)).mask
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
            ma_zero = np.ma.array(newimage, mask=(newimage==0)).mask
            newimage[ma_zero] = np.nan
            # write_fits(path_comb+'coadd_'+str(k), newheader, newimage)
            # tqdm.write(str(old_pixel_fov))
            # tqdm.write(str(new_pixel_fov))
            # tqdm.write(str(abs(newheader['CD1_1']*newheader['CD2_2'])))

            hyperimage.append(newimage)

        hyperimage = np.array(hyperimage)

        if cropedge:
            ## coadd.ref.fits here is 3D, different from that created by SWarp
            write_fits(path_tmp+'coadd.ref', newheader, hyperimage, wvl)
            
            reframe = improve(path_tmp+'coadd.ref')
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
        self.kfile = allist(kfile)

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
        [ref]
        Pereira-Santaella, Miguel, Almudena Alonso-Herrero, George H.
        Rieke, Luis Colina, Tanio Díaz-Santos, J.-D. T. Smith, Pablo G.
        Pérez-González, and Charles W. Engelbracht. “Local Luminous
        Infrared Galaxies. I. Spatially Resolved Observations with the
        Spitzer Infrared Spectrograph.” The Astrophysical Journal
        Supplement Series 188, no. 2 (June 1, 2010): 447.
        doi:10.1088/0067-0049/188/2/447.
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
    #     flist = allist(flist)
        
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
        flist = allist(flist)
        
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

class respect(improve):
    '''
    REstore SPECTra
    '''
    def __init__(self, verbose=False, tmpdir=None):
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
        # if tmpdir is None:
        #     path_tmp = os.getcwd()+'/tmp_rsp/'
        # else:
        #     path_tmp = tmpdir
        # if not os.path.exists(path_tmp):
        #     os.makedirs(path_tmp)
        
        # self.path_tmp = path_tmp

    def concat(self, flist, filOUT=None, comment=None,
               wsort=False, wrange=None):
        '''
        wsort=True can be used with wclean

        '''
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
                super().__init__(f)
                data.append(self.im)
                wave.append(self.wvl)
        ## Keep wavelengths in the given ranges (wrange)
        else:
            for f in flist:
                super().__init__(f)
                imin = closest(wmin, self.wvl[0])
                imax = closest(wmax, self.wvl[-1])
                iwi = 0
                iws = -1
                for i, w in enumerate(self.wvl):
                    if w<wrange[imin][0] and self.wvl[i+1]>wrange[imin][0]:
                        iwi = i+1
                    if w<wrange[imax][1] and self.wvl[i+1]>wrange[imax][0]:
                        iws = i
                data.append(self.im[iwi:iws])
                wave.append(self.wvl[iwi:iws])

        data = np.concatenate(data)
        wave = np.concatenate(wave)
        hdr = self.hdr
        ind = sorted(range(len(wave)), key=wave.__getitem__)
        wave = np.sort(wave)
        data = data[ind]
        self.wvl_concat = wave
        self.im_concat = data

        ## Write FITS file
        if filOUT is not None:
            write_fits(filOUT, hdr, data, wave, COMMENT=comment)

    def smooth(self, filIN, filUNC=None, BG_image=None, zerovalue=np.nan,
               wmin=None, wmax=None, lim_unc=1.e2, fltr_pn=None, cmin=5,
               filOUT=None):
        '''
        Remove spectral artifacts (Interpolate aberrant wavelengths)
        Anormaly if:
          abs(v - v_med) / unc > lim_unc

        ------ INPUT ------
        filIN               input spectral map (FITS)
        filUNC              input uncertainty map (FITS)
        filOUT              output smoothed spectral map (FITS)
        BG_image            background image used to generate unc map
        zerovalue           value used to replace zero value (Default:NaN)
        wmin                wavelength range to smooth (float)
        wmax                wavelength range to smooth (float)
        lim_unc             uncertainty dependant factor limit (positive float)
        fltr_pn             positive/negtive filter (Default: None)
                              'p' - smooth only positive aberrant
                              'n' - smooth only negtive aberrant
        cmin                minimum neighboring artifacts
        ------ OUTPUT ------
        im                  smoothed spectral map
        '''
        super().__init__(filIN)

        im = self.im
        wvl = self.wvl
        unc = self.uncert(filUNC=filUNC,BG_image=BG_image,zerovalue=zerovalue)

        if wmin is None:
            wmin = wvl[0]
        iwi = allist(wvl).index(wvl[closest(wvl,wmin)])
        if wmax is None:
            wmax = wvl[-1]
        iws = allist(wvl).index(wvl[closest(wvl,wmax)])

        if lim_unc<0:
            raise ValueError('lim_unc must be positive!')

        ## Scan every pixel/spectrum at each wavelength
        for w in trange(self.Nw, leave=False,
                        desc='<respect> smooth spectral map'):
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
            comment = "A <respect> smoothed spectral map"
            write_fits(filOUT, self.hdr, im, wvl,
                       COMMENT=comment)

        return im
            
    def mask(self):
        '''
        '''
        pass
        
class sextract(improve):
    '''
    AKARI/IRC spectroscopy slit coord extraction
    s means slit, spectral cube or SAV file

    ------ INPUT ------
    filOUT              output FITS file
    pathobs             path of IRC dataset
    parobs[0]           observation id
    parobs[1]           slit name
    parobs[2]           IRC N3 (long exp) frame (2MASS corrected; 90 deg rot)
    parobs[3]           NG(grism)/NP(prism)
    Nw                  num of wave
    Ny                  slit length
    Nx                  slit width
    ------ OUTPUT ------
    '''
    def __init__(self, pathobs=None, parobs=None, verbose=False):
        self.path = pathobs + parobs[0] + '/irc_specred_out_' + parobs[1]+'/'
        filIN = self.path + parobs[2]
        super().__init__(filIN)
        
        self.filSAV = self.path + parobs[0] + '.N3_' + parobs[3] + '.IRC_SPECRED_OUT'
        self.table = readsav(self.filSAV+savext, python_dict=True)['source_table']

        ## Slit width will be corrected during reprojection
        if parobs[1]=='Ns':
            self.slit_width = 3 # 5"/1.446" = 3.458 pix (Ns)
        elif parobs[1]=='Nh':
            self.slit_width = 2 # 3"/1.446" = 2.075 pix (Nh)

        if verbose==True:
            print('\n----------')
            print('Slit extracted from ')
            print('obs_id: {} \nslit: {}'.format(parobs[0], parobs[1]))
            print('----------\n')

    def rand_pointing(self, sigma=0.):
        '''
        Add pointing uncertainty to WCS

        ------ INPUT ------
        sigma               pointing accuracy (deg)
        ------ OUTPUT ------
        '''
        d_ro = abs(np.random.normal(0., sigma)) # N(0,sigma)
        d_phi = np.random.random() *2. * np.pi # U(0,2*pi)
        self.hdr['CRVAL1'] += d_ro * np.cos(d_phi)
        self.hdr['CRVAL2'] += d_ro * np.sin(d_phi)

        return d_ro, d_phi

    def spec_build(self, filOUT=None, write_unc=True,
                   Nx=0, Ny=32, Nsub=1, sig_pt=0.):
        '''
        Build the spectral cube/slit from spectra extracted by IDL pipeline
        (see IRC_SPEC_TOOL, plot_spec_with_image)

        ------ INPUT ------
        Nx                  number of (identical) pixels to fit slit width
                              Default: 0, 3 for Ns and 2 for Nh
        Ny                  number of pixels in spatial direction (Max=32)
                              Y axis in N3 frame (X axis in focal plane arrays)
        Nsub                number of subslits
        '''
        if Nx==0:
            Nx = self.slit_width
        ref_x = self.table['image_y'][0] # slit ref x
        ref_y = 512 - self.table['image_x'][0] # slit ref y

        ## Get slit coord from 2MASS corrected N3 frame
        ## Do NOT touch self.im (N3 frame, 2D) before this step
        self.crop(sizpix=(Nx, Ny), cenpix=(ref_x, ref_y))
        # self.hdr['CTYPE3'] = 'WAVE-TAB'
        self.hdr['CUNIT1'] = 'deg'
        self.hdr['CUNIT2'] = 'deg'
        self.hdr['BUNIT'] = 'MJy/sr'
        self.hdr['EQUINOX'] = 2000.0

        ## Add pointing unc
        self.rand_pointing(sig_pt)

        ## Read spec
        spec_arr = []
        for j in range(Ny):
            ## Ny/Nsub should be integer, or there will be shift
            ispec = math.floor(j / (math.ceil(Ny/Nsub)))
            # spec_arr.append(read_ascii(self.path+'spec'+str(ispec), '.spc', float))
            allslit = ascii.read(self.path+'spec'+str(ispec)+'.spc')
            subslit = []
            for k in allslit.keys():
                subslit.append(allslit[k])
            subslit = np.array(subslit)
            
            spec_arr.append(subslit)
        ## spec_arr.shape = (Ny,4,Nw)
        spec_arr = np.array(spec_arr)
        Nw = len(spec_arr[0,0,:])
        
        ## Broaden cube width
        cube = np.empty([Nw,Ny,Nx])
        unc = np.empty([Nw,Ny,Nx]) # Symmetric unc
        unc_N = np.empty([Nw,Ny,Nx]) # Asymmetric negtive
        unc_P = np.empty([Nw,Ny,Nx]) # Asymmetric positive
        wave = np.empty(Nw)
        for k in range(Nw):
            for j in range(Ny):
                for i in range(Nx):
                    cube[k][j][i] = spec_arr[j,1,k]
                    unc[k][j][i] = (spec_arr[j,3,k]-spec_arr[j,2,k])/2
                    unc_N[k][j][i] = (spec_arr[j,1,k]-spec_arr[j,2,k])
                    unc_P[k][j][i] = (spec_arr[j,3,k]-spec_arr[j,1,k])
            wave[k] = spec_arr[0,0,k]

        ## Save spec in wave ascending order$
        self.cube = cube[::-1]
        self.unc = unc[::-1]
        self.unc_N = unc_N[::-1]
        self.unc_P = unc_P[::-1]
        self.wvl = wave[::-1]

        if filOUT is not None:
            comment = "Assembled AKARI/IRC slit spectroscopy cube. "
            write_fits(filOUT, self.hdr, self.cube, self.wvl,
                       COMMENT=comment)

            if write_unc==True:
                uncom = "Assembled AKARI/IRC slit spec uncertainty cube. "
                write_fits(filOUT+'_unc', self.hdr, self.unc, self.wvl,
                           COMMENT=uncom)

                uncom_N = "Assembled AKARI/IRC slit spec uncertainty (N) cube. "
                write_fits(filOUT+'_unc_N', self.hdr, self.unc_N, self.wvl,
                           COMMENT=uncom)

                uncom_P = "Assembled AKARI/IRC slit spec uncertainty (P) cube. "
                write_fits(filOUT+'_unc_P', self.hdr, self.unc_P, self.wvl,
                           COMMENT=uncom)

        return self.cube

    def sav_build(self):
        '''
        Alternative extraction from SAV file
        Including wave calib, ?, etc. 
        (see IRC_SPEC_TOOL, plot_spec_with_image)
        '''
        filSAV = self.filSAV
        table = self.table
        ## Read SAV file
        image = readsav(filSAV+savext, python_dict=True)['specimage_n_wc']
        image = image[::-1] # -> ascending order
        noise = readsav(filSAV+savext, python_dict=True)['noisemap_n']
        noise = noise[::-1]
        wave = readsav(filSAV+savext, python_dict=True)['wave_array']
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
        return self.cube

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
                    print('WARNING: Left side fully covered! ')
                
                for j in range(Nw-i-1):
                    dw = wave[i+1+j] - wave[i]
                    if dw>0: # found the closest larger wave[i+1+j]
                        wmax = i+1+j
                        break
                if wmax==0:
                    print('WARNING: right side fully covered! ')

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
            edgeidx = np.ma.array(oldweight, 
                mask=np.logical_and(oldweight==0, newweight!=0)).mask
            if edgeidx.any():
                newimage[edgeidx] = edgeimage[edgeidx]

            oldweight = read_fits(path_tmp+'coadd.weight').data
            if np.sum(oldweight==0)!=0:
                SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE BILINEAR '+' old.fits', \
                    shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                edgeimage = read_fits(path_tmp+'coadd').data
                newweight = read_fits(path_tmp+'coadd.weight').data
                edgeidx = np.ma.array(oldweight, 
                    mask=np.logical_and(oldweight==0, newweight!=0)).mask
                if edgeidx.any():
                    newimage[edgeidx] = edgeimage[edgeidx]

                oldweight = read_fits(path_tmp+'coadd.weight').data
                if np.sum(oldweight==0)!=0:
                    SP.call('swarp '+swarp_opt+' -RESAMPLING_TYPE NEAREST '+' old.fits',
                            shell=True, cwd=path_tmp, stdout=devnull, stderr=SP.STDOUT)
                    edgeimage = read_fits(path_tmp+'coadd').data
                    newweight = read_fits(path_tmp+'coadd.weight').data
                    edgeidx = np.ma.array(oldweight, 
                        mask=np.logical_and(oldweight==0, newweight!=0)).mask
                    if edgeidx.any():
                        newimage[edgeidx] = edgeimage[edgeidx]

    ## SWarp is conserving surface brightness/pixel
    ## while the pixels size changes
    oldcdelt = get_pc(wcs=fixwcs(header=oldheader).wcs).cdelt
    refcdelt = get_pc(wcs=fixwcs(header=refheader).wcs).cdelt
    old_pixel_fov = abs(oldcdelt[0]*oldcdelt[1])
    new_pixel_fov = abs(refcdelt[0]*refcdelt[1])
    newimage = newimage * old_pixel_fov/new_pixel_fov
    ma_zero = np.ma.array(newimage, mask=(newimage==0)).mask
    newimage[ma_zero] = np.nan
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
                wsort=False, wrange=None):
    '''
    wsort=True can be used with wclean

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
            for i, w in enumerate(ds.wave):
                if w<wrange[imin][0] and ds.wave[i+1]>wrange[imin][0]:
                    iwi = i+1
                if w<wrange[imax][1] and ds.wave[i+1]>wrange[imax][0]:
                    iws = i
            data.append(ds.data[iwi:iws])
            wave.append(ds.wave[iwi:iws])

    data = np.concatenate(data)
    wave = np.concatenate(wave)
    hdr = ds.header    
    ind = sorted(range(len(wave)), key=wave.__getitem__)
    wave = np.sort(wave)
    data = data[ind]
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

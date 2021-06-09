#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging
# logging.disable(sys.maxsize)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/lib/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

## Local
sys.path.insert(0, testdir+'/..') ## laputan path
from laputan.inout import fitsext, read_fits, write_fits
from laputan.astrom import fixwcs
from laputan.imaging import (improve, islice, icrop, iswarp, iconvolve,
                              sextract, wclean, interfill, hextract, hswarp)

print('\n TEST improve ')
print('--------------')
imp = improve(datdir+'M83', verbose=True)
print('raw: ', imp.im[0,0,:])
imp.rand_norm(datdir+'M83_unc')
print('add sym rand: ', imp.im[0,0,:])
imp.rand_splitnorm([datdir+'M83_unc', datdir+'M83_unc'])
print('add split rand', imp.im[0,0,:])
imp.crop(outdir+'M83_imp_crop', sizval=(0.005,0.009), cenpix=(7,8))
print('See out/M83_imp_crop.fits [Done]')

print('\n TEST islice ')
print('-------------')
slc = islice(datdir+'M83', filSL=outdir+'M83_inv_sqrt',
             filUNC=datdir+'M83_unc', slicetype='inv_sq', postfix='_test_')
print(slc.filenames()[0])
if input('Clean slices (y/n): ')=='y':
    slc.clean()

print('\n TEST icrop ')
print('------------')
crp = icrop(datdir+'M83', filOUT=outdir+'M83_icrop',
            sizpix=(3,6), cenval=(204.2529675, -29.8656962),
            filUNC=[datdir+'M83_unc', datdir+'M83_unc'],
            dist='splitnorm', verbose=True)
print('See out/M83_icrop.fits [Done]')

print('\n TEST iswarp ')
print('-------------')
hdr_ref = fixwcs(datdir+'M82_09_L86'+fitsext).header
swp = iswarp(refheader=hdr_ref,
             center=None, pixscale=6.,
             tmpdir=outdir+'swp/')
swp.combine(datdir+'M82_08_L86', uncpdf='norm',
            filOUT=outdir+'M82_08_L86_rep', tmpdir=outdir+'swp/')
print('Reproject M82_08_L86 to M82_09_L86 (pixscale=6") [Done]')

swp_cube = iswarp((datdir+'M82_09_L86', datdir+'M82_04_SL1'),
                  center='9:55:51,69:40:45', pixscale=None,
                  tmpdir=outdir+'swp_cube/')
swp_cube.combine(datdir+'M82_04_SL1', uncpdf='norm',
                 filOUT=outdir+'M82_04_SL1_rep')
print('Reproject M82_04_SL1 to M82_09_L86 (recenter to M82 center) [Done]')

print('\n TEST iconvolve ')
print('----------------')
## See also idl/conv_prog.pro & idl/convolve_image.pro
convdir = outdir+'conv/' # see also idl/convolve_image.pro
if not os.path.exists(convdir):
    os.makedirs(convdir)
path_ker = datdir
path_idl = testdir+'/../idl/'
csv_ker = outdir+'kernelist' # see also idl/conv_prog.pro

irs_ker = []
psf = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
psf_ref = 'Gauss_06.0'
for p in psf:
    irs_ker.append(path_ker+'Kernel_HiRes_Gauss_0'+str(p)+'_to_'+psf_ref)
conv_cube = iconvolve(datdir+'M82_04_SL1',
                      kfile=irs_ker, klist=csv_ker,
                      convdir=convdir, filOUT=outdir+'M82_04_SL1'+'_conv')
conv_cube.do_conv(path_idl, verbose=False)
print('Convolve M82_04_SL1 [Done]')

irac_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
conv = iconvolve(datdir+'M82_IRAC4',
                 kfile=irac_ker, klist=csv_ker,
                 filOUT=outdir+'M82_IRAC4'+'_conv')
conv.do_conv(path_idl)
print('Convolve M82_IRAC4 [Done]')

print('\n TEST sextract ')
print('---------------')

print('\n TEST wclean ')
print('-------------')

print('\n TEST interfill ')
print('----------------')
data = read_fits(datdir+'IC10_SL2').data[0]
hdr = fixwcs(datdir+'IC10_SL2'+fitsext).header
newdata = interfill(data, axis=0)
write_fits(outdir+'IC10_fillgap', hdr, newdata)
print('interfill test [Done]')

print('\n TEST hextract ')
print('---------------')
hextract(datdir+'M82_08_L86', outdir+'M82_hextract', 20, 40, 32, 42)
print('hextract test [Done]')

print('\n TEST hswarp ')
print('-------------')
print('Reproject M82_08_L86 to M82_09_L86')
old = read_fits(datdir+'M82_08_L86')
oldimage = old.data
oldheader = old.header
refheader = read_fits(datdir+'M82_09_L86').header
hswp = hswarp(oldimage, oldheader, refheader, keepedge=True,
              tmpdir=outdir+'hswp/', verbose=True)
# print('hswarp image: ', hswp.image)
# print('hswarp image header: ', hswp.header)
print('See hswp/coadd.fits [Done]')

if input('Clean tmp files (y/n): ')=='y':
    swp.clean()
    swp_cube.clean()
    conv_cube.clean()
    conv.clean()
    conv.clean(outdir+'hswp/')


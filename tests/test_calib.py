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
import matplotlib.pyplot as plt

## Local
sys.path.insert(0, testdir+'/..') ## laputan path
from laputan.calib import (intercalib, photometry_profile, 
)

print('\n TEST intercalib ')
print('-----------------')
phots = 'IRAC1', 'IRAC4'
w_spec = np.arange(1,20,.1)
ic = intercalib()
ic.read_filter(phots, w_spec)
print('Wavelength center of the filters (IRAC1, IRAC4): ', ic.wcen)
print('spec off / broad band off: ', ic.specoff_ov_bboff)
phots = 'MIPS1'
ic = intercalib()
ic.read_filter(phots, w_spec)
print('Wavelength center of the filters (MIPS1): ', ic.wcen)
print('spec off / broad band off: ', ic.specoff_ov_bboff)

print('* via FITS file *')
c1 = intercalib(datdir+'M82')
sp1 = c1.synthetic_photometry(phots)
print('Fnu_filt = ', sp1.Fnu_filt.shape)
print('wcen = ', sp1.wcen)
print('smat = ', sp1.smat)

print('* via data array *')
x,y = 2,0
with fits.open(datdir+'M82.fits') as hdul:
    spec = hdul[0].data[:,y,x]
    wave = hdul[1].data
c2 = intercalib()
sp2 = c2.synthetic_photometry('IRAC3', wave, spec, verbose=False)
print('Fnu_filt = ', sp2.Fnu_filt)
print('wcen = ', sp2.wcen)
print('smat = ', sp2.smat)

print('\n TEST correct_spec ')
print('-------------------')
a1 = 1.5 # slope
b1 = -200. # offset
a2 = .5
b2 = 100.
wlim11=(None,5.2)
wlim12=(5.2,14.5)
wlim2=(5.2,7.6)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(8,9))
plt.subplots_adjust(left=.1, bottom=.05, \
                    right=.99, top=.95, wspace=0., hspace=.2)
ax1, ax2 = axes

print('* via FITS file *')

new_spec11 = c1.correct_spec(gain=a1, offset=b1, wlim=wlim11)
new_spec12 = c1.correct_spec(gain=a2, offset=b2, wlim=wlim12)
ax1.plot(c1.wvl, c1.im[:,y,x], c='k', label='y')
ax1.plot(c1.wvl, new_spec11[:,y,x], c='y',
         label='{}*y+{}, x=({},{})'.format(a1,b1,wlim11[0],wlim11[1]))
ax1.plot(c1.wvl, new_spec12[:,y,x], c='g',
         label='{}*y+{}, x=({},{})'.format(a2,b2,wlim12[0],wlim12[1]))
ax1.legend(loc='upper left')

print('* via data array *')
new_spec2 = c2.correct_spec(gain=a1, offset=b1, w_spec=wave, Fnu_spec=spec, wlim=wlim2)
ax2.plot(wave, spec, c='k', label='y')
ax2.plot(wave, new_spec2, c='y',
         label='{}*y+{}, x=({},{})'.format(a1,b1,wlim2[0],wlim2[1]))
ax2.legend(loc='upper left')

fig.savefig(outdir+'calib_correct_spec.png')
print('See ./out/calib_correct_spec.png [Done]')

print('\n TEST photometry_profile ')
print('-------------------------')
phot_lib = photometry_profile(None,
                              'IRAC1', 'IRAC2', 'IRAC3', 'IRAC4', 'MIPS1',)
                              # 'WISE1', 'WISE2', 'WISE3', 'WISE4',)
phot_lib.save(outdir+'calib_phot', transparent=1)
print('See ./out/calib_phot.png [Done]')

# plt.show()

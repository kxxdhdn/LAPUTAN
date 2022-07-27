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
# import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter

## rapyuta
from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot


## Input/output
filOUT = outdir+'psf_irs_hires.png'

sim_par_wave = [9.89, 19.51, 40.]
sim_par_fwhm = [3.94, 4.68, 10.1]
sim_per_wave = [9.89, 19.51, 40.]
sim_per_fwhm = [5.08, 5.24, 10.1]

# sim__wave = [9.89, 19.51] # SH
# sim__wave = [18.83, 37.14] # LH


wvl = np.linspace(9.89, 37.14, 1000)

## FWHM (arcsec)
fwhm_par = np.interp(wvl, sim_par_wave, sim_par_fwhm)
fwhm_per = np.interp(wvl, sim_per_wave, sim_per_fwhm)
fwhm_lam = np.sqrt(fwhm_par * fwhm_per)

## sigma (arcsec)
sigma_par = fwhm_par / (2. * np.sqrt(2.*np.log(2.)))
sigma_per = fwhm_per / (2. * np.sqrt(2.*np.log(2.)))
sigma_lam = np.sqrt(sigma_par * sigma_per)

colors = ['k','c','m']

i = closest(wvl, 20)
print(fwhm_par[i])
exit()

p = pplot(wvl, fwhm_lam, label='Geometric mean',
          xlim=(4.5, 44), ylim=(2.6, 10.5),
          xlog=1, ylog=0, clib=colors, lw=3, ls='--',
          xlabel=r'$\rm Wavelength,\ \lambda\ [\mu m]$',
          ylabel=r'$\rm FWHM\ [^{\prime\prime}]$',
          title=None, zorder=100, alpha=.8,
          figsize=(8,6), left=.1, bottom=.15, top=.95, right=.95,
          titlesize=20, xysize=20, tksize=20)
p.add_plot(wvl, fwhm_par, label='Parallel', lw=5, alpha=1)
p.add_plot(wvl, fwhm_per, label='Perpendicular', lw=5, alpha=1)
p.ax.legend(loc='upper left',fontsize=20, framealpha=0)
xtic = [5, 6, 7, 8, 9, 11, 13, 15, 20, 30, 40]
p.ax.set_xticks(xtic, minor=False) # major
p.ax.xaxis.set_major_formatter(ScalarFormatter()) # major

p.save(filOUT, transparent=False, figtight=True)

# p.show()

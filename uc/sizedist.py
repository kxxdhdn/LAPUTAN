#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Size distribution of interstellar grains via various dust models

"""
import sys, os
import numpy as np
from scipy.constants import c, physical_constants
from astropy.io import ascii
import matplotlib.transforms as mtransforms
from matplotlib.ticker import FuncFormatter

from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import plotool

## Print full numpy array
np.set_printoptions(threshold=sys.maxsize)

## Physical constants
h = physical_constants['Planck constant in eV s'][0]
# print(h)
hbar = physical_constants['reduced Planck constant in eV s'][0]
# print(hbar)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/lib/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

## INPUT
##-------
## dataset:
## /rmic
## /f_PAH_Z04, /f_gra_Z04, /f_sil_Z04
## /f_PAH_C11, /f_LamC_C11, /f_SamC_C11, /f_sil_C11
## /f_aSil_J13, /f_bHAC_J13, /f_esHAC_J13, /f_sHAC_J13, /f_vsHAC_J13
## /f_aSil_J17, /f_bHAC_J17, /f_esHAC_J17, /f_sHAC_J17, /f_vsHAC_J17
filin = datdir+'sizedist_G21hdr'
rmic = read_hdf5(filin, 'rmic')
f_PAH_Z04 = read_hdf5(filin, 'f_PAH_Z04')
f_gra_Z04 = read_hdf5(filin, 'f_gra_Z04')
f_sil_Z04 = read_hdf5(filin, 'f_sil_Z04')
f_PAH_C11 = read_hdf5(filin, 'f_PAH_C11')
f_LamC_C11 = read_hdf5(filin, 'f_LamC_C11')
f_SamC_C11 = read_hdf5(filin, 'f_SamC_C11')
f_sil_C11 = read_hdf5(filin, 'f_sil_C11')
f_aSil_J13 = read_hdf5(filin, 'f_aSil_J13')
f_bHAC_J13 = read_hdf5(filin, 'f_bHAC_J13')
f_esHAC_J13 = read_hdf5(filin, 'f_esHAC_J13')
f_sHAC_J13 = read_hdf5(filin, 'f_sHAC_J13')
f_vsHAC_J13 = read_hdf5(filin, 'f_vsHAC_J13')
f_aSil_J17 = read_hdf5(filin, 'f_aSil_J17')
f_bHAC_J17 = read_hdf5(filin, 'f_bHAC_J17')
f_esHAC_J17 = read_hdf5(filin, 'f_esHAC_J17')
f_sHAC_J17 = read_hdf5(filin, 'f_sHAC_J17')
f_vsHAC_J17 = read_hdf5(filin, 'f_vsHAC_J17')

varlist = [f_PAH_Z04, f_gra_Z04, f_sil_Z04,
           f_PAH_C11, f_LamC_C11, f_SamC_C11, f_sil_C11,
           f_aSil_J13, f_bHAC_J13, f_esHAC_J13, f_sHAC_J13, f_vsHAC_J13,
           f_aSil_J17, f_bHAC_J17, f_esHAC_J17, f_sHAC_J17, f_vsHAC_J17,
           ]
for var in varlist:
    var[var==0] = None

## MRN
a_gra_MRN = []
a_sil_MRN = []
for r in rmic:
    if r>0.005 and r<1:
        a_gra_MRN.append(r)
    if r>0.025 and r<0.25:
        a_sil_MRN.append(r)
a_gra_MRN = np.array(a_gra_MRN)
a_sil_MRN = np.array(a_sil_MRN)
f_gra_MRN = 10**-15.13 * a_gra_MRN**(-3.5)
f_sil_MRN = 10**-15.10 * a_sil_MRN**(-3.5)

## PLOT
##------
pl = plotool(1, 3, figsize=(12,4))
pl.trans = 'Axes'
lw1 = 2
pl.set_fig(wspace=0)
## Zubko2004
pl.set_ax(xlog=True, ylog=True, xlim=(3e-4,9e-1), ylim=(2e-9,9e-6),
          xtkform='log_sci', ytkform='log_sci',# xtktoggle=True,)
          # xlabel=r'Radius, a $\rm [\mu m]$',
          ylabel=r'Size distribution, $\rm a^4f(a) [10^{-21}\ cm^3/H]$',)
pl.plot(rmic, f_PAH_Z04*rmic**4*1e9, c='m', label='PAH')
pl.plot(rmic, f_gra_Z04*rmic**4*1e9, c='orange', label='Graphite')
pl.plot(rmic, f_sil_Z04*rmic**4*1e9, c='c', label='Silicates')
pl.plot(a_gra_MRN, f_gra_MRN*a_gra_MRN**4*1e9, c='grey', ls='--')
pl.plot(a_sil_MRN, f_sil_MRN*a_sil_MRN**4*1e9, c='grey', ls=':')
pl.ax.text(0.4, 0.5, 'MRN', c='grey', size=15, transform=pl.ax.transAxes)
pl.ax.text(0.85, 0.9, '(a)', style='italic', c='grey', size=15,
           transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='upper center', title='Zubko et al. (2004)')
pl.reset_handles()
pl.labels = []
## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

## Compiègne2011
pl.set_ax((0,1), xlog=True, ylog=True,xlim=(3e-4,9e-1), ylim=(2e-9,9e-6),
          xtkform='log_sci', ytkform='log_sci',
          # xtktoggle=True, ytktoggle=True,
          xlabel=r'Radius, a $\rm [\mu m]$',)
          # ylabel=r'Size distribution, ${\rm a}^4f({\rm a}) \rm [10^{-21}\ cm^3/H]$',)
pl.plot(rmic, f_PAH_C11*rmic**4*1e9, c='m', label='PAH')
pl.plot(rmic, f_LamC_C11*rmic**4*1e9, c='orange', label='a-C')
pl.plot(rmic, f_SamC_C11*rmic**4*1e9, c='orange', label='a-C')
pl.plot(rmic, f_sil_C11*rmic**4*1e9, c='c', label='Silicates')
pl.plot(a_gra_MRN, f_gra_MRN*a_gra_MRN**4*1e9, c='grey', ls='--')
pl.plot(a_sil_MRN, f_sil_MRN*a_sil_MRN**4*1e9, c='grey', ls=':')
pl.ax.text(0.4, 0.5, 'MRN', c='grey', size=15, transform=pl.ax.transAxes)
pl.ax.text(0.85, 0.9, '(b)', style='italic', c='grey', size=15,
           transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='upper center', title='Compiègne et al. (2011)')
pl.reset_handles()
pl.labels = []
## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

## Jones2017
pl.set_ax((0,2), xlog=True, ylog=True, xlim=(3e-4,9e-1), ylim=(2e-9,9e-6),
          xtkform='log_sci', ytkform='log_sci',)
          # xlabel=r'Radius, a $\rm [\mu m]$',
          # ylabel=r'Size distribution, ${\rm a}^4f({\rm a}) \rm [10^{-21}\ cm^3/H]$',)
pl.plot(rmic, f_esHAC_J17*rmic**4*1e9, c='m', label='a-C(:H)')
pl.plot(rmic, f_vsHAC_J17*rmic**4*1e9, c='m', label='a-C(:H)')
pl.plot(rmic, f_sHAC_J17*rmic**4*1e9, c='m', label='a-C(:H)')
pl.plot(rmic, f_bHAC_J17*rmic**4*1e9, c='orange', label='a-C/a-C(:H)')
pl.plot(rmic, f_aSil_J17*rmic**4*1e9, c='c', label='a-C/a-Sil')
pl.plot(a_gra_MRN, f_gra_MRN*a_gra_MRN**4*1e9, c='grey', ls='--')
pl.plot(a_sil_MRN, f_sil_MRN*a_sil_MRN**4*1e9, c='grey', ls=':')
pl.ax.text(0.4, 0.5, 'MRN', c='grey', size=15, transform=pl.ax.transAxes)
pl.ax.text(0.85, 0.9, '(c)', style='italic', c='grey', size=15,
           transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='upper center', title='Jones et al. (2017)')
pl.reset_handles()
## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

## LD01
# pl.set_ax((1,1), xlog=True, ylog=True, xlim=(1e-4,9e-1), ylim=(2e-9,9e-6),
#           xtkform='log_sci', ytkform='log_sci',
#           ytktoggle=True)
#           # xlabel=r'Radius, a $\rm [\mu m]$',
#           # ylabel=r'Size distribution, $\rm a^4f(a) [10^{-21}\ cm^3/H]$',)
# pl.plot(rmic, f_esHAC_J13*rmic**4*1e9, c='m', label='a-C(:H)')
# pl.plot(rmic, f_bHAC_J13*rmic**4*1e9, c='orange', label='a-C/a-C(:H)')
# pl.plot(rmic, f_aSil_J13*rmic**4*1e9, c='c', label='a-C/a-Sil')
# pl.plot(a_gra_MRN, f_gra_MRN*a_gra_MRN**4*1e9, c='grey', ls='--')
# pl.plot(a_sil_MRN, f_sil_MRN*a_sil_MRN**4*1e9, c='grey', ls=':')
# pl.ax.text(0.4, 0.5, 'MRN', c='grey', size=15, transform=pl.ax.transAxes)
# pl.ax.text(0.92, 0.9, '(c)', style='italic', c='grey', size=15,
#            transform=pl.ax.transAxes)
# pl.append_handles()
# pl.set_legend(loc='upper center', title='WD01')
# pl.reset_handles()
# pl.labels = []
# ## aesthetics
# for spine in pl.ax.spines.items():
#     pl.ax.spines[spine[0]].set_linewidth(2) # border width
# pl.ax.xaxis.set_ticks_position('both') # ticks
# pl.ax.yaxis.set_ticks_position('both')
# pl.ax.tick_params(axis='both', which='both', direction='in')
# pl.ax.tick_params(which='major', length=8, width=2)
# pl.ax.tick_params(which='minor', length=4, width=2)

pl.save(outdir+'sizedist.png', transparent=False, figtight=True)

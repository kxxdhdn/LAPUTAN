#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Cross-sections (Qabs & Qsca) of 

"""
import sys, os
import numpy as np
from scipy.constants import c, physical_constants
from astropy.io import ascii
import matplotlib.transforms as mtransforms
from matplotlib.ticker import FuncFormatter

from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.maths import gaussian
from rapyuta.plots import plotool

## Print full numpy array
np.set_printoptions(threshold=sys.maxsize)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/lib/DH21Ad/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

## INPUT
##-------
wave = ascii.read(datdir+'DH21_wave'+'.txt')[0]
wave = np.lib.recfunctions.structured_to_unstructured(np.array(wave))
Nw = len(wave)
aeff = ascii.read(datdir+'DH21_aeff'+'.txt')[0]
aeff = np.lib.recfunctions.structured_to_unstructured(np.array(aeff))
Nr = len(aeff)
aeff = np.array(aeff)

radius = 1e-1
jr = closest(aeff, radius)
print('radius = ', aeff[jr])
data = ascii.read(datdir+'q_DH21Ad_P0.20_Fe0.00_0.500'+'.dat',
                  format='no_header')['col'+str(jr)]
## dataset:
## Ncol = Nr = 169
## Nrow = Nw * 9 = 1129 * 9
## - Qext=Cext/(pi aeff^2)
##   Qext1-3, Qabs1-3, Qsca1-3
##   - 3 orientations
##     jori=1: k  ||  a, E perp a
##     jori=2: k perp a, E  ||  a
##     jori=3: k perp a, E perp a
##     ran: randomly-oriented
q_ext1 = np.array(data[0*Nw:1*Nw])
q_ext2 = np.array(data[1*Nw:2*Nw])
q_ext3 = np.array(data[2*Nw:3*Nw])
qext_ran = (q_ext1 + q_ext2 + q_ext3)/3.

q_abs1 = np.array(data[3*Nw:4*Nw])
q_abs2 = np.array(data[4*Nw:5*Nw])
q_abs3 = np.array(data[5*Nw:6*Nw])
qabs_ran = (q_abs1 + q_abs2 + q_abs3)/3.

q_sca1 = np.array(data[6*Nw:7*Nw])
q_sca2 = np.array(data[7*Nw:8*Nw])
q_sca3 = np.array(data[8*Nw:9*Nw])
qsca_ran = (q_sca1 + q_sca2 + q_sca3)/3.
        
## PLOT
##------
pl = plotool(1, 1, figsize=(12,6))
# pl.trans = 'Axes'
lw1 = 3
pl.set_fig(hspace=.5)
pl.set_ax(xlim=(2e-2,5e2), ylim=(2e-9,9),
          xlog=1, ylog=1, xtkform='log_sci', ytkform='log_sci',
          xlabel=r'Wavelength, $\rm \lambda\ [\mu m]$',
          ylabel=r'$\rm Q_{abs}\ &\ Q_{sca}$',
          xsize=20, ysize=20, xtksize=20, ytksize=20,)
pl.plot(wave, qext_ran, c='k', ls='--', label=r'$\rm Q_{ext}= Q_{sca}+ Q_{abs}$', zorder=100)
pl.plot(wave, qsca_ran, c='c', lw=lw1, label=r'$\rm Q_{sca}$')
pl.plot(wave, qabs_ran, c='m', lw=lw1, label=r'$\rm Q_{abs}$')
pl.plot([1e-3,wave[0]], [1,qext_ran[0]], c='k', ls='--', zorder=100)
pl.plot([1e-3,wave[0]], [1,qsca_ran[0]], c='c', lw=lw1)
pl.plot([1e-3,wave[0]], [1,qabs_ran[0]], c='m', lw=lw1)
# pl.ax.text(0.05, 0.9, '(a)', style='italic', c='grey', size=20,
#            transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='lower left', title=r'a = '+str(aeff[jr])+r'$\rm\ \mu m$',
              fontsize=20, title_fontsize=20)
## Regime limits
w1, w2 = 0, 0
flag1, flag2 = True, True
for iw in range(Nw):
    if abs(qext_ran[iw]-2)>1e-1 and flag1: # Qext = 2 (delta < 1e-1)
        w1 = wave[iw]
        iw1 = iw
        flag1 = False
    if wave[iw]>2*np.pi*radius*50 and flag2: # wave >> 2*pi*a^2 (> *50)
        w2 = wave[iw]
        iw2 = iw
        flag2 = False
    if not (flag1 or flag2):
        break
print('Lower lim of Mie regime (um):', w1)
print('Upper lim of Mie regime (um):', w2)
## extrap & fit
pl.plot([1e-3,1], [1,1], c='grey', lw=lw1, ls='--', zorder=20)
pl.ax.text(1.5, 1, r'$\rm Q_{abs} \simeq Q_{sca} \simeq 1$', c='grey', size=20)
wave_fit = np.logspace(np.log10(w2),np.log10(wave[-1]),100)
pl.plot(wave_fit, wave_fit**-2 * qabs_ran[iw2]/w2**-2,
        c='grey', lw=lw1, ls='--', zorder=20)
pl.ax.text(5e1, 1.2*qabs_ran[closest(wave, 5e1)], r'$\rm Q_{abs} \propto \lambda^{-2}$', c='grey', size=20)
pl.plot(wave_fit, wave_fit**-4 * qsca_ran[iw2]/w2**-4,
        c='grey', lw=lw1, ls='--', zorder=20)
pl.ax.text(5e1, 1.2*qsca_ran[closest(wave, 5e1)], r'$\rm Q_{sca} \propto \lambda^{-4}$', c='grey', size=20)
## shadow
for iw in range(Nw):
    ## Geometrical optics
    if wave[iw]<=w1:
        Nw_extrap = 100
        for i in range(Nw_extrap):
            dw = (w1-1e-2)/Nw_extrap
            wave_extrap = 1e-2 + i*dw
            alpha_exp = 0.5 * np.exp(-i/Nw_extrap)
            # print(wave_extrap, alpha_exp)
            pl.ax.axvspan(wave_extrap, wave_extrap+dw, color=(0,0,1), ec='None', alpha=alpha_exp)
    ## Mie regime
    elif wave[iw]<=w2:
        alpha_exp = 0.5 * np.exp((-iw+iw1)/500)
        mu = (np.log(w1)+np.log(w2))/2
        sig = (np.log(w2)-np.log(w1))*0.3
        alpha_gauss = 0.5 * gaussian(np.log(wave[iw]),mu,sig)/gaussian(mu,mu,sig)
        # print(wave[iw], alpha_exp, alpha_gauss)
        pl.ax.axvspan(wave[iw-1], wave[iw], color=(0,1,0), ec='None', alpha=alpha_gauss)
    ## Rayleigh regime
    else:
        alpha_exp = 0.5 * np.exp((-iw+iw2)/500)
        alpha_power = 1.5 * (np.log(iw)-np.log(iw2))**0.5/np.log(wave[-1])**0.5 + 0.1
        # print(wave[iw], alpha_exp, alpha_power)
        pl.ax.axvspan(wave[iw-1], wave[iw], color=(1,0,0), ec='None', alpha=alpha_power)
pl.ax.text(8e-3, 50, 'Geometrical optics', c=(0,0,1), fontsize=20)
pl.ax.text(8e-1, 50, 'Mie regime', c=(0,1,0), fontsize=20)
pl.ax.text(40, 50, 'Rayleigh regime', c=(1,0,0), fontsize=20)

## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

pl.save(outdir+'crossec.png', transparent=False, figtight=True)

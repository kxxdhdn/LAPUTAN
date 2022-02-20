#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Galaxy SED (Galliano et al. 2018)

"""
import sys, os
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot

## Print full numpy array
np.set_printoptions(threshold=sys.maxsize)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/lib/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

## INPUT
##-------
## dataset: 
## /Inupol, /fexc, /labline, /nuLnuAME,
## /nuLnuDIBs, /nuLnuERE, /nuLnuMS, /nuLnuOB,
## /nuLnudust, /nuLnuescgas, /nuLnuescstar,
## /nuLnuff, /nuLnusync, /nuLnutot, /polem,
## /polext, /sigERE, /sigpol, /vLvline, /w,
## /wline, /wmaxERE, /wpol
## used: 
## /labline, /nuLnuDIBs, /nuLnuMS, /nuLnuOB,
## /nuLnudust, /nuLnuescgas, /nuLnuescstar,
## /nuLnuff, /nuLnusync, /vLvline, /w, /wline
filin = datdir+'sed_ARAA'
w = read_hdf5(filin, 'w')
total = read_hdf5(filin, 'nuLnutot')
stars = read_hdf5(filin, 'nuLnuescstar')
gas = read_hdf5(filin, 'nuLnuescgas')
dust = read_hdf5(filin, 'nuLnudust')
labL = read_hdf5(filin, 'labline')
wavL = read_hdf5(filin, 'wline')
vLvline = read_hdf5(filin, 'vLvline')
Nline = len(wavL)
OB = read_hdf5(filin, 'nuLnuOB')
MS = read_hdf5(filin, 'nuLnuMS')
ff = read_hdf5(filin, 'nuLnuff')
sync = read_hdf5(filin, 'nuLnusync')
dibs = read_hdf5(filin, 'nuLnuDIBs')

## Extinction curve (D03)
filext = datdir+'Adust_D03'
lam = read_hdf5(filext, 'lambda (micron)')
Cext = read_hdf5(filext, 'C_extovH (cm^2ovH)')
lam0 = 0.5470
ind = closest(lam, lam0)
Cext0 = Cext[ind]
extCurve = interp1d(lam, Cext/Cext0,
                    kind='linear', fill_value='extrapolate')
## Extrapolate sync
Av = 10
for j, val in enumerate(sync):
    if val>0:
        ind = j
        break
syncFull = interp1d(w[ind:], sync[ind:]/np.exp( -extCurve(w[ind:])*Av/1.086 ),
                    kind='linear', fill_value='extrapolate')
sync0 = syncFull(w) * np.exp( -extCurve(w)*Av/1.086 )

## New total
lines = np.zeros(len(w))
for wi in wavL:
    ind = closest(w, wi)
    ## add lines to scgas
    lines[ind] = ff[ind] - gas[ind]
# print(ff-gas-lines)
total = dust + stars + gas + lines + sync0
## Ionization end
for j, g in enumerate(gas):
	if g>1:
		wion = w[j]
		print('Ionization end: ', w[j])
		break


## PLOT
##------
p = pplot(w, total, figsize=(15,10),
          c='k', lw=5, zorder=-.1,
          xlog=True, ylog=True, xlim=(1e-2,9e4), ylim=(5e2,8e9),
          xtkform='log_sci', ytkform='log_sci',
          xlabel=r'$\rm Wavelength, \lambda\ [\mu m]$',
          ylabel=r'$\rm Monochromatic luminosity, \nu L_\nu\ [L_\odot]$',
          tksize=20, xysize=20)
p.trans = 'Axes'
lw1 = 2
plt.rcParams['hatch.linewidth'] = lw1
## Ionizing photons
p.ax.fill_betweenx([0,8e9], wion, color='lightgrey', zorder=-1)
p.ax.text(1.2e-2, 3e9, 'Ionizing \nphotons', 
          c='dimgrey', fontsize=18, zorder=-1)
## dust
p.ax.fill_between(w, dust, color='lavenderblush', zorder=.1,
                  edgecolor='m', lw=lw1)
## stars (OB+MS)
p.ax.fill_between(w, stars, color='beige', zorder=.2)
p.ax.fill_between(w, OB+MS, stars, zorder=.2,
                  edgecolor='orange', lw=lw1, facecolor='None', hatch='--')
## gas (ff)
p.ax.fill_between(w, gas, color='lightcyan', zorder=.3)
p.ax.fill_between(w, ff, gas, zorder=.3,
                  edgecolor='darkcyan', lw=lw1, facecolor='None', hatch='--')
## line labels
for i in range(Nline):
    ind = closest(w, wavL[i])
    p.ax.text(wavL[i]*.93, total[ind]*1.5, labL[i], 
              c='darkcyan', fontsize=10, rotation=90)
## synchrotron
p.add_plot(w, sync0, c='limegreen', lw=lw1, zorder=.4)
## dibs
i1 = closest(w, 0.4)
i2 = closest(w, 1.4)
p.add_plot(w[i1:i2], dibs[i1:i2], c='b', lw=lw1, zorder=.4)
## labels
p.ax.text(5e1, 5e6, 'Dust', c='m', fontsize=35)
p.ax.text(7e-1, 1e7, 'Stars', c='orange', fontsize=35)
p.ax.text(5e-1, 5e4, 'Gas', c='darkcyan', fontsize=35)
p.ax.text(4e2, 5e8, 'Thermal \ndust emission', 
          c='m', fontsize=18)
p.ax.annotate('', xytext=(7e2,4e8),xy=(3e2,9e7),
              arrowprops=dict(fc='m',ec='m',width=lw1))
p.ax.text(2e-2, 1e9, 'Absorption \n& scattering', 
          c='orange', fontsize=18)
p.ax.annotate('', xytext=(5e-2,8e8),xy=(2e-1,2e8),
              arrowprops=dict(fc='orange',ec='orange',width=lw1))
p.ax.text(1e1, 1e4, 'Free-free', 
          c='darkcyan', fontsize=18)
p.ax.annotate('', xytext=(2e1,2e4),xy=(1e2,8e4),
              arrowprops=dict(fc='darkcyan',ec='darkcyan',width=lw1))
p.ax.text(4e1, 5e3, 'Escaping synchrotron', 
          c='limegreen', fontsize=18)
p.ax.annotate('', xytext=(2e2,4e3),xy=(5e2,2e3),
              arrowprops=dict(fc='limegreen',ec='limegreen',width=lw1))
p.ax.text(2e-1, 4e9, 'DIBs', 
          c='b', fontsize=18)
p.ax.annotate('', xytext=(2.5e-1,3e9),xy=(3.5e-1,5e8),
              arrowprops=dict(fc='b',ec='b',width=lw1))
p.ax.text(1e-1, 6e8, r'2175 $\rm \AA$', c='k', fontsize=14)
p.ax.vlines(0.2175, 1.5e8, 5.5e8, color='k', lw=.8)
p.ax.text(3.3, 5e9, 'Aromatic features', c='k', fontsize=14)
p.ax.hlines(4e9, 3.3, 17.0, color='k', lw=.8)
p.ax.vlines(3.3, 2.5e9, 4e9, color='k', lw=.8)
p.ax.vlines(6.2, 2.5e9, 4e9, color='k', lw=.8)
p.ax.vlines(7.7, 2.5e9, 4e9, color='k', lw=.8)
p.ax.vlines(8.6, 2.5e9, 4e9, color='k', lw=.8)
p.ax.vlines(11.3, 2.5e9, 4e9, color='k', lw=.8)
p.ax.vlines(12.7, 2.5e9, 4e9, color='k', lw=.8)
p.ax.vlines(17.0, 2.5e9, 4e9, color='k', lw=.8)
p.ax.text(1.5, 1.5e9, 'Aliphatic feature', c='k', fontsize=14)
p.ax.vlines(3.4, 2e8, 1.2e9, color='k', lw=.8)
p.ax.text(4, 9e8, 'Silicate features', c='k', fontsize=14)
p.ax.hlines(7e8, 9.7, 18.0, color='k', lw=.8)
p.ax.vlines(9.7, 5e8, 7e8, color='k', lw=.8)
p.ax.vlines(18.0, 5e8, 7e8, color='k', lw=.8)
## aesthetics
for spine in p.ax.spines.items():
    p.ax.spines[spine[0]].set_linewidth(2) # border width
p.ax.xaxis.set_ticks_position('both') # ticks
p.ax.yaxis.set_ticks_position('both')
p.ax.tick_params(axis='both', which='both', direction='in')
p.ax.tick_params(which='major', length=8, width=2)
p.ax.tick_params(which='minor', length=4, width=2)
# p.ax.set_axisbelow(False)
# for k, spine in p.ax.spines.items():  #ax.spines is a dictionary
#     spine.set_zorder(100)

p.save(outdir+'sed.png', transparent=False, figtight=True)

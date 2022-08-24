#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Dielectric functions of harmonic oscillator/a-C(:H)/DH21 Astrodust

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
## w(um)        n           k
dieps_BE = ascii.read(datdir+'zubko96_be'+'.dat')
w_BE = dieps_BE['w(um)']
n_BE = dieps_BE['n']
k_BE = dieps_BE['k']
Rem1_BE = n_BE**2 - k_BE**2 - 1 # eps1 = n^2 - k^2
Im_BE = 2*n_BE*k_BE # eps2 = 2nk
## dataset:
## E(eV)     Re(n)-1     Im(n)    Re(eps)-1   Im(eps)
dieps_Ad = ascii.read(datdir+'DH21Ad/'+'index_DH21Ad_P0.20_0.00_0.500'+'.txt')
w_Ad = h * c / dieps_Ad['E(eV)'] * 1e6
Rem1_Ad = dieps_Ad['Re(eps)-1']
Im_Ad = dieps_Ad['Im(eps)']
## electrons in dielectric
omega = np.linspace(0,3,100) # omega0
omega_p = 1 # omega0
gamma = 0.2 # omega0
Rem1_e = omega_p**2 * (1 - omega**2) / ((1 - omega**2)**2 + (gamma*omega)**2)
Im_e = omega_p**2 * gamma * omega / ((1 - omega**2)**2 + (gamma*omega)**2)
        
## PLOT
##------
pl = plotool(1, 3, figsize=(16,4))
pl.trans = 'Axes'
lw1 = 2
pl.set_fig(hspace=.5)
## e- in dielectric
pl.set_ax(xlim=(2,0),
          xlabel=r'Frequency, $\rm \omega/\omega_0$',
          ylabel=r'Dielectric function, $\epsilon_r(\omega)$',)
pl.plot(omega, Rem1_e, c='c', label=r'$\epsilon_1-1$')
pl.plot(omega, Im_e, c='m', label=r'$\epsilon_2$')
pl.ax.axhline(0, c='grey', ls='dashed')
pl.ax.axvline(1, c='grey')
pl.ax.text(0.05, 0.9, '(a)', style='italic', c='grey', size=15,
           transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='lower right')
## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

## a-C(:H) BE
pl.set_ax((0,1), xlog=True, ylog=True, xlim=(5e-2,8e2), ylim=(2e-2,2e2),
          xtkform='mylog', ytkform='mylog',
          xlabel=r'Wavelength, $\rm \lambda\ [\mu m]$',
          ylabel=r'Dielectric function, $\epsilon_r(\lambda)$',)
pl.plot(w_BE, Rem1_BE, c='c', label=r'$\epsilon_1-1$', subpos=(0,1))
pl.plot(w_BE, Im_BE, c='m', label=r'$\epsilon_2$', subpos=(0,1))
pl.ax.text(0.05, 0.9, '(b)', style='italic', c='grey', size=15,
           transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='lower right')
## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
# pl.ax.xaxis.set_major_formatter(
#     FuncFormatter(
#         lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)
#     )
# )
# pl.ax.yaxis.set_major_formatter(
#     FuncFormatter(
#         lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)
#     )
# )
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

## Astrodust
pl.set_ax((0,2), xlog=True, ylog=True, xlim=(5e-2,8e2), ylim=(2e-2,2e2),
          xtkform='mylog', ytkform='mylog',
          xlabel=r'Wavelength, $\rm \lambda\ [\mu m]$',
          ylabel=r'Dielectric function, $\epsilon_r(\lambda)$',)
pl.plot(w_Ad, Rem1_Ad, c='c', label=r'$\epsilon_1-1$', subpos=(0,2))
pl.plot(w_Ad, Im_Ad, c='m', label=r'$\epsilon_2$', subpos=(0,2))
pl.ax.text(0.05, 0.9, '(c)', style='italic', c='grey', size=15,
           transform=pl.ax.transAxes)
pl.append_handles()
pl.set_legend(loc='lower right')
## aesthetics
for spine in pl.ax.spines.items():
    pl.ax.spines[spine[0]].set_linewidth(2) # border width
pl.ax.xaxis.set_ticks_position('both') # ticks
pl.ax.yaxis.set_ticks_position('both')
# pl.ax.xaxis.set_major_formatter(
#     FuncFormatter(
#         lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)
#     )
# )
# pl.ax.yaxis.set_major_formatter(
#     FuncFormatter(
#         lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)
#     )
# )
pl.ax.tick_params(axis='both', which='both', direction='in')
pl.ax.tick_params(which='major', length=8, width=2)
pl.ax.tick_params(which='minor', length=4, width=2)

pl.save(outdir+'dieps.png', transparent=False, figtight=True)

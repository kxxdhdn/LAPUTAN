#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

MIR spectra (Galliano et al. 2018)

"""
import sys, os
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import ascii
import matplotlib.transforms as mtransforms

from rapyuta.arrays import closest
from rapyuta.inout import read_hdf5
from rapyuta.maths import gaussian
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
## FnuNGC1097, FnuNGC1569, FnucenA
## dFnuNGC1097, dFnuNGC1569, dFnucenA
## wNGC1097, wNGC1569, wcenA, rad, wmod
## LnuPAHi_max, LnuPAHi_min, LnuPAHn_max, LnuPAHn_min
filin = datdir+'specMIR_ARAA'
wNGC1097 = read_hdf5(filin, 'wNGC1097')
wNGC1569 = read_hdf5(filin, 'wNGC1569')
wcenA = read_hdf5(filin, 'wcenA')
FnuNGC1097 = read_hdf5(filin, 'FnuNGC1097')
FnuNGC1569 = read_hdf5(filin, 'FnuNGC1569')
FnucenA = read_hdf5(filin, 'FnucenA')
dFnuNGC1097 = read_hdf5(filin, 'dFnuNGC1097')
dFnuNGC1569 = read_hdf5(filin, 'dFnuNGC1569')
dFnucenA = read_hdf5(filin, 'dFnucenA')

for i, x in enumerate(FnuNGC1097):
    if x==0:
        FnuNGC1097[i] = np.nan

for i, x in enumerate(FnuNGC1569):
    if x==0:
        FnuNGC1569[i] = np.nan

for i, x in enumerate(FnucenA):
    if x==0:
        FnucenA[i] = np.nan

## Adust_D03
w_D03 = read_hdf5(datdir+'Adust_D03', 'lambda (micron)')
C0 = read_hdf5(datdir+'Adust_D03', 'C_extovH (cm^2ovH)')
iV = closest(w_D03, 0.5470) # normalized to V band
ext_D03 = C0/C0[iV]
Av = 10
Pabs = np.exp(-Av/1.086 * ext_D03)
# wabs = w_D03
wabs = np.linspace(2,32,1000)
fabs = interp1d(w_D03, Pabs, kind='linear', fill_value='extrapolate')
Pabs = fabs(wabs)

## ice/silicate data
# ice_CO = ascii.read(datdir+'ice_CO_15K'+'.dat')
# w_CO = ice_CO['col1']
# Av_CO = ice_CO['col2']
# print(w_CO)
# print(Av_CO)


## PLOT
##------
p = pplot(figsize=(20,12),
          # c='k', lw=5, zorder=-.1,
          xlog=True, ylog=True, xlim=(2.5,32), ylim=(1e-1,2e3),
          xtk=[3,5,7,10,15,20,30], ytkmi=np.logspace(-1,3,40), ytkform='mylog',
          xlabel=r'Wavelength, $\rm \lambda\ [\mu m]$',
          ylabel=r'Flux, $\rm F_\nu\ [MJy/sr]$',
          tksize=20, xysize=20)
p.trans = 'Axes'
lw1 = 2
## Cen A
p.add_plot(wcenA, FnucenA, dFnucenA,
           c='m', lw=lw1, ec='pink',
           label='Cen A')
ymax = FnucenA
wmax = wcenA
## NGC 1097
p.add_plot(wNGC1097, FnuNGC1097, dFnuNGC1097,
           c='orange', lw=lw1, ec='beige',
           label='NGC 1097')
## NGC 1569
p.add_plot(wNGC1569, FnuNGC1569/3, dFnuNGC1569/3,
           c='darkcyan', lw=lw1, ec='lightcyan',
           label='NGC 1569 \n'+r'($\rm F_\nu/3$)')
ymin = FnuNGC1569/3
wmin = wNGC1569

## lines
labL = [r'Br$\rm \alpha$', r'[ArII]', r'Pf$\rm \alpha$', r'[ArIII]',
        r'0-0 S(3)', r'[SIV]', r'0-0 S(2)', r'[NeII]',
        r'[NeIII]', r'0-0 S(1)', r'[SIII]', r'[OIV]', ]
wavL = [4.052, 6.985274, 7.459858, 8.99103,
        9.66491, 10.5105, 12.27861, 12.81355,
        15.555, 17.03483, 18.7129, 25.8903, ]
for i, w in enumerate(wavL):
    p.ax.vlines(w, ymin[closest(wmin,w)], 7e2, color='limegreen', lw=lw1, ls='--')
    p.ax.text(w*.985, 8e2, labL[i], c='limegreen', fontsize=15, rotation=90)

## bands
# p.ax.annotate('', xytext=(3.3,2e3), xy=(3.3,3e-1),
#               arrowprops=dict(arrowstyle='-|>',mutation_scale=30,
#                               lw=lw1,ls='--',color='b'))
p.ax.vlines(3.3, 3e-1, ymax[closest(wmax,3.3)], color='b', lw=lw1, ls='--') # 3.3
p.ax.text(3.2, 1.2e-1, 'C-H \nstretch \n3.3 '+r'$\rm \mu$m',
          c='b', fontsize=15)
# p.ax.annotate('', xytext=(3.4,1e-1), xy=(3.4,2e2),
#               arrowprops=dict(arrowstyle='-|>',mutation_scale=30,
#                               lw=lw1,ls='--',color='b'))
p.ax.vlines(3.4, ymin[closest(wmin,3.4)], 2e2, color='b', lw=lw1, ls='--') # 3.4
p.ax.text(3.3, 2.5e2, 'Aliphatic \nC-H \nstretch \n3.4 '+r'$\rm \mu$m',
          c='b', fontsize=15)
p.ax.vlines(6.2, 3e-1, ymax[closest(wmax,6.2)], color='b', lw=lw1, ls='--') # 6.2
p.ax.vlines(7.6, 3e-1, ymax[closest(wmax,7.6)], color='b', lw=lw1, ls='--') # 7.6
p.ax.hlines(3e-1, 6.2, 7.6, color='b', lw=lw1, ls='--')
p.ax.text(6.2, 1.2e-1, 'C-C \nstretch \n6.2, 7.6 '+r'$\rm \mu$m',
          c='b', fontsize=15)
p.ax.vlines(7.8, 3e-1, ymax[closest(wmax,7.8)], color='b', lw=lw1, ls='--') # 7.8
p.ax.vlines(8.6, 3e-1, ymax[closest(wmax,8.6)], color='b', lw=lw1, ls='--') # 8.6
p.ax.hlines(3e-1, 7.8, 8.6, color='b', lw=lw1, ls='--')
p.ax.text(7.8, 1.2e-1, 'C-H \nin-plane bending \n7.8, 8.6 '+r'$\rm \mu$m',
          c='b', fontsize=15)
p.ax.vlines(11.3, 3e-1, ymax[closest(wmax,11.3)], color='b', lw=lw1, ls='--') # 11.3
p.ax.vlines(12.7, 3e-1, ymax[closest(wmax,12.7)], color='b', lw=lw1, ls='--') # 12.7
p.ax.hlines(3e-1, 11.3, 12.7, color='b', lw=lw1, ls='--')
p.ax.text(11, 1.2e-1, 'C-H \nout-of-plane bending \n11.3, 12.7 '+r'$\rm \mu$m',
          c='b', fontsize=15)
# p.ax.annotate('', xytext=(17.0,2e3), xy=(17.0,3e-1),
#               arrowprops=dict(arrowstyle='-|>',mutation_scale=30,
#                               lw=lw1,ls='--',color='b'))
p.ax.vlines(17.0, 3e-1, ymax[closest(wmax,17.0)], color='b', lw=lw1, ls='--') # 17.0
p.ax.text(17, 1.2e-1, 'C-C-C \nbending \n17.0 '+r'$\rm \mu$m',
          c='b', fontsize=15)

## absorption (Gaussian)
trans = mtransforms.blended_transform_factory(p.ax.transData, p.ax.transAxes)
labE1 = [r'$\rm H_2$O'+' ice \n3.05 '+r'$\rm \mu$m',
         'CO ice \n4.7 '+r'$\rm \mu$m',
         r'$\rm H_2$O'+' ice \n6.02 '+r'$\rm \mu$m',
         r'C$\rm O_2$'+' ice \n15.2 '+r'$\rm \mu$m', ]
cenE1 = [3.05, 4.67, 6.02, 15.2]
sigE1 = [.27, .04, .49, .18] # Oberg2007, Fraser2004, Bisschop2007
for i in range(len(labE1)):
    Ngrad = 100 # color gradient number (3-sigma range)
    for j in range(Ngrad):
        left_edge = -3*sigE1[i]+6*sigE1[i]*j/Ngrad
        alpha_gauss = 0.5 * gaussian(left_edge,0,sigE1[i])/gaussian(0,0,sigE1[i])
        p.ax.axvspan(cenE1[i]+left_edge, cenE1[i]+left_edge+6*sigE1[i]/Ngrad,
                     color='deepskyblue', ec='None', alpha=alpha_gauss)
    p.ax.text(cenE1[i]*.96, 2.3e3, labE1[i], c='deepskyblue', fontsize=15)
labE2 = ['Silicates \n9.8 '+r'$\rm \mu$m',
         'Silicates \n18 '+r'$\rm \mu$m', ]
cenE2 = [9.8, 18.0]
sigE2 = [.7, .5] #
ampE2 = [2, 1]
for i in range(len(labE2)):
    Ngrad = 100 # color gradient number (3-sigma range)
    for j in range(Ngrad):
        left_edge = -3*sigE2[i]+6*sigE2[i]*j/Ngrad
        alpha_gauss = ampE2[i] * 0.5 * gaussian(left_edge,0,sigE2[i])/gaussian(0,0,sigE2[i])
        p.ax.axvspan(cenE2[i]+left_edge, cenE2[i]+left_edge+6*sigE2[i]/Ngrad,
                     color='grey', ec='None', alpha=alpha_gauss)
    p.ax.text(cenE2[i]*.96, 2.3e3, labE2[i], c='grey', fontsize=15)

## absorption (extinction curve)
# trans = mtransforms.blended_transform_factory(p.ax.transData, p.ax.transAxes)
# labE = [
#         # r'$\rm H_2$O'+' ice \n3.1 '+r'$\rm \mu$m',
#         # 'CO ice \n4.7 '+r'$\rm \mu$m',
#         # r'$\rm H_2$O'+' ice \n6.3 '+r'$\rm \mu$m',
#         # r'C$\rm O_2$'+' ice \n15.2 '+r'$\rm \mu$m',
#         'Silicates \n9.8 '+r'$\rm \mu$m',
#         'Silicates \n18 '+r'$\rm \mu$m', ]
# cenE = [
#         # 3.1, 4.7, 6.3, 15.2,
#         9.8, 18.0]
# Ngrad = len(wabs)
# for j in range(Ngrad-1):
#     for i in range(len(labE)):
#         if labE[i][:2]=='Si':
#             colorabs = 'grey'
#             if labE[i][:14]=='Silicates \n9.8':
#                 sigE = 2
#             else:
#                 sigE = 0.5
#         else:
#             colorabs = 'deepskyblue'
#             sigE = .3
#         p.ax.text(cenE[i]*.96, 2.3e3, labE[i], c=colorabs, fontsize=15)
#         imin = closest(wabs, cenE[i]-sigE, 'left')
#         imax = closest(wabs, cenE[i]+sigE, 'right')
#         amin = min(1-Pabs[imin:imax])
#         amax = max(1-Pabs[imin:imax])
#         if wabs[j]>cenE[i]-sigE and wabs[j]<cenE[i]+sigE:
#             alphabs = 0.5 * (1-Pabs[j] - amin)/(amax-amin)
#             print(amin, amax, alphabs)
#             p.ax.axvspan(wabs[j], wabs[j+1],
#                          color=colorabs, ec='None', alpha=alphabs)

## aesthetics
for spine in p.ax.spines.items():
    p.ax.spines[spine[0]].set_linewidth(2) # border width
p.ax.xaxis.set_ticks_position('both') # ticks
p.ax.yaxis.set_ticks_position('both')
p.ax.tick_params(axis='both', which='both', direction='in')
p.ax.tick_params(which='major', length=8, width=2)
p.ax.tick_params(which='minor', length=4, width=2)
p.ax.legend(loc='lower right',# bbox_to_anchor=(1,1),
            fontsize=20, framealpha=0,)

p.save(outdir+'specMIR.png', transparent=False, figtight=True)

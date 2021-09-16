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

## Local
sys.path.insert(0, testdir+'/..') ## rapyuta path
from rapyuta.inout import read_hdf5
from rapyuta.plots import plotool, pplot

## Data to be plotted
##--------------------

## Read the data file containing some uncertainties
ref = read_hdf5(datdir+"plotXY_D2G",name="Sample")
mulogOovH = read_hdf5(datdir+"plotXY_D2G", "Mean of 12+log(OovH)")
siglogOovH = read_hdf5(datdir+"plotXY_D2G", "Sigma of 12+log(OovH)")
gamlogOovH = read_hdf5(datdir+"plotXY_D2G", "Skewness of 12+log(OovH)")
mulnD2G = read_hdf5(datdir+"plotXY_D2G", "Mean of ln(dust-to-gas ratio)")
siglnD2G = read_hdf5(datdir+"plotXY_D2G", "Sigma of ln(dust-to-gas ratio)")
gamlnD2G = read_hdf5(datdir+"plotXY_D2G", "Skewness of ln(dust-to-gas ratio)")
rho = read_hdf5(datdir+"plotXY_D2G", "Correlation of 12+log(OovH) and ln(dust-to-gas ratio)")
mask = np.isfinite(mulnD2G)
mux0 = mulogOovH[mask]
sigx0 = siglogOovH[mask]
gamx0 = gamlogOovH[mask]
mulny0 = mulnD2G[mask]
siglny0 = siglnD2G[mask]
gamlny0 = gamlnD2G[mask]
rho0 = rho[mask]
Nsamp = np.size(mulogOovH)

print('\n TEST ellipses ')
print('---------------')
ep = plotool()
ep.figure(figsize=(9,5),nrows=1,ncols=2)
ep.set_fig(right=.99, wspace=.5)

## Ellipses
ep.set_ax(ylog=1,
          xlabel=r'Metallicity, $12+\log({\rm O/H})$',
          ylabel=r'Dust-to-HI mass ratio $M_{\rm dust}$/$M_{\rm HI}$')
ep.plot(mux0, mulny0, fmt='s', label='Symbol', yisln=True, 
        c='c', markersize=5, marker='*')
ep.eplot(mux0, mulny0, sigmax=sigx0, sigmay=siglny0, rho=rho0,
         yisln=True, ec='g',errinlegend='Error')
ep.set_legend(loc='upper left')

## Skewllipses
ep.set_ax((1,2), ylog=1,
          xlabel=r'Metallicity, $12+\log({\rm O/H})$',
          ylabel=r'Dust-to-HI mass ratio $M_{\rm dust}$/$M_{\rm HI}$')
ep.plot(mux0, mulny0, fmt='s', label='Symbol', yisln=True, 
        c='r', markersize=5, marker='*', alpha=.5, subpos=(1,2))
ep.eplot(mux0, mulny0, sigmax=sigx0, sigmay=siglny0, rho=rho0,
         gammax=gamx0, gammay=gamlny0,
         ec='m', elw=.6, els='dashed', efill=True, efillcolor='pink',
         yisln=True, errinlegend='Error', alpha=.5, subpos=(1,2))
ep.set_legend((1,2), loc='upper left')

ep.save(outdir+'ellipses')
print('See out/ellipses.png [Done]')


print('\n TEST pplot ')
print('------------')
## Ellipses (half)
pp = pplot(mux0[:int(Nsamp/2)], mulny0[:int(Nsamp/2)],
           sigmax=sigx0[:int(Nsamp/2)], sigmay=siglny0[:int(Nsamp/2)], rho=rho0[:int(Nsamp/2)],
           yisln=True, ec='g',errinlegend='Error I',
           figsize=(8,8), ylog=1, legend='center left', 
           right=.8, anchor=(1,.5),
           xlabel=r'Metallicity, $12+\log({\rm O/H})$',
           ylabel=r'Dust-to-HI mass ratio $M_{\rm dust}$/$M_{\rm HI}$',
           fmt='s', label='Symbol I', c='c', markersize=5, marker='*')

## Skewllipses (half)
pp.add_plot(mux0[int(Nsamp/2):], mulny0[int(Nsamp/2):],
            sigmax=sigx0[int(Nsamp/2):], sigmay=siglny0[int(Nsamp/2):], rho=rho0[int(Nsamp/2):],
            gammax=gamx0[int(Nsamp/2):], gammay=gamlny0[int(Nsamp/2):],
            ec='m', elw=.6, els='dashed', efill=True, efillcolor='pink',
            yisln=True, errinlegend='Error II', alpha=.5,
            fmt='s', label='Symbol II', c='r', markersize=5, marker='*')

pp.save(outdir+'test_pplot')
print('See out/test_pplot.png [Done]')

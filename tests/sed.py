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
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

## rapyuta
from rapyuta.inout import read_hdf5
from rapyuta.plots import pplot


## Input/output
filIN = datdir+'dustobs'
filOUT = outdir+'sed_shift.png'

## Read data
nuLnutot = read_hdf5(filIN, 'nuLnutot')
w0 = read_hdf5(filIN, 'w')

## Set redshift
# z = [0,.3,1,2,3,5,7,12]
z = [1,3,6]
colors = ['k','c','g','maroon']

p = pplot(w0, nuLnutot, label='Rest frame',
          xlim=(2e-2, 1e6), ylim=(1e1, 8e12),
          xlog=1, ylog=1, clib=colors, lw=2,
          xlabel=r'$\rm{Wavelength,}\ \lambda\ [\mu m]$',
          ylabel=r'$\rm{Monochromatic\ Luminosity,}\ \nu L_{\nu}\ [L_{\odot}]$',
          title=None,
          figsize=(12,6), left=.1, bottom=.15, top=.95, right=.95,
          titlesize=20, labelsize=20, ticksize=20)
p.ax.text(0.2, 2e12, r'$[H_{\alpha}]_{6563\AA}$',size=20,c='k')
p.ax.annotate('', xytext=(0.6563,1e12),xy=(0.6563,2e11),
              arrowprops=dict(facecolor='k',shrink=0,width=.1,
                              headwidth=5,headlength=10))
for i in range(len(z)):
    # p.add_plot(w0*(1+z[i]), nuLnutot, label='z='+str(z[i]))
    p.add_plot(w0*(1+z[i]), nuLnutot/(1+z[i])**2, label='z='+str(z[i]))

    lam = 0.6563*(1+z[i])
    if i==0:
      ytext = 3e10
      ybot = 2e10
      yupp = 4e9
    elif i==1:
      ytext = 3e9
      ybot = 2e9
      yupp = 4e8
    elif i==2:
      ytext = 2e12
      ybot = 1e12
      yupp = 2e11
    p.ax.text(lam*.5, ytext, str(lam)+r'$\mu m$',size=20,c=colors[i+1])
    p.ax.annotate('', xytext=(lam,ybot),xy=(lam,yupp),
                  arrowprops=dict(edgecolor=colors[i+1],facecolor=colors[i+1],
                                  shrink=0,width=.1,
                                  headwidth=5,headlength=10))

trans = mtransforms.blended_transform_factory(p.ax.transData, p.ax.transAxes)
# MRSwmin, MRSwmax = 4.88, 28.34
# PACSwmin, PACSwmax = 55, 210
# ALMAwmin, ALMAwmax = .3e3, 3.6e3
# MRSf = (w0-MRSwmin)*(w0-MRSwmax)
# ALMAf = (w0-ALMAwmin)*(w0-ALMAwmax)
# p.ax.fill_between(w0, 0, 1, where=MRSf<0,
#                   facecolor='b', alpha=0.2, transform=trans)
# p.ax.text(MRSwmin+1, 2.e10, 'JWST MIRI MRS')
# p.ax.hlines(5.e10, PACSwmin, PACSwmax)
# p.ax.vlines(PACSwmin, 4.e10, 6.e10)
# p.ax.vlines(PACSwmax, 4.e10, 6.e10, )
# p.ax.arrow(PACSwmax, 5.e10, PACSwmin-PACSwmax, 0)
# p.ax.text(PACSwmin, 2.e10, 'Herschel PACS')

# p.ax.fill_between(w0, 0, 1, where=ALMAf<0,
#                   facecolor='r', alpha=0.2, transform=trans)
# p.ax.text(ALMAwmin+300, 2.e10, 'ALMA bands')

## NIR
IRf = (w0-0.7)*(w0-2.5)
p.ax.fill_between(w0, 0, 1, where=IRf<0,
                  facecolor='r', alpha=0.2, transform=trans)
## MIR
IRf = (w0-2.5)*(w0-25)
p.ax.fill_between(w0, 0, 1, where=IRf<0,
                  facecolor='r', alpha=0.4, transform=trans)
## FIR
IRf = (w0-25)*(w0-300)
p.ax.fill_between(w0, 0, 1, where=IRf<0,
                  facecolor='r', alpha=0.6, transform=trans)
p.ax.legend(loc='upper right',# bbox_to_anchor=(1,1),
            fontsize=20, framealpha=0,)
p.save(filOUT, transparent=True)

# p.show()

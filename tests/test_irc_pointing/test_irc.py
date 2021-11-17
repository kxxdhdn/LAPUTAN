#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from rapyuta.inout import read_fits
from rapyuta.plots import pplot

slit = 'Ns'

## raw
name = 'raw'
ds = read_fits('3390001.1_'+slit,'3390001.1_'+slit+'_unc')
snr = np.mean(ds.data/ds.unc)
# p = pplot(ds.wave, ds.data[:,14,1], yerr=ds.unc[:,14,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
# p.save(name)
mask = (ds.data/ds.unc)[:,14,1] >0.1
p = pplot(ds.wave[mask], (ds.data/ds.unc)[:,14,1][mask], 
    figsize=((10,8)), left=.15, right=.99, bottom=.1, top=.98,
	c='k',label=name,ylog=1,#ylim=(5e1,2e4),
	title=None,
	xlabel=r'$\lambda\,(\mu m)$', ylabel='Signal-to-noise ratio',
    titlesize=20, labelsize=20, ticksize=20)

## pt3
name = 'pt3'
name = '2err'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
# p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
# p.save(name)
mask = (ds.data/ds.unc)[:,14,1] >0.1
p.add_plot(ds.wave[mask], (ds.data/ds.unc)[:,14,1][mask], c='c',label='pointing error')

p.ax.legend(loc='upper left', fontsize=20, framealpha=0)
p.save('test_irc_pt')
exit()

## splitnorm
name = 'splitnorm'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

## 2err
name = '2err'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

## 2errsup
name = '2errsup'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

## 2errfill
name = '2errfill'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

## 2errfillsup
name = '2errfillsup'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

## 2errswp
name = '2errswp'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)
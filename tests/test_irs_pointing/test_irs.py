#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from rapyuta.inout import read_fits
from rapyuta.plots import pplot

slit = 'Ns'

## raw
name = 'raw'
ds = read_fits('../lib/M83','../lib/M83_unc')
snr = np.mean(ds.data/ds.unc)
# p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
# p.save(name)
mask = (ds.data/ds.unc)[:,6,5] < 5000
p = pplot(ds.wave[mask], (ds.data/ds.unc)[:,6,5][mask], 
    figsize=((10,8)), left=.15, right=.99, bottom=.1, top=.98,
	c='k',label=name,ylog=1,ylim=(5e1,2e4),title=None,
	xlabel=r'$\lambda\,(\mu m)$', ylabel='Signal-to-noise ratio',
    titlesize=20, labelsize=20, ticksize=20)

## 0_2
name = 'pointing error'
ds = read_fits('irs_pt_0_2','irs_pt_0_2_unc')
snr = np.mean(ds.data/ds.unc)
# p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
# p.save(name)
p.add_plot(ds.wave, (ds.data/ds.unc)[:,6,5], c='c',label=name)

# ## norm
# name = 'random'
# ds = read_fits('irs_pt_norm','irs_pt_norm_unc')
# snr = np.mean(ds.data/ds.unc)
# # p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
# # p.save(name)
# mask = (ds.data/ds.unc)[:,6,5] < 5000
# p.add_plot(ds.wave[mask], (ds.data/ds.unc)[:,6,5][mask], c='m',label=name)

## 2err
# name = '2 errors'
# ds = read_fits('irs_pt_2err','irs_pt_2err_unc')
# snr = np.mean(ds.data/ds.unc)
# # p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
# # p.save(name)
# p.add_plot(ds.wave, (ds.data/ds.unc)[:,6,5], c='m',label=name)

p.ax.legend(loc='upper left', fontsize=20, framealpha=0)
p.save('test_irs_pt')
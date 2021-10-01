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
p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
p.save(name)

## 0_2
name = '0_2'
ds = read_fits('irs_pt_'+name,'irs_pt_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
p.save(name)

## 2err
name = 'norm'
ds = read_fits('irs_pt_'+name,'irs_pt_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
p.save(name)

## 2err
name = '2err'
ds = read_fits('irs_pt_'+name,'irs_pt_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,6,5], yerr=ds.unc[:,6,5],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.nanmean(ds.data/ds.unc)))
p.save(name)


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
p = pplot(ds.wave, ds.data[:,14,1], yerr=ds.unc[:,14,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

## pt3
name = 'pt3'
ds = read_fits('3390001.1_'+slit+'_'+name,'3390001.1_'+slit+'_'+name+'_unc')
snr = np.mean(ds.data/ds.unc)
p = pplot(ds.wave, ds.data[:,20,1], yerr=ds.unc[:,20,1],ec='r',c='k',label=name,legend=2,title='S/R={}'.format(np.mean(ds.data/ds.unc)))
p.save(name)

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
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
from astropy.io import fits
from astropy.wcs import WCS

## Local
sys.path.insert(0, testdir+'/..') ## rapyuta path
from rapyuta.inout import (fclean, read_fits, write_fits,
                            read_hdf5, write_hdf5,
                            read_ascii, write_ascii,
                            read_csv, write_csv, )

print('\n TEST read_fits ')
print('----------------')
get_fdata = read_fits(datdir+'M83', datdir+'M83_unc')
# print('header \n', get_fdata.header)
# print('header of TAB \n', get_fdata.header_w)
print('data (shape) \n', get_fdata.data.shape)
print('wave (shape) \n', get_fdata.wave.shape)
print('uncertainty (shape) \n', get_fdata.unc.shape)

print('\n TEST write_fits ')
print('----------------')
write_fits(outdir+'M83_copy',
           get_fdata.header, get_fdata.data, get_fdata.wave)
write_fits(outdir+'M83_irs_data_format_copy',
           get_fdata.header, get_fdata.data, get_fdata.wave, wmod=1)

print('\n TEST write_hdf5 ')
print('-----------------')
label = [['alpha', 'beta', 'gamma'], [1,2,3]]
data1 = np.arange(5, 23, 2).reshape((3,3))
# data2 = [1,2,3]
# data2 = [[1,2,3]]
data2 = np.arange(4).reshape((2,2))
write_hdf5(outdir+'test_inout', 'Label', label)
write_hdf5(outdir+'test_inout', 'Data', data1, append=True)
# write_hdf5(outdir+'test_inout', 'Data', data2, ind1=[2,5], append=True)
# write_hdf5(outdir+'test_inout', 'Data', data2, ind1=1, ind2=[0,3], append=True)
write_hdf5(outdir+'test_inout', 'Data', data2, ind1=[0,2], ind2=[1,3], append=True)
write_hdf5(outdir+'test_inout_h5grp', 'Grp data', data1, group='/test_grp', verbose=True)
write_hdf5(outdir+'test_inout_h5grp', 'Grp data', data2, group='/test_grp', 
    ind1=[0,2], ind2=[1,3], append=True, verbose=True)
print('See ./out [Done]')

print('\n TEST read_hdf5 ')
print('----------------')
get_label = read_hdf5(outdir+'test_inout', 'Label')
print('Label \n', get_label)
get_data = read_hdf5(outdir+'test_inout', 'Data')
print('Data \n', get_data)
get_grp = read_hdf5(outdir+'test_inout_h5grp', 'Grp data', group='test_grp')
print('Grp data \n', get_grp)

print('\n TEST write_ascii ')
print('------------------')
lab = ['col1', 'col2']
arr = np.arange(12).reshape((6,2))
arr_trans = np.arange(0, 12000, 1000).reshape((2,6))
write_ascii(outdir+'test_inout', comment='Empty file!')
write_ascii(outdir+'test_inout', lab, arr, comment='Data begins from line 3')
write_ascii(outdir+'test_inout', dset=arr_trans, trans=True, append=True)
# write_ascii(outdir+'ascii_csv', lab, arr, ascext='.csv') # not compatible
print('See ./out [Done]')

print('\n TEST read_ascii ')
print('-----------------')
get_arr = read_ascii(outdir+'test_inout', dtype=float)
print('col1 \n', get_arr['col1'])
print('col2 \n', get_arr['col2'])
get_adata = read_ascii(datdir+'filt_IRAC1', start_header=2)
print('col1 - Wave (microns) \n', get_adata['Wave'][:10])
print('col2 - Spectral Response (electrons/photon) \n', get_adata['col2'][:10])

print('\n TEST write_csv ')
print('----------------')
write_csv(outdir+'test_inout', lab, arr)
print('See ./out [Done]')

print('\n TEST read_csv ')
print('---------------')
get_cdata = read_csv(outdir+'test_inout', 'col2', 'col1')
print('col2: \n', get_cdata['col2'])
print('col1: \n', get_cdata['col1'])

if input('Test fclean (y/n): ')=='y':
    print('\n TEST fclean ')
    print('-------------')
    info = './out is empty! [Done]'
    fclean('out', info)

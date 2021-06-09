#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))

import numpy as np
import h5py as H5
h5ext = '.h5'
import subprocess as SP

def read_hdf5(file, *header):
	'''
	Read h5 file

	------ INPUT ------
	file                input h5 file
	header              labels of data to read
	------ OUTPUT ------
	dataset             data
	'''
	hf = H5.File(file+h5ext, 'r')
	dataset = []
	for hdr in header:
		data = hf.get(hdr)
		data = np.array(data)
		dataset.append(data)

	hf.close()

	return dataset


SP.call('synthetic_photometry', shell=True)

## Read the output
fortOUT = os.getcwd()+'/synthetic_photometry_output'

wcen, Fnu_filt, smat = read_hdf5(fortOUT, \
    'Central wavelength (microns)', \
    'Flux (x.Hz-1)', \
    'Standard deviation matrix')

print('Fnu_filt = ', Fnu_filt)
print('wcen = ', wcen)
print('smat = ', smat)

SP.call('rm -rf '+fortOUT+'.h5', shell=True, cwd=os.getcwd())

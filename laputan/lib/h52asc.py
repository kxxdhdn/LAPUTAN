#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy as np
# from astropy.table import Table
from astropy.io import ascii

## astylo
datdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, datdir+'/../..')
from playtime.inout import read_hdf5, write_ascii
ascext = '.txt'

## Path
h5dir = './data/'
filt = ['IRAC1', 'IRAC2', 'IRAC3', 'IRAC4',
        'MIPS1', 'MIPS2', 'MIPS3',
        'WISE1', 'WISE2', 'WISE3', 'WISE4',]

for f in filt:
	## Read
	reader = read_hdf5(h5dir+'filt_'+f,
		      'Filter wavelength (microns)',
		      'Filter transmission')

	## Write (astropy.ascii.write)
	# data = Table([wave, tran], names=['Wave', 'Spectral Response'])
	# ascii.write(data, 'filt_'+f+ascext, format='commented_header')
	## Write (astylo.iolib.write_ascii)
	comment = 'Average spectral response curve (electrons/photon - microns) for '+f+' array'
	write_ascii('filt_'+f, 
		        header=['Wave', 'Spectral_Response'],
		        dset=reader, trans=True, comment=comment)


print(">>> Coucou h52asc [done] <<<")

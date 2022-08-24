#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pathlib import Path
import numpy as np
# from astropy.table import Table
from astropy.io import ascii

## Local
try:
    root = Path(__file__).parent.absolute() / '../../..'
except NameError:
    root = Path().absolute() / '../../..'
sys.path.insert(0, root)

from rapyuta.inout import read_hdf5, write_ascii, ascext

## Path
filtdir = str(root / 'rapyuta/lib/filt')
filt = ['IRAC1', 'IRAC2', 'IRAC3', 'IRAC4',
        'MIPS1', 'MIPS2', 'MIPS3',
        'WISE1', 'WISE2', 'WISE3', 'WISE4',]

for f in filt:
    ## Read
    col1 = read_hdf5('filt_'+f, 'Filter wavelength (microns)')
    col2 = read_hdf5('filt_'+f, 'Filter transmission')
    data = np.array([col1, col2])#.reshape((len(col1),2))
    # print(data.shape)

    ## Write (astropy.ascii.write)
    # data = Table([wave, tran], names=['Wave', 'Spectral Response'])
    # ascii.write(data, 'filt_'+f+ascext, format='commented_header')
    ## Write (rapyuta.inout.write_ascii)
    comment = 'Average spectral response curve (electrons/photon - microns) of '+f
    write_ascii('filt_'+f, header=['Wave', 'Spectral_Response'],
    	        dset=data, trans=True, comment=comment)


print(">>> Coucou h52asc [done] <<<")
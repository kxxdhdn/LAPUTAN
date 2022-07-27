#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Convert ASCII to HDF5

"""
import os
import numpy as np

## rapyuta
from rapyuta.inout import write_hdf5
ascext = '.txt'

## Path
mroot = os.path.dirname(os.path.abspath(__file__))+'/'
datdir = mroot
filin = mroot+'kext_albedo_WD_MW_3.1A_60_D03.all'
# filin = mroot+'kext_albedo_WD_MW_4.0A_40_D03.all'
# filin = mroot+'kext_albedo_WD_MW_5.5A_30_D03.all'
filout = datdir+'../Adust_D03'
# filout = datdir+'Adust_MW_3.1A_60_D03'
# filout = datdir+'Adust_MW_4.0A_40_D03'
# filout = datdir+'Adust_MW_5.5A_30_D03'

## Read
dset = np.genfromtxt(filin+ascext, skip_header=80,
                     dtype={'names': ('lambda','albedo','cos','C_extovH','K_abs','cos2'),
                            'formats':  ('f8','f8','f8','f8','f8','f8')},
                     usecols=(0,1,2,3,4,5))

## Write
write_hdf5(filout, 'lambda (micron)', dset['lambda'])
write_hdf5(filout, 'albedo', dset['albedo'], append=True)
write_hdf5(filout, '<cos>', dset['cos'], append=True)
write_hdf5(filout, 'C_extovH (cm^2ovH)', dset['C_extovH'], append=True)
write_hdf5(filout, 'K_abs (cm^2ovg)', dset['K_abs'], append=True)
write_hdf5(filout, '<cos^2>', dset['cos2'], append=True)

print(">>> Coucou asc2h5 [Done] <<<")

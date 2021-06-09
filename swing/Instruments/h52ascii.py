#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Set dir
h5dir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/dat/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
	os.makedirs(outdir)

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

## Local
sys.path.insert(0, testdir+'/..') ## astylo path
from astylo.calib import (intercalib, photometry_profile)


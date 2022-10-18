#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

SST data processing and analyzing

    

"""

import os
import math
import numpy as np
from astropy.io import ascii
import subprocess as SP
import warnings

## Local
import rapyuta.utbox as UT
import rapyuta.inout as IO
import rapyuta.impro as IM
from inout import fitsext, h5ext
from .utils import *

##------------------------------------------------
##
##            Spitzer/IRS data cube
##
##------------------------------------------------


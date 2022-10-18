#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

List, Array, Tuple, Table, Etc.

    ramp, arrayize, listize, closest, zoom

"""

import warnings
import math
import numpy as np
from astropy.table import Table

## Local
import rapyuta.utbox as UT


##------------------------------------------------
##
##           (Copyright: F. Galliano)
##
##------------------------------------------------

def ramp(x0=0,x1=1,dx=None,dlnx=None,dlogx=None,N=None,log=None, \
         homogenize=True):
    '''
    GENERATE A REAL RAMP OF VALUES (unlike arange)

    This function generates a ramp between x0 and x1, in linear or log scale.
    There are two modes.
      1. The number of elements in the ramp can be enforced (through N). If log
    is True the steps are regular in logarithmic scale, otherwise (default), they
    are linear.
      2. The step between elements can be enforced. If dlnx is set, then it is
    used as the step in ln(x), if it is dlogx, then the step is in log10(x), 
    otherwise, dx is the step in x. If the chosen step is not chosen exaclty to
    go from x0 to x1, there are two choices.
            a) if homogenize is True (default), the actual size is the step is
         slightly modified to have evenly spaced values between x0 and x1.
           b) if homogenize is False, the N-1 first steps are excatly the entered
         value and the last one is smaller.
    
    Copyright: F. Galliano
    '''
    if ((dx == None and dlnx == None and dlogx == None and N == None) \
        or (dx != None and dlnx == None and dlogx == None and N != None)):
        UT.strike('ramp','you should select either step or N')

    # Step mode
    elif (N == None):
        if (log == None):
            if (dx == None and dlnx == None):
                log = True
                log10 = True
            elif (dx == None and dlogx == None):
                log = True
                log10 = False
            elif (dlnx == None and dlogx == None):
                log = False
                log10 = False
        if (log):
            if (log10):
                N = np.round( (np.log10(x1)-np.log10(x0)) / dlogx ) + 1
                if (not homogenize):
                    x = 10**(np.arange(N,dtype=float)*dlogx)*x0
                    if (np.max(x) < x1): x = np.append(x,x1)
                    return(x)
            else:
                N = np.round( (np.log(x1)-np.log(x0)) / dlnx ) + 1
                if (not homogenize):
                    x = np.exp(np.arange(N,dtype=float)*dlnx)*x0
                    if (np.max(x) < x1): x = np.append(x,x1)
                    return(x)
        else:
            N = np.round( (x1-x0) / dx ) + 1
            if (not homogenize):
                x = np.arange(N,dtype=float)*dx + x0
                if (np.max(x) < x1): x = np.append(x,x1)
                return(x)
            
    # Number mode
    elif (dx == None and dlnx == None and dlogx == None):
        if (log == None): log = False

    # Generate the ramp
    dindgen = np.append( np.arange(N-1,dtype=float) / (N-1), 1)
    if (log):
        x = np.exp(dindgen*(np.log(x1)-np.log(x0)) + np.log(x0))
    else:
        x = dindgen*(x1-x0) + x0
    return(x)

def arrayize(arr,N=None,default=None,NP=True,dtype=None):
    '''
    Transform scalar or list into a numpy array

    Returns a numpy array or a list (if NP=False) whether the input is a scalar, 
    a list or an array. If N is set to an integer and arr is a scalar, then the 
    output is a size N array with the value arr at each element. If arr is None, 
    then default is substituted.

    Copyright: F. Galliano
    '''
    islist = isinstance(arr,list)
    isdict = isinstance(arr,dict)
    isnumpy = ( (not np.isscalar(arr)) and not islist and not isdict and \
                not isinstance(arr,type(None)) )
    arrout = arr
    if (not islist and not isnumpy):
        if (arrout == None): arrout = default
        if (N != None):
            if (np.size(arrout) == 1):
                arrout = [arrout for i in np.arange(N)]
            elif (np.size(arrout) != N):
                UT.strike('arrayize', 'wrong size for default.',
                          cat='InputError')
        else:
            arrout = [arrout]
        if (NP): arrout = np.array(arrout,dtype=dtype)
    elif (islist and NP):
        arrout = np.array(arrout,dtype=dtype)
    Nout = len(arrout)
    if (N != None):
        if (N != Nout):
            UT.strike('arrayize', 'array is not of size N.',
                          cat='InputError')
        return(arrout)
    else:
        return(arrout,Nout)


##------------------------------------------------
##
##              Unsorted functions
##
##------------------------------------------------

def listize(arr):
    '''
    Convert any iterable to list object
    worked for int, float, string, tuple, ndarray, list, dict, set, etc.
    '''
    if arr is None:
        listout = None
    elif np.isscalar(arr):
        listout = [arr] # scalar (string, int, float, etc.)
    elif isinstance(arr, np.ndarray):
        listout = arr.tolist() # ndarray
    elif isinstance(arr, np.recarray):
        listout = arr.tolist() # recarray
    elif isinstance(arr, Table):
        colnam = arr.colnames
        listout = [dict(zip(colnam, row)) for row in arr]
    else:
        listout = list(arr) # others

    return listout
    
def closest(arr, val, side=None):
    '''
    Return index of element in the array closest to a given value.
    The input array can be unsorted with NaN values. 
    However, if there are repeating elements, 
    the smallest index of the same closest value will be returned.
    If side is defined, while there are no qualified value in array,
    InputError will be raised (strict criterion).

    ------ INPUT ------
    arr                 input array
    val                 target value
    side                nearest left or right side of target value (Default: None)
    ------ OUTPUT ------
    ind                 index of the closet value
    '''
    if side=='left':
        arr2list = [x if x<=val else np.nan for x in arr]
    elif side=='right':
        arr2list = [x if x>=val else np.nan for x in arr]
    else:
        arr2list = list(arr)

    ## The first element in min func must not be np.nan
    if np.isnan(arr2list[0]):
        arr2list[0] = np.inf
    ## The min func uses key as iterable to calculate the min value,
    ## then use the position of this value in key to display value in arr.
    ## The index func reobtain (the index of) that position.
    ## In this case, the input arr can be unsorted.
    ind = arr2list.index(min(arr2list, key=lambda x:abs(x-val)))

    if np.isinf(arr2list[ind]):
        warnings.warn('side condition was ignored. ')
        
        arr2list = list(arr)
        if np.isnan(arr2list[0]):
            arr2list[0] = np.inf
        ind =  arr2list.index(min(arr2list, key=lambda x:abs(x-val)))
    
    return ind

def zoom(coord, scale=1, zp=0, grid=None,
         integer=False, start=0):
    '''
    Convert old coordinates to zoomed coord, given the scaling factor.
    
    When new coord are integers (integer=True), each value is rounded down to the closest integer.

    In DS9, the lower left pixel (1,1) is centered at the image grid (1,1)
    Whereas the image grid non-NaN zeropoint is (0.5,0.5)

    See also scipy.ndimage.zoom
    ------ INPUT ------
    coord               old coordinates
    scale               positive scaling factor (Default: 1)
    zp                  zeropoint shift
    grid                quick config of grid coordinates
                          'DS9' - DS9 pixel grid is integer starting with 1
                          None - read keywords "integer" and "start" (Default)
    integer             coordinates are integers (Default: False)
    start               coordinates start with 0/1 (Default: 0)
    ------ OUTPUT ------
    newcoord            new coordinates
                          - list of all possible values if integer=True
    '''
    ## Check inputs
    if scale<=0:
        UT.strike('zoom', 'scaling factor must be positive.',
                  cat='InputError')

    ## Known grids
    if grid=='DS9':
        integer = True
        start = 1

    if integer:
        if scale>1:
            Nsub = math.ceil(scale)
        else:
            Nsub = math.ceil(1./scale)
        # newcoord0 = math.floor((coord-start - zp)/scale + start)
        newcoord0 = math.floor((coord-start - zp)/scale)
        if newcoord0>=0:
            newcoord0 += start # No 0 in the coord
        newcoord = [newcoord0+i for i in range(Nsub)]
    else:
        # newcoord0 = (coord-start - zp)/scale + start
        newcoord0 = (coord-start - zp)/scale
        if newcoord0>=0:
            newcoord0 += start # No 0 in the coord
        newcoord = listize(newcoord0)

    return newcoord

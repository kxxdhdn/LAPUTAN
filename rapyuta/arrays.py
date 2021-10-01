#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Arrays

    arrayize, listize, closest, ramp,
    pix2sup, sup2pix, 

"""

import math
import numpy as np
import warnings

## Local
from utilities import InputError


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
                raise InputError('<arrayize>',
                                 'wrong size for default.')
        else:
            arrout = [arrout]
        if (NP): arrout = np.array(arrout,dtype=dtype)
    elif (islist and NP):
        arrout = np.array(arrout,dtype=dtype)
    Nout = len(arrout)
    if (N != None):
        if (N != Nout):
            raise InputError('<arrayize>',
                             'array is not of size N.')
        return(arrout)
    else:
        return(arrout,Nout)
    
def listize(arr):
    '''
    Convert any iterable to list object
    worked for int, float, string, tuple, ndarray, list, dict, set, etc.
    '''
    if np.isscalar(arr):
        listout = [arr] # scalar (string, int, float, etc.)
    elif isinstance(arr, np.ndarray):
        listout = arr.tolist() # ndarray
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


def pix2sup(x, xscale=1, x0=None, origin=0):
    '''
    Convert pixel to super pixel coordinates given super pixel zero point

    ------ INPUT ------
    x                   pixel coordinates
    xscale              super pixel size (Default: 1)
    x0                  zero point of super pixel
    origin              zero point convention
                          1 if in FITS and Fortran standards
                          0 if in Numpy and C standards
    ------ OUTPUT ------
    xs                  super pixel coordinates
    '''
    if origin==1:
        if x0 is None:
            x0 = 1
            
        if x-x0>=0:
            xs = math.floor((x-x0)/xscale) + 1
        else:
            xs = math.floor((x-x0)/xscale)
            
    elif origin==0:
        if x0 is None:
            x0 = 0
            
        xs = math.floor((x-x0)/xscale)
    
    return xs

def sup2pix(xs, xscale=1, x0=None, Npix=None, origin=0):
    '''
    Convert super pixel to pixel coordinates given super pixel zero point

    ------ INPUT ------
    xs                  super pixel coordinates
    xscale              super pixel size (Default: 1)
    x0                  zero point of super pixel
    origin              zero point convention
                          1 if in FITS and Fortran standards
                          0 if in Numpy and C standards
    ------ OUTPUT ------
    xarr                pixel coordinates
    '''
    if origin==1:
        if x0 is None:
            x0 = 1
            
        x = x0 + (xs-1) * xscale

    elif origin==0:
        if x0 is None:
            x0 = 0
            
        x = x0 + xs * xscale

    if Npix is None:
        xarr = [x+i for i in range(xscale)]
    else:
        xarr = []
        for i in range(xscale):
            if x+i<=Npix:
                xarr.append(x+i)
    
    return xarr

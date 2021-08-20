#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Arrays

    allist, closest

"""

import numpy as np
import warnings

## Local
from utilities import InputError


def allist(allIN):
    '''
    Convert any iterable to list object
    worked for int, float, string, tuple, ndarray, list, dict, set, etc.
    '''
    if np.isscalar(allIN):
        listOUT = [allIN] # scalar (string, int, float, etc.)
    elif isinstance(allIN, np.ndarray):
        listOUT = allIN.tolist() # ndarray
    else:
        listOUT = list(allIN) # others

    return listOUT
    
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

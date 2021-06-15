#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Arrays

    allist, closest

"""

import numpy as np

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
    
def closest(arr, val):
    '''
    Return the index i corresponding to the closest arr[i] to val
    '''
    arr = list(arr)
    
    return arr.index(min(arr, key=lambda x:abs(x-val)))


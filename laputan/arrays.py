#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Arrays

    allist, closest

"""

import numpy as np

def allist(allIN):
    '''
    All (str, int, float, ndarray, list) input to list output
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


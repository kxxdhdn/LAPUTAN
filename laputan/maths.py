#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Maths

    f_lin, f_lin0, f_lin1, gaussian, gaussian2D,
    rms, nanrms, std, nanstd, nanavg, bsplinterpol

"""

import math
import numpy as np
import scipy.interpolate as interpolate

def f_lin(x, A, B):
    '''
    Y = A * x + B
    '''
    return A * x + B

def f_lin0(x, A):
    '''
    Y = A * x
    '''
    return A * x

def f_lin1(x, B):
    '''
    Y = x + B
    '''
    return x + B

def gaussian(x, mu, sigma):
    '''
    Normalized Gaussian function given variable x
    '''
    return 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2))

def gaussian2D(x1, x2, mu1, mu2, sig1, sig2, A=1.):
    '''
    2D Gaussian function given iid variables x & y (and amplitude A)
    '''
    return A * np.exp(- (x1 - mu1)**2 / (2 * sig1**2) - (x2 - mu2)**2 / (2 * sig2**2))

def rms(a, ddof=0):
    '''
    Calculate root mean square
    
    ------ INPUT ------
    a                   an array or a list
    ddof                Delta Degrees of Freedom
    ------ OUTPUT ------
    rms                 root mean square of a
    '''
    n = np.size(a) - ddof
    a = np.array(a)
    ms = np.sum(a*a) / n

    return np.sqrt(ms)

def nanrms(a, ddof=0):
    '''
    Calculate root mean square (NaNs treated as zeros)
    
    ------ INPUT ------
    a                   an array or a list
    ddof                Delta Degrees of Freedom
    ------ OUTPUT ------
    rms                 root mean square of a
    '''
    n = np.size(a) - ddof
    a = np.array(a)
    ms = np.nansum(a*a) / n

    return np.sqrt(ms)

def std(a, ddof=0):
    '''
    The same as np.std
    '''
    n = np.size(a) - ddof
    a = np.array(a)
    mu = np.mean(a)
    ms = np.sum((a - mu)**2) / n
    
    return np.sqrt(ms)

def nanstd(a, axis=None, weights=None, MaskedValue=np.nan):
    '''
    Weighted standard deviation with NaNs ignored (NaN convention diff from nanrms *)
    MaskedValue presents if all NaNs
    '''
    ma = np.ma.MaskedArray(a, mask=np.isnan(a))
    if weights is not None:
        wgt = np.ma.MaskedArray(weights, mask=np.isnan(a))
    else:
        wgt = weights

    mask_any = ma.mask.any(axis=axis)
    mask_all = ma.mask.all(axis=axis) # All NaNs

    avg = np.average(ma, axis=axis, weights=wgt) # NaNs ignored
    ## Extend avg shape to that of a
    if axis is None:
        avg_ext = np.tile(avg, a.shape)
    elif axis==0:
        avg_ext = np.repeat(avg[np.newaxis,:,:], a.shape[0], axis=0)
        # avg_ext = np.tile(avg, (a.shape[0],1,1)) # alternative
    else:
        avg_ext = np.repeat(avg, ma.shape[axis], axis=axis-1).reshape(ma.shape)

    ## Deviation array of masked a from avg along axis
    dev_ma = ma - avg_ext

    if weights is not None:
        ## Count nonzero weights
        wgt_nz = np.ma.masked_where(wgt==0, wgt)
        Nwgt = wgt_nz.count(axis=axis)
        # print('Number of nonzero weights along axis {}: \n'.format(axis), Nwgt)
        std = np.sqrt(Nwgt/(Nwgt-1) * np.average(dev_ma**2, axis=axis, weights=wgt))
    else:
        std = np.sqrt(np.average(dev_ma**2, axis=axis, weights=wgt))

    ## Convert output none value convention (Default: NaNs)
    if axis is not None:
        std = std.data
        std[mask_all] = MaskedValue

    return std

def nanavg(a, axis=None, weights=None, MaskedValue=np.nan):
    '''
    Numpy.average with NaNs ignored (weight normalisation included)
    or
    Numpy.nanmean with weights
    '''
    ## NaNs -> --
    ma = np.ma.MaskedArray(a, mask=np.isnan(a))
    if weights is not None:
        wgt = np.ma.MaskedArray(weights, mask=np.isnan(a))
    else:
        wgt = weights

    mask_any = ma.mask.any(axis=axis)
    mask_all = ma.mask.all(axis=axis)

    avg = np.average(ma, axis=axis, weights=wgt)

    ## Check weight sums
    # print('wgt sums: ', np.average(ma, axis=axis, weights=wgt, returned=True)[1])

    ## Convert output none value convention (Default: NaNs)
    if axis is not None:
        avg = avg.data
        avg[mask_all] = MaskedValue

    return avg

def bsplinterpol(x, y, x0):
    '''
    Monte-Carlo propagated error calculator
    
    ------ INPUT ------
    x                   in base x
    y                   in data y
    x0                  out base x
    ------ OUTPUT ------
    bspl(x0)            B-spline interpol out data
    '''
    mask = []
    for i, yi in enumerate(y):
        if np.isnan(yi)==1 or yi==0:
            mask.append(i)
    if len(x)-len(mask)>4: # number of knots (avoid all NaNs col)
        x = np.delete(x, mask)
        y = np.delete(y, mask)

    t, c, k = interpolate.splrep(x, y, s=0, k=4) # s - smooth
    # print('''\
    # t: {}
    # c: {}
    # k: {}
    # '''.format(t, c, k))
    bspl = interpolate.BSpline(t, c, k, extrapolate=False)
    
    return bspl(x0)

"""
------------------------------ MAIN (test) ------------------------------
"""
if __name__ == "__main__":

    pass


#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Mathematics

    f_lin, f_lin0, f_lin1, gaussian, gaussian2D,
    rms, nanrms, std, nanstd, nanavg,
    pix2sr, sr2arcsec2, rad2arcsec, hour2deg, deg2hour,
    icorr2ij, ij2icorr,
    bsplinterp

"""

import math
import numpy as np
import scipy.interpolate as interpolate
import warnings

## Local
import utbox as UT


## Python string formatting
## https://pyformat.info


##------------------------------------------------
##
##                Basic functions
##
##------------------------------------------------

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


##------------------------------------------------
##
##               Basic statistics
##
##------------------------------------------------

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


##------------------------------------------------
##
##               Unit conversions
##
##------------------------------------------------

def pix2sr(X, CDELT):
    '''
    X pixel = Y sr
    ------ INPUT ------
    X                   float or ndarray in pixel
    CDELT               float or ndarray (same dim if X is also ndarray)
    '''
    PFOV = abs(CDELT)
    
    return X * (PFOV * 2. * math.pi / 360.)**2.

def sr2arcsec2(X):
    '''
    X sr = Y arcsec^2
    '''
    return X * (360. * 3600. / (2. * math.pi))**2.

def rad2arcsec(X):
    '''
    X rad = Y arcsec
    '''
    return X * 360. * 3600. / (2. * math.pi)

def hour2deg(h, m, s, deg, arcmin, arcsec):
    '''
    Hour angle to degree conversion
    '''
    ra = (h + m/60. + s/3600.) * 360./24.
    if deg<0:
        dec = -(-deg + arcmin/60. + arcsec/3600.)
    else:
        dec = deg + arcmin/60. + arcsec/3600.

    return ra, dec

def deg2hour(ra, dec):
    '''
    Degree to hour angle conversion
    '''
    h = math.floor(ra * 24./360.)
    m = math.floor((ra*24./360. - h) * 60.)
    s = ((ra*24./360. - h)*60 - m) * 60.
    if dec<0:
        deg = math.ceil(dec) # negtive
        arcmin = -math.ceil((dec - deg) * 60.) # positive
        arcsec = -((dec - deg)*60. + arcmin) * 60. # positive
    else:
        deg = math.floor(dec)
        arcmin = math.floor((dec - deg) * 60.)
        arcsec = ((dec - deg)*60. - arcmin) * 60.
        
    # print('{:d}h{:d}m{:04.2f}s, {:d}d{:d}m{:04.2f}s'.format(h,m,s,deg,arcmin,arcsec))

    return h, m, s, deg, arcmin, arcsec


##------------------------------------------------
##
##              Unsorted statistics
##
##------------------------------------------------

def icorr2ij(Npar, icorr=None, upper=True):
    '''
    Get the pair of parameter indices and their correlation coefficient index
    
    ------ INPUT ------
    Npar                Number of correlated parameters
    icorr               Index of correlation coefficients (Default: None)
                          None - return all ij pairs
    upper               Upper or lower triangles (Default: True)
    ------ OUTPUT ------
    ij                  Index pair
    '''
    Ncorr = math.comb(Npar,2) # Python >=3.8
    # Ncorr = int( Npar*(Npar-1)/2 ) # alternative
    ij = np.empty((Ncorr,2), dtype=int)
    ic = 0
    for i in range(Npar):
        for j in range(Npar):
            if j>i:
                if upper:
                    ij[ic,:] = [i+1,j+1]
                else:
                    ij[ic,:] = [j+1,i+1]
                ic += 1

    if icorr is None:
        warnings.warn('<maths.icorr2ij> icorr is None! All icorr2ij pairs are returned. ')
        return ij
    else:
        return ij[icorr-1]

def ij2icorr(Npar, i, j, verbose=False):
    '''
    Get the correlation coefficient index given the pair of parameter indices
    
    ------ INPUT ------
    Npar                Number of correlated parameters
    i                   Index of the first parameter (>= 1)
    j                   Index of the second parameter (>= 1)
    verbose             Default: False
    ------ OUTPUT ------
    icorr               Indices of correlation coefficient
    '''
    icorr = np.empty((Npar,Npar), dtype=int)
    ic = 0
    for ii in range(Npar):
        for jj in range(Npar):
            if jj>ii:
                ic += 1
                icorr[jj,ii] = ic
                icorr[ii,jj] = ic
            elif ii==jj:
                icorr[ii,jj] = 0
                ## Return all ij2icorr pairs
    if verbose:
        if i<j:
            print('Upper triangle detected.')
        elif i>j:
            print('Lower triangle detected.')
        else:
            warnings.warn('<maths.ij2icorr> i equals to j! All ij2icorr pairs are returned. ')
    if i==j:
        return icorr
    else:
        return icorr[j-1,i-1]

    ## Opt.1
    # try:
    #     icorr = 0
    #     for ii in range(Npar):
    #         for jj in range(Npar):
    #             if j>i:
    #                 if jj>ii:
    #                     icorr += 1
    #                 if ii==i-1 and jj==j-1:
    #                     raise StopIteration
    #             elif i>j:
    #                 if ii>jj:
    #                     icorr += 1
    #                 if ii==i-1 and jj==j-1:
    #                     raise StopIteration
    #             else:
    #                 raise StopIteration
    # except StopIteration:
    #     return icorr

    ## Opt.2
    # return (i-1)*Npar - math.comb(i+1,2) + j


##------------------------------------------------
##
##              Numerical analysis
##
##------------------------------------------------

def bsplinterp(x, y, x0):
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

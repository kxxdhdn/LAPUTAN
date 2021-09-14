#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging
# logging.disable(sys.maxsize)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/dat/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
	os.makedirs(outdir)

import numpy as np
import matplotlib.pyplot as plt

## Local
sys.path.insert(0, testdir+'/..') ## rapyuta path
from rapyuta.maths import (f_lin, f_lin1, f_lin0,
	                          gaussian, gaussian2D,
                              rms, nanrms, std, nanstd,
                              nanavg, bsplinterpol, )

x0 = np.arange(-5., 5., .5)
rand = np.random.uniform(-1., 1., 20)
x = np.arange(-5., 5., .1)
x[10:20] = np.nan
a = np.arange(24, dtype=float).reshape((4,3,2))
a[(a>5) & (a<12)] = np.nan
wgt = np.arange(0, 1.2, .05).reshape((4,3,2))

print('\n TEST gaussian ')
print('---------------')
g = gaussian(x0, 0., 1.)
# print(g)
fig = plt.figure('TEST gaussian')
ax = fig.add_subplot()
ax.plot(x0, g, c='k')

print('\n TEST gaussian2D ')
print('-----------------')
g2 =gaussian2D(x0, x0*2, 0., 0., 1., 1.)
fig2 = plt.figure('TEST gaussian2d')
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(x0, x0*2, g2, s=.5)

print('\n TEST rms ')
print('----------')
print('rms (ddof=0) = ', rms(x0))
print('rms (ddof=1) = ', rms(x0, ddof=1))

print('\n TEST nanrms ')
print('-------------')
print('rms (ddof=0) = ', nanrms(x))
print('rms (ddof=1) = ', nanrms(x, ddof=1))

print('\n TEST std ')
print('----------')
print('numpy.std = ', np.std(x0, ddof=0))
print('std = ', std(x0))

print('\n TEST nanstd ')
print('-------------')
print('numpy.nanstd(a, axis=1) = \n', np.nanstd(a, 1))
print('nanstd(a, axis=1) = \n', nanstd(a, 1))
print('nanstd(a, axis=1, weights=wgt) = \n', nanstd(a, 1, wgt))
print('nanstd(a, axis=1, MaskedValue=100) = \n', nanstd(a, 1, MaskedValue=100), '\n')

# print('nanstd(a, axis=(0,2)) = \n', nanstd(a, (0,2)), '\n')

print('nanstd(a, axis=None) = ', nanstd(a))

print('\n TEST nanavg ')
print('-------------')
print('numpy.average(a, axis=1) = \n', np.average(a, 1))
print('nanavg(a, axis=1) = \n', nanavg(a, 1))
print('nanavg(a, axis=1, weights=wgt) = \n', nanavg(a, 1, wgt))
print('nanavg(a, axis=1 MaskedValue=100) = \n', nanavg(a, 1,  MaskedValue=100), '\n')

print('numpy.average(a, axis=(0,2)) = \n', np.average(a, (0,2)))
print('nanavg(a, axis=(0,2)) = \n', nanavg(a, (0,2)), '\n')

print('numpy.average(a, axis=None) = ', np.average(a))
print('nanavg(a, axis=None) = ', nanavg(a))

print('\n TEST bsplinterpol ')
print('-------------------')
x_interp = np.arange(-5., 5., .01)
g_interp = bsplinterpol(x0, g+.005*rand, x_interp)
ax.plot(x0, g+.005*rand, c='c')
ax.plot(x_interp, g_interp, c='r')
print('See figure TEST gaussian [Done]')

plt.show()

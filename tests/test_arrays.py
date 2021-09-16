#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging
# logging.disable(sys.maxsize)

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/lib/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
	os.makedirs(outdir)

import numpy as np

## Local
sys.path.insert(0, testdir+'/..') ## rapyuta path
from rapyuta.arrays import (listize, closest, )

a = np.arange(24, dtype=float).reshape(4,3,2)
b = [['a', 'b'], ['c', 'd']]
c = 'abcd'
d = 1.
x = np.arange(-5., 5., .1)
## x with NaNs
x[10:20] = np.nan
x[90:] = np.nan
## Unsorted x with repeating elements
x[60:70] = -2.345

print('\n TEST listize ')
print('-------------')
print('ndarray to list: ', listize(a), '\n')
print('string to list: ', listize(c), '\n')
print('int/float to list: ', listize(d), '\n')
print('list: ', listize(b))

print('\n TEST closest ')
print('--------------')
ind = closest(x, -2.35)
print('x[{0}]={1} is closest to -2.35'.format(ind, x[ind]))
ind = closest(x, -2.35, side='left')
print('x[{0}]={1} is the left nearest of -2.35'.format(ind, x[ind]))
ind = closest(x, -2.35, side='right')
print('x[{0}]={1} is the right nearest of -2.35'.format(ind, x[ind]))
ind = closest(x, -6, side='right')
print('x[{0}]={1} is the right nearest of -6'.format(ind, x[ind]))
ind = closest(x, 6, side='right')
print('x[{0}]={1} is the right nearest of 6'.format(ind, x[ind]))

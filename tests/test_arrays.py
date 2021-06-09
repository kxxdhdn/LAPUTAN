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
sys.path.insert(0, testdir+'/..') ## laputan path
from laputan.arrays import (allist, closest, )

a = np.arange(24, dtype=float).reshape(4,3,2)
b = [['a', 'b'], ['c', 'd']]
c = 'abcd'
d = 1.
x = np.arange(-5., 5., .1)
x[10:20] = np.nan


print('\n TEST allist ')
print('-------------')
print('ndarray to list: ', allist(a), '\n')
print('string to list: ', allist(c), '\n')
print('int/float to list: ', allist(d), '\n')
print('list: ', allist(b))

print('\n TEST closest ')
print('--------------')
ind = closest(x, -2.35)
print('x[{0}]={1} is closest to -2.35'.format(ind, x[ind]))

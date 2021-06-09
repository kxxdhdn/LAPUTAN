#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Plots

    plotool:
        figure, set_border, Cartesian, set_ax,
        plot, set_font, save, show
    pplot(plotool):
        add_plot

"""

from astropy import units as u
import numpy as np
from scipy import optimize
# import matplotlib as mpl
import matplotlib.colors as mplc
import matplotlib.pyplot as plt

# cmap = mpl.cm.viridis
# norm = mpl.colors.Normalize(vmin=0, vmax=1)
sizeXS, sizeS, sizeM = 4, 6, 8
sizeL, sizeXL, sizeXXL = 10, 12, 14

class plotool:
    '''
    PLOT tOOL
    '''
    def __init__(self, x=np.zeros(2), y=np.zeros(2),
                 xerr=None, yerr=None):
        
        # INPUTS
        self.x = x
        self.y = y
        self.xerr = xerr
        self.yerr = yerr

        self.ax = None

    def figure(self, figsize=None, figint=True,
               nrows=1, ncols=1):
        
        if figint==True:
            plt.ion()

        self.nrows = nrows
        self.ncols = ncols

        if nrows==1 and ncols==1:
            self.fig, self.ax = plt.subplots(nrows, ncols,
                figsize=figsize)
        else:
            self.fig, self.axes = plt.subplots(nrows, ncols,
                figsize=figsize)
            self.ax = self.axes[0,0]
        
    def set_border(self, left=None, bottom=None, right=None, top=None,
                   wspace=None, hspace=None):

        plt.subplots_adjust(left=left, bottom=bottom,
            right=right, top=top, wspace=wspace, hspace=hspace)

    def Cartesian(self, c=None, ls=None, lw=None,
                  ec=None, elw=None, lab=None):

        self.markers, self.caps, self.bars = self.ax.errorbar(
            self.x, self.y, yerr=self.yerr, xerr=self.xerr, c=c,
            ls=ls, lw=lw, ecolor=ec, elinewidth=elw, label=lab)

    def set_ax(self, xlog=False, ylog=False,
               basex=10,basey=10,nonposx='sym', nonposy='sym',
               xlim=(None,None), ylim=(None,None),
               xlab=None, ylab=None, legend=None, title=None):
        '''
        nonposx, nonposy: 'sym', 'mask', 'clip'
        '''

        if xlog==True:
            if nonposx=='sym':
                self.ax.set_xscale('symlog',base=basex)
            else:
                self.ax.set_xscale('log',base=basex,nonpositive=nonposx)
        if ylog==True:
            if nonposx=='sym':
                self.ax.set_yscale('symlog',base=basey)
            else:
                self.ax.set_yscale('log',base=basey,nonpositive=nonposy)
                
        self.ax.set_xlim(xlim[0], xlim[1])
        self.ax.set_ylim(ylim[0], ylim[1])

        # self.ax.set_xticks()
        # self.ax.set_yticks()
        # self.ax.set_xticklabels()
        # self.ax.set_yticklabels()
        
        self.ax.set_xlabel(xlab)
        self.ax.set_ylabel(ylab)
        self.ax.set_title(title)
        self.ax.legend(loc=legend)
        self.legend = legend
    
    def plot(self, nrow=1, ncol=1, x=None, y=None,
             xlim=(None, None), ylim=(None, None),
             xlog=False, ylog=False,
             basex=10,basey=10,nonposx='sym', nonposy='sym',
             c=None, ls=None, lw=None,
             ec=None, elw=None, xlab=None, ylab=None,
             lab=None, legend=None, title=None, mod='CA'):

        if x is None:
            x = self.x
        if y is None:
            y = self.y
            
        if mod=='CA':
            if self.nrows==1 and self.ncols==1:
                self.Cartesian(c, ls, lw, ec, elw, lab)
            else:
                self.ax = self.axes[nrow-1,ncol-1]
                self.Cartesian(c, ls, lw, ec, elw, lab)
            
            self.set_ax(xlog, ylog, basex, basey, nonposx, nonposy,
                        xlim, ylim, xlab, ylab, legend, title)
        
        else:
            print('*******************')
            print('Prochainement...')
            
            print('PL: polar')
            print('CL: cylindrical')
            print('SP: spherical')
            print('*******************')

    def set_font(self, fontsize=sizeM, subtitlesize=sizeM,
                 axesize=sizeS, xticksize=sizeS, yticksize=sizeS,
                 legendsize=sizeM, figtitlesize=sizeL):

        plt.rc('font', size=fontsize)            # general text
        plt.rc('axes', titlesize=subtitlesize)   # axes title
        plt.rc('axes', labelsize=axesize)        # x and y labels
        plt.rc('xtick', labelsize=xticksize)     # x tick
        plt.rc('ytick', labelsize=yticksize)     # y tick
        plt.rc('legend', fontsize=legendsize)    # legend
        plt.rc('figure', titlesize=figtitlesize) # figure title


    def save(self, savename=None, transparent=False):

        if savename is not None:
            self.fig.savefig(savename, transparent=transparent)
        else:
            print('WARNING: not saved! ')

    def show(self):

        plt.ioff()
        plt.show()

##-----------------------------------------------

##                    sub-classes

##-----------------------------------------------

class pplot(plotool):
    '''
    Plot single 2D curve
    '''
    def __init__(self, x=None, y=None, xerr=None, yerr=None,
                 xlim=(None, None), ylim=(None,None),
                 xlog=None, ylog=None,
                 basex=10,basey=10,nonposx='sym', nonposy='sym',
                 c=None, ls=None, lw=.5, ec='r', elw=.8,
                 xlab='X', ylab='Y', lab=None, legend=None, title='2D Curve',
                 figsize=None, figint=False,
                 left=.1, bottom=.1, right=.99, top=.9,
                 wspace=None, hspace=None, clib='base'):
        super().__init__(x, y, xerr, yerr)

        if x is None:
            x = self.x
        if y is None:
            y = self.y

        ## Named color lib
        if clib=='base':
            self.clib = list(mplc.BASE_COLORS) # 8 colors
        elif clib=='tableau':
            self.clib = list(mplc.TABLEAU_COLORS) # 10 colors
        elif clib=='ccs4' or clib=='x11':
            self.clib = list(mplc.CSS4_COLORS)
        elif clib=='xkcd':
            self.clib = list(mplc.XKCD_COLORS)
        else:
            raise ValueError('Unname colors! ')

        self.figure(figsize, figint)
        self.iplot = 0

        self.set_border(left=left, bottom=bottom,
            right=right, top=top, wspace=wspace, hspace=hspace)

        if c is None:
            c = self.clib[self.iplot]
        self.Cartesian(c, ls, lw, ec, elw, lab)

        self.set_ax(xlog, ylog, basex, basey, nonposx, nonposy,
                    xlim, ylim, xlab, ylab, legend, title)

        self.set_font()

    def add_plot(self, x=None, y=None, xerr=None, yerr=None,
                 c=None, ls=None, lw=.5, ec='r', elw=.8, lab=None):

        if x is None:
            x = self.x
        else:
            self.x = x
        if y is None:
            y = self.y
        else:
            self.y = y
            
        self.iplot += 1
        if c is None:
            c = self.clib[self.iplot]
            
        self.Cartesian(c, ls, lw, ec, elw, lab)

        self.set_ax(legend=self.legend)

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Plots

    plotool:
        figure, set_border, set_clib, set_ax,
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
import warnings

## Local
from utilities import merge_aliases

# cmap = mpl.cm.viridis
# norm = mpl.colors.Normalize(vmin=0, vmax=1)
sizeXS, sizeS, sizeM = 4, 6, 8
sizeL, sizeXL, sizeXXL = 10, 12, 14

##-----------------------------------------------
##
##            <plotool> based tools
##
##-----------------------------------------------

class plotool:
    '''
    plot Tool
    '''
    def __init__(self, x=np.zeros(2), y=np.zeros(2)):
        
        # INPUTS
        self.x = x
        self.y = y

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

    def set_clib(self, clib):
        if clib=='base':
            self.clib = list(mplc.BASE_COLORS) # 8 colors
        elif clib=='tableau':
            self.clib = list(mplc.TABLEAU_COLORS) # 10 colors
        elif clib=='ccs4' or clib=='x11':
            self.clib = list(mplc.CSS4_COLORS)
        elif clib=='xkcd':
            self.clib = list(mplc.XKCD_COLORS)
        else:
            self.clib = clib
        
    def set_ax(self, xlog=False, ylog=False,
               basex=10, basey=10, nonposx='sym', nonposy='sym',
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

        if xlim[0]!=None or xlim[1]!=None:
            self.ax.set_xlim(xlim[0], xlim[1])
        if ylim[0]!=None or ylim[1]!=None:
            self.ax.set_ylim(ylim[0], ylim[1])

        # self.ax.set_xticks()
        # self.ax.set_yticks()
        # self.ax.set_xticklabels()
        # self.ax.set_yticklabels()

        if xlab is not None:
            self.ax.set_xlabel(xlab)
        if ylab is not None:
            self.ax.set_ylabel(ylab)
        if title is not None:
            self.ax.set_title(title)
        if legend is not None:
            self.ax.legend(loc=legend)
        self.legend = legend
    
    def plot(self, nrow=1, ncol=1,
             x=None, y=None, xerr=None, yerr=None,
             fmt='', capsize=None, barsabove=False, # errorbar kw
             ecolor=None, ec=None, elinewidth=None, elw=None, # errorbar kw
             mod='CA', **kwargs):
        '''
        Like set_ax(), this is a clump operation.
        The idea is to all set in one command,
        while each single operation should also be valid.
        '''

        ## kw aliases
        ec = merge_aliases(None, ecolor=ecolor, ec=ec)
        elw = merge_aliases(None, elinewidth=elinewidth, elw=elw)
        
        if x is None:
            x = self.x
        else:
            self.x = x
        if y is None:
            y = self.y
        else:
            self.y = y

        ## CA: Cartesian using matplotlib.pyplot.errorbar
        if mod=='CA':
            if self.nrows!=1 or self.ncols!=1:
                self.ax = self.axes[nrow-1,ncol-1]
                
            self.markers, self.caps, self.bars = self.ax.errorbar(
                x=x, y=y, yerr=yerr, xerr=xerr,
                fmt=fmt, ecolor=ec, elinewidth=elw,
                capsize=capsize, barsabove=barsabove,
                **kwargs)
        
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
            warnings.warn('Not saved! ')

    def show(self):

        plt.ioff()
        plt.show()

class pplot(plotool):
    '''
    Uni-frame plot (1 row * 1 col)
    '''
    def __init__(self, x=None, y=None, xerr=None, yerr=None,
                 fmt='', capsize=None, barsabove=False, # errorbar kw
                 ecolor=None, ec=None, elinewidth=None, elw=None, # errorbar kw
                 figsize=None, figint=False, # figure kw
                 left=.1, bottom=.1, right=.99, top=.9, # set_border kw
                 wspace=None, hspace=None, # set_border kw
                 xlim=(None, None), ylim=(None,None), xlog=None, ylog=None, # set_ax kw
                 basex=10, basey=10, nonposx='sym', nonposy='sym', # set_ax kw
                 xlab='X', ylab='Y', legend=None, title='Untitled', # set_ax kw
                 clib='base', c=None, **kwargs):
        super().__init__(x, y)

        self.iplot = 0

        ## kw aliases
        ec = merge_aliases(None, ecolor=ecolor, ec=ec)
        elw = merge_aliases(None, elinewidth=elinewidth, elw=elw)

        ## Auto color
        self.set_clib(clib)
        if c is None:
            c = self.clib[self.iplot]

        ## Init figure
        self.figure(figsize, figint)

        ## set_border
        self.set_border(left=left, bottom=bottom,
            right=right, top=top, wspace=wspace, hspace=hspace)

        ## plot
        self.plot(x=x, y=y, xerr=xerr, yerr=yerr,
                  fmt=fmt, ec=ec, elw=elw, # errorbar kw
                  capsize=capsize, barsabove=barsabove, # errorbar kw
                  c=c, **kwargs)

        ## set_ax
        self.set_ax(xlog, ylog, basex, basey, nonposx, nonposy,
                    xlim, ylim, xlab, ylab, legend, title)
        
        self.set_font()

    def add_plot(self, x=None, y=None, xerr=None, yerr=None,
                 fmt='', capsize=None, barsabove=False, # errorbar opt
                 ecolor=None, ec=None, elinewidth=None, elw=None, # errorbar opt
                 c=None, **kwargs):

        self.iplot += 1

        ## kw aliases
        ec = merge_aliases(None, ecolor=ecolor, ec=ec)
        elw = merge_aliases(None, elinewidth=elinewidth, elw=elw)

        ## Auto color
        if self.iplot==len(self.clib):
            self.iplot = 0
        if c is None:
            c = self.clib[self.iplot]
        
        if x is None:
            x = self.x
        else:
            self.x = x
        if y is None:
            y = self.y
        else:
            self.y = y

        ## plot
        self.plot(x=x, y=y, xerr=xerr, yerr=yerr,
                  fmt=fmt, ec=ec, elw=elw, # errorbar kw
                  capsize=capsize, barsabove=barsabove, # errorbar kw
                  c=c, **kwargs)

        ## Add legend
        self.ax.legend(loc=self.legend)

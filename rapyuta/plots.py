#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Plots

    ellipse
    plotool:
        figure, set_figure, set_clib, set_ax, set_legend, 
        plot, save, show
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
##                 Plot utilities
##
##-----------------------------------------------

# def ellipse(xmean=None, ymean=None, xstdev=None, ystdev=None, rho=None,
#             xlim=(None,None), ylim=(None,None), Npt=300, xlog=False, ylog=False):
#     '''
#     Uncertainty ellipse

#     Function used to plot uncertainty ellipses
#     (or 1-sigma contour of a bivariate normal distribution).
#     Inspired by F. Galliano's Python library.
#     '''
#     x = AR.ramp(x0=xmean-xstdev*(1-1.E-5),x1=xmean+xstdev*(1-1.E-5),N=Npt)
    
#     c1 = rho * (x-xmean)/xstdev
#     c2 = NP.sqrt( (1-rho**2) * (1-(x-xmean)**2/xstdev**2) )
#     y1 = ystdev * ( c1 - c2 ) + ymean
#     y2 = ystdev * ( c1 + c2 ) + ymean
#     xplot = NP.concatenate((x,x[::-1],[x[0]]))
#     yplot = NP.concatenate((y1,y2[::-1],[y1[0]]))
#     if (xlog): xplot = NP.exp(xplot)
#     if (ylog): yplot = NP.exp(yplot)
#     if (xmin != None): xplot[xplot < xmin] = xmin
#     if (xmax != None): xplot[xplot > xmax] = xmax
#     if (ymin != None): yplot[yplot < ymin] = ymin
#     if (ymax != None): yplot[yplot > ymax] = ymax
#     return(xplot,yplot)

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
        
        if figint:
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

    def set_fig(self, left=None, right=None,
                bottom=None, top=None,
                wspace=None, hspace=None,
                title=None, tfsize=20):

        self.fig.subplots_adjust(left=left, right=right,
            bottom=bottom, top=top, wspace=wspace, hspace=hspace)
        
        if title is not None:
            self.fig.suptitle(title,size=tfsize)

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
        
    def set_ax(self, subax=(1,1), # ax = axes[subax[0]-1,subax[1]-1]
               xlog=False, ylog=False, # ax.set_xscale
               basex=10, basey=10, nonposx='clip', nonposy='clip', # ax.set_xscale
               xlim=(None,None), ylim=(None,None), #ax.set_xlim
               xtickfsize=10, ytickfsize=10, # ax.xaxis.set_tick_params(labelsize=)
               xlabel=None, ylabel=None, xfsize=10, yfsize=10, # ax.set_xlabel
               title=None, tfsize=10, # ax.set_title (subplot title)
               ):
        '''
        nonposx, nonposy: 'sym', 'mask', 'clip'
        '''
        if self.nrows!=1 or self.ncols!=1:
            self.ax = self.axes[subax[0]-1,subax[1]-1]

        if xlog:
            if nonposx=='sym':
                self.ax.set_xscale('symlog',base=basex)
            else:
                self.ax.set_xscale('log',base=basex,nonpositive=nonposx)
        if ylog:
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

        self.ax.xaxis.set_tick_params(labelsize=xtickfsize)
        self.ax.yaxis.set_tick_params(labelsize=ytickfsize)
        
        if xlabel is not None:
            self.ax.set_xlabel(xlabel,size=xfsize)
        if ylabel is not None:
            self.ax.set_ylabel(ylabel,size=yfsize)
        if title is not None:
            self.ax.set_title(title,size=tfsize)
        
    def set_legend(self, subax=None,
                   loc='center right', fontsize=10,
                   anchor=None, shrinkx=1., shrinky=1.,
                   **kwargs):
        '''
        - bbox_to_anchor rules: (1,1) correspond to upper right of the axis
                 bbox_to_anchor = (1,1)
        .--------.
        |        |
        |  axis  |
        |        |
        .--------.
        
        - lengend loc is relative to the bbox_to_anchor as follows:
                           |
        lower right        |       lower left
        ------------bbox_to_anchor-----------
        upper right        |       upper left
                           |
        '''
        if subax is None:
            self.fig.subplots_adjust(right=shrinkx,
                                     top=shrinky)
            
            self.fig.legend(loc=loc, fontsize=fontsize, bbox_to_anchor=anchor,
                            **kwargs)
        else:
            if self.nrows!=1 or self.ncols!=1:
                self.ax = self.axes[subax[2]-1,subax[1]-1]
            
            # shrink current axis
            box = self.ax.get_position()
            self.ax.set_position([box.x0, box.y0, box.width*shrinkx, box.height*shrinky])

            self.ax.legend(loc=loc, fontsize=fontsize, bbox_to_anchor=anchor,
                           **kwargs)            
        
    def plot(self, subax=(1,1),
             x=None, y=None, xerr=None, yerr=None,
             fmt='', capsize=None, barsabove=False, # errorbar kw
             ecolor=None, ec=None, elinewidth=None, elw=None, # errorbar kw
             mod='CA', **kwargs):
        '''
        Like set_ax(), this is a clump operation.
        The idea is to all set in one command,
        while each single operation should also be valid.
        '''
        if self.nrows!=1 or self.ncols!=1:
            self.ax = self.axes[subax[0]-1,subax[1]-1]

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
    def __init__(self, x=None, y=None, xerr=None, yerr=None, # errorbar kw
                 fmt='', capsize=None, barsabove=False, # errorbar kw
                 ecolor=None, ec=None, elinewidth=None, elw=None, # errorbar kw
                 figsize=None, figint=False, # figure kw
                 left=.1, bottom=.1, right=.99, top=.9, # set_fig kw
                 wspace=.1, hspace=.1, # set_fig kw
                 title='Untitled', titlesize=None, # set_fig kw
                 xlog=None, ylog=None, # set_ax kw
                 basex=10, basey=10, nonposx='clip', nonposy='clip', # set_ax kw
                 xlim=(None, None), ylim=(None,None), # set_ax kw
                 ticksize=None, # set_ax kw
                 xlabel='X', ylabel='Y', labelsize=None, # set_ax kw
                 legend=None, legendsize=None, anchor=None, # set_legend kw
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

        ## set_fig
        self.set_fig(left=left, bottom=bottom, right=right, top=top,
            wspace=wspace, hspace=hspace, title=title, tfsize=titlesize)

        ## plot
        self.plot(x=x, y=y, xerr=xerr, yerr=yerr,
                  fmt=fmt, ec=ec, elw=elw, # errorbar kw
                  capsize=capsize, barsabove=barsabove, # errorbar kw
                  c=c, **kwargs)

        ## set_ax
        self.set_ax((1,1), xlog, ylog, basex, basey, nonposx, nonposy, xlim, ylim,
                    ticksize, ticksize, xlabel, ylabel, labelsize, labelsize)

        if legend is not None:
            self.set_legend(loc=legend, fontsize=legendsize, anchor=anchor)
        self.legend = legend
        self.legendsize = legendsize
        self.anchor = anchor

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
        if self.legend is not None:
            self.set_legend(loc=self.legend,
                fontsize=self.legendsize, anchor=self.anchor)

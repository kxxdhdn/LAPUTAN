#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Visualization

    pplot(plotutils):
        add_plot

"""

import warnings
import numpy as np
# import matplotlib as mpl
import matplotlib.patches as mpatches

## Local
import rapyuta.utbox as UT
import rapyuta.latte as LA
from .utils import *


##------------------------------------------------
##
##            <plotutils> based tools
##
##------------------------------------------------

class pplot(plotutils):
    '''
    Uni-frame plot (1 row * 1 col)

    ------ INPUT ------
    plotutils.plot(**kwargs)
    plotutils.eplot(**kwargs)
    '''
    def __init__(self, x=None, y=None, yerr=None, xerr=None,
                 figsize=(8,6), figint=False,
                 ## errorbar kw
                 fmt='', capsize=None, barsabove=False,
                 ## eplot kw
                 edgecolor='grey', ecolor='grey', ec='grey',
                 elinewidth=None, elw=None, elinestyle=None, els=None,
                 mask=None, xmin=None, xmax=None, ymin=None, ymax=None,
                 sigmax=None, sigmay=None, rho=None,
                 gammax=None, gammay=None,xisln=False, yisln=False,
                 efillcolor=None, efc=None, efill=None, ehatch=None,
                 errinlegend=None, alpha=1,
                 ## set_fig kw
                 left=.15, bottom=.1, right=.95, top=.9,
                 wspace=None, hspace=None, title=None, titlesize=20,
                 ## set_ax kw
                 xlog=None, ylog=None,
                 basex=10, basey=10, nonposx='clip', nonposy='clip',
                 xlim=(None, None), ylim=(None,None),
                 tk=None, tkmi=None, tkform=None, tksize=None, # same xy
                 xtk=None, xtkmi=None, xtkform=None, xtksize=15, # x
                 ytk=None, ytkmi=None, ytkform=None, ytksize=15, # y
                 xlabel='X', ylabel='Y', xsize=None, ysize=None, xysize=15,
                 ## set_legend kw
                 loc=None, legendsize=15, legendalpha=1,
                 anchor=None, figtight=False,
                 ## Other kw
                 autocol='base', c=None, **kwargs):
        super().__init__(figsize=figsize, figint=figint, x=x, y=y)

        self.pltid = 0

        ## Keyword aliases
        ## Note that ec and ecolor are shared by plot (errorbar) and eplot (ellipse)
        ec = UT.merge_aliases('grey', edgecolor=edgecolor, ecolor=ecolor, ec=ec) # Default: 'grey'
        elw = UT.merge_aliases(None, elinewidth=elinewidth, elw=elw)
        els = UT.merge_aliases(None, elinestyle=elinestyle, els=els)
        efc = UT.merge_aliases(None, efillcolor=efillcolor, efc=efc)

        ## Auto color
        self.set_autocol(autocol)
        if c=='auto':
            c = self.autocol[self.pltid]

        ## set_fig
        self.set_fig(left=left, bottom=bottom, right=right, top=top,
            wspace=wspace, hspace=hspace, title=title, tsize=titlesize)

        ## set_ax
        if tk is not None:
            xtk = tk
            ytk = tk
        if tkmi is not None:
            xtkmi = tkmi
            ytkmi = tkmi
        if tkform is not None:
            xtkform = tkform
            ytkform = tkform
        if tksize is not None:
            xtksize = tksize
            ytksize = tksize
        if xysize is not None:
            xsize = xysize
            ysize = xysize
        self.set_ax(xlog=xlog, ylog=ylog,
                    basex=basex, basey=basey, nonposx=nonposx, nonposy=nonposy,
                    xlim=xlim, ylim=ylim,
                    xtk=xtk, xtkmi=xtkmi, xtkform=xtkform, xtksize=xtksize,
                    ytk=ytk, ytkmi=ytkmi, ytkform=ytkform, ytksize=ytksize,
                    xlabel=xlabel, xsize=xsize, ylabel=ylabel, ysize=ysize)

        ## plot
        ell = (rho is not None)
        if (not ell):
            
            self.plot(x=x, y=y, yerr=yerr, xerr=xerr,
                      xisln=xisln, yisln=yisln,
                      fmt=fmt, ec=ec, elw=elw, els=els, # errorbar kw
                      capsize=capsize, barsabove=barsabove, # errorbar kw
                      c=c, alpha=alpha, **kwargs)
            self.append_handles()
            
        else:

            self.plot(x=x, y=y, xisln=xisln, yisln=yisln,
                      fmt=fmt, c=c, alpha=alpha, **kwargs)
            self.append_handles()
            ## Ellipse errors
            self.eplot(x=x, y=y, mask=mask,
                       xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                       xisln=xisln, yisln=yisln,
                       sigmax=sigmax, sigmay=sigmay, rho=rho,
                       gammax=gammax, gammay=gammay,
                       ec=ec, elw=elw, els=els,
                       efill=efill, efc=efc, ehatch=ehatch,
                       errinlegend=errinlegend, alpha=alpha)
            self.append_handles()
            ## Replace line by patch in legend
            epatch = mpatches.Patch(ec=ec, lw=elw, ls=elw,
                                    fill=efill, fc=efc, alpha=alpha,
                                    hatch=ehatch, label=errinlegend)
            self.handles[-1] = epatch

        ## set_legend
        if loc is not None:
            self.set_legend(loc=loc, fontsize=legendsize,
                            bbox_to_anchor=anchor, figtight=figtight,
                            framealpha=legendalpha)
        self.loc = loc
        self.legendsize = legendsize
        self.anchor = anchor
        self.figtight = figtight
        self.legendalpha = legendalpha

    def add_plot(self, x=None, y=None, yerr=None, xerr=None,
                 ## errorbar kw
                 fmt='', capsize=None, barsabove=False,
                 ## eplot kw
                 edgecolor='grey', ecolor='grey', ec='grey',
                 elinewidth=None, elw=None, elinestyle=None, els=None,
                 mask=None, xmin=None, xmax=None, ymin=None, ymax=None,
                 sigmax=None, sigmay=None, rho=None,
                 gammax=None, gammay=None,xisln=False, yisln=False,
                 efillcolor=None, efc=None, efill=None, ehatch=None,
                 errinlegend=None, alpha=1,
                 ## set_legend kw (only when addlegend=True)
                 addlegend=False, loc=None, legendsize=15,
                 legendalpha=1, anchor=None, figtight=False,
                 ## Other (errorbar) kw
                 c=None, label=None, **kwargs):
        '''
        ------ INPUT ------
        plotutils.plot(**kwargs)
        plotutils.eplot(**kwargs)
        '''
        self.pltid += 1

        ## Keyword aliases
        ## Note that ec and ecolor are shared by plot (errorbar) and eplot (ellipse)
        ec = UT.merge_aliases('grey', edgecolor=edgecolor, ecolor=ecolor, ec=ec) # Default: 'grey'
        elw = UT.merge_aliases(None, elinewidth=elinewidth, elw=elw)
        els = UT.merge_aliases(None, elinestyle=elinestyle, els=els)
        efc = UT.merge_aliases(None, efillcolor=efillcolor, efc=efc)

        ## Auto color
        if self.pltid==len(self.autocol):
            self.pltid = 0
        if c=='auto':
            c = self.autocol[self.pltid]
        
        if x is None:
            x = self.x
        else:
            self.x = x
        if y is None:
            y = self.y
        else:
            self.y = y

        ## plot
        if addlegend:
            ## The reset kw of append_handles must NOT used with
            ## reset_handles at the same time!
            # self.reset_handles()
            self.loc = loc
            self.legendsize = legendsize
            self.anchor = anchor
            self.figtight = figtight
            self.legendalpha = legendalpha
            
        ell = (rho is not None)
        if (not ell):
            
            self.plot(x=x, y=y, yerr=yerr, xerr=xerr,
                      xisln=xisln, yisln=yisln,
                      fmt=fmt, ec=ec, elw=elw, els=els, # errorbar kw
                      capsize=capsize, barsabove=barsabove, # errorbar kw
                      c=c, alpha=alpha, label=label, **kwargs)
            self.append_handles(reset=False)
            
        else:

            self.plot(x=x, y=y, xisln=xisln, yisln=yisln,
                      fmt=fmt, c=c, alpha=alpha,
                      label=label, **kwargs)
            self.append_handles(reset=False)
            ## Ellipse errors
            self.eplot(x=x, y=y, mask=mask,
                       xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                       xisln=xisln, yisln=yisln,
                       sigmax=sigmax, sigmay=sigmay, rho=rho,
                       gammax=gammax, gammay=gammay,
                       ec=ec, elw=elw, els=els,
                       efill=efill, efc=efc, ehatch=ehatch,
                       errinlegend=errinlegend, alpha=alpha)
            self.append_handles()
            ## Replace line by patch in legend
            epatch = mpatches.Patch(ec=ec, lw=elw, ls=els,
                                    fill=efill, fc=efc, alpha=alpha,
                                    hatch=ehatch, label=errinlegend)
            self.handles[-1] = epatch

        ## set_legend
        if (self.loc is not None) and (label is not None):
            self.set_legend(loc=self.loc, fontsize=self.legendsize,
                            bbox_to_anchor=self.anchor, figtight=self.figtight,
                            framealpha=self.legendalpha)

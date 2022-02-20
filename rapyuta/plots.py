#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Plots

    plotrange, Btau, Ctau, ellipse, SUE
    plotool:
        set_clib, set_fig, set_ax, 
        reset_handles, append_handles, get_handles, set_legend, 
        plot, eplot, save, show, close
    pplot(plotool):
        add_plot, add_legend

"""

import warnings
from astropy import units as u
import numpy as np
from scipy import optimize
# import matplotlib as mpl
from matplotlib.ticker import (
    NullFormatter, ScalarFormatter, LogFormatter,
    LogFormatterExponent, LogFormatterSciNotation,
    PercentFormatter, FuncFormatter
)
import matplotlib.colors as mplc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

## Local
from utilities import InputError, merge_aliases
from arrays import arrayize, ramp

# cmap = mpl.cm.viridis
# norm = mpl.colors.Normalize(vmin=0, vmax=1)


##------------------------------
## Automatic plot range setting
##------------------------------
def plotrange(x,y,xran=None,yran=None,xlog=False,ylog=False,mask=None, \
              errx=None,erry=None,xisln=False,yisln=False):
    '''
    Automatically sets the x and y ranges for (X,Y) plots, based on the entered
    data set.

    Copyright: F. Galliano
    '''

    # Check the imput
    N = np.size(x)
    if (np.size(y) != N):
        UT.strike('plotrange','x and y should have the same size.')
    xran = arrayize(xran,N=2)
    yran = arrayize(yran,N=2)

    # X error bar settings
    if errx is not None:
        if (np.isscalar(errx)): errx = np.array([errx])
        sex = np.shape(errx)
        if (len(sex) == 2):
            if (sex != (2,N) ): UT.strike('plotrange','wrong size for errx.')
        elif (len(sex) == 1):
            if (sex != (N,) ): UT.strike('plotrange','wrong size for errx.')
            errx = np.array([errx,errx])
    else:
        errx = np.zeros((2,N))
        
    # Y error bar settings
    if erry is not None:
        if (np.isscalar(erry)): erry = np.array([erry])
        sey = np.shape(erry)
        if (len(sey) == 2):
            if (sey != (2,N) ): UT.strike('plotrange','wrong size for erry.')
        elif (len(sey) == 1):
            if (sey != (N,) ): UT.strike('plotrange','wrong size for erry.')
            erry = np.array([erry,erry])
    else:
        erry = np.zeros((2,N))
    
    # Homogenize the arrays and account for errors
    xlow = np.array(x,dtype=float).flatten() - errx[0,:]
    xhigh = xlow + errx[1,:]
    ylow = np.array(y,dtype=float).flatten() - erry[0,:]
    yhigh = ylow + erry[1,:]

    # Lin/Log
    if (xisln): xlow, xhigh = np.exp(xlow), np.exp(xhigh)
    if (yisln): ylow, yhigh = np.exp(ylow), np.exp(yhigh)
    
    # Mask
    mask = arrayize(mask,default=True,N=N)
    if (xlog): mask = ( mask & (xlow > 0) & (xhigh > 0) )
    if (ylog): mask = ( mask & (ylow > 0) & (yhigh > 0) )
    if (xran[0] != None): mask = ( mask & (xlow >= xran[0]) )
    if (xran[1] != None): mask = ( mask & (xhigh <= xran[1]) )
    if (yran[0] != None): mask = ( mask & (ylow >= yran[0]) )
    if (yran[1] != None): mask = ( mask & (yhigh <= yran[1]) )
            
    # Plain range
    xran = np.array([ np.min(xlow[mask]), np.max(xhigh[mask]) ])
    yran = np.array([ np.min(ylow[mask]), np.max(yhigh[mask]) ])

    # Add aesthetical margins
    fracmarg = 0.03
    if (not xlog):
        dxr = xran[1] - xran[0]
        xran += ( dxr*(-fracmarg), dxr*fracmarg )        
    else:
        dxr = xran[1] / xran[0]
        xran *= ( dxr**(-fracmarg), dxr**fracmarg )
    if (not ylog):
        dyr = yran[1] - yran[0]
        yran += ( dyr*(-fracmarg), dyr*fracmarg )
    else:
        dyr = yran[1] / yran[0]
        yran *= ( dyr**(-fracmarg), dyr**fracmarg )

    # Output
    return(xran,yran)

##-----------------------------------------------
##
##    Plotting functions for ellipses and SUEs
##
##           (Copyright: F. Galliano)
##
##-----------------------------------------------

## Function for SUEs
def Btau(tau):
    return( (1-2/np.pi)*(tau-1)**2 + tau )

## Function for SUEs
def Ctau(tau):
    return( ( (4/np.pi-1)*tau**2 + (3-8/np.pi)*tau + 4/np.pi-1 ) * (tau-1) )

## Ellipse
def ellipse(xmean=None,ymean=None,xstdev=None,ystdev=None,rho=None, \
            xmin=None,xmax=None,ymin=None,ymax=None,Npt=300, \
            xisln=False,yisln=False):
    """
    UNCERTAINTY ELLIPSES

    Function to plot uncertainty ellipses (or 1 sigma contour of a bivariate
    normal distribution). The parameters are the means (xmean,ymean), the 
    standard deviations (xstdev,ystdev) and the correlation coefficients (rho).
    The optional bounds (xmin,xmax,ymin,ymax) have the effect of truncating the
    ellipses in case there is a range of parameter space that is forbidden.

    It is important to notice that the xisln/yisln parameters are not related to
    the log settings of the axes where we plot the ellipse, but are here to 
    indicate that the moments of the variable to plot correspond to the natural 
    logarithm (ln) of the variable we want to display. For instance, for 
    displaying the ellipses of (x,y) where, for x, the moments are those of lnx,
    we would write:
        ellipse(xmean=mean_of_lnx,ymean=mean_of_y,xstdev=stdev_of_lnx, \
                ystdev=stdev_of_y,rho=correl_coeff_of_lnx_and_y,xisln=True)
    """
    x = ramp(x0=xmean-xstdev*(1-1.E-5),x1=xmean+xstdev*(1-1.E-5),N=Npt)
    c1 = rho * (x-xmean)/xstdev
    c2 = np.sqrt( (1-rho**2) * (1-(x-xmean)**2/xstdev**2) )
    y1 = ystdev * ( c1 - c2 ) + ymean
    y2 = ystdev * ( c1 + c2 ) + ymean
    xplot = np.concatenate((x,x[::-1],[x[0]]))
    yplot = np.concatenate((y1,y2[::-1],[y1[0]]))
    if (xisln): xplot = np.exp(xplot)
    if (yisln): yplot = np.exp(yplot)
    if (xmin != None): xplot[xplot < xmin] = xmin
    if (xmax != None): xplot[xplot > xmax] = xmax
    if (ymin != None): yplot[yplot < ymin] = ymin
    if (ymax != None): yplot[yplot > ymax] = ymax
    return(xplot,yplot)


## SUEs (1 sigma contour of a bivariate split-normal distribution)
def SUE(xmean=None,ymean=None,xstdev=None,ystdev=None,rho=None, \
        xskew=None,yskew=None,xmin=None,xmax=None,ymin=None,ymax=None, \
        Npt=300,xisln=False,yisln=False):
    """
    SKEWED UNCERTAINTY ELLIPSES (SUE)

    Function to plot uncertainty SUEs (or 1 sigma contour of a bivariate
    split-normal distribution). The parameters are the means (xmean,ymean), the 
    standard deviations (xstdev,ystdev), the skewnesses (xskew,yskew) and the 
    correlation coefficients (rho). The optional bounds (xmin,xmax,ymin,ymax) 
    have the effect of truncating the SUEs in case there is a range of 
    parameter space that is forbidden.

    It is important to notice that the xisln/yisln parameters are not related to
    the log settings of the axes where we plot the SUE, but are here to 
    indicate that the moments of the variable to plot correspond to the natural 
    logarithm (ln) of the variable we want to display. For instance, for 
    displaying the ellipses of (x,y) where, for x, the moments are those of lnx,
    we would write:
        SUE(xmean=mean_of_lnx,ymean=mean_of_y,xstdev=stdev_of_lnx, \
            ystdev=stdev_of_y,xskew=skewness_of_lnx,yskew=skewness_of_y, \
            rho=correl_coeff_of_lnx_and_y,xisln=True)
    """

    # Rotation angle
    theta = 1./2 * np.arctan( 2*rho*xstdev*ystdev / (xstdev**2-ystdev**2) )

    # Numerically solve for taux and tauy (tau=1.D2 ==> skew=0.99)
    taugrid = ramp(N=10000,x0=1.E-2,x1=1.E2,log=True)
    Ax = np.sqrt(np.pi/2) \
       * ( (np.cos(theta))**3*xskew*xstdev**3 \
         + (np.sin(theta))**3*yskew*ystdev**3 ) \
       / ( (np.sin(theta))**6 + (np.cos(theta))**6 ) \
       * ( ( (np.cos(theta))**2 - (np.sin(theta))**2 ) \
         / ( (np.cos(theta))**2*xstdev**2 \
           - (np.sin(theta))**2*ystdev**2 ) )**1.5
    Ay = np.sqrt(np.pi/2) \
       * ( (np.cos(theta))**3*yskew*ystdev**3 \
         - (np.sin(theta))**3*xskew*xstdev**3 ) \
       / ( (np.cos(theta))**6 + (np.sin(theta))**6 ) \
       * ( ( (np.cos(theta))**2 - (np.sin(theta))**2 ) \
         / ( (np.cos(theta))**2*ystdev**2 \
           - (np.sin(theta))**2*xstdev**2 ) )**1.5
    taux = np.exp(np.interp(Ax,Ctau(taugrid)/(Btau(taugrid))**1.5, \
                            np.log(taugrid)))
    tauy = np.exp(np.interp(Ay,Ctau(taugrid)/(Btau(taugrid))**1.5, \
                            np.log(taugrid)))
    if (not np.isfinite(taux) or taux > 1.E2): taux = 1.E2
    if (not np.isfinite(tauy) or tauy > 1.E2): tauy = 1.E2
 
    # Rest of the parameters
    lambdax = np.sqrt( ( (np.cos(theta))**2*xstdev**2 \
                       - (np.sin(theta))**2*ystdev**2 ) \
                     / ( (np.cos(theta))**2 - (np.sin(theta))**2 ) / Btau(taux) )
    lambday = np.sqrt( ( (np.cos(theta))**2*ystdev**2 \
                       - (np.sin(theta))**2*xstdev**2 ) \
                     / ( (np.cos(theta))**2 - (np.sin(theta))**2 ) / Btau(tauy) )
    x0 = xmean - np.sqrt(2/np.pi) * ( np.cos(theta)*lambdax*(taux-1) \
                                    - np.sin(theta)*lambday*(tauy-1) )
    y0 = ymean - np.sqrt(2/np.pi) * ( np.sin(theta)*lambdax*(taux-1) \
                                    + np.cos(theta)*lambday*(tauy-1) )

    # Draw the SUE
    matrot = np.array([ [ np.cos(theta), -np.sin(theta) ], \
                        [ np.sin(theta), np.cos(theta) ] ])
    xell_ax1 = np.zeros(2)
    yell_ax1 = np.zeros(2)
    xell_ax2 = np.zeros(2)
    yell_ax2 = np.zeros(2)
    for k in np.arange(4):
        if (k == 0):
            xell_sub = ramp(N=Npt,x0=-lambdax,x1=0) + x0
            rx = 1-(xell_sub-x0)**2/lambdax**2
            yell_sub = np.zeros(Npt)
            yell_sub[rx >= 0] = -lambday * np.sqrt(rx[rx >= 0]) + y0
            yell_sub[rx < 0] = np.nan
        elif (k == 1):
            xell_sub = ramp(N=Npt,x0=0,x1=lambdax*taux) + x0
            rx = 1-(xell_sub-x0)**2/lambdax**2/taux**2
            yell_sub = np.zeros(Npt)
            yell_sub[rx >= 0] = -lambday * np.sqrt(rx[rx >= 0]) + y0
            yell_sub[rx < 0] = np.nan
        elif (k == 2):
            xell_sub = (ramp(N=Npt,x0=0,x1=lambdax*taux))[::-1] + x0
            rx = 1-(xell_sub-x0)**2/lambdax**2/taux**2
            yell_sub = np.zeros(Npt)
            yell_sub[rx >= 0] = lambday*tauy * np.sqrt(rx[rx >= 0]) + y0
            yell_sub[rx < 0] = np.nan
        elif (k == 3):
            xell_sub = (ramp(N=Npt,x0=-lambdax,x1=0))[::-1] + x0
            rx = 1-(xell_sub-x0)**2/lambdax**2
            yell_sub = np.zeros(Npt)
            yell_sub[rx >= 0] = lambday*tauy * np.sqrt(rx[rx >= 0]) + y0
            yell_sub[rx < 0] = np.nan

        # Add the limit case (half ellipse)
        mask = np.logical_and(np.isfinite(yell_sub),np.isfinite(xell_sub))
        xell_sub = xell_sub[mask]
        yell_sub = yell_sub[mask]
        Nsub = np.count_nonzero(mask)

        # Rotate the ellipse
        for j in np.arange(Nsub):
            vecell = np.matmul(matrot, \
                               np.array([xell_sub[j]-x0,yell_sub[j]-y0]))
            xell_sub[j] = vecell[0] + x0
            yell_sub[j] = vecell[1] + y0
        if (k == 0):
            xell = xell_sub
            yell = yell_sub
        else:
            xell = np.concatenate((xell,xell_sub))
            yell = np.concatenate((yell,yell_sub))
    xplot = np.concatenate((xell,[xell[0]]))
    yplot = np.concatenate((yell,[yell[0]]))
                  
    # Logs and limits
    if (xisln):
        xplot = np.exp(xplot)
        x0 = np.exp(x0)
    if (yisln):
        yplot = np.exp(yplot)
        y0 = np.exp(y0)
    if (xmin != None):
        xplot[xplot < xmin] = xmin
        if (x0 < xmin): x0 = xmin
    if (xmax != None):
        xplot[xplot > xmax] = xmax
        if (x0 > xmax): x0 = xmax
    if (ymin != None):
        yplot[yplot < ymin] = ymin
        if (y0 < ymin): y0 = ymin
    if (ymax != None):
        yplot[yplot > ymax] = ymax
        if (y0 > ymax): y0 = ymax
    return(xplot,yplot,x0,y0)

##-----------------------------------------------
##
##            <plotool> based tools
##
##-----------------------------------------------

class plotool:
    '''
    plot Tool

    ------ INPUT ------
    nrows,ncols         Default: 1,1
    figint              Interactive figure (Default: True)

    plt.subplots(**kwargs)
    '''
    def __init__(self, nrows=1, ncols=1,
                 x=np.zeros(2), y=np.zeros(2), figint=True, **kwargs):
        '''
        ------ self ------
        figid               ID number of current figure (Default: 0)
        horder              ID number of current handle in self.handles (Default: 0)
        handles             Current handle list (Default: [])
        labels              Current label list (Default: [])
        legend              Current legend (Default: None)
        trans               transformation of coordinate system
                              'Data'
                              'Axes' (default)
                              'Figure'
                              'Display'
        '''
        ## INPUTS
        self.figid = 0
        self.horder = 0
        self.handles = []
        self.labels = []
        self.legend = None
        self.trans = 'Axes'
        self.x = x
        self.y = y

        if figint:
            plt.ion()

        self.nrows = nrows
        self.ncols = ncols

        self.fig, self.axes = plt.subplots(nrows, ncols, squeeze=False,
                                           num=self.figid, **kwargs)
        ## De-squeeze/Inflate (if squeeze=True)
        # if nrows==1 and ncols==1:
        #     self.axes = np.array(self.axes)[np.newaxis,np.newaxis]
        # elif nrows==1:
        #     self.axes = self.axes[np.newaxis,:]
        # elif ncols==1:
        #     self.axes = self.axes[:,np.newaxis]

        self.ax = self.axes[0,0]

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
        
    def set_fig(self, left=None, right=None, bottom=None, top=None,
                wspace=None, hspace=None, title=None, tsize=None):

        self.fig.subplots_adjust(left=left, right=right,
            bottom=bottom, top=top, wspace=wspace, hspace=hspace)
        
        if title is not None:
            self.fig.suptitle(title, size=tsize)
        
    def set_ax(self, subpos=None, # ax = axes[subpos[0],subpos[1]]
               xlog=False, ylog=False, # ax.set_xscale
               basex=10, basey=10, nonposx='clip', nonposy='clip', # ax.set_xscale
               xlim=(None,None), ylim=(None,None), #ax.set_xlim
               xtk=None, xtkmi=None, xtkform=None, xtksize=None, # ax.xaxis.set_tick_params
               ytk=None, ytkmi=None, ytkform=None, ytksize=None, 
               xlabel=None, xsize=None, ylabel=None, ysize=None, # ax.set_xlabel
               xtktoggle=False, ytktoggle=False,
               title=None, tsize=None, # ax.set_title (subplot title)
               **kwargs):
        '''
        nonposx, nonposy: 'sym', 'mask', 'clip'

        ------ INPUT ------
        subpos              identify subplot (Default: None - current self.ax)
        xtk,ytk             list of major tick labels
        xtkmi,ytkmi         list of minor tick labels
        xtkform,ytkform     formatter of (major) tick labels
        xtksize,ytksize     size of (major) tick labels
        xsize,ysize         size of axes labels
        xtktoggle           toggle x tick and label to top
        ytktoggle           toggle y tick and label to right

        self.ax.tick_params(**kwargs):
        axis                {'x', 'y', 'both' (default)}
        reset               Default: False
        which               {'major' (dafault), 'minor', 'both'}
        direction           {'in', 'out' (default), 'inout'}
        length              tick length
        width               tick width
        '''
        if self.nrows!=1 or self.ncols!=1:
            if subpos is not None:
                self.ax = self.axes[subpos[0],subpos[1]]
            ## else: keep current self.ax (Default: (0,0))

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

        ## x ticks
        if xtkform=='log':
            xformat = LogFormatter()
        elif xtkform=='log_exp':
            xformat = LogFormatterExponent()
        elif xtkform=='log_sci':
            xformat = LogFormatterSciNotation()
        elif xtkform=='mylog':
            xformat = FuncFormatter(
                lambda a,pos: (
                    '{{:.{:1d}f}}'.format(int(np.maximum(-np.ma.log10(np.array([a])).filled(0),0)))
                ).format(a)
            )
        elif xtkform=='pct':
            xformat = PercentFormatter()
        else:
            xformat = ScalarFormatter()
        if xtk is not None:
            self.ax.set_xticks(xtk, minor=False) # major
        if xtkmi is not None:
            self.ax.set_xticks(xtkmi, minor=True) # minor
        self.ax.xaxis.set_major_formatter(xformat) # major
        self.ax.xaxis.set_minor_formatter(NullFormatter()) # minor
        self.ax.xaxis.set_tick_params(labelsize=xtksize)
        if xtktoggle:
            self.ax.xaxis.tick_top()
            self.ax.xaxis.set_label_position('top')
        ## y ticks
        if ytkform=='log':
            yformat = LogFormatter()
        elif ytkform=='log_exp':
            yformat = LogFormatterExponent()
        elif ytkform=='log_sci':
            yformat = LogFormatterSciNotation()
        elif ytkform=='mylog':
            yformat = FuncFormatter(
                lambda a,pos: (
                    '{{:.{:1d}f}}'.format(int(np.maximum(-np.ma.log10(np.array([a])).filled(0),0)))
                ).format(a)
            )
        elif ytkform=='pct':
            yformat = PercentFormatter()
        else:
            yformat = ScalarFormatter()
        if ytk is not None:
            self.ax.set_yticks(ytk, minor=False) # major
        if ytkmi is not None:
            self.ax.set_yticks(ytkmi, minor=True) # minor
        self.ax.yaxis.set_major_formatter(yformat) # major
        self.ax.yaxis.set_minor_formatter(NullFormatter()) # minor
        self.ax.yaxis.set_tick_params(labelsize=ytksize)
        if ytktoggle:
            self.ax.yaxis.tick_right()
            self.ax.yaxis.set_label_position('right')
        ## both ticks
        self.ax.tick_params(**kwargs)
        
        if xlabel is not None:
            self.ax.set_xlabel(xlabel, size=xsize)
        if ylabel is not None:
            self.ax.set_ylabel(ylabel, size=ysize)
        if title is not None:
            self.ax.set_title(title, size=tsize)

    def reset_handles(self):
        '''
        Reset handles (before plot with new legend)
        '''
        if self.trans=='Figure':
            self.fig.add_artist(self.legend)
        elif self.trans=='Axes' or self.trans=='Data':
            self.ax.add_artist(self.legend)
        self.handles = []
        self.horder = 0

        return self.handles

    def append_handles(self):
        '''
        Append currently added handle (after plot)
        '''
        handles, labels = self.ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if label not in self.labels:
                self.handles.append(handle)
                self.labels.append(label)
                self.horder += 1

        return self.handles
        
    def get_handles(self):
        '''
        Get non-repeated handles (after plot)
        '''
        handles, labels = self.ax.get_legend_handles_labels()
        handlist, labelist = [], []
        for handle, label in zip(handles, labels):
            if label not in labelist:
                handlist.append(handle)
                labelist.append(label)
        self.horder = len(handlist)
        self.handles = handlist
        self.labels = labelist

        return self.handles
        
    def set_legend(self, subpos=None,
                   left=None, right=None, bottom=None, top=None, figtight=False,
                   handles=None, **kwargs):
        '''
        - bbox_to_anchor rules: (1,1) correspond to upper right of the axis

                    bbox_to_anchor = (1,1)
          .--------.
          |        |
          |  Axes  |
          | (bbox) |
          |        |
          .--------.
        
        - lengend loc is relative to the reference point (bbox_to_anchor) as follows:

                               |
          lower right (2)      |       lower left (1)
          ---------------bbox_to_anchor---------------
          upper right (3)      |       upper left (4)
                               |
        
        ------ INPUT ------
        subpos              identify subplot (Default: None - current self.ax)
        x0,y0               bottom left of current Axes (Default: None)
        shrinkx,shrinky     current Axes size (Default: 1)
        figtight            ignore Axes settings (Default: False)
        handles             Default: None - self.handles

        self.fig/ax.legend(**kwargs):
        title               title of legend box
        loc                 relative position of legend box
        bbox_to_anchor      reference point of legend box
        bbox_transform      reference frame of coordinates
                              self.fig.transFigure
                              self.ax.transAxes (default)
        ------ OUTPUT ------
        self.legend         registrated with self.fig/ax.add_artist(self.legend)
        '''
        if handles is None:
            handles = self.handles

        if self.trans=='Figure':
            self.fig.subplots_adjust(left=left, right=right,
                                     bottom=bottom, top=top)
            
            self.legend = self.fig.legend(handles=handles, **kwargs)
        elif self.trans=='Axes' or self.trans=='Data':
            if self.nrows!=1 or self.ncols!=1:
                if subpos is not None:
                    self.ax = self.axes[subpos[0],subpos[1]]
                ## else: keep current self.ax (Default: (0,0))
            
            ## Set current axes
            bbox = self.ax.get_position()
            if left is None:
                x0 = bbox.x0
            else:
                x0 = bbox.x0 + left*(bbox.x1-bbox.x0)
            if right is None:
                x1 = bbox.x1
            else:
                x1 = bbox.x0 + right*(bbox.x1-bbox.x0)
            if bottom is None:
                y0 = bbox.y0
            else:
                y0 = bbox.y0 + bottom*(bbox.y1-bbox.y0)
            if top is None:
                y1 = bbox.y1
            else:
                y1 = bbox.y0 + top*(bbox.y1-bbox.y0)
            self.ax.set_position([x0, y0, x1-x0, y1-y0])

            self.legend = self.ax.legend(handles=handles, **kwargs)
        else:
            raise InputError('<plotool.set_legend>',
                             'Non-recognized transformation!')

        if figtight:
            self.fig.tight_layout()

        return self.legend
        
    def plot(self, x=None, y=None, yerr=None, xerr=None, xisln=False, yisln=False,
             fmt='', capsize=None, barsabove=False, # errorbar kw
             ecolor=None, ec=None, elinewidth=None, elw=None, # errorbar kw
             subpos=None, mod='CA', **kwargs):
        '''
        Like set_ax(), this is a clump operation.
        The idea is to all set in one command,
        while each single operation should also be valid.

        ------ INPUT ------
        subpos              identify subplot (Default: None - current self.ax)

        self.ax.errorbar(**kwargs)
        '''
        if self.nrows!=1 or self.ncols!=1:
            if subpos is not None:
                self.ax = self.axes[subpos[0],subpos[1]]
            ## else: keep current self.ax (Default: (0,0))

        ## Keyword aliases
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

        ## Log inputs
        xp, xperr, yp, yperr = None, None, None, None
        if (xisln):
            if x is not None: xp = np.exp(x)
            if xerr is not None: xperr = x * (1. - np.exp(-xerr)) ## suppose xmin <-> lnxmin
        else:
            if x is not None: xp = x
            if xerr is not None: xperr = xerr
        if (yisln):
            if y is not None: yp = np.exp(y)
            if yerr is not None: yperr = y * (1. - np.exp(-yerr))
        else:
            if y is not None: yp = y
            if yerr is not None: yperr = yerr

        ## CA: Cartesian using matplotlib.pyplot.errorbar
        if mod=='CA':
                
            self.markers, self.caps, self.bars = self.ax.errorbar(
                x=xp, y=yp, yerr=yperr, xerr=xperr,
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

    def eplot(self, x=None, y=None, mask=None,
              xmin=None, xmax=None, ymin=None, ymax=None,
              xisln=False, yisln=False,
              ## Uncertainty kw
              sigmax=None, sigmay=None, rho=None,
              gammax=None, gammay=None,
              ## Bar/ellipse keywords
              edgecolor=None, ecolor=None, ec=None,
              elinewidth=None, elw=None, elinestyle=None, els=None,
              efillcolor=None, efc=None, efill=False, ehatch=None,
              errinlegend=None, alpha=1,
              ## Other kw
              subpos=None, **kwargs):
        '''
        DISPLAY ERROR BARS/ELLIPSES/SUES
        
        x and y are (N,) shape arrays. sigmax and sigmay are either (N,) or (2,N)
        shape arrays. If xisln is set, then it is assumed that all the x related
        quantities are moments of lnx (x is mu(lnx), sigmax is stdev(lnx), etc.).
        This is independent of the xlog setting.

        If rho is an (N,) shape array, an ellipse
        is drawn, instead. If gammax or gammay are not None, an "asymmetric 
        ellipse" or SUE (1 sigma contour of a bivariate split-normal 
        distribution) is drawn. 

        ------ INPUT ------
        subpos              identify subplot (Default: None - current self.ax)

        self.plot(**kwargs)
        '''
        if self.nrows!=1 or self.ncols!=1:
            if subpos is not None:
                self.ax = self.axes[subpos[0],subpos[1]]
            ## else: keep current self.ax (Default: (0,0))
            
        ## kw aliases
        ec = merge_aliases(None, edgecolor=edgecolor, ecolor=ecolor, ec=ec)
        elw = merge_aliases(None, elinewidth=elinewidth, elw=elw)
        els = merge_aliases(None, elinestyle=elinestyle, els=els)
        efc = merge_aliases(None, efillcolor=efillcolor, efc=efc)
        
        ## Central values
        xp, N = arrayize(x)
        yp = arrayize(y,N=N)
        if (N != np.size(yp)):
            raise InputError('<plotool.eplot>',
                             'x and y must have the same size.')
            
        ## Ellipse/SUE
        ell = (rho is not None)
        skewll = ( ( (gammax is not None) or (gammay is not None) ) \
                   and ell )

        ## X error bar settings
        uncx = False
        if sigmax is not None:
            uncx = True
            if (np.isscalar(sigmax)): sigmax = np.array([sigmax])
            sex = np.shape(sigmax)
            if (len(sex) == 2):
                if (sex != (2,N) ):
                    raise InputError('<plotool.eplot>',
                                     'wrong size for errx.')
            elif (len(sex) == 1):
                if (sex != (N,) ):
                    raise InputError('<plotool.eplot>',
                                     'wrong size for errx.')
                sigmax = np.array([sigmax,sigmax])

        ## Y error bar settings
        uncy = False
        if sigmay is not None:
            uncy = True
            if (np.isscalar(sigmay)): sigmay = np.array([sigmay])
            sey = np.shape(sigmay)
            if (len(sey) == 2):
                if (sey != (2,N) ):
                    raise InputError('<plotool.eplot>',
                                     'wrong size for erry.')
            elif (len(sey) == 1):
                if (sey != (N,) ):
                    raise InputError('<plotool.eplot>',
                                     'wrong size for erry.')
                sigmay = np.array([sigmay,sigmay])
        
        if (ell):
            ## Setting the error bars and the potential lower/upper limits
            if (not uncx or not uncy):
                raise InputError('<plotool.eplot>',
                                 'for ellipses, both x and y errors must be set.')
            xcenell = xp.copy()
            xerr = sigmax[0,:]
            if (xisln):
                xcen = np.exp(xp)
            else:
                xcen = xp.copy()
            ycenell = yp.copy()
            yerr = sigmay[0,:]
            if (yisln):
                ycen = np.exp(yp)
            else:
                ycen = yp.copy()
                
            ## Mask
            if mask is None:
                mask = np.ones(N,bool)
            Nsub = np.count_nonzero(mask)
            xcen = xcen[mask]
            ycen = ycen[mask]
            xcenell = xcenell[mask]
            ycenell = ycenell[mask]
            xerr = xerr[mask]
            yerr = yerr[mask]
            rho = arrayize(rho,N=N)
            rho = rho[mask]
            if (skewll):
                gammax = arrayize(gammax,N=N)
                gammay = arrayize(gammay,N=N)
                gammax = gammax[mask]
                gammay = gammay[mask]

            ## Lines
            elw = arrayize(elw,N=N)
            elw = elw[mask]
            els = arrayize(els,N=N)
            els = els[mask]

            ## Colors
            efill = arrayize(efill,N=N)
            if ec is None:
                ec = np.array(['b' for i in range(N)])
                efc = np.array(['lightblue' for i in range(N)])
            else:
                ec = arrayize(ec,N=N)
                efc = arrayize(efc,N=N)
            ec = ec[mask]
            efc = efc[mask]

            ## Hatching the ellipses
            if (ehatch):
                ehatch = np.array(['back slanted' for i in np.arange(N)])
            elif (not ehatch):
                ehatch = np.array(['' for i in np.arange(N)])
            else:
                ehatch = arrayize(ehatch,default='',N=N)
            ehatch = ehatch[mask]
            nothatched = (ehatch == '')
            efill = np.logical_or(efill,(ehatch != ''))
            
            ## Draw the error bars
            if (not skewll):

                for i in np.arange(Nsub):
                    if i==0:
                        elegend = errinlegend
                    else:
                        elegend = None
                    xell, yell = ellipse(xmean=xcenell[i],ymean=ycenell[i],
                                         xstdev=xerr[i],ystdev=yerr[i],
                                         rho=rho[i],xmin=xmin,xmax=xmax,
                                         ymin=ymin,ymax=ymax,xisln=xisln,yisln=yisln)
                    self.plot(xell, yell,
                              color=ec[i], linewidth=elw[i], linestyle=els[i],
                              label=elegend, alpha=alpha, **kwargs)
                    if (efill[i]):
                        self.ax.fill(xell,yell,hatch=ehatch[i],fill=nothatched[i],
                                     color=efc[i],alpha=alpha)

            elif(skewll):

                for i in np.arange(Nsub):
                    if i==0:
                        elegend = errinlegend
                    else:
                        elegend = None
                    xell, yell, x0, y0 \
                        = SUE(xmean=xcenell[i],ymean=ycenell[i],
                              xstdev=xerr[i],ystdev=yerr[i],
                              xskew=gammax[i],yskew=gammay[i],
                              rho=rho[i],xmin=xmin,xmax=xmax,
                              ymin=ymin,ymax=ymax,
                              xisln=xisln,yisln=yisln)
                    xcen[i] = x0
                    ycen[i] = y0
                    self.plot(xell, yell,
                              color=ec[i], linewidth=elw[i], linestyle=els[i],
                              label=elegend, alpha=alpha, **kwargs)
                    if (efill[i]):
                        self.ax.fill(xell,yell,hatch=ehatch[i],fill=nothatched[i],
                                     color=efc[i],alpha=alpha)
            
    def save(self, savename=None, transparent=False,
             figtight=False, close=True, **kwargs):
        '''

        ------ INPUT ------
        self.fig.savefig(**kwargs)
        '''
        ## Use self.fig.set_size_inches to define ranges
        if figtight:
            bbox = 'tight'
        else:
            bbox = None
            
        if savename is not None:
            self.fig.savefig(savename, transparent=transparent,
                             bbox_inches=bbox, **kwargs)
        else:
            warnings.warn('Not saved! ')

        if (close): plt.close(self.figid)

    def show(self):

        plt.ioff()
        self.fig.show(self.fig)

    def close(self):

        plt.close(self.figid)

class pplot(plotool):
    '''
    Uni-frame plot (1 row * 1 col)

    ------ INPUT ------
    plotool.plot(**kwargs)
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
                 clib='base', c=None, **kwargs):
        super().__init__(figsize=figsize, figint=figint, x=x, y=y)

        self.pltid = 0

        ## Keyword aliases
        ## Note that ec and ecolor are shared by plot (errorbar) and eplot (ellipse)
        ec = merge_aliases('grey', edgecolor=edgecolor, ecolor=ecolor, ec=ec) # Default: 'grey'
        elw = merge_aliases(None, elinewidth=elinewidth, elw=elw)
        els = merge_aliases(None, elinestyle=elinestyle, els=els)
        efc = merge_aliases(None, efillcolor=efillcolor, efc=efc)

        ## Auto color
        self.set_clib(clib)
        if c is None:
            c = self.clib[self.pltid]

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
        plotool.plot(**kwargs)
        '''
        self.pltid += 1

        ## Keyword aliases
        ## Note that ec and ecolor are shared by plot (errorbar) and eplot (ellipse)
        ec = merge_aliases('grey', edgecolor=edgecolor, ecolor=ecolor, ec=ec) # Default: 'grey'
        elw = merge_aliases(None, elinewidth=elinewidth, elw=elw)
        els = merge_aliases(None, elinestyle=elinestyle, els=els)
        efc = merge_aliases(None, efillcolor=efillcolor, efc=efc)

        ## Auto color
        if self.pltid==len(self.clib):
            self.pltid = 0
        if c is None:
            c = self.clib[self.pltid]
        
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
            self.reset_handles()
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
            self.append_handles()
            
        else:

            self.plot(x=x, y=y, xisln=xisln, yisln=yisln,
                      fmt=fmt, c=c, alpha=alpha,
                      label=label, **kwargs)
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
            epatch = mpatches.Patch(ec=ec, lw=elw, ls=els,
                                    fill=efill, fc=efc, alpha=alpha,
                                    hatch=ehatch, label=errinlegend)
            self.handles[-1] = epatch

        ## set_legend
        if (self.loc is not None) and (label is not None):
            self.set_legend(loc=self.loc, fontsize=self.legendsize,
                            bbox_to_anchor=self.anchor, figtight=self.figtight,
                            framealpha=self.legendalpha)

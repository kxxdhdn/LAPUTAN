#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""

Input & Output

    fclean, write_fits, read_fits, write_hdf5, read_hdf5,
    write_ascii, read_ascii, write_csv, read_csv

"""

import sys, logging
logging.disable(sys.maxsize)

import subprocess as SP
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import h5py as H5
import csv

global fitsext, h5ext, ascext, csvext, savext
fitsext = '.fits'
h5ext = '.h5'
ascext = '.txt'
csvext = '.csv'
savext = '.sav'

def fclean(file, *alert):
    '''
    Clean folder/files
    '''
    SP.call('rm -rf '+file, shell=True)
    for text in alert:
        print(text)

def write_fits(file, header, data, wave=None, wmod=0, **hdrl):
    '''
    Write fits file

    ------ INPUT ------
    file                output FITS filename
    header              header of primary HDU
    data                data in primary HDU
    wave                data in table 1 (ndarray. Default: None)
    wmod                wave table format (0 - Image; 1 - BinTable. Default: 0)
    ------ OUTPUT ------
    '''
    for key, value in hdrl.items():
        header[key] = value
    primary_hdu = fits.PrimaryHDU(header=header, data=data)
    hdul = fits.HDUList(primary_hdu)
    
    ## Add table
    if wave is not None:
        ## Convert wave format
        if isinstance(wave, fits.fitsrec.FITS_rec):
            if wmod==0:
                wave = wave[0][0][:,0]
        else:
            Nw = len(wave)
            if wmod==1:
                wave = np.array(wave).reshape((Nw,1))
                col = fits.Column(array=[wave], format=str(Nw)+'E',
                                  name='WAVE-TAB', unit='um', dim='(1,{})'.format(Nw))
                tab = fits.BinTableHDU.from_columns([col], name='WCS-TAB ')
                wave = tab.data
        ## Create table
        if wmod==0:
            hdu = fits.ImageHDU(data=wave, name='WAVE-TAB')
        elif wmod==1:
            hdu = fits.BinTableHDU(data=wave, name='WCS-TAB ')

        hdul.append(hdu)

    hdul.writeto(file+fitsext, overwrite=True)
    
def read_fits(file, file_unc=None, wmod=0):
    '''
    Read fits file (auto detect dim)

    ------ INPUT ------
    file                input FITS filename
    file_unc            input uncertainty file
    wmod                output wave mode (Default: 0)
                          0 - 1darray; 
                          1 - FITS_rec.
    ------ OUTPUT ------
    ds                  output dataset
      header              header of primary HDU
      data                data in primary HDU
      header_w            header of W-TAB
      wave                data in table 1 (None if does not exist)
      unc                 uncertainty array
    '''
    ## Initialize output object
    ds = type('', (), {})()
    ds.header_w = None
    ds.wave = None
    ds.unc = None
    
    ## Read header & data
    with fits.open(file+fitsext) as hdul:
        ds.HDUL = hdul
        hdr = hdul[0].header
        ds.data = hdul[0].data
        ds.header = hdr
        try:
            ds.wcs = WCS(hdr)
        except:
            ds.wcs = WCS(None)
            # print('Not supported WCS (probably with 3D SIP distortion), use astrom.fixwcs')

        ## Read wavelength
        if len(hdul)==2:
            ds.header_w = hdul[1].header
            wave = hdul[1].data

            if isinstance(hdul[1], fits.BinTableHDU):
                if wmod==0:
                    wave = wave[0][0][:,0] ## Convert FITS_rec to 1darray
            elif isinstance(hdul[1], fits.ImageHDU):
                Nw = len(wave)
                if wmod==1:
                    wave = np.array(wave).reshape((Nw,1))
                    col = fits.Column(array=[wave], format=str(Nw)+'E',
                                      name='WAVE-TAB', unit='um', dim='(1,{})'.format(Nw))
                    tab = fits.BinTableHDU.from_columns([col], name='WCS-TAB ')
                    wave = tab.data
            
            ds.wave = wave
    
    if file_unc is not None:
        ## Read uncertainty data
        with fits.open(file_unc+fitsext) as hdul:
            ds.unc = hdul[0].data

    return ds

def write_hdf5(file, name, data, group='/',
               ind1=None, ind2=None, ind3=None, ind4=None,
               ind5=None, ind6=None, ind7=None,
               append=False, verbose=False):
    '''
    Write dataset into a h5 file (a single name/data_array per time, dim < 7)
    Inspired by SwING inout library

    ------ INPUT ------
    file                output h5 filename
    name                name of the dataset (len <= 80)
    data                dataset (dim < 7)
    group               name of the group (Default: '/')
    indx                array index ([idimx_inf,idimx_sup] or idimx if data is scalar, Default: None)
    append              True: if not overwrite (Default: False)
    verbose             courtesy notification (Default: False)
    ------ OUTPUT ------
    '''
    ## Preliminaries
    ##---------------
    if (append):
        hf = H5.File(file+h5ext, 'a')
    else:
        hf = H5.File(file+h5ext, 'w')
    
    h5 = hf.require_group(group) # Default: group = '/'

    ## Check if it is a subarray
    subarr = ( (ind1 is not None) or (ind2 is not None) or (ind3 is not None) or \
               (ind4 is not None) or (ind5 is not None) or (ind6 is not None) )

    ## Convert all data to array format
    darr = np.array(data)
    sharr = darr.shape
    if darr.dtype.kind == 'U':
        ## create_dataset does not support lists of UTF-8 yet
        ## "these strings are supposed to store only ASCII-encoded text"
        ## See http://docs.h5py.org/en/stable/strings.html
        asciiList = [n.encode('ascii','ignore') for n in darr.flatten()]
        darr = np.array(asciiList).reshape(sharr)

    ## Write dataset (dset)
    ##----------------------
    
    ## Case 0: check if the dataset already exists
    try:
        dset = h5[name]
        hasdset = True
    except:
        hasdset = False
    if (hasdset):
        createdset = False
        ## No indx input
        if (not subarr):
            del h5[name]
            createdset = True
    else:
        createdset = True
        if (subarr):
            h5.close()
            raise ValueError('Dataset does not exist')

    ## Case 1: no dataset OR no subarr
    if (createdset):
        dset = h5.create_dataset(name, data=darr)

    ## Case 2: sub-array filling
    if (subarr):
        ## Dimension of the array
        Ndim = len(sharr)

        ## Indices indx[idimx_inf,idimx_sup]
        if (Ndim > 0):
            if ind1 is None:
                ind1 = [0,sharr[0]]
            elif (np.size(ind1) == 1):
                ind1 = [ind1, ind1+1]
        if (Ndim > 1):
            if ind2 is None:
                ind2 = [0,sharr[1]]
            elif (np.size(ind2) == 1):
                ind2 = [ind2, ind2+1]
        if (Ndim > 2):
            if ind3 is None:
                ind3 = [0,sharr[2]]
            elif (np.size(ind3) == 1):
                ind3 = [ind3, ind3+1]
        if (Ndim > 3):
            if ind4 is None:
                ind4 = [0,sharr[3]]
            elif (np.size(ind4) == 1):
                ind4 = [ind4, ind4+1]
        if (Ndim > 4):
            if ind5 is None:
                ind5 = [0,sharr[4]]
            elif (np.size(ind5) == 1):
                ind5 = [ind5, ind5+1]
        if (Ndim > 5):
            if ind6 is None:
                ind6 = [0,sharr[5]]
            elif (np.size(ind6) == 1):
                ind6 = [ind6, ind6+1]
        if (Ndim > 6):
            if ind7 is None:
                ind7 = [0,sharr[6]]
            elif (np.size(ind7) == 1):
                ind7 = [ind7, ind7+1]

        ## Write the sub-array if indx are set
        if Ndim == 0:
            dset = arr[0]
        elif Ndim == 1:
            dset[ind1[0]:ind1[1]] = np.reshape( darr, (ind1[1]-ind1[0],) )
        elif Ndim == 2:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0]) )
        elif Ndim == 3:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0],ind3[0]-ind3[1]) )
        elif Ndim == 4:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1],ind4[0]:ind4[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1]) )
        elif Ndim == 5:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                 ind4[0]:ind4[1],ind5[0]:ind5[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1], \
                                   ind5[0]-ind5[1]) )
        elif Ndim == 6:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                 ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1], \
                                   ind5[0]-ind5[1],ind6[0]-ind6[1]) )
        elif Ndim == 7:
            dset[ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                 ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1],ind7[0]:ind7[1]] = \
                np.reshape( darr, (ind1[1]-ind1[0],ind2[1]-ind2[0], \
                                   ind3[0]-ind3[1],ind4[0]-ind4[1], \
                                   ind5[0]-ind5[1],ind6[0]-ind6[1],ind7[0]-ind7[1]) )

    hf.flush()
    hf.close()
    
    ## Courtesy notification
    ##-----------------------
    if (verbose):
        if (append):
            if (subarr):
                print('[write_hdf5] Dataset {}/{} in the file:'.format(group,name))
                print('    ', file+h5ext, ' has been modified.')
            else:
                print('[write_hdf5] Dataset {}/{} has been added in the file:'.format(group,name))
                print('    ', file+h5ext, '.')
        else:
            print('[write_hdf5] Dataset {}/{} has been written in the new file:'.format(group,name))
            print('    ', file+h5ext, '.')

def read_hdf5(file, name, group='/',
              ind1=None, ind2=None, ind3=None, ind4=None,
              ind5=None, ind6=None, ind7=None):
    '''
    Read h5 file (a single name/data_array per time, dim < 7)
    Inspired by SwING inout library

    ------ INPUT ------
    file                input h5 filename
    name                name of the dataset (len <= 80)
    group               name of the group (Default: '/')
    indx                array index ([idimx_inf,idimx_sup] or idimx if data is scalar, Default: None)
    ------ OUTPUT ------
    dset                dataset
    '''
    ## Preliminaries
    ##---------------
    hf = H5.File(file+h5ext, 'r')
    h5 = hf.require_group(group) # Default: group = '/'

    Ndim = h5[name].ndim
    shapIN = h5[name].shape

    ## Read dataset (dset)
    ##----------------------
    
    ## Indices indx[idimx_inf,idimx_sup]
    if (Ndim > 0):
        if ind1 is None:
            ind1 = [0,shapIN[0]]
        elif (np.size(ind1) == 1):
            ind1 = [ind1, ind1+1]
    if (Ndim > 1):
        if ind2 is None:
            ind2 = [0,shapIN[1]]
        elif (np.size(ind2) == 1):
            ind2 = [ind2, ind2+1]
    if (Ndim > 2):
        if ind3 is None:
            ind3 = [0,shapIN[2]]
        elif (np.size(ind3) == 1):
            ind3 = [ind3, ind3+1]
    if (Ndim > 3):
        if ind4 is None:
            ind4 = [0,shapIN[3]]
        elif (np.size(ind4) == 1):
            ind4 = [ind4, ind4+1]
    if (Ndim > 4):
        if ind5 is None:
            ind5 = [0,shapIN[4]]
        elif (np.size(ind5) == 1):
            ind5 = [ind5, ind5+1]
    if (Ndim > 5):
        if ind6 is None:
            ind6 = [0,shapIN[5]]
        elif (np.size(ind6) == 1):
            ind6 = [ind6, ind6+1]
    if (Ndim > 6):
        if ind7 is None:
            ind7 = [0,shapIN[6]]
        elif (np.size(ind7) == 1):
            ind7 = [ind7, ind7+1]
        
    # Read the array or the sub-array if some indx are set
    if Ndim==0:
        dset = h5[name]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==1:
        dset = h5[name][ind1[0]:ind1[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==2:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==3:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==4:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==5:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1],ind5[0]:ind5[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==6:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')
    elif Ndim==7:
        dset = h5[name][ind1[0]:ind1[1],ind2[0]:ind2[1],ind3[0]:ind3[1], \
                        ind4[0]:ind4[1],ind5[0]:ind5[1],ind6[0]:ind6[1], \
                        ind7[0]:ind7[1]]
        if (str(dset.dtype).startswith('|S')):
            dset = np.array(dset, dtype='unicode')

    hf.close()

    return dset

def write_ascii(file, header=None, dset=None, trans=False,
                ascext=ascext, append=False, comment=None):
    '''
    Write ASCII file
    Supported format: commented_header
    See also: astropy.io.ascii.write, numpy.savetxt
    
    ------ INPUT ------
    file                input ASCII filename
    header              data header
    dset                dataset (1darray[nrow], 2darray[nrow,ncol] or list)
    trans               transpose dset (Default: False)
    ascext              ASCII file suffix (Default: '.txt')
    append              True: if not overwrite (Default: False)
    comment             single line comment on top
    ------ OUTPUT ------
    '''
    if append==True:
        mod = 'a'
    else:
        mod = 'w'

    with open(file+ascext, mod) as f:
        pass
        ## Comment on the top
        ##--------------------
        if comment is not None:
            f.write('## '+comment+'\n')
        
        ## Number of columns/rows
        ncol = 0
        nrow = 0
        ## Dataset
        ##---------
        rows = ''
        if dset is not None:
            dset = np.array(dset)
            ## Transpose dset
            if trans==True:
                Ny, Nx = dset.shape
                darr = []
                for x in range(Nx):
                    for y in range(Ny):
                        darr.append(dset[y,x])
                darr = np.array(darr).reshape((Nx,Ny))
            else:
                darr = dset
            ## Update number of columns/rows
            Ndim = darr.ndim
            if Ndim==1:
                ncol = 1
                nrow = len(darr)
            else:
                nrow, ncol = darr.shape
            ## Convert dataset array to string
            if Ndim==1 or Ndim==2:
                if nrow<500:
                    rows += str(darr).replace('[','').replace(']','')+'\n'
                else:
                    ## Python 500 row limit
                    for i in range(nrow//500+1):
                        rows += str(darr[500*i:500*(i+1)]).replace('[','').replace(']','')+'\n'
            else:
                raise ValueError(
                    'Expected 1D or 2D array, got {}D array instead!'.format(Ndim))
        ## Header
        ##--------
        hrow = ''
        if header is not None:
            if isinstance(header, str):
                header = [header]
            else:
                header = list(header)
                
            for i, h in enumerate(header):
                if i==0:
                    hrow += h
                else:
                    hrow += ' '+h
            hrow += '\n'
            ## Update number of columns/rows
            ncol = len(header)
            nrow += 1
        ## Column width control
        ##----------------------
        ## Table in a single string
        if nrow==0:
            fwrows = ''
        else:
            cells = hrow.split()
            cells.extend(rows.split())
            cells = np.array(cells).reshape((-1,ncol))
            col_width = []
            for x in range(ncol):
               col_width.append(len(max(cells[:,x], key=len)))
            ## Table with fixed column widths in a single string
            if hrow=='':
                fwrows = '' # no header
            else:
                fwrows = '# '
            for y in range(nrow):
                for x in range(ncol):
                    c = cells[y,x]
                    ## '{message:{fill}{align}{width}}'.format(
                    ##     message=c, fill=' ', align='>', width=col_width+1)
                    if y==0 and x==0 and header is not None:
                        ## Due to hashtag
                        fwrows += ''.ljust(col_width[x]-len(c))+c
                    else:
                        fwrows += ''.ljust(col_width[x]+2-len(c))+c
                    if x==ncol-1:
                        fwrows += '\n'

        f.write(fwrows)

def read_ascii(file, ascext=ascext, dtype=str, start_header=-1):
    '''
    Read ASCII file
    Supported format: commented_header
    See also: astropy.io.ascii.read, numpy.genfromtxt

    ------ INPUT ------
    file                input ASCII filename
    ascext              ASCII file suffix (Default: '.txt')
    dtype               data type (Default: 'str')
    start_header        line number of the header line (Default: -1)
    ------ OUTPUT ------
    dset                dataset (dict)
    '''
    with open(file+ascext, 'r') as f:
        ## f.read() -> str | f.readlines() -> list
        header = None
        darr = []
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            # print(line)
            if line[0]=='#':
                ## By default, header is taken from the last commented line before data
                if isinstance(start_header, int)==True:
                    if start_header==-1:
                        header = line.split()[1:]
                    elif start_header>0 and i==start_header-1:
                        header = line.split()[1:]
            else:
                line = list(map(dtype, line.split()))
                row = []
                for cell in line:
                    row.append(cell)
                darr.append(row)

    darr = np.array(darr)

    ## Number of columns/rows
    nrow, ncol = darr.shape

    ## Default header = ['col1','col2',...]
    header_default = []
    for x in range(ncol):
        header_default.append('col'+str(x+1))

    ## Dict out
    dset = {}
    for x, h in enumerate(header_default): 
        dset[h] = darr[:,x]
        if header is not None:
            dset[header[x]] = darr[:,x]

    return dset

def write_csv(file, header, dset, trans=False, append=False):
    '''
    Read fits file
    See also: astropy.io.ascii.write

    ------ INPUT ------
    file                output csv filename
    header              data labels in list('label1', 'label2', ...)
    dset                dataset (1darray[nrow], 2darray[nrow,ncol] or list)
    trans               transpose dset (Default: False)
    append              True: if not overwrite (Default: False)
    ------ OUTPUT ------
    '''
    if append==True:
        mod = 'a'
    else:
        mod = 'w'

    with open(file+csvext, mod, newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)

        ## Header
        ##--------
        writer.writeheader()

        ## Dataset
        ##---------
        dset = np.array(dset)
        ## Transpose dset
        if trans==True:
            Ny, Nx = dset.shape
            darr = []
            for x in range(Nx):
                for y in range(Ny):
                    darr.append(dset[y,x])
            darr = np.array(darr).reshape((Nx,Ny))
        else:
            darr = dset
        ## Number of columns/rows
        Ndim = darr.ndim
        if Ndim==1:
            ncol = 1
            nrow = darr.shape[0]
        else:
            nrow, ncol = darr.shape
        ## Write
        for y in range(nrow):
            ## Init dict
            rows = {h: [] for h in header}
            ## Write dict
            for x in range(ncol):
                rows[header[x]] = darr[y,x]
            ## Write csv row
            writer.writerow(rows)
            
def read_csv(file, *header):
    '''
    Read csv file

    ------ INPUT ------
    file                input csv filename
    header              labels of data to read
    ------ OUTPUT ------
    dset                dataset (dict)
    '''
    with open(file+csvext, 'r', newline='') as csvfile:
        dset = {}
        for h in header:
            csvfile.seek(0) # reset pointer
            reader = csv.DictReader(csvfile) # pointer at end
            col = []
            for row in reader:
                col.append(row[h])
            data = np.array(col)
            dset[h] = data

    return dset

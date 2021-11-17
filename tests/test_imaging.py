#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, logging, warnings
# logging.disable(sys.maxsize)
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", message="Skipping SYSTEM_VARIABLE record")

## Set dir
testdir = os.path.dirname(os.path.abspath(__file__))
datdir = testdir+'/lib/'
outdir = testdir+'/out/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

## Local
sys.path.insert(0, testdir+'/..') ## rapyuta path
from rapyuta.inout import fitsext, read_fits, write_fits
from rapyuta.astrom import fixwcs
from rapyuta.imaging import ( improve, Jy_per_pix_to_MJy_per_sr,
                              islice, icrop, irebin, igroupixel, ismooth,
                              imontage, iswarp, iconvolve, cupid,
                              wclean, interfill, hextract, hswarp,
                              concatenate )

print('\n TEST improve ')
print('--------------')
imp = improve(datdir+'M83', verbose=True)
print('raw: ', imp.im[0,:,3])
imp.rand_norm(datdir+'M83_unc')
print('add sym rand: ', imp.im[0,:,3])
imp.rand_splitnorm([datdir+'M83_unc', datdir+'M83_unc'])
print('add split rand', imp.im[0,:,3])
imp.rand_pointing(sigma=0.2, fill='near') # Spitzer IRS  pointing accuracy
print('add pointing rand: ', imp.im[0,:,3])
imp.crop(outdir+'M83_imp_crop', sizval=(0.005,0.009), cenpix=(7,8))
print('See out/M83_imp_crop.fits [Done]')

## IRS pointing error estimates
# Nmc = 20
# mcimage = []
# for j in range(Nmc+1):
#     imp = improve(datdir+'M83', verbose=True)
#     if j==0:
#         write_fits(outdir+'irs_pt',
#                    imp.hdr, imp.im, imp.wvl)
#     else:
#         imp.rand_norm(datdir+'M83_unc')
#         # imp.rand_pointing(sigma=0.2)
#         mcimage.append(imp.im)
# if Nmc>1:
#     mcimage = np.array(mcimage)
#     unc = np.nanstd(mcimage,axis=0)
#     write_fits(outdir+'irs_pt_unc',
#                imp.hdr, unc, imp.wvl)

imp1 = improve(datdir+'M83', verbose=True)
## Non-square pixel is not compatible with DS9
imp1.rebin(pixscale=(7,3), total=False, extrapol=True,
           filOUT=outdir+'M83_imp_rebin')
print('See out/M83_imp_rebin.fits [Done]')

imp2 = improve(datdir+'M83')
plt.figure()
plt.plot(imp2.wvl, imp2.im[:,0,0], c='k', label='raw')
imp2.smooth(6, wstart=None)
plt.plot(imp2.wvl, imp2.im[:,0,0], c='g', alpha=.5, label='smoothed')
plt.savefig(outdir+'M83_imp_smooth')
plt.close()
print('See out/M83_imp_smooth.fits [Done]')

imp3 = improve(datdir+'M82')
plt.figure()
data = read_fits(datdir+'M82').data[:,0,2]
wvl = read_fits(datdir+'M82').wave
unc = read_fits(datdir+'M82_unc').data[:,0,2]
plt.errorbar(wvl, data, unc, c='k', ecolor='r', label='raw')
data = imp3.artifact(filUNC=datdir+'M82_unc',
                     wmin=2.5, wmax=5, lim_unc=1.e2, cmin=6, fltr_pn='n')[:,0,2]
plt.errorbar(wvl, data, c='g', alpha=.5, label='Artifacts removed')
plt.legend(loc='best')
plt.xlim(2,6)
plt.ylim(-10,20)
plt.savefig(outdir+'M82_SL_imp_artifact')
plt.close()
print('See out/M82_SL_imp_artifact [Done]')

print('\n TEST Jy_per_pix_to_MJy_per_sr ')
print('-------------------------------')
Jy_per_pix_to_MJy_per_sr(datdir+'M82_IRAC1_DP', 
                         filOUT=outdir+'M82_IRAC1_unit')
print('See out/M82_IRAC1_unit.fits [Done]')

print('\n TEST islice ')
print('-------------')
slc = islice(datdir+'M83', filSL=outdir+'M83_inv_sqrt',
             filUNC=datdir+'M83_unc', slicetype='inv_sq', postfix='')
# if input('Clean slices (y/n): ')=='y':
slc.clean()

print('\n TEST icrop ')
print('------------')
crp = icrop(datdir+'M83', filOUT=outdir+'M83_icrop',
            sizpix=(3,6), cenval=(204.2529675, -29.8656962),
            filUNC=[datdir+'M83_unc', datdir+'M83_unc'],
            dist='splitnorm', verbose=True)
print('See out/M83_icrop.fits [Done]')

print('\n TEST irebin ')
print('-------------')
rbn = irebin(datdir+'M83', filOUT=outdir+'M83_irebin', verbose=True,
             pixscale=1, total=False, extrapol=False)
print('See out/M83_irebin.fits [Done]')

print('\n TEST igroupixel ')
print('-----------------')
grp = igroupixel(datdir+'M83', xscale=3, yscale=4,
                 filOUT=outdir+'M83_igroupixel')
print('See out/M83_igroupixel.fits [Done]')

print('\n TEST ismooth ')
print('--------------')
newave = read_fits(datdir+'M83').wave[:80] - 0.5
smt = ismooth(datdir+'M83', smooth=6, wgrid=newave, wstart=None,
              filOUT=outdir+'M83_ismooth')
print('See out/M83_ismooth.fits [Done]')

print('\n TEST imontage ')
print('---------------')
ds = read_fits(datdir+'M82_04_SL1')
hdr = fixwcs(datdir+'M82_04_SL1'+fitsext).header
write_fits(outdir+'sip', hdr, ds.data[0])

mtg = imontage('exact', tmpdir=outdir+'mtg/')
hdr_ref = fixwcs(datdir+'M82_template'+fitsext).header
mtg.reproject_mc(outdir+'sip', refheader=hdr_ref, filOUT=outdir+'sip_mtgrep')
hdr_nosip = hdr.copy()
for kw in hdr.keys():
    if ('A_' in kw) and (not 'PA' in kw) and (not 'RA' in kw):
        # print(kw)
        del hdr_nosip[kw]
    if 'B_' in kw:
        # print(kw)
        del hdr_nosip[kw]
    if 'AP_' in kw:
        # print(kw)
        del hdr_nosip[kw]
    if 'BP_' in kw:
        # print(kw)
        del hdr_nosip[kw]
write_fits(outdir+'nosip', hdr_nosip, ds.data[0])
mtg.reproject_mc(outdir+'nosip', refheader=hdr_ref, filOUT=outdir+'nosip_mtgrep')
print('SIP kw can be removed since they are encoded in CD matrix! [Done]')
mtg.coadd((datdir+'M82_04_SL1',datdir+'M82_06N_SL1'),
          refheader=hdr_ref, dist='norm', sig_pt=.2, fill_pt='near', Nmc=2,
          filOUT=outdir+'M82_SL1_mtgcoadd')
print('Coadd M82_04_SL1 & M82_06N_SL1 to M82_09_L86 [Done]')

print('\n TEST iswarp ')
print('-------------')
hdr_ref = fixwcs(datdir+'M82_09_L86'+fitsext).header
swp = iswarp(refheader=hdr_ref, tmpdir=outdir+'swp/')
# swp.combine(datdir+'M82_08_L86', dist='norm',
#             filOUT=outdir+'M82_08_L86_rep', tmpdir=outdir+'swp/')
swp.combine_mc(datdir+'M82_08_L86', dist='norm', sig_pt=.2, Nmc=2,
               filOUT=outdir+'M82_08_L86_rep', tmpdir=outdir+'swp/')
print('Reproject M82_08_L86 to M82_09_L86 [Done]\n')
swp_cube = iswarp((datdir+'M82_09_L86', datdir+'M82_04_SL1'),
                  center='9:55:51,69:40:45', pixscale=6.,
                  tmpdir=outdir+'swp_cube/')
swp_cube.combine(datdir+'M82_04_SL1', dist='norm',
                 filOUT=outdir+'M82_04_SL1_rep')
print('Reproject M82_04_SL1 (pixscale=6" recentered to M82 center) [Done]')
swp_coadd = iswarp((datdir+'M82_09_L86', datdir+'M82_04_SL1'),
                 refheader=hdr_ref, tmpdir=outdir+'swp_cube/')
swp_coadd.combine_mc((datdir+'M82_04_SL1',datdir+'M82_06N_SL1'),
                   dist='norm', sig_pt=.2, Nmc=2,
                   keepedge=True, cropedge=True, filOUT=outdir+'M82_SL1_swp_coadd')
print('Coadd M82_04_SL1 & M82_06N_SL1 to M82_09_L86 [Done]')

print('\n TEST iconvolve ')
print('----------------')
## See also idl/conv_prog.pro & idl/convolve_image.pro
convdir = outdir+'conv/' # see also idl/convolve_image.pro
if not os.path.exists(convdir):
    os.makedirs(convdir)
path_ker = datdir
path_idl = testdir+'/../idl/'
csv_ker = outdir+'kernelist' # see also idl/conv_prog.pro

irs_ker = []
psf = [2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.]
psf_ref = 'Gauss_06.0'
for p in psf:
    irs_ker.append(path_ker+'Kernel_HiRes_Gauss_0'+str(p)+'_to_'+psf_ref)
conv_cube = iconvolve(datdir+'M82_04_SL1',
                      psf=psf, kfile=irs_ker, klist=csv_ker,
                      dist='norm', sig_pt=.2,
                      convdir=convdir, filOUT=outdir+'M82_04_SL1'+'_conv')
conv_cube.do_conv(path_idl, verbose=False)
print('Convolve M82_04_SL1 [Done]')

irac_ker = path_ker+'Kernel_HiRes_IRAC_8.0_to_Gauss_06.0'
conv = iconvolve(datdir+'M82_IRAC4',
                 kfile=irac_ker, klist=csv_ker,
                 filOUT=outdir+'M82_IRAC4'+'_conv')
conv.do_conv(path_idl)
print('Convolve M82_IRAC4 [Done]')

print('\n TEST cupid ')
print('------------')
parobs = [
    ['3390001.1','NG','F011100297_N002','Ns',24,2,5],
    ['3390001.1','NG','F011100297_N002','Nh',28,7,3],
]
buildir = outdir+'/cubuild/'
if not os.path.exists(buildir):
    os.makedirs(buildir)
fits_irc = []
for obs in parobs:
    fits_irc.append(obs[0]+'_'+obs[3])
## MC add pointing unc
Nmc = 2
for i in range(len(parobs)):
    for j in range(Nmc+1):
        if j==0:
            verbose = False
        else:
            verbose = False
        
        cup = cupid(datdir, obsid=parobs[i][0], slit=parobs[i][3],
                    spec=parobs[i][1], imref=parobs[i][2], verbose=verbose)
        if j==0:
            cup.spec_build(buildir+fits_irc[i],
                           filRAW=buildir+'raw', tmpdir=buildir, fiLOG=buildir+'build_history',
                           Nx=parobs[i][6], Ny=parobs[i][4], Nsub=parobs[i][5],
                           pixscale=1, wmax=4.3, supix=True)
        else:
            cup.spec_build(buildir+fits_irc[i]+'_'+str(j), tmpdir=buildir,
                           Nx=parobs[i][6], Ny=parobs[i][4], Nsub=parobs[i][5],
                           pixscale=1, wmax=4.3, supix=True,
                           dist='splitnorm', sig_pt=3, fill_pt='med', swarp=False)
    mcimage = []
    for j in range(Nmc+1):
        if j==0:
            hdr = read_fits(buildir+fits_irc[i]).header
            wvl = read_fits(buildir+fits_irc[i]).wave
            data = read_fits(buildir+fits_irc[i]).data
            write_fits(outdir+fits_irc[i], hdr, data, wvl)
        else:
            mcimage.append(read_fits(buildir+fits_irc[i]+'_'+str(j)).data)
    if Nmc>1:
        mcimage = np.array(mcimage)
        unc = np.nanstd(mcimage, axis=0)
        write_fits(outdir+fits_irc[i]+'_unc', hdr, unc, wvl)

print('\n TEST wclean ')
print('-------------')

print('\n TEST interfill ')
print('----------------')
data = read_fits(datdir+'IC10_SL2').data[0]
hdr = fixwcs(datdir+'IC10_SL2'+fitsext).header
newdata = interfill(data, axis=0)
write_fits(outdir+'IC10_fillgap', hdr, newdata)
print('interfill test [Done]')

print('\n TEST hextract ')
print('---------------')
hextract(datdir+'M82_08_L86', outdir+'M82_hextract', 20, 40, 32, 42)
print('hextract test [Done]')

print('\n TEST hswarp ')
print('-------------')
print('Reproject M82_08_L86 to M82_09_L86')
old = read_fits(datdir+'M82_08_L86')
oldimage = old.data
oldheader = old.header
refheader = read_fits(datdir+'M82_09_L86').header
hswp = hswarp(oldimage, oldheader, refheader, keepedge=True,
              tmpdir=outdir+'hswp/', verbose=False)
# print('hswarp image: ', hswp.data)
# print('hswarp image header: ', hswp.header)
print('See hswp/coadd.fits [Done]')

print('\n TEST concatenate ')
print('------------------')
flist = [datdir+'M82_SL1', datdir+'M82_SL2']
concat = concatenate(flist, filOUT=outdir+'M82_SL_concat',
                     keepfrag=False, cropedge=True)
wvl = concat.wave
data = concat.data[:,0,0]
uncl = [f+'_unc' for f in flist]
unc = concatenate(uncl, filOUT=outdir+'M82_SL_concat_unc').data[:,0,0]
plt.figure()
plt.errorbar(wvl, data, unc, c='k', ecolor='r')
plt.savefig(outdir+'M82_SL_concat')
plt.close()
print('See out/M82_SL_concat.png [Done]')

if input('Clean tmp files (y/n): ')=='y':
    mtg.clean()
    swp.clean()
    swp_cube.clean()
    conv_cube.clean()
    conv.clean()
    conv.clean(outdir+'hswp/')


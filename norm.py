# ************************************************
# SCRIPT TO NORMALIZE THE FINAL SPECTRUM
#
# Roberta Razera - July / 2021
# ************************************************
import matplotlib
matplotlib.use('Agg')
import astropy.io.fits as pyfits
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
import math
import os
from matplotlib import colors
matplotlib.rc('font',**{'family':'sans-serif'})
matplotlib.rc('text', usetex=True)
matplotlib.rc('xtick', labelsize=9) 
matplotlib.rc('ytick', labelsize=9) 
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib.colors import LogNorm
import pylab
import pandas as pd
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import matplotlib.artist as artists
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from datetime import datetime
import utils
import pandas as pd

save_name = 'apStar-r12-2M17265563-2813558.dat'
mypath = '/home/robertarazera/Documents/2M17265563-2813558/PLOTS/'
obj = '2M17265563-2813558'

# Normalization of the synthetic spectra
print('----------------------------------')
print('--- Starting the normalization ---')
print('----------------------------------')

path1 = '/home/robertarazera/Downloads/ESPECTROS_APOGEE/'
syn_spec = path1 + save_name
print ("ESPECTRO:")
print syn_spec
syn = np.loadtxt(syn_spec)
syn_wave = syn[:,1] #dont forget to see if this columns are as expected: wave, flux
syn_flux = syn[:,2]

frac = 0.02
it = 35
delta = 0.5

x = utils.fit_cont(syn_wave, syn_flux, frac = frac, it = it, delta = delta)

lowess_x = list(zip(*x))[0]
lowess_y = list(zip(*x))[1]

old_flux = []
for i in lowess_x:
    flux_f = syn_flux[syn_wave == i]
    old_flux = np.append(old_flux,flux_f)

flux_norm = old_flux/lowess_y

data = np.c_[lowess_x,flux_norm]

# - Saving the normalized spectrum
save_norm = '2M17265563-2813558_norm.dat'
add='0.5D_'
np.savetxt(mypath+add+save_norm,data,fmt='%.4f %.5f',delimiter=' ',header='Wave Flux_norm')
saved = mypath+add+save_norm
print('Normalized spectrum saved as: %s' % (saved))

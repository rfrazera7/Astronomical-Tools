# ************************************************
# SCRIPT TO READ FITS, DERIVE AIR WAVELENGTH, FLUX
# AND NORMALIZE THE FINAL SPECTRUM
#
# Roberta Razera - May / 2021
# ************************************************
from astropy.io import fits
import numpy as np
import utils
import pandas as pd

elem1 = "apStar-r12-"
column_names = ["name", "Na", "Mg", "Al", "Si", "P", "S", "K", "Ca", "Ti", "V", "Cr", "Mn", "Co", "Ni", "Cu"]
df = pd.read_csv("/home/robertarazera/Downloads/abundancias_parametros_apogee.csv", names=column_names)
name = df.name.to_list()

for i in range(len(name)):
	if i > 0:
		new_name = name[i]
		star = elem1 + new_name + ".fits"
		path = ('/home/robertarazera/Downloads/ESPECTROS_APOGEE/')
		spec = fits.open('/home/robertarazera/Downloads/ESPECTROS_APOGEE/' + star)
		
		# Defining function to convert to wave_air and save data
		def write_dat(flux = None):
			print('-----------------------------------------------------------------')
			print("------------------ NOW LET'S CHECK THE HEADER: ------------------")
			print('-----------------------------------------------------------------')
			print(repr(spec[1].header))
			print('-----------------------------------------------------------------')
			print("Check to see if the names in the header are correct!")
			crval = spec[1].header['CRVAL1']
			cdelt = spec[1].header['CDELT1']
			naxis = spec[1].header['NAXIS1']
			print ('NAXIS1 = '), naxis, ('='), spec[1].header['NAXIS1'] # Number of points
			print ('CDELT1 = '), cdelt, ('='), spec[1].header['CDELT1'] # Delta lambda
			print ('CRVAL1 = '), crval, ('='), spec[1].header['CRVAL1'] # First lambda, where the spectrum starts
			print('-----------------------------------------------------------------')
			fac = crval + (cdelt * naxis)
			print("Now let's get wavelength in A (linear), then transform it to AIR:")
			print('-----------------------------------------------------------------')
			wave = 10**(np.linspace(crval, fac, num=naxis))
			AIR = wave/(1.0+5.792105E-2/(238.0185E0-(1.E4/wave)**2)+1.67917E-3/(57.362E0-(1.E4/wave)**2)) 
			print ("AIR= "), AIR
			print('-----------------------------------------------------------------')
			data = np.c_[wave, AIR, flux]
			print('------------------------ Saving data... -------------------------')
			print('-----------------------------------------------------------------')
			print data
			print('##############################################################################################################################################')
			saved = elem1 + new_name + ".dat"
			np.savetxt(path+saved, data, delimiter=' ', fmt='%.5f %.5f %.4f', header='wave_vac, wave_air, flux')
			print ("-> -> -> -> -> -> ->"), ('Data saved as: %s' % (saved)), ("at: %s" % (path)), ("<- <- <- <- <- <- <-")
			print('##############################################################################################################################################')
			# Normalization of the synthetic spectra
			print('-----------------------------------------------------------------')
			print('------------------- Starting the normalization ------------------')
			print('-----------------------------------------------------------------')
			syn_spec = path + saved
			syn = np.loadtxt(syn_spec)
			syn_wave = syn[:,1]
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

			save_dir = path + 'NORMALIZADOS/'
			save_norm = elem1 + new_name + '-norm.dat'
			np.savetxt(save_dir+save_norm, data, fmt='%.4f %.5f',delimiter=' ',header='Wave Flux_norm')
			saved2 = save_dir+save_norm
			print('Normalized spectrum saved as: %s' % (saved2))
			print(' (     (        )  (    (        )      (      ')
			print(' )\ )  )\ )  ( /(  )\ ) )\ )  ( /(      )\ )   ')
			print('(()/( (()/(  )\())(()/((()/(  )\()) (  (()/(   ')
			print(' /(_)) /(_))((_)\  /(_))/(_))((_)\  )\  /(_))  ')
			print('(_))_|(_))   _((_)(_)) (_))   _((_)((_)(_))_   ')
			print('| |_  |_ _| | \| ||_ _|/ __| | || || __||   \  ')
			print('| __|  | |  | .` | | | \__ \ | __ || _| | |) | ')
			print('|_|   |___| |_|\_||___||___/ |_||_||___||___/  ')
		print('									')
		print('									')
		print('									')
		print('----------------------------------------------------------------')
		print ("We are going to work with |"), star, ("|")
		print('----------------------------------------------------------------')
		spec.info()
		print('----------------------------------------------------------------')
		if len(spec[1].data.shape) == 1:
			flux = spec[1].data
			print('									')
			print ('This is your flux: '), flux
			print('									')			
			write_dat(flux=flux)
		else:
			flux = spec[1].data[0]	
			print('									')
			print ('This is your flux: '), flux
			print('									')
			write_dat(flux=flux)






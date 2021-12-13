# Minimal demonstration of using AstroPy specutils to read a plot reduced
# spectra from Liverpool Telescope SPRAT level 2 reduced FITS files

from glob import glob
import os
import argparse
import numpy as np
import specutils as sp
import statistics as st

from astropy.io import fits
from astropy import units as u
from astropy.visualization import quantity_support
import astropy.wcs as fitswcs
from astropy.wcs.docstrings import ORIGIN

from spectres import spectres
from scipy.signal import find_peaks
from matplotlib import pyplot as plt

quantity_support()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a SPRAT spectrum')

    parser.add_argument('-d', dest='dir_str', help='Directory containing the reduced SPRAT "_2.fits" files to display. Directory must be at same level as sprat_splot.py script. Mandatory.', default = '2021zby')
    parser.add_argument('-t', dest='datime_str', help='Date of the observation within the night. Also stored within the name of the educed SPRAT "_2.fits" files. Required format is in YYYYMMDD. Mandatory.', default = '20211008')
    parser.add_argument('-x', dest='extn_str', choices=['2','SPEC_NONSS','3','SPEC_SS','4','NORMFLUX','5','FLUX'], default='5', help='Multi-extention to display (default: 5, FLUX). Mandatory.')
    parser.add_argument('-z', dest='redshift_str', default='none', help='Redshift of host used in plotting the rest wavelength axis (default: none). Optional.')
    parser.add_argument('-b', dest='binfactor_str', default='none', help='Bin factor used to rebin and smooth noisy spectra (default: none). Optional.')
    parser.add_argument('-v', dest='verbose', action='store_true', help='Turn on verbose mode (default: Off).')

    args = parser.parse_args()

    # Parse the string names of the extensions into extension numbers
    if args.extn_str in ['2','SPEC_NONSS']: 
      extension = 2
      unitName = "adu"
    elif args.extn_str in ['3','SPEC_SS']: 
      extension = 3
      unitName = "adu"
    elif args.extn_str in ['4','NORMFLUX']: 
      # Relative flux normalized to 5500A is dimensionless
      extension = 4
      unitName = "Normalised"
    elif args.extn_str in ['5','FLUX']: 
      extension = 5
      unitName = "erg/s/cm2/A"
    
if args.verbose:
  print (args)

# Convert the host redshift into a number if one is provided
try:
  host_z = float(args.redshift_str)
except:
  host_z = 'none'

# Convert the bin factor into a number if one is provided
try:
  bf = float(args.binfactor_str)
except:
  bf = 'none'

# Read in and perfrom a naive median stack of data to get specdata
obj = args.dir_str
date = int(args.datime_str)
data = []
os.chdir(f'{obj}')
for fn in glob(f"v_e_{date}*_2.fits"):
    f=fits.open(fn)
    data.append(f[extension].data)
data = np.array(data)
specdata = np.median(data, axis=0)[0]

# Spectra are stored as 2D NAXIS1 x 1 arrays. Read and convert to a 1D NAXIS vector. 
specheader = f[extension].header
wcsheader = dict()
for key in ("CDELT1", "CRVAL1", "CUNIT1", "CTYPE1", "CRPIX1"):
    wcsheader[key] = specheader[key]

# Grab some exposure information from the header
telescope =   '# TELESCOPE = ' + specheader["TELESCOP"]
instrument =  '# INSTRUMENT = ' + specheader["INSTRUME"]
user =        '# OWNER = ' + specheader["USERID"]
proposal =    '# PROPOSAL = ' + specheader["PROPID"]
object =      '# OBJECT = ' + specheader["CAT-NAME"]
ra =          '# RA = ' + specheader["CAT-RA"]
dec =         '# DEC = ' + specheader["CAT-DEC"]
dateObs =     '# DATE-OBS = ' + specheader["DATE-OBS"]
mjdObs =      '# MJD-OBS = ' + str(specheader["MJD"])
exptime =     '# EXPTIME-OBS = ' + str(specheader["EXPTIME"]) + 's'
airmass =     '# AIRMASS-OBS = ' + str(specheader["AIRMASS"])
seeing =      '# SEEING-OBS = ' + str(specheader["ESTSEE"])
colData =     '# wl fl'
colType =     f'# [A] [{unitName}]'

# For debug purposes
#for item in specheader:
#  print(item +' : '+str(specheader[item]))

flux = specdata * u.Unit("adu")
# create WCS Wavelength Calibration from fits header
my_wcs = fitswcs.WCS(wcsheader)

# Make Spectrum1D object 
sp1d = sp.Spectrum1D(flux=flux,wcs=my_wcs)

# Extract the wavelength and flux values as separate xy axes lists
wl_list = sp1d.spectral_axis.value
flux_list = sp1d.flux.value

# Find the median flux around 6000A region and produce a normalised flux list using this median
median_wl_list_idx = np.argsort(wl_list)[len(wl_list)//2]
mid_range_fluxes = sorted(flux_list[median_wl_list_idx-25:median_wl_list_idx+26])
median_mid_flux = st.median(mid_range_fluxes)
flux_list_norm = flux_list / median_mid_flux
# Idenitfy the position of any outlier peaks and troughs that are significantly bigger than their neighbours and remove them from the spectrum
outlier_peaks_list_ids, _ = find_peaks(flux_list_norm, threshold=0.5)
outlier_troughs_list_ids, _ = find_peaks(-flux_list_norm, threshold=0.5)
all_outlier_list_ids = sorted([*outlier_peaks_list_ids, *outlier_troughs_list_ids])
wl_list = np.delete(wl_list, all_outlier_list_ids)
flux_list = np.delete(flux_list, all_outlier_list_ids)

# Perform any aspectral smoothing if a binning factor has been provided
if bf != 'none':
  wav = wl_list
  flux = flux_list
  wav_bin = np.linspace(wav[0], wav[-1], int(len(wav)/bf))
  flux_bin = spectres(wav_bin, wav, flux)
  wav_bin = wav_bin
else:
  wav_bin = wl_list
  flux_bin = flux_list


ascii_spec = np.column_stack((wav_bin, flux_bin))
np.savetxt(f"{obj}_{date}_SPRAT.txt",ascii_spec, fmt="%f %g")

# Create text file of spectrum
with open(f"{obj}_{date}_SPRAT.txt",'r') as f:
  with open('tempfile.txt','w') as f2: 
    f2.writelines([
        telescope + '\n', 
        instrument + '\n', 
        user + '\n', 
        proposal + '\n', 
        object + '\n', 
        ra + '\n', 
        dec + '\n', 
        dateObs + '\n', 
        mjdObs + '\n', 
        exptime + '\n', 
        airmass + '\n', 
        seeing + '\n', 
        colData + '\n', 
        colType + '\n'])
    f2.write(f.read())
os.rename('tempfile.txt',f"{obj}_{date}_SPRAT.txt")

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.plot(wav_bin,flux_bin)
ax1.set_xlabel(f"Observed Wavelength [A]")
ax1.set_ylabel(f"Flux [{unitName}]")

if host_z != 'none':

  ax2 = ax1.twiny()
  ax2.set_xticks(ax1.get_xticks())
  ax2.set_xbound(ax1.get_xbound())
  ax2.set_xticklabels([ int(round(x / (host_z + 1),0)) for x in ax1.get_xticks()])
  
  ax2.set_xlabel(f"Rest Wavelength [A]")

plt.savefig(f"{obj}_{date}_SPRAT.png")
plt.show()

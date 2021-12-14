# Minimal demonstration of using AstroPy specutils to read a plot reduced
# spectra from Liverpool Telescope SPRAT level 2 reduced FITS files

import os
import argparse
import pandas as pd

from astrorapid.classify import Classify

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run AstroRAPID on an object from the ZTF forced Photometry')

    parser.add_argument('-n', dest='filename_str', help='Name of the ZTF forced photometry .txt file to fit with RAPID. File must be inside candidates direcotry. Mandatory.', default = 'ZTF21aaqytjr')
    parser.add_argument('-z', dest='redshift_str', help='Spectroscopic or photometric redshift of the object or host. Optional.', default = '0.003319')
    parser.add_argument('-ra', dest='ra_deg_str', help='Right Ascension of the object in degrees. Mandatory', default = '186.423654	')
    parser.add_argument('-dec', dest='dec_deg_str', help='Declination of the object in degrees. Mandatory', default = '7.228395')
    parser.add_argument('-mwebv', dest='ebv_str', help='Total galactic visual extinction (mwebv) in the line-of-sight of the object. Mandatory', default = '0.2851')
    parser.add_argument('-v', dest='verbose', action='store_true', help='Turn on verbose mode (default: Off).')

    args = parser.parse_args()

    # check if file exists
    objid = args.filename_str
    filepath = os.path.dirname(os.path.abspath(__file__)) + f'/candidates/{objid}.txt'
    if os.path.isfile(filepath) : 
      print(f'{objid}.txt file found!')
    else:
      print(f'{objid}.txt file not found!')
      quit()
    
if args.verbose:
  print(args)

# Convert the ra, dec, dust reddening and redshift into a number if one is provided
ra = float(args.ra_deg_str)
dec = float(args.dec_deg_str)
mwebv = float(args.ebv_str)
try:
  redshift = float(args.redshift_str)
except:
  redshift = 0.0

# Read in forced photometry file and drop comment lines
data = []
with open(filepath) as fp:
  for line in fp:
      if not line.lstrip().startswith('#'):
        line=line.strip()
        line=line.strip("\n")
        data.append(line)

# turn photometry into a semi dataframe
headers = data[0].replace(" ", "").split(",")
data.pop(0)
def floatify(val):
    try:
        return float(val)
    except Exception as e:
        return val
semi_df = [[floatify(x) for x in e.split()] for e in data]

# convert semi dataframe into a full pandas dataframe
full_df = pd.DataFrame.from_records(semi_df)
full_df.columns = headers

# change working directory to the outputs folder
outpath = os.getcwd() + '/outputs/'
os.chdir(outpath)

# extract from the dataframe and create the mjd, flux, fluxerr and passband lists that rapid needs for making the lightcurve
mjd = list(full_df["jd"] - 2400000.5)
flux = list(full_df["forcediffimflux"])
fluxerr = list(full_df["forcediffimfluxunc"])
filters = list(full_df["filter"])
passband = [f.replace("ZTF_", "") for f in filters]

# use the flux and fluxerr list to create the list of photo flags (0 - non detection, 4096 - normal detection, 6144 - trigger detection)
photflag = []
trigger = 0
for idx, val in enumerate(flux):
  # if flux is not a 3-sig detection then flag as a non detection
  # else flag as a normal detection
  if val < 3*fluxerr[idx]:
    photflag.append(0)
  else:
    # but if flux is the first 5-sig detection then flag as the trigger detection
    if val >= 5*fluxerr[idx] and trigger == 0:
      trigger = 1
      photflag.append(6144)
    else:
      photflag.append(4096)

# Light curve information for RAPID
light_curve_info = (mjd, flux, fluxerr, passband, photflag, ra, dec, objid, redshift, mwebv)

# Classification
light_curve_list = [light_curve_info, ]  # Add more light curves to be classified to this list

classification = Classify(known_redshift=True)
predictions = classification.get_predictions(light_curve_list)
#print(predictions)

# Plot classifications vs time
classification.plot_light_curves_and_classifications()
classification.plot_classification_animation()
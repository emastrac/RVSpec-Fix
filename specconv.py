5#! /usr/bin/env python3

###############################################
#                                             #
#   ***   SPECTRUM CONVERTER TEMPLATE   ***   #
#    ***          for RAVESPAN         ***    #
#     ***      by Bogumil Pilecki     ***     #
#      ***       configured by       ***      #
#       ***    Everett Mastracci    ***       #
###############################################

#################################################################
# This file is part of the RaveSpan program
# Copyright (C) 2024 Bogumil Pilecki
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################
#
# RaveSpan program includes the following files:
# * ravespan.py  - main RaveSpan GUI
# * librvc.py    - radial velocity curve plot
# * libspec.py   - spectrum related operations
# * libdata.py   - information about data points
# * libanal.py   - analysis window, velocity measurements
# * libdial.py   - various additional dialog windows
# * libcommpn.py - a bunch of small, useful functions 
# * bin/rvspec       - executable python script to run RaveSpan
# * work/specconv.py - a template program for spectra conversion to RaveSpan format
#
#################################################################


from astropy.io import fits as pyfits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
import sys
import numpy as np
import os
from scipy import interpolate


#####################################################
def get_fileholder(fname):
    try:
        pyf = pyfits.open(fname)
        return pyf
    except:
        print("IO error: no file or file cannot be read.")
        sys.exit(0)

args = sys.argv[1:]

if len(args) < 2:
    good_input = False
    while not good_input:
        print("Please enter the file mode: (object/template)")
        mode = input()
        if not (mode == "object" or mode == "template"):
            print("Bad input. Try again:")
        else:
            good_input = True
    print("Please enter the FITS file name(s) to process:")
    args = input().split(' ')
else:
    mode = args[0]
    args = args[1:]

# ITERATE THROUGH FILES GIVEN AS ARGUMENTS:
# e.g.
# spec_convert objectX_2017_02_12.fits objectX_2017_02_13.fits objectX_2017_02_14.fits
# 
for ia,fname in enumerate(args):
    print("%d/%d"%(ia+1, len(args)))
    print(" * file:", fname)
    
    # OPEN FILE
    pyf = get_fileholder(fname)

    # GET SPECTRUM
    spec = pyf[0].data
    NaN_count = 0
    for i in range(0, len(spec)):
        if np.isnan(spec[i]):
            spec[i] = np.median(spec[~np.isnan(spec)])
            NaN_count = NaN_count + 1
    print("WARNING! Replaced " + str(NaN_count) + " NaN values.")
    if isinstance(spec[0], np.ndarray):
        lambda_0 = spec[3][0]*10
        unique_wavelengths, unique_flux_indices = np.unique(spec[3], return_index=True)
        wavelengths_resamp = np.linspace(spec[3][0], spec[3][spec[3].size - 1], spec[3].size)
        diff = np.diff(wavelengths_resamp)
        spec = interpolate.PchipInterpolator(unique_wavelengths, spec[4][unique_flux_indices])(wavelengths_resamp)
        
    if mode == "template":
        spec = np.divide(spec, np.median(spec))

    # GET HEADER
    pyf.verify("fix")
    pheader = pyf[0].header
    # GET ALL THE NECESSARY DATA FROM THE HEADER
    # ...
    # JD/HJD/BJD = ...      // time of the observation 
    # LAMBDA_O = ...        // wavelength of the beginning of the spectrum
    # RESOLUTION = ...      // spectrum resolution
    # BARY_CORR = ...       // barycentric correction (in km/s) 
    # INSTRUMENT = ...      // instrument name
    # ...
    lambda_0 = 0
    # resolution = np.diff(wavelengths_resamp)[0]*10
    if "CDELT1" in pheader.keys():
        resolution = float(pheader["CDELT1"])
    elif "DELTA_WL" in pheader.keys():
        resolution = float(pheader["DELTA_WL"])
    elif "SP-RESOL" in pheader.keys():
        resolution = float(pheader["SP-RESOL"])
    else:
        print(f"Resolution not found in header for file {fname}. Please enter the resolution in angstroms per pixel:")
        resolution = float(input())

    if "CRPIX1" in pheader.keys() and "CRVAL1" in pheader.keys():
        spec = spec[int(pheader["CRPIX1"]) - 1:]
        lambda_0 = pheader["CRVAL1"]
    elif "WAVELENG" in pheader.keys() and "REFPIXEL" in pheader.keys():
        lambda_0 = pheader["WAVELENG"] - (pheader["REFPIXEL"]*resolution)
    else:
        print(f"Wavelength data not found in header for file {fname}. Please enter the lowest wavelength in the spectrum in angstroms:")
        lambda_0 = float(input())

    if "OBJECT" in pheader.keys():
        object_name = pheader["OBJECT"]
    else:
        print(f"Object name not found in header for file {fname}. Please enter it:")
        object_name = input()

    if mode == "object":
        obs_time = 0
        barycorr = 0
        instrument = ""
        obs_date_keywords = ["JD", "DATE-OBS", "TIME-OBS", "TIME-END", "DATE-END"]
        date_found = False
        for key in obs_date_keywords:
            if key in pheader.keys() and not date_found:
                if key == "JD": 
                    obs_time = pheader["JD"]
                else:
                    obs_time = Time(pheader[key], format="fits", scale="utc").jd
                date_found = True
        if not date_found:
            print(f"Observation date not found in header for file {fname}. Please enter the date in Julian Date (JD) format:")
            obs_time = float(input())
        if "VHELIO" in pheader.keys():
            barycorr = float(pheader["VHELIO"])
        elif "LATITUDE" in pheader.keys() and "LONGITUD" in pheader.keys() and "RA" in pheader.keys() and "DEC" in pheader.keys():
            sc = SkyCoord(pheader["RA"], pheader["DEC"], unit=(u.hourangle, u.deg))
            try:
                obs_site = EarthLocation.from_geodetic(pheader["LATITUDE"], pheader["LONGITUD"])
            except:
                print("Error querying EarthLocation. Please enter the barycentric correction manually in km/s:")
                barycorr = float(input())
            else:
                barycorr = sc.radial_velocity_correction("barycentric", Time(obs_time, format = "jd"), obs_site).value/1000
        else:
            print(f"Barycentric correction data not found in header for file {fname}. Please enter the value in kilometers per second:")
            barycorr = float(input())

        if "INSTRUME" in pheader.keys():
            instrument = pheader["INSTRUME"]
        else:
            print(f"Instrument not found in header for file {fname}. Please enter the name of the instrument used for the observation:")
            instrument = input()

    # TRANSFORM THE SPECTRUM IF NECESSARY
    # RaveSpan only reads uniformly sampled spectra,
    # if your spectra have non-uniform sampling,
    # you have to resample it.
    # ...
    # spec = spectrum_conversion(spec, ...)
    # ...  


    # PREPARE THE FILE NAME
    # it has 5 segments:
    # * obs_time  -  time of the observation
    # * lambda_0 - wavelength of the beginning of the spectrum
    # * resolution - spectrum resolution
    # * barycorr - barycentric correction (in km/s)
    # * instrument - instrument name
    # e.g.
    # 7234.12334_4501.12345431_0.01555_-2.123_HARPS
    if mode == "object":
        filename = '%.5f_'%obs_time+str(lambda_0)+'_'+str(resolution)+'_'+'%.3f'%(barycorr)+'_'+instrument
    else:
        filename = object_name+'_'+str(lambda_0)+'_'+str(resolution)


    # SAVE SPECTRUM AS A BINARY FILE
    if mode == "object":
        path = "specdb" + '/' + object_name
    else:
        path = "templates"
    os.makedirs(path, exist_ok=True)
    np.array(spec, np.float32).tofile(path + '/' + filename)
    open("objects/" + object_name + ".obj", 'a').close()

    # THEN COPY IT TO THE CORRESPONDING DIRECTORY IN specdb/
    # e.g. /home/user/ravespan/specdb/object_name/7234.12334_4501.12345431_0.01555_-2.123_HARPS

    if mode == "object":
        print("   ---->", "objects/" + object_name + ".obj")
    print("   ---->", path + '/' + filename)

###############
### THE END ###
###############

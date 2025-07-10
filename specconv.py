#! /usr/bin/env python3

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
import sys
from numpy import *
import os


#####################################################
def get_fileholder(fname):
    try:
        pyf = pyfits.open(fname)
        return pyf
    except:
        print("IO error: no file or file cannot be read.")
        sys.exit(0)

args = sys.argv[1:]
args = ["dao.fits"]

if len(args) == 0:
    print("No filename given.")
    sys.exit()

# ITERATE THROUGH FILES GIVEN AS ARGUMENTS:
# e.g.
# spec_convert objectX_2017_02_12.fits objectX_2017_02_13.fits objectX_2017_02_14.fits
# 
for ia,fname in enumerate(args):
    print("%d/%d"%(ia+1, len(args)))
    print(" * file:", fname)
    
    # OPEN FILE
    pyf = get_fileholder(fname)


    # GET HEADER
    pheader = pyf[0].header
    # GET ALL THE NECESSARY DATA FROM THE HEADER
    # ...
    # JD/HJD/BJD = ...      // time of the observation 
    # LAMBDA_O = ...        // wavelength of the beginning of the spectrum
    # RESOLUTION = ...      // spectrum resolution
    # BARY_CORR = ...       // barycentric correction (in km/s) 
    # INSTRUMENT = ...      // instrument name
    # ...
    obs_time = 0
    lambda_0 = 0
    resolution = 0
    barycorr = 0
    instrument = ""
    obs_date_keywords = ["DATE-OBS", "TIME-OBS", "TIME-END", "DATE-END"]
    lambda_0_keywords = ["WAVELENG"]
    date_found = False
    lambda_0_found = False
    for key in obs_date_keywords:
        if key in pheader.keys():
            obs_time = Time(pheader[key], format="fits", scale="utc").jd
            date_found = True
    if not date_found:
        print(f"Observation date not found in header for file {fname}. Please enter the date in Julian Date (JD) format:")
        obs_time = float(input())

    if "WAVELENG" in pheader.keys() and "REFPIXEL" in pheader.keys() and "DELTA_WL" in pheader.keys():
        lambda_0 = pheader["WAVELENG"] - (pheader["REFPIXEL"]*pheader["DELTA_WL"])
    else:
        print(f"Wavelength data not found in header for file {fname}. Please enter the lowest wavelength in the spectrum in angstroms:")
        lambda_0 = float(input())

    if "DELTA_WL" in pheader.keys():
        resolution = pheader["DELTA_WL"]
    else:
        print(f"Resolution not found in header for file {fname}. Please enter the resolution in angstroms per pixel:")
        resolution = float(input())

    if "VHELIO" in pheader.keys():
        barycorr = pheader["VHELIO"]
    else:
        print(f"Barycentric correction data not found in header for file {fname}. Please enter the value in kilometers per second:")
        barycorr = float(input())

    if "INSTRUME" in pheader.keys():
        instrument = pheader["INSTRUME"]
    else:
        print(f"Instrument not found in header for file {fname}. Please enter the name of the instrument used for the observation:")
        instrument = input()

    # GET SPECTRUM
    spec = pyf[0].data
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
    filename = '%.5f_'%obs_time+str(lambda_0)+'_'+str(resolution)+'_'+'%.3f'%(barycorr)+'_'+instrument


    # SAVE SPECTRUM AS A BINARY FILE
    if '_' in fname:
        object_name = fname.split('_')[0]
        path = "specdb" + '/' + object_name + '/'
    else:
        path = "specdb" + '/' + fname.split(".")[0] + '/'
    os.makedirs(path, exist_ok=True)
    array(spec, float32).tofile(path + filename)

    # THEN COPY IT TO THE CORRESPONDING DIRECTORY IN specdb/
    # e.g. /home/user/ravespan/specdb/object_name/7234.12334_4501.12345431_0.01555_-2.123_HARPS

    print("   ---->", path + filename)



###############
### THE END ###
###############




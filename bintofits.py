# Binary file to FITS file converter

import numpy as np
from astropy.io import fits

filename = "TEMP_3800_0.02"
path = "templates/"
spec = np.fromfile(path + filename, dtype=np.float32)

hdu = fits.PrimaryHDU(data=spec)
hdr = hdu.header
keywords = filename.split('_')
hdr["JD"] = float(keywords[0])
hdr["WAVELENG"] = float(keywords[1]) + (spec.size*float(keywords[2]))/2
hdr["SP-RESOL"] = float(keywords[2])
hdr["VHELIO"] = float(keywords[3])
hdr["INSTRUME"] = keywords[4]

hdu.writeto(filename + ".fits", overwrite=True)
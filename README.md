# RVSpec-Fix
Library fixes for RaveSpan to allow the program to run.
You can find the download for RaveSpan [here](https://users.camk.edu.pl/pilecki/ravespan/?page=download).

## Installation
### Windows
1. If you have a version of Python 3.9 installed, skip to step 3. Otherwise, download [Python 3.9.13](https://www.python.org/downloads/release/python-3913/) (or any version of Python 3.9).
2. Install the version of Python you just downloaded and make note of the installation location.
3. Download RaveSpan from [Dr. Pilecki's website](https://users.camk.edu.pl/pilecki/ravespan/?page=download).
4. Optionally, download the template library slib from the same page.
5. Go to the [Code](https://github.com/emastrac/RVSpec-Fix) section of this repository, and download the files from the green Code dropdown menu.
6. Extract `ravespan_p3.tgz` and `RVSpec-Fix-main.zip` (it does not matter where).
7. If you would like to use my FITS file conversion script, replace `ravespan_p3/work/specconv.py` with the version in `RVSpec-Fix-main`.
8. Drag the rest of the `.py` files from `RVSpec-Fix-main` into `ravespan_p3/rvspeclib` and replace the existing library files.
9. Open a command prompt and `pip install` `matplotlib`, `numpy`, `scipy`, and `PyQt5`.
* If you had another version of Python installed already, make sure to do this with `\path\to\your\Python\Python39\Scripts\pip.exe install ...` instead.
10. Do `cd \path\to\your\ravespan_p3` and run `python3 setup.py install --user`
* Again, if you had another version of Python installed already, do `\path\to\your\Python\Python39\python.exe setup.py install --user`
* Verify that the patch worked correctly by navigating to `%APPDATA$\Python\Python39\site-packages\rvspeclib` and make sure that the Python files say Edited by Everett Mastracci at the top (except for `__init__.py` and `libdata.py`).
11. If you downloaded the extra template library, extract it into `ravespan_p3/work`.
12. You may now move the work directory wherever you like.
13. To launch RaveSpan, run `\path\to\your\Python\Python39\python.exe %APPDATA$\Python\Python39\Scripts\rvspec` from your work directory. You may want to create a batch file to do this for you.

### Linux
Will create guide if people want it. Mostly the same process as Windows.

## FITS File Conversion
RaveSpan cannot parse FITS files, it needs to convert them into raw data first. To do this, run `python3 specconv.py yourmode yourfile.fits yourfile2.fits ...` in your work directory. The two mode are object and template. I have set up the conversion script to attempt to automatically pick up header data to format the file title the way RaveSpan requires, but your mileage may vary. Ensure that the file output follows the form `OBSERVATIONDATE_INITIALWAVELENGTH_RESOLUTION_BARYCORRECTION_INSTRUMENT` for objects and `OBJECTNAME_INITIALWAVELENGTH_RESOLUTION` for templates.
* Will probably require you to `pip install astropy`

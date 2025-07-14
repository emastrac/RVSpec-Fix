"""
######################################################################
######                   RaveSpan functions                     ######
######################################################################
"""
#################################################################
# This file is part of the RaveSpan program
# Copyright (C) 2024 Bogumil Pilecki
#
# Edited by Everett Mastracci
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

RV_ERR_DEF = 0.2

import inspect
from scipy.interpolate import splrep, splev
from scipy.optimize import leastsq
from scipy.linalg import svd

from PyQt5 import QtWidgets, QtCore

import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import cm, ticker, rcParams

from math import log, sqrt, cos, sin, tan, atan, pi
import numpy as np
from numpy.random import randn
import sys, os, time
from mpfit.mpfit import mpfit

dpi = 2.*pi
vlight = 299792.458 # km/s

params = {'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10}
rcParams.update(params)


class GenericThread(QtCore.QThread):
    def __init__(self, function, *args, **kwargs):
        QtCore.QThread.__init__(self)
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def __del__(self):
        self.wait()

    def run(self):
        self.function(*self.args,**self.kwargs)
        return


class MyNavigationToolbar_old(NavigationToolbar):
    def __init__(self, plotCanvas, parent=None, orientation='horizontal', coordinates=True):
        if orientation == "vertical":
            coordinates = False
        NavigationToolbar.__init__(self, plotCanvas, parent, coordinates)
        self.setIconSize(QtCore.QSize(25,25))
        act = self.actions()
        if len(act)>10-int(coordinates==False):
            self.removeAction(act[8])
        self.removeAction(act[7])
        self.removeAction(act[1])
        self.removeAction(act[2])
        if orientation == "vertical":
            self.setOrientation(QtCore.Qt.Vertical)
            self.setFixedWidth(28)
            self.removeAction(act[6])
            self.removeAction(act[3])
        else:
            self.setFixedHeight(28)
            self.setMaximumWidth(370)
            self.removeAction(act[6])   
            self.removeAction(act[3])   

class MyNavigationToolbar(NavigationToolbar):
    def __init__(self, plotCanvas, parent=None, orientation='horizontal',
                 acts=[0,4,5,9], widgets=None, coordinates=False, icon_size=25,
                 widget_size=20, fixed_size=28, max_size=370):

        NavigationToolbar.__init__(self, plotCanvas, parent, coordinates)
        self.setIconSize(QtCore.QSize(icon_size,icon_size))

        act = self.actions()
        for i in range(10):
            self.removeAction(act[i])

        if widgets is not None:
            for w in widgets:
                w.setFixedWidth(widget_size)
                self.addWidget(w)  

        for a in acts:
            self.addAction(act[a])      

        if orientation[:1].lower() == "h":
            self.setOrientation(QtCore.Qt.Horizontal)
            self.setFixedHeight(fixed_size)
            self.setMaximumWidth(max_size)
        else:
            self.setOrientation(QtCore.Qt.Vertical)
            self.setFixedWidth(fixed_size)
            

def set_QSpin(parent, val,rmin=0.0,rmax=1.0,step=0.1,decimals=None,suffix=None):
    obj = QtWidgets.QDoubleSpinBox(parent)
    if decimals is not None:
        obj.setDecimals(decimals)    
    obj.setRange(rmin,rmax)
    obj.setValue(val)
    obj.setSingleStep(step)
    if suffix is not None:
        obj.setSuffix(suffix)        
    return obj

def set_QIntSpin(parent, val,rmin=0,rmax=9,step=1,suffix=None):
    obj = QtWidgets.QSpinBox(parent)
    obj.setRange(rmin,rmax)
    obj.setValue(val)
    obj.setSingleStep(step)
    if suffix is not None:
        obj.setSuffix(suffix)
    return obj


####################################################################


def ensure_dir_exists(tdir):
    if not os.path.exists(tdir):
        ndir = "/".join(tdir.split('/')[:-1])
        ensure_dir_exists(ndir)
        try:
            os.mkdir(tdir, 0o755)
        except:
            print("Cannot create dir:", tdir)
            sys.exit(1)  

def np_fix_hjd(hjd):
    hjd[hjd<10000.]  += 2450000.
    hjd[hjd<100000.] += 2400000.
    return hjd


def fix_hjd(hjd):
    if hjd<10000.:
        hjd += 2450000.
    elif hjd<100000.:
        hjd += 2400000.
    return hjd

full_hjd = fix_hjd

rms = lambda x: (sum(x**2)/len(x))**0.5

def ymdhms2jd(*d):
    """ gregorian date to julian date according to: http://www.tondering.dk/claus/cal/julperiod.php """
    defdate = [2000, 1, 1, 12.0, 0.0, 0.0]
    if type(d) is tuple:
        fulldate = list(d) + defdate[len(d):]
    else:
        fulldate = defdate
    year, month, day, hour, minute, second = fulldate    
    a = (14-month)//12
    y = year + 4800 - a
    m = month + 12*a - 3
    jd = day + (153*m+2)//5 + 365*y + y//4 - y//100 + y//400 - 32045
    frac = (hour - 12. + (minute + second/60.)/60.)/24.
    return jd+frac

def jd2ymd(jd, hms=False):
    """ julian date to gregorian date according to: http://www.tondering.dk/claus/cal/julperiod.php
        -> error: a=jd+32044.5; c=a-(); """
    a = jd + 32044.5
    b = (4*a + 3) // 146097
    c = a - (146097*b)//4
    d = (4*c + 3) // 1461
    e = c - (1461*d)//4
    m = (5*e + 2) // 153

    day = e - (153*m + 2)//5 + 1
    month = m+3 - 12*(m//10)
    year = 100*b + d-4800 + (m//10)

    if hms:
        rest = (day*24)%24.0
        h = int(rest)
        rest = (rest-h)*60
        m = int(rest)
        rest = (rest-m)*60
        s = int(rest)
        return (year,month,day,h,m,s)
    return (year,month,day)

fitbase = lambda x, a: a
fitzeros  = lambda x, *p: 0.0*x
fitgauss = lambda x, v, sig, h: h * np.exp(-0.5 * ((x - v)/sig)**2)
def fitrotational(x, v, vrot, area, limbdark):
    c0 = 1.-limbdark/3.
    c1 = 2. * (1.-limbdark) / (pi * c0)
    c2 = limbdark / (2. * c0)
    xpar = 1. - ((x-v)/vrot)**2
    xpar[xpar<0.0] = 0.0
    rotprof = area * (c1*np.sqrt(xpar) + c2*xpar) / vrot
    return rotprof

gauss = lambda x, mu, sig: np.exp(-0.5 * ((x - mu)/sig)**2) / (sqrt(2.*pi)*sig)


def fit_orbital(x, hjd0, per, v0, k, ecc, aop):
    if ecc>0.0:
        x = x - per * (aop - pi/2.)/dpi   
        ph = (x - hjd0) / per % 1.0 
        v = get_true_anomaly(0.0, 0.0, ecc, M=ph*dpi)%dpi
        rv_shape = (np.cos(v+aop) + ecc*np.cos(aop)) 
    else:
        ph = (x-hjd0)/per%1.0
        rv_shape = - np.sin(dpi*ph)
    return v0 + k*rv_shape

def fit_starpuls(x, hjd0, per, a0, *terms):
    puls = a0
    nterms = len(terms)
    if nterms>0:
        if nterms%2==0:
            print("Error: bad number of parameters.")
            raise ValueError
        ph = (x-hjd0)/per%1.0*dpi
        a1 = terms[0]
        puls += a1*np.cos(ph)
        terms = np.array(terms[1:]).reshape(nterms//2,2)
        for j, (ai,bi) in enumerate(terms):
            i = 2+j
            puls += ai*np.cos(i*ph) + bi*np.sin(i*ph)
    return puls
    

fit_orbital_1 = lambda x, per,hjd0,ecc,aop,v0,k: fit_orbital(x, hjd0, per, v0, k, ecc, aop)
fit_orbital_2 = lambda x, per,hjd0,ecc,aop,v0,k: fit_orbital(x, hjd0, per, v0, -k, ecc, aop)
fit_thirdbody_A = lambda x, per,hjd0,ecc,aop,k: fit_orbital(x, hjd0, per, 0.0, k, ecc, aop)
fit_thirdbody_B = lambda x, per,hjd0,ecc,aop,k: fit_orbital(x, hjd0, per, 0.0, -k, ecc, aop)


def get_axis_lims(data,marg=0.04,minspan=0.0,vmin=None,vmax=None):
    if len(data) == 0:
        dmin, dmax = 0.0, 1.0
    else:
        dmin, dmax = min(data), max(data)
    dspan = dmax - dmin
    if dspan < minspan:
        spancor = (minspan-dspan)/2.
        dmin -= spancor
        dmax += spancor
    dmarg = dspan*marg
    if vmin is not None:
        dmin = min(dmin, vmin)
    if vmax is not None:
        dmax = max(dmax, vmax)
    return dmin-dmarg, dmax+dmarg

def get_up_std(x,x0=0.0,nmin=2):
    xp = x[x>x0]-x0
    lp = len(xp)
    if lp < nmin:
        return None
    return sqrt(sum(xp**2)/lp)


def create_verts_from_errs(x,y,err):
    n = len(x)
    verts = np.c_[x,y+err,x,y-err].reshape(n,2,2)
    return verts


def find_vertex(x1,x2,x3,y1,y2,y3):
    d = (x1 - x2) * (x1 - x3) * (x2 - x3)
    a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / d
    b = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / d
    c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / d
    xv = -b / (2.*a)
    yv = c - b*b/(4*a)
    return xv,yv

def last_not_None(alist):
    """ returns last non-None value from the list """
    last = None
    llen = len(alist)
    for i,el in enumerate(reversed(alist)):
        if el is not None:
            last = llen - i - 1
            break
    return last


def convolve_with_rotational_profile(spec, f0, res, vsini, limbdark=0.6, step_size=250):
    """ Rotational profile from:
        Gray, D. F. 1992, The observation and analysis of stellar photospheres,
        2nd ed. (Cambridge University Press)
        
        Spectrum is cut into N pieces each of a length of step_size angstroms and rotational profiles is applied to each part using constant (average) wavelength. It is faster than calculation of rotational profile for every wavelength in a spectrum. """
        
    ntotal = len(spec)
    newspec = np.ones(ntotal) # ones is not defined
    nstep = int(step_size/res) + 1

    c0 = 1.-limbdark/3.
    c1 = 2. * (1.-limbdark) / (pi * c0)
    c2 = limbdark / (2. * c0)

    for i,j in zip(range(0,ntotal-nstep,nstep), range(nstep,ntotal,nstep)):
        if j >= ntotal-nstep: j = ntotal
        l0 = f0 + res*(i+j)/2.
        prof_radius = l0*vsini/vlight 
        nran = int(prof_radius/res)+1
        xprof = np.arange(-nran, nran+1)*res/prof_radius # arrange is not defined
        xpar = 1.-xprof**2
        xpar[xpar<0.0] = 0.0
        rotprof = (c1*np.sqrt(xpar) + c2*xpar) / prof_radius # previously used stdlib sqrt
        rotprof /= sum(rotprof)
        nmin = max(0,i-nran)
        nmax = min(ntotal,j+nran)
        imin = max(0,i)
        jmax = min(ntotal,j)
        newspec[imin:jmax] = np.convolve(spec[nmin:nmax], rotprof, 'same')[imin-nmin:jmax-nmin] # convolve is not defined

    return np.array(newspec,dtype=np.float32) # array, float32 not defined



###################################################################################        
###################################################################################


def get_ecc_anomaly(P,t,e, tp=0.0, M=None):
    if M is None: M = get_mean_anomaly(P,t, tp)
    E = np.empty_like(M)
    for i,Mi in enumerate(M):
        E[i]=newton(lambda xE,yM,ye: xE-ye*sin(xE)-yM, Mi, lambda xE,yM,ye: 1.0-ye*cos(xE), args=(Mi,e))
    return E

def get_mean_anomaly(P,t,tp=0.0):
    return 2.0*pi*(t-tp)/P

def get_true_anomaly(P,t,e, tp=0.0, M=None):
    E=get_ecc_anomaly(P,t,e, tp,M)
    v = 2.*np.arctan(sqrt((1.+e)/(1.-e))*np.tan(E/2.))
    return v    

    
###############################           NUMERICAL METHODS             ##############################################
# Netwon-Raphson method taken from scipy.optimize.minpack.py
def newton(func, x0, fprime=None, args=(), tol=1.48e-8, maxiter=50):
    """Given a function of a single variable and a starting point,
    find a nearby zero using Newton-Raphson.

    fprime is the derivative of the function.  If not given, the
    Secant method is used.
    """
    if fprime is not None:
        p0 = x0
        for iter in range(maxiter):
            myargs = (p0,)+args
            fval = func(*myargs)
            fpval = fprime(*myargs)
            if fpval == 0:
                print("Warning: zero-derivative encountered.")
                return p0
            p = p0 - func(*myargs)/fprime(*myargs)
            if abs(p-p0) < tol:
                return p
            p0 = p
    else: # Secant method
        p0 = x0
        p1 = x0*(1+1e-4)
        q0 = func(*((p0,)+args))
        q1 = func(*((p1,)+args))
        for iter in range(maxiter):
            if q1 == q0:
                if p1 != p0:
                    print("Tolerance of %s reached" % (p1-p0))
                return (p1+p0)/2.0
            else:
                p = p1 - q1*(p1-p0)/(q1-q0)
            if abs(p-p1) < tol:
                return p
            p0 = p1
            q0 = q1
            p1 = p
            q1 = func(*((p1,)+args))
    raise RuntimeError("Failed to converge after %d iterations, value is %s" % (maxiter,p))
    
    
    
    

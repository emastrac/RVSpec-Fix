"""
######################################################################
######                       RV ANALYSIS                        ######
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

from __future__ import print_function
from PyQt5 import QtWidgets, QtCore, QtGui
from .libcommon import *


class SPECANAL(QtWidgets.QMainWindow):
    """ Analysis Window - CCF, TODCOR, Broadening Function profiles and velocity
        determination. """
    def __init__(self, parent):
        super(SPECANAL, self).__init__(parent)
        self.setParent(parent)
        self.main = parent

        self.setWindowTitle("RV Analysis")
        
        self.main_frame = QtWidgets.QWidget()

        self.dpi = 90
        bgcolor = [v / 255.0 for v in self.palette().color(QtGui.QPalette.Window).getRgb()]
        self.fig = Figure((6., 4.0), facecolor=bgcolor, edgecolor=bgcolor, linewidth=-1, dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        
        self.method2ax = {"ccf": "1d", "bf": "1d", "todcor": "2d"}
        self.ax, self.daxes = {}, ['1d', '2d']
        self.ax['1d'] = self.fig.add_axes([0.11,0.11, 0.88 ,0.88])
        self.ax['2d'] = self.fig.add_axes([0.12,0.09, 0.89,0.89])

        for kax in self.daxes:
            self.ax[kax].set_zorder(-1)         
            self.ax[kax].set_visible(False)     
            self.ax[kax].grid(True)

        self.ax['1d'].set_xlabel("RV [km/s]", fontsize=12)
        self.ax['2d'].set_xlabel("RV1 [km/s]")
        self.ax['2d'].set_ylabel("RV2 [km/s]")

        self.bfsmooth_min = 0.5  
        self.bfsmooth_max = 9.95  
        self.bfsmooth_sig = 3.0  
        self.bfsmooth_step = 0.05 
        self.bf_slider = None

        self.bf_slider = QtWidgets.QSlider(QtCore.Qt.Vertical)
        self.bf_slider.setRange(int(self.bfsmooth_min/self.bfsmooth_step), int(self.bfsmooth_max/self.bfsmooth_step))
        self.bf_slider.setValue(int(self.bfsmooth_sig/self.bfsmooth_step))
        self.bf_slider.setTracking(True)
        self.bf_slider.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.bf_slider.setFixedWidth(30)
        self.bf_slider.setEnabled(False)
        self.bf_slider.setToolTip("Set Broadening Function smoothing level.\nWorks only if BF method results are present.")
        self.bf_slider.valueChanged.connect(self.bf_slider_event)
        self.bf_value_label=QtWidgets.QLabel('{0:.2f}'.format(self.bfsmooth_sig),self.main_frame)
        font = self.bf_value_label.font()
        font.setPointSize(10)
        self.bf_value_label.setFont(font)
        self.bf_value_label.setToolTip("Broadening Function - sigma of smoothing gaussian [km/s]")
        
        self.styles = {0: 'b-', 1: 'r-', 2: 'g-', 3: 'm--'}
        
        self.current_method = None
        self.Nv = self.main.Nv
        
        self.plot_1d = None
        self.vmark_1d = [None]*self.Nv
        self.fit_1D_plot = [None]*self.Nv
        
        self.def_fitpars = {'fakegauss': [0.]*3, 'fakerota': [0.]*4}
        self.def_fitpars['gauss'] = [0., 5., 0.6]
        self.def_fitpars['rota']  = [0., 5., 0.6, 0.6]

        self.todcor = None
        self.todcor_mark = [None, None]
        self.todcor_fit_contour = None
        
        self.create_fitfunc_dictionaries()

        self.mpl_toolbar = MyNavigationToolbar(self.canvas,self.main_frame, 'vertical')
        
        lab_vfixed = QtWidgets.QLabel("Fixed:")
        self.vfixed = []
        for i in range(self.Nv):
            self.vfixed += [QtWidgets.QCheckBox(str(i+1))]
            self.vfixed[i].setToolTip('Fix velocity #%d'%(i+1))
            self.vfixed[i].clicked.connect(self.canvas.setFocus)

        self.cb_mark_rvs = QtWidgets.QCheckBox("RVs")
        self.cb_mark_rvs.setChecked(True)
        self.cb_mark_rvs.setToolTip('Mark components radial velocities.')
        self.cb_mark_rvs.clicked.connect(self.canvas.setFocus)

        self.cb_show_fits = QtWidgets.QCheckBox("Fits")
        self.cb_show_fits.setChecked(True)
        self.cb_show_fits.setToolTip('Show fit curves.')
        self.cb_show_fits.clicked.connect(self.canvas.setFocus)

        self.cb_mark_rvs.stateChanged.connect(self.show_velos)
        self.cb_show_fits.stateChanged.connect(self.show_fits)
        
        self.cb_manual = QtWidgets.QCheckBox("Manual")
        self.cb_manual.setToolTip('Allow manual setting of velocities.')
        self.cb_manual.clicked.connect(self.canvas.setFocus)
        
        self.hbox = QtWidgets.QHBoxLayout()
        v1box = QtWidgets.QVBoxLayout()
        v1box.addWidget(self.mpl_toolbar)
        v1box.addWidget(self.bf_slider)
        v1box.addWidget(self.bf_value_label)
        self.hbox.addLayout(v1box)
        v2box = QtWidgets.QVBoxLayout()
        v2box.addWidget(self.canvas)
        h1box = QtWidgets.QHBoxLayout()
        h1box.addSpacing(15)
        h1box.addWidget(lab_vfixed)
        for i in range(self.Nv):
            h1box.addWidget(self.vfixed[i])
        h1box.addStretch()
        h1box.addWidget(self.cb_mark_rvs)
        h1box.addWidget(self.cb_show_fits)
        h1box.addWidget(self.cb_manual)
        v2box.addLayout(h1box)
        self.hbox.addLayout(v2box)
        
        self.main_frame.setLayout(self.hbox)
        self.setCentralWidget(self.main_frame)   

        #self.resize(534,330)

        self.canvas.mpl_connect('button_press_event', self._mpressed)
        self.canvas.mpl_connect('key_press_event', self._kpressed)
        self.canvas.setFocusPolicy(QtCore.Qt.WheelFocus)

    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()

    def _mpressed(self, event):
        self.activateWindow()
        if event.inaxes:
            if event.button == 1 and self.current_method in ["ccf", "bf"]:
                metpow = self.data[np.argmin(abs(self.xdata-event.xdata))]
                print("v = %.3f   (%s power = %.3f)"%(event.xdata, self.current_method, metpow) )
            if self.cb_manual.isChecked():
                if self.current_method in ["ccf","bf"] and event.button in [2,3]:
                    mi = (event.button-1)//2
                    fitenabled = (event.key!="shift")
                    fixother = (event.key=="control")

                    self.setvelo_1D_click(event.xdata, mi, fitenabled=fitenabled, fixother=fixother)
                    
                elif self.current_method == "todcor" and event.button == 3:
                    vx = event.xdata
                    vy = event.ydata
                    if self.main.cb_fit.isChecked():
                        vx, vy = self.click_fit_todcor(vx, vy)
                        print("clicked: %.3f, %.3f"%(event.xdata, event.ydata) )
                    self.main.curve.set_selected_velocity( 0, vx, draw=False)
                    self.main.curve.set_selected_velocity( 1, vy )
                    self.plot_2D_marks((vx,vy), draw=True)
                    print("v1=%.3f"%vx)
                    print("v2=%.3f"%vy)

    def _kpressed(self,event):
        if self.cb_manual.isChecked() and event.inaxes:
            if self.current_method in ["ccf","bf"]:
                if event.key in ['1', '!', 'q']:
                    fitenabled = event.key in ['1', 'q']
                    fixother = (event.key == 'q')
                    res = self.setvelo_1D_click(event.xdata, 0, fitenabled=fitenabled, fixother=fixother)
                elif event.key in ['2', '@', 'w']:
                    fitenabled = event.key in ['2', 'w']
                    fixother = (event.key == 'w')
                    res = self.setvelo_1D_click(event.xdata, 1, fitenabled=fitenabled, fixother=fixother)
                elif event.key in ['3', '#', 'e']:
                    fitenabled = event.key in ['3', 'e']
                    fixother = (event.key == 'e')
                    res = self.setvelo_1D_click(event.xdata, 2, fitenabled=fitenabled, fixother=fixother)
                elif event.key in ['4', '$', 'r']:
                    fitenabled = event.key in ['4', 'r']
                    fixother = (event.key == 'r')
                    res = self.setvelo_1D_click(event.xdata, 3, fitenabled=fitenabled, fixother=fixother)
                elif event.key == "f":
                    pass
                elif event.key == "c":
                    pass

    def bf_slider_event(self, value):
        if self.current_method == "bf":
            self.bfsmooth_sig = value*self.bfsmooth_step
            sigma = self.bfsmooth_sig
            self.bf_value_label.setText('{0:.2f}'.format(sigma))
            self.data = self.gauss_smooth(self.bfdata, sigma)
            self.update_bf()

    def show_bf_smooth_slider(self):
        self.bf_slider.setEnabled(True)
        self.bf_slider.setToolTip("Broadening Function smoothing level.")


    def hide_bf_smooth_slider(self):
        self.bf_slider.setEnabled(False)
        self.bf_slider.setToolTip("Broadening Function smoothing level.\nNOT ACTIVE with other methods !")

    def get_fixed_v(self, force_ifixother=None):
        fixed_v = np.ones(self.Nv, bool)
        
        if force_ifixother is not None:
            fixed_v[force_ifixother] = False
        else:
            for i in range(self.Nv):
                if not self.vfixed[i].isChecked():
                    fixed_v[i] = False

        return fixed_v


    def write_velos(self, velos, velerrs=None, mi=None):
        """ Writes velocities to the terminal. """
        velos = self.make_vlist(velos,mi)
        if velerrs is None:
            velerrs = [None]*len(velos)
        else:
            velerrs = self.make_vlist(velerrs,mi)
        for i,(v,ve) in enumerate(zip(velos,velerrs)):
            v_str = ""
            if v is not None and type(v) is not str:
                v_str += "v%d = %7.3f"%(i+1, v)        
            if ve is not None:
                v_str += "   v%d_err = %.3f"%(i+1, ve)
            print(v_str)

    def setvelo_1D_click(self, v, vind, fitenabled=True, fixother=False):
        """ vind - component number """
        if fitenabled and self.main.cb_fit.isChecked():
            xi = np.argmin( np.abs(self.xdata-v) )
            v,ve = self.fit_1D_profile(self.xdata, self.data, xi, vind, fixother=fixother) 
        else:
            tmp = ['keep']*self.Nv
            tmp[vind] = v
            v = tmp
            ve = None

        self.set_rvc_1D_marks(v, ve, vind=vind, draw=True)
        self.write_velos(v,ve,vind)

    def update_1D_mark(self, num, val):
        if self.vmark_1d[num] is None:
            self.vmark_1d[num], = self.curax.plot([val]*2,[-1.,2.1],\
                        self.styles[num], visible=self.cb_mark_rvs.isChecked())
        else:
            self.vmark_1d[num].set_xdata([val]*2)
            self.vmark_1d[num].set_ydata([-1.,2.1])

    def make_vlist(self, v, vind):
        """ if v is a scalar, create a v list => always return a list ! """

        vlist = v
        if np.isscalar(v):
            vlist = [None]*self.Nv
            vlist[vind] = v
        return vlist 

    def set_rvc_1D_marks(self, velos, verrors=None, setvelos=True, vind=None, draw=False):
        velos = self.make_vlist(velos, vind)
        if setvelos:
            self.main.curve.set_selected_velocities(velos, verrors)
        for i,v in enumerate(velos):
            if v is None or type(v) is str: continue
            self.update_1D_mark(i, v)
        if draw:
            self.curax.set_xlim(self.xmin, self.xmax)
            self.curax.set_ylim(self.ymin, self.ymax)
            self.canvas.draw()

    def clean_1D_marks(self, mask=None, remove=False):
        if mask is None: mask = [1]*self.Nv
        for i in range(self.Nv):
            if mask[i] and self.vmark_1d[i] is not None:
                if remove:
                    self.vmark_1d[i].remove()
                    self.vmark_1d[i] = None
                else:
                    self.vmark_1d[i].set_xdata([])
                    self.vmark_1d[i].set_ydata([])

    def plot_1D_fit(self, num, x,y):
        if self.fit_1D_plot[num] is None:
            self.fit_1D_plot[num], = self.curax.plot(x,y, self.styles[num], zorder=2, \
                                    visible=self.cb_show_fits.isChecked())
        else:
            self.fit_1D_plot[num].set_xdata(x)
            self.fit_1D_plot[num].set_ydata(y)

    def clean_1D_fits(self, mask=None, remove=False):
        if mask is None:
            mask = [1]*self.Nv
        elif type(mask) is int:
            mask = [1 if mask==i else 0 for i in range(self.Nv)]  
        for i in range(self.Nv):
            if mask[i] and self.fit_1D_plot[i] is not None:
                if remove:
                    self.fit_1D_plot[i].remove()
                    self.fit_1D_plot[i] = None
                else:
                    self.fit_1D_plot[i].set_xdata([])
                    self.fit_1D_plot[i].set_ydata([])

    def get_fit_vals(self,fpoly,xmin,xmax,n):
        fxvals = np.linspace(xmin, xmax, n, endpoint=True)
        return fxvals, np.polyval(fpoly, fxvals)



    def create_fitfunc_dictionaries(self):
        """ Creates description for predefined functions later used in multifunction fitting method. """
        self.dfuncs = {"base": fitbase, "gauss": fitgauss, "fakegauss": fitzeros, \
                       "rota": fitrotational, "fakerota": fitzeros}
        self.dfpars = {"base": ["a"], "gauss": ("v","sig","h"), "fakegauss": "gauss" , \
                       "rota": ["v","vrot","h","ld"], "fakerota": "rota"}
        self.dfpar_lims  =  {"base": [[-0.1,0.5]], "gauss": [[None,None],[0.1,None],[0.0,1.]], \
                             "rota": [[None,None],[0.1,None],[0.0,20.0],[0.0,1.0]]}
        
        for k in self.dfuncs.keys():
            if type(self.dfpars[k]) is str:
                self.dfpars[k] = self.dfpars[self.dfpars[k]]
        
        for k in self.dfuncs.keys():
            if k not in self.dfpar_lims:
                self.dfpar_lims[k] = [[None,None] for p in self.dfpars[k]]
        
        self.dfpar_islim = {}
        self.dfnpars = {}
        for k in self.dfuncs.keys():
            npars = len(self.dfpars[k])
            self.dfnpars[k] = npars
            pl = self.dfpar_lims[k]
            self.dfpar_islim[k] = [[pl[i][0] is not None, pl[i][1] is not None]  for i in np.arange(npars)]

    def nfun_eval(self, x, pars, funcs):
        funval, shi = 0.0, 0
        for fkey in funcs:
            npar = self.dfnpars[fkey]
            fun  = self.dfuncs[fkey]
            funval += fun(x, *pars[shi:shi+npar])
            shi += npar
        return funval       

    def mperrfun_nfun(self,pars,x,y,funcs,fjac=None):
        funval = self.nfun_eval(x, pars, funcs)
        dev = (y - funval)/0.01
        return [0, dev] 

    def mpfit_nfun(self,x,y,funcs,pars,fixpars):
        order, limited, limits = [], [], []

        n = 0
        for i,fkey in enumerate(funcs):
            npar = self.dfnpars[fkey]
            n += npar
            okeys = [self.dfpars[fkey][j]+str(i) for j in range(npar)]
            order.extend(okeys)
            limited += self.dfpar_islim[fkey]
            limits  += self.dfpar_lims[fkey]
          
        pari = [{'value':pars[i], 'fixed':fixpars[i], 'limited': limited[i], \
                 'limits': limits[i]} for i in range(n)]
        fa = {'x':x, 'y':y, 'funcs': funcs}
        m = mpfit( self.mperrfun_nfun,  functkw=fa, parinfo=pari, iterfunct=None)
        if (m.status <= 0): 
            print('error message = ', m.errmsg)
        return m.params, m.perror


    def fit_n_funcs(self, velos0, fun="rota", fakefun="fakerota", fixed_v=None):
        n = len(velos0)
        np = self.dfnpars[fun]
        if n == 0 or velos0.count(None)==n:
            return [], [], []
        funcs = ["base"]
        pars, fixed = [0.0], [0]
        for i,v in enumerate(velos0):
            if v is None:
                funcs += [fakefun]
                pars += self.def_fitpars[fakefun]
                fixed += [1] * np
            else:
                funcs += [fun]
                pars += self.def_fitpars[fun]
                pars[3+np*i] -= 0.1*i
                pars[2+np*i] *= self.res    
                
                if v is not None:
                    pars[1+np*i] = velos0[i]
                fixed += [0] * np
                if fixed_v is not None:
                    if fixed_v[i]: fixed[1+np*i] = 1 


        pars, epars = self.mpfit_nfun(self.xdata, self.data, funcs, pars, fixed)

        profiles, velos, velerrs = [], [], []
        usefun = [fakefun, fun]
        for i in range(n):
            if velos0[i] is None:
                profiles += [None]
                velos += [None]
            else:
                funcs = ["base"] + [usefun[int(i==j)] for j in range(n)]
                xvec = self.xdata
                yvec = self.nfun_eval(xvec, pars, funcs)
                v = pars[1 + i*np]
                ve = epars[1 + i*np]
                profiles += [[xvec, yvec]]
                velos += [v]
                velerrs += [ve]
                if fun == "rota":
                    print("v%drot: %5.2f"%(i+1, pars[2+i*np]), end=' ')
                    print("| S%d = %5.2f"%(i+1, pars[3+i*np]) )
                elif fun == "gauss":
                    sig = pars[2+i*np]
                    h = pars[3+i*np]
                    print("v%dfwhm/2: %5.2f"%(i+1, 2.355*sig/2.), end=' ')
                    print("| S%d = %5.2f"%(i+1, sig*h*sqrt(2*pi)) )

        return profiles, velos, velerrs
   
    def fit_polynomial(self, x, y, method="3ppoly", npoints=1000):
        fpoly = np.polyfit(x,y,4)
        xmin, xmax = np.min(x), np.max(x)
        fxvals, fvals = self.get_fit_vals(fpoly, xmin, xmax, npoints)
        imax = np.argmax(fvals)
        if imax==0 or imax==len(fvals)-1:
            vfmax = fxvals[imax]
        elif method == "hires":
            span = (xmax - xmin)/10.
            xcen = fxvals[imax]
            fxv,fv = self.get_fit_vals(fpoly, xcen-span, xcen+span, npoints)
            vfmax = fxvals[np.argmax(fvals)]
        else:
            vfmax,fmax = find_vertex(*np.r_[fxvals[imax-1:imax+2], fvals[imax-1:imax+2]])
        return fxvals, fvals, vfmax

    def fit_1D_profile(self, x, y, im, vind, fixother=False):
        """ convenience function for fit_1d_profiles if only 1 velocity is to be set """
        imax2 = [None]*self.Nv
        imax2[vind] = im
        ifixother = None
        if fixother: ifixother = vind
        velos, velerrs = self.fit_1D_profiles(x, y, imax2, ifixother=ifixother)
        return velos, velerrs

    def fit_1D_profiles(self, x, y, imax2, ifixother=None, nmin=10, npar=7):
        """ Returns improved velocities for maxima given in ixmax2 list.
        It may be of length 1 or 2. If element of the list is None than it will be igored. """
        method = self.main.preferences.get_fit_function(self.current_method)
        velos = []
        if method == "poly":
            for i,im in enumerate(imax2):
                if im is None:
                    velos += ['keep']
                    continue
                n = max(nmin, int(npar/self.res))
                s = slice( max(0, im-n), min(x.size, im+n+1))
                fxvals, fvals, vfmax = self.fit_polynomial(x[s], y[s])
                self.plot_1D_fit(i, fxvals, fvals)
                velos += [vfmax]
            velerrs = [None]*len(velos)
        elif method in ["gauss", "rota"]:
            velos0 = [None]*self.Nv
            notfixed = None
            fixed_v = self.get_fixed_v(ifixother)
            for i,im in enumerate(imax2):
                if im is None or fixed_v[i]:
                    velos0[i] = self.main.curve.get_selected_velocity(i)    
                else:
                    velos0[i] = x[im]
            
            profiles, velos, velerrs = self.fit_n_funcs(velos0, fun=method,
                                            fakefun="fake"+method, fixed_v = fixed_v)
            for i, prof in enumerate(profiles):
                if prof is None:
                    self.clean_1D_fits(i)
                else:
                    (fxvals, fvals) = prof
                    self.plot_1D_fit(i, fxvals, fvals)
        else:
            print("Method", method, "not implemented.")
            return [0]        
        return velos, velerrs

    def filter_out_too_close(self, max_list, v_values):
        """ Clear list of detected velocities """
        cmax_list = max_list.copy()
        
        for imax in max_list:
            if imax in cmax_list:
                diff = abs(v_values[cmax_list]-v_values[imax])
                cond = (diff==0.0) | (diff>5*self.res)
                cmax_list = cmax_list[cond]
        
        return cmax_list

    def auto_find_1D_max(self, xmin, xmax, nmin=10, npar=7, detection_limit={"ccf":0.2, "bf":0.01}):
        """ Find up to N maxima, if fitting is enabled -> get improved coordinates from fitted functions. """
        cond = (self.xdata>xmin) & (self.xdata<xmax)
        cxdata, cdata = self.xdata[cond], self.data[cond]
        di = np.diff( np.sign( np.diff(cdata) ) )
        maxcond = (di < 0) & (cdata[1:-1]>detection_limit[self.current_method])
        imaxs = np.where(maxcond)[0]+1

        argAll = np.argsort(cdata[imaxs])[::-1]
        imaxAll = imaxs[argAll]
        
        if len(imaxAll) == 0:
            return [],[]

        imaxClear = self.filter_out_too_close(imaxAll, cxdata)

        n_auto_v = self.main.preferences.get_n_auto_vel()

        imaxN = imaxClear[:n_auto_v]
        
        if self.main.preferences.get_v_from_model():
            mvelos = self.main.curve.get_model_for_current_hjd()
            new_imaxN = []
            for mv in mvelos[:2]:
                ai = np.argmin(abs(cxdata - mv))
                new_imaxN.append(ai)
            imaxN = np.array(new_imaxN,int)
        elif self.main.preferences.get_ID_from_model():
            velos = cxdata[imaxN]
            mvelos = self.main.curve.get_model_for_current_hjd()
            new_imaxN = []
            for mv in mvelos[:len(velos)]:
                ai = np.argmin(abs(velos - mv))
                velos[ai] = -99999.
                new_imaxN.append(imaxN[ai])
            imaxN = np.array(new_imaxN,int)

        if self.main.cb_fit.isChecked():
            velos,velerrs = self.fit_1D_profiles(cxdata, cdata, imaxN)
        else:
            velos = cxdata[imaxN]
            velerrs = [None]*len(velos)

        self.write_velos(velos, velerrs)
        return velos,velerrs


    def set_current_method(self, method):
        self.clean_plots_for(method)
        
        self.current_method = method
        self.current_axes = self.method2ax[method]
        self.curax = self.ax[self.current_axes]
        for kax in self.daxes:
            if kax == self.current_axes:
                self.curax.set_zorder(1)
                self.curax.set_visible(True)
            else:
                self.ax[kax].set_zorder(-1)
                self.ax[kax].set_visible(False)


    def plot_cached_data(self, data, vleft, res, method):
        xdata = np.arange(len(data))*res + vleft
        if method == "ccf":
            self.plot_ccf(xdata, data, res, from_cache=True)
        elif method == "bf":
            self.plot_bf(xdata, data, res, from_cache=True)

    def clean_1D_plot(self):
        if self.plot_1d is not None :
            self.plot_1d.remove()
            self.plot_1d = None

    def clean_plots_for(self, method):
        if method != "bf":
            self.hide_bf_smooth_slider()
        met_axes = self.method2ax[method]
        if self.current_method is not None:
            if self.current_axes == "1d":
                if met_axes == "1d":
                    self.clean_1D_marks()
                    self.clean_1D_fits()
                else: 
                    self.clean_1D_plot()
                    self.clean_1D_fits(remove=True)
                    self.clean_1D_marks(remove=True)
                    newwidth = self.height()+65
                    self.resize( newwidth, self.height())
            elif self.current_axes == "2d":
                if met_axes == "2d":
                    self.clean_2D_marks()
                    self.clean_2D_fit()
                else: 
                    newwidth = self.height()*534./330.
                    self.resize( newwidth, self.height())
                    self.clean_2D_image()
        elif met_axes == "2d":
            newwidth = self.height()+65
            self.resize( newwidth, self.height())

    def plot_1D_data(self, xdata, data):
        if self.plot_1d is None:
            self.plot_1d = self.ax["1d"].plot(xdata, data, 'k-')[0]
        else:
            self.plot_1d.set_xdata(xdata)
            self.plot_1d.set_ydata(data)

    def get_1D_vrange(self, minran=20.):
        v0 = self.main.get_v0_from_field()

        ampl = self.main.get_model_amplitude()
        mul,add = self.main.preferences.get_vrange_pars()

        ran = max( ampl*mul+add, minran)
        vmin, vmax = v0 - ran, v0 + ran
        return vmin, vmax

    def get_bf_yrange(self, data, data_mask=None):
        if data_mask is not None:
            data = data[data_mask]
        
        rmin, rmax = np.nanmin(data), np.nanmax(data)
        yspan = rmax-rmin
        rmax += 0.1*yspan
        rmin -= 0.05*yspan
        return rmin, rmax

    def find_and_set_velos(self, method, xmin, xmax, from_cache=False):
        if from_cache:
            velos = self.main.curve.get_selected_velocities()
            self.set_rvc_1D_marks(velos, setvelos=False)
        else: 
            if self.current_method == "bf":
                ydata = self.bfdata
            else:
                ydata = self.data
            self.main.curve.cache_anal_1D_data(self.xdata, ydata, method)
            if self.main.cb_auto_max.isChecked():
                velos, verrors = self.auto_find_1D_max(xmin, xmax)
                self.set_rvc_1D_marks(velos,verrors)

    """
    ##############################################################
    #                      PLOT 1D METHOD
    ###############################################################
    """

    def plot_1D_method(self, method, xdata, data, from_cache=False, plot_anal=True, yran=None):
        """ plots data for given 1D-method
            plot_anal - if False: only analysis is done"""
        self.xdata = xdata
        self.data = data
        self.xmin, self.xmax = self.get_1D_vrange()
        if yran is None:
            self.ymin, self.ymax = 0., max(1.,max(data)) 
        elif yran == "bf":
            data_mask = (xdata>self.xmin) & (xdata<self.xmax)
            self.ymin, self.ymax = self.get_bf_yrange(data, data_mask)
        else:
            self.ymin, self.ymax = yran
                
        self.plot_1D_data(xdata, data)

        self.find_and_set_velos(method, self.xmin, self.xmax, from_cache=from_cache)

        self.curax.set_ylabel(method.upper(), fontsize=11)
        self.curax.set_xlim(self.xmin, self.xmax)
        self.curax.set_ylim(self.ymin, self.ymax)

        if plot_anal:
            self.canvas.draw()
               

    def plot_ccf(self, xdata, data, res, from_cache=False, plot_anal=True):
        """ Plots CCF data """
        self.set_current_method("ccf")
        self.setWindowTitle("RV Analysis (CCF)")
        self.res = res

        self.plot_1D_method("ccf", xdata, data, from_cache, plot_anal)


    """
    ##############################################################
    #                       BF
    ###############################################################
    """
            
    def plot_bf(self, xdata, data, res, v0shift=0.0, from_cache=False, plot_anal=True):  
        """ Plots Broadening Function data """
        self.set_current_method("bf")
        self.setWindowTitle("RV Analysis (BF)")

        self.res = res
        self.bfdata = data
        data = self.gauss_smooth(self.bfdata, self.bfsmooth_sig)
        
        self.plot_1D_method("bf", xdata, data, from_cache, plot_anal, yran="bf")
        self.show_bf_smooth_slider()
        

    def gauss_smooth(self, data, sigma):
        """Smoothing gaussian - sigma=0.4248*FWHM"""
        sig_in_res_unit = sigma/self.res
        nran = int(6.*sig_in_res_unit)
        nran = max(nran, 5)
        nran = min(nran, len(data)//2)
        gaussprof = gauss(np.arange(-nran, nran+1), 0., sig_in_res_unit)/self.res
        sdata = np.convolve(self.bfdata, gaussprof, 'same')

        return sdata

    def update_bf(self):
        self.plot_1d.set_xdata(self.xdata)
        self.plot_1d.set_ydata(self.data)
        cond = (self.xdata>self.xmin) & (self.xdata<self.xmax)
        cdata = self.data[cond]
        dmin, dmax = min(cdata), max(cdata)
        yspan = dmax - dmin
        self.ymax = min(max( self.ymax, dmax+0.05*yspan), dmax+0.3*yspan)
        self.ymin = max(min( self.ymin, dmin-0.025*yspan), dmin-0.05*yspan)
        self.curax.set_ylim(self.ymin,self.ymax)
        self.canvas.draw()




    """
    ###############################################################
    #                       2D  methods
    ###############################################################
    """

    def plot_todcor(self, data, res, v0shift=0.0, plot_anal=True):
        self.clean_plots_for("todcor")

        self.set_current_method("todcor")
        self.setWindowTitle("RV Analysis (TODCOR)")

        self.data = data
        self.res = res
        self.vshift = v0shift
        ndata = len(data)
        nran = ndata//2
        
        if self.todcor is None:
            self.todcor = self.curax.imshow(data, origin="lower",cmap=cm.gist_heat)
        else:
            self.todcor.set_data(data)

        if self.main.cb_auto_max.isChecked():
            velos = self.auto_find_2D_max()
            for i,v in enumerate(velos):
                self.main.curve.set_selected_velocity( i, v, draw=i) 
            self.plot_2D_marks(velos)
        else:
            self.clean_2D_marks()
        vmin = self.vshift - res*nran
        vmax = self.vshift + res*nran
        self.todcor.set_extent((vmin,vmax,vmin,vmax))
        dspan = np.ptp(data) # ptp was removed from the ndarray class in NumPy 2.0
        self.todcor.set_clim(data.min()+0.0*dspan,data.max())
        self.curax.axis((vmin,vmax,vmin,vmax))
        if plot_anal:
            self.canvas.draw() 


    def plot_2D_marks(self, vals, draw=False):
        if len(vals) == 0: return
        if self.todcor_mark[0] is None:
            mark_rvs = self.cb_mark_rvs.isChecked()
            self.todcor_mark[0], = self.curax.plot([vals[0]]*2,[-500,500],'b-', visible=mark_rvs)
            self.todcor_mark[1], = self.curax.plot([-500,500],[vals[1]]*2,'r-', visible=mark_rvs)
        else:
            self.todcor_mark[0].set_xdata([vals[0]]*2)
            self.todcor_mark[0].set_ydata([-500,500])
            self.todcor_mark[1].set_xdata([-500,500])
            self.todcor_mark[1].set_ydata([vals[1]]*2)
            
        if draw:
            self.canvas.draw()

    def clean_2D_marks(self, remove=False):
        for i in range(2):
            if self.todcor_mark[i] is not None:
                if remove:
                    self.todcor_mark[i].remove()
                    self.todcor_mark[i] = None
                else:
                    self.todcor_mark[i].set_xdata([])
                    self.todcor_mark[i].set_ydata([])

    def plot_2D_fit(self, X, Y, Z):
        if self.todcor_fit_contour is not None:
            self.clean_2D_fit()
        lims = self.curax.axis()
        self.todcor_fit_contour = self.curax.contour(Y,X,Z, cmap=cm.Greys)
        if not self.cb_show_fits.isChecked():
            for cc in self.todcor_fit_contour.collections:
                cc.set_visible(False)
        self.curax.axis(lims)

    def clean_2D_fit(self):
        if self.todcor_fit_contour is not None:
            for cc in self.todcor_fit_contour.collections:
                cc.remove()
            self.todcor_fit_contour = None

    def todcor_profile(self,p,x,y):
        prof = 0.0
        for i in range(16):
            xpow = x**(i//3)
            ypow = y**(i%3)
            prof += p[i]*xpow*ypow
        return prof

    def get_2dfit_velos(self, pars, xmin, xmax, ymin, ymax, n, fitcont=False):
        fx = np.linspace(xmin,xmax,n)
        fy = np.linspace(ymin,ymax,n)
        X,Y = np.meshgrid(fx,fy) # meshgrid is not defined
        x,y = np.ravel(X),np.ravel(Y) # ravel is not defined
        fvals = self.todcor_profile(pars, x,y)
        if fitcont:
            Z=fvals.reshape(X.shape)
            self.plot_2D_fit(X,Y,Z)
        imax = np.argmax(fvals)
        return x[imax], y[imax]

    def fit_2D_profile(self, fdata, nx0,ny0):
        lx,ly = self.data.shape
        nx0 -= lx//2
        ny0 -= ly//2

        lx,ly = fdata.shape
        pars = np.array([1.]+[0.]*15)
        errfun = lambda p,x,y: np.ravel(self.todcor_profile(p,x,y) - fdata) # ravel is not defined
        x = np.arange(nx0, nx0+lx)*self.res + self.vshift
        y = np.arange(ny0, ny0+ly)*self.res + self.vshift
        Y,X = np.meshgrid(y,x) # meshgrid is not defined  

        pars, success = leastsq(errfun, pars, args=(X,Y))
        print("%d params of 2D poly, values from %g to %g"%(len(pars), min(pars), max(pars)))

        xmin, xmax = nx0*self.res + self.vshift, (nx0+lx)*self.res + self.vshift
        ymin, ymax = ny0*self.res + self.vshift, (ny0+ly)*self.res + self.vshift
        xcen, ycen = self.get_2dfit_velos(pars, xmin, xmax, ymin, ymax, 300, fitcont=True)
        xspan = (xmax - xmin)/10.
        yspan = (ymax - ymin)/10.
        xbest, ybest = self.get_2dfit_velos(pars, xcen-xspan, xcen+xspan, ycen-yspan, ycen+yspan, 500)
        velos2D = [xbest, ybest]
        return velos2D

    def click_fit_todcor(self, x, y, nmin=8, npar=4):
        lx,ly = self.data.shape
        xran, yran = lx//2, ly//2
        x -= self.vshift
        y -= self.vshift
        y, x = x, y
        nx = int(x/self.res+0.5) + xran
        ny = int(y/self.res+0.5) + yran
        print("ANAL:", self.data[nx,ny])

        n = max (nmin, int(npar/self.res))
        nxmin = max(0,nx-n)
        nymin = max(0,ny-n)
        sx = slice(nxmin, min(lx,nx+n+1))
        sy = slice(nymin, min(ly,ny+n+1))
        velos = np.array(self.fit_2D_profile(self.data[sx,sy], nxmin, nymin))
        return velos[::-1]
 
    def auto_find_2D_max(self, nmin=8, npar=4, detection_limit=0.97):
        lx,ly = self.data.shape
        xran, yran = lx//2, ly//2
        num = self.data.argmax()
        nx = num//ly
        ny = num%lx
        if self.data[nx,ny]<0.97: return []
        
        if self.main.cb_fit.isChecked():
            n = max (nmin, int(npar/self.res))
            nxmin=max(0,nx-n)
            nxmax=min(lx,nx+n+1)
            nymin=max(0,ny-n)
            nymax=min(ly,ny+n+1)
            velos = self.fit_2D_profile(self.data[nxmin:nxmax,nymin:nymax], nxmin, nymin)
        else:
            nxabs = nx - xran
            nyabs = ny - yran
            velos = np.array([nxabs*self.res,nyabs*self.res]) + self.vshift
        
        for i,v in enumerate(velos):
            print("v%d=%.3f"%(i+1,v) )
        return velos[::-1]

    def clean_2D_image(self):
        pass
    

    """
    ##############################################################
    """

    def set_velos_visibility(self, value):
        self.analysis.show_velos(self.cb_show_velos.isChecked())

    def set_fits_visibility(self, value):
        self.analysis.show_fits(self.cb_show_fits.isChecked())
        
    def show_velos(self, status=True):
        """ Mark velocities in RV Analysis window. """
        if self.current_method in ["ccf","bf"]:
            marks = self.vmark_1d
        elif self.current_method == "todcor":
            marks = self.todcor_mark
        else:
            return
        for i in range(2):
            if marks[i] is not None:
                marks[i].set_visible(status)
        self.canvas.draw()

    def show_fits(self, status=True):
        """ Plot CCF and TODCOR fits in RV Analysis window. """
        if self.current_method in ["ccf","bf"]:
            for i in range(2):
                if self.fit_1D_plot[i] is not None:
                    self.fit_1D_plot[i].set_visible(status)
        elif self.current_method == "todcor":
            if self.todcor_fit_contour is not None:
                for cc in self.todcor_fit_contour.collections:
                    cc.set_visible(status)
        self.canvas.draw()
    


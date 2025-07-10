"""
#############################################################################
###########                     SPECTRUM                     ################
#############################################################################
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


from PyQt5 import QtWidgets, QtCore, QtGui
from .libcommon import *
from numpy.fft import rfft, irfft,fft,ifft
from numpy import max as npmax

class SPECTRUM(QtWidgets.QMainWindow):
    """ Spectrum View Window: spectrum, template, zoom panel """
    def __init__(self, parent):
        super(SPECTRUM, self).__init__(parent)
        self.setParent(parent)
        self.main = parent
        self.analysis = self.main.analysis
        self.sd_dialog = SDDialog(self)
        
        self.setWindowTitle("Spectrum")
        self.init_vals()
        self.load_masks()

        self.main_frame = QtWidgets.QWidget()

        self.dpi = 80
        bgcolor = [v / 255.0 for v in self.palette().color(QtGui.QPalette.Window).getRgb()]
        self.fig = Figure((7.8, 4.7), facecolor=bgcolor, edgecolor=bgcolor, linewidth=-1, dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.ax = self.fig.add_axes([0.06,0.08,0.92 ,0.9])
        
        self.spec_plot, = self.ax.plot([2000,9999],[1,1],'r-',alpha=0.5)
        self.tpl_plot = [None,None]
        self.tpl_plot[0], = self.ax.plot([2000,9999],[1,1],'b-', alpha=0.4)
        self.tpl_plot[1], = self.ax.plot([2000,9999],[1,1],'m-', alpha=0.4, visible=False)
        self.mask_areas = self.ax.fill(ec='k',lw=0.5, visible=True)
        self.zoom_mark, = self.ax.plot([5000,5000],[0.,3.],'k:', visible=False)
        self.ax.set_xlabel("wavelength")

        self.zoomfig = Figure((1.0, 1.0), facecolor=bgcolor, edgecolor=bgcolor, linewidth=-1, dpi=self.dpi)
        self.zoomcanvas = FigureCanvas(self.zoomfig)
        self.zoomcanvas.setParent(self)
        
        self.axzoom = self.zoomfig.add_axes([0.02,0.02,0.96 ,0.96])
        self.axzoom.xaxis.set_ticklabels([])
        self.axzoom.yaxis.set_ticklabels([])
        self.zoomspec, = self.axzoom.plot([4995,5005],[1,1],'r-',alpha=0.55)
        self.zoomtpl, = self.axzoom.plot([4995,5005],[1,1],'b-',alpha=0.33, visible=False)
        self.zoommark, = self.axzoom.plot([5000,5000],[0.,1.2],'k:')
        self.zoomtext = self.axzoom.text(0.61,0.1,"5000 A", size=10, transform=self.axzoom.transAxes)
        self.axzoom.axis(xmin=4995, xmax=5005, ymin=0.0, ymax=1.2)

        self.load_lines()

        self.canvas.draw()
        self.zoomcanvas.draw()

        self.cb_show_spec = QtWidgets.QCheckBox("Spectrum")
        self.cb_show_spec.setToolTip("Show selected spectrum.")
        self.spec_plot.set_visible(self.cb_show_spec.isChecked())   
        self.cb_show_tpl = []
        for i in range(2):
            self.cb_show_tpl.append(QtWidgets.QCheckBox("T%d"%(i+1)))
            self.cb_show_tpl[i].setToolTip('Show template %d in the Spectrum window.'%(i+1))
            state = (i==0)
            self.cb_show_tpl[i].setChecked(state)
            self.tpl_plot[i].set_visible(state)
            self.cb_show_tpl[i].stateChanged.connect(self.sync_show_tpl)
        
        self.cb_normalization = QtWidgets.QCheckBox("Normalized")
        self.cb_normalization.setToolTip("Toggle normalization of selected spectrum.")
        self.cb_normalization.setChecked(True)
        self.cb_zoom_lock = QtWidgets.QCheckBox("Zoom lock")
        self.cb_zoom_lock.setToolTip("Lock zoom window at selected wavelength.")
        self.cb_show_spec.stateChanged.connect(self.sync_show_spec)
        self.cb_normalization.stateChanged.connect(self.sync_normalization)
        self.cb_zoom_lock.stateChanged.connect(self.set_zoommark_visibility)

        self.cb_zoom_stick = QtWidgets.QCheckBox("Stick")
        self.cb_zoom_stick.setChecked(True)
        self.cb_zoom_stick.setToolTip("Stick to spectrum minima - center zoom window\nat spectrum minimum.")
        self.cb_zoom_stretch = QtWidgets.QCheckBox("Stretch")
        self.cb_zoom_stretch.setToolTip("Stretch spectrum vertically in zoom window\nto fill available Y-space.")
        self.cb_zoom_stretch.stateChanged.connect(self.on_zoom_change)
        zoomrange_lab = QtWidgets.QLabel("Zoom:")
        self.zoomrange = set_QSpin(self, 6.,rmin=1.,rmax=99.99,step=0.5)
        self.zoomrange.setToolTip("Set zoom range in angstroms around the pointed wavelength.")
        self.zoomrange.valueChanged.connect(self.on_zoom_change)

        self.cb_update_limits = QtWidgets.QCheckBox("")
        self.cb_update_limits.setChecked(True)
        self.cb_update_limits.setToolTip("Update limits on any change in the plot.")

        self.mpl_toolbar = MyNavigationToolbar(self.canvas, self, widgets=[self.cb_update_limits], coordinates=False)

        gridlay = QtWidgets.QGridLayout()
        gridlay.addWidget(self.cb_show_spec,1,1,1,2)
        for i in range(2):
            gridlay.addWidget(self.cb_show_tpl[i],1,3+i)
        gridlay.addWidget(self.cb_normalization,1,5,1,2)
        gridlay.addWidget(self.cb_zoom_lock,1,7,1,3)
        gridlay.addWidget(self.zoomcanvas,1,9,2,3)

        gridlay.addWidget(self.mpl_toolbar,2,1,1,3)
        gridlay.addWidget(self.cb_zoom_stick,2,4)
        gridlay.addWidget(self.cb_zoom_stretch,2,5,1,2)
        gridlay.addWidget(zoomrange_lab,2,7)
        gridlay.addWidget(self.zoomrange,2,8,1,1)
        
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addLayout(gridlay)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

        self.canvas.mpl_connect('button_press_event', self._mpressed)
        self.canvas.mpl_connect('motion_notify_event', self._mmotion)
        self.canvas.setFocusPolicy(QtCore.Qt.WheelFocus)

    def show_or_activate(self):
        if not self.plotting_enabled:
            self.plotting_enabled = True
            self.canvas.draw()
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()

    def init_vals(self):
        self.xmin, self.xmax = 2000, 9999
        self.ymin, self.ymax = 0, 1.2

        self.plotting_enabled = False
        
        self.dbase = {"tpl": {}, "xtpl": {}, "2tpl_xc": {}, "svd": {}}
        self.dbind = {"tpl": [], "xtpl": [], "2tpl_xc": [], "svd": []}
        self.dbmax = {"tpl":  4, "xtpl":  2, "2tpl_xc":  1, "svd":  2, "spec": 50, "norm_spec": 50, "xspec": 4}
        self.new_object_clear()
        
        self.tpl_specs = [None]*2
        self.xtpls = [None]*2
        self.tpl_f0 = [0,0]
        self.tpl_res = [0,0]
        self.tpl_rot = [0,0]


    def new_object_clear(self):
        self.dbase["spec"], self.dbind["spec"] = {}, []
        self.dbase["xspec"], self.dbind["xspec"] = {}, []
        self.dbase["norm_spec"], self.dbind["norm_spec"] = {}, []
        self.spec = None        
        self.norm_spec = None   
        self.xspec = None       
        self.f0, self.res, self.speclen = None, None, None
        self.zoom = None
        self.sd_res = None


    def load_lines(self):
        self.lines = {
            3933.7: ("Ca","ionised Calcium"), 3968.5: ("Ca","ionised Calcium"),
            3969.5: ("H-eps",),
            4101.7: ("H-delta",), \
            4307.9: ("Fe",), \
            4314.2: ("CH","CH molecule"), \
            4340.5: ("H-gamma",), \
            4383.6: ("Fe",), 4668.1: ("Fe",), \
            4471. : ("HeI",),
            4861.3: ("H-beta",), \
            4957.6: ("Fe",), \
            5183.6: ("Mg",), 5172.7: ("Mg",), 5168.9: ("Mg",), 5167.3: ("Mg",), \
            5269.6: ("Fe",), \
            5460.7: ("Hg",), \
            5875.6: ("He",), \
            5890.0: ("Na",), 5895.9: ("Na",), \
            6276. : ("atmO2->","atmospheric O2"), 6287. : ("<-atmO2","atmospheric O2"), \
            6562.8: ("H-alpha",), \
            6867.2: ("atmO2->","atmospheric O2"), 6884. : ("<-atmO2","atmospheric O2"), \
            7593.7: ("atmO2->","atmospheric O2"), 7621. : ("<-atmO2","atmospheric O2") \
            }
        for wlength in self.lines.keys():
            data = self.lines[wlength]
            self.ax.plot([wlength]*2, [0.0,0.05], '-', c='b')
            self.axzoom.text(wlength, 1.12, data[0], size=9, horizontalalignment='center', \
                                                             verticalalignment='top')

    def setChecked_ifnotNone(self, obj, state):
        if obj is not None:
            obj.setChecked(state)

    def sync_show_spec(self, state):
        """ EVENT. """
        self.setChecked_ifnotNone(self.main.cb_show_spec, state)
        self.show_spec(state)

    def sync_show_tpl(self, state):
        """ EVENT. """
        for i in range(2):
            state = self.cb_show_tpl[i].isChecked()
            if self.main.cb_show_tpl[i] is not None:
                stateR = self.main.cb_show_tpl[i].isChecked()
                if state != stateR:
                    self.main.cb_show_tpl[i].setChecked(state)
                print("SYNCING TPLS - SPEC")
            self.tpl_plot[i].set_visible(state)
        print("PLOTTTING TPL")
        self.canvas.draw()

    def sync_normalization(self, state):
        self.main.cb_normalization.setChecked(state)
        self.plot_data()   

    def remove_from_db(self, idb, key):
        """ remove an object from a database """
        if key in self.dbind[idb]:
            del(self.dbase[idb][key])
            self.dbind[idb].remove(key)

    def add_to_db(self, idb, key, value):
        """ add and object to a database, keep them clean """
        self.dbase[idb][key] = value
        self.dbind[idb].append(key)
        self.db_maintenance(idb)
        

    def db_maintenance(self, idb):
        """ if database 'idb' exceeds size limit: remove the oldest items """
        overflow = len(self.dbind[idb]) - self.dbmax[idb]
        if overflow > 0:
            for i in range(overflow):
                db_id = self.dbind[idb][0]
                del(self.dbase[idb][db_id])
                self.dbind[idb] = self.dbind[idb][1:]
    
    def clear_db_area(self, area):
        self.dbase[area], self.dbind[area] = {}, []
    
    def clear_db_norm(self):
        self.clear_db_area("norm_spec")
    
    def _mpressed(self, event):
        if self.xspec is None:
            return None
        if event.inaxes:
            if event.button > 1:
                if self.cb_zoom_lock.isChecked():
                    self.plot_zoom_window(event.xdata, event.ydata)
                    self.canvas.draw()
                else:
                    self.cb_zoom_lock.setChecked(True)

    def _mmotion(self, event):
        if self.xspec is None or \
           self.cb_zoom_lock.isChecked():
            return None            
        if event.inaxes and event.button is None:
            try:
                self.plot_zoom_window(event.xdata, event.ydata)
            except TypeError:
                print("WARNING!!!")
                print("_mmotion: event", event, event.inaxes, event.button)
                print("plot_zoom_window: xev, self.f0", event.xdata, self.f0)

    def on_zoom_change(self, val):
        self.draw_zoom_window()

    def get_range_indices(self, x, ran):
        nmax = len(self.xspec)
        i0 = int((x-ran-self.f0)/self.res+0.5)
        i1 = int((x+ran-self.f0)/self.res+0.5)
        if i1 < 0:
            return 0, 0
        elif i0 > nmax:
            return nmax, nmax
        return max(i0, 0), min(i1, nmax)


    def plot_zoom_window(self, xev, yev):
        """ plot zoom window for given mouse pointer coords """
        ran = self.zoomrange.value()/2.
        
        if xev < self.xspec[0] - ran or xev > self.xspec[-1] + ran:
            return
            
        if self.cb_zoom_stick.isChecked():
            si0, si1 = self.get_range_indices(xev, 3. * ran)  
            if si1 - si0 > 0:
                icenter = np.argmin(self.norm_spec[si0:si1])+si0
            else:
                return 
        else:
            nmax = len(self.xspec)
            icenter = int((xev-self.f0)/self.res+0.5)
            icenter = max(min(icenter,nmax-1), 0)
       
        xcent = self.xspec[icenter]

        self.draw_zoom_window(xcent, ran)
        

    def draw_zoom_window(self, xcent=None, ran=None):
        """ direct zoom window drawing """
        
        if xcent is None:
            if self.zoom is not None:
                xcent = self.zoom
            else:
                return
        else:
            self.zoom = xcent
            
        if ran is None:
            ran = self.zoomrange.value()/2.

        i0, i1 = self.get_range_indices(xcent, ran)
        if i1-i0<1.5*ran:
            self.zoomspec.set_xdata([])
            self.zoomspec.set_ydata([])
            self.zoomcanvas.draw()
            return

        sli = slice(i0,i1)
        xdata = self.xspec[sli]
        ydata = self.norm_spec[sli]
        if self.cb_zoom_stretch.isChecked():
            ydata = 1.-(1.-ydata)*0.99/max(1.-ydata)
        self.zoomspec.set_xdata(xdata)
        self.zoomspec.set_ydata(ydata)
        
        self.zoomtext.set_text("%.2f"%xcent)
        mark_xdata = [xcent]*2
        self.zoommark.set_xdata(mark_xdata) 
        self.zoom_mark.set_xdata(mark_xdata) 
        
        self.axzoom.axis(xmin=self.f0+i0*self.res, xmax=self.f0+i1*self.res, \
                         ymax=0.2+max(1., np.average(ydata)))
        self.zoomcanvas.draw()                
    

    def set_zoommark_visibility(self,value):
        if value:
            self.zoom_mark.set_visible(True)
        else:
            self.zoom_mark.set_visible(False)
        self.canvas.draw()

    def load_spectrum_file(self, path, ftype):
        if ftype == 'txt':
            data = np.loadtxt(path, unpack=True, dtype=np.float32)
            if len(data)>1: spec = data[1]
            else:           spec = data[0]
        else:
            spec = np.fromfile(path,dtype=np.float32)
        return spec


    """####################### 
       #    LOAD SPECTRUM    #
       #######################"""

    def load_spec(self, path, f0, res, ftype, plot_data=True):
        """ Main loading spectrum method. Loads spectrum either from file or the databse """
        self.path = path
        
        spec, xspec, norm_spec = self._load_any_spec(path, f0, res, ftype)
        self.spec = spec
        self.xspec = xspec
        self.norm_spec = norm_spec
        
        self.f0, self.res, self.speclen = f0, res, len(spec)

        if plot_data:
            self.draw_zoom_window()
            self.plot_data()
            

    def _load_any_spec(self, path, f0, res, ftype):
        """ loads spectrum, either from file or the databse """
        if path in self.dbase["spec"]:
            spec = self.dbase["spec"][path]
        else:
            spec = self.load_spectrum_file(path, ftype)
            spec /= np.median(spec)
            self.add_to_db("spec", path, spec)
        speclen = len(spec)
        xind = (f0,res,speclen)
        if xind in self.dbase["xspec"]:
            xspec = self.dbase["xspec"][xind]
        else:
            xspec = np.linspace(f0, f0+speclen*res, speclen, endpoint=False)
            self.add_to_db("xspec", xind, xspec)
        if path in self.dbase["norm_spec"]:
            norm_spec = self.dbase["norm_spec"][path]
        else:
            norm_spec = self.create_norm_spec_from(xspec, spec)    
            self.add_to_db("norm_spec", path, norm_spec)
        
        return spec, xspec, norm_spec


    def update_norm_spec(self):
        if self.spec is None:
            return
        self.norm_spec = self.create_norm_spec()    
        self.add_to_db("norm_spec", self.path, self.norm_spec)


    def get_signal_to_noise_ratio(self, path, f0, res, ftype):
        
        spec, xspec, norm_spec = self._load_any_spec(path, f0, res, ftype)
        
        snrat = self.calc_sn_ratio(norm_spec)
        
        return snrat
        
        
    def calc_sn_ratio(self, norm_spec=None, nmin = 10):
        if norm_spec is None:
            norm_spec = self.norm_spec
        up_std = get_up_std(norm_spec, 1.0, nmin)
        if up_std is None:
            return 0.0
        return 1./up_std
        


    def get_weighted_mode(self,hn):
        n = len(hn)
        imode = np.argmax(hn)
        imc,nmc = float(imode*hn[imode]), float(hn[imode])
        if imode > 0:
            iml,nml = float((imode-1)*hn[imode-1]), float(hn[imode-1])
        if imode < n-1:
            imr,nmr = float((imode+1)*hn[imode+1]), float(hn[imode+1])

        if imode == 0:
            iwmode = (0.5 + (imr+imc)/(nmr+nmc))
        elif imode == n-1:
            iwmode = (0.5 + (iml+imc)/(nml+nmc))
        else:
            iwmode = (0.5 + (iml+imc+imr)/(nml+nmc+nmr))
        return iwmode


    def normalize_spectrum(self, xspec, spec, binfill=1000, nhist=25, norm_at=None, mode='stretch', clean=True):
        """ Returns normalized spectrum.
        norm_at - forces normalization at given level (default is to normalize at mode value)
        clean   - iterative cleaning of hot pixels (use with caution with high norm_at values) """

        if len(xspec) == 0:
            return np.empty(0)

        if clean:
            spec = spec.copy()
        
        xslen = len(xspec)
        n = xslen // binfill
        bxvals = np.empty(n+2, np.float64) 
        bvals = np.empty(n+2, np.float64)  
        icents = np.ones(n, np.float64)    

        bins = np.linspace(min(xspec), max(xspec), n+1, endpoint=True)  
        ibins = (xslen * np.arange(n+1)) // n     
        
        islices = [slice(ibins[i-1],ibins[i]) for i in range(1,n+1)]    

        bxvals[0]  = bins[0]
        bxvals[-1] = bins[-1]
        for i in range(0,n):
            bxvals[i+1] = (bins[i]+bins[i+1])/2.
            if norm_at is None:
                ss = spec[islices[i]]
                hn,hv = np.histogram(ss,nhist)
                icent = self.get_weighted_mode(hn) / float(nhist)   
                icents[i] = (icent>0.49 and icent) or 0.49    

        if norm_at is None:
            sicents = icents[icents>0.49]
            ncent = (len(sicents)>0 and np.average(sicents)) or 0.5
            print("Normalization at mode level %.3f level"%ncent)
        else:
            ncent = norm_at
            print("Normalization forced to %.3f level."%ncent)

        for i in range(0,n):
            ss = spec[islices[i]]
            len_ss = len(ss)

            if min(ss) < 0.0:
                ss[ss<0.0] = 0.0

            am = np.argsort(ss)[int((len_ss-1)*ncent)]
            ssm = ss[am]
            if i>0 and ssm < 0.5*bvals[i]:
                bvals[i+1] = bvals[i]
            elif ssm<0.01:
                bvals[i+1] = 0.01
            else:
                bvals[i+1] = ssm 

            iclean = clean
            while (iclean):
                up_std = get_up_std(ss-ssm)
                if up_std is None: break
                cond = ss > ssm + 3.*up_std
                if any(cond):
                    ss[cond] = ssm
                else:
                    iclean = False
          
        bvals[0]=bvals[1]
        bvals[-1]=bvals[-2]    

        try:
            spl = splrep(bxvals, bvals)
            if mode == 'stretch':
                norm_spec = np.float32(spec / splev(xspec, spl))
            else:
                norm_spec = np.float32(spec + 1.0 - splev(xspec, spl))

        except TypeError:
            print("Normalization ERROR ! Not normalized spectrum will be used instead.")
            norm_spec = spec
        
        norm_spec[norm_spec<=0.0] = 1.0 
        
        return norm_spec

    def create_norm_spec(self):
        """ Creates normalized spectrum in the class object."""
        binfill, nhist, norm_at, clean = self.main.preferences.get_norm_params()
        norm_spec = self.normalize_spectrum(self.xspec, self.spec, binfill, nhist, norm_at=norm_at, clean=clean)
        return norm_spec

    def create_norm_spec_from(self, xspec, spec):
        """ Creates normalized spectrum in the class object."""
        binfill, nhist, norm_at, clean  = self.main.preferences.get_norm_params()
        norm_spec = self.normalize_spectrum(xspec, spec, binfill, nhist, norm_at=norm_at, clean=clean)
        return norm_spec

    def get_spectrum(self, mode = "normalized"):
        """ Returns a pair: wavelengths array and spectrum array (normalized or not). """
        if mode == "normalized":
            if self.norm_spec is not None:
                return self.xspec, self.norm_spec
        else:
            if self.spec is not None:
                return self.xspec, self.spec
        return None, None

    def get_template(self, num, binw=1):
        """ Returns a pair: wavelengths array and template spectrum array. """
        xtpl, stpl = self.xtpls[num], self.tpl_specs[num]
        if binw>1:
            print(" -> binning enabled:", binw)
            stpl = self.bin_ydata(stpl, binw)
            xtpl = self.create_binned_xdata(binw, self.tpl_f0[num], self.tpl_res[num], len(xtpl))
            
        return xtpl, stpl


    def load_masks(self):
        """ LOAD MASKS """
        self.masks = {"full": [[0000,10000]]}
        mask_fname = None
        for fname in ["masks.list", "masks.ranges", "spectrum.masks"]:
            if os.path.exists(self.main.datadir + "/" + fname):
                mask_fname = fname
                break
        if mask_fname is None:
            print("Warning: file '%s/masks.list' not found."%self.main.datadir)
            return
        with open(self.main.datadir + "/" + mask_fname) as smfile:
            for line in smfile:
                ldata = line.split()
                submasks = "".join(ldata[1:]).split(',')
                mask = [[int(x) for x in sm.split('-')] for sm in submasks]
                self.masks[ldata[0]] = mask
    
    def select_mask(self, name):
        """ SELECT MASK """
        if name not in self.masks:
            print("Warning: mask selected but doesn't exist in the database!")
        self.use_mask = self.masks[name]
        for ma in self.mask_areas:   ma.remove()
        
        fill_mask = []
        for ran in self.use_mask:
            fill_mask.append([ran[0],ran[1],ran[1],ran[0]])
            fill_mask.append([0,0,9,9])
        self.mask_areas = self.ax.fill(*fill_mask,fc=(1.0,1,0.75), ec='k', lw=0.5, zorder=-1)
        self.ax.axis(ymin=self.ymin, ymax=self.ymax, xmin=self.xmin, xmax=self.xmax)
        self.canvas.draw()


    def create_binned_xdata(self, binw, lambda_0, res, dlen):
        lambda_0 += res*(binw/2.-0.5)
        res *= binw
        blen = dlen // binw
        xdata = np.linspace(lambda_0, lambda_0+blen*res, blen, endpoint=False)
        return xdata
    
    def bin_xdata(self, data, binw):
        if binw < 2 or len(data) == 0:
            return data
        res = (data[-1]-data[0])/(len(data)-1)
        shift = res*(binw/2.-0.5)
        nbin = len(data)//binw
        nmiss = len(data)%binw
        eshift = res*nmiss + shift
        bdata = np.linspace(data[0]+shift, data[-1]-eshift, nbin)
        return bdata

    def bin_ydata(self, data, binw):
        if binw < 2:
            return data
        blen = len(data)//binw
        bdata = np.zeros(blen)
        for i in range(binw):
            bdata += data[i:blen*binw:binw]
        print("Ydata:", blen)
        return bdata/float(binw)

    def clear_spectrum(self):
        self.xspec = None
        self.spec_plot.set_xdata([2000,9999])
        self.spec_plot.set_ydata([1,1])
        self.axzoom.axis(xmin=2000, xmax=9999, ymin=0.0, ymax=1.2)
        self.canvas.draw()
    
    def get_mask_range(self):
        return self.use_mask[0][0], self.use_mask[-1][-1] 

    def subtract_templates(self, xspec, spec, tpls_to_sub, res=None, debug=False):
        """ subtract selected templates """

        logxspec = np.log(xspec)
        logxtpl = []
        tspec = []
        
        for i, tpl in enumerate(tpls_to_sub):
            if tpl:
                v_comp = self.main.curve.get_selected_velocity(i, barycorr=True, force=True)
                if v_comp is None:
                    print("Warning! Velocity for %d component not known."%i)
                    continue
                lxt = np.log(self.xtpls[i]) + v_comp/vlight
                logxtpl.append(lxt)
                tspec.append(self.tpl_specs[i])
    
        if len(logxtpl)==0:
            return xspec, spec
    
        if any(map(lambda x: len(x)<10, [logxspec]+logxtpl)):
            print("Warning! Spectrum too short.")
            return xspec, spec

        if res is None:
            res = logxspec[1] - logxspec[0]
        
        xmin, xmax = self.common_range([logxspec]+logxtpl)
        common_xspec = np.arange(xmin, xmax+res, res)

        cond = (logxspec>xmin-3.*res) & (logxspec<xmax+3.*res)
        spec_spline = splrep(logxspec[cond], spec[cond])
        i_spec = splev(common_xspec, spec_spline)

        for lxt, lst in zip(logxtpl, tspec):
            cond = (lxt>xmin-3.*res) & (lxt<xmax+3.*res)
            spec2sub_spline = splrep(lxt[cond], lst[cond])
            i_spec2sub = splev(common_xspec, spec2sub_spline)
            i_spec -= i_spec2sub-1.0

        return np.exp(common_xspec), i_spec # exp is not defined
    
    
    def plot_data(self):
        """ PLOT SPECTRUM DATA """
        if self.xspec is None or len(self.xspec)==0:
            return

        xdata = self.xspec.copy()
        norm_on = self.cb_normalization.isChecked()
        ydata = self.norm_spec if norm_on else self.spec

        tpls_to_subtract = self.sd_dialog.get_to_subtract()
        if self.main.get_sdmode() and any(tpls_to_subtract):
            xdata, ydata = self.subtract_templates(xdata, ydata, tpls_to_subtract, res = self.res/vlight)

        minres = self.main.preferences.get_plot_minres()
        bin = int(minres/self.res)
        if bin>1:
            xdata = self.bin_xdata(xdata, bin) 
            ydata = self.bin_ydata(ydata, bin) 

        if self.main.preferences.get_mask_only_view():
            xmin, xmax = self.get_mask_range()
            cond = (xdata>xmin) & (xdata<xmax)
            if sum(cond)>10:
                xdata = xdata[cond]
                ydata = ydata[cond]

        xmin,xmax = xdata[0],xdata[-1] 
 
        self.spec_plot.set_xdata(xdata)
        self.spec_plot.set_ydata(ydata)

        ymin,ymax = 0, max(ydata)
        ys10 = 0.1 * (ymax - ymin)
        self.ymin = 0
        self.ymax = 1.5 if norm_on else ymax+ys10*0.8
        self.xmin = xmin
        self.xmax = xmax
        if self.cb_update_limits.isChecked():
            self.ax.axis(ymin=self.ymin, ymax=self.ymax, xmin=self.xmin, xmax=self.xmax)
        
        if self.plotting_enabled:
            self.canvas.draw()

    def apply_gaussian_profile(self, tnum, sig, step_size=250):
        spec = self.tpl_specs[tnum]
        res = self.tpl_res[tnum]
        f0 = self.tpl_f0[tnum]
        ntotal = len(self.tpl_specs[tnum])
        newspec = np.ones(ntotal)
        nstep = int(step_size/res) + 1
        for i,j in zip(range(0,ntotal-nstep,nstep), range(nstep,ntotal,nstep)):
            if j >= ntotal-nstep: j = ntotal
            l0 = f0 + res*(i+j)/2
            prof_radius = l0*3.*sig/vlight 
            nran = int(prof_radius/res)+1
            xprof = np.arange(-nran, nran+1)
            rotprof = gauss(xprof,0.,sig)
            nmin = max(0,i-nran)
            nmax = min(ntotal,j+nran)
            imin = max(0,i)
            jmax = min(ntotal,j)
            newspec[imin:jmax] = np.convolve(spec[nmin:nmax], rotprof, 'same')[imin-nmin:jmax-nmin] # convolve is not defined
        return newspec

    def apply_rotational_profile(self, tnum, vsini, step_size=250):
        spec = self.tpl_specs[tnum]
        res = self.tpl_res[tnum]
        f0 = self.tpl_f0[tnum]
        limbdark = self.main.preferences.get_tpl_limbdark()
        newspec = convolve_with_rotational_profile(spec, f0, res, vsini, limbdark, step_size)
        return newspec


    """####################### 
       #    LOAD TEMPLATE    #
       #######################"""        

    def load_template(self, tpl_file, tnum=0, para=False, vsini=0):
        """
        tnum - only 0 for CCF, 0 or 1 for TODCOR
        para - is paratemplate """

        i = tnum
        
        slib_path = self.main.slibdir + "/"+tpl_file
        tpl_path = self.main.tpl_dir + "/"+tpl_file
        
        if para:
            self.tpl_specs[i] = np.fromfile(slib_path,dtype=np.float32)
        elif tpl_file in self.dbase["tpl"]:
            self.tpl_specs[i] = self.dbase["tpl"][tpl_file]
            self.dbind["tpl"].remove(tpl_file)
            self.dbind["tpl"].append(tpl_file)
        elif os.path.exists(tpl_path):
            self.tpl_specs[i] = np.fromfile(tpl_path,dtype=np.float32)
            self.add_to_db("tpl", tpl_file, self.tpl_specs[i])
            if os.path.exists(tpl_path+".info"):
                with open(tpl_path+".info") as finfo:
                    ftxt = finfo.read()
                    print(tpl_file+":")
                    print(ftxt)
        else:
            return -1
        
        tpllen = len(self.tpl_specs[i])
        vals = tpl_file.split('_')
        
        shift = 0
        
        if vals[-1] == "SPECIAL":
            fname = self.main.tpl_dir + "/"+tpl_file.replace("_SPECIAL", ".data")
            with open(fname,'r') as fdata:
                l1 = fdata.readline().split()
                l2 = fdata.readline().split()
                if len(l1)<2 or len(l2)<2:
                    print("Error reading template data file: "+fname)
                    raise IOError
                vals[-1] = l1[1]
                vals.append(l2[1])
        elif para:
            shift = 2
            
        self.tpl_f0[i] = float(vals[shift+1])
        self.tpl_res[i] = float(vals[shift+2])
        
        if para and vsini > 0:
            self.tpl_specs[i]=self.apply_rotational_profile( i, vsini)

        xtpl_ind = (self.tpl_f0[i], self.tpl_res[i], tpllen)
        if xtpl_ind in self.dbase["xtpl"]:
            self.xtpls[i] = self.dbase["xtpl"][xtpl_ind]
            self.dbind["xtpl"].remove(xtpl_ind)
            self.dbind["xtpl"].append(xtpl_ind)
        else:
            self.xtpls[i] = np.linspace(self.tpl_f0[i], self.tpl_f0[i]+tpllen*self.tpl_res[i], \
                                     tpllen, endpoint=False)
            self.add_to_db("xtpl", xtpl_ind, self.xtpls[i])

        if self.xspec is None:
            self.xmin, self.xmax = self.xtpls[i][0], self.xtpls[i][-1]
            if self.xtpls[1-i] is not None and self.cb_show_tpl[1-i].isChecked():
                self.xmin = min(self.xmin, self.xtpls[1-i][0])
                self.xmax = max(self.xmax, self.xtpls[1-i][-1])

        self.plot_template(i)

        return 0

    def get_template_res(self, num):
        return self.tpl_res[num]

    def enable_plotting(self, enabled):
        self.plotting_enabled = enabled

    def plot_template(self, i):
        """ PLOT TEMPLATE """
        self.tpl_plot[i].set_xdata(self.xtpls[i])
        self.tpl_plot[i].set_ydata(self.tpl_specs[i])
        if self.cb_show_tpl[i].isChecked() and self.plotting_enabled:
            #print("PLOT TEMPLATE")
            self.ax.axis(ymin=self.ymin, ymax=self.ymax, xmin=self.xmin, xmax=self.xmax)
            self.canvas.draw()

    def show_mask(self, status=True):
        """ EVENT. Shows mask in RV SPEC window. """
        for ma in self.mask_areas:
            ma.set_visible(status)
        self.canvas.draw()

    def show_spec(self, status=True):
        """ EVENT. Shows spectrum in RV SPEC window. """
        self.spec_plot.set_visible(status)
        #print("PLOTTING SPECTRUM")
        self.canvas.draw()
    
    def switch_norm_spec(self, status=True):
        """ EVENT. Switches between normalized and raw spectra. """
        self.plot_data()


    def template_log_and_mask(self, i_tpl, res, slogxspec, method=None):

        logxtpl = np.log(self.xtpls[i_tpl])
        cond = (logxtpl>slogxspec[0]-777./vlight) & (logxtpl<slogxspec[-1]+777./vlight)
        spl = splrep(logxtpl[cond], self.tpl_specs[i_tpl][cond])
        sltpl = np.float32(splev(slogxspec, spl))

        if method == "TODCOR":
            logxtpl2 = np.log(self.xtpls[1-i_tpl])
            cond = (logxtpl2>slogxspec[0]-777./vlight) & (logxtpl2<slogxspec[-1]+777./vlight)
            spl = splrep(logxtpl2[cond], self.tpl_specs[1-i_tpl][cond])
            sltpl2 = np.float32(splev(slogxspec, spl))
        else:    
            sltpl2 = None

        return sltpl, sltpl2

    def apply_binning(self, xdata, ydata, res):
        minres = self.main.preferences.get_minres()
        bin_size = int(minres/res)
        if bin_size > 1:
            xdata = self.bin_xdata(xdata, bin_size)
            ydata = self.bin_ydata(ydata, bin_size) 
            print("Binning data for calculations:", bin_size)
        return xdata, ydata
    
    
    """##############################
       #     APPLY LOG AND MASK     #
       ##############################"""   
       
    def _get_log_spline_xspec_from_mask(self, xdata, ydata, with_tpl, method, mask_mode, vmargin, res):
        """ internal function used in apply_log_and_mask """
        slogxspec = np.empty(0)

        minx = min(xdata)
        maxx = max(xdata)

        if type(with_tpl) is tuple:
            i_tpl, other_tpl = with_tpl
            minx = max(minx, self.xtpls[i_tpl][0])  
            maxx = min(maxx, self.xtpls[i_tpl][-1]) 
            if method == "TODCOR":
                minx = max(minx, self.xtpls[other_tpl][0])
                maxx = min(maxx, self.xtpls[other_tpl][-1])

        print("Mask:", end=" ")
        if mask_mode == 'minmax':
            rleft = max(minx, self.use_mask[0][0])      
            rright = min(maxx, self.use_mask[-1][-1])   
            vleft = log(rleft)-vmargin/vlight
            vright = log(rright)+vmargin/vlight
            slogxspec = np.arange(vleft, vright, res)
        else:
            for ran in self.use_mask:
                rleft = max(minx,ran[0])
                rright = min(maxx,ran[1])
                print(rleft,"-", rright, end=' ')
                if rleft < rright:
                    slogxspec = np.r_[slogxspec, np.arange(log(rleft), log(rright), res)]
                else:
                    print("(excluded!)", end=' ')
        print()

        return slogxspec
        

    def _sd_mode_tpl_subtraction(self, slspec, slogxspec, other_tpl, write_files=False):
        """ internal function used in apply_log_and_mask """
        vrel_other = self.main.curve.get_selected_velocity(other_tpl, v0_relative=True)
        if vrel_other is not None:
            logxtpl_2sub = np.log(self.xtpls[other_tpl]) + vrel_other/vlight
            tpl_2sub = self.tpl_specs[other_tpl]-1.0
            spl = splrep(logxtpl_2sub, tpl_2sub)
            slspec_2sub = np.float32(splev(slogxspec, spl, ext=1)) 
            slspec -= slspec_2sub
            if write_files:
                np.savetxt("sd_org_2sub",zip(logxtpl_2sub, tpl_2sub)) # savetxt is not defined
                np.savetxt("sd_int_2sub",zip(slogxspec, slspec_2sub))
                #np.savetxt("sd_org",zip(logxspec[np.linalg.cond], self.norm_spec[np.linalg.cond])) # ???????
                np.savetxt("sd_int",zip(slogxspec, slspec+slspec_2sub))
                np.savetxt("sd_sub",zip(slogxspec, slspec))

    def _apply_log_and_mask_for(self, xdata, ydata, resolution, method=None, mask_mode='cut', vmargin=0.0, with_tpl=True, sd_mode=False, vshift=0.0, allow_binning=False, write_files=False):
        """ Applies logarithm and mask for given spectra. """

        if method in ["CCF", "BF"]:
            i_tpl = 1-int(self.main.t1_radio.isChecked())
        else:
            i_tpl = 0 

        other_tpl = 1-i_tpl

        res = resolution/vlight

        if allow_binning:
            xdata, ydata = self.apply_binning(xdata, ydata, resolution)

        logxspec = np.log(xdata)
        logxspec -= vshift/vlight

        slogxspec = self._get_log_spline_xspec_from_mask(xdata, ydata, with_tpl and (i_tpl, other_tpl), method, mask_mode, vmargin, res)

        if len(slogxspec) < 50:
            print("Warning ! Mask is too narrow or outside spectrum.")
            return None, None, False

        cond = (logxspec>slogxspec[0]-777./vlight) & (logxspec<slogxspec[-1]+777./vlight)   
        spl = splrep(logxspec[cond], ydata[cond])
        slspec = np.float32(splev(slogxspec, spl))

        if sd_mode:
            self._sd_mode_tpl_subtraction(slspec, slogxspec, other_tpl, write_files=write_files)
        
        if with_tpl:
            template_specs = self.template_log_and_mask(i_tpl, res, slogxspec, method=method)
        else:
            template_specs = (None, None)
                    
        star_spec = (slogxspec, slspec)       

        return star_spec, template_specs, True
        
        
        
    def apply_log_and_mask(self, resolution, method=None, use_norm=False, mask_mode='cut', vmargin=0.0, with_tpl=True, sd_mode=False, vshift=0.0, allow_binning=False):
        """ Normal use apply_log_and_mask function. """

        xdata = self.xspec
        if use_norm:
            ydata = self.norm_spec
        else:
            ydata = self.spec

        star_spec, templ_specs, is_ok = self._apply_log_and_mask_for(xdata, ydata, resolution, method, mask_mode, vmargin, with_tpl, sd_mode, vshift, allow_binning)
        
        if not is_ok:
            return False
        
        self.masked_tpl_spec, self.masked_tpl2_spec = templ_specs
        self.masked_spec = star_spec[1]     

        return True   



    """#################################################
       #     RADIAL VELOCITY DETERMINATION METHODS     #
       #################################################"""        

                
    def prepare_xdata(self, n, res, v0=0.0):
        """ n must be odd """
        if n%2 == 0:
            raise ValueError("n must be odd, but now is %d"%n)
        dmin, dmax = -n//2+1, n//2+1
        if True:
            xdata = v0 + res * np.arange(dmin, dmax)
        else:
            xdata = v0 + res * np.linspace(dmin, dmax, n, endpoint=False)
        
        return xdata


    """##################### 
       #     CALCULATE     #
       #####################"""        

    def calculate(self, method, resolution, use_norm, sd_mode=False, plot_anal=True):
        """ Prepare data for calculations, do calculations with selected method,
            send results to the Analysis Window. """
        if self.spec is None:
            msgBox = QtWidgets.QMessageBox.warning(self, "No spectrum", \
            "Cannot calculate as there is no spectrum loaded.", \
            QtWidgets.QMessageBox.Ok)
            return
        v0 = self.main.get_v0_from_field()
        barycorr = self.main.curve.get_currpoint_barycorr()
        vshift = v0 - barycorr
        is_ok = self.apply_log_and_mask(resolution,method, use_norm, sd_mode=sd_mode, vshift=vshift, allow_binning=True)
        if not is_ok:
            msgBox = QtWidgets.QMessageBox.warning(self, "Mask error", \
            "Mask is too narrow or outside spectrum.", QtWidgets.QMessageBox.Ok)
            return
         
        mul,add = self.main.preferences.get_vrange_pars()
        ran = self.main.get_model_amplitude()*mul + add
        
        n = int(2.*ran/resolution)

        if method == "CCF":
            xcorr = self.cross_correlation(self.masked_spec, self.masked_tpl_spec)
            xdata = self.prepare_xdata(len(xcorr), resolution, v0)
            self.analysis.plot_ccf(xdata, xcorr, resolution, plot_anal = plot_anal)
        
        elif method == "TODCOR":
            x2dcorr = self.todcorr(n)
            self.analysis.plot_todcor(x2dcorr, resolution, v0, plot_anal = plot_anal)
                
        elif method == "BF":
            nsingular = self.main.preferences.get_num_singular_values()
            n = min(nsingular,max(101,n)) 
            n += (n+1)%2   
            xdata = self.prepare_xdata(n, resolution, v0)
            bfdata = self.broadening_function(n, resolution)
            self.analysis.plot_bf(xdata, bfdata, resolution, plot_anal = plot_anal)
            
        else:
            msgBox = QtWidgets.QMessageBox.warning(self, "Not implemented", \
            "The method is not yet implemented.", QtWidgets.QMessageBox.Ok)
            return 

        if plot_anal:
            self.analysis.show()


    """###################### 
       #     CCF METHOD     #
       ######################"""
   
    def cross_correlation(self, masked_spec, masked_tpl_spec):
        """ Cross-correlation method """
        ccf = irfft(rfft(self.masked_spec) * np.conjugate(rfft(self.masked_tpl_spec))) # conjugate is not defined
        xn = len(ccf)
        ccf = np.r_[ccf[-2000:], ccf[:2001]]
        rsms = rms(self.masked_spec)
        rsmts = rms(self.masked_tpl_spec)
        ccf /= xn*rsms*rsmts
        return 1.-(1.-ccf)/(1.-min(ccf))
   
    def xcorr(self,s1,s2):
       xc = irfft(rfft(s1) * np.conjugate(rfft(s2))) # conjugate is not defined
       return xc


    """####################### 
       #    TODCOR METHOD    #
       #######################"""
    
    def todcorr(self,n):
        """ TODCOR method """
        nran = n//2
        nran2 = nran*2+1
        slen = len(self.masked_spec)
        rms_mspec = rms(self.masked_spec)
        rms_mtpl1 = rms(self.masked_tpl_spec)
        rms_mtpl2 = rms(self.masked_tpl2_spec)
        
        print("TODCOR:", rms_mspec, rms_mtpl1, rms_mtpl2)
        
        a = self.xcorr(self.masked_spec,self.masked_tpl_spec)/(slen*rms_mspec*rms_mtpl1)
        b = self.xcorr(self.masked_spec,self.masked_tpl2_spec)/(slen*rms_mspec*rms_mtpl2)
        d = self.xcorr(self.masked_tpl_spec,self.masked_tpl2_spec)/(slen*rms_mtpl1*rms_mtpl2)
                
        x2dcorr = np.empty((nran2,nran2))
        ii = np.arange(-nran,nran+1,dtype=int)

        lii2 = len(ii)//2
        ai = a[ii]  
        aa = ai**2
        for i in ii:
            di = d[ii-i]
            x2dcorr[nran+i] = ( b[i]**2 + aa - (2*b[i]*ai*di) )/(1-di**2)

        x2dcorr = x2dcorr**0.5
        return x2dcorr 


    """############################## 
       # BROADENING FUNCTION METHOD #
       ##############################"""
    
    def broadening_function(self, n, res):
        """ Broadening Function (written by PK, optimized and corrected by BP) """
        lenspec = len(self.masked_spec)
        lefttrim = n//2
        righttrim = lenspec-n//2 - lenspec%2
        objspec = 1.-self.masked_spec[lefttrim:righttrim]
        lenspec = len(objspec)

        tpl_hash = sum(self.masked_tpl_spec)    
        db_ind = (n, tpl_hash)

        if db_ind in self.dbase["svd"].keys():
            U,singvals,V = self.dbase["svd"][db_ind]
        else:
            print("* Creating new SVD matrices (n = %d)."%n)
            tplspec = self.masked_tpl_spec.copy()
            lentpl = len(tplspec)
            if lentpl%2:
                tplspec = tplspec[:-1]
            tplspec = 1.-tplspec              
            indx,indy = np.meshgrid(np.arange(n,dtype=np.int32),np.arange(lenspec,dtype=np.int32))
            ind = n-1-indx+indy
            des_mat = tplspec[ind]  
            U,singvals,V = svd(des_mat, full_matrices=False, compute_uv=True, overwrite_a=True)              
            self.add_to_db("svd", db_ind, (U,singvals,V) )
        
        bfunc = np.dot(np.transpose(V), np.dot(np.diag(1.0/singvals), np.dot(np.transpose(U),objspec)) )


        return bfunc*10.


    def save_analplot(self, xdata, data, prefix, format = "%.4f"):
        """ NOT IN USE """
        cachepath = self.main.cachedir + "/"+self.main.data["id"]
        ensure_dir_exists(cachepath)
        hjd = self.main.curve.get_current_hjd()
        datapath = cachepath+"_"+prefix+"_"+str(hjd)
        if format == "b":
            data.tofile(datapath+"_bin")
        else:
            data.tofile(datapath+"_txt",sep="\n",format = format)




    """############################## 
       #   SPECTRAL DISENTANGLING   #
       ##############################"""

    def spectra_separation(self, resolution):

        sd_from_tpl = self.main.get_sdmode()
        self.sd_res = resolution/vlight
        
        dirpath, indices, flist, f0, res, baryc, ftype, vmask1, vmask2, v1, v2 = self.main.curve.get_filelist_and_stuff() 
        if dirpath is None:
            return
        self.v0 = self.main.get_v0_from_field()
        niter = self.sd_dialog.get_sd_niter()
        calc_per_iter = 2. if sd_from_tpl else 3.

        self.sd_cost_load = 1.0
        self.sd_cost = 3.0
        self.sd_unit = (self.sd_cost_load + calc_per_iter*niter*self.sd_cost)
        self.sd_pbar_max = len(flist)*self.sd_unit

        sd_xspecdb, sd_specdb, sd_rvs, sd_snrat = self.load_files(dirpath, indices, flist, f0, res, baryc,
                                                        ftype, vmask1&vmask2, self.v0, v1, v2)

        if sd_specdb is None:
            print("No files for analysis.")
            return None
        elif len(sd_specdb)<3:
            print("Only %d file(s) selected for analysis... aborting."%len(sd_specdb))
            return None

        if sd_from_tpl:
            self.specAxy = [self.sd_xtpl, self.sd_tpl-0.5]
            first_calc_mode = self.sd_dialog.get_sd_avgmet()
            print("SD from TPL")
        else:
            i_widest = np.argmax(map(len,sd_xspecdb))
            sd_x = sd_xspecdb[i_widest].copy()
            self.specAxy = [sd_x, 0.5*np.ones(len(sd_x))]
            first_calc_mode = 'maximum'
            print("SD NOT from TPL (maximum function is forced as a method for the first substep)")
            print(np.exp(sd_x[0]), np.exp(sd_x[-1])) # exp not defined

        calc_mode = self.sd_dialog.get_sd_avgmet()

        try:
            for i in range(niter):
                self.specBxy = self.calc_sspectrum(self.specAxy, sd_xspecdb, sd_specdb, sd_rvs, sd_snrat,
                                                   mode='calc_sec',
                                                   calc_mode = first_calc_mode if i==0 else calc_mode,
                                                   write_files=False if i==0 else False)
                self.specAxy = self.calc_sspectrum(self.specBxy, sd_xspecdb, sd_specdb, sd_rvs, sd_snrat,
                                                   calc_mode = calc_mode, debug=False if i==0 else False)

            if not sd_from_tpl:
                self.specBxy = self.calc_sspectrum(self.specAxy, sd_xspecdb, sd_specdb, sd_rvs, sd_snrat,
                                                   mode='calc_sec', calc_mode = calc_mode)
                
        except MemoryError:
            print("Not enough memory to continue operation.")
            return None

        self.sd_pbar_val = self.sd_pbar_max
        self.update_sd_progress_bar()

        self.main.set_sdmode(True)

        print("\nSPECTRAL DISENTANGLING HAS FINISHED")

    
        sdid, sdname = self.create_sd_spectra([self.specAxy, self.specBxy])
        for i in range(2):
            self.main.add_template(sdname[i])
            self.main.set_template(sdid[i], sdname[i], i)
        
        self.sd_dialog.update_res_text()


    def update_sd_progress_bar(self):
        val = 100.*self.sd_pbar_val/self.sd_pbar_max
        if val > 99.5:
            val = 100
        elif val < 0.5:
            val = 0
        self.sd_dialog.update_progress_bar(val)


    def load_files(self, dirpath, indices, filelist, f0, res, baryc, ftype, vmask, v0, v1, v2):
        """ Load all spectra with velocities measured. """
        
        nfiles = len(filelist)
        for key in ["spec", "norm_spec"]:
            if self.dbmax[key]<nfiles:
                self.dbmax[key] = nfiles

        use_norm = self.cb_normalization.isChecked()
        
        sd_xspecdb = []
        sd_specdb = []
        sd_rvs = []
        sd_snrat = []
        self.sd_tpl = None
        
        print("Total # of files:",nfiles)
        for i in range(nfiles):
            self.sd_pbar_val = i*self.sd_cost_load
            self.update_sd_progress_bar()
            
            if vmask[i] != 1:
                print("--- File ignored (a measurement is disabled):", filelist[i])
                self.sd_pbar_max -= self.sd_unit
                continue
                        
            rvtmp = [v1[i], v2[i]]

            print("+++ File loaded:", filelist[i])
            spec, xspec, norm_spec = self._load_any_spec(dirpath+'/'+filelist[i], f0[i], res[i], ftype[i])

            snrat = self.calc_sn_ratio(norm_spec)
            ic = indices[i]
            self.main.curve.snrat[ic] = snrat 

            vshift = v0 - baryc[i]
            lmargin = abs(min(rvtmp)-v0)
            rmargin = abs(max(rvtmp)-v0)
            margin = abs(rvtmp[1]-rvtmp[0])
            print("### v0:", self.v0, " # baryc:", baryc[i], " # vshift=", vshift)
            print("### margin:", margin)

            if use_norm:
                ydata = norm_spec
            else:
                ydata = spec            
            star_spec, templ_specs, is_ok = self._apply_log_and_mask_for(xspec, ydata, res[i], mask_mode='minmax', vmargin=margin, with_tpl=(self.sd_tpl is None), vshift=vshift, allow_binning=True)

            if not is_ok:
                msgBox = QtWidgets.QMessageBox.warning(self, "Mask error", 
                "Mask is too narrow or outside spectrum.", QtWidgets.QMessageBox.Ok)
                # print("--- skipped: Mask is too narrow or outside spectrum.")
                continue

            if self.sd_tpl is None:
                self.sd_xtpl = star_spec[0]     
                self.sd_tpl = templ_specs[0]    

            
            sd_xspecdb += [star_spec[0]]        
            sd_specdb += [star_spec[1]]         
            sd_rvs += [rvtmp]
            sd_snrat += [snrat]

        for key in ["spec", "norm_spec"]:
            self.dbmax[key] = 50

        sd_rvs = np.array(sd_rvs) # array is not defined
        sd_snrat = np.array(sd_snrat)
        sd_snrat[sd_snrat<1.0] = 1.0    
        
        self.main.dataview.update_col('snrat')

        return sd_xspecdb, sd_specdb, sd_rvs, sd_snrat


    def common_range(self, xs):
        x0 = [x[0] for x in xs]
        x1 = [x[-1] for x in xs]
        return max(x0), min(x1)
    

    def subtract_spectrum(self, spec_xy, spec2sub_xy, veldiff=0.0, res=None, debug=False):
        """ veldiff is the v_shift of spec2sub_xy in reference of spec_xy """
        xspec, spec = spec_xy
        xspec2sub_t, spec2sub = spec2sub_xy

        if len(xspec)<10 or len(xspec2sub_t)<10:
            print("Warning! Spectrum too short.")
            return None
        if res is None:
            res = (xspec[-1] - xspec[0])/float(len(xspec)-1)
            print("Resolution is not given and is calculated!")

        xspec2sub = xspec2sub_t.copy() + veldiff/vlight
        
        xmin, xmax = self.common_range([xspec, xspec2sub])
        common_xspec = np.arange(xmin, xmax+self.sd_res, res)

        cond = (xspec>xmin-3.*res) & (xspec<xmax+3.*res)
        spec_spline = splrep(xspec[cond], spec[cond])
        i_spec = splev(common_xspec, spec_spline)

        cond = (xspec2sub>xmin-3.*res) & (xspec2sub<xmax+3.*res)
        spec2sub_spline = splrep(xspec2sub[cond], spec2sub[cond])
        i_spec2sub = splev(common_xspec, spec2sub_spline)
        
        subtracted = i_spec - i_spec2sub

        if debug:
            np.savetxt("sdio",zip(common_xspec,i_spec)) # savetxt is not defined
            np.savetxt("sdis",zip(common_xspec,i_spec2sub))
            np.savetxt("sdid",zip(common_xspec, subtracted))
        
        return common_xspec, subtracted

    def calc_sspectrum(self, xy_specX, sd_xspecdb, sd_specdb, sd_rvs, sd_snrat, mode='calc_primary', calc_mode='median', write_files=False, debug=False):
        """ Calculate single star spectra """

        rvs_2sub = sd_rvs[:,0].copy()  
        rvs_2keep = sd_rvs[:,1].copy()   
        if mode == "calc_primary":
            rvs_2sub, rvs_2keep = rvs_2keep, rvs_2sub      

        all_common_xsp = []
        all_cleared_sp = []

        print("\n *** starting spectral disentangling step in mode: %s\n"%mode)
                
        for i, (x_spec,spec,rv_2keep, rv_2sub) in enumerate(zip(sd_xspecdb, sd_specdb, rvs_2keep, rvs_2sub)):

            self.sd_pbar_val += self.sd_cost
            self.update_sd_progress_bar()

            vshift = rv_2sub - self.v0
            print("::: r2sub=", rv_2sub, " : r2keep=", rv_2keep, " ::: vshift=", vshift)
            common_xspec, spec_subother = self.subtract_spectrum((x_spec,spec), xy_specX, veldiff = vshift, res = self.sd_res, debug=debug if i==0 else False)

            common_xspec -= (rv_2keep-self.v0)/vlight

            all_common_xsp += [common_xspec]
            all_cleared_sp += [spec_subother]
            print("preparing spectrum %d"%i)
            
            if write_files:
                np.savetxt("sdspec_%02d"%i,zip(common_xspec,spec_subother)) # savetxt is not defined
        
        final_spectra = None
        xmin,xmax = self.common_range(all_common_xsp)
        multi_common_xspec = np.arange(xmin, xmax+self.sd_res, self.sd_res)
        for xsp,sp in zip(all_common_xsp, all_cleared_sp):
            multi_spline = splrep(xsp, sp)
            final_subspec = splev(multi_common_xspec, multi_spline)
            if final_spectra is None:
                final_spectra = final_subspec
            else:
                final_spectra = np.c_[final_spectra, final_subspec] # c_ is not defined

        print("\nCombining spectra using", calc_mode)
        calc_fun = {'average': np.mean, 'median': np.median, 'maximum': npmax, 'weighted average': np.average}

        if calc_mode not in calc_fun:
            print("Warning! '%s' mode is not available: 'average' mode is used instead.")
            calc_mode = 'average'

        if 'weighted' in calc_mode:
            avg_spec = calc_fun[calc_mode](final_spectra, weights=sd_snrat, axis=1)
        else:
            avg_spec = calc_fun[calc_mode](final_spectra, axis=1)
        
        med_spec = np.median(avg_spec)
        
        if write_files:
            np.savetxt("sdspec_med",zip(multi_common_xspec, avg_spec)) # savetxt is not defined

        print("\n --> sd step finished\n")

        out_spec = avg_spec-med_spec
        return [multi_common_xspec, out_spec+1.0]


    def create_sd_spectra(self, sdspecs):
        sdid_l, sdname_l = [], []
        for i in range(2):
            name = str(self.main.obj_id_field.text())
            if name.strip() == "":
                name = "nameless"
            if '_' in name: 
                name = name.replace('_', '.')
            sdid = name + "-sd" + str(i+1)

            xspec = np.exp(sdspecs[i][0])  # exp is not defined
            
            
            f0 = xspec[0]
            res = xspec[1] - xspec[0]
            
            new_x = np.arange(f0,xspec[-1],res)
            
            spl = splrep(xspec, sdspecs[i][1])
            new_sp = splev(new_x, spl)

            binfill, nhist, norm_at, clean = self.main.preferences.get_norm_params()
            new_norm = self.normalize_spectrum(new_x, new_sp, binfill, nhist, mode='shift', clean=False)
 
            sdname = sdid + "_SPECIAL"
            sddata = sdid + ".data"
            sdpath =  self.main.tpl_dir + '/' + sdname
            np.float32(new_norm).tofile(sdpath)

            sddpath =  self.main.tpl_dir + '/' + sddata
            with open(sddpath,'w') as fsdd:
                fsdd.write("wavelength0= %f\n"%f0)
                fsdd.write("resolution= %.12f\n"%res)           

            print("Spectrum of component %d written as: %s"%(i+1,sdid))
            
            sdid_l.append(sdid)
            sdname_l.append(sdname)
        return sdid_l, sdname_l


    def get_current_spec(self):
        xdata = self.spec_plot.get_xdata()
        ydata = self.spec_plot.get_ydata()
        return xdata, ydata
        

"""
#########################################
#   Sepctral Disentangling DIALOG BOX   #
#########################################
"""

class SDDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(SDDialog, self).__init__(parent)
        self.spectrum = parent
        self.main = self.spectrum.main

        grid = QtWidgets.QGridLayout()

        self.sd_started = False
        self.sd_time_clicked = 0.0
        
        self.disent_butt = QtWidgets.QPushButton("Spectral\ndisentangling")
        self.disent_butt.setToolTip("Warning! Long execution - double click to run.")
        self.disent_butt.clicked.connect(self.spectra_disentangling)
        self.cb_sdmode = QtWidgets.QCheckBox("SD mode")
        self.cb_sdmode.setStatusTip('Enable Spectral Disentangling mode')
        self.cb_sdmode.stateChanged.connect(self.sd_mode_enabled)

        lab_avgmet = QtWidgets.QLabel("Method:")
        self.combo_avg = QtWidgets.QComboBox()
        self.combo_avg.addItems(["median", "average","maximum",'weighted average'])

        lab_niter = QtWidgets.QLabel("# iter:")
        self.spin_niter = set_QIntSpin(self,1,rmin=1,rmax=5,step=1)
        self.spin_niter.setToolTip("Number of subiterations for each SD step")
        
        # SUBTRACT T1 / SUBTRACT T2
        self.cb_sub_t1 = QtWidgets.QCheckBox("Subtr. T1")
        self.cb_sub_t2 = QtWidgets.QCheckBox("Subtr. T2")
        self.cb_sub_t1.stateChanged.connect(self.subtract_template)
        self.cb_sub_t2.stateChanged.connect(self.subtract_template)

        # CALCULATE ALL V1 and V2
        pb_calc_v1 = QtWidgets.QPushButton("Calc v1")
        pb_calc_v2 = QtWidgets.QPushButton("Calc v2")
        pb_calc_v1.setEnabled(False)
        pb_calc_v2.setEnabled(False)
        pb_calc_v1.clicked.connect(self.main.calc_all_v1)
        pb_calc_v2.clicked.connect(self.main.calc_all_v2)

        # SAVE T1 / SAVE T2
        pb_save_t1 = QtWidgets.QPushButton("Save T1")
        pb_save_t2 = QtWidgets.QPushButton("Save T2")
        pb_save_t1.setToolTip("Saves T1 spectrum (ie. of component 1)\nusing normalization given by L2/L1.")
        pb_save_t2.setToolTip("Saves T2 spectrum (ie. of component 2)\nusing normalization given by L2/L1.")
        pb_save_t1.clicked.connect(self.saveSDtemplate1)
        pb_save_t2.clicked.connect(self.saveSDtemplate2)

        # SAVE CURRENT
        pb_save = QtWidgets.QPushButton("Save\ncurrently\nplotted")
        pb_save.clicked.connect(self.saveSpec_SD)

        # CLOSE DIALOG
        button_box = QtWidgets.QDialogButtonBox()
        close_button = button_box.addButton(QtWidgets.QDialogButtonBox.Close)
        close_button.clicked.connect(self.hide)
        close_button.setAutoDefault(False)

        ### LAYOUT / GROUPS
        group_t1 = QtWidgets.QGroupBox("Component 1:")
        group_t2 = QtWidgets.QGroupBox("Component 2:")


        maxcol = 9
        row = 1

        grid.addWidget(self.disent_butt, row,1,1,3)
        grid.addWidget(lab_avgmet, row,4,1,2)
        grid.addWidget(self.combo_avg, row,6,1,2)        
        grid.addWidget(self.cb_sdmode, row,8,1,2)

        row += 1

        grid.addWidget(lab_niter, row,1)
        grid.addWidget(self.spin_niter, row,2)

        # PROGRESS BAR
        self.sd_progress = QtWidgets.QProgressBar(self)
        self.sd_progress.setToolTip("Shows progess of Spectral Disentangling analysis.")        
        grid.addWidget(self.sd_progress, row,3,1,maxcol-2)

        row += 1

        lab_lumrat_0 = QtWidgets.QLabel("L2/L1 =")
        self.spin_lumrat_b = set_QSpin(self,1.0,rmin=0.001,rmax=99.999,step=0.01,decimals=3)
        lab_lumrat_1 = QtWidgets.QLabel("+")
        self.spin_lumrat_a = set_QSpin(self,0.0,rmin=-9.999,rmax=99.999,step=0.001,decimals=4)
        lab_lumrat_2 = QtWidgets.QLabel("<html>* (&lambda; -</html>")
        self.spin_lumrat_s = set_QSpin(self,4500.,rmin=0.0,rmax=9999.9,step=1.0,decimals=1)
        lab_lumrat_3 = QtWidgets.QLabel(") / 1000")

        lumrat_hbox = QtWidgets.QHBoxLayout()

        lumrat_hbox.addWidget(lab_lumrat_0)
        lumrat_hbox.addWidget(self.spin_lumrat_b)
        lumrat_hbox.addWidget(lab_lumrat_1)
        lumrat_hbox.addWidget(self.spin_lumrat_a)
        lumrat_hbox.addWidget(lab_lumrat_2)
        lumrat_hbox.addWidget(self.spin_lumrat_s)
        lumrat_hbox.addWidget(lab_lumrat_3)
        grid.addLayout(lumrat_hbox, row, 1, 1, maxcol)

        row += 1

        self.lab_sdres = [None, None]
        for i in range(2):
            self.lab_sdres[i] = QtWidgets.QLabel("Res T%d: ?"%(i+1))
            grid.addWidget(self.lab_sdres[i], row,4+i*3,1,3)
        
        lab_binsize = QtWidgets.QLabel("Bin:")
        self.spin_binsize = set_QIntSpin(self,1,rmin=1,rmax=32,step=1)
        self.spin_binsize.valueChanged.connect(self.binning_changed)

        grid.addWidget(lab_binsize, row,1)
        grid.addWidget(self.spin_binsize, row,2)

        row += 1

        gopt_vbox = QtWidgets.QVBoxLayout()
        gopt_vbox.addWidget(pb_save)
        gopt_vbox.addWidget(button_box)
        grid.addLayout(gopt_vbox, row,8,2,2)

        grid.addWidget(group_t1,row,1,1,maxcol-2)
        gt1_hbox = QtWidgets.QHBoxLayout(group_t1)
        gt1_hbox.addWidget(self.cb_sub_t1)
        gt1_hbox.addWidget(pb_calc_v1)
        gt1_hbox.addWidget(pb_save_t1)

        row += 1

        grid.addWidget(group_t2,row,1,1,maxcol-2)
        gt2_hbox = QtWidgets.QHBoxLayout(group_t2)
        gt2_hbox.addWidget(self.cb_sub_t2)
        gt2_hbox.addWidget(pb_calc_v2)
        gt2_hbox.addWidget(pb_save_t2)
        
        self.setLayout(grid)
        self.setWindowTitle("Spectral Disentangling")

    ###############
    ### METHODS ###
    ###############

    def spectra_disentangling(self):
        if self.sd_started:
            return
        if time.time() - self.sd_time_clicked < 1.0:
            self.sd_started = True
            print("... and started.")
            self.spectrum.spectra_separation(self.main.get_resolution())
            self.sd_started = False
        else:
            print("Spectra disentangling enabled... (click again within 1 second to start)")
            self.sd_time_clicked = time.time()
            def gt_fun():
                org_text = self.disent_butt.text()
                for x in np.arange(1.0,0.0,-0.1):
                    self.disent_butt.setText("Click again\nin %.1f second"%x)
                    time.sleep(0.1)
                self.disent_butt.setText(org_text)
            self.genericThread = GenericThread(gt_fun)
            self.genericThread.start()

    def sd_mode_enabled(self, status):
        if not self.sd_started and status:
            msgBox = QtWidgets.QMessageBox.warning(self, "SD mode enabled", "SD mode is now enabled. In this mode if you select T1 and run velocity analysis T2 will be subtracted from each spectrum and only velocity of the first component will be extracted and vice verse if you select T2. Moreover when you start SD analysis, the process will start from a template (in other case it starts from a flat spectrum uniformly equal to 1.", QtWidgets.QMessageBox.Ok)

    def binning_changed(self, value):
        self.update_res_text(value)

    def update_res_text(self, value=None):
        if value is None:
            value = self.spin_binsize.value()
        for i in range(2):
            sdres = self.spectrum.get_template_res(i)
            if sdres is None or sdres == 0.0:
                self.lab_sdres[i].setText("Res T%d: ?"%(i+1))
            else:
                self.lab_sdres[i].setText("Res T%d: %.4f"%(i+1,sdres*value))        

    def update_progress_bar(self, progress=0.0):
        self.sd_progress.setValue(progress)

    def get_sd_avgmet(self):
        return str(self.combo_avg.currentText())

    def get_sd_niter(self):
        return self.spin_niter.value()

    def saveSpec_SD(self):
        xspec, spec = self.spectrum.get_current_spec()
        self.main.save_spectrum(xspec, spec)

    def get_binning(self):
        return self.spin_binsize.value()

    def get_l2l1_params(self):
        return [self.spin_lumrat_b.value(), self.spin_lumrat_a.value(), self.spin_lumrat_s.value()]

    def saveSDtemplate1(self):
        binsize = self.get_binning()
        l2l1params = self.get_l2l1_params()
        self.main.saveTemplate(0, l2l1=l2l1params, binning=binsize)

    def saveSDtemplate2(self):
        binsize = self.get_binning()
        l2l1params = self.get_l2l1_params()
        self.main.saveTemplate(1, l2l1=l2l1params, binning=binsize)

    def subtract_template(self):
        """ EVENT """
        self.spectrum.plot_data()

    def get_to_subtract(self):
        sub1 = self.cb_sub_t1.isChecked()
        sub2 = self.cb_sub_t2.isChecked()
        return sub1, sub2
        
    def show_or_activate(self):
        if self.isHidden():
            self.show()
        else:
            self.activateWindow()
            


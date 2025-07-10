"""
#############################################################################
###########                     RV CURVE                     ################
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


from __future__ import print_function

from PyQt5 import QtWidgets, QtCore, QtGui
from .libcommon import *


class RV_CURVE(QtWidgets.QMainWindow):
    """ Radial Velocity Curve View Window - observations and model """
    def __init__(self, parent):
        super(RV_CURVE, self).__init__(parent)
        self.setParent(parent)
        self.main = parent
        self.spectrum = self.main.spectrum
        self.analysis = self.main.analysis

        self.setWindowTitle("RV Curve")
        self.main_frame = QtWidgets.QWidget()
        self.left_clicked = 0.0
        
        self.Nv = self.main.Nv
        self.init_vals()
        
        self.colors = {0: ['k+', None, None, None], \
                       1: ['bo', 'ro', 'go', 'mo'], \
                       2: ['wo', 'ws', 'wv', 'w^']}
        self.colors_txt = {"active": 1, "inactive": 2, "empty": 0}
        self.tcolors = {}
        for key in self.colors_txt:
            self.tcolors[key] = self.colors[self.colors_txt[key]]
            
        self.mstyle = ['k-', 'k--', 'k:', 'k.']

        self.dpi = 80
        bgcolor = [v / 255.0 for v in self.palette().color(QtGui.QPalette.Window).getRgb()]
        self.fig = Figure((8.0, 5.0), facecolor=bgcolor, edgecolor=bgcolor, linewidth=-1, dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        
        self.axres = self.fig.add_axes([0.09,0.85,0.88,0.13])
        self.base_res, = self.axres.plot([-0.1,1.1],[0,0], 'k:', zorder=-2)
        self.v_res, self.ve_res, self.voff_res = [], [], []
        self.v_plot, self.v_off = [], []
        self.v_model = []
        self.v_sel, self.vres_sel = [], []

        # Fixes error on startup
        self.axres.xaxis.set_tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        self.ax = self.fig.add_axes([0.09,0.09,0.88,0.76], sharex=self.axres)
        self.sel_plot, = self.ax.plot([],[], 'o', c=(1,1,0), mec='k', ms=9, mew=0.5, zorder=10)
        self.rsel_plot, = self.axres.plot([],[], 'o', c=(1,1,0), mec='k', ms=8, mew=0.5, zorder=10)
        self.base_plot, = self.ax.plot([0,1],[0,0], 'k:', zorder=-10)
        self.lc_plot, = self.ax.plot([],[],'.', c=(0.8,0.8,0.8), ms=3, zorder=-9, visible=True)

        marker_size = 4.5
        for i in range(self.Nv):
            self.v_plot   += [self.ax.plot([],[], self.tcolors["active"][i], mec='k', ms=marker_size, mew=0.5, zorder=9-i)[0]]
            self.v_off    += [self.ax.plot([],[], self.tcolors["inactive"][i], mec='k', ms=marker_size, zorder=5-i)[0]]
            self.v_model  += [self.ax.plot([],[], self.mstyle[i], zorder=-2-i)[0]]
            self.v_sel    += [self.ax.plot([],[], self.tcolors["active"][i], ms=marker_size, mew=0.5, mec='k', zorder=19-i)[0]]

        residual_ms = 5
        for i in range(2):
            self.v_res    += [self.axres.plot([],[], self.tcolors["active"][i], # redundant marker definition
                                              ms=residual_ms, zorder=9-i)[0]]
            self.ve_res   += [self.axres.vlines([0],[0],[0], self.tcolors["active"][i][0],
                                                zorder=9-i, alpha=0.5, visible=False)]
            self.voff_res += [self.axres.plot([],[], self.tcolors["active"][i], # redundant marker definition
                                              ms=residual_ms, zorder=5-i)[0]]            
            self.vres_sel += [self.axres.plot([],[], self.tcolors["active"][i], # redundant marker definition
                                              ms=residual_ms, zorder=19-i)[0]]

        self.std_txt = self.ax.text(0.985, 0.035,'x.xx km/s', fontsize=10, ha='right', va='center',
                                        transform=self.ax.transAxes)
        self.std1_txt = self.axres.text(0.985, 0.93,'x.xx km/s', fontsize=10, ha='right', va='top',
                                        transform=self.axres.transAxes, color='blue')
        self.std2_txt = self.axres.text(0.985, 0.02,'x.xx km/s', fontsize=10, ha='right', va='bottom',
                                        transform=self.axres.transAxes, color='red')
            
        xmajor_formatter = ticker.FormatStrFormatter('%.1f')
        self.ax.xaxis.set_major_formatter(xmajor_formatter)
        ymajor_formatter = ticker.FormatStrFormatter('%.4g')
        self.ax.yaxis.set_major_formatter(ymajor_formatter)

        self.ax.set_ylabel("v [km/s]")
        self.ax.set_xlabel("HJD - 2450000")

        self.axres.axis(ymin=-0.81, ymax=0.81, xmin=0, xmax=1)
        self.ax.axis(ymin=-50, ymax=50, xmin=0, xmax=1)
        self.canvas.draw()

        self.cb_update_limits = QtWidgets.QCheckBox("")
        self.cb_update_limits.setChecked(True)
        self.cb_update_limits.setToolTip("Update limits on any change in the plot.")
        
        self.mpl_toolbar = MyNavigationToolbar(self.canvas, self, widgets=[self.cb_update_limits], coordinates=False)
        self.t_label = QtWidgets.QLabel("T=?")
        self.v_label = QtWidgets.QLabel("v=?")
        self.butt_phased = QtWidgets.QPushButton("Phased")
        self.butt_phased.setToolTip("Toggle between raw and phased RV curve.")
        self.cb_show_model = QtWidgets.QCheckBox("Model")
        self.cb_show_model.setToolTip("Show model RV curve.")
        self.cb_show_model.setChecked(True)
        self.butt_phased.clicked.connect(self.sync_phased)
        self.cb_show_model.stateChanged.connect(self.sync_show_model)

        next_butt = QtWidgets.QPushButton(">",self)
        prev_butt = QtWidgets.QPushButton("<",self)
        next_butt.setToolTip("Select next data point.")
        prev_butt.setToolTip("Select previous data point.")
        next_butt.setMaximumWidth(self.main.smallbutton_width)
        prev_butt.setMaximumWidth(self.main.smallbutton_width)
        prev_butt.clicked.connect(self.select_prev_datapoint)
        next_butt.clicked.connect(self.select_next_datapoint)
        prev_butt.setStatusTip('Go to previous spectrum.')
        next_butt.setStatusTip('Go to next spectrum.')
        
        hbox1 = QtWidgets.QHBoxLayout()
        hbox1.addWidget(self.mpl_toolbar)
        hbox1.addWidget(self.cb_show_model)
        hbox1.addWidget(self.butt_phased)
        hbox1.addWidget(prev_butt)
        hbox1.addWidget(next_butt)
        hbox1.addWidget(self.t_label)
        hbox1.addWidget(self.v_label)
        
        vbox = QtWidgets.QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addLayout(hbox1)
        vbox.setStretch(0,1)
                
        self.main_frame.setLayout(vbox)          
        self.setCentralWidget(self.main_frame)   

        self.canvas.mpl_connect('motion_notify_event', self._mmotion)
        self.canvas.mpl_connect('button_press_event', self._mpressed)
        self.canvas.mpl_connect('key_press_event', self._kpressed)
        self.canvas.setFocusPolicy(QtCore.Qt.WheelFocus)

    def init_vals(self):
        """ Set some initial values and create lists and arrays. """
        self.obj_fid = None
        self.curr_dpoint = None
        self.anal_1D_data = []
        self.dirpath = None
        self.v_instr_filter = np.ones(0,bool)
        self.data_sigma = None              
        self.error_scaling = 1.0     
        self.spec_files,self.f0, self.res = [], [], []
        self.baryc, self.instr, self.ftype = [], [], []
        self.snrat, self.dates = [], []
        self.instrument_set = set()
        self.lc_hjd, self.lc_data = np.empty(0), np.empty(0)
        self.hjds, self.xdata = np.empty(0), np.empty(0)
        self.v, self.ve, self.vmask = [], [], []
        for i in range(self.Nv):
            self.v  += [np.empty(0)]
            self.ve += [np.empty(0)]
            self.vmask += [np.empty(0,bool)]


    def data_collection(self, key=None):
        data_aliases = {"hjd": self.hjds, "rv1": self.v[0], "err1": self.ve[0],
                   "rv2": self.v[1], "err2": self.ve[1], "rv3": self.v[2], "err3": self.ve[2],
                   "rv4": self.v[3], "err4": self.ve[3],
                   "mask1": self.vmask[0], "mask2": self.vmask[1], "mask3": self.vmask[2], "mask4": self.vmask[3],
                   "instr": self.instr, 'snrat': self.snrat, "date": self.dates,
                   "w0": self.f0, "wres": self.res, "fname": self.spec_files}
        if key is None:
            return data_aliases.keys()
        else:
            return data_aliases[key]


    def setChecked_ifnotNone(self, obj, state):
        if obj is not None:
            obj.setChecked(state)

    def sync_show_model(self, state):
        self.setChecked_ifnotNone(self.main.cb_show_model, state)
        self.plot_data(stats=False)
        self.canvas.setFocus()

    def sync_phased(self):
        if self.curve_is_phased():
            self.butt_phased.setText("Raw RVC")
            self.setWindowTitle("RV Curve - Raw")
            self.setChecked_ifnotNone(self.main.rv_raw_radio, True)
        else:
            self.butt_phased.setText("Phased")
            self.setWindowTitle("RV Curve - Phased")
            self.setChecked_ifnotNone(self.main.rv_phased_radio, True)
        self.plot_data(stats=False)
        self.canvas.setFocus()

    def show_err(self, status=True):
        """ EVENT. Shows errors in RV Curve window (residuals subplot). """
        for i_ve in self.ve_res:
            i_ve.set_visible(status)
        self.canvas.draw()

    def rescale_errors(self):
        if self.data_sigma is None:
            print(" * Warning! The standard deviation is not yet calculated for this curve. Please fit a model and try again.")
            return 
        errsum = 0.0
        errnum = 0
        for i in range(2):
            vsel = self.ve[i][self.vmask[i]==1]
            errnum += len(vsel)
            errsum += sum(vsel)
        avgerr = errsum/errnum
        self.error_scaling = self.data_sigma/avgerr
        print("Error scaling factor:", self.error_scaling)
        self.main.dataview.scaling_factor_txt.setTitle("esf=%.3g"%self.error_scaling)
        
        self.main.dataview.update_cols("err1", "err2")
        self.plot_data()

    def calculate_ith_snrat(self, ic, update_DV=False):
        path = self.dirpath+'/'+self.spec_files[ic]
        snrat = self.spectrum.get_signal_to_noise_ratio(path,  self.f0[ic], self.res[ic], self.ftype[ic])
        self.snrat[ic] = snrat
        
        if update_DV:
            self.main.dataview.update_col('snrat')
            

    def calculate_snrat(self):
        sortind = self.get_sorted_indices()
        lhjds = len(sortind)
        for i in range(lhjds):
            self.calculate_ith_snrat(sortind[i])
        
        self.main.dataview.update_col('snrat')


    def show_model_is_checked(self):
        return self.cb_show_model.isChecked()

    def curve_is_phased(self):
        return self.butt_phased.text() == "Phased"

    def get_axis_span(self):
        x0,x1,y0,y1 = self.ax.axis()
        xspan, yspan = x1-x0, y1-y0
        return xspan, yspan

    def cache_anal_1D_data(self, xdata, data, method):
        vleft = xdata[0]
        res = xdata[1]-vleft
        self.anal_1D_data[self.curr_dpoint] = (data, vleft, res, method)

    def swap_pointed_velos(self, event):
        dpoint = np.argmin(abs(self.xdata-event.xdata))
        res = self.swap_velos(dpoint)
        if res is not None:
            self.plot_data()

    def get_pointed_point(self, x, y):
        lendata = len(self.xdata)
        xscale = (self.xmax-self.xmin)*self.height()/self.width()
        yscale = self.ymax-self.ymin
        xdist2 = ((self.xdata-x)/xscale)**2
        alldist = np.empty(0)
        for i in range(self.Nv):
            ydist2 = ((self.ydata[i]-y)/yscale)**2
            dist2 = xdist2+ydist2
            dist2[self.vmask[i]==0] = 999.
            alldist = np.r_[alldist, dist2]
        ixiy = np.argmin(alldist)
        dpoint = ixiy%lendata
        iv = ixiy//lendata
        return dpoint, iv

    def delete_pointed(self, event):
        dpoint, iv = self.get_pointed_point(event.xdata, event.ydata)
        res = self.delete_datapoint(dpoint,iv)
        if res is not None:
            self.plot_data()        

    def toggle_pointed(self, event):
        dpoint, iv = self.get_pointed_point(event.xdata, event.ydata)
        res = self.toggle_datapoint(dpoint,iv)
        if res is not None:
            self.plot_data()

    def toggle_primary(self):
        self.toggle_ith(0)
        return True

    def toggle_secondary(self):
        self.toggle_ith(1)
        return True

    def toggle_ith(self, i):
        cond = (self.vmask[i]>0) & self.v_instr_filter
        self.vmask[i][cond] = any(self.vmask[i][cond]==1)+1
        return True

    def _mmotion(self, event):
        if event.inaxes:
            if self.curve_is_phased():
                self.t_label.setText("ph=%.3f"%event.xdata)
            else:
                self.t_label.setText("T=%.2f"%event.xdata)
            self.v_label.setText("v=%.2f"%event.ydata)

    def _mpressed(self, event):
        self.activateWindow()
        if len(self.xdata)==0: return None
        if event.inaxes:
            if event.button==3 or event.button==1 and (time.time()-self.left_clicked)<0.25:
                indices = np.where(self.v_instr_filter)[0]
                filtered_data = self.xdata[self.v_instr_filter]
                fdata_ind = np.argmin(abs(filtered_data-event.xdata))
                dpoint = indices[fdata_ind]
                self.select_datapoint(dpoint)
            elif event.button==2:
                self.toggle_pointed(event)
            elif event.button==1: self.left_clicked=time.time()

    def _kpressed(self,event):
        if len(self.xdata)==0: return None
        if event.key=='s':
            res = self.swap_velos()
        elif event.key=='e' and event.inaxes:
            res = self.swap_pointed_velos(event)
        elif event.key=='S':
            res = self.swap_all_velos()
        elif event.key == 'd':
            res = self.delete_datapoint()
        elif event.key=='f' and event.inaxes:
            res = self.delete_pointed(event)
        elif event.key=='x':
            res = self.toggle_datapoint()
        elif event.key=='c' and event.inaxes:
            res = self.toggle_pointed(event)
        elif event.key in ['b', '1']:
            res = self.toggle_primary()
        elif event.key in ['r', '2']:
            res = self.toggle_secondary()
        else:
            return
        if res is not None:
            self.plot_data()

    def select_datapoint(self, ic, draw=True):
        if len(self.xdata)==0: return None
        if ic is None or ic == self.curr_dpoint: return
        self.curr_dpoint = ic
        self.main.set_instrument(self.instr[ic])
        self.main.set_date(self.hjds[ic])
        self.main.dataview.select_point(ic)
        if draw:
            self.plot_data()
        self.spectrum.load_spec(self.dirpath+'/'+self.spec_files[ic], self.f0[ic], self.res[ic], self.ftype[ic])
        snrat = self.spectrum.calc_sn_ratio()   
        self.snrat[ic] = snrat
        self.main.dataview.update_col('snrat')
        
        if draw and self.anal_1D_data[ic] is not None:
            data, vleft, res, method = self.anal_1D_data[ic]
            self.analysis.plot_cached_data(data, vleft, res, method)
        print("* spectrum HJD = %f, INSTRUMENT: %s "%(self.hjds[ic],self.instr[ic]), end=' ')
        if self.ftype[ic] != "bin":
            print("TYPE:", self.ftype[ic])
        else:
            print()
        if self.spectrum.cb_show_spec.isChecked():
            self.spectrum.show()        
        self.canvas.setFocus()


    def get_filelist_and_stuff(self):
        if len(self.hjds) == 0:
            return [None]*10
        indices = np.where(self.v_instr_filter)[0]
        filt_file = [self.spec_files[i] for i in indices]
        filt_ftype = [self.ftype[i] for i in indices]
        filt_f0 = self.f0[self.v_instr_filter]
        filt_res = self.res[self.v_instr_filter]
        filt_baryc = self.baryc[self.v_instr_filter]
        filt_v0 = self.v[0][self.v_instr_filter]
        filt_v1 = self.v[1][self.v_instr_filter]
        filt_vmask0 = self.vmask[0][self.v_instr_filter]
        filt_vmask1 = self.vmask[1][self.v_instr_filter]
        return self.dirpath, indices, filt_file, filt_f0, filt_res, filt_baryc, filt_ftype, filt_vmask0, filt_vmask1, filt_v0, filt_v1
        

    def calculate_all(self, method, resolution, use_norm, sd_mode=False):
        sortind = self.get_sorted_indices()
        lhjds = len(sortind)
        for i in range(lhjds):
            self.select_datapoint(sortind[i], draw=False)
            plot_analysis = self.main.cb_ca_update.isChecked() or i == lhjds-1
            self.spectrum.calculate(method, resolution, use_norm, sd_mode=sd_mode, plot_anal = plot_analysis)
            QtWidgets.QApplication.processEvents()

    def get_sorted_indices(self, getcurrent=False):
        indices = np.where(self.v_instr_filter)[0]
        filtered_hjds = self.hjds[self.v_instr_filter]
            
        if self.curve_is_phased():
            per,hjd0 = self.main.get_ephem_from_fields()
            vals = (filtered_hjds-hjd0)/per%1.0
            filtered_sortind = np.argsort(vals)
        else:
            filtered_sortind = np.argsort(filtered_hjds)

        sortind = indices[filtered_sortind]

        if getcurrent:
            current = self.curr_dpoint
            if current is not None:
                if current in sortind:
                    current = list(sortind).index(current)
                else:
                    current = None
            return sortind, current
        else:
            return sortind

    def get_datapoint(self, which):
        ndata_in_use = sum(self.v_instr_filter)
        if ndata_in_use == 0: return None
        sortind, sicurr = self.get_sorted_indices(getcurrent=True)
        if which == "prev":
            if sicurr in [None, 0]: return sortind[-1]
            return sortind[sicurr - 1]
        elif which == "next" or True:
            if sicurr in [None, ndata_in_use-1]: return sortind[0]
            return sortind[sicurr + 1]                    

    def get_prev_datapoint(self):
        return self.get_datapoint("prev")

    def get_next_datapoint(self):
        return self.get_datapoint("next")
        
        
    def select_prev_datapoint(self):
        prev_dp = self.get_prev_datapoint()
        if prev_dp is not None:
            self.select_datapoint(prev_dp)

    def select_next_datapoint(self):
        next_dp = self.get_next_datapoint()
        if next_dp is not None:
            self.select_datapoint(next_dp)

    def get_current_hjd(self):
        return self.hjds[self.curr_dpoint]

    def clean_files(self):
        self.init_vals()
        self.clear_selected()

    def update_instruments(self):
        avg_res, avg_f0, avg_n = {}, {}, {}
        for f0,res,instr in zip(self.f0,self.res, self.instr):
            if instr in avg_res:
                avg_res[instr] += res
                avg_f0[instr] += f0
                avg_n[instr] += 1
            else:
                avg_res[instr] = res
                avg_f0[instr] = f0
                avg_n[instr] = 1
        self.instrument_set = avg_res.keys()

        for instr in self.instrument_set:
            avg_res[instr] /= avg_n[instr]
            avg_f0[instr]  /= avg_n[instr]

        self.main.instruments.update_table(self.instrument_set, avg_res, avg_f0, avg_n)

    def new_object(self):
        self.clean_files()
        self.spectrum.new_object_clear()
        self.update_instruments()
        self.plot_data()

    def get_files(self, dirpath, obj_fid):
        """ GET FILES - reads all spectrum files from the object directory
            and interprets additional information from file names. """
        self.clean_files()

        self.dirpath = str(dirpath)
        if not self.dirpath.startswith('/'):
            if self.dirpath[0] == '.':
                self.dirpath = self.main.specdir + self.dirpath[1:] 
            else:
                self.dirpath = self.main.specdir + '/' + self.dirpath
        if not os.path.isdir(self.dirpath):
            msgBox = QtWidgets.QMessageBox.warning(self, "Data Error", "No directory with data for the selected object is found.", QtWidgets.QMessageBox.Ok)
            return
        
        self.obj_fid = obj_fid
        spfiles = sorted(os.listdir(self.dirpath))
        self.hjds=[]
        for spfile in spfiles:
            if spfile[0] not in '0123456789' or spfile.count('_')<4:
                continue
            
            self.spec_files += [spfile]           
            vals = spfile.split('_')
            self.hjds += [vals[0]]
            self.f0 += [vals[1]]
            self.res += [vals[2]]
            self.baryc += [vals[3]]
            rest = vals[4].split('.')   
            self.instr += [rest[0]]     
            if len(rest)>1: self.ftype+=[rest[1]]   
            else:           self.ftype+=['bin']     
            
        ndata = len(self.hjds)
        self.hjds = np.array([fix_hjd(float(x)) for x in self.hjds])
        self.dates = ["%d.%02d.%02d %02d:%02d"%jd2ymd(jd, hms=True)[:-1] for jd in self.hjds]
        for i in range(self.Nv):
            self.v[i]  = np.zeros(ndata)
            self.ve[i] = np.ones(ndata)*0.2
            self.vmask[i] = np.zeros(ndata,int)

        self.f0 = np.array(self.f0,float)
        self.res = np.array(self.res,float)
        self.baryc = np.array(self.baryc,float)
        self.snrat = np.zeros(ndata, float)
        self.anal_1D_data = [None]*ndata

        self.update_instruments()

        self.load_cached_rv()

        self.main.dataview.fill_table()

        self.load_lc()

        self.plot_data()
        self.show()

    def collect_data(self):
        collection = [self.hjds,self.v,self.ve,self.vmask,self.instr]
        return collection

    def toggle_datapoint(self, idp=None, num='all', plot_data=False):
        if idp is None:
            if self.curr_dpoint is None: return None
            idp = self.curr_dpoint
        
        if idp >= len(self.vmask[0]):
            return
        
        for i in range(self.Nv):
            if num in [i, 'all'] or (i<=2 and num == "both"):
                if self.vmask[i][idp] > 0:
                    self.vmask[i][idp] = 3 - self.vmask[i][idp] 
        self.main.dataview.update_row(idp)
        if plot_data:
            self.plot_data()
        return idp

    def delete_datapoint(self, idp=None, num='all', plot_data=False):
        if idp is None:
            if self.curr_dpoint is None: return None
            idp = self.curr_dpoint
        
        for i in range(self.Nv):
            if num in [i, 'all'] or (i<=2 and num == "both"):
                if num != 'all' and self.vmask[i][idp] == 0 and abs(self.v[i][idp])<999.9:
                    self.vmask[i][idp] = 2 
                else:
                    self.vmask[i][idp] = 0 

        self.main.dataview.update_row(idp)
        if plot_data:
            self.plot_data()
        return idp

    def swap_all_velos(self):
        if len(self.hjds) == 0: return None

        for vlist in [self.v, self.ve, self.vmask]:
            vlist[0], vlist[1] = vlist[1], vlist[0]
        return True

    def swap_velos(self,i=None, plot_data=False):
        if i is None:
            if self.curr_dpoint is None:
                return None
            else:
                i = self.curr_dpoint

        for vlist in [self.v, self.ve, self.vmask]:
            vlist[0][i], vlist[1][i] = vlist[1][i], vlist[0][i]

        self.main.dataview.update_row(i)
        if plot_data:
            self.plot_data()
        return i

    def get_values(self, key, scaled_errors=True):
        """ dictionary-like access to normal class attributes """
        if key[:3] == "err" and scaled_errors:
            scaled_errs={"err1": self.ve[0]*self.error_scaling, "err2": self.ve[1]*self.error_scaling,
                         "err3": self.ve[2]*self.error_scaling, "err4": self.ve[3]*self.error_scaling}
            return scaled_errs[key]
        return self.data_collection(key)

    def smartsave_rvcurve(self, filepath=None, cols=None, ignore=None):
        """ Save RV curve for the current object.
            cols - save only selected columns (default is: hjds+rv1+rv2+instr) """
        colsfmt = {'hj': "%10.5f ", 'rv': " %8.3f", 'er': " %6.3f",
                   'in': "  %s", 'da': "  %s", 'fn': "  %s"}
        dmask = {0: '!', 1: ' ', 2: 'x'}
        if filepath is None:
            if self.obj_fid is None:
                return
            ensure_dir_exists(self.main.rvs_dir)
            filepath = self.main.rvs_dir + "/"+self.obj_fid+".rv"

        if cols is None:
            cols = ['hjd','rv1','err1']
            for i in range(1,self.Nv):
                if any(self.vmask[i]):
                    cols.extend(['rv%d'%(i+1), 'err%d'%(i+1)])
            cols += ['instr']
        
        if self.main.VERB:
            print("Writing RVC with the following columns:", cols)

        with open(filepath,'w') as frvc:
            frvc.write("# "+" ".join(cols)+"\n")
            for j,hjd in enumerate(self.hjds):
                if ignore == '2d' and \
                   self.vmask[0][j] != 1 and self.vmask[1][j] != 1:
                    continue
                line = ""
                for cid in cols:
                    if cid[:2] in colsfmt:
                        cfmt = colsfmt[cid[:2]]
                    else:
                        cfmt = " %8.3f"

                    if cid == 'hjd':
                        line += cfmt%(hjd - 2450000.)
                    else:
                        line += cfmt%self.get_values(cid)[j]
                        if cid[:2] == "rv":
                            i = int(cid[2:])-1
                            line += dmask[self.vmask[i][j]]
                line += "\n"
                frvc.write(line)


    def save_rv_curve(self, filepath=None, ignore=None):
        """ Save RV curve for the current object. """
        if filepath is None:
            if self.obj_fid is None:
                return
            ensure_dir_exists(self.main.rvs_dir)
            filepath = self.main.rvs_dir + "/"+self.obj_fid+".rv"
        dmask = {0: '!', 1: ' ', 2: 'x'}
        nwrite = 2
        for i in range(2,self.Nv):
            if any(self.vmask[i]): nwrite = i+1
            
        with open(filepath,'w') as frvc:
            frvc.write("# version 1.03\n")
            for j,hjd in enumerate(self.hjds):
                if ignore == '2d' and \
                   self.vmask[0][j] != 1 and self.vmask[1][j] != 1:
                    continue
                line = "%10.5f "%(hjd - 2450000)
                for i in range(nwrite):
                    line += " %8.3f%1s"%(self.v[i][j],
                                          dmask[self.vmask[i][j]])
                    line += " %6.3f "%self.ve[i][j]
                line += " %s\n"%self.instr[j]
                frvc.write(line)


    def load_cached_rv(self, errors=False):
        """ Load RV curve that was already calculated and saved for the current object.
        Starting from version 1.03 errors are written and read. """
        filepath = self.main.rvs_dir + "/"+self.obj_fid+".rv"
        if os.path.exists(filepath):
            fver = 0.0 
            with open(filepath) as frv:
                c5hjds = np.array([float("%.5f"%x) for x in self.hjds])
                for line in frv:
                    if line[0] == '#':
                        fver = float(line.split()[-1])
                        if fver >= 1.03:
                            errors = True
                        continue
                    ldata = line.split()
                    ihjd = float(ldata[0])+2450000.
                    ldata = ldata[1:]   
                    
                    if ldata[-1][0].isalpha():
                        iinstr = ldata[-1]
                        ldata = ldata[:-1]
                    else:
                        iinstr = None
                    
                    nlvel = len(ldata)
                    mul = 1
                    if errors:
                        mul = 2
                        nlvel //= mul
                        ive = np.empty(self.Nv)
                    iv = np.empty(self.Nv)     
                    ivm = np.ones(self.Nv, int)
                    
                    for i in range(self.Nv):
                        if i>=nlvel:
                            iv[i], ivm[i] = 0.0, 0
                            if errors:
                               ive[i] = 0.0 
                        else:
                            vi_mask = ldata[mul*i][-1]
                            if vi_mask in ['!', 'x']:
                                ivm[i] = 0 if vi_mask=='!' else 2
                                iv[i] = float(ldata[mul*i][:-1])
                            else:
                                iv[i] = float(ldata[mul*i])
                            if errors:
                                ive[i] = float(ldata[mul*i+1])
                    
                    
                    it = np.where(c5hjds==ihjd)[0]
                    if len(it)>0:
                        instrums = [self.instr[k] for k in it]
                        if iinstr is None:      
                            l = 0
                        elif iinstr in instrums:
                            l = instrums.index(iinstr)
                        else:
                            continue
                        j = it[l]
                        self.hjds[j] = ihjd
                        for i in range(self.Nv):
                            self.v[i][j] = iv[i]
                            self.vmask[i][j] = ivm[i]
                            if errors:
                                self.ve[i][j] = ive[i]


    def load_lc(self, lc_path=None):
        if lc_path is None:
            lc_path = self.main.lc_dir + "/"+self.obj_fid+".lc"
        if os.path.exists(lc_path):
            self.lc_hjd, self.lc_data = np.loadtxt(lc_path, usecols=(0,1), unpack=True)
            corr_hjd_list = list(map(lambda x: fix_hjd(float(x)), self.lc_hjd))
            self.lc_hjd = np.array(corr_hjd_list)
        else:
            self.lc_hjd, self.lc_data = np.empty(0), np.empty(0)
                

    def show_light_curve(self, state, plot_data=False):
        self.lc_plot.set_visible(state)
        if plot_data:
            self.plot_data(stats=False)

    def get_currpoint_barycorr(self):
        return self.baryc[self.curr_dpoint]

    def get_selected_velocity(self, num, v0_relative=False, barycorr=False, force=False):
        val = None
        if self.vmask[num][self.curr_dpoint]==1 or force and self.vmask[num][self.curr_dpoint]==2:
            val = self.v[num][self.curr_dpoint]
            if v0_relative:
                val -= self.main.get_v0_from_field()
            if barycorr:
                val -= self.get_currpoint_barycorr()
        return val

    def get_selected_velocities(self):
        velos = [self.get_selected_velocity(0), self.get_selected_velocity(1)]
        return velos

    def get_velocity(self, i, num):
        val = None
        if self.vmask[num][i]==1:
            val = self.v[num][i]
        return val

    def get_velocities(self, i):
        velos = [self.get_velocity(i, 0), self.get_velocity(i, 1)]
        return velos


    def get_vel_data(self, ipoint=None, num=None):
        """ num is a component number, ie. num=0 for primary and 1 for secondary
            ipoint tells for which observation point we want data, if None - for current """
        if ipoint is None:
            ipoint = self.curr_dpoint
        v = [self.v[0][ipoint], self.v[1][ipoint]]
        ve = [self.ve[0][ipoint], self.ve[1][ipoint]]
        vmask = [self.vmask[0][ipoint], self.vmask[1][ipoint]]

        if num is not None:
            return v[num], ve[num], vmask[num]
        return v, ve, vmask


    def get_data_for_id(self, key, ipoint=None):
        if ipoint is not None and ipoint>len(self.hjds):
            return None


        if ipoint is None:
            if key not in self.data_collection():
                return None
            return self.get_values(key)
        else:
            if key not in self.data_collection():
                return None
            return self.get_values(key)[ipoint]

    def get_data_for_ids(self, keys, ipoint=None):
        if ipoint is not None and ipoint>len(self.hjds):
            return None
        data = []
        for key in keys:
            data += [self.get_data_for_id(key, ipoint)]
        return data

    def set_selected_velocity(self, num, val, val_e=None, keep=False, draw=True):
        vX = self.v[num]
        veX = self.ve[num]
        vXmask = self.vmask[num]

        if val is None:
            if not keep:
                vX[self.curr_dpoint] = 0.0
                vXmask[self.curr_dpoint] = 0
        elif type(val) is str and val == 'keep':
            pass
        else:
            vX[self.curr_dpoint] = val 
            if val_e is not None and val_e>0.0009:
                veX[self.curr_dpoint] = val_e
            if vXmask[self.curr_dpoint] == 0:
                vXmask[self.curr_dpoint] = 1
        
        if num<2:
            self.main.dataview.update_point(self.curr_dpoint, num)
        
        if draw:
            self.plot_data()

    def set_selected_velocities(self, velos, verrors=None, draw=True):
        if verrors is None:
            verrors = [None]*len(velos)
        for i, (v, ve) in enumerate(zip(velos,verrors)):
            if v is not None:
                self.set_selected_velocity(i, v, ve, draw=False)
        if draw:
            self.plot_data()

    def clear_selected(self, draw=False):
        self.sel_plot.set_xdata([])
        self.sel_plot.set_ydata([])
        self.rsel_plot.set_xdata([])
        self.rsel_plot.set_ydata([])
        if draw:
            self.canvas.draw() 

    def save_rvc_interactive(self, cols=None):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save RV curve as...', self.main.datadir)[0]
        if not filepath:
            return
        self.smartsave_rvcurve(filepath, cols=cols, ignore='2d')
        

    def save_velocity_data(self,filepath, hjds, vel_cols, write_header=False):
        with open(filepath,'w') as frvc:
            for j,hjd in enumerate(hjds):
                line = "%10.5f "%hjd
                for vc in vel_cols:
                    line += " %8.3f"%vc[j]
                frvc.write(line+"\n")

    def save_residual_data(self, filepath):
        per,hjd0,ecc,aop,v0,k1,k2=self.main.get_pars_from_fields(hjdmode="short")

        xdata = self.get_xdata_from_hjds()
        rvs = self.get_model(self.get_xdata_from_hjds(), mode='all')
        v_instr = self.v_instr_filter

        v_on = [None]*self.Nv
        ydata = [None]*self.Nv
        for i in range(2):
            v_on[i] = self.vmask[i]==1
            #ydata[i] = self.v[i] - rvs[i]
            ydata[i] = self.v[i] + self.v_shifts - rvs[i] # TEST

        iv_on = v_on[0] & v_on[1] & v_instr
        self.save_velocity_data(filepath, xdata[iv_on], [ydata[0][iv_on], ydata[1][iv_on]])

    def save_residual(self):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save residual RVC as...', self.main.datadir)[0]
        if not filepath:
            return  
        self.save_residual_data(filepath)

    def save_orbital_data(self, filepath, use_mask=True):
        per,hjd0,ecc,aop,v0,k1,k2=self.main.get_pars_from_fields(hjdmode="short")

        xdata = self.get_xdata_from_hjds()

        rvs = self.get_model(self.get_xdata_from_hjds(), mode="additional")

        v_instr = self.v_instr_filter

        v_on = [None]*self.Nv
        ydata = [None]*self.Nv

        for i in range(2):
            v_on[i] = self.vmask[i]==1
            ydata[i] = self.v[i] - rvs[i] # ADD V_SHIFTS

        if use_mask:
            iv_on = v_on[0] & v_on[1] & v_instr
        else:
            iv_on = v_instr
            for i in range(2):
                ydata[i][~v_on[i]] = np.nan
                
        self.save_velocity_data(filepath, xdata[iv_on], [ydata[0][iv_on], ydata[1][iv_on]])

    def save_orbital(self):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save orbital RVC as...', self.main.datadir)[0]
        if not filepath:
            return  
        self.save_orbital_data(filepath, use_mask=False)

    def save_3rd_body_data(self, filepath):
        per,hjd0,ecc,aop,v0,k1,k2=self.main.get_pars_from_fields(hjdmode="short")

        xdata = self.get_xdata_from_hjds()
        rvs = self.get_model(self.get_xdata_from_hjds(), exclude='3rd_body')

        v_instr = self.v_instr_filter
        v_on = [None]*self.Nv
        ydata = [None]*self.Nv

        for i in range(2):
            v_on[i] = self.vmask[i]==1
            ydata[i] = self.v[i] - rvs[i] # ADD V_SHIFTS

        iv_on = v_on[0] & v_on[1] & v_instr
        self.save_velocity_data(filepath, xdata[iv_on], [ydata[0][iv_on], ydata[1][iv_on]])

    def save_3rd_body(self):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save 3rd body RVC as...', self.main.datadir)[0]
        if not filepath:
            return  
        self.save_3rd_body_data(filepath)

    def save_puls_data(self, filepath, num=0):
        per,hjd0,ecc,aop,v0,k1,k2=self.main.get_pars_from_fields(hjdmode="short")

        xdata = self.get_xdata_from_hjds()
        rvs = self.get_model(self.get_xdata_from_hjds(), exclude='puls%d'%(num+1))

        v_instr = self.v_instr_filter

        iv_on = (self.vmask[num]==1) & v_instr
        ydata = self.v[num] - rvs[num] # ADD V_SHIFTS

        self.save_velocity_data(filepath, xdata[iv_on], [ydata[iv_on]])

    def save_puls1(self):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save pulsational (A) RVC as...', self.main.datadir)[0]
        if not filepath:
            return  
        self.save_puls_data(filepath, num=0)

    def save_puls2(self):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save pulsational (B) RVC as...', self.main.datadir)[0]
        if not filepath:
            return  
        self.save_puls_data(filepath, num=1)



    def multivelo_eval(self, x, pars, funcs_and_stuff, vmask=None):
        """ x data may be a single vector or a list of arrays separate for each velocity """

        n_pars = len(pars)

        if type(x) is not list:
            x = [x]*self.Nv

        velos = []        
        for i, ith_velo_stuff in enumerate(funcs_and_stuff):
            if ith_velo_stuff is not None:
                v = np.zeros_like(x[i])
                for fun,param_window in ith_velo_stuff:
                    param_window = np.r_[param_window, np.zeros(n_pars-len(param_window), bool)]
                    v += fun(x[i], *pars[param_window])
                velos += [v]
        return velos

    def mperrfun_multivelo(self, pars, mx, my, my_e, funcs_and_stuff, fjac=None):
        funvals = self.multivelo_eval(mx, pars, funcs_and_stuff)
        dev = np.empty(0)
        for y, y_e, funval in zip(my, my_e, funvals):
            idev = (y - funval)/y_e
            dev = np.r_[dev,idev]
        return [0, dev] 

    def mpfit_multivelo(self, data, funcs_and_stuff, pars_and_stuff):
        """ 'my' is a list of velocity-vectors, one for each component """
        mx, my, my_e = data
        pars, fixpars, limited, limits = pars_and_stuff
        n = len(pars)
          
        pari = [{'value':pars[i], 'fixed':fixpars[i], 'limited': limited[i], \
                 'limits': limits[i]} for i in range(n)]
        fa = {'mx':mx, 'my':my, 'my_e':my_e, 'funcs_and_stuff': funcs_and_stuff}
        m = mpfit( self.mperrfun_multivelo,  functkw=fa, parinfo=pari, iterfunct=None)
        if (m.status <= 0): 
            print('error message = ', m.errmsg)
        return m.params, m.perror

    def prepare_pars_and_funcs_and_stuff(self, for_fit=False, mode="all", exclude=[]):
        """ prepare sets of functions, parameters, etc. for given variabilities """
        if type(exclude) not in [list, tuple, set]:
            exclude = [exclude]

        model_mods_org, model_opts = self.main.fitdial.get_rvmods_and_opts()        
        model_mods = set(model_mods_org).difference(exclude)
        
        ith_velo_stuff = [[] for i in range(self.Nv)]

        limited, limits = [], []
        n_params = 0
        pars = np.empty(n_params)
        fixpars = np.empty(n_params)
        porder_list = []

        hjdmode = "short"

        if not mode == "additional":
            params, param_mask = self.main.get_pardict_from_fields(hjdmode=hjdmode, withmask=True)
            
            porder = ["per","t0","ecc","aop","v0","k1","k2"]
            pars = np.array([params[k] for k in porder])
            fixpars = np.array([not param_mask[k] for k in porder], int)
            plim, plimited = self.main.get_parlims_for_pars(porder)
            limits.extend(plim)
            limited.extend(plimited)
            porder_list.append(["basic",porder])

            fun = fit_orbital_1
            param_window = np.array([1,1,1,1,1,1,0], dtype=bool)
            funcs_and_stuff = fun, param_window
            ith_velo_stuff[0] = [funcs_and_stuff]

            fun = fit_orbital_2
            param_window = np.array([1,1,1,1,1,0,1], dtype=bool)
            funcs_and_stuff = fun, param_window
            ith_velo_stuff[1] = [funcs_and_stuff]            

            n_params = len(pars)

        if not mode == "basic" and "3rd_body" in model_mods:
            params, param_mask = self.main.fitdial.get_parvalues_for("3rd_body", withmask=True)
            porder = ['per', 't0', 'ecc', 'aop', 'k12', 'k3']
            pars = np.r_[pars, np.array([params[k] for k in porder])]
            fixpars = np.r_[fixpars, np.array([not param_mask[k] for k in porder], int)]
            plim, plimited = self.main.fitdial.get_parlims_for_pars(porder)
            limits.extend(plim)
            limited.extend(plimited)
            porder_list.append(["3rd_body",porder])
            param_window = np.r_[np.zeros(n_params,bool), np.array([1,1,1,1,0,1], dtype=bool)]
            ith_velo_stuff[2] = [[fit_thirdbody_B, param_window]]
            if model_opts['tb_mode'] == "12+3":
                param_window = np.r_[np.zeros(n_params,bool), np.array([1,1,1,1,1,0], dtype=bool)]
                ith_velo_stuff[0].append([fit_thirdbody_A, param_window])
                ith_velo_stuff[1].append([fit_thirdbody_A, param_window])
            elif model_opts['tb_mode'] == "1+3":
                param_window = np.r_[np.zeros(n_params,bool), np.array([1,1,1,1,1,0], dtype=bool)]
                ith_velo_stuff[0].append([fit_thirdbody_A, param_window])
            elif model_opts['tb_mode'] == "2+3":
                param_window = np.r_[np.zeros(n_params,bool), np.array([1,1,1,1,1,0], dtype=bool)]
                ith_velo_stuff[1].append([fit_thirdbody_A, param_window])
            else:
                print("Error: option '%s' for 'tb_mode' not implemented."%model_opts['tb_mode'])
                raise ValueError
        
            n_params = len(pars)        

        if not mode == "basic":
            for i,rvmod in enumerate(["puls1", "puls2"]):
                if rvmod in model_mods:
                    params, param_mask = self.main.fitdial.get_parvalues_for(rvmod, withmask=True)
                    nharm = self.main.fitdial.get_puls_nharm(i)
                    if nharm>0:
                        pw_enabled = 2+nharm*2
                        porder = ["t0", "per", "a0", "a1"]
                    else:
                        pw_enabled = 3
                        porder = ["t0", "per", "a0"]
                        param_mask["t0"] = 0
                    
                    for j in range(1, nharm):
                        for c in ["a","b"]:
                            key = "%s%d"%(c,j+1)
                            porder.append(key)
                    
                    pars = np.r_[pars, np.array([params[k] for k in porder])]
                    fixpars = np.r_[fixpars, np.array([not param_mask[k] for k in porder], int)]
                    plim, plimited = self.main.fitdial.get_parlims_for_pars(porder)
                    limits.extend(plim)
                    limited.extend(plimited)
                    porder_list.append([rvmod, porder])
                    param_window = np.r_[np.zeros(n_params,bool), np.ones(pw_enabled, dtype=bool)]
                    ith_velo_stuff[i].append([fit_starpuls, param_window])
                    n_params = len(pars)

                
        if for_fit:
            pars_and_stuff = [pars, fixpars, limited, limits]
            return pars_and_stuff, ith_velo_stuff, porder_list
        else:
            pars_and_stuff = pars
            return pars_and_stuff, ith_velo_stuff

    def get_model_for_x(self, x):
        """ NOT IN USE """
        pars, funcs_and_stuff = self.prepare_pars_and_funcs_and_stuff()
        ys = self.multivelo_eval(x*np.ones(1), pars, funcs_and_stuff)
        y = np.array([iy[0] for iy in ys])
        return y

    def get_model_for_current_hjd(self):
        """ get model values for the current observation point """
        pars, funcs_and_stuff = self.prepare_pars_and_funcs_and_stuff()
        hjd = self.get_current_hjd()-2450000
        ys = self.multivelo_eval(hjd*np.ones(1), pars, funcs_and_stuff)
        y = np.array([iy[0] for iy in ys])
        return y


    def get_model(self, model_x, mode="all", exclude=[]):
        """ Get parameteres and functions for  """
        if mode == "basic":
            pars, funcs_and_stuff = self.prepare_pars_and_funcs_and_stuff(mode="basic", exclude=exclude)
        elif mode == "additional":
            pars, funcs_and_stuff = self.prepare_pars_and_funcs_and_stuff(mode="additional", exclude=exclude)
        else: 
            pars, funcs_and_stuff = self.prepare_pars_and_funcs_and_stuff(exclude=exclude)
        model_y = self.multivelo_eval(model_x, pars, funcs_and_stuff)
        
        return model_y


    def fit_multi_model(self):
        if len(self.v[0]) == 0:
            return None
        pars_and_stuff, funcs_and_stuff, porder_list = self.prepare_pars_and_funcs_and_stuff(for_fit=True)
        fixpars = pars_and_stuff[1]
        if sum(fixpars) == len(fixpars):
            print("No parameters to fit. Ignoring fitting process.")
            return


        autoswap, autodisable = self.main.preferences.get_autoadjust_pars()

        if autoswap or autodisable:
            v0 = self.main.get_v0_from_field()
            d1v1,d2v1 = abs(self.v[0]-self.rvm1), abs(self.v[1]-self.rvm1)
            d1v2,d2v2 = abs(self.v[0]-self.rvm2), abs(self.v[1]-self.rvm2)
            diff = d1v1+d2v2 - (d1v2+d2v1)
            d1v0, d2v0 = self.v[0]-v0, self.v[1]-v0
            for i in range(len(diff)):
                if autodisable and d1v0[i]*d2v0[i]>0:
                    if diff[i]<0: self.vmask[1][i]=2
                    else:      self.vmask[0][i]=2
                if autoswap and diff[i]>0:
                    self.swap_velos(i)

        m, x, y, ye = [None]*self.Nv, [None]*self.Nv, [None]*self.Nv, [None]*self.Nv

        instr_shifts = self.main.instruments.collect_colvalues('shift')
        self.v_weights = np.zeros(len(self.instr), float)
        self.v_shifts = np.zeros(len(self.instr), float)
        for i, vins in enumerate(self.instr):
            self.v_shifts[i] = instr_shifts[vins]/1000.

        mi_sum = 0.0
        for vi in range(self.Nv):
            m[vi] = (self.vmask[vi] == 1) & self.v_instr_filter
            x[vi] = self.hjds[m[vi]]-2450000
            y[vi] = (self.v[vi] + self.v_shifts)[m[vi]]
            ye[vi] = self.ve[vi][m[vi]]*self.error_scaling
            
            mi_sum += sum(m[vi])

        if mi_sum < sum(1.0 - fixpars):
            print("ERROR: number of parameters must not exceed data")
            return None

        data = [x, y, ye]
        
        p, pe = self.mpfit_multivelo(data, funcs_and_stuff, pars_and_stuff)

        if pe is None:
            pe = np.zeros(len(p))

        pars_in_order = porder_list[:]
        i = 0   
                
        for ig, group in enumerate(porder_list):
            print(group[0]+":")
            pars_in_order[ig].append([])
            for par in group[1]:
                print("%7s = %10.7g +- %7.3g"%(par, p[i], pe[i]), "" if fixpars[i] else "*")
                pars_in_order[ig][2].append(p[i])
                i += 1
            
            if group[0] == "basic":
                print("basic (calculated):")
                self.print_basic_calculated(p, pe, uncertainties=True)
                
        print
        return pars_in_order

    def print_basic_calculated(self, pars, epars, uncertainties=False):
        """ Print basic information like masses and semi-major axes """
        
        per,ecc,k1,k2 = pars[np.array([0, 2, 5, 6], dtype=int)]
        
        if uncertainties:
            try:
                from uncertainties import ufloat
                from uncertainties.unumpy import sqrt, sin
            except ImportError:
                uncertainties = False # i dont even know what is happening here
                
        if uncertainties:
            e_per,e_ecc,e_k1,e_k2 = epars[np.array([0, 2, 5, 6], dtype=int)]
            per = ufloat(per, e_per)
            ecc = ufloat(ecc, e_ecc)
            k1 = ufloat(k1, e_k1)
            k2 = ufloat(k2, e_k2)

        par1 = 86400./2./pi/695500. 
        par2 = 2*pi* 6.67384*1.9891e10 
        
        A1 = par1 * k1 * per * sqrt(1-ecc**2)
        A2 = par1 * k2 * per * sqrt(1-ecc**2)
        A = A1 + A2
        M1 = (k1+k2)**2*k2*(per*86400)*(1-ecc**2)**1.5 / par2
        M2 = (k1+k2)**2*k1*(per*86400)*(1-ecc**2)**1.5 / par2
        M = M1+M2
        q = M2/M1
        
        if uncertainties:
            print("  Asin(i) = %15s (%s + %s)"%(A.format('.5g'), A1.format('.5g'), A2.format('.5g')))
            print(" Msin3(i) = %15s (%s + %s)"%(M.format('.5g'),M1.format('.5g'),M2.format('.5g')))
            print("  q = %s"%q.format('.5g'))
        else:
            print("  Asin(i) = %.5g (%.4g + %.4g)"%(A, A1, A2))
            print(" Msin3(i) = %.5g (%.4g + %.4g)"%(M,M1,M2))
            print("  q = %.5g"%q)


    def clear_model(self):
        """ clean all models """
        for i in range(self.Nv):
            self.v_model[i].set_xdata([])
            self.v_model[i].set_ydata([])

    def get_xdata_from_hjds(self):
        return self.hjds - 2450000

    def decrease_npoints(self, x, y, limit=2.0, niter=5):
        for i in range(niter):
            n = len(x)
            if n < 200:
                break
            diff = abs((y[1:-1]-y[:-2])-(y[2:]-y[1:-1]))
            keep = diff > limit             
            keep = np.r_[[True],keep,[True]]   
            all_even = np.c_[np.ones(n//2+1,bool),np.zeros(n//2+1,bool)].ravel()[:n]   
            cond = keep | all_even
            ntrue = sum(cond)
            if ntrue<100 or ntrue>int(0.9*n):
                break
            x = x[cond]
            y = y[cond]
        return x, y

    def plot_data(self, stats=True, raw_model_maxcycles=20, update_limits=False):
        """ PLOT CURVE DATA (RAW AND PHASED), WITH OR WITHOUT MODEL """


        per,hjd0,ecc,aop,v0,k1,k2=self.main.get_pars_from_fields(hjdmode="short")

        if self.curve_is_phased():
            self.ax.set_xlabel("phase (P=%.3g)"%per)
            self.xdata = (self.get_xdata_from_hjds()-hjd0)/per%1.0
            self.xmin, self.xmax = -0.04, 1.04
        else:
            self.ax.set_xlabel("HJD - 2450000")
            self.xdata = self.get_xdata_from_hjds()
            self.xmin, self.xmax = get_axis_lims(self.xdata, minspan=per)

        rawmode_model_is_bad = (not self.curve_is_phased()) and \
                           (self.xmax-self.xmin)/per > raw_model_maxcycles
        showmodel_ok = self.show_model_is_checked() and not rawmode_model_is_bad
        if showmodel_ok:
            if self.curve_is_phased():
                mod_hjds = np.linspace(hjd0,hjd0+per,100,endpoint=True)
                mod_xdata = (mod_hjds-hjd0)/per
                mod_xdata = [mod_xdata] * self.Nv
                rv1,rv2 = self.get_model(mod_hjds, mode="basic")[:2]
            else:
                mxmin = self.xmin 
                mxmax = self.xmax 
                step = per/50.
                mod_xdata = np.arange(mxmin, mxmax+per/100., step)
                mod_xdata = [mod_xdata] * self.Nv

                affected = self.main.fitdial.get_enabled_periods()
                for i in range(2):
                    istep = step
                    for sp in affected.get(i,[]):
                        if sp/10. < istep:
                            istep = max(sp/10., step/100.)
               
                    mod_xdata[i] = np.arange(mxmin, mxmax+per/100., istep)
                
                mod_hjds = mod_xdata
                rvs = self.get_model(mod_hjds)
                rv1 = rvs[0]
                rv2 = rvs[1]
            
            xd,yd = self.decrease_npoints(mod_xdata[0], rv1)
            self.v_model[0].set_xdata(xd)
            self.v_model[0].set_ydata(yd)
            xd,yd = self.decrease_npoints(mod_xdata[1], rv2)
            self.v_model[1].set_xdata(xd)
            self.v_model[1].set_ydata(yd)
            mod_rvlim = np.array([min(rv1), max(rv1), min(rv2), max(rv2)])
        else:
            self.clear_model()
            mod_rvlim = np.ones(1)*v0

        
        rv = self.get_model(self.get_xdata_from_hjds())
        self.rvm1, self.rvm2 = rv[:2]   

        if self.curve_is_phased():
            rv_advars = self.get_model(self.get_xdata_from_hjds(), mode="additional")

        active_instr = self.main.instruments.collect_active()
        instr_shifts = self.main.instruments.collect_colvalues('shift')
        self.v_instr_filter = np.zeros(len(self.instr),bool)
        self.v_shifts = np.zeros(len(self.instr), float)

        for i, vins in enumerate(self.instr):
            self.v_instr_filter[i] = vins in active_instr
            self.v_shifts[i] = instr_shifts[vins]/1000.

        yres_on = [None]*2      

        v_on = [None]*self.Nv
        v_off = [None]*self.Nv
        self.ydata = [None]*self.Nv
        for i in range(self.Nv):
            v_on[i] = self.vmask[i]==1
            self.ydata[i] = self.v[i] + self.v_shifts    
            
            if i<2 and self.curve_is_phased():
                self.ydata[i] -= rv_advars[i]
            
            self.ydata[i][self.vmask[i]==0] = v0
            if i == 0:
                v_off[i] = ~v_on[i]
            else:
                v_off[i] = self.vmask[i]>1      
            if i<2:
                yres = self.v[i] - rv[i] + self.v_shifts    
                on_and_instr = v_on[i] & self.v_instr_filter
                yres_on[i] = yres[on_and_instr]
                xres_on = self.xdata[on_and_instr]
                eres_on = (self.ve[i]*self.error_scaling)[on_and_instr]
                
                self.v_res[i].set_xdata(xres_on)
                self.v_res[i].set_ydata(yres_on[i])
                
                verts = create_verts_from_errs(xres_on, yres_on[i], eres_on)
                self.ve_res[i].set_verts(verts)

        npoints = sum(map(len, yres_on))
        if stats and npoints>1:
            rstd = sqrt( sum( [sum(dsqr) for dsqr in [iy**2 for iy in yres_on]]) / npoints )
            
            y2sum = np.asarray([sum(dsqr) for dsqr in [iy**2 for iy in yres_on]])
            ylens = np.asarray([len(dset) for dset in yres_on])

            istd = np.sqrt(y2sum / ylens)
            
            self.data_sigma = rstd
            print("sig = %.1f m/s (%.1f m/s + %.1f m/s)"%(rstd*1000.,
                                                   istd[0]*1000., istd[1]*1000.))

            if rstd<1.0:
                self.std_txt.set_text("%.0f m/s"%(rstd*1000.))
                self.std1_txt.set_text("%.0f m/s"%(istd[0]*1000.))
                self.std2_txt.set_text("%.0f m/s"%(istd[1]*1000.))
            else:
                self.std_txt.set_text("%.2f km/s"%(rstd))
                self.std1_txt.set_text("%.2f km/s"%(istd[0]))
                self.std2_txt.set_text("%.2f km/s"%(istd[1]))
        
        self.base_res.set_xdata([self.xmin,self.xmax])
        self.base_plot.set_xdata([self.xmin,self.xmax])
        self.base_plot.set_ydata([v0]*2)

        for i in range(self.Nv):
            iv_on = v_on[i] & self.v_instr_filter
            iv_off = v_off[i] & self.v_instr_filter
            self.v_plot[i].set_xdata(self.xdata[iv_on])
            self.v_plot[i].set_ydata(self.ydata[i][iv_on])
            self.v_off[i].set_xdata(self.xdata[iv_off])
            self.v_off[i].set_ydata(self.ydata[i][iv_off])

        if self.curr_dpoint is not None:
            ic = self.curr_dpoint
            sx, sy = [], []     
            rsx, rsy = [], []   
            
            for i in range(self.Nv-1,-1,-1):
                if self.vmask[i][ic] > 0:
                    sx += [self.xdata[ic]]
                    sy += [self.ydata[i][ic]]
                    self.v_sel[i].set_xdata([self.xdata[ic]])
                    self.v_sel[i].set_ydata([self.ydata[i][ic]])
                    colmark = self.colors[self.vmask[i][ic]][i]
                    self.v_sel[i].set_color(colmark[0])
                    self.v_sel[i].set_marker(colmark[1])
                else:
                    self.v_sel[i].set_xdata([])
                    self.v_sel[i].set_ydata([])
            
            for i in range(2-1,-1,-1):
                if self.vmask[i][ic]==1:
                    yres = self.v[i][ic] - rv[i][ic] + self.v_shifts[ic]
                    rsx += [self.xdata[ic]]
                    rsy += [yres]
                    self.vres_sel[i].set_xdata([self.xdata[ic]])
                    self.vres_sel[i].set_ydata([yres])
                else:
                    self.vres_sel[i].set_xdata([])
                    self.vres_sel[i].set_ydata([])                            
            
            if len(sx)==0:
                sx = [self.xdata[ic]]
                sy = [v0]

            self.sel_plot.set_xdata(sx)
            self.sel_plot.set_ydata(sy)
            self.rsel_plot.set_xdata(rsx)
            self.rsel_plot.set_ydata(rsy)
        
        minv, maxv = [], []
        for yd in self.ydata:
            if len(yd) > 0:
                minv.append(min(yd))
                maxv.append(max(yd))
        obs_rvlim = np.array(minv + maxv)
        self.ymin, self.ymax = get_axis_lims(np.r_[obs_rvlim,mod_rvlim], minspan=20.0)

        if len(self.lc_data)>0:
            lc_ph = (self.lc_hjd - 2450000 - hjd0) / per % 1.0
            lc_val = self.lc_data

            scale = -(self.ymax-self.ymin)/2.0
            avg = np.average(lc_val)
            if len(lc_val)>1:
                lc_span = np.ptp(lc_val)
                if lc_span < 0.05: lc_span = 0.05
                lc_val = scale*(lc_val-avg)/lc_span+v0
        
            lc_x, lc_y = np.empty(0), np.empty(0)
            if self.curve_is_phased():
                lc_x = lc_ph
                lc_y = lc_val
            elif not rawmode_model_is_bad:
                nperiods = int((self.xmax-self.xmin)/per)+1
                ph0 = (self.xmin-hjd0)/per
                xobs = self.xmin + (lc_ph-ph0)%1.0*per
                for i in range(nperiods):
                    lc_x = np.r_[lc_x, xobs+i*per]
                    lc_y = np.r_[lc_y, lc_val]

            self.lc_plot.set_xdata(lc_x)
            self.lc_plot.set_ydata(lc_y)
        else:
            self.lc_plot.set_xdata(np.empty(0))
            self.lc_plot.set_ydata(np.empty(0))

        plotres_ran = self.main.preferences.get_res_range()
        if self.cb_update_limits.isChecked():
            self.axres.axis(xmin=self.xmin, xmax=self.xmax, ymin=-plotres_ran, ymax=plotres_ran)
            self.ax.axis(ymin=self.ymin, ymax=self.ymax)

        self.canvas.draw()
        

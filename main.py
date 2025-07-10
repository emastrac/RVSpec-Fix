#! /usr/bin/env python3

#################################################################
#                RV & Spectral Analysis - RaveSpan              #
#                      by Bogumil Pilecki                       #
#                  & Edited by Everett Mastracci                #
version="1.5.0_edit"                                            #
reldate="10 July 2025"                                          #
# Copyright (C) 2010-2023 Bogumil Pilecki                       #
# Send comments and wishes to:                                  #
# E-mail: pilecki AT camk DOT edu DOT pl                        #
#################################################################
#                                                               #
# If you publish results or data obtained using this software   #
# please cite:                                                  #
# Pilecki et al., 2017, ApJ, 842, 110                           #
#                                                               #
# you may also consider citing one of these papers with         #
# the usage examples as well:                                   #
# Pilecki et al. 2013, MNRAS, 436, 953                          #
# Pilecki et al. 2015, ApJ, 806, 29                             #
# Pilecki et al. 2018, ApJ, 862, 43                             #
# Pilecki et al. 2018, ApJ, 868, 30                             #
#                                                               #
#################################################################
#
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
# This file is part of the RaveSpan program.
# RaveSpan program includes the following files:
# * main.py      - main RaveSpan GUI
# * librvc.py    - radial velocity curve plot
# * libspec.py   - spectrum related operations
# * libdata.py   - information about data points
# * libanal.py   - analysis window, velocity measurements
# * libdial.py   - various additional dialog windows
# * libcommpn.py - a bunch of small, useful functions 
# * bin/rvspec       - executable python script to run RaveSpan
# * work/specconv.py - a template program for spectra conversion to RaveSpan format
#
# RaveSpan makes use of mpfit.py program (details inside the file).
#
#################################################################

from .librvc  import *
from .libspec import *
from .libanal import *
from .libdial import *
from .libdata import *

class RV_MENU(QtWidgets.QMainWindow):
    """ Main window """
    def __init__(self, objfile=None, datadir='data', slibdir='data/slib', verbose=False, parent=None):
        QtWidgets.QMainWindow.__init__(self, parent)
        self._initvals(datadir,slibdir,verbose)

        if sys.platform == "darwin":
            self.smallbutton_width = 50
        else:
            self.smallbutton_width = 30

        self.Nv = 4
        self.analysis = SPECANAL(self)
        self.spectrum = SPECTRUM(self)
        self.curve = RV_CURVE(self)
        self.dataview = DATA_VIEW(self)

        self.setWindowTitle('RaveSpan - Main')
        
        screen = QtWidgets.QDesktopWidget().screenGeometry()
        self.gui_menu()
        self.statusBar()

        self.main_widget = QtWidgets.QWidget(self)

        ### MAIN LAYOUT ###
        self.main_layout = QtWidgets.QVBoxLayout(self.main_widget)
        self.gui_object_data()
        self.gui_rv_calculations()
        self.gui_show_options() 
        self.gui_bottom_buttons()
        
        self.preferences.load_preferences()
        
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.adjustSize()
        self.move(self.geometry().left(),self.geometry().top())
        self.curve.move(self.geometry().right()+3,self.geometry().top())
        self.curve.adjustSize()
        self.analysis.move(self.geometry().right()+3,self.curve.geometry().bottom()+67)
        # self.spectrum.move(self.geometry().right()+3,screen.height()-self.spectrum.height()) 
        # Causes window to draw off screen on Windows


    """
    ##############################################################
    ###                    EVENT FUNCTIONS                     ###
    ##############################################################
    """

    def closeEvent(self, event):
        """ Executed when X button in the left/right top window corner is clicked. """
        event.accept()
               

    def _initvals(self, datadir, slibdir, verbose):
        self.VERB = verbose
        self.ca_started = False
        self.ca_time_clicked = 0
        self.filepath = None
        self.data = {}
        self.tpl_files = {}
        self.paratemplates = []
        self.instr_mask_dict = {}
        self.show_rv = False
        self.show_spec = False
        self.plot_on_pchange = True
        self.tpl_activation_disabled = False
        self.mask_activation_disabled = False
        self.datadir = datadir
        self.slibdir = slibdir
        self.tpl_buttons_enabled = True
        self.tpl_dir = self.datadir + '/templates'
        self.obj_dir = self.datadir + '/objects'
        self.specdir = self.datadir + '/specdb'
        self.rvs_dir = self.datadir + '/rvs'
        self.lc_dir  = self.datadir + '/lcs'
        self.cachedir = self.datadir + '/cache'
        self.RVMethods = ["ccf","todcor","bf"]
        self.default_data_values()


    def get_pars_from_fields(self, hjdmode):
        """ OLD, but used in some functions,
            ought to be replaced by get_pardict_from_fields.
            returns:
                [per,hjd0,ecc,aop,v0,k1,k2]
            hjdmodes:
              - short    5421.1234
              - full     2455421.1234 """
        per  = self.obj_spin["per"].value()
        hjd0 = self.obj_spin["t0"].value()
        if hjdmode == "full":
            hjd0 = full_hjd(hjd0)
        ecc  = self.obj_spin["ecc"].value()
        aop  = self.obj_spin["aop"].value()
        v0 = self.obj_spin["v0"].value()
        k1 = self.obj_spin["k1"].value()
        k2 = self.obj_spin["k2"].value()
        return [per, hjd0, ecc, aop, v0, k1, k2]

    def get_ephem_from_fields(self, fixhjd = True):
        per, hjd0 = self.obj_spin["per"].value(), self.obj_spin["t0"].value()
        if fixhjd:
            hjd0 = fix_hjd(hjd0)
        return per, hjd0

    def get_v0_from_field(self):
        return self.obj_spin["v0"].value()


    def get_model_amplitude(self):
        k1,k2 = self.obj_spin["k1"].value(), self.obj_spin["k2"].value()
        ecc, aop = self.obj_spin["ecc"].value(), self.obj_spin["aop"].value()
        amplitude = max(k1,k2) * (1.0 + ecc*cos(aop))
        return amplitude
        

    def get_pardict_from_fields(self, hjdmode, withmask=False):
        """ returns a dictionary of parameters """
        pars = ['per', 't0', 'ecc', 'aop', 'v0', 'k1', 'k2']
        pd = {p: self.obj_spin[p].value() for p in pars}
        if hjdmode == "full":
            pd['t0'] = full_hjd(pd['t0'])

        if withmask:
            pm = {p: self.obj_cb[p].isChecked() for p in pars}
            return pd, pm
        return pd
    
    def get_parlims_for_pars(self, parlist):
        dplims = {"ecc": [0.0,0.99], "aop": [-0.1,2*pi+0.1],
                  "per": [0.001,99999.], 't0': [-49999,49999.],
                  "k1": [0.0, 999.], "k2": [-999.0, 999.], "v0": [-999.,999.]}

        parlims = []
        parlimited = []
        for p in parlist:
            parlims += [dplims[p]]
            parlimited += [[1,1]]
            
        return parlims, parlimited
    

    def set_pars_from_fitting(self, per, hjd0, ecc, aop, v0, k1, k2):
        """ [per,hjd0,ecc,aop,v0,k1,k2] """
        for name,value in [["per", per], ["t0", hjd0], ["ecc", ecc], ["aop", aop], ["v0", v0], ["k1", k1], ["k2", k2]]:
            self.obj_spin[name].blockSignals(True)
            self.obj_spin[name].setValue(value)
            self.obj_spin[name].blockSignals(False)


    def fit_params(self):
        pars_in_order = self.curve.fit_multi_model()
        if pars_in_order is None:   return
        
        self.plot_on_pchange = False
        for group in pars_in_order:
            if group[0] == "basic":
                self.set_pars_from_fitting(*group[2])
            elif group[0] == "3rd_body":
                self.fitdial.set_parvalues_for("3rd_body", dict(zip(*group[1:])))
            elif group[0] in ["puls1", "puls2"]:
                self.fitdial.set_parvalues_for(group[0], dict(zip(*group[1:])))
        self.plot_on_pchange = True        
        self.curve.plot_data()

    def add_masks(self, obj=None, clear=False):
        if obj is None: obj = self.mask_combo
        if clear: obj.clear()
        masks = [*self.spectrum.masks.keys()]
        for mask in masks:
            obj.addItem(mask)
        self.spectrum.select_mask(masks[0])

    def get_paratemplate_index(self, tid):
        for i,ptobj in enumerate(self.paratemplates):
            if tid == ptobj['id']:
                return i
        return None

    def set_paratemplate(self, tid, tname, rotation, tnum=0):
        self.tpl_activation_disabled = True
        pti = self.add_paratemplate(tid, tname, rotation, tnum)
        self.tpl_combos[tnum].setCurrentIndex(pti)
        self.spectrum.load_template(tname, tnum, para=True, vsini=rotation)
        self.tpl_activation_disabled = False

    def add_paratemplate(self, tid, tname, rotation, tnum=0):
        paratpl_obj = {'id': tid, 'fname': tname, 'rotation': rotation}
        if paratpl_obj in self.paratemplates:
            pti = self.tpl_combos[tnum].findText(tid)
            return pti
        self.paratemplates.insert(0, paratpl_obj)
        other_selected = self.tpl_combos[1-tnum].currentText()
        if len(self.paratemplates)>5:
            other_selected == self.paratemplates[5]['id']
            self.paratemplates[4]=self.paratemplates[5]
        self.paratemplates=self.paratemplates[:5]
        for i in range(2):
            ptt = self.tpl_combos[i].currentText()
            self.tpl_combos[i].clear()
            for pto in self.paratemplates:
                self.tpl_combos[i].addItem(pto['id'])
            for key in self.tpl_files:
                self.tpl_combos[i].addItem(key)
            self.tpl_combos[i].insertSeparator(len(self.paratemplates))
            pti = self.tpl_combos[i].findText(ptt)
            self.tpl_combos[i].setCurrentIndex(pti)
        return 0

    def set_mask(self, mask):
        self.mask_activation_disabled = True
        pmi = self.mask_combo.findText(mask)
        self.mask_combo.setCurrentIndex(pmi)
        self.spectrum.select_mask(mask)
        self.mask_activation_disabled = False

    def set_template(self, tid, tname, tnum=0):
        self.tpl_activation_disabled = True
        pti = self.tpl_combos[tnum].findText(tid)
        self.tpl_combos[tnum].setCurrentIndex(pti)
        self.spectrum.load_template(tname, tnum)
        self.tpl_activation_disabled = False


    def add_template(self, tname, tnum=0):
        """ add interpolated template to combo box """
        tdata = tname.split('_')
        if tdata[0] in self.tpl_files:  
            self.spectrum.remove_from_db("tpl", tname)
            return
        self.tpl_files[tdata[0]] = tname
        for i in range(2):
            self.tpl_combos[i].addItem(tdata[0])
        

    def add_templates(self, tnum=0):
        if not os.path.exists(self.tpl_dir):
            if tnum==0:
                print("Warning: directory '%s' doesn't exist."%self.tpl_dir)
            return

        obj = self.tpl_combos[tnum]
        tpls = sorted(os.listdir(self.tpl_dir))
        loaded = False
        for tpl in tpls:
            if tpl.endswith(".info"):
                continue
            tdata = tpl.split('_')
            if len(tdata)>2 or len(tdata)==2 and tdata[-1]=="SPECIAL":
                self.tpl_files[tdata[0]] = str(tpl)
                obj.addItem(tdata[0])
                if not loaded:
                    self.spectrum.load_template(self.tpl_files[tdata[0]], tnum)
                    loaded = True


    def refreshTemplates(self):
        """ EVENT """
        for i in range(2):
            ptt = self.tpl_combos[i].currentText()
            self.tpl_combos[i].clear()
            self.add_templates(i)
            
            pti = self.tpl_combos[i].findText(ptt)
            if pti>=0:
                self.tpl_combos[i].setCurrentIndex(pti)
    

    def show_or_activate(self, window):
        if window.isHidden():
            window.show()
        else:
            window.activateWindow()

    def disable_tpl_buttons(self):
        self.tpl_buttons_enabled = False
        for i in range(2):
            self.tpl_buttons[i].setEnabled(False)
            self.tpl_buttons[i].setToolTip("Sorry, library with template spectra not found !")

    def template_pressed(self, tnum=0):
        self.tdialog.select_template(tnum)
        self.show_or_activate(self.tdialog)

    def template_activated(self, stext, tnum=0):
        if self.tpl_activation_disabled: return
        if stext in self.tpl_files:
            state = self.spectrum.load_template(self.tpl_files[stext],tnum)
            if state<0:
                print("\nError reading template:", stext, end=" ")
                print("=> Template list refresh was forced")
                msgBox = QtWidgets.QMessageBox.warning(self, "Template error", "Error reading template: " + \
                                          stext + "\nTemplate list refresh was forced.", QtWidgets.QMessageBox.Ok)
                self.refreshTemplates()
        else:
            pti = self.get_paratemplate_index(stext)
            if pti is not None:
                self.spectrum.load_template(self.paratemplates[pti]['fname'], tnum, para=True)
        self.sd_dialog.update_res_text()

    def mask_activated(self, mtext):
        if self.mask_activation_disabled: return
        mtext = str(mtext)
        if mtext.strip() != "":
            self.spectrum.select_mask(mtext)
            cur_instr = str(self.instr_name.text()).strip()
            if cur_instr != "":
                self.instr_mask_dict[cur_instr] = mtext

    def resolution_changed(self, value):
        pass

    def method_activated(self, value):
        tpl2_on = (value == "TODCOR")
        if not tpl2_on:
            num = 1-int(self.t1_radio.isChecked())
            self.tdialog.tmpl[num].toggle()
    

    def get_method(self):
        method = str(self.method_combo.currentText())
        return method

    def set_sdmode(self, state):
        self.sd_dialog.cb_sdmode.setChecked(state)

    def get_sdmode(self):
        return self.sd_dialog.cb_sdmode.isChecked()

    def get_resolution(self):
        return self.r_spin.value()

    def get_options_for_vcalc(self):
        resolution = self.r_spin.value()
        method = self.get_method()
        use_norm = self.cb_normalization.isChecked()
        sd_mode = self.get_sdmode()
        return resolution, method, use_norm, sd_mode  

    def calc_all_v1(self):
        pass

    def calc_all_v2(self):
        pass

    def calc_one_pressed(self):
        if self.curve.curr_dpoint is None:
            self.curve.select_next_datapoint()
        resolution, method, use_norm, sd_mode = self.get_options_for_vcalc()
        self.spectrum.calculate(method, resolution, use_norm, sd_mode=sd_mode)

    def calc_all_pressed(self):
        if self.ca_started:
            return
        if time.time() - self.ca_time_clicked < 1.0:
            self.ca_started = True
            self.calc_all_butt.setEnabled(False)
            print("... and started.")
            resolution, method, use_norm, sd_mode = self.get_options_for_vcalc()
            self.curve.calculate_all(method, resolution, use_norm, sd_mode=sd_mode)
            self.calc_all_butt.setEnabled(True)
            self.ca_started = False
        else:
            print("Calcualte All enabled... (click again within 1 second to start)")
            self.ca_time_clicked = time.time()
            def gt_fun():
                for x in np.arange(1.0,0.0,-0.1):
                    self.calc_all_butt.setText("%.1f"%x)
                    time.sleep(0.1)
                self.calc_all_butt.setText("Calc All")
            self.genericThread = GenericThread(gt_fun)
            self.genericThread.start()


    def do_show_rv(self):
        self.show_or_activate(self.curve)

    def do_show_spec(self):
        self.spectrum.show_or_activate()

    def do_show_anal(self):
        self.show_or_activate(self.analysis)

    def do_show_data(self):
        self.show_or_activate(self.dataview)

    def sync_normalization(self, value):
        self.spectrum.cb_normalization.setChecked(value)

    def ca_update_changed(self,state):
        self.preferences.set_ca_update(state)

    def toggle_phased(self):
        if not self.curve.curve_is_phased():
            self.curve.butt_phased.setText("Phased")
            self.curve.plot_data(stats=False)

    def toggle_raw(self):
        if self.curve.curve_is_phased():
            self.curve.butt_phased.setText("Raw RVC")
            self.curve.plot_data(stats=False)

    def q_adjusted(self, val):
        k2v = val * self.obj_spin["k1"].value()
        self.obj_spin["k2"].setValue(k2v)

    def per_changed(self):
        val = self.obj_spin["per"].value()
        self.obj_spin["per"].setSingleStep(10e-6*val**2)
        
    def aop_changed(self, val):
        print("AOP changed")
        if  val < 0.0 or val > dpi:
            self.obj_spin["aop"].blockSignals(True)
            self.obj_spin["aop"].setValue(val%dpi)
            self.obj_spin["aop"].blockSignals(False)
        self.param_changed()
        
    def param_changed(self):
        if self.plot_on_pchange:
            self.curve.plot_data()

    def set_instrument(self, instr):
        self.instr_name.setEnabled(True)
        self.instr_name.setText(instr)
        if self.cb_instr_mask.isChecked():
            print(instr)
            print(self.instr_mask_dict)
            instr = str(instr).strip()
            if instr in self.instr_mask_dict.keys():
                mask = self.instr_mask_dict[instr]
                print(mask + " selected")
                self.set_mask(mask)

    def set_date(self, jd):
        self.date.setEnabled(True)
        year,month,day,h,m,s = jd2ymd(jd, hms=True)
        txt = "%d.%02d.%02d %02d:%02d"%(year,month,day,h,m)
        self.date.setText(txt)
        

    def set_fields(self, data=None):
        """ Set the parameter fields (main and additional)."""
        self.instr_name.setEnabled(False)
        self.instr_name.setText('not set')
        self.date.setText('not set')
        self.date.setEnabled(False)
        self.plot_on_pchange = False
        print("Plotting disabled")
        
        self.fitdial.clear_fields()
        
        if data is None:
            data = self.data

        self.obj_path_field.setText("")

        for key,value in data["main_pars"].items():
            if key in ["per","t0","v0","k1","k2","ecc","aop"]:
                self.obj_spin[key].blockSignals(True)
                self.obj_spin[key].setValue(value)
                self.obj_spin[key].blockSignals(False)
        
        for key,value in data["general"].items():
            if key == 'id':     
                self.obj_id_field.setText(value)
            elif key == 'data':
                self.obj_path_field.setText(value)
            
        self.fitdial.set_fields(data)
        self.param_changed()
        
        self.plot_on_pchange = True


    def id_changed(self,value):
        self.tdialog.set_interp_name(value)


    def default_data_values(self, set_id="new_object", set_datapath=None):
        self.data = {}
        for key in ["general", "main_pars", "rv_method", "3rd_body", "puls1", "puls2"]:
            self.data[key] = {}
            
        if set_id:  self.data["general"]["id"] = set_id
        if set_datapath:  self.data["general"]["data"] = set_datapath
        self.data['main_pars']["per"] = 1.0
        self.data['main_pars']["t0"] = 5000.
        self.data['main_pars']["v0"] = 0.	
        self.data['main_pars']["k1"] = 150.
        self.data['main_pars']["k2"] = 90.
        self.data['main_pars']["ecc"] = 0.
        self.data['main_pars']["aop"] = 0.
        

    def newObject(self):
        self.default_data_values()
        self.setWindowTitle('RaveSpan - Main')
        self.set_fields()
        self.curve.new_object()

        self.spectrum.clear_spectrum()
        self.set_sdmode(False)


    def openObject(self):
        ##PYQT5: getOpenFileName is now the new name for getOpenFileNameAndFilter, which returns a tuple (filename, selected filter)
        filepath = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', self.obj_dir)[0]
        if filepath:
            self.loadObject(filepath)        



    def loadObject(self, filepath=None):
        if not os.path.exists(filepath):
            msgBox = QtWidgets.QMessageBox.warning(self, "File error", "File %s does not exist."%filepath, QtWidgets.QMessageBox.Ok)
            return
        if not filepath.endswith(".obj"):
            msgBox = QtWidgets.QMessageBox.warning(self, "File error", "This is not an object file.", QtWidgets.QMessageBox.Ok)
            return
        obj_fid = os.path.basename(filepath)[:-4]

        self.default_data_values(set_id=False)

        self.data["general"]['objpath'] = filepath
        
        with open(filepath) as fobj:
            group = "main_pars"    
                                   
            for line in fobj:
                line = line.strip()
                if line == '':
                    continue
                if line[0]=='[' and line[-1]==']':
                    group = line[1:-1]
                elif '=' in line:
                    lsplit = line.split()
                    if len(lsplit) != 3:
                        continue
                    key, eq, val = lsplit
                    if group == 'main_pars' and key in ["per","t0","v0","k1","k2","ecc","aop"]:
                        self.data[group][key] = float(val)
                    elif group == "general" and key in ["id","data"]:
                        self.data[group][key] = val
                    elif group == "rv_method":
                        if key == 'method':
                            mi = self.method_combo.findText(val)
                            if mi>=0: self.method_combo.setCurrentIndex(mi)
                        elif key in ['templ1','templ2']:
                            ti = int(key[-1:])-1
                            mi = self.tpl_combos[ti].findText(val)
                            if mi>=0: self.tpl_combos[ti].setCurrentIndex(mi)
                        elif key == 'mask':
                            mi = self.mask_combo.findText(val)
                            if mi>=0: self.mask_combo.setCurrentIndex(mi)
                        elif key == 'res':
                            self.r_spin.setValue(float(val))
                        self.data[group][key] = val
                    elif group in self.fitdial.rvmod_order:
                        if key in self.fitdial.rvmod_parorder[group] or \
                                    key[0] in ['a','b'] and key[1].isdigit():
                            self.data[group][key] = float(val)
                        elif key in ["nharm"]:
                            self.data[group][key] = int(val)
                        elif key in ["enabled"]:
                            self.data[group][key] = eval(val)
                        elif key in ["type"]:
                            self.data[group][key] = val
                else:
                    lsplit = line.split()
                    if len(lsplit) == 2:
                        if lsplit[0] in ["per","t0","v0","k1","k2","ecc","aop"]:
                            self.data[group][lsplit[0]] = float(lsplit[1])
                        elif lsplit[0] in ["id"]:
                            self.data[group][lsplit[0]] = lsplit[1]
                    
                                    
        if "id" not in self.data["general"]:
            self.data["general"]["id"] = obj_fid
        if "data" not in self.data["general"]:
            self.data["general"]['data'] = obj_fid

        self.setWindowTitle('RaveSpan - '+obj_fid)
        self.set_fields()
        self.id_changed(self.data["general"]["id"])

        self.curve.get_files(self.data["general"]['data'], obj_fid)

        self.spectrum.clear_spectrum()
        self.set_sdmode(False)
    
    
    
    def saveObject(self):
        if 'objpath' not in self.data["general"] or \
           self.data["general"]['objpath'] is None:
            self.saveAsObject()
            return
        oid = self.obj_id_field.text()
        odatapath = self.obj_path_field.text()
        odata = self.get_pars_from_fields(hjdmode="short")
        
        with open(self.data["general"]['objpath'],'w') as fobj:
            fobj.write("\n[general]\n")
            fobj.write("id     = %s\n"%str(oid))
            fobj.write("data   = %s\n"%str(odatapath))
            fobj.write("\n[main_pars]\n")
            for key,val in zip(["per","t0","ecc","aop","v0","k1","k2"], odata):
                fobj.write("%-6s = %s\n"%(key,val))
            fobj.write("\n[rv_method]\n")
            fobj.write("method = %s\n"%str(self.method_combo.currentText()))
            fobj.write("templ1 = %s\n"%str(self.tpl_combos[0].currentText()))
            fobj.write("templ2 = %s\n"%str(self.tpl_combos[1].currentText()))
            fobj.write("mask   = %s\n"%str(self.mask_combo.currentText()))            
            fobj.write("res    = %s\n"%str(self.r_spin.value()))
            self.fitdial.write_rv_mod_data(fobj)

        self.curve.save_rv_curve()

    def saveAsObject(self):
        filepath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save as...', self.obj_dir)[0]
        if not filepath:
            return  
        if not filepath.endswith(".obj"):
            filepath += ".obj"
        self.data["general"]['objpath'] = filepath
        self.saveObject()

    def browse_data(self):
        datapath = QtWidgets.QFileDialog.getExistingDirectory(self, 'Open file', self.specdir, QtWidgets.QFileDialog.ShowDirsOnly)
        if not datapath:
            return  

        spfiles = os.listdir(datapath)
        ngood = 0
        for spfile in spfiles:
            if spfile[0] in '0123456789' and spfile.count('_')>3:
                ngood+=1

        if ngood == 0:
            msgBox = QtWidgets.QMessageBox.warning(self, "Data Error", "No data in the selected directory is found.", QtWidgets.QMessageBox.Ok)
        
        datapath = datapath.replace(os.path.realpath(self.specdir),'')
        if datapath[0] == '/':
            datapath = datapath[1:]
        self.data["general"]['data'] = datapath
        self.obj_path_field.setText(datapath)
        self.curve.get_files(datapath, self.data["general"]['id'])


    def save_spectrum(self, xspec, spec):
        if spec is None: return
        specpath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save spectrum as...', self.datadir)[0]
        if specpath == "": return
        
        with open(specpath,'w') as fspec:
            for x,s in np.c_[xspec,spec]: # c_ is not defined
                fspec.write("%10.4f %8.4f\n"%(x,s))
            print("File written to:", specpath)

    def saveSpectrum(self):
        mode = "normalized" if self.cb_normalization.isChecked() else "raw"
        xspec,spec = self.spectrum.get_spectrum(mode)
        self.save_spectrum(xspec, spec)
    
    def saveTemplate1(self):
        self.saveTemplate(0)

    def saveTemplate2(self):
        self.saveTemplate(1)

    def saveTemplate(self, num, l2l1=1.0, binning=1):
        xspec, spec = self.spectrum.get_template(num, binw=binning)
        if spec is None: return
        specpath = QtWidgets.QFileDialog.getSaveFileName(self, 'Save template as...', self.datadir)[0]
        if specpath == "": return
        
        with open(specpath,'w') as fspec:
            if l2l1 != 1.0:
                if l2l1[1] == 0.0:
                    fl2l1 = l2l1[0]
                else:
                    fl2l1 = l2l1[0] + l2l1[1]*(xspec-l2l1[2])/1000.
                    lev = 0.04
                    if sum(fl2l1<lev) > 0:
                        fl2l1[fl2l1<lev] = lev
                        print("\n*** WARNING! Values of normalization function lower than 0.04 were converted to 0.04\n")
                    print("Some values of normalization function across spectrum printed below:")
                    print(fl2l1[::10000])
                if num == 0:
                    renorm_spec = (spec-1.)*(fl2l1+1.) + 1.
                else:
                    renorm_spec = (spec-1.)*(fl2l1+1.)/fl2l1 + 1.
            else:
                renorm_spec = spec
            for x,s in np.c_[xspec,renorm_spec]: # c_ is not defined
                fspec.write("%10.4f %8.4f\n"%(x,s))
            print("File written to:", specpath)
        
    
    def rescaleErrors(self):
        self.curve.rescale_errors()

    def about(self):
        about_text = "RaveSpan %s by B. Pilecki\n"%version
        about_text+= "(RV & spectrum analyzer)\n\n"
        about_text+= "RaveSpan is a simple and powerful tool \nfor analyzing spectra and radial velocity \ncurves.\n\n"
        about_text+= "Idea, graphical design, code writing, \n"
        about_text+= "CCF, TODCOR & BF implementations,\n"
        about_text+= "and Spectral Disentangling mode by:\n"
        about_text+= "\t\tBogumil Pilecki\n\n"
        about_text+= "Send comments, suggestions \nand report bugs to:\n"
        about_text+= "\tpilecki@camk.edu.pl\n\n"      
        about_text+= "Great thanks to:\n"
        about_text+= "\tMarek Gorski for help\n"
        about_text+= "\twith implementation of:\n"
        about_text+= "\tCCF & TODCOR,\n"
        about_text+= "and to:\n"
        about_text+= "\tSlavek Rucinski\n"
        about_text+= "\t& Theo Pribulla for help\n"
        about_text+= "\twith Broadening Function.\n\n"
        about_text+= "version: %s (%s)\n"%(version,reldate)
        msgBox = QtWidgets.QMessageBox.information(self, "About RaveSpan", about_text, QtWidgets.QMessageBox.Ok)
        


    """
    #################################################################
    ###                        GUI                                ###
    #################################################################
    """

    def create_menu(self, menu_obj, menu_list):
        """ Creates menu for given menu_list:
            [TYPE, (TEXT, FUNCTION (, TIP))]"""
        for menu_item in menu_list:
            item_type = menu_item[0]
            if item_type == 'a':
                item_obj, item_fun = menu_item[1:3]
                taction = QtWidgets.QAction(item_obj, self)
                if len(menu_item)>3:
                    taction.setStatusTip(menu_item[3])
                taction.triggered.connect(item_fun)
                menu_obj.addAction(taction)
            elif item_type == 's':
                menu_obj.addSeparator()

    def gui_menu(self):
        icon_path = "icons/"

        menubar = self.menuBar()
        mbFile = menubar.addMenu('&File')

        filemenu_order = ["new", "open", "save", "saveas", "--sep", "savespec", "savetpl1", "savetpl2", "--sep", "quit"]
        filemenu_items = {"new": ['gtk-new.png', "New", None, 'Create new object', self.newObject],
            "open": ['gtk-open.png', "Open", "Ctrl+O", 'Open object', self.openObject],
            'save': ['gtk-save.png', "Save", "Ctrl+S", 'Save object', self.saveObject],
          'saveas': ['gtk-save-as.png', "Save as", None, 'Save object as...', self.saveAsObject],
        'savespec': ['gtk-save-as.png', "Save spectrum", None, 'Save spectrum as...', self.saveSpectrum],
        'savetpl1': ['gtk-save-as.png', "Save template 1", None, 'Save template 1 as...', self.saveTemplate1],
        'savetpl2': ['gtk-save-as.png', "Save template 2", None, 'Save template 2 as...', self.saveTemplate2],
        'quit': ['gtk-quit.png', "Quit", "Ctrl+Q", 'Quit application', QCoreApplication.instance().quit]
        }
                        
        for fmitem in filemenu_order:
            if fmitem == '--sep':
                mbFile.addSeparator()
            else:
                icon, name, shortcut, statustip, eventfun = filemenu_items[fmitem]
                obj = QtWidgets.QAction(QtGui.QIcon(icon_path+icon), name, self)
                if shortcut is not None:
                    obj.setShortcut(shortcut)
                obj.setStatusTip(statustip)
                obj.triggered.connect(eventfun)
                mbFile.addAction(obj)

        self.instruments = InstrumentDialog(self)

        mbView = menubar.addMenu('&View')
        view_menu_list = [['a', "Viev RV", self.do_show_rv],
                          ['a', "View Spectrum", self.do_show_spec],
                          ['a', "View Analysis", self.do_show_anal],
                          ['a', "View Data", self.do_show_data]]
        self.create_menu(mbView, view_menu_list)

        self.sd_dialog = self.spectrum.sd_dialog
        self.editmask = EditMaskDialog(self)
        self.preferences = PreferencesDialog(self)
        self.preferences.move(self.geometry().left()+400,self.geometry().top())
 
        ### TOOLS
        mbTools = menubar.addMenu('&Tools')
        tools_menu_list = [
                      ['a', "Refresh templates", self.refreshTemplates, "Refresh Template lists"],
                      ['a', "Rescale errors", self.rescaleErrors, "Rescale errors according to actual scatter"],
                      ['a', "Edit mask file", self.editmask.show_or_activate, "Create/modify predefined masks"],
                      ['a', "Instrument setup", self.instruments.show_or_activate, "Show data gathered with selected instruments"],
                      ['s'],
                      ['a', "Spectral disentangling", self.sd_dialog.show_or_activate, "Dialog window for Spectral Disentangling"],
                      ['s'],
                      ['a', "Preferences", self.preferences.show_or_activate, "Customize RaveSpan behaviour"]]
        
        self.create_menu(mbTools, tools_menu_list)

        self.helpdial = HelpDialog(self)

        ### OPTIONS
        mbOptions = menubar.addMenu('&Info')
        options_menu_list = [
                      ['a', "About", self.about],
                      ['a', "Help", self.helpdial.show_or_activate]]
        self.create_menu(mbOptions, options_menu_list)


    def gui_object_data(self):

        self.obj_spin = {}
        self.obj_cb = {}

        object_group = QtWidgets.QGroupBox("Object:")
        self.main_layout.addWidget(object_group)
        self.object_layout = QtWidgets.QVBoxLayout(object_group)
        
        self.id_grid = QtWidgets.QGridLayout()
        self.object_layout.addLayout(self.id_grid)
        
        # ID
        idlabel = QtWidgets.QLabel('ID:', self.main_widget)
        self.obj_id_field = QtWidgets.QLineEdit("", self.main_widget)
        pb_browse_obj = QtWidgets.QPushButton("&Browse")
        pb_browse_obj.clicked.connect(self.openObject)
        self.obj_id_field.textChanged[str].connect(self.id_changed)
        self.id_grid.addWidget(idlabel,1,1)
        self.id_grid.addWidget(self.obj_id_field,1,2)
        self.id_grid.addWidget(pb_browse_obj,1,3)

        # PATH
        pathlabel = QtWidgets.QLabel('data:', self.main_widget)
        self.obj_path_field = QtWidgets.QPushButton("", self.main_widget)
        self.obj_path_field.clicked.connect(self.browse_data)
        self.id_grid.addWidget(pathlabel,2,1)
        self.id_grid.addWidget(self.obj_path_field,2,2,1,2)

        # EPHEMERIS - PERIOD
        self.ephem_grid = QtWidgets.QGridLayout()
        self.object_layout.addLayout(self.ephem_grid)
        self.obj_cb['per'] = QtWidgets.QCheckBox("P:", self.main_widget)
        self.obj_spin["per"] = set_QSpin(self,1.0,rmin=0.00001,rmax=99999.,step=0.1,decimals=8)
        self.obj_spin["per"].valueChanged.connect(self.param_changed)
        self.obj_spin["per"].editingFinished.connect(self.per_changed)
        self.ephem_grid.addWidget(self.obj_cb['per'],1,1)
        self.ephem_grid.addWidget(self.obj_spin["per"],1,2,1,2)

        # EPHEMERIS - HJD0
        self.obj_cb['t0'] = QtWidgets.QCheckBox('T0:', self.main_widget)
        self.obj_spin["t0"] = set_QSpin(self,5000.0,rmin=0.0,rmax=9999.99999,step=0.1,decimals=6)
        self.obj_spin["t0"].valueChanged.connect(self.param_changed)
        self.ephem_grid.addWidget(self.obj_cb['t0'],2,1)
        self.ephem_grid.addWidget(self.obj_spin["t0"],2,2,1,2)

        self.obj_cb['ecc'] = QtWidgets.QCheckBox('ecc:', self.main_widget)
        self.obj_spin["ecc"] = set_QSpin(self,0.0,rmin=0.0,rmax=0.98,step=0.01,decimals=6)
        self.obj_spin["ecc"].valueChanged.connect(self.param_changed)
        self.ephem_grid.addWidget(self.obj_cb['ecc'],3,1)
        self.ephem_grid.addWidget(self.obj_spin["ecc"],3,2)

        self.obj_cb['aop'] = QtWidgets.QCheckBox('aop:', self.main_widget)
        self.obj_spin["aop"] = set_QSpin(self,0.0,rmin=-0.1,rmax=dpi+0.1,step=0.1,decimals=5)
        self.obj_spin["aop"].valueChanged.connect(self.aop_changed)
        self.ephem_grid.addWidget(self.obj_cb['aop'],4,1)
        self.ephem_grid.addWidget(self.obj_spin["aop"],4,2)


        # ORBIT
        self.obj_cb['v0'] = QtWidgets.QCheckBox('V0:', self.main_widget)
        self.obj_spin["v0"] = set_QSpin(self,0.0,rmin=-999,rmax=999,step=0.5, decimals=3)
        self.obj_spin["v0"].valueChanged.connect(self.param_changed)
        self.ephem_grid.addWidget(self.obj_cb['v0'],1,4)
        self.ephem_grid.addWidget(self.obj_spin["v0"],1,5)


        self.obj_cb['k1'] = QtWidgets.QCheckBox('K1:', self.main_widget)
        self.obj_cb['k2'] = QtWidgets.QCheckBox('K2:', self.main_widget)
        self.obj_spin["k1"] = set_QSpin(self,100.0,rmin=0.0,rmax=999.,step=0.2, decimals=4)
        self.obj_spin["k2"] = set_QSpin(self,50.0,rmin=-999.0,rmax=999.,step=0.2, decimals=4)
        self.obj_spin["k1"].valueChanged.connect(self.param_changed)
        self.obj_spin["k2"].valueChanged.connect(self.param_changed)
        self.ephem_grid.addWidget(self.obj_cb['k1'],2,4)
        self.ephem_grid.addWidget(self.obj_spin["k1"],2,5)
        self.ephem_grid.addWidget(self.obj_cb['k2'],3,4)
        self.ephem_grid.addWidget(self.obj_spin["k2"],3,5)

        self.obj_cb['per'].setChecked(0)
        self.obj_cb['t0'].setChecked(0)
        self.obj_cb['aop'].setChecked(0)
        self.obj_cb['ecc'].setChecked(0)
        self.obj_cb['v0'].setChecked(1)
        self.obj_cb['k1'].setChecked(1)
        self.obj_cb['k2'].setChecked(1)

        pb_fit_more = QtWidgets.QPushButton('+')
        pb_fit_more.setStatusTip("Additional fitting parameters.")
        pb_fit_more.setMaximumWidth(self.smallbutton_width)

        self.fitdial = FittingDialog(self)
        pb_fit_more.clicked.connect(self.fitdial.show_or_activate)

        pb_fit_params = QtWidgets.QPushButton("Fit")
        pb_fit_params.setStatusTip("Fit selected parameters to the data.")
        pb_fit_params.clicked.connect(self.fit_params)
        self.ephem_grid.addWidget(pb_fit_more,4,4)
        self.ephem_grid.addWidget(pb_fit_params,4,5)



    def gui_rv_calculations(self):
        """ GUI RV Calculations Group Box """
        self.rv_calc_group = QtWidgets.QGroupBox("RV calculations:")
        self.main_layout.addWidget(self.rv_calc_group)
        self.rv_calc_grid = QtWidgets.QGridLayout(self.rv_calc_group)
        self.rv_calc_grid.setColumnStretch(6,1)
        self.grid_row = 0
        self.gui_method_selection()     # 1 - Method LAYOUT
        self.gui_template_selection()   # 2 - Template LAYOUT
        self.gui_mask_selection()       # 3 - Mask LAYOUT
        self.gui_calc_resolution()      # 4 - Resolution LAYOUT
        self.gui_calc_butts()           # 5 - Calculate Buttons LAYOUT
        self.gui_calc_opts()

    def gui_method_selection(self):    
        self.method_combo = QtWidgets.QComboBox(self)
        button_group = QtWidgets.QButtonGroup()
        self.t1_radio = QtWidgets.QRadioButton("T1")
        self.t1_radio.setStatusTip("Use template 1 for this method.")
        self.t1_radio.toggle()
        self.t2_radio = QtWidgets.QRadioButton("T2")
        self.t2_radio.setStatusTip("Use template 2 for this method.")
        button_group.addButton(self.t1_radio)
        button_group.addButton(self.t2_radio)
        self.method_combo.setStatusTip("Select analysis method to use.")
        for method in self.RVMethods:
            self.method_combo.addItem(method.upper())

        self.method_combo.currentIndexChanged[str].connect(self.method_activated)
        
        self.rv_calc_grid.addWidget(self.method_combo,self.grid_row,1,1,2)
        self.rv_calc_grid.addWidget(self.t1_radio,self.grid_row+1,1,1,1)
        self.rv_calc_grid.addWidget(self.t2_radio,self.grid_row+1,2,1,1)

    def gui_template_selection(self):

        self.tpl_buttons = [None,None]
        self.tpl_combos = [None,None]
        for i in range(2):
            self.tpl_buttons[i] = QtWidgets.QPushButton('T%d'%(i+1),self)
            self.tpl_buttons[i].setStatusTip('Select template %d by parameters.'%(i+1))
            self.tpl_buttons[i].setMaximumWidth(self.smallbutton_width)
            
            self.tpl_combos[i] = QtWidgets.QComboBox(self)
            self.add_templates(i)
            self.tpl_combos[i].setStatusTip('Select template %d by name.'%(i+1))

        self.tdialog = TemplateDialog(self)

        self.tpl_buttons[0].clicked.connect(lambda: self.template_pressed(0))
        self.tpl_buttons[1].clicked.connect(lambda: self.template_pressed(1))
        self.tpl_combos[0].currentIndexChanged[str].connect(lambda s: self.template_activated(s,0))
        self.tpl_combos[1].currentIndexChanged[str].connect(lambda s: self.template_activated(s,1))
        
        self.rv_calc_grid.addWidget(self.tpl_buttons[0],self.grid_row,3)
        self.rv_calc_grid.addWidget(self.tpl_buttons[1],self.grid_row+1,3)
        self.rv_calc_grid.addWidget(self.tpl_combos[0],self.grid_row,4,1,4)
        self.rv_calc_grid.addWidget(self.tpl_combos[1],self.grid_row+1,4,1,4)
        self.grid_row += 2

    def gui_mask_selection(self):
        d_label = QtWidgets.QLabel('Date:', self.main_widget)
        d_label.setAlignment(QtCore.Qt.AlignLeft|QtCore.Qt.AlignCenter)
        self.date = QtWidgets.QLabel('not set', self.main_widget)
        self.date.setEnabled(False)
        self.rv_calc_grid.addWidget(d_label,self.grid_row,1,1,1)
        self.rv_calc_grid.addWidget(self.date,self.grid_row,2,1,3)

        i_label = QtWidgets.QLabel('Instr:', self.main_widget)
        i_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignCenter)
        self.instr_name = QtWidgets.QLabel('not set', self.main_widget)
        self.instr_name.setEnabled(False)
        self.rv_calc_grid.addWidget(i_label,self.grid_row,5,1,1)
        self.rv_calc_grid.addWidget(self.instr_name,self.grid_row,6,1,2)

        m_label = QtWidgets.QLabel('Mask:', self.main_widget)
        m_label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignCenter)
        self.mask_combo = QtWidgets.QComboBox(self)
        self.add_masks(self.mask_combo)
        self.mask_combo.currentIndexChanged[str].connect(self.mask_activated)
        self.rv_calc_grid.addWidget(m_label,self.grid_row+1,4,1,2)
        self.rv_calc_grid.addWidget(self.mask_combo,self.grid_row+1,6,1,2)
        self.grid_row += 2

    def gui_calc_resolution(self):
        r_label = QtWidgets.QLabel('Resolution:', self.main_widget)
        self.r_spin = set_QSpin(self,1.0,rmin=0.1,rmax=19.9,step=0.1,suffix=" km/s")
        self.rv_calc_grid.addWidget(r_label,self.grid_row,4,1,2)
        self.rv_calc_grid.addWidget(self.r_spin,self.grid_row,6,1,2)
        self.grid_row += 1

    def gui_calc_butts(self):
        calc_one_butt = QtWidgets.QPushButton("Calc &One",self)
        next_butt = QtWidgets.QPushButton(">",self)
        prev_butt = QtWidgets.QPushButton("<",self)
        next_butt.setMaximumWidth(self.smallbutton_width)
        prev_butt.setMaximumWidth(self.smallbutton_width)
        self.calc_all_butt = QtWidgets.QPushButton("Calc All",self)
        calc_one_butt.clicked.connect(self.calc_one_pressed)
        prev_butt.clicked.connect(self.curve.select_prev_datapoint)
        next_butt.clicked.connect(self.curve.select_next_datapoint)
        self.calc_all_butt.clicked.connect(self.calc_all_pressed)
        
        calc_one_butt.setStatusTip('Calculate RVs for selected point only.')
        prev_butt.setStatusTip('Go to previous spectrum.')
        next_butt.setStatusTip('Go to next spectrum.')
        self.calc_all_butt.setStatusTip('Calculate RVs for all points.')
        self.calc_all_butt.setToolTip('Warning! Long execution - double click to run.')
        self.rv_calc_grid.addWidget(prev_butt,self.grid_row,4,1,1)
        self.rv_calc_grid.addWidget(calc_one_butt,self.grid_row,6,1,2)
        self.rv_calc_grid.addWidget(next_butt,self.grid_row,5,1,1)
        self.rv_calc_grid.addWidget(self.calc_all_butt,self.grid_row+1,6,1,2)

        disent_butt = QtWidgets.QPushButton("SD",self)
        disent_butt.setMaximumWidth(66)
        disent_butt.setStatusTip("Show spectral disentangling window.")
        disent_butt.clicked.connect(self.sd_dialog.show_or_activate)
        self.rv_calc_grid.addWidget(disent_butt,self.grid_row+1,1,1,2)

        self.grid_row += 2


    def gui_calc_opts(self):
        self.cb_instr_mask=QtWidgets.QCheckBox("Instr-Mask")
        self.cb_auto_max=QtWidgets.QCheckBox("Auto&Max")
        self.cb_fit=QtWidgets.QCheckBox("&Fit")
        self.cb_auto_max.toggle()
        self.cb_fit.toggle()
        self.cb_instr_mask.setStatusTip('Link mask to the instrument.')
        self.cb_auto_max.setStatusTip('Find maxima automatically.')
        self.cb_fit.setStatusTip('Get velocities from fitted functions.')
        self.rv_calc_grid.addWidget(self.cb_instr_mask,3,1,1,2)
        self.rv_calc_grid.addWidget(self.cb_auto_max,4,1,1,2)
        self.rv_calc_grid.addWidget(self.cb_fit,5,1,1,2)

    def gui_show_options(self):
        """ GUI SHOW OPTIONS """

        show_options_group = QtWidgets.QGroupBox("Options:")
        self.main_layout.addWidget(show_options_group)

        self.show_grid = QtWidgets.QGridLayout(show_options_group)

        self.pb_show_rv = QtWidgets.QPushButton('Show RV')
        self.pb_show_rv.setStatusTip("Show RV Curve window.")
        self.pb_show_rv.clicked.connect(self.do_show_rv)

        self.show_grid.addWidget(self.pb_show_rv,1,1,1,2)

        self.rv_raw_radio = None
        self.rv_phased_radio = None
        self.cb_show_model = None
        self.cb_model_id = None

        self.pb_show_data= QtWidgets.QPushButton('Show Data')
        self.pb_show_data.setStatusTip("Show Data View window.")
        self.pb_show_data.clicked.connect(self.do_show_data)

        self.show_grid.addWidget(self.pb_show_data,2,1,1,2)

        self.pb_show_spec = QtWidgets.QPushButton('Spectrum')
        self.pb_show_spec.setStatusTip("Show Spectrum window.")
        self.pb_show_spec.clicked.connect(self.do_show_spec)

        self.show_grid.addWidget(self.pb_show_spec,1,3,1,2)

        # SPECTRUM
        self.cb_show_spec = None
        
        self.cb_show_tpl = [None,None]

        self.cb_normalization=QtWidgets.QCheckBox("Normalized")
        self.cb_normalization.setChecked(self.spectrum.cb_normalization.isChecked())
        self.cb_normalization.setStatusTip("Toggle normalization of selected spectrum.")
        self.cb_normalization.stateChanged.connect(self.sync_normalization)

        self.show_grid.addWidget(self.cb_normalization,2,3,1,2)

        self.pb_show_anal = QtWidgets.QPushButton('Analysis')
        self.pb_show_anal.setStatusTip('Show Analysis window.')
        self.pb_show_anal.clicked.connect(self.do_show_anal)
        self.show_grid.addWidget(self.pb_show_anal,1,5,1,2)

        self.cb_ca_update=QtWidgets.QCheckBox('CA update')  
        self.cb_ca_update.setStatusTip('Update Analysis while calculating all points.')
        self.cb_ca_update.stateChanged.connect(self.ca_update_changed)

        self.show_grid.addWidget(self.cb_ca_update,2,5,1,2)


    def gui_bottom_buttons(self):
        """ GUI BOTTOM BUTTONS """

        save_butt = QtWidgets.QPushButton("&Save",self)
        save_butt.setStatusTip('Save object data and RV curve.')
        save_butt.clicked.connect(self.saveObject)
        self.show_grid.addWidget(save_butt,3,1,1,2)

        prefs_butt = QtWidgets.QPushButton("&Prefs",self)
        prefs_butt.setStatusTip('Edit preferences.')
        prefs_butt.clicked.connect(self.preferences.show_or_activate)
        self.show_grid.addWidget(prefs_butt,3,3,1,2)

        quit_butt = QtWidgets.QPushButton("&Quit",self)
        quit_butt.setStatusTip('Quit the application.')
        quit_butt.clicked.connect(QCoreApplication.instance().quit)
        self.show_grid.addWidget(quit_butt,3,5,1,2)



if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)  # start Qt application
    menu = RV_MENU()                    # create Matplotlib widget
    menu.show()

    sys.exit(app.exec_())     # exit with the same return code of Qt application
    
    

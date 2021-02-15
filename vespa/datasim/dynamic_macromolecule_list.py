# Python modules

# Third party modules
import numpy as np
import wx
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT

# Our modules
import vespa.common.minf_parabolic_info as minf
import vespa.common.wx_gravy.common_dialogs as common_dialogs

from vespa.datasim.util_datasim import calc_lw



##### Dynamic Metabolite List Class ###########################################

class DynamicMacromoleculeList(object):

    def __init__(self, grid_parent, tab, GridSizer, datasim):
        
        self.tab         = tab
        self.grid_parent = grid_parent
        self.datasim     = datasim
        self.list_lines  = []

        self.minppm = datasim.pts2ppm(datasim.dims[0])
        self.maxppm = datasim.pts2ppm(0)

        # We follow the wx CamelCaps naming convention for this wx object.
        self.GridSizer = GridSizer
        self.set_new_values()


    def __str__(self):
        """ width is in PPM and we calculate width_hz and width_damp values """

        lines = []
        vals = self.get_values()
        for i,line in enumerate(vals):
            check = str(line['checkbox'])
            ppm   = "{:.4f}".format(str(line['ppm']))
            area  = "{:.4f}".format(str(line['area']))
            phase = "{:.4f}".format(str(line['phase']))
            width = "{:.4f}".format(str(line['width']))
            width_hz   = "{:.4f}".format(str(self.width_ppm2hz(width)))
            width_damp = "{:.4f}".format(str(self.width_ppm2damp(width)))

            lines.append("Line "+str(i)+
                         ", Check = "+check+
                         ", PPM = "+ppm+
                         ", Area = "+area+
                         ", Phase [deg] = "+phase+
                         ", Width [ppm] = "+width+
                         ", Width [Hz] = "+width_hz+
                         ", Damp [sec] = "+width_damp)

        return '\n'.join(lines)


    def width_ppm2hz(self, val):
        """ val here is in [ppm], convert to [hz]"""
        return self.datasim.ppm2hz(val, rel=True)

    def width_ppm2damp(self, val):
        """ val here is in [ppm], convert to [sec]"""
        width_hz = self.datasim.ppm2hz(val, rel=True)
        damp = self._optimize_damp_given_hz(width_hz, tb=100000.0) # Pure Lorentz
        return damp

    def width_hz2ppm(self, val):
        """ val here is in [hz], convert to [ppm]"""
        return self.datasim.hz2ppm(val, rel=True)

    def width_damp2ppm(self, val):
        """ val here is in [sec], convert to [ppm]"""
        ta = val
        tb = 100000.0                       # Pure Lorentz
        hz = calc_lw(ta, tb)
        width = self.width_hz2ppm(hz)
        return width

    def _optimize_damp_given_hz(self, lw, ta=-1.0, tb=-1.0):

        info = {'ta': ta, 'tb': tb, 'orig_lw': lw }

        def _lw_function(val, info):
            ta = val if info["ta"] == -1 else info["ta"]
            tb = val if info["tb"] == -1 else info["tb"]
            width_hz = calc_lw(ta, tb)
            return np.abs(info["orig_lw"] - width_hz)

        # Call parabolic interpolation, Brent's method 1D minimization routine
        val_a = 0.00001  # lower bound
        val_b = 0.06     # "some" point in the middle
        val_c = 2000.0      # upper bound
        res, maxit = minf.minf_parabolic_info(val_a, val_b, val_c, _lw_function, info)

        return res


    def set_new_values(self):

        self.minppm = self.datasim.pts2ppm(self.datasim.dims[0]-1)
        self.maxppm = self.datasim.pts2ppm(0)

        # Set up default values from datasim
        self.checks = self.datasim.mmol_flags
        self.ppms   = self.datasim.mmol_ppms
        self.areas  = self.datasim.mmol_areas
        self.phases = self.datasim.mmol_phases
        self.widths = self.datasim.mmol_widths

        # error check PPM range, mainly for Resolution change and HLSVD import
        if not ((self.ppms > self.minppm).all() and (self.ppms < self.maxppm).all()):
            msg = 'Warning: Imported PPM values are outside current min/max PPM range, setting to boundary values.'
            common_dialogs.message(msg, "TabDatasim - Import Analysis Hlsvd Results File", common_dialogs.E_OK)
            self.ppms = np.clip(self.ppms, self.minppm, self.maxppm)
            self.datasim.mmol_ppms = self.ppms

        self.remove_all_rows()

        self.grid_parent.Freeze()
        for check, ppm, area, phase, width in zip(self.checks, self.ppms, self.areas, self.phases, self.widths):
            width_hz   = self.width_ppm2hz(width)
            width_damp = self.width_ppm2damp(width)
            self.add_row(check, ppm, area, phase, width, width_hz, width_damp)
        self.grid_parent.Thaw()
        self.grid_parent.Fit()
        self.grid_parent.Layout()
        
    
    def add_row(self, _check, _ppm, _area, _phase, _width, _width_hz, _width_damp):
        """Adds a row to the end of the list."""

        # create widgets to go into the line
        list_line = { }

        checkbox = wx.CheckBox(self.grid_parent, -1, '')
        
        ppm   = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        area  = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        phase = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        width = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        widhz = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        damp  = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        
        # keep a copy of panel and widgets to access later
        line = { "checkbox":checkbox, "ppm":ppm, "area":area, "phase":phase, "width":width, "widhz":widhz, "damp":damp}

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["checkbox"], 0, wx.ALIGN_CENTER_VERTICAL)
        for key in ("ppm","area","phase","width","widhz","damp"):
            self.GridSizer.Add(line[key], 0, wx.ALIGN_CENTER_VERTICAL)

        # Configure the controls I just created

        checkbox.SetValue(_check)

        # All of the floatspins have the same size. 
        fs_size = wx.Size(90, -1)

        # Note. On these Spin and FloatSpin widgets, if the value you want to
        #    set is outside the wxGlade standard range, you should make the 
        #    call to reset the range first and then set the value you want.
        ppm.SetDigits(3)
        ppm.SetIncrement(0.05)
        ppm.SetRange(self.minppm,self.maxppm)
        ppm.SetValue(_ppm)
        ppm.SetMinSize(fs_size)

        area.SetDigits(5)
        area.SetIncrement(0.05)
        area.SetRange(0.00001,10000.0)
        area.SetValue(_area)
        area.SetMinSize(fs_size)

        phase.SetDigits(2)
        phase.SetIncrement(1.0)
        phase.SetRange(-360.0,360.0)
        phase.SetValue(_phase)
        phase.SetMinSize(fs_size)

        width.SetDigits(7)
        width.SetIncrement(0.05)
        width.SetRange(0.00001,10000.0)
        width.SetValue(_width)
        width.SetMinSize(fs_size)

        widhz.SetDigits(7)
        widhz.SetIncrement(5.0)
        widhz.SetRange(0.00001,10000.0)
        widhz.SetValue(_width_hz)
        widhz.SetMinSize(fs_size)

        damp.SetDigits(7)
        damp.SetIncrement(0.01)
        damp.SetRange(-10000.0,10000.0)
        damp.SetValue(_width_damp)
        damp.SetMinSize(fs_size)

        self.list_lines.append(line)

        self.grid_parent.Bind(wx.EVT_CHECKBOX, self.update_macromolecules, checkbox)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_macromolecules, ppm)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_macromolecules, area)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_macromolecules, phase)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_macromol_width, width)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_macromol_widhz, widhz)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_macromol_damp,  damp)

        self.tab.PanelMacromoleculeLines.Layout()
        self.tab.Layout()

        
    def remove_checked_rows(self):

        checklist = []
        for i, line in enumerate(self.list_lines):
            if line["checkbox"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later indices in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self.list_lines[i].values()):
                if hasattr(item, "Destroy"):                    # It's a wx control
                    item.Destroy()
            del self.list_lines[i]
            
        self.GridSizer.Layout()
        self.tab.PanelBaselineSignals.Layout()


    def remove_all_rows(self):
        self.select_all()
        self.remove_checked_rows()

        
    def get_values(self):
        return [self.get_line_values(line) for line in self.list_lines]


    def get_line_values(self, line):
        return { "checkbox"   : line["checkbox"].GetValue(),
                 "ppm"        : line["ppm"].GetValue(),
                 "area"       : line["area"].GetValue(),
                 "phase"      : line["phase"].GetValue(),
                 "width"      : line["width"].GetValue(),
                 "width_hz"   : line["widhz"].GetValue(),
                 "width_damp" : line["damp"].GetValue()}

    def select_all(self):
        for line in self.list_lines:
            line["checkbox"].SetValue(True)


    def deselect_all(self):
        for line in self.list_lines:
            line["checkbox"].SetValue(False)

    def update_macromol_width(self, event):
        val_ppm  = event.GetEventObject().GetValue()
        val_hz   = self.width_ppm2hz(val_ppm)
        val_damp = self.width_ppm2damp(val_ppm)
        for line in self.list_lines:
            if line["width"] == event.GetEventObject():
                line["widhz"].SetValue(val_hz)
                line["damp"].SetValue(val_damp)
        self.update_macromolecules(event)

    def update_macromol_widhz(self, event):
        val_hz   = event.GetEventObject().GetValue()
        val_ppm  = self.width_hz2ppm(val_hz)
        val_damp = self.width_ppm2damp(val_ppm)
        for line in self.list_lines:
            if line["widhz"] == event.GetEventObject():
                line["width"].SetValue(val_ppm)
                line["damp"].SetValue(val_damp)
        self.update_macromolecules(event)

    def update_macromol_damp(self, event):
        val_damp = event.GetEventObject().GetValue()
        val_ppm  = self.width_damp2ppm(val_damp)
        val_hz   = self.width_ppm2hz(val_ppm)
        for line in self.list_lines:
            if line["damp"] == event.GetEventObject():
                line["width"].SetValue(val_ppm)
                line["widhz"].SetValue(val_hz)
        self.update_macromolecules(event)

    def update_macromolecules(self, event):
        self.tab.on_dynamic_macromolecule_list(None)
   
   
    def get_macromolecule_settings(self):
        
        checks, ppms, areas, phases, widths, widhzs, damps  = [], [], [], [], [], [], []
        vals = self.get_values()

        for i,line in enumerate(vals):
            checks.append(True) if line['checkbox'] else checks.append(False)
            ppms.append(line['ppm'])   
            areas.append(line['area'])
            phases.append(line['phase'])
            widths.append(line['width'])
            widhzs.append(line['width_hz'])
            damps.append(line['width_damp'])

        checks = np.array(checks)
        ppms   = np.array(ppms)
        areas  = np.array(areas)
        phases = np.array(phases)
        widths = np.array(widths)
        widhzs = np.array(widhzs)
        damps  = np.array(damps)

        return checks, ppms, areas, phases, widths, widhzs, damps

    

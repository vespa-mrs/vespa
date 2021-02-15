# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.common.util.ppm as util_ppm

from vespa.analysis.constants import GisoDefaultFixedT2
from vespa.analysis.constants import GisoDynMetabolite

from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY
from vespa.common.wx_gravy.widgets.floatspin_multiplier.floatspin_multiplier_base import FloatSpinMultiplier




class DynamicListGisoMetabolite(object):

    def __init__(self, _inner_notebook, tab, GridSizer, dataset, external_event_handler):
        
        self.tab             = tab
        self.dataset         = dataset
        self._inner_notebook = _inner_notebook
        
        self.external_event_handler = external_event_handler
        
        self._list_lines      = []

        # We follow the wx CamelCaps naming convention for this wx object.
        self.GridSizer = GridSizer


    def event_handler(self, event):
        self.external_event_handler(event)


    def get_values(self, only_checked=True):
        checks = []
        names = []
        area_scales = []
        peak_ppms = []
        search_ppms = []
        db_ppms = []
        fix_t2  = []
        search_ph0 = []
        
        lines = [self._get_line_values(line) for line in self._list_lines]

        for i,line in enumerate(lines):
            if not only_checked or line['checkbox']:
                checks.append(line['checkbox'])
                names.append(self.full_names[i])
                area_scales.append(line['scale'])   
                peak_ppms.append(line['shift'])
                search_ppms.append(line['width'])
                db_ppms.append(line['dbppm'])
                fix_t2.append(line['fixt2'])
                search_ph0.append(line['phas0'])
                
        if only_checked:
            return names, area_scales, peak_ppms, search_ppms, db_ppms, fix_t2, search_ph0
        else:
            return checks, names, area_scales, peak_ppms, search_ppms, db_ppms, fix_t2, search_ph0


    def select_all(self): 
        for line in self._list_lines: 
            line["checkbox"].SetValue(True) 
    
    
    def set_new_values(self, respect_current=True, preset=False):
        """
        Tells this control to grab new values from the block. If
        respect_current is True, the state of existing selections will 
        be applied to the new values if possible.
        """    
        if respect_current:
            tmp = self.get_values(False)
            prev_checks = tmp[0]
            prev_names = tmp[1] 
            prev_area_scales = tmp[2] 
            prev_peak_ppms = tmp[3]
            prev_search_ppms = tmp[4]
            prev_db_ppms = tmp[5]
            prev_fix_t2 = tmp[6]
            prev_search_ph0 = tmp[7]

        prior = self.tab.block.set.prior
        self.full_names = sorted(prior.basis_set_names)

        # Set up default values from prior or calculations
        # We sort the list of metab names.
        
        area_scales     = [1.0]    * len(self.full_names)
        checks          = [False]  * len(self.full_names)
        search_ppms     = [0.10]   * len(self.full_names) 
        peak_ppms       = [prior.basis_set[name].peak_ppm for name in self.full_names]
        db_ppms         = [False]  * len(self.full_names)
        fix_t2          = [1000.0] * len(self.full_names)
        search_ph0      = [0.0]     * len(self.full_names)
        all_ppms        = [prior.basis_set[name].all_ppms for name in self.full_names]
        
        if preset:
            for i, name in enumerate(list(self.tab.block.set.prior_list)):
                if name in self.full_names:
                    index = self.full_names.index(name)
                    checks[index]      = True
                    area_scales[index] = self.tab.block.set.prior_area_scale[i]
                    peak_ppms[index]   = self.tab.block.set.prior_peak_ppm[i]
                    search_ppms[index] = self.tab.block.set.prior_search_ppm[i]
                    db_ppms[index]     = self.tab.block.set.prior_db_ppm[i]
                    fix_t2[index]      = self.tab.block.set.prior_fix_t2[i]
                    search_ph0[index]  = self.tab.block.set.prior_search_ph0[i]


        # if loading a new prior file on an existing one then we need to
        # maintain existing user set parameters if metabolites in both sets
        if respect_current:
            if not len(prev_checks):
                # we are initializing for first time, take values from block
                prev_names       = list(self.tab.block.set.prior_list)
                prev_checks      = [True] * len(prev_names)
                prev_area_scales = list(self.tab.block.set.prior_area_scale)
                prev_peak_ppms   = list(self.tab.block.set.prior_peak_ppm)
                prev_search_ppms = list(self.tab.block.set.prior_search_ppm)
                prev_db_ppms     = self.tab.block.set.prior_db_ppm
                prev_fix_t2      = self.tab.block.set.prior_fix_t2
                prev_search_ph0  = list(self.tab.block.set.prior_search_ph0)
                
            for i, name in enumerate(prev_names):
                if name in self.full_names:
                    index = self.full_names.index(name)
                    checks[index]      = prev_checks[i]
                    area_scales[index] = prev_area_scales[i]
                    peak_ppms[index]   = prev_peak_ppms[i]
                    search_ppms[index] = prev_search_ppms[i]
                    db_ppms[index]     = prev_db_ppms[i]
                    fix_t2[index]      = prev_fix_t2[i]
                    search_ph0[index]  = prev_search_ph0[i]

        self._remove_all_rows()

        for arguments in zip(self.full_names, checks, area_scales, peak_ppms, search_ppms, db_ppms, fix_t2, search_ph0):
            self._add_row(*arguments)


    ##################################################################
    #
    #      Internal use/private functions    
    #
    
    def _add_row(self, full_name, flag, mult, pkppm, wid, dbflag, t2val, ph0val):
        '''Adds a row to the end of the list.'''
        # helper calcs
        maxppm = self.dataset.pts2ppm(0)
        minppm = self.dataset.pts2ppm(self.dataset.spectral_dims[0]-1)

        self.GridSizer.SetRows(self.GridSizer.GetRows() + 1)

        # create widgets to go into the line
        list_line = { }

        checkbox = wx.CheckBox(self._inner_notebook, -1, full_name)
        
        scale = FloatSpinMultiplier(self._inner_notebook, 
                                    increment=1.25, 
                                    digits=5, 
                                    style=wx.SP_ARROW_KEYS|wx.SP_WRAP|wx.TE_PROCESS_ENTER, 
                                    agwStyle=FS_LEFT
                                   )
        shift = FloatSpin(self._inner_notebook, agwStyle=FS_LEFT)
        width = FloatSpin(self._inner_notebook, agwStyle=FS_LEFT)
        dbppm = wx.CheckBox(self._inner_notebook, -1, '')
        fixt2 = FloatSpin(self._inner_notebook, agwStyle=FS_LEFT)
        phas0 = FloatSpin(self._inner_notebook, agwStyle=FS_LEFT)
        
        # keep a copy of panel and widgets to access later
        line = { "checkbox" : checkbox, 
                 "scale"    : scale, 
                 "shift"    : shift,
                 "width"    : width,
                 "dbppm"    : dbppm,
                 "fixt2"    : fixt2,
                 "phas0"    : phas0  }

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["checkbox"], 0, wx.ALIGN_CENTER_VERTICAL)
        for key in ("scale", "shift", "width"):
            self.GridSizer.Add(line[key], 0, wx.ALIGN_CENTER_VERTICAL)
        self.GridSizer.Add(line["dbppm"], 0, wx.ALIGN_CENTER_VERTICAL|wx.ALIGN_CENTER_HORIZONTAL)
        for key in ("fixt2", "phas0"):
            self.GridSizer.Add(line[key], 0, wx.ALIGN_CENTER_VERTICAL)
        
        # Configure the controls I just created
        checkbox.SetValue(flag)

        # All of the floatspins have the same size. 
        floatspin_size = wx.Size(70, -1)

        # Note. On these Spin and FloatSpin widgets, if the value you want to
        #    set is outside the wxGlade standard range, you should make the 
        #    call to reset the range first and then set the value you want.
        scale.multiplier = GisoDynMetabolite.AREA_SCALE_MULT
        scale.SetDigits(GisoDynMetabolite.AREA_SCALE_DIGITS)
        scale.SetIncrement(GisoDynMetabolite.AREA_SCALE_INCR)
        scale.SetRange(GisoDynMetabolite.AREA_SCALE_MIN,GisoDynMetabolite.AREA_SCALE_MAX)
        scale.SetValue(mult)
        scale.SetMinSize(floatspin_size)
        
        shift.SetDigits(GisoDynMetabolite.SEARCH_CENTER_DIGITS)
        shift.SetIncrement(GisoDynMetabolite.SEARCH_CENTER_INCR)
        shift.SetRange(minppm,maxppm)
        shift.SetValue(pkppm)
        shift.SetMinSize(floatspin_size) 

        width.SetDigits(GisoDynMetabolite.SEARCH_WIDTH_DIGITS)
        width.SetIncrement(GisoDynMetabolite.SEARCH_WIDTH_INCR)
        width.SetRange(GisoDynMetabolite.SEARCH_WIDTH_MIN,GisoDynMetabolite.SEARCH_WIDTH_MAX)
        width.SetValue(wid)
        width.SetMinSize(floatspin_size) 

        dbppm.SetValue(dbflag)

        fixt2.SetDigits(GisoDefaultFixedT2.DIGITS)
        fixt2.SetIncrement(GisoDefaultFixedT2.INCR)
        fixt2.SetRange(GisoDefaultFixedT2.MIN,GisoDefaultFixedT2.MAX)
        fixt2.SetValue(t2val)
        fixt2.SetMinSize(floatspin_size) 

        phas0.SetDigits(GisoDynMetabolite.SEARCH_PHASE0_DIGITS)
        phas0.SetIncrement(GisoDynMetabolite.SEARCH_PHASE0_INCR)
        phas0.SetRange(GisoDynMetabolite.SEARCH_PHASE0_MIN,GisoDynMetabolite.SEARCH_PHASE0_MAX)
        phas0.SetValue(ph0val)
        phas0.SetMinSize(floatspin_size) 


        self._list_lines.append(line)

        self._inner_notebook.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)
        self._inner_notebook.Bind(EVT_FLOATSPIN,   self.event_handler, scale)
        self._inner_notebook.Bind(EVT_FLOATSPIN,   self.event_handler, shift)
        self._inner_notebook.Bind(EVT_FLOATSPIN,   self.event_handler, width)
        self._inner_notebook.Bind(wx.EVT_CHECKBOX, self.event_handler, dbppm)
        self._inner_notebook.Bind(EVT_FLOATSPIN,   self.event_handler, fixt2)
        self._inner_notebook.Bind(EVT_FLOATSPIN,   self.event_handler, phas0)
        
        self.tab.PanelMetabolite.Layout()
        
        
    def _remove_all_rows(self):
        for line in self._list_lines:
            for control in list(line.values()):
                control.Destroy()
        
        self._list_lines = [ ]
            
        # There's always one row that contains the headings.
        self.GridSizer.SetRows(1)
        self.GridSizer.Layout()

        
    def _get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "checkbox"     : line["checkbox"].GetValue(),
                 "scale"        : line["scale"].GetValue(),
                 "shift"        : line["shift"].GetValue(),
                 "width"        : line["width"].GetValue(),
                 "dbppm"        : line["dbppm"].GetValue(),
                 "fixt2"        : line["fixt2"].GetValue(),
                 "phas0"        : line["phas0"].GetValue()
               }


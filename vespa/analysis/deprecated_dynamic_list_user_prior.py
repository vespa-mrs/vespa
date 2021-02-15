# Python modules

# 3rd party modules
import wx

# Our modules
import vespa.common.util.ppm as util_ppm
import vespa.common.wx_gravy.util as wx_util

from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY


#------------------------------------------------------------------------------

class DynamicListUserPrior(object):

    def __init__(self, PanelLines, PanelPrior, GridSizer, dataset, external_event_handler):
        
        self._dataset = dataset
        self._PanelLines = PanelLines
        self._PanelPrior = PanelPrior
        
        self.external_event_handler = external_event_handler
        
        # We follow the wx CamelCaps naming convention for this wx object.
        self._GridSizer = GridSizer
                
        self._list_lines = []


    @property
    def lines(self):
        return [self._get_line_values(line) for line in self._list_lines]

        
    def set_new_values(self, previous=None):
        if not previous:
            previous = self._dataset.user_prior.spectrum.get_rows()
            
        self.select_all()
        self.remove_checked_rows()

        for row_vals in previous:
            self.add_row(row_vals)

    
    def add_row(self, row_vals, update=False):
        '''
        Adds a row to the end of the list. 
        
        '''
        maxppm = self._dataset.pts2ppm(0)
        minppm = self._dataset.pts2ppm(self._dataset.spectral_dims[0]-1)

        self._GridSizer.SetRows(self._GridSizer.GetRows() + 1)

        # create widgets to go into the line
        list_line = { }

        checkbox = wx.CheckBox(self._PanelLines)
        
        value_ppm    = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        value_area   = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        value_phase  = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        value_lwhz   = FloatSpin(self._PanelLines, agwStyle=FS_LEFT)
        
        # keep a copy of panel and widgets to access later
        line = { "check"        : checkbox, 
                 "value_ppm"    : value_ppm, 
                 "value_area"   : value_area, 
                 "value_phase"  : value_phase,
                 "value_lwhz"   : value_lwhz, 
               }

        # Add the controls to the grid sizer
        self._GridSizer.Add(line["check"], 0, wx.ALIGN_CENTER_VERTICAL)
        for key in ("value_ppm", "value_area", "value_phase", "value_lwhz"):
            self._GridSizer.Add(line[key], 0, wx.EXPAND)

        # Configure the controls I just created

        # All of the floatspins have the same size. 
        floatspin_size = wx.Size(70, -1)

        # Note. On these Spin and FloatSpin widgets, if the value you want to
        #    set is outside the wxGlade standard range, you should make the 
        #    call to reset the range first and then set the value you want.

        wx_util.configure_spin(value_ppm,  70, 2, 0.05,(minppm,maxppm))
        wx_util.configure_spin(value_area, 70, 3, 0.1, (0.001,100000.0))
        wx_util.configure_spin(value_phase,70, 1, 5.0, (-360,360))
        wx_util.configure_spin(value_lwhz, 70, 2, 1.0, (0.001,10000.0))

        checkbox.SetValue(row_vals[0])
        value_ppm.SetValue(row_vals[1])
        value_area.SetValue(row_vals[2])
        value_phase.SetValue(row_vals[3])
        value_lwhz.SetValue(row_vals[4])

        self._list_lines.append(line)

        # only need to update if a metabolite is added/removed from basis
        self._PanelLines.Bind(wx.EVT_CHECKBOX, self.event_handler, checkbox)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_ppm)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_area)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_phase)
        self._PanelLines.Bind(EVT_FLOATSPIN, self.event_handler, value_lwhz)

        self._PanelPrior.Layout()  

        if update:
            self.event_handler()
        
        
    def remove_checked_rows(self):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self._list_lines):
            if line["check"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self._list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    # It's a wx control
                    item.Destroy()
                
            del self._list_lines[i]
            
        # Reduce the # of rows in the grid sizer
        rows = self._GridSizer.GetRows()
        self._GridSizer.SetRows(rows - len(checklist))
        self._GridSizer.Layout()
        self._PanelPrior.Layout()        
        self.event_handler()


    def select_all(self):
        for line in self._list_lines:
            line["check"].SetValue(True)


    def deselect_all(self):
        for line in self._list_lines:
            line["check"].SetValue(False)
            
            
    def event_handler(self, event=None):
        self.external_event_handler(event)


################   Private Methods

    def _get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "check"        : line["check"].GetValue(),
                 "value_ppm"    : line["value_ppm"].GetValue(),
                 "value_area"   : line["value_area"].GetValue(),
                 "value_phase"  : line["value_phase"].GetValue(),
                 "value_lwhz"   : line["value_lwhz"].GetValue(),
                 "limit_ppm"    : 0.0,
                 "limit_area"   : 0.0,
                 "limit_phase"  : 0.0,
                 "limit_lwhz"   : 0.0
               }



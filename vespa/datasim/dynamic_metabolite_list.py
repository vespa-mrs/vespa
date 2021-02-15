# Python modules

# 3rd party modules
import wx

# Our modules
from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY
from vespa.common.wx_gravy.widgets.floatspin_multiplier.floatspin_multiplier_base import FloatSpinMultiplier




##### Dynamic Metabolite List Class ###########################################

class DynamicMetaboliteList(object):

    def __init__(self, grid_parent, tab, GridSizer, datasim):
        
        self.tab         = tab
        self.grid_parent = grid_parent
        self.datasim     = datasim
        self.list_lines  = []

        # We follow the wx CamelCaps naming convention for this wx object.
        self.GridSizer = GridSizer

        self.set_new_values()


    def __str__(self):
        return self.__unicode__()
    
    
    def __unicode__(self):
        lines = [ ]
        vals = self.get_values()
        for i,line in enumerate(vals):
            check = str(line['checkbox'])
            sname = self.full_names[i]
            scale = "{:.4f}".format(str(line['scale']))
            decay = "{:.4f}".format(str(line['decay']))

            lines.append("Line "+str(i)+", Check = "+check+", Metab = "+sname+", Scale = "+ppm+", T2 (Ta) Decay [sec] = "+width)

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)    
    

    def set_new_values(self):

        # Set up default values from prior or calculations
        self.full_names  = list(self.datasim.names)
        self.mets_scales = self.datasim.mets_scales
        self.mets_flags  = self.datasim.mets_flags
        self.mets_decays = self.datasim.mets_decays
            
        self.remove_all_rows()

        for full_name, flag, mult, decay in \
            zip(self.full_names, self.mets_flags, self.mets_scales, self.mets_decays):
            self.add_row(full_name, flag, mult, decay)
        
    
    def add_row(self, full_name, flag, mult, decayval):
        """
        Adds a row to the end of the list. 
        """

        # create widgets to go into the line
        list_line = { }

        checkbox = wx.CheckBox(self.grid_parent, -1, full_name)
        
        scale = FloatSpinMultiplier(self.grid_parent, increment=1.25, digits=5, 
                                     style=wx.SP_ARROW_KEYS|wx.SP_WRAP|wx.TE_PROCESS_ENTER, 
                                     agwStyle=FS_LEFT)
        decay = FloatSpin(self.grid_parent, agwStyle=FS_LEFT)
        
        # keep a copy of panel and widgets to access later
        line = { "checkbox" : checkbox, "scale" : scale, "decay" : decay}

        # Add the controls to the grid sizer
        self.GridSizer.Add(line["checkbox"], 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT , 10)
        for key in ("scale", "decay"):
            self.GridSizer.Add(line[key], 0, wx.ALIGN_CENTER_VERTICAL)

        # Configure the controls I just created

        checkbox.SetValue(flag)

        # All of the floatspins have the same size. 
        floatspin_size = wx.Size(90, -1)

        # Note. On these Spin and FloatSpin widgets, if the value you want to
        #    set is outside the wxGlade standard range, you should make the 
        #    call to reset the range first and then set the value you want.
        scale.multiplier = 1.25
        scale.SetDigits(5)
        scale.SetIncrement(1.0)
        scale.SetRange(0.00001,100000.0)
        scale.SetValue(mult)
        scale.SetMinSize(floatspin_size)
        
        decay.SetDigits(5)
        decay.SetIncrement(0.05)
        decay.SetRange(0.00001,10000.0)
        decay.SetValue(decayval)
        decay.SetMinSize(floatspin_size) 

        self.list_lines.append(line)

        self.grid_parent.Bind(wx.EVT_CHECKBOX, self.update_metabolites, checkbox)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_metabolites, scale)
        self.grid_parent.Bind(EVT_FLOATSPIN,   self.update_metabolites, decay)

        self.tab.Layout()
        
        
    def remove_checked_rows(self):
        # gather indices of all checked boxes
        checklist = []
        
        for i, line in enumerate(self.list_lines):
            if line["checkbox"].GetValue():
                checklist.append(i)
        
        # remove in reverse order so we don't invalidate later
        # indices by removing the ones preceding them in the list
        checklist.reverse()
        
        for i in checklist:
            # Each line is a dict of controls + mixture info
            for item in list(self.list_lines[i].values()):
                if hasattr(item, "Destroy"):
                    # It's a wx control
                    item.Destroy()
                
            del self.list_lines[i]
            
        self.GridSizer.Layout()


    def remove_all_rows(self):
        self.select_all()
        self.remove_checked_rows()

        
    def get_values(self):
        return [self.get_line_values(line) for line in self.list_lines]


    def get_line_values(self, line):
        # Returns a dict containing  values of the controls in the line.
        return { "checkbox"     : line["checkbox"].GetValue(),
                 "scale"        : line["scale"].GetValue(),
                 "decay"        : line["decay"].GetValue()
               }


    def select_all(self):
        for line in self.list_lines:
            line["checkbox"].SetValue(True)


    def deselect_all(self):
        for line in self.list_lines:
            line["checkbox"].SetValue(False)


    def update_metabolites(self, event):
        self.tab.on_dynamic_metabolite_list(None)
   
   
    def get_metabolite_settings(self):
        
        checks = []
        scales = []
        decays = []
        
        vals = self.get_values()

        for i,line in enumerate(vals):
            if line['checkbox']:
                checks.append(True)
            else:
                checks.append(False)
            scales.append(line['scale'])   
            decays.append(line['decay'])

        return checks, scales, decays


    def get_last_settings(self):
        
        checks = []
        names  = []
        scales = []
        decays = []
        
        vals = self.get_values()

        for i,line in enumerate(vals):
            checks.append(line['checkbox'])
            names.append(self.full_names[i])
            scales.append(line['scale'])  
            decays.append(line['decay']) 

        return checks, names, scales, decays
    
    
    
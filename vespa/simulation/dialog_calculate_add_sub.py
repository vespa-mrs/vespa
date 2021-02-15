# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.simulation.util_simulation_config as util_simulation_config
import vespa.simulation.auto_gui.calculate_add_sub as calculate_add_sub
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.wx_gravy.util as wx_util



#------------------------------------------------------------------------------
# Note. GUI Architecture/Style
#
# Many of the GUI components in Vespa are designed using the WxGlade 
# application to speed up development times. The GUI components are designed
# interactively and users can preview the resultant window/panel/dialog, but
# while event functions can be specified, only stub functions with those 
# names are created. The WxGlade files (with *.wxg extensions) are stored in
# the 'wxglade' subdirectory. The ouput of their code generation are stored
# in the 'auto_gui' subdirectory. 
#
# To used these GUI classes, each one is inherited into a unique 'vespa'
# class, where program specific initialization and other functionality are
# written.  Also, the original stub functions for widget event handlers are 
# overloaded to provide program specific event handling.
#------------------------------------------------------------------------------


class DialogCalculateAddSub(calculate_add_sub.MyDialog):
    """
    This dialog is used to create a new tab from the active tab. The user 
    selects which Loop of the Experiment contains the Off and On states of an
    edited pulse sequence simulation. This loop must be of at least dimension
    two, but can be more. The assumption is that the first entry in this Loop
    contains the Off state and is followed by the on state in the next entry.
    All other entries in this Loop are ignored. On selecting a Loop, this
    dialog checks only to be sure that it has at least a dimensionality of two.

    
    """
    def __init__(self, parent, experiment):
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        calculate_add_sub.MyDialog.__init__(self, parent)
        
        self.parent       = parent
        self.experiment   = experiment
        self.index        = None
        self.comment      = ''
        
        #------------------------------
        # Initialize widget controls
        
        self._initialize_controls()

        self.Layout()
        self.Fit()
        
        
        # Under Linux, the focus is elsewhere (?) unless I manually
        # set it to the list
        self.ChoiceLoop.SetFocus()
        


    ##### Event Handlers ######################################################

    def on_loop(self, event):

        self.index = self.ChoiceLoop.GetSelection()

        if len(self.experiment.dims[self.index]) < 2:

            for i,dim in enumerate(self.experiment.dims):
                if len(dim) >= 2:
                    self.index = i
                    self.ChoiceLoop.SetSelection(i)
                    break
    

    def on_ok(self,event):
        
        self.comment = self.TextComment.GetValue()
        # we were successful, close dialog
        self.EndModal(wx.ID_OK)


    def on_cancel(self, event):
        # Nothing to do here since changes are only written on "OK"
        self.EndModal(wx.ID_CANCEL)




    ##### Internal helper functions  ##########################################

        
    def _initialize_controls(self):

        #----------------------------------------------------------------------
        # set up Experiment Loop list controls

        # All we need from the pulse seq is loop labels.
        loop_labels = self.experiment.pulse_sequence.loop_labels

        self.ChoiceLoop.Clear()
        
        for i, dim in enumerate(self.experiment.dims):
            if dim == mrs_experiment.DEFAULT_LOOP:
                # This is an empty loop.
                self.ChoiceLoop.Append("Empty Loop")
            else:
                label = "Loop%d - %s" % (i + 1, loop_labels[i])
                self.ChoiceLoop.Append(label)

        for i,dim in enumerate(self.experiment.dims):
            if len(dim) >= 2:
                self.index = i
                self.ChoiceLoop.SetSelection(i)
                break

        cmt = "This Add-Sub Experiment was derived from existing Experiment ...\n" \
              "Name = "+self.experiment.name+"\n" \
              "UUID = "+self.experiment.id+"\n"
              
        self.TextComment.SetValue(cmt)


        # We add the OK & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOk, self.ButtonCancel = \
                    wx_util.add_ok_cancel(self, self.LabelOkCancelPlaceholder,
                                          self.on_ok)





        


# Python modules

import xml.etree.cElementTree as ElementTree

# 3rd party modules
import wx
import wx.html

# Our modules
import vespa.common.auto_gui.experiment_browser as gui_experiment_browser
import vespa.common.util.xml_ as util_xml
import vespa.common.wx_gravy.util as wx_util


class DialogExperimentBrowser(gui_experiment_browser.MyDialog):
    """Displays a dialog box that contains a list of experiments. The
    list can be filtered by certain criteria and a preview of each 
    experiment is given. The user is expected to choose a dialog from the 
    list.
    
    When the dialog closes, dialog.selected_experiment_id is set to 
    the id of the selected experiment, or None if the user didn't choose an
    experiment.
    """
    def __init__(self, parent, db):
        self.selected_experiment_id = None
        
        if not parent:
            parent = wx.GetApp().GetTopWindow()

        gui_experiment_browser.MyDialog.__init__(self, parent)
        
        self.db = db
        
        html_sizer = self.LabelHtml.GetContainingSizer()
        parent = self.LabelHtml.GetParent()
        self.LabelHtml.Destroy()
        
        self.html_ctrl = wx.html.HtmlWindow(parent)
        
        html_sizer.Add(self.html_ctrl, 1, wx.EXPAND|wx.ALIGN_TOP)
        
        self.html_ctrl.SetBackgroundColour(self.GetBackgroundColour())
        
        # We add the Open & Cancel buttons dynamically so that they're in the 
        # right order under OS X, GTK, Windows, etc.
        self.ButtonOpen, self.ButtonCancel = \
            wx_util.add_ok_cancel(self, self.LabelOpenCancelPlaceholder,
                                  self.on_open, ok_text="Open")

        # The default option is to display all nuclei & all b0s. 
        self.ComboIsotope.Append("any")
        self.ComboIsotope.SetSelection(0)
        self.ComboB0.Append("any")
        self.ComboB0.SetSelection(0)
        
        self.populate_list()
        
        self.SetSize( (600, 400) )
        
        self.Layout()
        
        self.Center()
        
        # I build sorted, unique lists of the isotopes and populate the 
        # combobox with the list. I only do this once per
        #invocation of the dialog box.
        isotopes = set()
        
        for experiment in self.experiments:
            isotopes.add(experiment.isotope)
            
        self.ComboIsotope.AppendItems(sorted(isotopes))
        
        b0s = [str(b0) for b0 in self.db.fetch_b0_bins(True)]
        self.ComboB0.AppendItems(b0s)
        
        # Under OS X & Linux, the focus is elsewhere (?) unless I manually
        # set it to the list 
        self.ListExperiments.SetFocus()   

        
    def on_open(self, event):
        index = self.ListExperiments.GetSelection()
        
        if index != wx.NOT_FOUND:
            self.selected_experiment_id = self.experiments[index].id
        
            self.Close()        
        #else:
            # On the Mac, a click on an empty listbox row doesn't generate
            # a click event but it does set the current selection to 
            # wx.NOT_FOUND (ouch). Therefore, it's possible (on the Mac,
            # at least), to get this dialog in a state where the open 
            # button is enabled but nothing is selected.
        
        
    def on_filter_change(self, event):
        self.populate_list()

        
    def on_list_click(self, event):
        index = self.ListExperiments.GetSelection()
        
        if index != wx.NOT_FOUND:
            experiment = self.experiments[index]
        else:
            experiment = None

        self.item_selected(experiment)
        
        
    def on_list_double_click(self, event):
        index = self.ListExperiments.GetSelection()
        
        if index != wx.NOT_FOUND:
            self.selected_experiment_id = self.experiments[index].id
        
            self.Close()        
        #else:
            # On the Mac, a click on an empty listbox row doesn't generate
            # a click event but it does set the current selection to 
            # wx.NOT_FOUND (ouch). Therefore, it's possible (on the Mac,
            # at least), to get this dialog in a state where the open 
            # button is enabled but nothing is selected.
                
        
    def item_selected(self, experiment=None):
        """To be called when an item is selected"""
        
        # Build a small HTML document using ElementTree. Using ElementTree is
        # a little more cumbersome than assembling ordinary strings, but 
        # it eliminates bugs related to improperly escaped HTML.
        html = ElementTree.Element("html")
        body = ElementTree.SubElement(html, "body")
        if experiment:
            p = ElementTree.SubElement(body, "p")
            b = util_xml.TextSubElement(p, "b", "Comment: ")
            b.tail = experiment.comment

            metabolites = [metabolite.name for metabolite in experiment.metabolites]
            p = ElementTree.SubElement(body, "p")
            b = util_xml.TextSubElement(p, "b", "Metabolites: ")
            b.tail = ", ".join(metabolites)

        html = ElementTree.tostring(html)
        self.html_ctrl.SetPage(html)
        
        # I can set the HTML control's background color during init, but
        # the control forgets that info every time I call SetPage() so I
        # have to reset it here.
        self.html_ctrl.SetBackgroundColour(self.GetBackgroundColour())

        
    def populate_list(self):
        # Clear the current selection's preview
        self.item_selected(None)

        # Figure out if the user selected any filtering criteria
        isotope = self.ComboIsotope.GetValue()
        b0 = self.ComboB0.GetValue()
        
        # If "any" is selected, we don't pass any criteria to the query
        if isotope == "any":
            isotope = None
        b0 = None if (b0 == "any") else b0
            
        self.experiments = self.db.fetch_experiment_previews(b0, isotope)
        
        # Populate the listbox
        names = [experiment.name for experiment in self.experiments]
        self.ListExperiments.Set(names)
        
        self.ButtonOpen.Enable(bool(names))

        if names:
            # Select the first item
            self.ListExperiments.SetSelection(0)
            self.item_selected(self.experiments[0])
                        
        
        
        

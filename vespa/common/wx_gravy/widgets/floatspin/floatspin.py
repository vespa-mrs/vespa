
"""\
FloatSpin objects

@copyright: 2016-2021     Brian J. Soher
@license: MIT (see LICENSE.txt) - THIS PROGRAM COMES WITH NO WARRANTY
"""

import wx
from edit_windows import ManagedBase, EditStylesMixin
import time
import common, compat, misc
import new_properties as np
from collections import OrderedDict

from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY



class EditFloatSpin(ManagedBase, EditStylesMixin):
    """  Class to handle FloatSpin objects  """
    WX_CLASS = 'FloatSpin'
    _PROPERTIES = ["Widget", "range", "value", "increment", "digits", "extrastyle", "style"]
    PROPERTIES = ManagedBase.PROPERTIES + _PROPERTIES + ManagedBase.EXTRA_PROPERTIES

    def __init__(self, name, parent, index):
        ManagedBase.__init__(self, name, parent, index)
        EditStylesMixin.__init__(self)

        # initialise instance properties
        self.range      = np.FloatRangePropertyA( "0.0, 100.0" )
        self.value      = np.SpinDoublePropertyA(0.0, val_range=(0.0,100.0), immediate=True, default_value=0.0)
        self.increment  = np.SpinDoublePropertyA(1.0, val_range=(0.0,100.0), immediate=True, default_value=1.0)
        self.digits     = np.SpinPropertyA(3, val_range=(0,20), immediate=True, default_value=3)
        self.extrastyle = WidgetExtraStyleProperty()
        

    def create_widget(self):
        mi,ma  = self.properties["range"].get_tuple()
        value  = self.properties["value"].get()
        incr   = self.properties["increment"].get()
        digits = self.properties["digits"].get()
        
        self.widget = FloatSpin(self.parent_window.widget, wx.ID_ANY, min_val=mi, max_val=ma,
                                value=value, increment=incr, digits=digits)
    

    def finish_widget_creation(self, level, sel_marker_parent=None):
        ManagedBase.finish_widget_creation(self, level, sel_marker_parent)
        self.widget.Bind(wx.EVT_CHILD_FOCUS, self._on_set_focus)
        self.widget.Bind(wx.EVT_SET_FOCUS, self._on_set_focus)
        self.widget.Bind(wx.EVT_SPIN, self.on_set_focus)


    def _on_set_focus(self, event):
        # within a short time window, we ignore focus events as these seem due losing focus
        if not misc.focused_widget is self and time.time()-misc.focused_time > 0.05:
            # don't set focused_widget during event, as this may cause crashes; use delay instead
            misc.set_focused_widget(self, delayed=True)
        event.Skip()


    def _properties_changed(self, modified, actions):  # from SpinCtrlDouble

        if not modified or "range" in modified and self.widget:
            mi,ma = self.properties["range"].get_tuple()
            self.widget.SetRange(mi, ma)
            self.properties["value"].set_range(mi,ma)
            self.properties["increment"].set_range(mi,ma)

        if not modified or "increment" in modified and self.widget:
            self.widget.SetIncrement(self.increment)

        if not modified or "digits" in modified and self.widget:
            self.widget.SetDigits(self.digits)

        if not modified or "value" in modified or "range" in modified:
            # check that value is inside range
            value_p = self.properties["value"]
            if value_p.is_active():
                mi,ma = self.properties["range"].get_tuple()
                value = value_p.get()
                if value<mi:
                    value_p.set(mi)
                    value = mi
                elif value>ma:
                    value_p.set(ma)
                    value = ma
                if self.widget:
                    self.widget.SetValue(value)

        # FS style changes do not need any action here
        
        EditStylesMixin._properties_changed(self, modified, actions)
        ManagedBase._properties_changed(self, modified, actions)
        

# end of class EditFloatSpin


class WidgetExtraStyleProperty(np.WidgetStyleProperty):
    # for Extra-Styles widget/settings for FloatSpin; e.g. FS_LEFT, FS_RIGHT, FS_CENTER,FS_READONLY

    FLAG_DESCRIPTION = OrderedDict()
    FLAG_DESCRIPTION['Extra-Styles'] = ['FS_LEFT', 'FS_CENTRE', 'FS_RIGHT', 'FS_READONLY']
    FLAG_NAMES  = sum( list(FLAG_DESCRIPTION.values()), [] )

    EXCLUDES = {'FS_LEFT'   : set(['FS_CENTRE', 'FS_RIGHT']),
                'FS_CENTRE' : set(['FS_LEFT',   'FS_RIGHT']),
                'FS_RIGHT'  : set(['FS_LEFT',   'FS_CENTRE']) }

    STYLE_DEFS = {
        'FS_LEFT':     { 'desc': _('Same as wxTE_LEFT for wxTextCtrl: the numerical text is left aligned.'),
                         'exclude': 'FS_CENTRE|FS_RIGHT'  },
        'FS_CENTRE':   { 'desc': _('Same as wxTE_CENTER for wxTextCtrl: the numerical text is center aligned.'),
                         'exclude': 'FS_LEFT|FS_RIGHT'  },
        'FS_RIGHT':    { 'desc': _('Same as wxTE_RIGHT for wxTextCtrl: the numerical text is right aligned.'),
                         'exclude': 'FS_LEFT|FS_CENTRE'  },
        'FS_READONLY': { 'desc': _('The numericalS text is read only.') },
        }

    def _ensure_values(self):
        """ from the base classs - but needs some rethink for this way of using it """
        self._values = []  # the associated flag values
        widget_writer = self.owner.widget_writer

        for name in self._names:
            wx_name = widget_writer.cn_f(name)
            if not wx_name:  # cn_f() returns an empty string if the given styles are not supported
                self._values.append(None)
                continue
            try:
                self._values.append( self.owner.wxname2attr(wx_name) )  # FIXME bjs - this always fails
            except:
                self._values.append(None)                               # but this is not objectionable as a fallback
                continue

    def set_owner(self, owner, attname):
        """ style information is taken from class constants above """
        np._CheckListProperty.set_owner(self, owner, attname)
        self.style_defs = self.STYLE_DEFS
        self.styles     = self.FLAG_DESCRIPTION
        self._names     = sum( list(self.styles.values()), [] )
        self._values    = None
        self.set('FS_LEFT')
        self.default_value = set(set(['FS_LEFT']))
        self.set("")
        self.modified = False

    def create_editor(self, panel, sizer):
        """ similar to base class but simpler """
        self._choices = [] # the checkboxes
        self._ensure_values()

        tooltips = self._create_tooltip_text()

        static_box = wx.StaticBox(panel, -1, _("Extra-Style"), style=wx.FULL_REPAINT_ON_RESIZE)
        box_sizer = wx.StaticBoxSizer(static_box, wx.VERTICAL)
        for name, flag_value in zip(self._names, self._values):
            if name in self.STYLE_DEFS:
                style_def = self.STYLE_DEFS[name]
            else:
                style_def = {}
            if "obsolete" in style_def or "rename_to" in style_def:
                self._choices.append(None)
                continue
            checkbox = wx.CheckBox(panel, -1, name)

            if name in tooltips:
                compat.SetToolTip( checkbox, tooltips[name] )

            self._choices.append(checkbox)
            box_sizer.Add(checkbox)

        sizer.Add(box_sizer, 0, wx.ALL | wx.EXPAND, 5) 

        self.update_display(True)
        for checkbox in self._choices:
            if checkbox is not None:
                checkbox.Bind(wx.EVT_CHECKBOX, self.on_checkbox)     


def builder(parent, index):
    """ factory function for EditFloatSpin objects.  """
    name = parent.toplevel_parent.get_next_contained_name('float_spin_%d')
    with parent.frozen():
        editor = EditFloatSpin(name, parent, index)
        editor.properties["style"].set_to_default()
        editor.check_defaults()
        if parent.widget: editor.create()
    return editor


def xml_builder(parser, base, name, parent, index):
    """ factory function to build EditFloatSpin objects from a XML file """
    editor = EditFloatSpin( name, parent, index)
    editor.properties["value"].set_active(False)
    return editor


def initialize():
    """ initialization function for the module: returns a wxBitmapButton to be added to the main palette. """
    import os
    common.widget_classes['EditFloatSpin'] = EditFloatSpin
    common.widgets['EditFloatSpin'] = builder
    common.widgets_from_xml['EditFloatSpin'] = xml_builder
    return common.make_object_button('EditFloatSpin', os.path.join(os.path.dirname(__file__), 'floatspin.png'))

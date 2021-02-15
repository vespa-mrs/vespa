#!/usr/bin/env python

# Python modules


# 3rd party modules
#import wx.aui as aui
import wx.lib.agw.aui as aui              #?? Not anymore in wxPython 4.0.6 ??


# NB. bjs 2020-03-23 - maybe deprecate seeing as how the only aui I use
#                      is set up directly without calling this class?


class AUIManager(aui.AuiManager):
    """
    AUI Manager class
    """

    #----------------------------------------------------------------------
    def __init__(self, managed_window):
        """Constructor"""
        aui.AuiManager.__init__(self)
        self.SetManagedWindow(managed_window)
# Python modules


# 3rd party modules
import wx
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??

# Our modules
import vespa.common.wx_gravy.util as wx_util



"""This module contains two classes that are supersets (subclasses) of
wx Notebooks.

VespaNotebook subclasses wx.Notebook and VespaAuiNotebook subclasses
(surprise!) wx.AuiNotebook.

Both Vespa classes add the same thing, which are the properties and
methods in _NotebookHelper. Don't instantiate _NotebookHelper directly;
that's not what it's there for.
"""

class _NotebookHelper(object):
    def __init__(self):
        pass

    @property
    def tabs(self):
        """Returns a (potentially empty) list of this notebook's tabs. If this
        notebook has a welcome tab, it's not included in the list."""
        return [self.GetPage(i) for i in range(self.GetPageCount())
                                if self.GetPageText(i) != self.WELCOME_TAB_TITLE]


    @property
    def active_tab(self):
        """Returns the currently active page/tab, or None if the current tab
        is the welcome tab.
        """
        page = None
        if self.is_welcome_tab_open:
            # Nothing to do; no other page can be open when the welcome
            # page is open.
            pass
        else:
            # GetSelection() returns the index of the current page.
            # Occasionally there isn't an active page.
            index = self.GetSelection()
            if index != -1:
                # Occasionally GetSelection() returns an index that's == GetPageCount(), which
                # raises an error when we try to use it since the index is 0-based. We use min()
                # to ensure we don't try to reference a page that doesn't exist.
                index = min(index, self.GetPageCount() - 1)
                page = self.GetPage(index)

        return page


    @property
    def is_welcome_tab_open(self):
        return self.GetPageText(0) == self.WELCOME_TAB_TITLE


    def activate_tab(self, tab=None, index=-1):
        """Given a tab or an index (but not both!), activates the associated
        tab."""
        if tab:
            index = self.get_tab_index(tab)

        # The wx docs I have (2.8.9.2) state that SetSelection() is
        # deprecated on wx.Notebooks in favor of ChangeSelection().
        # Unfortunately, wx.AuiNotebook doesn't provide ChangeSelection()
        # at all, only SetSelection(). The two are not equivalent.
        # SetSelection() generates a "page changed" event while
        # ChangeSelection() does not.
        self.SetSelection(index)


    def get_tab_index(self, find_this_tab):
        """Given a tab, returns the index of that tab in the notebook,
        or -1 if it's not found."""
        found = False
        for i, tab in enumerate(self.tabs):
            if tab is find_this_tab:
                found = True
                break

        return i if found else -1


    def get_tab_text(self, tab=None, index=-1):
        """Given a tab or an index (but not both!), returns the text
        on that tab."""
        if tab:
            index = self.get_tab_index(tab)

        return self.GetPageText(index) if index != -1 else ""


    def remove_tab(self, tab=None, index=-1):
        """Given a tab or an index (but not both!), removes that tab."""
        if tab:
            index = self.get_tab_index(tab)

        if index != -1:
            self.DeletePage(index)


    def get_tab_by_label(self, label):
        for i in range(self.GetPageCount()):
            if label == self.GetPageText(i):
                return self.GetPage(i)

        return None

    def get_tab_index_by_label(self, label):
        for i in range(self.GetPageCount()):
            if label == self.GetPageText(i):
                return i

        return None


    def show_welcome_tab(self):
        """Creates & displays the welcome page"""
        html_control = wx.html.HtmlWindow(self, -1, wx.DefaultPosition,
                                          wx.Size(400, 300))
        html_control.SetPage(self.WELCOME_TAB_TEXT)

        self.AddPage(html_control, self.WELCOME_TAB_TITLE, False)


class VespaAuiNotebook(aui.AuiNotebook, _NotebookHelper):
    DEFAULT_STYLE = wx.BORDER_NONE

    DEFAULT_AGW_STYLE = aui.AUI_NB_MIDDLE_CLICK_CLOSE    |       \
                        aui.AUI_NB_CLOSE_ON_ACTIVE_TAB   |       \
                        aui.AUI_NB_SCROLL_BUTTONS        |       \
                        aui.AUI_NB_TAB_MOVE              |       \
                        aui.AUI_NB_TAB_SPLIT             |       \
                        aui.AUI_NB_TOP 

    # Subclasses may override this but probably don't need to.
    WELCOME_TAB_TITLE = "Welcome Info"

    # Subclasses with a welcome tab definitely want to override this.
    # Subclasses without a welcome tab can just ignore it.
    WELCOME_TAB_TEXT = """
    <html><body>
    <h1>Welcome to Vespa</h1>
    </body></html>
    """

    def __init__(self, parent, style=0, agw_style=0):
        aui.AuiNotebook.__init__(self, parent)
        _NotebookHelper.__init__(self)

        if style:
            # Use the caller's style
            pass
        else:
            style = self.DEFAULT_STYLE

        if agw_style:
            # Use the caller's style
            pass
        else:
            agw_style = self.DEFAULT_AGW_STYLE

        self.SetWindowStyleFlag(style)
        
        self.SetAGWWindowStyleFlag(agw_style)
        
        # The name of the stock tab art changed somewhere between wx 2.8 and 3.0.
        if hasattr(aui, 'AuiDefaultTabArt'):
            tab_art_provider = aui.AuiDefaultTabArt
        else:
            tab_art_provider = aui.AuiGenericTabArt
            
        self.SetArtProvider(tab_art_provider())


class VespaNotebook(wx.Notebook, _NotebookHelper):
    def __init__(self, parent, style=0):
        wx.Notebook.__init__(self, parent)
        _NotebookHelper.__init__(self)

        if style:
            # Use the caller's style
            pass
        else:
            style = self.DEFAULT_STYLE

        self.SetWindowStyleFlag(style)

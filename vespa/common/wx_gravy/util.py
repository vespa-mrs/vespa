# Python modules

import os
import tempfile
import webbrowser

# 3rd party modules
import wx
#import wx.aui as aui
import wx.lib.agw.aui as aui        # NB. wx.aui version throws odd wxWidgets exception on Close/Exit  ?? Not anymore in wxPython 4.0.6 ??

# Our modules
import vespa.common.util.misc as util_misc

# wx.GetStockLabel() returns the label associated with a wx stock button
# id (like wx.ID_CANCEL). There's no function to do the reverse (map a label
# to an id) so _STOCK_ITEM_IDS fills that void. It's keyed by properly
# capitalized strings and also by lower-case versions of those strings.
# Not all wx stock ids are represented, only ones I thought we might use.
_STOCK_ITEM_IDS = {
                    "About"             : wx.ID_ABOUT,
                    "Add"               : wx.ID_ADD,
                    "Cancel"            : wx.ID_CANCEL,
                    "Clear"             : wx.ID_CLEAR,
                    "Close"             : wx.ID_CLOSE,
                    "Copy"              : wx.ID_COPY,
                    "Delete"            : wx.ID_DELETE,
                    "Edit"              : wx.ID_EDIT,
                    "Replace"           : wx.ID_REPLACE,
                    "New"               : wx.ID_NEW,
                    "No"                : wx.ID_NO,
                    "OK"                : wx.ID_OK,
                    "Open"              : wx.ID_OPEN,
                    "Preferences"       : wx.ID_PREFERENCES,
                    "Print"             : wx.ID_PRINT,
                    "Properties"        : wx.ID_PROPERTIES,
                    "Quit"              : wx.ID_EXIT,
                    "Refresh"           : wx.ID_REFRESH,
                    "Remove"            : wx.ID_REMOVE,
                    "Revert to Saved"   : wx.ID_REVERT_TO_SAVED,
                    "Save"              : wx.ID_SAVE,
                    "Save As..."        : wx.ID_SAVEAS,
                    "Select all"        : wx.ID_SELECTALL,
                    "Stop"              : wx.ID_STOP,
                    "Undelete"          : wx.ID_UNDELETE,
                    "Undo"              : wx.ID_UNDO,
                    "Yes"               : wx.ID_YES,
                  }

# Add lower case keys.
_STOCK_ITEM_IDS.update(dict([(key.lower(), value) for key, value
                                                  in _STOCK_ITEM_IDS.items()]))


def add_ok_cancel(dialog, placeholder, ok_event_handler=None,
                  cancel_event_handler=None, ok_text="OK"):
    """Adds OK & Cancel buttons to the dialog in the lower right. The
    placeholder parameter should be a control located in a sizer
    approximately where you want the OK & Cancel buttons. This function
    destroys the placeholder control.

    The newly-created OK and Cancel buttons are returned.

    If ok_text is something other than the default, the "OK" button gets
    that text instead. If the button text (stripped of '&' characters) is in
    this module's _STOCK_ITEM_IDS dict, it will be given the appropriate wx
    id.

    This function is useful because OK and Cancel buttons appear in a
    different order on Windows than on OS X & GTK, and getting them in the
    right place takes more code than one might expect.
    """
    parent = placeholder.GetParent()

    # Using the proper ids (e.g. wx.ID_CANCEL) ensures that these buttons
    # get the appropriate artwork under GTK.
    # Also, a button with id == ID_CANCEL is automatically bound to
    # self.Close()  | wx.ALIGN_CENTER_VERTICAL
    id_ = _STOCK_ITEM_IDS.get(ok_text.lower().replace("&", ""), wx.ID_ANY)

    ok = wx.Button(parent, id_, ok_text)
    cancel = wx.Button(parent, wx.ID_CANCEL, "Cancel")

    # The placeholder tells me where to add the OK & Cancel. Once I've used
    # it, I get rid of it.
    containing_sizer = placeholder.GetContainingSizer()
    parent = placeholder.GetParent()
    placeholder.Destroy()

    # The OK and Cancel buttons appear in a different order on Windows
    # than on OS X & GTK, but the wx.StdDialogButtonSizer handles that
    # for me.
    button_sizer = wx.StdDialogButtonSizer()

    button_sizer.AddButton(ok)
    button_sizer.AddButton(cancel)

    button_sizer.SetAffirmativeButton(ok)
    ok.SetDefault()

    if ok_event_handler:
        dialog.Bind(wx.EVT_BUTTON, ok_event_handler, ok)
    if cancel_event_handler:
        dialog.Bind(wx.EVT_BUTTON, cancel_event_handler, cancel)

    button_sizer.Realize()

    containing_sizer.Add(button_sizer, 1,
                         wx.BOTTOM | wx.ALIGN_BOTTOM | wx.ALIGN_RIGHT,
                         10)
    containing_sizer.Fit(dialog)

    # Under wx, creation order = tab order. Since I created the OK button
    # first, that's before Cancel in the tab order. That's not a problem
    # as long as OK is to the left of cancel but if it isn't, I have to
    # correct the tab order.
    ok_left, _ = ok.GetPosition()
    cancel_left, _ = cancel.GetPosition()
    if ok_left > cancel_left:
        # Ooops, OK is to the right of cancel
        ok.MoveAfterInTabOrder(cancel)

    return (ok, cancel)


def configure_spin(control, width, digits=None, increment=None, min_max=None):
    """A convenience function for setting common attributes of spin
    controls including regular integer spin controls, floatspins and
    floatspin multipliers. There is no return value.

    Note that integer spin controls will complain if you attempt to set
    digits or increment.

    The width supplied should be the width that makes the control look best on
    Windows. The standard control fonts are a little bigger under OS X and GTK.
    As a result, the controls on those platforms have to be a little wider to
    display the same amount of text. To correct for that, we multiply the
    width by a platform-dependent fudge factor.
    """
    # We start by saving the current value so that we can restore it at the
    # end of this code in case the interim changes alter it.
    val = control.GetValue()

    # Widen controls under non-Windows platforms
    if "__WXMAC__" in wx.PlatformInfo:
        width = int(width * 1.2)
    elif "__WXGTK__" in wx.PlatformInfo:
        width = int(width * 1.1)

    control.SetMinSize((width, -1))
    control.SetSize((width, -1))

    if digits is not None:
        control.SetDigits(digits)

    if increment is not None:
        control.SetIncrement(increment)

    if min_max is not None:
        control.SetRange(min_max[0], min_max[1])

    control.SetValue(val)


def copy_to_clipboard(some_text):
    """Copies some_text to the clipboard."""
    wx.TheClipboard.Open()
    data = wx.TextDataObject()
    data.SetText(some_text)
    wx.TheClipboard.SetData(data)
    wx.TheClipboard.Close()


def display_text_as_file(the_text):
    """Given a string (Unicode or ASCII), writes the string to disk and
    displays that file in the operating system's default text editor.
    """
    # We build strings using \n as the newline inidicator, but occasionally we
    # mix those strings with text files read from disk (e.g. a VASF params
    # file). In that case, the string will have mixed newline indicators
    # which confuses the heck out of most text editors.
    # We fix that here by ensuring that the newlines are consistent.
    # Note that Notepad (the default text viewer for many under Windows)
    # can't handle unadorned "\n" as a new line.
    newline = "\r\n" if wx.Platform == "__WXMSW__" else "\n"
    the_text = util_misc.normalize_newlines(the_text, newline)

    the_text = the_text.encode("utf-8")

    fd, filename = tempfile.mkstemp(".txt")
    os.write(fd, the_text)
    os.close(fd)

    display_file(filename)

    # It would be polite to remove the temp file I just created, but
    # since the call to webbrowser.open() is asynchronous, there's a good
    # chance that file removal will happen before the OS has had a
    # chance to open it. We leave it behind as trash that will hopefully
    # be cleaned up by the OS on the next reboot. Apologies for the mess!


def display_file(filename):
    """Given a filename, this function asks the operating system to open
    the file in the default application associated with the filetype. This
    isn't guaranteed to work for all file types, but has worked reliably
    for us with PDFs, PNGs and text files.
    """
    # Display the text file with the webbrowser module. As of Python 2.6, the
    # webbrowser documentation warns, "[O]n some platforms, trying to open a
    # filename using this function, may work and start the operating system's
    # associated program. However, this is neither supported nor portable."
    # Despite the warning, this works for us under OS X, Windows, Ubuntu,
    # Kubuntu and Fedora.
    #
    # Relative paths don't seem to work well here, so we ensure the path is
    # fully qualified. Also note that we must pass a URI rather than just a
    # raw path, hence the addition of "file://".
    filename = os.path.realpath(filename)
    webbrowser.open("file://" + filename)


def get_selected_item_indices(list_ctrl):
    """Returns a (possibly empty) list of the items selected in a wx.ListCtrl
    or wx.Listbox. Note that the latter class already has the method
    GetSelections(), so you don't really need to call this function for
    a wx.Listbox. There's no similar method on the wx.ListCtrl class,
    hence this function.
    """
    is_listbox = isinstance(list_ctrl, wx.ListBox)

    if is_listbox:
        indices = list_ctrl.GetSelections()
    else:
        indices = [ ]
        # wxPython doesn't provide a simple function for getting the index
        # of all selected items in a ListCtrl
        index = list_ctrl.GetFirstSelected()
        while index != -1:
            indices.append(index)
            index = list_ctrl.GetNextSelected(index)

    # Under Windows & GTK, the list returned by listbox.GetSelections() is
    # always sorted smallest to largest. Under OS X, it's not sorted and
    # might be e.g. (2, 1, 0) or even (1, 0, 2). I don't know how list
    # controls behave, but sorting is cheap so we always do it here to
    # eliminate surprises.
    return sorted(indices)


def is_select_all(event):
    """Given the event generated by EVT_KEY_DOWN, returns True if the key
    combo indicates select all (Ctrl+a, or Cmd+a on the Mac). Returns False
    otherwise."""
    return ((event.GetKeyCode() == ord("A")) and (event.GetModifiers() == wx.MOD_CMD))


def is_wx_floatspin_ok():
    """Returns True if the wxPython version is sufficiently recent to provide
    the floatspin features we require. False otherwise.

    "Sufficiently recent" means version 2.8.11.0 or greater (including 2.9).
    """
    version = wx.__version__.strip()
    version = version.split('.')
    version = [int(element) for element in version if element.strip()]
    return (version[1] > 8) or ((version[1] == 8) and (version[2] >= 11))


def select_list_ctrl_items(list_ctrl, item_indices=None, select=True):
    """Given a wx.ListCtrl or a wx.ListBox and the index of one or more
    items therein, selects or deselects the item(s). item_indices can be a
    None, a single index or a list of them. When left at its default (None),
    all items in the list are affected.
    """
    is_listbox = isinstance(list_ctrl, wx.ListBox)

    if item_indices is None:
        # operate on all items
        if is_listbox:
            item_indices = list(range(list_ctrl.GetCount()))
        else:
            item_indices = list(range(list_ctrl.GetItemCount()))
    else:
        # Is item_indices a list or a single item?
        if not util_misc.is_iterable(item_indices):
            # Bartender, make it a tuple!
            item_indices = (item_indices, )

    for item_index in item_indices:
        if is_listbox:
            if select:
                list_ctrl.Select(item_index)
            else:
                list_ctrl.Deselect(item_index)
        else:
            list_ctrl.Select(item_index, select)


def send_close_to_active_tab_OLD_wx300(notebook):
    """
    This was original code that worked with wxPython 2.8.x that (maybe) did not
    have a CallAfter() thing, or maybe because it was C++ wrapped that made it 
    OK?
    
    Given an AUI Notebook, simulates a click on the active tab's close icon 
    (which usually a circle with an X in it, but varies under GTK versus OS X 
    etc.). This allows us to programmatically (attempt to) close a tab exactly
    the same way the user does. That is all the same event calls are made.
    
    """
    # I close the tab by mimicing the event that gets sent when the
    # user clicks on the tab's close button. The event's id has to
    # be that of the AuiTabCtrl. I figured this out by capturing a
    # real tab click event and examining its properties.
    tab_control = None
    for tab_control in notebook.GetChildren():
        if isinstance(tab_control, aui.AuiTabCtrl):
            break

    event = aui.AuiNotebookEvent(aui.EVT_AUINOTEBOOK_BUTTON.typeId, tab_control.Id)
    event.SetInt(aui.AUI_BUTTON_CLOSE)
    event.SetEventObject(tab_control)
    notebook.GetEventHandler().ProcessEvent(event)




    
    
def send_close_to_active_tab(notebook):
    """
    Given an AUI Notebook, simulates a click on the active tab's close icon 
    (which usually a circle with an X in it, but varies under GTK versus OS X 
    etc.). This allows us to programmatically (attempt to) close a tab exactly
    the same way the user does. That is all the same event calls are made.
    
    As of update to vespa 0.10.0 using wxPython 4.0.4, we had to amend this.
    
    Derived from aui.OnTabButton(self, event) which handles the 
    EVT_AUINOTEBOOK_BUTTON event for class AuiNotebook in wx.lib.aui.auibook.py
    I had to do this in here because it was calling a CallAfter() event to do
    the actual Tab destroy and this bollixed up the code that follows.
    
    
    """
    # I close the tab by mimicing the event that gets sent when the
    # user clicks on the tab's close button. The event's id has to
    # be that of the AuiTabCtrl. I figured this out by capturing a
    # real tab click event and examining its properties.
    tabs = None
    for tabs in notebook.GetChildren():
        if isinstance(tabs, aui.AuiTabCtrl):
            break

#    event = aui.AuiNotebookEvent(aui.EVT_AUINOTEBOOK_BUTTON.typeId, tabs.Id)
#    event.SetInt(aui.AUI_BUTTON_CLOSE)
#    event.SetEventObject(tabs)

    # if the close button is to the right, use the active
    # page selection to determine which page to close
    selection = tabs.GetActivePage()

    if selection == -1 or not tabs.GetEnabled(selection):
        return

    close_wnd = tabs.GetWindowFromIdx(selection)

    # ask owner if it's ok to close the tab
    e = aui.AuiNotebookEvent(aui.wxEVT_COMMAND_AUINOTEBOOK_PAGE_CLOSE, notebook.GetId())
    idx = notebook._tabs.GetIdxFromWindow(close_wnd)
    e.SetSelection(idx)
    e.SetOldSelection(0)
    e.SetEventObject(notebook)
    notebook.GetEventHandler().ProcessEvent(e)
    if not e.IsAllowed():
        return

    if repr(close_wnd.__class__).find("AuiMDIChildFrame") >= 0:
        close_wnd.Close()

    else:
        main_idx = notebook._tabs.GetIdxFromWindow(close_wnd)
        notebook.DeletePage(main_idx)

    # notify owner that the tab has been closed
    e2 = aui.AuiNotebookEvent(aui.wxEVT_COMMAND_AUINOTEBOOK_PAGE_CLOSED, notebook.GetId())
    e2.SetSelection(idx)
    e2.SetEventObject(notebook)
    notebook.GetEventHandler().ProcessEvent(e2)

    if notebook.GetPageCount() == 0:
        mgr = notebook.GetAuiManager()
        win = mgr.GetManagedWindow()
        win.SendSizeEvent()






def set_font_to_monospace(control):
    """Surprisingly, this function sets the font of the given control to
    the default monospace font.
    """
    font = control.GetFont()
    point_size = font.GetPointSize()
    font = wx.Font(point_size, wx.FONTFAMILY_TELETYPE, wx.FONTSTYLE_NORMAL,
                   wx.FONTWEIGHT_NORMAL)
    control.SetFont(font)


def show_wx_inspector(window):
    from wx.lib.inspection import InspectionTool

    if not InspectionTool().initialized:
        InspectionTool().Init()

    # Find a widget to be selected in the tree.  Use either the one under the
    # cursor, if any, or this frame.
    selected_window = wx.FindWindowAtPointer()
    if not selected_window:
        selected_window = window
    InspectionTool().Show(selected_window, True)


def wrap_label(label, layout_control=None, fudge_factor=.97):
    """Given a label (wx.StaticText) control, wraps the text therein so
    that it fits inside its sizer.

    The layout_control param is typically the parent (e.g. a dialog, panel,
    etc.) of the label control. If layout_control is not None, this function
    will call layout_control.Layout() once it's done wrapping. This is
    important when the label has gotten taller or shorter (added/removed a
    line) as a result of wrapping.

    fudge_factor exists because sizing the label to the exact width of the
    sizer doesn't render correctly, so we size the label to the
    sizer's width * fudge_factor. Adjust if necessary.
    """
    # In wx, labels (wxStaticText) don't wrap automatically. They do, however,
    # have a .Wrap() method that requires an explicit pixel width. The correct
    # width to use is that of the sizer, but the sizer's width seems to
    # include it's left & right borders which are not meant to be drawn on.
    # So I get the sizer's width, shrink it a little and pass that to .Wrap().

    width = label.GetContainingSizer().GetSize().width

    # wx implements wrapping by adding newlines to the text. That's fine
    # if you only want to wrap once but it makes a mess if you try to
    # wrap the same text multiple times. Here we remove the newlines.
    label_text = label.GetLabel()
    label_text = util_misc.normalize_newlines(label_text, " ")
    label.SetLabel(label_text)

    label.Wrap(int(width * fudge_factor))

    if layout_control:
        layout_control.Layout()

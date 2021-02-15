# Python modules


# 3rd party modules
import wx

# Our modules
import vespa.common.util.misc as util_misc

# Some shortcut styles for message boxes
Q_OK_CANCEL = wx.ICON_QUESTION | wx.OK | wx.CANCEL
Q_YES_NO = wx.ICON_QUESTION | wx.YES_NO
I_OK = wx.ICON_INFORMATION | wx.OK
X_OK = wx.ICON_EXCLAMATION | wx.OK
E_OK = wx.ICON_ERROR | wx.OK


def message(message_text, title=None, style=wx.ICON_INFORMATION | wx.OK):
    """Displays a message box. 
    
    The title defaults to the title of the App's top level window.
    
    The param style should be a message style ORed with a button style. 

    Message styles are: wx.ICON_EXCLAMATION, wx.ICON_ERROR, wx.ICON_QUESTION
    and wx.ICON_INFORMATION.

    The button styles are: wx.OK, wx.CANCEL and wx.YES_NO. 
    You can add wx.YES_DEFAULT or wx.NO_DEFAULT to make Yes or No the
    default choice.

    NOTE: the return code is one of wx.OK, wx.CANCEL, wx.YES or wx.NO.
    
    This differs from wx.MessageDialog.ShowModal() which returns wx.ID_OK,
    wx.ID_CANCEL, etc.
    """
    if not title:
        # An app may call this to display a message box before it has 
        # created any windows, so we can't assume a top window always exists.
        top_window = wx.GetApp().GetTopWindow()
        title = top_window.GetTitle() if top_window else "Vespa"

    return wx.MessageBox(message_text, title, style)


def pickdir(message="", default_path=""):
    """Prompts the user to select a directory and returns the selected
    directory as a string. If the user hits cancel, the function returns an
    empty string.
    
    If the message param is empty (the default), the message is a standard
    one chosen by wx.

    Under OS X and Windows, this dialog remembers the directory where it was 
    last invoked, and if the default_path param is empty, it will start in 
    that same directory. Under GTK a reasonable default is used.
    """    
    if not message:
        message = wx.DirSelectorPromptStr
    
    # Under OS X and Windows, wx automatically remembers the last dir that 
    # the file selector used. Ubuntu forgets, however. This might be a 
    # function of Ubuntu, Linux, GTK, I don't know.
    if (not default_path) and ("__WXGTK__" in wx.PlatformInfo):
        default_path = util_misc.get_documents_dir()
        
    return wx.DirSelector(message, default_path)


def pickfile(message="", filetype_filter="All files (*.*)|*.*", 
             enforce_filter=False, default_path="", default_file="",
             multiple=False, new_file=False):
    """
    Prompts the user to select a single file to be opened and returns the 
    selected file as a string. If the user hits cancel, the function returns 
    an empty string.
    
    If the message param is empty, the message is a standard one chosen by wx
    (which is "Select a file" under OS X, Windows & Linux).

    If you supply a filetype_filter, it should look something like this:
       "JPEG files (*.jpg)|*.jpg"
    One can supply multiple file types like so:
       "JPEG files (*.jpg)|*.jpg|Bitmap files (*.bmp)|*.bmp"
       
    If you supply a filetype_filter and enforce_filter is False (the default), 
    this function will append "|All files (*.*)|*.*" to the filter to allow 
    users to select a file with any extension. Otherwise, the user is limited
    to selecting a file that has one of the extensions in the filter.
    
    Under OS X and Windows, this dialog remembers the directory where it was 
    last invoked, and if the default_path param is empty, it will start in 
    that same directory. Under GTK a reasonable default is used.

    When multiple=False (the default), the function returns a string ("" if 
    the user chose cancel). When multiple=True, a sorted list (possibly empty)
    is returned. 
    
    When new_file=False (the default), the user must select a file that 
    already exists. When new_file=True, the user may type in the name of a
    file that does not yet exist.
    """  
    # I tried using the wx.CHANGE_DIR flag so that wx would handle remembering
    # where the last file was opened, but it didn't work on Linux so we had to 
    # roll our own solution.
    flags = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
    if new_file:
        flags = wx.FD_OPEN
    
    if multiple:
        return_value = [ ]
        flags |= wx.FD_MULTIPLE
    else:
        return_value = ""
    
    if not message:
        message = wx.FileSelectorPromptStr
    
    if not enforce_filter and not filetype_filter.endswith("*.*"):
        # Ensure the user can select all file types.
        if not filetype_filter.endswith("|"):
            filetype_filter += "|"
        filetype_filter += "All files (*.*)|*.*"
        
    # Under OS X and Windows, wx automatically remembers the last dir that 
    # the file selector used. Ubuntu forgets, however. This might be a 
    # function of Ubuntu, Linux, GTK, I don't know.
    if (not default_path) and ("__WXGTK__" in wx.PlatformInfo):
        default_path = util_misc.get_documents_dir()
        
    dialog = wx.FileDialog(wx.GetApp().GetTopWindow(), message, default_path, 
                           default_file, filetype_filter, flags)

    if dialog.ShowModal() == wx.ID_OK:
        if multiple:
            return_value = sorted(dialog.GetPaths())
        else:
            return_value = dialog.GetPath()

    return return_value
    
             
def save_as(message="Save file", filetype_filter="", default_path="", 
             default_filename=""):
    """Prompts the user to indicate the name of a file to be saved and 
    returns the name as a string. If the user hits cancel, the function 
    returns an empty string.
    
    The filename_filter param is only a suggestion. The user is free to 
    choose any name and extension she chooses. If you supply a 
    filename_filter, it should look something like this:
       "JPEG files (*.jpg)|*.jpg"
       
    Under OS X and Windows, this dialog remembers the directory where it was 
    last invoked, and if the default_path param is empty, it will start in 
    that same directory. Under GTK a reasonable default is used.    
    """    
    flags = wx.FD_SAVE | wx.FD_CHANGE_DIR | wx.FD_OVERWRITE_PROMPT 
    
    # Under OS X and Windows, wx automatically remembers the last dir that 
    # the file selector used. Ubuntu forgets, however. This might be a 
    # function of Ubuntu, Linux, GTK, I don't know.
    if (not default_path) and ("__WXGTK__" in wx.PlatformInfo):
        default_path = util_misc.get_documents_dir()
        
    # Under Windows, OS X & Ubuntu 9, with wxPython 2.8, the 
    # default_extension parameter to wx.FileSelector() is ignored. 
    default_extension = ""
    return wx.FileSelector(message, default_path, default_filename, 
                           default_extension, filetype_filter, flags)

